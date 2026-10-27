"""``ttm-bad-periods``: transient bad periods per tetrode for one experiment's recording.

Scans the recording as the sorter sees it (bandpass 300 to 6000 Hz plus global
common median reference) in 0.1 s bins, builds 1 s features, and writes into
``<session>/bad_periods/``:

  bin_stats_0p1s.npz      per-bin sums (reused unless ``--rescan``)
  features_1s.npz         rms, correlations and their running baselines
  bad_periods.csv         one row per tetrode and interval (``ttm.bad_periods``)
  bad_periods_summary.csv intervals, flagged seconds and percent per tetrode
  bad_periods.png         the intervals over the session, one row per tetrode

``--binary`` reads an already materialized preprocessed binary instead of
filtering the store, which is much faster. Detector thresholds are the
library defaults unless given.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from ttm import bad_periods as bp
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered


def _recording(exp, binary):
    import spikeinterface as si

    if binary is not None:
        return si.load(binary)
    from ttm.sorting import preprocess_for_sorting

    return preprocess_for_sorting(si.read_zarr(str(exp.store_path())), cmr="global")


def _plot(table: pd.DataFrame, names: dict, duration_s: float, out_png, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"loud": "C3", "common": "C0", "global": "k"}
    gids = sorted(names)
    fig, ax = plt.subplots(figsize=(14, 0.3 * len(gids) + 1.2))
    for row, g in enumerate(gids):
        sub = table[table["group"] == g]
        for kind, c in colors.items():
            s = sub[sub["kind"].str.split("+").str[0] == kind]
            ax.broken_barh(list(zip(s["start_s"] / 3600, (s["end_s"] - s["start_s"]) / 3600)),
                           (row - 0.4, 0.8), color=c, lw=0)
    ax.set_yticks(range(len(gids)), [names[g] for g in gids])
    ax.invert_yaxis()
    ax.set_xlim(0, duration_s / 3600)
    ax.set_xlabel("session time (h)   [red: loud wire, blue: shared signal, black: most tetrodes at once]")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


def main() -> None:
    line_buffered()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--binary", default=None, help="materialized preprocessed binary to scan instead of the store")
    ap.add_argument("--rescan", action="store_true", help="recompute bin_stats_0p1s.npz")
    ap.add_argument("--n-jobs", type=int, default=48)
    ap.add_argument("--wire-loud", type=float, default=None)
    ap.add_argument("--common-corr-excess", type=float, default=None)
    a = ap.parse_args()
    exp = experiment_from_args(a)
    out = exp.bad_periods_dir
    out.mkdir(parents=True, exist_ok=True)

    stats_path = out / "bin_stats_0p1s.npz"
    if stats_path.exists() and not a.rescan:
        stats = dict(np.load(stats_path))
        print(f"scan: reuse {stats_path.name}")
    else:
        rec = _recording(exp, a.binary)
        t0 = time.time()
        stats = bp.bin_stats(rec, bin_s=0.1, n_jobs=a.n_jobs)
        np.savez(stats_path, **stats)
        print(f"scan: {stats['ss'].shape[0]} bins of 0.1 s in {time.time() - t0:.0f} s -> {stats_path.name}")

    feats = bp.bin_features(stats, bin_s=1.0)
    base = bp.baselines(feats)
    np.savez(out / "features_1s.npz", **feats, **base)
    names = {}
    for g, name in zip(feats["groups"], _group_names(a, exp, feats)):
        names[int(g)] = name
    kw = {k: v for k, v in dict(wire_loud=a.wire_loud, common_corr_excess=a.common_corr_excess).items()
          if v is not None}
    table = bp.detect_bad_periods(feats, base, group_names=names, **kw)
    table.to_csv(out / "bad_periods.csv", index=False)

    duration = float(feats["t"][-1] + feats["bin_s"])
    valid = bp.valid_seconds(table, duration, feats["group_ids"])
    summary = pd.DataFrame({
        "group": [int(g) for g in feats["group_ids"]],
        "tetrode": [names[int(g)] for g in feats["group_ids"]],
        "intervals": [int((table["group"] == g).sum()) for g in feats["group_ids"]],
        "flagged_s": [round(duration - valid[int(g)], 1) for g in feats["group_ids"]],
        "flagged_pct": [round(100 * (duration - valid[int(g)]) / duration, 3) for g in feats["group_ids"]],
        "longest_s": [float((table[table["group"] == g].eval("end_s - start_s")).max()) if (table["group"] == g).any()
                      else 0.0 for g in feats["group_ids"]],
    })
    summary.to_csv(out / "bad_periods_summary.csv", index=False)
    _plot(table, names, duration, out / "bad_periods.png", f"{exp.subject}/{exp.experiment} bad periods")
    print(summary.to_string(index=False))
    print(f"{len(table)} intervals -> {out / 'bad_periods.csv'}")


def _group_names(a, exp, feats) -> list[str]:
    """Tetrode names per channel from the store's ``tetrode`` property, else ``TT<group + 1>``."""
    try:
        import spikeinterface as si

        rec = si.load(a.binary) if a.binary is not None else si.read_zarr(str(exp.store_path()))
        if "tetrode" in rec.get_property_keys():
            return [str(v) for v in rec.get_property("tetrode")]
    except Exception:
        pass
    return [f"TT{int(g) + 1}" for g in feats["groups"]]


if __name__ == "__main__":
    main()
