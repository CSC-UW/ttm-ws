"""``ttm-bad-periods``: transient bad periods per tetrode for one experiment's recording.

Scans the recording as the sorter sees it (bandpass 300 to 6000 Hz plus the
aligned median reference, the ``preprocess_for_sorting`` defaults) in 0.1 s
bins, builds 1 s features, and writes into ``<session>/bad_periods/``:

  bin_stats_0p1s.npz      per-bin sums, reused when made from the same input
                          (source, reference, channels left out) unless ``--rescan``
  bad_periods_params.json the input, the thresholds and the package versions
  features_1s.npz         rms, correlations and their running baselines
  bad_periods.csv         one row per tetrode and interval (``ttm.bad_periods``)
  bad_periods_summary.csv intervals, flagged seconds and percent per tetrode
  bad_periods.png         the intervals over the session, one row per tetrode

``--binary`` reads an already materialized preprocessed binary instead of
filtering the store, which is much faster (``ttm-sort binary`` writes one).
``--cmr global`` scans with the plain median of all channels instead, as the
tables made before 2026-11-03 were. The manifest's ``badChannels`` and
every channel of its ``badTetrodes`` are removed before the reference (and
from a ``--binary`` that still holds them) unless ``--keep-bad-channels``.
Detector thresholds are the library defaults unless given.
"""

from __future__ import annotations

import argparse
import importlib.metadata as md
import json
import time

import numpy as np
import pandas as pd

from ttm import bad_periods as bp
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered


def _bad_channels(exp, a) -> list:
    import spikeinterface as si

    from ttm.sorting import bad_channel_ids

    if a.keep_bad_channels:
        return []
    return bad_channel_ids(si.read_zarr(str(exp.store_path())), exp.badChannels, exp.badTetrodes)


def _recording(exp, binary, drop, cmr="global_aligned"):
    import spikeinterface as si

    if binary is not None:
        rec = si.load(binary)
        present = [c for c in rec.channel_ids if str(c) in {str(d) for d in drop}]
        return rec.remove_channels(present) if present else rec
    from ttm.sorting import preprocess_for_sorting

    return preprocess_for_sorting(si.read_zarr(str(exp.store_path())), cmr=cmr, bad_channel_ids=drop)


def scan_key(exp, binary, cmr, drop) -> dict:
    """What the scan was made from; a cached scan is reused only when this matches."""
    source = str(binary) if binary is not None else str(exp.store_path())
    return {"source": source, "cmr": "binary" if binary is not None else cmr,
            "dropped": sorted(str(c) for c in drop)}


def cached_key(stats: dict) -> dict | None:
    if "scan_key" not in stats:
        return None
    return json.loads(str(stats["scan_key"]))


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
    ap.add_argument("--cmr", choices=["global_aligned", "global"], default="global_aligned",
                    help="reference applied to the store (ignored with --binary)")
    ap.add_argument("--n-jobs", type=int, default=48)
    ap.add_argument("--wire-loud", type=float, default=None)
    ap.add_argument("--common-corr-excess", type=float, default=None)
    ap.add_argument("--keep-bad-channels", action="store_true",
                    help="keep the manifest's badChannels and badTetrodes in the reference and the scan")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    drop = _bad_channels(exp, a)
    print(f"bad channels left out: {[str(c) for c in drop] or 'none'}")
    out = exp.bad_periods_dir
    out.mkdir(parents=True, exist_ok=True)

    stats_path = out / "bin_stats_0p1s.npz"
    key = scan_key(exp, a.binary, a.cmr, drop)
    stats = None
    if stats_path.exists() and not a.rescan:
        cached = dict(np.load(stats_path))
        if cached_key(cached) == key:
            stats = cached
            print(f"scan: reuse {stats_path.name}")
        else:
            print(f"scan: {stats_path.name} was made from {cached_key(cached)}, not {key}; rescanning")
    if stats is None:
        rec = _recording(exp, a.binary, drop, a.cmr)
        t0 = time.time()
        stats = bp.bin_stats(rec, bin_s=0.1, n_jobs=a.n_jobs)
        np.savez(stats_path, **stats, scan_key=json.dumps(key))
        print(f"scan: {stats['ss'].shape[0]} bins of 0.1 s in {time.time() - t0:.0f} s -> {stats_path.name}")
    stats.pop("scan_key", None)

    feats = bp.bin_features(stats, bin_s=1.0)
    base = bp.baselines(feats)
    np.savez(out / "features_1s.npz", **feats, **base)
    names = {}
    for g, name in zip(feats["groups"], _group_names(a, exp, feats, drop)):
        names[int(g)] = name
    kw = {k: v for k, v in dict(wire_loud=a.wire_loud, common_corr_excess=a.common_corr_excess).items()
          if v is not None}
    table = bp.detect_bad_periods(feats, base, group_names=names, **kw)
    table.to_csv(out / "bad_periods.csv", index=False)
    versions = {}
    for name in ("ttm", "spikeinterface", "numpy"):
        try:
            versions[name] = md.version(name)
        except md.PackageNotFoundError:
            versions[name] = None
    (out / "bad_periods_params.json").write_text(json.dumps({
        "scan": key, "thresholds": kw or "library defaults", "versions": versions,
        "written": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2) + "\n")

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


def _group_names(a, exp, feats, drop=()) -> list[str]:
    """Tetrode names per channel from the store's ``tetrode`` property, else ``TT<group + 1>``."""
    try:
        import spikeinterface as si

        rec = si.load(a.binary) if a.binary is not None else si.read_zarr(str(exp.store_path()))
        present = [c for c in rec.channel_ids if str(c) in {str(d) for d in drop}]
        if present:
            rec = rec.remove_channels(present)
        if "tetrode" in rec.get_property_keys():
            return [str(v) for v in rec.get_property("tetrode")]
    except Exception:
        pass
    return [f"TT{int(g) + 1}" for g in feats["groups"]]


if __name__ == "__main__":
    main()
