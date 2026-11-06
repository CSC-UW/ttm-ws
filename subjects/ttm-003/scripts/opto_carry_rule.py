"""What the carry's in-pulse rule changes: two carry runs of the same windows, with and without the pulse table.

The runs (``ttm-sort carry --run A --pulses ...`` and ``ttm-sort carry --run B --bank A/bank0.npz --whitening
A/whitening.npz``, same ``--windows``) start from one bank and whitening, so their first window has the same
fits and differs only by the rule; later windows differ also through the template updates. Per unit and
window, the accepted spikes in each run: in all, 1 to 10 ms after onset, 15 ms after onset to the offset, and
in the 100 ms after the offset, and whether the unit is driven. Where the fits are the same, the spikes the
rule adds and drops, with their median amplitude, matched-filter snr and cosine, against those of the spikes
both runs accept, by position relative to the pulse.

With ``--stretch NAME`` (an ``opto_units.py`` stretch), each driven unit of the run with the table is matched to
the unit of that stretch's sort that shares most of its spikes within 10 samples, over the windows both
cover.

Writes ``results/opto_carry_rule_<E>.csv`` (one row per unit and window),
``results/opto_carry_rule_<E>_added.csv`` (the comparison of added and accepted spikes) and, with
``--stretch``, ``results/opto_carry_rule_<E>_match.csv``.

Usage: python subjects/ttm-003/scripts/opto_carry_rule.py --subject TTM-003 --experiment TTM-003-48H-2
           --with-run pulsetest --without-run nopulsetest --pulses subjects/ttm-003/results/laser_pulses_<session>.csv
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from ttm import carry as cr
from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
WINDOW_S = 1800.0
POST_S = 0.1


def latency(f: np.ndarray, on: np.ndarray, off: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    j = np.searchsorted(on, f, "right") - 1
    ok = j >= 0
    lat = np.full(f.size, np.inf)
    plen = np.full(f.size, np.nan)
    lat[ok] = (f[ok] - on[j[ok]]) / FS
    plen[ok] = (off[j[ok]] - on[j[ok]]) / FS
    return lat, plen


def zones(lat: np.ndarray, plen: np.ndarray) -> dict[str, np.ndarray]:
    return {"onset_1_10ms": (lat >= 0.001) & (lat < 0.010), "pulse_15ms_off": (lat >= 0.015) & (lat < plen),
            "post_0_100ms": (lat >= plen) & (lat < plen + POST_S),
            "outside": ~((lat >= 0) & (lat < plen + POST_S))}


def match_units(exp, stretch: str, d_with: pathlib.Path, df: pd.DataFrame) -> None:
    info = json.loads((RESULTS / f"opto_units_{exp.experiment}_{stretch}.json").read_text())
    fa = int(info["start_h"] * 3600 * FS)
    fb = int(info["end_h"] * 3600 * FS)
    sorting = cr.sort_window(None, exp.experiment_dir / "opto_units" / stretch, cr.CarryParams())
    kept = set(pd.read_csv(RESULTS / f"opto_units_{exp.experiment}_{stretch}.csv").unit)
    rows = []
    for u in sorted(set(df.loc[df.driven, "unit"])):
        ws = sorted(df.loc[df.unit == u, "window"])
        parts = []
        for w in ws:
            z = np.load(d_with / f"spikes_w{w:02d}.npz")
            parts.append(z["frame"][(z["unit"] == u) & (z["reason"] == 0)])
        fr = np.concatenate(parts)
        lo, hi = max(fa, int(ws[0] * WINDOW_S * FS)), min(fb, int((ws[-1] + 1) * WINDOW_S * FS))
        fr = fr[(fr >= lo) & (fr < hi)]
        for v in sorting.unit_ids:
            if int(v) not in kept:
                continue
            g = sorting.get_unit_spike_train(v).astype(np.int64) + fa
            g = g[(g >= lo) & (g < hi)]
            if g.size == 0:
                continue
            k = np.searchsorted(g, fr)
            near = np.minimum(np.abs(fr - g[np.clip(k - 1, 0, g.size - 1)]), np.abs(g[np.clip(k, 0, g.size - 1)] - fr)) <= 10
            shared = int(near.sum())
            if shared:
                rows.append({"carry_unit": int(u), "stretch": stretch, "sort_unit": int(v), "carry_spikes": int(fr.size),
                             "sort_spikes": int(g.size), "shared": shared, "share_of_carry": shared / fr.size,
                             "share_of_sort": shared / g.size})
    m = pd.DataFrame(rows).sort_values(["carry_unit", "shared"], ascending=[True, False])
    m.round(4).to_csv(RESULTS / f"opto_carry_rule_{exp.experiment}_match.csv", index=False)
    print(m.groupby("carry_unit").head(3).round(3).to_string(index=False))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--with-run", required=True)
    ap.add_argument("--without-run", required=True)
    ap.add_argument("--pulses", required=True)
    ap.add_argument("--stretch", default=None, help="an opto_units.py stretch to match the driven units to")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    d_with = exp.sorting_dir(a.with_run, subdir="sortings_carry") / "carry"
    d_without = exp.sorting_dir(a.without_run, subdir="sortings_carry") / "carry"
    p = pd.read_csv(a.pulses).sort_values("onset_s")
    on, off = (p.onset_s.to_numpy() * FS).astype(np.int64), (p.offset_s.to_numpy() * FS).astype(np.int64)
    windows = sorted(int(f.stem[-2:]) for f in d_with.glob("spikes_w*.npz") if (d_without / f.name).exists())
    rows, added = [], []
    for w in windows:
        sp, sn = np.load(d_with / f"spikes_w{w:02d}.npz"), np.load(d_without / f"spikes_w{w:02d}.npz")
        driven = np.load(d_with / f"bank_w{w:02d}.npz")["driven"].tolist()
        same = np.array_equal(sp["frame"], sn["frame"]) and np.array_equal(sp["unit"], sn["unit"])
        a0, a1 = w * WINDOW_S * FS, (w + 1) * WINDOW_S * FS
        n_pulses = int(((on >= a0) & (on < a1)).sum())
        acc_w, acc_wo = sp["reason"] == 0, sn["reason"] == 0
        print(f"w{w:02d}: {n_pulses} pulses, driven {driven}, same fits {same}, "
              f"accepted {acc_w.sum():,} with the table and {acc_wo.sum():,} without")
        for u in sorted(set(sp["unit"][acc_w].tolist()) | set(sn["unit"][acc_wo].tolist())):
            r = {"window": w, "unit": int(u), "driven": u in driven, "pulses": n_pulses, "same_fits": same}
            for tag, s, acc in (("with", sp, acc_w), ("without", sn, acc_wo)):
                f = s["frame"][acc & (s["unit"] == u)]
                z = zones(*latency(f, on, off))
                r[f"n_{tag}"] = int(f.size)
                r.update({f"{k}_{tag}": int(v.sum()) for k, v in z.items() if k != "outside"})
            if same:
                m = sp["unit"] == u
                r["added"] = int((m & acc_w & ~acc_wo).sum())
                r["dropped"] = int((m & ~acc_w & acc_wo).sum())
                if u in driven:
                    z = zones(*latency(sp["frame"], on, off))
                    for kind, sel in (("both", m & acc_w & acc_wo), ("added", m & acc_w & ~acc_wo)):
                        for k, zm in z.items():
                            q = sel & zm
                            if q.any():
                                added.append({"window": w, "unit": int(u), "spikes": kind, "zone": k, "n": int(q.sum()),
                                              "a_median": float(np.median(sp["a"][q])),
                                              "snr_median": float(np.median(sp["snr"][q])),
                                              "cos_pl_median": float(np.median(sp["cos_pl"][q]))})
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"opto_carry_rule_{exp.experiment}.csv", index=False)
    ad = pd.DataFrame(added)
    ad.round(3).to_csv(RESULTS / f"opto_carry_rule_{exp.experiment}_added.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(df[df.driven].to_string(index=False))
    print(ad.round(2).to_string(index=False))
    if a.stretch:
        match_units(exp, a.stretch, d_with, df)
    o = df[~df.driven & (df.n_without >= 500)].copy()
    o["change"] = (o.n_with - o.n_without) / o.n_without
    big = o[o.n_without >= 10000]
    print(f"other units with 500 spikes or more: change {o.change.min():+.4f} to {o.change.max():+.4f}; "
          f"with 10,000 or more: largest |change| {big.change.abs().max():.4f}")


if __name__ == "__main__":
    main()
