"""One source seen on two tetrodes: a post-hoc duplicate rule for per-tetrode matching, scored on the
hybrid without re-running a carry.

Per window, for each pair of units on different tetrodes with 200 or more accepted spikes: the counts
of one unit's spikes at each lag from the other's within 2 ms, the excess within 0.5 ms over the mean
count at 1 to 2 ms (each unit's spike time sits on its own trough, so one source's copies can be
offset by a few samples), the lag of the peak, the share of the excess within one sample of the peak, and
the correlation of the two fits' snr across the pairs at the peak (within one sample). A pair is
flagged when the excess is ``--excess`` or more of the smaller unit's spikes, the share at the peak
is ``--sharp`` or more and the correlation ``--corr`` or more: the two units then fire together to
the sample with amplitudes that scale together, which one source on two tetrodes does and two cells
firing together need not. Within a flagged pair, of each two fits at the peak lag (within one sample)
the one with the lower snr is removed.

Ground-truth units are injected at times independent of the real units', so the hybrid measures
what the rule removes from units that do not share a source; it cannot measure what the rule costs
two real cells that fire together.

Prints the flagged pairs and, over the arm's windows, pooled recall and precision of the small (trough
6 MAD or less) and large ground-truth units before and after (window means, the units that stop
after window 5 left out from window 6 on), per ground-truth unit, and the real units that lose the
most.

Run: .venv/bin/python investigations/carry-bakeoff/xtet_dedup.py ARM [ARM ...] [--real WINDOW]
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402

D = 60
CORE = 15


def lags(fu: np.ndarray, fv: np.ndarray) -> np.ndarray:
    """Every lag of sorted ``fv`` from ``fu`` within ``D`` samples."""
    lo, hi = np.searchsorted(fv, fu - D), np.searchsorted(fv, fu + D + 1)
    n = hi - lo
    if n.sum() == 0:
        return np.zeros(0, np.int64)
    start = np.repeat(np.cumsum(n) - n, n)
    j = np.repeat(lo, n) + np.arange(n.sum()) - start
    return fv[j] - np.repeat(fu, n)


def flagged_pairs(fr, un, snr, grp, *, min_n=200, excess_min=0.05, sharp_min=0.5, corr_min=0.3):
    units = [u for u, n in zip(*np.unique(un, return_counts=True)) if n >= min_n and int(u) in grp]
    sorted_ = {}
    for u in units:
        o = np.argsort(fr[un == u], kind="stable")
        sorted_[u] = (fr[un == u][o], snr[un == u][o], np.flatnonzero(un == u)[o])
    out = []
    for i, u in enumerate(units):
        fu, su, iu = sorted_[u]
        for v in units[i + 1:]:
            if grp[int(u)] == grp[int(v)]:
                continue
            fv, sv, iv = sorted_[v]
            d = lags(fu, fv)
            if d.size == 0:
                continue
            c = np.bincount(d + D, minlength=2 * D + 1).astype(np.float64)
            base = np.r_[c[:D - 30], c[D + 31:]].mean()
            core = c[D - CORE:D + CORE + 1] - base
            excess = core.sum()
            if excess <= 0 or excess / min(fu.size, fv.size) < excess_min:
                continue
            pk = int(core.argmax())
            sharp = core[max(pk - 1, 0):pk + 2].sum() / excess
            lag = pk - CORE
            j = np.clip(np.searchsorted(fv, fu + lag), 1, fv.size - 1)
            pick = np.where(np.abs(fv[j - 1] - fu - lag) <= np.abs(fv[j] - fu - lag), j - 1, j)
            m = np.abs(fv[pick] - fu - lag) <= 1
            if m.sum() < 30:
                continue
            r = float(np.corrcoef(su[m], sv[pick[m]])[0, 1])
            row = {"u": int(u), "v": int(v), "tet_u": grp[int(u)] + 1, "tet_v": grp[int(v)] + 1,
                   "n_u": fu.size, "n_v": fv.size, "excess": excess / min(fu.size, fv.size), "lag": lag,
                   "sharp": sharp, "r": r, "pairs": int(m.sum()),
                   "flag": bool(sharp >= sharp_min and r >= corr_min)}
            if row["flag"]:
                lose_u = su[m] < sv[pick[m]]
                row["drop"] = np.r_[iu[m][lose_u], iv[pick[m]][~lose_u]]
            out.append(row)
    return out


def dedup(sp, grp, **kw):
    acc = np.flatnonzero(sp["reason"] == 0)
    pairs = flagged_pairs(sp["frame"][acc], sp["unit"][acc], sp["snr"][acc], grp, **kw)
    keep = np.ones(sp["frame"].size, bool)
    for p in pairs:
        if p["flag"]:
            keep[acc[p["drop"]]] = False
    return (sp["reason"] == 0) & keep, pairs


def score_gt(fr, un, gt, w):
    a, b = w * bk.WIN, (w + 1) * bk.WIN
    mw = (gt["frame"] >= a) & (gt["frame"] < b)
    gfr, gun = gt["frame"][mw], gt["unit"][mw]
    rows = []
    for i in range(len(gt["group"])):
        mine = np.sort(fr[un == 1000 + i])
        gi = gfr[gun == i]
        rows.append({"window": w, "gt": i, "n_true": gi.size, "n_unit": mine.size,
                     "found": int(near(gi, mine).sum()), "hits": int(near(mine, np.sort(gi)).sum())})
    return rows


def run_hybrid(name, kw):
    d = bk.BAKE / "arms" / name
    gt = np.load(bk.BAKE / "gt.npz")
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    small = {p["gt"] for p in plan if p["trough_mad"] <= 6}
    leaves = {p["gt"] for p in plan if p["kind"] == "leaves"}
    b0 = bk.Bank.load(bk.BAKE / "bank_start.npz")
    grp = dict(zip(b0.ids.tolist(), b0.group.tolist()))
    rows, real, flags = [], [], []
    w = 0
    print(f"== {name}")
    while (d / f"spikes_w{w:02d}.npz").exists():
        sp = dict(np.load(d / f"spikes_w{w:02d}.npz"))
        keep, pairs = dedup(sp, grp, **kw)
        acc = sp["reason"] == 0
        flags += [{k: v for k, v in p.items() if k != "drop"} | {"window": w, "removed": p.get("drop", []).__len__()}
                  for p in pairs]
        for tag, m in (("before", acc), ("after", keep)):
            for r in score_gt(sp["frame"][m], sp["unit"][m], gt, w):
                rows.append({**r, "rule": tag})
        for u in np.unique(sp["unit"][acc]):
            if u < 1000:
                real.append({"window": w, "unit": int(u), "tetrode": grp[int(u)] + 1,
                             "before": int((acc & (sp["unit"] == u)).sum()),
                             "after": int((keep & (sp["unit"] == u)).sum())})
        w += 1
    fl = pd.DataFrame(flags)
    fl["pair"] = [f"u{r.u} TT{r.tet_u} / u{r.v} TT{r.tet_v}" for r in fl.itertuples()]
    print("pairs with excess 0.05 or more in a window: windows flagged / windows seen; medians of excess, lag, "
          "share at the peak, snr correlation; spikes removed")
    t = fl.groupby("pair").agg(flagged=("flag", "sum"), seen=("flag", "size"), excess=("excess", "median"),
                               lag=("lag", "median"), sharp=("sharp", "median"), r=("r", "median"),
                               removed=("removed", "sum"))
    print(t.sort_values(["flagged", "removed"], ascending=False).round(2).to_string())
    df = pd.DataFrame(rows)
    df = df[~df["gt"].isin(leaves) | (df.window < 6)]
    df["size"] = np.where(df["gt"].isin(small), "small", "large")
    pooled = df.groupby(["rule", "size", "window"]).agg(found=("found", "sum"), n_true=("n_true", "sum"),
                                                        hits=("hits", "sum"), n_unit=("n_unit", "sum"))
    pooled["recall"] = pooled.found / pooled.n_true
    pooled["precision"] = pooled.hits / pooled.n_unit.clip(lower=1)
    print(pooled.groupby(["size", "rule"])[["recall", "precision"]].mean().round(3).to_string())
    per = df.groupby(["gt", "rule"]).agg(found=("found", "sum"), n_true=("n_true", "sum"), hits=("hits", "sum"),
                                         n_unit=("n_unit", "sum"))
    per["recall"] = per.found / per.n_true
    per["precision"] = per.hits / per.n_unit.clip(lower=1)
    t = per[["recall", "precision"]].unstack("rule").round(3)
    print("per ground-truth unit, all windows pooled")
    print(t[(t["precision"]["after"] != t["precision"]["before"]) | (t["recall"]["after"] != t["recall"]["before"])]
          .to_string())
    rl = pd.DataFrame(real).groupby(["unit", "tetrode"])[["before", "after"]].sum()
    rl["lost"] = 1 - rl.after / rl.before
    print("real units losing 5 % or more")
    print(rl[rl.lost >= 0.05].sort_values("lost", ascending=False).round(3).to_string())


def run_real(name, w, kw):
    d = bk.BAKE_ROOT / "real48" / name
    sp = dict(np.load(d / f"spikes_w{w:02d}.npz"))
    b = bk.Bank.load(d / f"bank_w{w:02d}.npz")
    grp = dict(zip(b.ids.tolist(), b.group.tolist()))
    keep, pairs = dedup(sp, grp, **kw)
    acc = sp["reason"] == 0
    print(f"== {name}, window {w}: flagged pairs")
    for p in sorted(pairs, key=lambda p: -p["excess"]):
        if p["flag"]:
            print(f"  u{p['u']} TT{p['tet_u']} / u{p['v']} TT{p['tet_v']}: {p['excess']:.2f}, {p['lag']:+d}, "
                  f"{p['sharp']:.2f}, {p['r']:.2f}, {p['pairs']}; removed {p['drop'].size}")
    print("not flagged with excess 0.05 or more: "
          + "; ".join(f"u{p['u']}/u{p['v']} sharp {p['sharp']:.2f} r {p['r']:.2f}" for p in pairs if not p["flag"]))
    rows = [{"unit": int(u), "tetrode": grp[int(u)] + 1, "before": int((acc & (sp["unit"] == u)).sum()),
             "after": int((keep & (sp["unit"] == u)).sum())} for u in np.unique(sp["unit"][acc])]
    rl = pd.DataFrame(rows).set_index(["unit", "tetrode"])
    rl["lost"] = 1 - rl.after / rl.before
    print(rl[rl.lost >= 0.02].sort_values("lost", ascending=False).round(3).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arms", nargs="+")
    ap.add_argument("--real", type=int, default=None, help="score a real48 arm's window instead")
    ap.add_argument("--excess", type=float, default=0.05)
    ap.add_argument("--sharp", type=float, default=0.5)
    ap.add_argument("--corr", type=float, default=0.3)
    a = ap.parse_args()
    kw = {"excess_min": a.excess, "sharp_min": a.sharp, "corr_min": a.corr}
    pd.set_option("display.width", 200)
    for name in a.arms:
        if a.real is None:
            run_hybrid(name, kw)
        else:
            run_real(name, a.real, kw)


if __name__ == "__main__":
    main()
