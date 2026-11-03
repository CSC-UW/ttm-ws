"""What the per-tetrode carry's false positives on ground-truth units are, against the joint carry's.

On the aligned-reference hybrid with the re-sorted seed bank (``BAKE_REFERENCE=aligned_seed``), arms
``W_s45_c50_lnb_b2`` (per tetrode) and ``W_joint_masked_s45_c50_lnb_b2`` (joint), windows 0 to 11
(the units that stop after window 5 left out from window 6). Each accepted spike of a ground-truth
unit is a hit when a ground-truth spike of the unit is within 10 samples. A false positive is
classed by what is near it: a ground-truth spike of the unit 11 to 45 samples away (a second fit of
the same spike), a ground-truth spike of another ground-truth unit of the tetrode within 10 samples,
an accepted spike of a real unit of the tetrode within 10 samples (in the same arm), or none of
these. Also whether the other arm accepted a spike of the unit within 10 samples, and whether the
same arm accepted a larger spike (by snr) on another tetrode within 3 ms. Per arm and class the
medians of snr, cos_pl, the fit amplitude and the peeled trough, against the hits'.

Writes ``results/fp_anatomy.csv`` (per unit, arm and class: count and medians).

Run: BAKE_REFERENCE=aligned_seed .venv/bin/python investigations/carry-bakeoff/fp_anatomy.py
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

HERE = pathlib.Path(__file__).parent
ARMS = {"per tetrode": "W_s45_c50_lnb_b2", "joint": "W_joint_masked_s45_c50_lnb_b2"}
TOL = 10
XT = int(0.003 * bk.FS)


def near(ref, other, lo=0, hi=TOL):
    """Whether each of ``ref`` has an element of sorted ``other`` at a distance in [lo, hi]."""
    if other.size == 0:
        return np.zeros(ref.size, bool)
    out = np.zeros(ref.size, bool)
    for sgn in (-1, 1):
        t = ref + sgn * lo
        i = np.searchsorted(other, t - hi if sgn < 0 else t, side="left")
        j = np.searchsorted(other, t if sgn < 0 else t + hi, side="right")
        out |= j > i
    return out


def larger_elsewhere(fr, snr, g, all_fr, all_snr, all_g):
    o = np.argsort(all_fr)
    af, asn, ag = all_fr[o], all_snr[o], all_g[o]
    lo = np.searchsorted(af, fr - XT)
    hi = np.searchsorted(af, fr + XT, side="right")
    out = np.zeros(fr.size, bool)
    for k in range(fr.size):
        s = slice(lo[k], hi[k])
        out[k] = bool(((ag[s] != g) & (asn[s] > snr[k])).any())
    return out


def main():
    assert bk.REFERENCE == "aligned_seed", "run with BAKE_REFERENCE=aligned_seed"
    gt = np.load(bk.BAKE_ROOT / "v2_aligned" / "gt.npz")
    plan = json.loads((bk.BAKE_ROOT / "v2" / "gt_plan.json").read_text())
    leaves = {p["gt"] for p in plan if p["kind"] == "leaves"}
    b0 = bk.Bank.load(bk.BAKE / "bank_start.npz")
    rows = []
    for w in range(12):
        sp = {k: np.load(bk.BAKE / "arms" / v / f"spikes_w{w:02d}.npz") for k, v in ARMS.items()}
        acc = {k: s["reason"] == 0 for k, s in sp.items()}
        ugrp = {}
        for k, s in sp.items():
            bank = bk.Bank.load(bk.BAKE / "arms" / ARMS[k] / f"bank_w{w:02d}.npz")
            ugrp[k] = dict(zip(bank.ids.tolist(), bank.group.tolist()))
        a, b = w * bk.WIN, (w + 1) * bk.WIN
        mw = (gt["frame"] >= a) & (gt["frame"] < b)
        gfr, gun = gt["frame"][mw], gt["unit"][mw]
        for i, g in enumerate(gt["group"]):
            if i in leaves and w >= 6:
                continue
            u, g = 1000 + i, int(g)
            own = np.sort(gfr[gun == i])
            others = np.sort(gfr[(gun != i) & np.isin(gun, np.flatnonzero(gt["group"] == g))])
            for k, s in sp.items():
                ok = acc[k]
                m = ok & (s["unit"] == u)
                fr = s["frame"][m]
                o = np.argsort(fr)
                fr = fr[o]
                feats = {f: s[f][m][o] for f in ("snr", "cos_pl", "a", "trough")}
                ag = np.array([ugrp[k].get(int(x), -1) for x in s["unit"][ok]])
                real = np.sort(s["frame"][ok][(ag == g) & (s["unit"][ok] < 1000)])
                other_arm = [kk for kk in sp if kk != k][0]
                mo = acc[other_arm] & (sp[other_arm]["unit"] == u)
                fo = np.sort(sp[other_arm]["frame"][mo])
                hit = near(fr, own)
                cls = np.full(fr.size, "none", dtype=object)
                cls[near(fr, real)] = "real unit"
                cls[near(fr, others)] = "other ground truth"
                cls[near(fr, own, 11, 45)] = "second fit"
                cls[hit] = "hit"
                shared = near(fr, fo)
                xt = larger_elsewhere(fr, feats["snr"], g, s["frame"][ok], s["snr"][ok], ag)
                for c in np.unique(cls):
                    for sh in (True, False):
                        mm = (cls == c) & (shared == sh)
                        if not mm.any():
                            continue
                        rows.append({"window": w, "gt": i, "tetrode": f"TT{g + 1}", "arm": k, "cls": c,
                                     "in_other_arm": sh, "n": int(mm.sum()), "xtet_larger": float(xt[mm].mean()),
                                     **{f"{f}_med": float(np.nanmedian(v[mm])) for f, v in feats.items()}})
    df = pd.DataFrame(rows)
    df.round(4).to_csv(HERE / "results" / "fp_anatomy.csv", index=False)
    pd.set_option("display.width", 250)

    def agg(d):
        n = d.n.sum()
        return pd.Series({"n": n, "xtet_larger": (d.xtet_larger * d.n).sum() / n,
                          **{c: (d[c] * d.n).sum() / n for c in ("snr_med", "cos_pl_med", "a_med", "trough_med")}})

    print("all ground-truth units, per arm, class and whether the other arm has the spike "
          "(medians are spike-weighted means of per-window medians)")
    print(df.groupby(["arm", "cls", "in_other_arm"]).apply(agg, include_groups=False).round(3).to_string())
    print("\nfalse positives per unit (per tetrode / joint), by class")
    fp = df[df.cls != "hit"].pivot_table(index=["gt", "tetrode"], columns=["arm", "cls"], values="n", aggfunc="sum",
                                         fill_value=0)
    print(fp.to_string())


if __name__ == "__main__":
    main()
