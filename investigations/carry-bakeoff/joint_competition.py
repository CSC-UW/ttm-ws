"""Does joint matching lose ground-truth spikes that coincide with events on other tetrodes?

circus-omp (spikeinterface 0.105.2, default ``vicinity=2``) takes, in each pass over a 1 s chunk,
only the best template score over the whole probe within 3 ms either way as a candidate; later
passes can take the next best once the first is fitted, until five passes in a row add no fit inside
the amplitude bounds. Per tetrode runs that loop on each tetrode alone.

On the hybrid, windows 0 to 11: each ground-truth spike is classed by the accepted spikes of the
per-tetrode arm (``W_s45_c50_lnb_b2_dt2``) on other tetrodes within 3 ms either way: none; only
smaller ones (snr below the ground-truth unit's median snr in the window); one or more larger. Per
class and ground-truth size (trough 6 MAD or less, more): the spikes' share, and recall in each arm.
The units that stop after window 5 are left out from window 6 on.

Run: .venv/bin/python investigations/carry-bakeoff/joint_competition.py [ARM ...]
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402

REF = "W_s45_c50_lnb_b2_dt2"
ARMS = ["W_s45_c50_lnb_b2_dt2", "W_joint_masked_s45_c50_lnb_b2", "W_joint_masked_june"]
V = 90


def classify(gi, ref_fr, ref_snr, own_snr):
    """0: no other-tetrode spike within V samples; 1: only smaller ones; 2: a larger one."""
    out = np.zeros(gi.size, np.int8)
    lo, hi = np.searchsorted(ref_fr, gi - V), np.searchsorted(ref_fr, gi + V + 1)
    for k in np.flatnonzero(hi > lo):
        out[k] = 2 if ref_snr[lo[k]:hi[k]].max() > own_snr else 1
    return out


def main():
    arms = sys.argv[1:] or ARMS
    gt = np.load(bk.BAKE / "gt.npz")
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    small = {p["gt"] for p in plan if p["trough_mad"] <= 6}
    leaves = {p["gt"] for p in plan if p["kind"] == "leaves"}
    b0 = bk.Bank.load(bk.BAKE / "bank_start.npz")
    grp = dict(zip(b0.ids.tolist(), b0.group.tolist()))
    rows = []
    for w in range(12):
        a, b = w * bk.WIN, (w + 1) * bk.WIN
        ref = dict(np.load(bk.BAKE / "arms" / REF / f"spikes_w{w:02d}.npz"))
        racc = ref["reason"] == 0
        r_fr, r_un, r_snr = ref["frame"][racc], ref["unit"][racc], ref["snr"][racc]
        r_g = np.array([grp[int(u)] for u in r_un])
        found = {}
        for arm in arms:
            sp = dict(np.load(bk.BAKE / "arms" / arm / f"spikes_w{w:02d}.npz"))
            acc = sp["reason"] == 0
            found[arm] = (sp["frame"][acc], sp["unit"][acc])
        mw = (gt["frame"] >= a) & (gt["frame"] < b)
        for i in range(len(gt["group"])):
            if i in leaves and w >= 6:
                continue
            g = int(gt["group"][i])
            gi = gt["frame"][mw & (gt["unit"] == i)]
            if gi.size == 0:
                continue
            own = r_snr[r_un == 1000 + i]
            own_snr = float(np.median(own)) if own.size else 5.0
            m = r_g != g
            o = np.argsort(r_fr[m])
            cls = classify(gi, r_fr[m][o], r_snr[m][o], own_snr)
            row = {"window": w, "gt": i, "size": "small" if i in small else "large"}
            for arm in arms:
                fr, un = found[arm]
                hit = near(gi, np.sort(fr[un == 1000 + i]))
                for c in range(3):
                    row[f"n{c}"] = int((cls == c).sum())
                    row[f"{arm}|{c}"] = int(hit[cls == c].sum())
            rows.append(row)
    df = pd.DataFrame(rows)
    names = {0: "no other-tetrode spike within 3 ms", 1: "only smaller ones", 2: "a larger one"}
    out = []
    for size, d in df.groupby("size"):
        tot = sum(d[f"n{c}"].sum() for c in range(3))
        for c in range(3):
            n = d[f"n{c}"].sum()
            r = {"size": size, "class": names[c], "share": n / tot}
            for arm in arms:
                r[arm] = d[f"{arm}|{c}"].sum() / max(n, 1)
            out.append(r)
    pd.set_option("display.width", 220)
    print("recall of ground-truth spikes by what the per-tetrode arm accepted on other tetrodes within 3 ms")
    print(pd.DataFrame(out).round(3).to_string(index=False))
    per = []
    for i, d in df.groupby("gt"):
        r = {"gt": i, "tetrode": f"TT{int(gt['group'][i]) + 1}", "size": d["size"].iloc[0]}
        for c in (0, 2):
            n = d[f"n{c}"].sum()
            for arm in arms:
                r[f"{arm.split('_')[1]}{arm.split('_')[-1]}|{c}"] = d[f"{arm}|{c}"].sum() / max(n, 1)
        per.append(r)
    print("per ground-truth unit, recall with no other-tetrode spike (|0) and with a larger one (|2)")
    print(pd.DataFrame(per).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
