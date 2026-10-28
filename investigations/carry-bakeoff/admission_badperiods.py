"""Stage A, window 1: per-tetrode rules with each tetrode's own bad periods dropped, against the joint runs.

Fits inside a bad period of their own tetrode (``bad_periods.csv``) are dropped; recall and precision are
counted on ground-truth spikes outside the bad periods of their tetrode, for every run, so a masked run
is not charged for the spikes its mask removed.

Run: .venv/bin/python investigations/carry-bakeoff/admission_badperiods.py [--window 1]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402
from admission import load, refractory  # noqa: E402
from score import near  # noqa: E402


def in_own_bad(frames, groups):
    out = np.zeros(frames.size, bool)
    for g, iv in bk.bad_intervals().items():
        if g < 0 or not iv:
            continue
        m = groups == g
        st, en = np.array(iv).T
        out[m] = bk.in_intervals(frames[m], st, en)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=1)
    a = ap.parse_args()
    w = a.window
    gt = np.load(bk.BAKE / "gt.npz")
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    m = (gt["frame"] >= w * bk.WIN) & (gt["frame"] < (w + 1) * bk.WIN)
    gfr, gun = gt["frame"][m], gt["unit"][m]
    gbad = in_own_bad(gfr, gt["group"][gun])
    sizes = {"small": [p["gt"] for p in plan if p["trough_mad"] <= 6],
             "large": [p["gt"] for p in plan if p["trough_mad"] > 6]}
    runs = []
    for name in ("A_joint", "A_joint_masked", "A_pt"):
        z = load(name, w)
        runs.append((name, z, z["reason"] == 0))
    z = load("A_pt_open", w)
    base = z["reason"] != 1
    for k in (4.0, 4.5):
        keep = refractory(z, base & (z["snr"] >= k) & (z["a"] >= 0.5) & (z["cos_pl"] >= 0.6))
        runs.append((f"pt snr>={k} a>=0.5 cos_pl>=0.6 +refr", z, keep))
    rows = []
    for name, z, keep in runs:
        for drop in (False, True):
            kk = keep & ~in_own_bad(z["frame"], z["group"]) if drop else keep
            fr, un = z["frame"][kk], z["unit"][kk]
            row = {"run": name, "own bad periods dropped": drop}
            for cname, mem in sizes.items():
                tp = ntrue = found = nunit = 0
                for i in mem:
                    gi = gfr[(gun == i) & ~gbad]
                    mine = np.sort(fr[un == 1000 + i])
                    mine_ok = mine[~in_own_bad(mine, np.full(mine.size, gt["group"][i]))]
                    found += near(gi, mine).sum()
                    ntrue += gi.size
                    tp += near(mine_ok, np.sort(gi)).sum()
                    nunit += mine_ok.size
                row[f"recall_{cname}"] = found / max(ntrue, 1)
                row[f"precision_{cname}"] = tp / max(nunit, 1)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(bk.BAKE / f"admission_badperiods_w{w:02d}.csv", index=False)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
