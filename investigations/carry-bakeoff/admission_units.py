"""Stage A per ground-truth unit: recall and precision of the joint runs against per-tetrode rules, one window.

Run: .venv/bin/python investigations/carry-bakeoff/admission_units.py [--window 0]
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

RULES = {"s4 a0.5 cos0.6": (4.0, 0.5, 0.6), "s4.5 a0.5 cos0.6": (4.5, 0.5, 0.6)}


def per_unit(z, keep, gt, w):
    a, b = w * bk.WIN, (w + 1) * bk.WIN
    m = (gt["frame"] >= a) & (gt["frame"] < b)
    out = {}
    fr, un = z["frame"][keep], z["unit"][keep]
    for i in range(len(gt["group"])):
        gi = gt["frame"][m & (gt["unit"] == i)]
        mine = np.sort(fr[un == 1000 + i])
        out[i] = (near(gi, mine).mean(), near(mine, np.sort(gi)).mean() if mine.size else np.nan)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=0)
    a = ap.parse_args()
    gt = np.load(bk.BAKE / "gt.npz")
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    res = {}
    for name in ("A_joint", "A_joint_masked"):
        z = load(name, a.window)
        res[name] = per_unit(z, z["reason"] == 0, gt, a.window)
    z = load("A_pt_open", a.window)
    for rn, (k, r, c) in RULES.items():
        keep = (z["reason"] != 1) & (z["snr"] >= k) & (z["a"] >= r) & (z["cos_pl"] >= c)
        res[rn] = per_unit(z, refractory(z, keep), gt, a.window)
    rows = []
    for p in plan:
        row = {"gt": p["gt"], "tetrode": p["tetrode"], "mad": p["trough_mad"], "kind": p["kind"]}
        for name, d in res.items():
            row[f"{name} R"], row[f"{name} P"] = d[p["gt"]]
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(bk.BAKE / f"admission_units_w{a.window:02d}.csv", index=False)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.2f}"))


if __name__ == "__main__":
    main()
