"""The two matchers on the hybrid with the global median and with the aligned median reference.

Per matcher (``W_s45_c50_lnb_b2`` per tetrode, ``W_joint_masked_s45_c50_lnb_b2`` joint, same rules
and rebuild) and reference: ``v2`` (global median, mp_baseline's seed bank), ``v2_aligned``
(aligned median, the same seed bank) and ``v2_aligned_seed`` (aligned median, the seed bank
re-sorted on it, ``seed_aligned.py``). The ground truth is the same in all three. Means over the 12
windows of ``compare_carries.pooled`` (small: trough 6 MAD or less): recall, precision and F1, the
half-size in-pulse recall; per ground-truth unit, precision and recall over the windows; and the
real units' accepted spikes per tetrode summed over the windows, for the tetrodes with few units.

Writes ``results/compare_reference.csv`` (the pooled rows) and
``results/compare_reference_units.csv`` (per ground-truth unit).

Run: .venv/bin/python investigations/carry-bakeoff/compare_reference.py (after score.py on each)
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from compare_carries import pooled  # noqa: E402

HERE = pathlib.Path(__file__).parent
REFS = {"global": "v2", "aligned": "v2_aligned", "aligned, re-seeded": "v2_aligned_seed"}
ARMS = {"per tetrode": "W_s45_c50_lnb_b2", "joint": "W_joint_masked_s45_c50_lnb_b2"}
FEW = ("TT1", "TT2", "TT3", "TT13", "TT14", "TT15", "TT16", "TT9")


def f1(r, p):
    return 2 * r * p / (r + p) if r + p else 0.0


def main():
    plan = json.loads((bk.BAKE_ROOT / "v2" / "gt_plan.json").read_text())
    small = [p["gt"] for p in plan if p["trough_mad"] <= 6]
    large = [p["gt"] for p in plan if p["trough_mad"] > 6]
    leaves = {p["gt"] for p in plan if p["kind"] == "leaves"}
    rows, units, real = [], [], {}
    for ref, sub in REFS.items():
        for arm, name in ARMS.items():
            d = bk.BAKE_ROOT / sub / "arms" / name
            g = pd.read_csv(d / "score_gt.csv")
            g = g[~(g["gt"].isin(leaves) & (g.window >= 6))]
            per_w = [(pooled(gw, small), pooled(gw, large)) for _, gw in g.groupby("window")]
            rs = np.mean([s.recall for s, _ in per_w]); ps = np.mean([s.precision for s, _ in per_w])
            rl = np.mean([l_.recall for _, l_ in per_w]); pl = np.mean([l_.precision for _, l_ in per_w])
            half = g.recall_laser_half.dropna()
            rows.append({"reference": ref, "matcher": arm, "small_recall": rs, "small_precision": ps,
                         "small_f1": f1(rs, ps), "large_recall": rl, "large_precision": pl,
                         "large_f1": f1(rl, pl), "laser_half": half.mean()})
            for u, gu in g.groupby("gt"):
                found = (gu.recall * gu.n_true).sum() / gu.n_true.sum()
                hits = (gu.precision.fillna(0) * gu.n_unit).sum() / max(gu.n_unit.sum(), 1)
                units.append({"reference": ref, "matcher": arm, "gt": u, "tetrode": gu.tetrode.iloc[0],
                              "recall": found, "precision": hits})
            r = pd.read_csv(d / "score_real.csv")
            real[(ref, arm)] = r.groupby("tetrode").n.sum()
    df = pd.DataFrame(rows).round(3)
    du = pd.DataFrame(units).round(3)
    df.to_csv(HERE / "results" / "compare_reference.csv", index=False)
    du.to_csv(HERE / "results" / "compare_reference_units.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.to_string(index=False))
    print("\nper ground-truth unit, precision (recall)")
    piv = du.assign(v=du.precision.map("{:.2f}".format) + " (" + du.recall.map("{:.2f}".format) + ")")
    print(piv.pivot_table(index=["gt", "tetrode"], columns=["matcher", "reference"], values="v", aggfunc="first")
          .to_string())
    print("\nreal units' accepted spikes over the 12 windows, tetrodes with few units")
    rt = pd.DataFrame({f"{a}, {r}": s for (r, a), s in real.items()}).reindex(list(FEW)).fillna(0)
    print(rt.astype(int).to_string())


if __name__ == "__main__":
    main()
