"""Whether ground-truth units' false positives sit next to large units' accepted spikes more often
than their hits do (a sign that large-unit subtraction residuals produce them).

Run: .venv/bin/python investigations/template-superres/fp_near_large.py w03_w08
"""
import sys

import numpy as np
import pandas as pd

OUT = ("/tmp/claude-295375638/-Users-gfindlay-ad-wisc-edu-projects-ttm-ws/"
       "e980bd61-79fb-4043-a37f-61959b628a69/scratchpad/superres/")
ARM = ("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix/carry_bakeoff/"
       "v2_aligned_seed/arms/W_s45_c50_lnb_b2/")

d = pd.read_parquet(OUT + f"fits_{sys.argv[1]}.parquet")
tr = d.groupby(["window", "unit"]).trough_mad.first()
rows = []
for w in sorted(d.window.unique()):
    sp = np.load(ARM + f"spikes_w{w:02d}.npz")
    acc = sp["reason"] == 0
    for g in sorted(d[d.window == w].group.unique()):
        x = d[(d.window == w) & (d.group == g) & (d.reason == 0) & (d.unit >= 1000)]
        units = d[(d.window == w) & (d.group == g)].unit.unique()
        for u_s in x.unit.unique():
            big = [u for u in units if u != u_s and tr[(w, u)] >= 9]
            if not big:
                continue
            f_big = np.sort(sp["frame"][acc & np.isin(sp["unit"], big)])
            y = x[x.unit == u_s]
            j = np.clip(np.searchsorted(f_big, y.frame.to_numpy()), 1, max(f_big.size - 1, 1))
            dist = np.minimum(np.abs(f_big[j - 1] - y.frame.to_numpy()), np.abs(f_big[j] - y.frame.to_numpy()))
            for lab, m in (("hit", y.hit.astype(bool).to_numpy()), ("fp", ~y.hit.astype(bool).to_numpy())):
                rows.append({"window": w, "group": g, "unit": u_s, "small": tr[(w, u_s)] <= 6, "label": lab,
                             "n": int(m.sum()), "near_60": float(np.mean(dist[m] <= 60)) if m.any() else np.nan})
t = pd.DataFrame(rows)
print(t.pivot_table(index=["window", "group", "unit", "small"], columns="label", values=["n", "near_60"]).round(3).to_string())
