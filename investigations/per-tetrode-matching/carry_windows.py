"""Window by window, what the carry-forward gives the units of a few tetrodes, joint against per tetrode.

Same runs and scores as `carry_compare.py`. Per tetrode asked for, run and
1800 s window: the clean-time spikes of its single units, their median trough
amplitude on the deepest wire (MAD), and their median cosine to the unit's
window-local template.

Run: .venv/bin/python investigations/per-tetrode-matching/carry_windows.py \
         joint=mp_baseline pertet=mp_baseline_pertet_w14 \
         joint_masked=mp_baseline_masked_w14 pertet_masked=mp_baseline_masked_pertet_w14 \
         --tetrodes 1 6 9 13 14 16
"""
import argparse

import numpy as np
import spikeinterface as si
from ttm import bad_periods as bp

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
SORTING = "assembled_reseed_c12"
FS = 30000.0
WIN = int(1800 * FS)
N_WIN = 14
ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+", help="name=directory under the session's sortings")
ap.add_argument("--tetrodes", type=int, nargs="+", default=[1, 9, 13, 14, 16])
a = ap.parse_args()
runs = dict(r.split("=") for r in a.runs)
table = bp.read_table(TABLE)

out = {}
for name, d in runs.items():
    s = si.load(f"{S}/{d}/{SORTING}")
    sc = np.load(f"{S}/{d}/purity_scores_{SORTING}_w14.npz")
    mua = (np.asarray(s.get_property("is_mua")).astype(bool) if "is_mua" in s.get_property_keys()
           else np.zeros(len(s.unit_ids), bool))
    grp = np.asarray(s.get_property("group")).astype(int)
    for t in a.tetrodes:
        f, am, r = [], [], []
        for u in np.asarray(s.unit_ids)[(grp == t - 1) & ~mua]:
            fu = sc[f"frames_{u}"]
            k = (fu < N_WIN * WIN) & ~bp.bad_mask(table, t - 1, fu, FS)
            f.append(fu[k])
            am.append(sc[f"amps_{u}"][k])
            r.append(sc[f"rs_{u}"][k])
        f, am, r = (np.concatenate(x) if x else np.zeros(0) for x in (f, am, r))
        w = (f // WIN).astype(int)
        out[(name, t)] = [(int((w == i).sum()), np.median(np.abs(am[w == i])) if (w == i).any() else np.nan,
                           np.nanmedian(r[w == i]) if (w == i).any() else np.nan) for i in range(N_WIN)]

for t in a.tetrodes:
    print(f"\nTT{t}: clean spikes / median trough (MAD) / median cosine per window")
    print("  window " + "".join(f"{n:>26s}" for n in runs))
    for i in range(N_WIN):
        print(f"  {i:6d} " + "".join(f"{out[(n, t)][i][0]:>10d} {out[(n, t)][i][1]:6.1f} {out[(n, t)][i][2]:7.3f}"
                                     for n in runs))
