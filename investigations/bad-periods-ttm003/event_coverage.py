"""TT5 to TT12 around any stretch: event coverage and spikes per 60 s, baseline against a masked run.

Like `window1_coverage.py`, for the full 48 h sortings: per 60 s from T0 to
T1, the baseline's detected events at 10 MAD or more on TT5 to TT12, the
share claimed within 0.5 ms by a single unit of the same tetrode in each
sorting, the spikes of each on those tetrodes, and how many tetrodes have a
bad period in that minute.

Run: .venv/bin/python investigations/bad-periods-ttm003/event_coverage.py 96300 98100 \
         [--sorting assembled_reseed_c12] [--masked mp_baseline_masked]
"""
import argparse
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _wobble_eval import TOL_MS, _by_tet, _within_tol
from ttm import bad_periods as bp

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
FS = 30000.0
ap = argparse.ArgumentParser()
ap.add_argument("t0", type=float)
ap.add_argument("t1", type=float)
ap.add_argument("--sorting", default="assembled_reseed_c12")
ap.add_argument("--masked", default="mp_baseline_masked")
a = ap.parse_args()
runs = {"baseline": "mp_baseline", "masked": a.masked}

z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
ps, pg, am = z["peak_sample"], z["peak_group"], z["amp_mad"]
good = (pg >= 4) & (pg <= 11)
m = (ps >= a.t0 * FS) & (ps < a.t1 * FS) & (am >= 10) & good
ps, pg = ps[m], pg[m]
table = bp.read_table(TABLE)
edges = np.arange(a.t0, a.t1 + 1, 60)
b = np.digitize(ps / FS, edges) - 1
tol = int(TOL_MS * 1e-3 * FS)
cl, sp = {}, {}
for name, d in runs.items():
    s = si.load(f"{S}/{d}/{a.sorting}")
    if "is_mua" in s.get_property_keys():
        s = s.select_units(np.asarray(s.unit_ids)[~np.asarray(s.get_property("is_mua")).astype(bool)])
    bt = _by_tet(s)
    cl[name] = _within_tol(ps, pg, bt, tol)
    allsp = np.concatenate([bt[g] for g in range(4, 12) if g in bt])
    sp[name] = np.histogram(allsp / FS, edges)[0]
print("t0      n_ev  flagged_tt  cov_base cov_mask | spikes base mask (TT5-TT12)")
for i in range(len(edges) - 1):
    k = b == i
    t = table[(table.start_s < edges[i + 1]) & (table.end_s > edges[i])]
    print(f"{edges[i]:7.0f} {k.sum():6d}  {t.group.nunique():3d}        "
          + "  ".join(f"{cl[n][k].mean() * 100:6.1f}" if k.any() else "   nan" for n in runs)
          + "  | " + " ".join(f"{sp[n][i]:7d}" for n in runs))
