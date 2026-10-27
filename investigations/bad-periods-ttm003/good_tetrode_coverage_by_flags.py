"""TT5 to TT12: coverage of their clean events by the number of other tetrodes flagged in that second.

The baseline's detected events at 10 MAD or more on TT5 to TT12 in the first 14
windows, outside their own tetrode's flagged time, grouped by how many other
tetrodes have a bad period in the same second; per group, the share claimed
within 0.5 ms by any unit of the same tetrode in the baseline, the masked run
and the template-seed control.

Run: .venv/bin/python investigations/bad-periods-ttm003/good_tetrode_coverage_by_flags.py
"""
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _wobble_eval import TOL_MS, _by_tet, _within_tol
from ttm import bad_periods as bp

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
FS = 30000.0; NS = 14 * 1800
runs = {"baseline": "mp_baseline", "masked": "mp_baseline_masked_w14", "seed1": "mp_baseline_seed1_w14"}
z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
ps, pg, am = z["peak_sample"], z["peak_group"], z["amp_mad"]
m = (ps < NS * FS) & (am >= 10) & (pg >= 4) & (pg <= 11)
ps, pg = ps[m], pg[m]
table = bp.read_table(TABLE)
flag = np.zeros((16, NS), bool)
for r in table.itertuples():
    a, b = int(np.floor(r.start_s)), int(np.ceil(r.end_s))
    if a < NS:
        flag[r.group, a:min(b, NS)] = True
sec = (ps // FS).astype(int)
own = flag[pg, sec]
k = flag.sum(0)[sec] - own
ps, pg, k = ps[~own], pg[~own], k[~own]
tol = int(TOL_MS * 1e-3 * FS)
cl = {n: _within_tol(ps, pg, _by_tet(si.load(f"{S}/{d}/assembled_reseed_c12")), tol) for n, d in runs.items()}
bins = [(0, 0), (1, 1), (2, 2), (3, 4), (5, 8)]
for lo, hi in bins:
    i = (k >= lo) & (k <= hi)
    print(f"{lo}-{hi} other tetrodes flagged: {i.sum():8d} events, {len(np.unique(ps[i] // FS)):6d} s | "
          + "  ".join(f"{n} {cl[n][i].mean() * 100:5.1f}" for n in runs))
# single-tetrode flags outside the 2287-2700 s event
out = (ps < 2280 * FS) | (ps > 2710 * FS)
for lo, hi in [(1, 1), (2, 4)]:
    i = (k >= lo) & (k <= hi) & out
    print(f"outside 2280-2710 s, {lo}-{hi} flagged: {i.sum():7d} events | " + "  ".join(f"{n} {cl[n][i].mean() * 100:5.1f}" for n in runs))
