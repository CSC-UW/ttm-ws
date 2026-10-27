"""Window 1 (1800 to 3600 s) on TT5 to TT12, which have no flagged time then.

Per 60 s: the baseline's detected events at 10 MAD or more on TT5 to TT12, the
share claimed within 0.5 ms by any unit of the same tetrode in the baseline,
the masked run and the template-seed control, and the spikes of each sorting on
those tetrodes. The 2287 to 2696 s event on the even-channel tetrodes falls in
this window.

Run: .venv/bin/python investigations/bad-periods-ttm003/window1_coverage.py
"""
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _wobble_eval import TOL_MS, _by_tet, _within_tol

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
FS = 30000.0
runs = {"baseline": "mp_baseline", "masked": "mp_baseline_masked_w14", "seed1": "mp_baseline_seed1_w14"}
z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
ps, pg, am = z["peak_sample"], z["peak_group"], z["amp_mad"]
good = (pg >= 4) & (pg <= 11)
m = (ps >= 1800 * FS) & (ps < 3600 * FS) & (am >= 10) & good
ps, pg = ps[m], pg[m]
tol = int(TOL_MS * 1e-3 * FS)
edges = np.arange(1800, 3601, 60)
b = np.digitize(ps / FS, edges) - 1
cl, sp = {}, {}
for name, d in runs.items():
    s = si.load(f"{S}/{d}/assembled_reseed_c12")
    bt = _by_tet(s)
    cl[name] = _within_tol(ps, pg, bt, tol)
    allsp = np.concatenate([bt[g] for g in range(4, 12) if g in bt])
    sp[name] = np.histogram(allsp / FS, edges)[0]
print("t0     n_ev  cov_base cov_mask cov_seed1 | spikes base mask seed1 (TT5-TT12)")
for i in range(len(edges) - 1):
    k = b == i
    print(f"{edges[i]:5.0f} {k.sum():6d}  " + "  ".join(f"{cl[n][k].mean()*100:6.1f}" for n in runs)
          + "  | " + " ".join(f"{sp[n][i]:7d}" for n in runs))
