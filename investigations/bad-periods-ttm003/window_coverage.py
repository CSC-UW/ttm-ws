"""Clean-time coverage per tetrode and 1800 s window, for any number of 48 h sortings.

The share of the baseline's detected events at 10 MAD or more in each
tetrode's clean time (outside its rows of the bad-period table) that a single
unit of the same tetrode claims within 0.5 ms, per carry-forward window from
W0 to W1 for the tetrodes asked for, then over the 48 h for all of them. With
the masked run and an unmasked control (another template seed), a difference
between masked and baseline can be set against the one between two unmasked
runs.

Run: .venv/bin/python investigations/bad-periods-ttm003/window_coverage.py \
         baseline=mp_baseline masked=mp_baseline_masked control=mp_baseline_seed1 \
         [--tetrodes 5 15] [--windows 0 96] [--sorting assembled_reseed_c12]
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
WIN = int(1800 * FS)
ap = argparse.ArgumentParser()
ap.add_argument("runs", nargs="+", help="name=directory under the session's sortings")
ap.add_argument("--tetrodes", type=int, nargs="*", default=[])
ap.add_argument("--windows", type=int, nargs=2, default=(0, 96))
ap.add_argument("--sorting", default="assembled_reseed_c12")
a = ap.parse_args()
runs = dict(r.split("=") for r in a.runs)
print(" ".join(sys.argv[1:]))

z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
m = z["amp_mad"] >= 10
ps, pg = z["peak_sample"][m], z["peak_group"][m]
table = bp.read_table(TABLE)
clean = np.ones(len(ps), bool)
for g in np.unique(pg):
    k = pg == g
    clean[k] = ~bp.bad_mask(table, int(g), ps[k], FS)
ps, pg = ps[clean], pg[clean]
w = ps // WIN
tol = int(TOL_MS * 1e-3 * FS)
cl = {}
for name, d in runs.items():
    s = si.load(f"{S}/{d}/{a.sorting}")
    if "is_mua" in s.get_property_keys():
        s = s.select_units(np.asarray(s.unit_ids)[~np.asarray(s.get_property("is_mua")).astype(bool)])
    cl[name] = _within_tol(ps, pg, _by_tet(s), tol)


def pct(k):
    return "  ".join(f"{cl[n][k].mean() * 100:6.1f}" if k.any() else "     -" for n in runs)


head = "  ".join(f"{n[:6]:>6s}" for n in runs)
for t in a.tetrodes:
    print(f"\nTT{t}: window  n_events  {head}")
    for i in range(*a.windows):
        k = (pg == t - 1) & (w == i)
        print(f"  {i:3d} {k.sum():8d}  {pct(k)}")
print(f"\nover the 48 h: tetrode  n_events  {head}")
for g in range(16):
    k = pg == g
    print(f"  TT{g + 1:<3d} {k.sum():9d}  {pct(k)}")
print(f"  all   {len(ps):9d}  {pct(np.ones(len(ps), bool))}")
