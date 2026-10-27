"""Does an artifact on some tetrodes stop circus-omp on the others?

Templates like the baseline's in window 1: every unit of
`mp_baseline/assembled_reseed_c12` with 100 spikes or more in window 0, built
from its window-0 spikes the way the carry-forward re-estimates them. Then
circus-omp as the carry runs it (amplitudes 0.8 and up, 1 s chunks) on
T0 to T1 s of the baseline binary in three ways: unmasked with all units,
masked (`mask_bad_periods` with the session's table) with all units, and
unmasked with only the units of TT5 to TT12, which have no flagged time there.
Last, where in each 1 s chunk the TT5 to TT12 spikes fall, unmasked and masked.

Run: .venv/bin/python investigations/bad-periods-ttm003/omp_coupling.py 2400 2460
"""
import sys
import time

import numpy as np
import spikeinterface as si

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _mp_common import build_templates_object, run_matching
from ttm import bad_periods as bp

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
FS = 30000.0
AMP = {"amplitudes": [0.8, float("inf")]}
GOOD = range(4, 12)
T0, T1 = float(sys.argv[1]), float(sys.argv[2])

rec = si.load(f"{S}/mp_baseline/binary")
masked = bp.mask_bad_periods(rec, bp.read_table(TABLE))
base = si.load(f"{S}/mp_baseline/assembled_reseed_c12")
w0 = int(1800 * FS)
ids = [u for u in base.unit_ids if (base.get_unit_spike_train(u) < w0).sum() >= 100]
s0 = base.select_units(ids).frame_slice(0, w0)
win0 = rec.frame_slice(0, w0)
win0.reset_times()
s0.register_recording(win0)
tmpl, _ = build_templates_object(s0, win0, with_snr=False, n_jobs=32)
ch_group = dict(zip(rec.channel_ids, rec.get_property("group")))
mask = tmpl.sparsity.mask
ugroup = {u: ch_group[rec.channel_ids[np.flatnonzero(mask[i])[0]]] for i, u in enumerate(tmpl.unit_ids)}
good_ids = [u for u in tmpl.unit_ids if ugroup[u] in GOOD]
print(f"{len(ids)} units with >= 100 window-0 spikes, {len(good_ids)} on TT5-TT12; {T0:.0f} to {T1:.0f} s")


def match(r, t):
    seg = r.frame_slice(int(T0 * FS), int(T1 * FS))
    seg.reset_times()
    s, _ = run_matching(seg, t, method_kwargs=AMP, n_jobs=16)
    return s


for name, r, t in [("unmasked, all units", rec, tmpl), ("masked, all units", masked, tmpl),
                   ("unmasked, TT5-TT12 units only", rec, tmpl.select_units(good_ids))]:
    t0 = time.time()
    s = match(r, t)
    n_good = sum(len(s.get_unit_spike_train(u)) for u in s.unit_ids if ugroup[u] in GOOD)
    n_other = sum(len(s.get_unit_spike_train(u)) for u in s.unit_ids if ugroup[u] not in GOOD)
    print(f"{name:32s} TT5-TT12 {n_good:8d}   other tetrodes {n_other:8d}   ({time.time() - t0:.0f} s)", flush=True)

for name, r in [("unmasked, all units", rec), ("masked, all units", masked)]:
    s = match(r, tmpl)
    f = np.concatenate([s.get_unit_spike_train(u) for u in s.unit_ids if ugroup[u] in GOOD])
    print(f"{name:32s} TT5-TT12 spikes by tenth of the 1 s chunk: "
          f"{np.bincount((f % int(FS)) // int(FS / 10), minlength=10).tolist()}", flush=True)
