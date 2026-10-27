"""How stable is TT15's three-unit seed sort? TT15's 1800 s seed window re-sorted 14 ways.

MountainSort5 scheme 2 as `ttm.tracking.sort_chunk` runs it in `mp_baseline.py`,
on TT15's four channels only: unmasked and masked (`mask_bad_periods` with the
session's table) with whitening seeds 42 (the pipeline's) and 1 to 4, and
unmasked with the window moved by 5, 10, 30 and 60 s. Each cluster is matched
to the baseline seed sort of TT15 by spike times within 0.3 ms over 60 to
1800 s.

Run: .venv/bin/python investigations/bad-periods-ttm003/tt15_seed_stability.py OUT_DIR
"""
import pathlib
import sys

import numpy as np
import spikeinterface as si
from joblib import Parallel, delayed
from ttm import bad_periods as bp
from ttm.tracking import sort_chunk, to_int_numpy_sorting

S = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix")
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
OUT = pathlib.Path(sys.argv[1])
FS = 30000.0
GROUP = 14
TOL = 9  # 0.3 ms


def run(name, masked, wseed, start_s):
    rec = si.load(S / "mp_baseline" / "binary")
    if masked:
        rec = bp.mask_bad_periods(rec, bp.read_table(TABLE))
    ch = rec.channel_ids[rec.get_property("group") == GROUP]
    r = rec.select_channels(ch).frame_slice(int(start_s * FS), int((start_s + 1800) * FS))
    r.reset_times()
    s = to_int_numpy_sorting(sort_chunk(r, OUT / name, whitening_seed=wseed, sort_n_jobs=1))
    return name, start_s, {int(u): s.get_unit_spike_train(u) + int(start_s * FS) for u in s.unit_ids}


def agree(a, b):
    """Fraction of a's spikes with a b spike within TOL samples."""
    i = np.searchsorted(b, a)
    d = np.minimum(np.abs(a - b[np.clip(i, 0, len(b) - 1)]), np.abs(a - b[np.clip(i - 1, 0, len(b) - 1)]))
    return float((d <= TOL).mean())


variants = [("unmasked_s42", False, 42, 0.0), ("masked_s42", True, 42, 0.0)]
variants += [(f"unmasked_s{w}", False, w, 0.0) for w in (1, 2, 3, 4)]
variants += [(f"masked_s{w}", True, w, 0.0) for w in (1, 2, 3, 4)]
variants += [(f"unmasked_t{t}", False, 42, float(t)) for t in (5, 10, 30, 60)]
res = Parallel(n_jobs=len(variants))(delayed(run)(*v) for v in variants)

ref = si.read_npz_sorting(S / "mp_baseline" / "seed_sort" / "by_group" / str(GROUP) / "sorter_output" / "firings.npz")
base = {int(u): ref.get_unit_spike_train(u) for u in ref.unit_ids}
lo, hi = 60 * FS, 1800 * FS  # common span of all variants
base = {u: t[(t >= lo) & (t < hi)] for u, t in base.items()}
print("reference (baseline seed sort, TT15):", {u: len(t) for u, t in base.items()})
for name, start_s, units in res:
    units = {u: t[(t >= lo) & (t < hi)] for u, t in units.items()}
    print(f"\n{name}: {len(units)} units")
    for u, t in sorted(units.items()):
        best = max(base, key=lambda b: agree(base[b], t) if len(t) else 0)
        print(f"  unit {u}: {len(t):7d} spikes; covers {agree(base[best], t):.3f} of ref {best}, "
              f"{agree(t, base[best]):.3f} of its spikes in ref {best}")
