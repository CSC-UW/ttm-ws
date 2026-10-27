"""What per-tetrode circus-omp adds over joint circus-omp on a stretch: shape, size, refractory, and a noise null.

Same templates and matcher settings as `compare_stretches.py`, unmasked. Spikes
of the two arms are split into found by both (same unit within 3 samples),
joint only and per tetrode only; per tetrode and category: the median cosine
of each spike's snippet to its unit's template on the unit's four wires, the
share below 0.6, and the median amplitude scale (`_mp_common.per_spike_fit`).
Per unit and arm: spikes with an interval under 1.5 ms to the unit's previous
spike. Last, both arms on the same stretch with every channel replaced by
white noise at its clean noise level (`mask_bad_periods` with every tetrode
flagged): spikes there are fits to noise.

Run: .venv/bin/python investigations/per-tetrode-matching/extra_spikes.py T0 T1
"""
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _mp_common import build_templates_object, per_spike_fit
from ttm import bad_periods as bp
from ttm.matching import find_spikes_per_group

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
FS = 30000.0
WIN = int(1800 * FS)
AMP = {"amplitudes": [0.8, float("inf")]}
JOB = {"n_jobs": 16, "chunk_duration": "1s", "progress_bar": False}
RP = int(1.5e-3 * FS)
T0, T1 = float(sys.argv[1]), float(sys.argv[2])
f0, f1 = int(T0 * FS), int(T1 * FS)

rec = si.load(f"{S}/mp_baseline/binary")
base = si.load(f"{S}/mp_baseline/assembled_reseed_c12")
w = max(f0 // WIN - 1, 0)
lo, hi = w * WIN, (w + 1) * WIN
ids = [u for u in base.unit_ids if ((st := base.get_unit_spike_train(u)) >= lo).sum() - (st >= hi).sum() >= 100]
sw = base.select_units(ids).frame_slice(lo, hi)
rw = rec.frame_slice(lo, hi)
rw.reset_times()
sw.register_recording(rw)
tmpl, _ = build_templates_object(sw, rw, with_snr=False, n_jobs=32)
groups = np.asarray(rec.get_property("group")).astype(int)
ugroup = np.array([groups[np.flatnonzero(m)[0]] for m in tmpl.sparsity.mask])
print(f"{T0:.0f} to {T1:.0f} s, templates from window {w} ({len(ids)} units)")


def match(r):
    seg = r.frame_slice(f0, f1)
    seg.reset_times()
    j = find_spikes_from_templates(seg, tmpl, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    p = find_spikes_per_group(seg, tmpl, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    return seg, j, p


def in_other(x, y):
    hit = np.zeros(len(x), bool)
    for u in np.unique(x["cluster_index"]):
        mx = x["cluster_index"] == u
        sy = np.sort(y["sample_index"][y["cluster_index"] == u])
        if len(sy) == 0:
            continue
        sx = x["sample_index"][mx]
        i = np.clip(np.searchsorted(sy, sx), 1, len(sy) - 1) if len(sy) > 1 else np.zeros(len(sx), int)
        d = np.minimum(np.abs(sy[i] - sx), np.abs(sy[np.maximum(i - 1, 0)] - sx))
        hit[mx] = d <= 3
    return hit


seg, j, p = match(rec)
aj, rj = per_spike_fit(j, tmpl, seg)
app, rp = per_spike_fit(p, tmpl, seg)
j_both, p_both = in_other(j, p), in_other(p, j)
rows = []
for g in range(16):
    for name, sel, a, r in [("both", (ugroup[p["cluster_index"]] == g) & p_both, app, rp),
                            ("joint only", (ugroup[j["cluster_index"]] == g) & ~j_both, aj, rj),
                            ("per tetrode only", (ugroup[p["cluster_index"]] == g) & ~p_both, app, rp)]:
        if sel.sum():
            rows.append({"tetrode": f"TT{g + 1}", "category": name, "n": int(sel.sum()),
                         "cos_median": round(float(np.nanmedian(r[sel])), 3),
                         "cos_below_0.6_pct": round(float(np.mean(r[sel] < 0.6) * 100), 1),
                         "amp_median": round(float(np.nanmedian(a[sel])), 3)})
pd.set_option("display.width", 200)
print("\nspikes by category, per tetrode")
print(pd.DataFrame(rows).to_string(index=False))


def rp_violations(sp):
    out = {}
    for u in np.unique(sp["cluster_index"]):
        st = np.sort(sp["sample_index"][sp["cluster_index"] == u])
        out[int(u)] = (len(st), int((np.diff(st) < RP).sum()))
    return out


vj, vp = rp_violations(j), rp_violations(p)
print("\nper unit: spikes and intervals under 1.5 ms, joint against per tetrode")
print("unit  tetrode   joint_n  joint_rp   pertet_n  pertet_rp")
for u in sorted(set(vj) | set(vp)):
    nj, cj = vj.get(u, (0, 0))
    n_p, cp = vp.get(u, (0, 0))
    print(f"{tmpl.unit_ids[u]:4}  TT{ugroup[u] + 1:<5d} {nj:8d} {cj:9d}  {n_p:9d} {cp:10d}")

flag_all = pd.DataFrame({"group": range(16), "tetrode": [f"TT{g + 1}" for g in range(16)],
                         "start_s": T0 - 5.0, "end_s": T1 + 5.0, "kind": "null", "channels": "",
                         "max_loud": 0.0, "max_corr_excess": 0.0, "flagged_s": T1 - T0})
null = bp.mask_bad_periods(rec, flag_all, noise_levels=bp.clean_noise_levels(
    rec, [(int(a * FS), int(b * FS), np.flatnonzero(groups == g).tolist())
          for g, a, b in bp.read_table(TABLE)[["group", "start_s", "end_s"]].itertuples(index=False)]))
_, nj, np_ = match(null)
print(f"\nwhite noise at the clean noise level, {T1 - T0:.0f} s: joint {len(nj)} spikes, per tetrode {len(np_)}")
for g in range(16):
    a_, b_ = (ugroup[nj["cluster_index"]] == g).sum(), (ugroup[np_["cluster_index"]] == g).sum()
    if a_ or b_:
        print(f"  TT{g + 1:<4d} joint {a_:6d}  per tetrode {b_:6d}")
