"""Are the spikes only per-tetrode circus-omp finds small fits, or events shared with other tetrodes?

Same templates and settings as `compare_stretches.py`, unmasked. Per spike,
the matched-filter signal-to-noise: the projection of its snippet on its unit's
unit-norm template over the noise level of the unit's wires. Per tetrode class:
its percentiles for joint spikes, per-tetrode spikes joint also has (same unit
within 3 samples) and per-tetrode spikes joint lacks; per-tetrode spikes and
their share of intervals under 1.5 ms above floors of that signal-to-noise;
and the share of spikes with a spike of the same arm on another tetrode
within 0.2 ms.

Run: .venv/bin/python investigations/per-tetrode-matching/snr_and_coincidence.py T0 T1
"""
import sys

import numpy as np
import spikeinterface as si
from spikeinterface.core import get_noise_levels
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _mp_common import build_templates_object, per_spike_fit
from ttm.matching import find_spikes_per_group

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
FS = 30000.0
WIN = int(1800 * FS)
RP = int(1.5e-3 * FS)
JOB = {"n_jobs": 16, "chunk_duration": "1s", "progress_bar": False}
AMP = {"amplitudes": [0.8, float("inf")]}
CLASSES = {"TT4": [3], "TT5-TT12": list(range(4, 12)), "few-unit": [0, 1, 2, 12, 13, 14, 15], "TT9": [8]}
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
ug = np.array([groups[np.flatnonzero(m)[0]] for m in tmpl.sparsity.mask])
dense = tmpl.get_dense_templates()
tnorm = np.array([np.linalg.norm(dense[i][:, tmpl.sparsity.mask[i]]) for i in range(len(ug))])
seg = rec.frame_slice(f0, f1)
seg.reset_times()
nl = get_noise_levels(seg, return_in_uV=False)
sig = np.array([np.mean(nl[tmpl.sparsity.mask[i]]) for i in range(len(ug))])
print(f"{T0:.0f} to {T1:.0f} s, templates from window {w}; template ||t|| over noise per unit: "
      f"{dict(zip(tmpl.unit_ids.tolist(), np.round(tnorm / sig, 1).tolist()))}")

j = find_spikes_from_templates(seg, tmpl, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
p = find_spikes_per_group(seg, tmpl, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)


def snr(sp):
    a, _ = per_spike_fit(sp, tmpl, seg)
    ci = sp["cluster_index"]
    return a * tnorm[ci] / sig[ci]


def in_other(x, y):
    hit = np.zeros(len(x), bool)
    for u in np.unique(x["cluster_index"]):
        mx = x["cluster_index"] == u
        sy, sx = np.sort(y["sample_index"][y["cluster_index"] == u]), x["sample_index"][mx]
        if len(sy) > 1:
            i = np.clip(np.searchsorted(sy, sx), 1, len(sy) - 1)
            hit[mx] = np.minimum(np.abs(sy[i] - sx), np.abs(sy[i - 1] - sx)) <= 3
    return hit


def short_isi(sp):
    out = np.zeros(len(sp), bool)
    for u in np.unique(sp["cluster_index"]):
        idx = np.flatnonzero(sp["cluster_index"] == u)
        o = idx[np.argsort(sp["sample_index"][idx])]
        out[o[1:]] = np.diff(sp["sample_index"][o]) < RP
    return out


def coincident(sp, tol=6):
    out = np.zeros(len(sp), bool)
    g_sp = ug[sp["cluster_index"]]
    for g in range(16):
        m = g_sp == g
        st = np.sort(sp["sample_index"][g_sp != g])
        if m.any() and len(st) > 1:
            x = sp["sample_index"][m]
            i = np.clip(np.searchsorted(st, x), 1, len(st) - 1)
            out[m] = np.minimum(np.abs(st[i] - x), np.abs(st[i - 1] - x)) <= tol
    return out


sj, s_p = snr(j), snr(p)
pb, vj, vp, cj, cp = in_other(p, j), short_isi(j), short_isi(p), coincident(j), coincident(p)


def q(v):
    return np.round(np.nanpercentile(v, [10, 50, 90]), 1).tolist()


print("\nsignal-to-noise, 10th/50th/90th percentile: joint | per tetrode, also joint | per tetrode only")
for n, gs in CLASSES.items():
    mj, mp = np.isin(ug[j["cluster_index"]], gs), np.isin(ug[p["cluster_index"]], gs)
    print(f"  {n:9s} {q(sj[mj])} | {q(s_p[mp & pb])} | {q(s_p[mp & ~pb])}")
print("\nper-tetrode spikes at or above a signal-to-noise floor: spikes, intervals under 1.5 ms (joint first)")
for n, gs in CLASSES.items():
    mj, mp = np.isin(ug[j["cluster_index"]], gs), np.isin(ug[p["cluster_index"]], gs)
    row = f"  {n:9s} joint {mj.sum():6d} {vj[mj].mean() * 100:5.2f}% |"
    for th in [0, 6, 8, 10, 12]:
        k = mp & (s_p >= th)
        row += f" >={th}: {k.sum():6d} {vp[k].mean() * 100:5.2f}% |"
    print(row)
print("\nshare with a spike on another tetrode within 0.2 ms: joint | per tetrode, also joint | per tetrode only")
for n, gs in CLASSES.items():
    mj, mp = np.isin(ug[j["cluster_index"]], gs), np.isin(ug[p["cluster_index"]], gs)
    print(f"  {n:9s} {cj[mj].mean() * 100:5.1f}% | {cp[mp & pb].mean() * 100:5.1f}% | {cp[mp & ~pb].mean() * 100:5.1f}%")
print("\nper tetrode only: intervals under 1.5 ms among those coincident with another tetrode and the rest")
for n, gs in CLASSES.items():
    mp = np.isin(ug[p["cluster_index"]], gs) & ~pb
    print(f"  {n:9s} coincident {vp[mp & cp].mean() * 100:5.1f}% (n={int((mp & cp).sum())}), "
          f"not {vp[mp & ~cp].mean() * 100:5.1f}% (n={int((mp & ~cp).sum())})")
