"""Per-tetrode circus-omp settings against joint circus-omp on a clean stretch: recall and precision proxies.

Same templates as `compare_stretches.py`. Arms: joint with June's settings
(amplitudes 0.8 and up, max_failures 5), and per tetrode with each amplitude
floor and max_failures given. Per arm and tetrode class (TT4; TT5 to TT12;
the few-unit tetrodes TT1 to TT3 and TT13 to TT16; TT9 alone, its units being
small): spikes, the share of the baseline's detected events at 10 MAD or more
in clean time claimed by a unit of the same tetrode within 0.5 ms, the share
of spikes with an interval under 1.5 ms to the unit's previous spike, and the
median cosine of each spike's snippet to its unit's template.

Run: .venv/bin/python investigations/per-tetrode-matching/settings_sweep.py T0 T1 \
         [--floors 0.8 0.9 1.0] [--failures 5 2 1]
"""
import argparse
import sys
import time

import numpy as np
import spikeinterface as si
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _mp_common import build_templates_object, per_spike_fit
from _wobble_eval import TOL_MS
from ttm import bad_periods as bp
from ttm.matching import find_spikes_per_group

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
FS = 30000.0
WIN = int(1800 * FS)
JOB = {"n_jobs": 16, "chunk_duration": "1s", "progress_bar": False}
RP = int(1.5e-3 * FS)
CLASSES = {"TT4": [3], "TT5-TT12": list(range(4, 12)), "few-unit": [0, 1, 2, 12, 13, 14, 15], "TT9": [8]}
ap = argparse.ArgumentParser()
ap.add_argument("t0", type=float)
ap.add_argument("t1", type=float)
ap.add_argument("--floors", type=float, nargs="+", default=[0.8, 0.9, 1.0])
ap.add_argument("--failures", type=int, nargs="+", default=[5, 2, 1])
a = ap.parse_args()
f0, f1 = int(a.t0 * FS), int(a.t1 * FS)

rec = si.load(f"{S}/mp_baseline/binary")
table = bp.read_table(TABLE)
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
seg = rec.frame_slice(f0, f1)
seg.reset_times()

z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
k = (z["peak_sample"] >= f0) & (z["peak_sample"] < f1) & (z["amp_mad"] >= 10)
ps, pg = z["peak_sample"][k] - f0, z["peak_group"][k].astype(int)
clean = np.ones(len(ps), bool)
for g in range(16):
    m = pg == g
    clean[m] = ~bp.bad_mask(table, g, ps[m] + f0, FS)
ps, pg = ps[clean], pg[clean]
tol = int(TOL_MS * 1e-3 * FS)


def claimed(sp):
    out = np.zeros(len(ps), bool)
    sg = ugroup[sp["cluster_index"]]
    for g in range(16):
        st = np.sort(sp["sample_index"][sg == g])
        m = pg == g
        if len(st) > 1 and m.any():
            i = np.clip(np.searchsorted(st, ps[m]), 1, len(st) - 1)
            out[m] = np.minimum(np.abs(st[i] - ps[m]), np.abs(st[i - 1] - ps[m])) <= tol
    return out


def short_isi(sp):
    out = np.zeros(len(sp), bool)
    for u in np.unique(sp["cluster_index"]):
        idx = np.flatnonzero(sp["cluster_index"] == u)
        order = idx[np.argsort(sp["sample_index"][idx])]
        out[order[1:]] = np.diff(sp["sample_index"][order]) < RP
    return out


arms = [("joint", 0.8, 5)] + [("per tetrode", fl, mf) for fl in a.floors for mf in a.failures]
print(f"{a.t0:.0f} to {a.t1:.0f} s, templates from window {w} ({len(ids)} units); {len(ps)} clean events >= 10 MAD")
print(f"{'arm':28s} {'class':9s} {'spikes':>7s} {'claimed%':>9s} {'isi<1.5ms%':>11s} {'cos_med':>8s} {'wall s':>7s}")
for how, fl, mf in arms:
    kw = {"amplitudes": [fl, float("inf")], "max_failures": mf}
    t = time.time()
    fn = find_spikes_from_templates if how == "joint" else find_spikes_per_group
    sp = fn(seg, tmpl, method="circus-omp", method_kwargs=kw, job_kwargs=JOB)
    dt = time.time() - t
    c, v = claimed(sp), short_isi(sp)
    _, r = per_spike_fit(sp, tmpl, seg)
    sg = ugroup[sp["cluster_index"]]
    for name, gs in CLASSES.items():
        m, me = np.isin(sg, gs), np.isin(pg, gs)
        print(f"{how + f', amp>={fl}, mf={mf}':28s} {name:9s} {m.sum():7d} {c[me].mean() * 100:9.1f} "
              f"{v[m].mean() * 100:11.2f} {np.nanmedian(r[m]):8.3f} {dt:7.1f}")
