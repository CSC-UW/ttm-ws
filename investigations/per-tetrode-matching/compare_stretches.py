"""circus-omp on all 64 channels at once against circus-omp per tetrode, on stretches of TTM-003 48 h.

Templates as the carry-forward would use them on the stretch: every unit of
`mp_baseline/assembled_reseed_c12` with 100 spikes or more in the 1800 s window
before the stretch's window, built from those spikes. Each stretch is matched
as the carry runs it (amplitudes 0.8 and up, 1 s chunks, 16 jobs) jointly
(`find_spikes_from_templates`) and per tetrode (`ttm.matching.find_spikes_per_group`),
on the unmasked binary and, with --masked, also through `mask_bad_periods`.
Per arm: spikes on TT5 to TT12 and on the other tetrodes, the share of the
baseline's detected events at 10 MAD or more in clean time that a unit of the
same tetrode claims within 0.5 ms (per tetrode and pooled), where in the 1 s
chunk the TT5 to TT12 spikes fall, and wall time. Then joint against per
tetrode spike for spike (same unit within 3 samples).

Run: .venv/bin/python investigations/per-tetrode-matching/compare_stretches.py T0 T1 [--masked]
"""
import argparse
import sys
import time

import numpy as np
import spikeinterface as si
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _mp_common import build_templates_object
from _wobble_eval import TOL_MS
from ttm import bad_periods as bp
from ttm.matching import find_spikes_per_group

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
FS = 30000.0
WIN = int(1800 * FS)
AMP = {"amplitudes": [0.8, float("inf")]}
JOB = {"n_jobs": 16, "chunk_duration": "1s", "progress_bar": False}
GOOD = range(4, 12)
ap = argparse.ArgumentParser()
ap.add_argument("t0", type=float)
ap.add_argument("t1", type=float)
ap.add_argument("--masked", action="store_true")
a = ap.parse_args()
f0, f1 = int(a.t0 * FS), int(a.t1 * FS)

rec = si.load(f"{S}/mp_baseline/binary")
table = bp.read_table(TABLE)
base = si.load(f"{S}/mp_baseline/assembled_reseed_c12")
w = f0 // WIN - 1
lo, hi = max(w, 0) * WIN, (max(w, 0) + 1) * WIN
ids = [u for u in base.unit_ids if ((st := base.get_unit_spike_train(u)) >= lo).sum() - (st >= hi).sum() >= 100]
sw = base.select_units(ids).frame_slice(lo, hi)
rw = rec.frame_slice(lo, hi)
rw.reset_times()
sw.register_recording(rw)
tmpl, _ = build_templates_object(sw, rw, with_snr=False, n_jobs=32)
groups = np.asarray(rec.get_property("group")).astype(int)
ugroup = np.array([groups[np.flatnonzero(m)[0]] for m in tmpl.sparsity.mask])
flagged = table[(table.start_s < a.t1) & (table.end_s > a.t0)]
print(f"{a.t0:.0f} to {a.t1:.0f} s; templates from window {max(w, 0)} ({len(ids)} units, "
      f"{int(np.isin(ugroup, GOOD).sum())} on TT5-TT12); tetrodes with bad periods here: "
      f"{sorted(set(flagged.tetrode), key=lambda s: int(s[2:]))}")

z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
k = (z["peak_sample"] >= f0) & (z["peak_sample"] < f1) & (z["amp_mad"] >= 10)
ps, pg = z["peak_sample"][k] - f0, z["peak_group"][k].astype(int)
clean = np.ones(len(ps), bool)
for g in range(16):
    m = pg == g
    clean[m] = ~bp.bad_mask(table, g, ps[m] + f0, FS)
ps, pg = ps[clean], pg[clean]
tol = int(TOL_MS * 1e-3 * FS)


def claimed(spikes):
    out = np.zeros(len(ps), bool)
    sg = ugroup[spikes["cluster_index"]]
    for g in range(16):
        st = np.sort(spikes["sample_index"][sg == g])
        m = pg == g
        if len(st) == 0 or not m.any():
            continue
        i = np.clip(np.searchsorted(st, ps[m]), 1, len(st) - 1) if len(st) > 1 else np.zeros(m.sum(), int)
        d = np.minimum(np.abs(st[i] - ps[m]), np.abs(st[np.maximum(i - 1, 0)] - ps[m]))
        out[m] = d <= tol
    return out


arms = [("joint", "unmasked"), ("per tetrode", "unmasked")]
if a.masked:
    arms += [("joint", "masked"), ("per tetrode", "masked")]
recs = {"unmasked": rec}
if a.masked:
    recs["masked"] = bp.mask_bad_periods(rec, table)
res = {}
for how, which in arms:
    seg = recs[which].frame_slice(f0, f1)
    seg.reset_times()
    t = time.time()
    if how == "joint":
        sp = find_spikes_from_templates(seg, tmpl, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    else:
        sp = find_spikes_per_group(seg, tmpl, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    res[(how, which)] = (sp, time.time() - t)

print(f"\n{'arm':24s} {'TT5-TT12':>9s} {'other':>9s} {'claimed % (10 MAD, clean)':>27s} {'wall s':>7s}")
for (how, which), (sp, dt) in res.items():
    good = np.isin(ugroup[sp["cluster_index"]], GOOD)
    c = claimed(sp)
    print(f"{how + ', ' + which:24s} {good.sum():9d} {(~good).sum():9d} {c.mean() * 100:27.1f} {dt:7.1f}")

print(f"\nclaimed % per tetrode ({len(ps)} clean events at 10 MAD or more)")
print("tetrode  n_events  " + "  ".join(f"{h[:5]}/{w_[:4]}" for h, w_ in res))
cl = {key: claimed(sp) for key, (sp, _) in res.items()}
for g in range(16):
    m = pg == g
    if m.sum():
        print(f"  TT{g + 1:<4d} {m.sum():8d}  " + "  ".join(f"{cl[key][m].mean() * 100:10.1f}" for key in res))

print("\nTT5-TT12 spikes by tenth of the 1 s chunk")
for key, (sp, _) in res.items():
    f = sp["sample_index"][np.isin(ugroup[sp["cluster_index"]], GOOD)]
    print(f"  {key[0] + ', ' + key[1]:24s} {np.bincount((f % int(FS)) // int(FS / 10), minlength=10).tolist()}")


def agreement(x, y):
    both = 0
    for u in range(len(tmpl.unit_ids)):
        sx = np.sort(x["sample_index"][x["cluster_index"] == u])
        sy = np.sort(y["sample_index"][y["cluster_index"] == u])
        if len(sx) and len(sy):
            i = np.clip(np.searchsorted(sy, sx), 1, max(len(sy) - 1, 1)) if len(sy) > 1 else np.zeros(len(sx), int)
            d = np.minimum(np.abs(sy[i] - sx), np.abs(sy[np.maximum(i - 1, 0)] - sx))
            both += int((d <= 3).sum())
    return both


for which in recs:
    x, y = res[("joint", which)][0], res[("per tetrode", which)][0]
    both = agreement(x, y)
    print(f"\n{which}: joint {len(x)} spikes, per tetrode {len(y)}; in both (same unit, 3 samples) {both}; "
          f"joint only {len(x) - both}, per tetrode only {len(y) - both}")
    for g in range(16):
        xg, yg = x[ugroup[x["cluster_index"]] == g], y[ugroup[y["cluster_index"]] == g]
        if len(xg) or len(yg):
            b = agreement(xg, yg)
            print(f"  TT{g + 1:<4d} joint {len(xg):7d}  per tetrode {len(yg):7d}  both {b:7d}")
