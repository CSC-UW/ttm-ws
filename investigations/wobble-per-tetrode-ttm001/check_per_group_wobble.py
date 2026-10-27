"""Does `ttm.matching.find_spikes_per_group` run wobble correctly?

`find_spikes_per_group` is tested with circus-omp only
(`packages/ttm/tests/test_matching.py`). Here, with wobble at one absolute
threshold (0.35 times the median template energy of the whole bank):

1. Synthetic, as in the tests: three tetrodes of four channels, two units
   each, Gaussian noise. Per group against wobble run separately on each
   tetrode's channels with that tetrode's templates and the same threshold;
   then an artifact on one tetrode, which must leave the others' spikes
   unchanged.
2. TTM-001: 20 s of the stretch with the window-55 bank. Per group with 16
   jobs against the separate per-tetrode runs, and against per group with
   one job.

Equality is of the (sample, unit) pairs.

Run: .venv/bin/python investigations/wobble-per-tetrode-ttm001/check_per_group_wobble.py T0
"""
import sys

import numpy as np
import spikeinterface as si
from _stretch import BINARY, FS, JOB, bank_templates, method_kwargs
from spikeinterface.core import Templates
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates
from ttm.matching import find_spikes_per_group, group_templates, template_groups

T, NB, N_TET = 90, 30, 3
SYN_JOB = {"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False}


def key(sp):
    return sorted(zip(sp["sample_index"].tolist(), sp["cluster_index"].tolist()))


def by_hand(rec, tmpl, mk, job):
    groups = np.asarray(rec.get_property("group")).astype(int)
    ug = template_groups(tmpl, groups)
    out = []
    for g in np.unique(ug):
        ch = np.flatnonzero(groups == g)
        units = np.flatnonzero(ug == g)
        s = find_spikes_from_templates(rec.select_channels(np.asarray(rec.channel_ids)[ch]),
                                       group_templates(tmpl, units, ch), method="wobble", method_kwargs=mk,
                                       job_kwargs=job)
        s["cluster_index"] = units[s["cluster_index"]]
        s["channel_index"] = ch[s["channel_index"]]
        out.append(s)
    return np.concatenate(out)


def synthetic(artifact_on=None, seed=0):
    t = np.arange(T) - NB
    dense = np.zeros((2 * N_TET, T, 4 * N_TET), np.float32)
    mask = np.zeros((2 * N_TET, 4 * N_TET), bool)
    for u in range(2 * N_TET):
        amp = -120.0 - 30 * (u % 2)
        base = amp * np.exp(-(t**2) / 8.0) - 0.4 * amp * np.exp(-((t - 6) ** 2) / 20.0)
        g, c = u // 2, u % 4
        dense[u, :, 4 * g + c] = base
        dense[u, :, 4 * g + (c + 1) % 4] = 0.4 * base
        mask[u, 4 * g: 4 * g + 4] = True
    tmpl = Templates(templates_array=np.stack([dense[u][:, mask[u]] for u in range(len(mask))]),
                     sampling_frequency=FS, nbefore=NB, is_in_uV=False, sparsity_mask=mask,
                     channel_ids=np.array([f"CH{c}" for c in range(4 * N_TET)]), unit_ids=np.arange(2 * N_TET))
    rng = np.random.default_rng(seed)
    n = int(6 * FS)
    x = rng.normal(0, 10.0, (n, 4 * N_TET)).astype(np.float32)
    for u in range(len(dense)):
        for t0 in np.sort(rng.choice(np.arange(200, n - 200, 120), size=90, replace=False)):
            x[t0 - NB: t0 - NB + T] += dense[u]
    if artifact_on is not None:
        a = np.arange(n // 3, 2 * n // 3)
        burst = 400.0 * np.sin(2 * np.pi * 900.0 * a / FS) * (rng.random(len(a)) < 0.6)
        x[a, 4 * artifact_on: 4 * artifact_on + 4] += burst[:, None].astype(np.float32)
    rec = si.NumpyRecording([x], FS, channel_ids=np.array([f"CH{c}" for c in range(4 * N_TET)]))
    rec.set_property("group", np.repeat(np.arange(N_TET), 4))
    return rec, tmpl


rec, tmpl = synthetic()
mk = method_kwargs("wobble", tmpl)
pg = find_spikes_per_group(rec, tmpl, method="wobble", method_kwargs=mk, job_kwargs=SYN_JOB)
hand = by_hand(rec, tmpl, mk, SYN_JOB)
joint = find_spikes_from_templates(rec, tmpl, method="wobble", method_kwargs=mk, job_kwargs=SYN_JOB)
print(f"synthetic, threshold {mk['parameters']['threshold']:.0f}: per group {len(pg)} spikes, separate per "
      f"tetrode {len(hand)}, joint {len(joint)}; per group equals separate: {key(pg) == key(hand)}; channel "
      f"indices equal: {sorted(pg['channel_index'].tolist()) == sorted(hand['channel_index'].tolist())}")
noisy, _ = synthetic(artifact_on=0)
pn = find_spikes_per_group(noisy, tmpl, method="wobble", method_kwargs=mk, job_kwargs=SYN_JOB)
jn = find_spikes_from_templates(noisy, tmpl, method="wobble", method_kwargs=mk, job_kwargs=SYN_JOB)
others = lambda s: s[s["cluster_index"] >= 2]
print(f"artifact on tetrode 0: other tetrodes' spikes unchanged per group: {key(others(pg)) == key(others(pn))} "
      f"({len(others(pg))} and {len(others(pn))}); jointly: {key(others(joint)) == key(others(jn))} "
      f"({len(others(joint))} and {len(others(jn))})")

t0 = float(sys.argv[1])
full = si.load(BINARY)
seg = full.frame_slice(int(t0 * FS), int((t0 + 20) * FS))
seg.reset_times()
tmpl, ugroup = bank_templates(full)
mk = method_kwargs("wobble", tmpl)
pg = find_spikes_per_group(seg, tmpl, method="wobble", method_kwargs=mk, job_kwargs=JOB)
hand = by_hand(seg, tmpl, mk, JOB)
serial = find_spikes_per_group(seg, tmpl, method="wobble", method_kwargs=mk,
                               job_kwargs={"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False})
print(f"\nTTM-001 {t0:.0f} to {t0 + 20:.0f} s, {len(tmpl.unit_ids)} units, threshold "
      f"{mk['parameters']['threshold']:.0f}: per group {len(pg)} spikes, separate per tetrode {len(hand)}, "
      f"per group with one job {len(serial)}; per group equals separate: {key(pg) == key(hand)}; "
      f"16 jobs equal one: {key(pg) == key(serial)}")
print("spikes per tetrode, per group: " + str(np.bincount(ugroup[pg["cluster_index"]], minlength=16).tolist()))
