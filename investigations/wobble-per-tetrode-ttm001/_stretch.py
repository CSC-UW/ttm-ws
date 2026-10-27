"""Shared pieces of the joint against per-tetrode matcher test on TTM-001 (paths literal, as run on tononi-2).

The bank is the one entering window 55 of June's carry-forward, rebuilt from
June's spike trains by `investigations/wobble-fullsize-ttm001/01_span_bank.py`
(`bank_w55.pkl`, 80 units, every tetrode with 2 to 9 units). Matcher settings
are June's: circus-omp with amplitudes 0.8 and up; wobble with an objective
threshold of 0.35 times the bank's median template energy (one absolute
threshold, the same for every tetrode when run per tetrode), then a cosine
gate of 0.60 on each spike's raw snippet against its template. 1 s chunks,
16 jobs. Joint is `find_spikes_from_templates`, per tetrode is
`ttm.matching.find_spikes_per_group`.
"""
import pickle
import sys
import time

import numpy as np
import spikeinterface as si
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _mp_common import (
    per_spike_fit,
    templates_from_dense,
    tsq_median,
    wobble_method_kwargs,
)
from ttm.matching import find_spikes_per_group

si.set_global_job_kwargs(progress_bar=False)

SESSION = "/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52"
TE = f"{SESSION}/sortings_seed42_pcafix/track_eval"
BINARY = f"{TE}/mp_long_s2000_d170000/binary"
BANK = f"{TE}/wobble_span/bank_w55.pkl"
TEMPLATE_AMP = f"{TE}/wobble_span/template_amp_w55.npy"
COVERAGE = f"{TE}/mp_long_s2000_d170000/spike_coverage.npz"
TABLE = f"{SESSION}/bad_periods_crop_s2000_d170000/bad_periods.csv"
OUT = f"{TE}/wobble_pertet"
FS = 30000.0
NBEFORE = 30
AMP = {"amplitudes": [0.8, float("inf")]}
WOBBLE_FACTOR = 0.35
GATE_R = 0.60
JOB = {"n_jobs": 16, "chunk_duration": "1s", "progress_bar": False}


def bank_templates(rec):
    """(Templates of the window-55 bank, tetrode of each template) in the recording's channel order."""
    with open(BANK, "rb") as f:
        b = pickle.load(f)
    ids = sorted(b["bank"])
    groups = np.asarray(rec.get_property("group")).astype(int)
    dense = np.stack([b["bank"][u] for u in ids]).astype(np.float32)
    mask = np.stack([groups == b["group"][u] for u in ids])
    tmpl = templates_from_dense(dense, mask, NBEFORE, np.asarray(ids), rec.channel_ids)
    return tmpl, np.array([b["group"][u] for u in ids])


def method_kwargs(matcher, tmpl):
    if matcher == "circus-omp":
        return AMP
    return wobble_method_kwargs(tmpl, threshold=WOBBLE_FACTOR * tsq_median(tmpl))


def match(seg, tmpl, matcher, per_tetrode, job=JOB):
    """(spikes, wall seconds) of one matcher on ``seg``, jointly or per tetrode."""
    find = find_spikes_per_group if per_tetrode else find_spikes_from_templates
    t = time.time()
    sp = find(seg, tmpl, method=matcher, method_kwargs=method_kwargs(matcher, tmpl), job_kwargs=job)
    return sp, time.time() - t


def gate(sp, tmpl, seg):
    """(spikes passing the cosine gate, a, r of the kept spikes, wall seconds of the gate)."""
    t = time.time()
    a, r = per_spike_fit(sp, tmpl, seg)
    keep = r >= GATE_R
    return sp[keep], a[keep], r[keep], time.time() - t
