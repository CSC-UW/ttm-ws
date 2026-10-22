"""Tests for ``ttm.capture``: the shape-based capture of low-amplitude spikes.

Synthetic tetrode traces with two known cells (one large, one small, different wire
profiles) check that the excised noise estimate ignores spikes, that detection finds every
planted trough, that a half-amplitude copy of a template scores cosine 1 and is assigned
to it, and that the refractory rule keeps the better of two colliding events. The last
test runs the module on 60 s of a real TTM-003 session and is skipped without the data.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from ttm import capture as cp

FS = 30000.0


def _template(amp, profile, width=3.0):
    t = np.arange(cp.N_SAMP) - cp.NBEFORE
    w = -amp * np.exp(-(t ** 2) / (2 * width ** 2)) + 0.3 * amp * np.exp(-((t - 12) ** 2) / 50.0)
    return (w[:, None] * np.asarray(profile)[None, :]).astype(np.float32)


def _synthetic(seconds=4.0, seed=0):
    rng = np.random.default_rng(seed)
    n = int(seconds * FS)
    x = rng.normal(0, 10.0, size=(n, 4)).astype(np.float32)
    ta = _template(200.0, [1.0, 0.5, 0.2, 0.1])
    tb = _template(90.0, [0.1, 0.2, 0.6, 1.0])
    fa = np.sort(rng.choice(np.arange(200, n - 200, 90), 300, replace=False))
    fb = np.sort(rng.choice(np.arange(245, n - 200, 90), 300, replace=False))
    scale_a = rng.uniform(0.4, 1.0, fa.size)  # burst attenuation down to 40 %
    for f, s in zip(fa, scale_a):
        x[f - cp.NBEFORE:f + cp.NAFTER] += s * ta
    for f in fb:
        x[f - cp.NBEFORE:f + cp.NAFTER] += tb
    return x, ta, tb, fa, fb


def test_excised_noise_ignores_spikes():
    x, ta, tb, fa, fb = _synthetic()
    plain = cp.noise_mad(x)
    excised, frac = cp.noise_mad_excised(x)
    assert 0 < frac < 0.5
    assert np.all(excised <= plain + 1e-6)
    assert np.allclose(excised, 10.0, atol=0.6)


def test_detect_finds_planted_troughs():
    x, ta, tb, fa, fb = _synthetic()
    noise, _ = cp.noise_mad_excised(x)
    frames, amp = cp.detect_troughs(x, noise, threshold_mad=3.5)
    for planted in (fa, fb):
        j = np.searchsorted(frames, planted)
        near = np.minimum(np.abs(frames[np.clip(j - 1, 0, frames.size - 1)] - planted),
                          np.abs(frames[np.clip(j, 0, frames.size - 1)] - planted))
        assert (near <= 3).mean() > 0.97
    assert np.all(amp >= 3.5)


def test_half_amplitude_copy_scores_cosine_one_and_is_assigned():
    x, ta, tb, fa, fb = _synthetic()
    templates = np.stack([ta, tb])
    clean = np.zeros((3000, 4), np.float32)
    clean[1000 - cp.NBEFORE:1000 + cp.NAFTER] = 0.5 * ta
    clean[2000 - cp.NBEFORE:2000 + cp.NAFTER] = tb
    r = cp.cosine_to_templates(clean, np.array([1000, 2000]), templates, s_max=0)
    assert r[0, 0] == pytest.approx(1.0, abs=1e-5)
    assert r[1, 1] == pytest.approx(1.0, abs=1e-5)
    assert r[0, 1] < 0.8 and r[1, 0] < 0.8
    a = cp.assign_events(np.array([1000, 2000]), np.array([10.0, 6.0]), r, gate_r=0.6)
    assert list(a.label) == [0, 1]
    assert list(a.reason) == ["accepted", "accepted"]


def test_assignment_recovers_attenuated_spikes_on_noisy_traces():
    x, ta, tb, fa, fb = _synthetic()
    noise, _ = cp.noise_mad_excised(x)
    frames, amp = cp.detect_troughs(x, noise, threshold_mad=3.5)
    trains = {0: fa[::2], 1: fb[::2]}  # templates from half the spikes
    templates, ids = cp.templates_from_trains(x, trains)
    r = cp.cosine_to_templates(x, frames, templates)
    a = cp.assign_events(frames, amp, r, gate_r=0.6)
    for planted, lab in ((fa, 0), (fb, 1)):
        got = a.train(lab)
        j = np.searchsorted(got, planted)
        near = np.minimum(np.abs(got[np.clip(j - 1, 0, got.size - 1)] - planted),
                          np.abs(got[np.clip(j, 0, got.size - 1)] - planted))
        assert (near <= 3).mean() > 0.95, lab
    pool = a.pool
    assert pool.size < 0.2 * frames.size


def test_refractory_rule_keeps_better_match():
    frames = np.array([100, 120, 500, 2000, 2010, 2020])
    score = np.array([0.9, 0.95, 0.8, 0.7, 0.99, 0.8])
    keep = cp.resolve_refractory(frames, score, refr_frames=45)
    assert list(keep) == [False, True, True, False, True, False]


def test_ccg_verdict_and_amplitude_split_on_one_cell_versus_two():
    rng = np.random.default_rng(1)
    n_frames = int(600 * FS)
    one = np.sort(rng.choice(n_frames, 20000, replace=False))
    one = one[np.r_[True, np.diff(one) > 90]]  # a 3 ms refractory period
    amp = rng.uniform(5, 15, one.size)
    split = cp.amplitude_split(one, amp, n_frames)
    assert split["ratio"] < 0.2
    other = np.sort(rng.choice(n_frames, 15000, replace=False))
    v = cp.ccg_verdict(one, other)
    assert v["verdict"] == "distinct"
    assert cp.ccg_verdict(one[::2], one[1::2])["verdict"] == "duplicate"
    c, nv = cp.rp_contamination(one, n_frames)
    assert nv == 0 and c == 0.0


@pytest.mark.requires_nfs
def test_tt6_in_pulse_spikes_are_captured_on_real_data():
    import spikeinterface as si
    from ttm.subject import find_subject

    try:
        exp = find_subject("TTM-003").experiment("TTM-003-OPTO1")
        out = exp.sorting_dir("asis", subdir="sortings_opto")
    except Exception:
        pytest.skip("TTM-003 manifests or config not available")
    if not (out / "preprocessed").exists():
        pytest.skip("OPTO1 preprocessed recording not mounted")
    rec = si.load(str(out / "preprocessed"))
    sorting = si.load(str(out / "sorting"))
    pulses = pathlib.Path(exp.manifest_path).parents[1] / "results" / f"stimulation_pulses_{exp.session}.csv"
    import pandas as pd
    pp = pd.read_csv(pulses)
    pp = pp[~pp.glitch]
    a, b = int(26 * 60 * FS), int(29 * 60 * FS)   # train 6, 9.8 mW
    x = cp.load_tetrode_traces(rec, 5, a, b)
    noise, _ = cp.noise_mad_excised(x)
    frames, amp = cp.detect_troughs(x, noise)
    groups = np.asarray(sorting.get_property("group")).astype(int)
    trains = {int(u): sorting.get_unit_spike_train(u).astype(np.int64) - a
              for u, g in zip(sorting.unit_ids, groups) if g == 5}
    trains = {u: t[(t >= 0) & (t < b - a)] for u, t in trains.items()}
    templates, ids = cp.templates_from_trains(x, trains)
    r = cp.cosine_to_templates(x, frames, templates)
    asg = cp.assign_events(frames, amp, r)
    cell = asg.train(ids.index(21)) / FS + a / FS
    sel = pp[(pp.onset_s > a / FS) & (pp.onset_s < b / FS)]
    dur = float(np.median(sel.offset_s - sel.onset_s))
    in_pulse = sum(((cell >= on + 0.01) & (cell < on + dur)).sum() for on in sel.onset_s) / (len(sel) * (dur - 0.01))
    assert in_pulse > 80, in_pulse
