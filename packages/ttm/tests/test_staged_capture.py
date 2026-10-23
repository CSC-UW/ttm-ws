"""Tests for ``ttm.staged_capture``: small-spike capture on the residual of a full-size sorting.

Synthetic tetrode traces hold a large cell whose spikes the "full-size sorting" has, and a
second cell that fires full-size spikes outside the pulses (also sorted) and half-size spikes
inside them (not sorted). Subtracting the sorted spikes must leave the half-size spikes as
residual events, the noise-corrected gate must offer them to the second cell's template, the
light test must admit the second cell and not the first, and acceptance must respect each
unit's refractory period against its sorted spikes.
"""

from __future__ import annotations

import numpy as np

from ttm import capture as cp
from ttm import staged_capture as sc

FS = 30000.0


def _template(amp, profile, width=3.0):
    t = np.arange(cp.N_SAMP) - cp.NBEFORE
    w = -amp * np.exp(-(t ** 2) / (2 * width ** 2)) + 0.3 * amp * np.exp(-((t - 12) ** 2) / 50.0)
    return (w[:, None] * np.asarray(profile)[None, :]).astype(np.float32)


def _scene(seed=0):
    rng = np.random.default_rng(seed)
    n = int(30 * FS)
    x = rng.normal(0, 10.0, size=(n, 4)).astype(np.float32)
    ta = _template(200.0, [1.0, 0.5, 0.2, 0.1])
    tb = _template(120.0, [0.1, 0.2, 0.6, 1.0])
    onsets = np.arange(2.0, 28.0, 1.0)
    offsets = onsets + 0.2
    in_pulse = lambda f: ((f / FS)[:, None] >= onsets[None, :] + 0.01).any(axis=1) & \
        ((f / FS)[:, None] < offsets[None, :]).any(axis=1)
    fa = np.sort(rng.choice(np.arange(300, n - 300, 120), 150, replace=False))
    fa = fa[~in_pulse(fa)]                       # the large cell is silenced by the light
    fb_out = np.sort(rng.choice(np.arange(345, n - 300, 120), 60, replace=False))
    fb_out = fb_out[~in_pulse(fb_out)]
    fb_in = np.concatenate([np.arange(int((o + 0.012) * FS), int((o + 0.19) * FS), 240) for o in onsets])
    for f in fa:
        x[f - cp.NBEFORE:f + cp.NAFTER] += ta
    for f in fb_out:
        x[f - cp.NBEFORE:f + cp.NAFTER] += tb
    for f in fb_in:
        x[f - cp.NBEFORE:f + cp.NAFTER] += 0.5 * tb
    return x, np.stack([ta, tb]), fa, fb_out, fb_in, onsets, offsets


def test_residual_keeps_only_unsorted_spikes():
    x, T, fa, fb_out, fb_in, on, off = _scene()
    frames = np.concatenate([fa, fb_out])
    labels = np.concatenate([np.zeros(fa.size, int), np.ones(fb_out.size, int)])
    order = np.argsort(frames)
    res, a = sc.subtract_sorted(x, frames[order], labels[order], T)
    assert np.nanmedian(a) > 0.9 and np.nanmedian(a) < 1.1
    noise, _ = cp.noise_mad_excised(x)
    ev, amp = sc.residual_events(res, noise, frames[order], threshold_mad=4.0)
    hit = sc.near(fb_in, ev, 3)
    assert hit.mean() > 0.95
    assert sc.near(ev, np.sort(fb_in), 3).mean() > 0.8


def test_gate_offers_half_size_spikes_and_light_admits_only_driven_unit():
    x, T, fa, fb_out, fb_in, on, off = _scene()
    frames = np.concatenate([fa, fb_out])
    labels = np.concatenate([np.zeros(fa.size, int), np.ones(fb_out.size, int)])
    order = np.argsort(frames)
    res, _ = sc.subtract_sorted(x, frames[order], labels[order], T)
    noise, _ = cp.noise_mad_excised(x)
    ev, amp = sc.residual_events(res, noise, frames[order], threshold_mad=4.0)
    off_ = sc.offer_events(res, ev, amp, T, noise, gate_rc=0.8)
    got_b = off_.frames[off_.offered & (off_.best == 1)]
    assert sc.near(fb_in, got_b, 3).mean() > 0.9
    trains = {0: np.sort(fa), 1: np.sort(fb_out)}
    resp = {}
    for i in (0, 1):
        new = off_.frames[off_.offered & (off_.best == i)]
        spk = np.sort(np.concatenate([trains[i], new])) / FS
        resp[i] = sc.light_response(spk, on, off)
    assert sc.admit(resp[1], min_ratio=2.0, max_p=1e-3, min_in_hz=5.0)
    assert not sc.admit(resp[0], min_ratio=2.0, max_p=1e-3, min_in_hz=5.0)


def test_accept_respects_host_refractory_and_eligibility():
    off_ = sc.Offer(frames=np.array([100, 130, 1000, 1020, 5000]), amp_mad=np.full(5, 5.0),
                    best=np.array([0, 0, 0, 0, 1]), r=np.array([0.9, 0.8, 0.7, 0.95, 0.9]),
                    rc=np.full(5, 0.9), a=np.full(5, 0.5), offered=np.ones(5, bool),
                    pass_index=np.zeros(5, int))
    host = {7: np.array([120]), 8: np.zeros(0, np.int64)}
    keep = sc.accept_offered(off_, np.array([True, True, True, True, False]), host, [7, 8], refr_frames=45)
    assert keep.tolist() == [False, False, False, True, False]


def test_governing_train_uses_most_recent_and_first_before_any():
    k = sc.governing_train(np.array([0.0, 5.0, 10.0, 11.0, 99.0]), np.array([10.0, 50.0]))
    assert k.tolist() == [0, 0, 0, 0, 1]
