import dataclasses

import numpy as np
import pytest

pytest.importorskip("spikeinterface.sortingcomponents")
import spikeinterface as si
from ttm import carry as cr

FS = 30000.0
P = cr.CarryParams()


def _template(seed=0, n_ch=4):
    rng = np.random.default_rng(seed)
    t = np.arange(P.nbefore + P.nafter) - P.nbefore
    shape = -np.exp(-(t / 4.0) ** 2) + 0.4 * np.exp(-((t - 12) / 8.0) ** 2)
    return (shape[:, None] * rng.uniform(0.5, 1.5, size=n_ch)[None, :] * 12).astype(np.float32)


def _tetrode(n=300_000, frames=(), amps=(), T=None, seed=1):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, 4)).astype(np.float32)
    if len(frames):
        cr.add_spikes(x, np.asarray(frames, np.int64), np.zeros(len(frames), np.int64),
                      np.asarray(amps, np.float32), T[None], P.nbefore)
    return x


def _run(x, frames, amps, T, *, pulses=None, bad=(), params=P):
    frames = np.asarray(frames, np.int64)
    sigma = cr.mf_noise(x, T[None])
    return cr.tetrode_window(x, 0, frames, np.zeros(frames.size, np.int64), np.asarray(amps, np.float64), T[None],
                             sigma, params, fs=FS, pulses=pulses, bad=bad)


def test_admission_floors_bad_periods_and_refractory() -> None:
    T = _template()
    real = np.arange(1000, 290_000, 1500)
    x = _tetrode(frames=real, amps=np.ones(real.size), T=T)
    noise_fits = np.arange(1700, 290_000, 15000)
    twins = real[:10] + 20
    frames = np.r_[real, noise_fits, twins]
    amps = np.r_[np.ones(real.size), np.full(noise_fits.size, 0.6), np.full(twins.size, 0.1)]
    o = np.argsort(frames, kind="stable")
    r = _run(x, frames[o], amps[o], T, bad=[(100_000, 110_000)])
    reason = np.empty_like(r["reason"])
    reason[o] = r["reason"]
    n = real.size
    in_bad = (real >= 100_000) & (real < 110_000)
    assert (reason[:n][in_bad] == cr.REASON_BAD_PERIOD).all()
    assert (reason[:n][~in_bad] == cr.REASON_ACCEPTED).mean() > 0.98
    assert (reason[n:n + noise_fits.size] != cr.REASON_ACCEPTED).all()
    assert (reason[n + noise_fits.size:] != cr.REASON_ACCEPTED).all()


def test_refractory_keeps_the_higher_snr_spike() -> None:
    T = _template()
    real = np.arange(1000, 290_000, 1500)
    x = _tetrode(frames=real, amps=np.ones(real.size), T=T)
    early = real[:20] - 30
    frames = np.sort(np.r_[real, early])
    amps = np.ones(frames.size)
    r = _run(x, frames, amps, T)
    acc = frames[r["reason"] == cr.REASON_ACCEPTED]
    assert np.isin(real[:20], acc).all()
    assert not np.isin(early, acc).any()
    assert (r["reason"][np.isin(frames, early)] != cr.REASON_ACCEPTED).all()


def test_rebuild_blends_toward_the_window_mean() -> None:
    T = _template()
    T_true = T * 1.3
    real = np.arange(1000, 290_000, 500)
    x = _tetrode(frames=real, amps=np.ones(real.size), T=T_true)
    r = _run(x, real, np.full(real.size, 1.3), T)
    assert r["n_used"][0] == P.rebuild_max_spikes
    expected = 0.8 * T + 0.2 * T_true
    np.testing.assert_allclose(r["newT"][0], expected, atol=0.25)
    frozen = _run(x, real, np.full(real.size, 1.3), T, params=dataclasses.replace(P, blend=0.0))
    np.testing.assert_allclose(frozen["newT"][0], T, atol=1e-5)


def test_relaxed_rule_inside_pulses_only_for_units_that_fire_through_them() -> None:
    T = _template()
    on = np.arange(10_000, 290_000, 30_000)
    off = on + int(0.207 * FS)
    outside = np.arange(2000, 290_000, 3000)
    outside = outside[~cr.in_intervals(outside, on, off + 3000)]
    inside = np.concatenate([np.arange(o + 600, f, 200) for o, f in zip(on, off)])
    small = 0.45
    x = _tetrode(frames=np.r_[outside, inside], amps=np.r_[np.ones(outside.size), np.full(inside.size, small)], T=T)
    frames = np.sort(np.r_[outside, inside])
    amps = np.where(np.isin(frames, inside), small, 1.0)
    without = _run(x, frames, amps, T)
    with_pulses = _run(x, frames, amps, T, pulses=(on, off))
    is_in = np.isin(frames, inside)
    assert (without["reason"][is_in] == cr.REASON_ACCEPTED).mean() < 0.2
    assert with_pulses["driven"] == [0]
    assert (with_pulses["reason"][is_in] == cr.REASON_ACCEPTED).mean() > 0.9


def _binary(tmp_path, n=4 * 60_000):
    T0, T1 = _template(0), _template(3)
    rng = np.random.default_rng(5)
    x = rng.normal(size=(n, 8)).astype(np.float32) * 4
    f0 = np.arange(500, n - 500, 700)
    f1 = np.arange(850, n - 500, 1100)
    for f, T, ch in ((f0, T0, slice(0, 4)), (f1, T1, slice(4, 8))):
        y = np.zeros((n, 4), np.float32)
        cr.add_spikes(y, f, np.zeros(f.size, np.int64), np.full(f.size, 3.0, np.float32), T[None], P.nbefore)
        x[:, ch] += y
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    rec.set_property("group", np.repeat([0, 1], 4))
    rec = rec.save(folder=tmp_path / "binary", format="binary")
    return tmp_path / "binary", rec, cr.Bank(np.array([0, 1]), np.array([0, 1]), np.stack([T0, T1]) * 3)


def test_carry_writes_each_window_and_resumes(tmp_path) -> None:
    path, rec, bank = _binary(tmp_path)
    params = dataclasses.replace(P, window_s=2.0, whiten_s=2.0, reseed_every=0, rebuild_min_accepted=20,
                                 rebuild_min_selected=10, rel_min_lead=10)
    W = cr.whitening_matrix(rec, params)
    groups = np.asarray(rec.get_property("group")).astype(int)
    bank_w = bank.whitened(W, groups)
    out = tmp_path / "carry"
    cr.run_carry(path, W, bank_w, out, params=params, windows=range(0, 2), n_jobs=1, workers=1)
    first = {f.name: np.load(f)["frame"] for f in out.glob("spikes_w*.npz")}
    assert sorted(first) == ["spikes_w00.npz", "spikes_w01.npz"]
    acc = np.load(out / "spikes_w00.npz")["reason"] == 0
    assert acc.sum() > 100
    mtime = (out / "spikes_w00.npz").stat().st_mtime_ns
    cr.run_carry(path, W, bank_w, out, params=params, windows=range(0, 4), n_jobs=1, workers=1)
    assert (out / "spikes_w00.npz").stat().st_mtime_ns == mtime
    assert (out / "spikes_w03.npz").exists()
    full = tmp_path / "carry_full"
    cr.run_carry(path, W, bank_w, full, params=params, windows=range(0, 4), n_jobs=1, workers=1)
    for w in range(4):
        a, b = np.load(out / f"spikes_w{w:02d}.npz"), np.load(full / f"spikes_w{w:02d}.npz")
        np.testing.assert_array_equal(a["frame"], b["frame"])
        np.testing.assert_array_equal(a["reason"], b["reason"])
    s = cr.export_sorting(out, tmp_path / "sorting", FS)
    assert s.get_num_units() == 2 and list(s.get_property("group")) == [0, 1]


def test_run_info_refuses_other_parameters(tmp_path) -> None:
    path = tmp_path / "carry_params.json"
    cr.write_run_info(path, P, binary="b")
    cr.write_run_info(path, P, binary="b")
    with pytest.raises(ValueError, match="different parameters"):
        cr.write_run_info(path, dataclasses.replace(P, blend=0.5))
