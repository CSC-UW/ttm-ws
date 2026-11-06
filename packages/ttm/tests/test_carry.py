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
    assert with_pulses["sustained"] == [0]
    assert (with_pulses["reason"][is_in] == cr.REASON_ACCEPTED).mean() > 0.9


def _sustained_case(post_ms=(2, 4, 20, 60)):
    """A unit firing half-size spikes through 207 ms pulses, with half-size fits at ``post_ms`` after
    each offset, and full-size spikes elsewhere."""
    T = _template()
    on = np.arange(10_000, 290_000, 30_000)
    off = on + int(0.207 * FS)
    outside = np.arange(2000, 290_000, 3000)
    outside = outside[~cr.in_intervals(outside, on, off + 3000)]
    inside = np.concatenate([np.arange(o + 600, f, 200) for o, f in zip(on, off)])
    after = {ms: off + int(ms * 1e-3 * FS) for ms in post_ms}
    small = np.r_[inside, *after.values()]
    x = _tetrode(frames=np.r_[outside, small], amps=np.r_[np.ones(outside.size), np.full(small.size, 0.45)], T=T)
    frames = np.sort(np.r_[outside, small])
    amps = np.where(np.isin(frames, small), 0.45, 1.0)
    return T, x, frames, amps, (on, off), inside, after


def test_relaxed_span_ends_pulse_relax_post_ms_after_the_offset() -> None:
    T, x, frames, amps, pulses, inside, after = _sustained_case()
    r = _run(x, frames, amps, T, pulses=pulses)
    assert r["sustained"] == [0]
    acc = r["reason"] == cr.REASON_ACCEPTED
    assert acc[np.isin(frames, inside)].mean() > 0.9
    assert acc[np.isin(frames, after[2])].mean() > 0.9 and acc[np.isin(frames, after[4])].mean() > 0.9
    assert not acc[np.isin(frames, after[20])].any() and not acc[np.isin(frames, after[60])].any()
    wide = _run(x, frames, amps, T, pulses=pulses, params=dataclasses.replace(P, pulse_relax_post_ms=100.0))
    assert (wide["reason"][np.isin(frames, after[60])] == cr.REASON_ACCEPTED).mean() > 0.9


def test_without_the_pulse_rule_pulse_times_still_stay_out_of_the_floor() -> None:
    T, x, frames, amps, pulses, inside, _ = _sustained_case(post_ms=())
    off_rule = dataclasses.replace(P, pulse_relax=False)
    r = _run(x, frames, amps, T, pulses=pulses, params=off_rule)
    assert r["sustained"] == []
    assert (r["reason"][np.isin(frames, inside)] == cr.REASON_ACCEPTED).mean() < 0.2
    on = np.arange(10_000, 290_000, 30_000)
    off = on + int(0.207 * FS)
    outside = np.arange(2000, 290_000, 6000)
    outside = outside[~cr.in_intervals(outside, on, off + 3000)]
    inside = np.concatenate([np.arange(o + 600, f, 450) for o, f in zip(on, off)])
    frames = np.sort(np.r_[outside, inside])
    amps = np.where(np.isin(frames, inside), 0.45, 1.0)
    assert np.isin(frames, inside).mean() > 0.5

    def admit(pulses):
        return cr.admit(frames, np.zeros(frames.size, np.int64), amps, np.full(frames.size, 10.0),
                        np.ones(frames.size), 1, off_rule, fs=FS, a=0, n_frames=300_000, pulses=pulses)[0]

    assert (admit((on, off))[np.isin(frames, inside)] == cr.REASON_FLOOR).all()
    assert (admit(None)[np.isin(frames, inside)] == cr.REASON_ACCEPTED).all()


def test_pulse_rule_takes_no_fit_inside_a_bad_period() -> None:
    T, x, frames, amps, pulses, inside, _ = _sustained_case(post_ms=())
    on, off = pulses
    bad = [(int(on[2]), int(off[2]))]
    r = _run(x, frames, amps, T, pulses=pulses, bad=bad)
    in_bad = cr.in_intervals(frames, np.array([bad[0][0]]), np.array([bad[0][1]]))
    assert r["sustained"] == [0]
    assert not (r["reason"][in_bad] == cr.REASON_ACCEPTED).any()
    assert (r["reason"][np.isin(frames, inside) & ~in_bad] == cr.REASON_ACCEPTED).mean() > 0.9


def test_run_info_reads_the_single_pulse_span_of_older_runs(tmp_path) -> None:
    import json

    path = tmp_path / "carry_params.json"
    old = P.to_dict()
    for k in ("pulse_exclude_post_ms", "pulse_relax", "pulse_relax_post_ms"):
        old.pop(k)
    old["pulse_post_ms"] = 100.0
    old = {(f"driven_{k[10:]}" if k.startswith("sustained_") else k): v for k, v in old.items()}
    assert "driven_start_ms" in old
    path.write_text(json.dumps({"params": old, "inputs": {}}))
    cr.write_run_info(path, dataclasses.replace(P, pulse_relax_post_ms=100.0))
    with pytest.raises(ValueError, match="different parameters"):
        cr.write_run_info(path, P)
    assert cr.record_readmit(path, P) == {"pulse_relax_post_ms": (100.0, 5.0)}
    info = json.loads(path.read_text())
    assert info["params"] == P.to_dict() and info["carried_with"] == old
    cr.write_run_info(path, P)
    with pytest.raises(ValueError, match="cannot change"):
        cr.record_readmit(path, dataclasses.replace(P, pulse_exclude_post_ms=50.0))
    no_noise = {k: v for k, v in old.items() if k != "seed_noise_chunk_s"}
    assert cr._params_dict(no_noise)["seed_noise_chunk_s"] is None


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


def test_readmit_replays_the_carry_and_applies_new_pulse_parameters(tmp_path) -> None:
    path, rec, bank = _binary(tmp_path)
    params = dataclasses.replace(P, window_s=2.0, whiten_s=2.0, reseed_every=0, rebuild_min_accepted=20,
                                 rebuild_min_selected=10, rel_min_lead=10, sustained_min_spikes=5)
    W = cr.whitening_matrix(rec, params)
    groups = np.asarray(rec.get_property("group")).astype(int)
    on = np.arange(3_000, rec.get_num_frames() - 10_000, 9_000)
    pulses = (on, on + 3_000)
    out = tmp_path / "carry"
    cr.run_carry(path, W, bank.whitened(W, groups), out, params=params, windows=range(4), pulses=pulses,
                 bad_periods={0: [(20_000, 25_000)]}, n_jobs=1, workers=1)
    saved = {w: np.load(out / f"spikes_w{w:02d}.npz")["reason"] for w in range(4)}
    for w in range(4):
        z = np.load(out / f"spikes_w{w:02d}.npz")
        r, _ = cr.readmit_window(z, cr.Bank.load(out / f"bankin_w{w:02d}.npz"), params, fs=FS,
                                 a=w * 60_000, n_frames=60_000, pulses=pulses, bad_periods={0: [(20_000, 25_000)]})
        np.testing.assert_array_equal(r, saved[w])
    cr.readmit_carry(out, dataclasses.replace(params, pulse_relax=False), fs=FS, n_frames=rec.get_num_frames(),
                     pulses=pulses, bad_periods={0: [(20_000, 25_000)]})
    for w in range(4):
        z = np.load(out / f"spikes_w{w:02d}.npz")
        np.testing.assert_array_equal(z["reason_carry"], saved[w])
        assert z["frame"].size == saved[w].size


def test_channel_noise_is_the_mean_chunk_mad_clear_of_bad_periods() -> None:
    rng = np.random.default_rng(3)
    sd = np.array([1.0, 2.0, 3.0, 4.0, 1.5, 1.5, 1.5, 1.5])
    x = (rng.normal(size=(60_000, 8)) * sd).astype(np.float32)
    x[30_000:33_000, :4] += 200 * rng.normal(size=(3_000, 4)).astype(np.float32)
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    rec.set_property("group", np.repeat([0, 1], 4))
    clean = cr.channel_noise(rec, {0: [(30_000, 33_000)]}, n_jobs=2)
    np.testing.assert_allclose(clean, sd, rtol=0.03)
    loud = cr.channel_noise(rec, None, n_jobs=2)
    assert (loud[:4] > clean[:4] * 1.05).all()
    np.testing.assert_allclose(loud[4:], clean[4:])
    assert cr.window_bad_periods({0: [(10, 50), (90, 200)], 1: [(300, 400)]}, 40, 100) == {0: [(0, 10), (50, 60)]}


def test_seed_snr_divides_by_the_given_noise() -> None:
    T = _template()
    f = np.arange(1000, 290_000, 1500)
    x = _tetrode(frames=f, amps=np.full(f.size, 6.0), T=T)
    from probeinterface import Probe

    rec = si.NumpyRecording(x, sampling_frequency=FS)
    probe = Probe(ndim=2)
    probe.set_contacts(positions=np.array([[0, 0], [20, 0], [0, 20], [20, 20]], float), shapes="circle",
                       shape_params={"radius": 5})
    probe.set_device_channel_indices(np.arange(4))
    rec.set_probe(probe)
    rec.set_property("group", np.zeros(4, int))
    sorting = si.NumpySorting.from_samples_and_labels([f + P.nbefore], [np.zeros(f.size, int)], FS)
    sorting.set_property("group", np.array([0]))
    _, snr1 = cr.unit_templates(sorting, rec, n_jobs=1, noise=np.ones(4))
    _, snr2 = cr.unit_templates(sorting, rec, n_jobs=1, noise=np.full(4, 2.0))
    _, again = cr.unit_templates(sorting, rec, n_jobs=1, noise=np.ones(4))
    np.testing.assert_allclose(snr2, snr1 / 2, rtol=1e-5)
    np.testing.assert_array_equal(snr1, again)


def test_run_info_refuses_other_parameters(tmp_path) -> None:
    path = tmp_path / "carry_params.json"
    cr.write_run_info(path, P, binary="b")
    cr.write_run_info(path, P, binary="b")
    with pytest.raises(ValueError, match="different parameters"):
        cr.write_run_info(path, dataclasses.replace(P, blend=0.5))


@pytest.mark.parametrize("drop", [["5"], ["0", "1", "2", "3"]])
def test_carry_with_a_bad_wire_or_a_bad_tetrode_removed(tmp_path, drop) -> None:
    path, rec, bank = _binary(tmp_path)
    params = dataclasses.replace(P, window_s=2.0, whiten_s=2.0, reseed_every=0, rebuild_min_accepted=20,
                                 rebuild_min_selected=10, rel_min_lead=10)
    sub = cr.open_sorting_recording(path, drop)
    assert sub.get_num_channels() == 8 - len(drop)
    groups = np.asarray(sub.get_property("group")).astype(int)
    keep = [str(c) not in drop for c in rec.channel_ids]
    T = np.zeros_like(bank.T)
    for k, g in enumerate(bank.group):
        ch = np.flatnonzero(np.asarray(keep)[4 * g:4 * g + 4])
        T[k, :, :ch.size] = bank.T[k][:, ch]
    W = cr.whitening_matrix(sub, params)
    assert W.shape == (sub.get_num_channels(),) * 2
    for g in np.unique(groups):
        ch = np.flatnonzero(groups == g)
        blk = W[np.ix_(ch, ch)]
        assert blk.shape == (ch.size, ch.size)
    bank_w = cr.Bank(bank.ids, bank.group, T).whitened(W, groups)
    out = tmp_path / "carry"
    cr.run_carry(path, W, bank_w, out, params=params, windows=range(0, 2), drop_channel_ids=drop, n_jobs=1,
                 workers=1)
    z = np.load(out / "spikes_w00.npz")
    units = set(z["unit"][z["reason"] == 0].tolist())
    last = cr.Bank.load(out / "bank_w01.npz")
    if len(drop) == 4:
        assert units == {1} and last.ids.tolist() == [1]
    else:
        assert units == {0, 1}
        k = int(np.flatnonzero(last.ids == 1)[0])
        assert np.abs(last.T[k, :, 3]).max() == 0
        assert np.abs(last.T[k, :, :3]).max() > 0


def test_bad_channel_ids_resolve_channels_and_tetrodes() -> None:
    from ttm.sorting import bad_channel_ids

    rec = si.NumpyRecording(np.zeros((100, 8), np.float32), sampling_frequency=FS)
    rec.set_property("group", np.repeat([0, 1], 4))
    assert [str(c) for c in bad_channel_ids(rec, ["2"], ["TT2"])] == ["2", "4", "5", "6", "7"]
    with pytest.raises(ValueError, match="CH99"):
        bad_channel_ids(rec, ["CH99"])
    assert bad_channel_ids(rec, ["CH99"], missing="ignore") == []
