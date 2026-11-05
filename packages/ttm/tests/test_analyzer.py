import numpy as np
import pytest

pytest.importorskip("spikeinterface.sortingcomponents")
import spikeinterface as si
from ttm import analyzer as an

FS = 30000.0


def _sorting(trains: dict, groups: list[int]):
    ids = list(trains)
    fr = np.concatenate([np.asarray(trains[u], np.int64) for u in ids])
    lab = np.concatenate([np.full(len(trains[u]), u) for u in ids])
    o = np.argsort(fr, kind="stable")
    s = si.NumpySorting.from_samples_and_labels([fr[o]], [lab[o]], sampling_frequency=FS, unit_ids=np.asarray(ids))
    s.set_property("group", np.asarray(groups))
    return s


def test_unit_periods_leave_out_each_tetrodes_bad_intervals() -> None:
    s = _sorting({0: [10, 20], 7: [30]}, [0, 1])
    p = an.unit_periods(s, 1000, {1: [(500, 600), (100, 200)]})
    got = {(int(r["unit_index"]), int(r["start_sample_index"]), int(r["end_sample_index"])) for r in p}
    assert got == {(0, 0, 1000), (1, 0, 100), (1, 200, 500), (1, 600, 1000)}


def test_with_multi_unit_adds_one_unit_per_tetrode() -> None:
    s = _sorting({0: [10, 50], 1: [20], 2011: [30, 40]}, [0, 0, 1])
    out = an.with_multi_unit(s, {0: np.array([5, 70]), 1: np.array([60]), 2: np.zeros(0, np.int64)})
    assert list(out.unit_ids) == [0, 1, 2011, 10001, 10002]
    assert list(out.get_property("group")) == [0, 0, 1, 0, 1]
    assert list(out.get_property("is_mua")) == [False, False, False, True, True]
    assert list(out.get_unit_spike_train(10001)) == [5, 70]
    assert list(out.get_unit_spike_train(2011)) == [30, 40]
    with pytest.raises(ValueError):
        an.with_multi_unit(s, {0: np.array([5])}, unit_id_base=2010)


def _clusters(n_own=400, n_other=900, seed=0):
    rng = np.random.default_rng(seed)
    own = rng.normal(0, 1, (n_own, 6)) + 2.5
    a = rng.normal(0, 1.3, (n_other // 2, 6))
    b = rng.normal(0, 0.8, (n_other - n_other // 2, 6)) - 2.0
    X = np.vstack([own, a, b])
    lab = np.r_[np.zeros(n_own, int), np.ones(n_other // 2, int), np.full(n_other - n_other // 2, 2)]
    return X, lab


def test_pc_metrics_match_spikeinterface_definitions() -> None:
    from spikeinterface.metrics.quality.pca_metrics import (
        d_prime_metric,
        mahalanobis_metrics,
        simplified_silhouette_score,
    )

    X, lab = _clusters()
    iso, lr = an.mahalanobis_all(X[lab == 0], X[lab != 0])
    iso_si, lr_si = mahalanobis_metrics(X, lab, 0)
    assert iso == pytest.approx(iso_si, rel=1e-9)
    assert lr == pytest.approx(lr_si, rel=1e-6)
    assert an.d_prime_all(X, lab == 0) == pytest.approx(abs(d_prime_metric(X, lab, 0)), rel=1e-6)
    assert an.d_prime_all(X, lab == 1) == pytest.approx(abs(d_prime_metric(X, lab, 1)), rel=1e-6)
    assert an.silhouette_all(X, lab, 0) == pytest.approx(simplified_silhouette_score(X, lab, 0), rel=1e-9)


def test_isolation_distance_is_missing_when_the_unit_outnumbers_the_rest() -> None:
    X, lab = _clusters(n_own=1000, n_other=900)
    iso, lr = an.mahalanobis_all(X[lab == 0], X[lab != 0])
    assert np.isnan(iso)
    assert np.isfinite(lr)


def _save_recording(path, x, groups):
    from probeinterface import Probe

    rec = si.NumpyRecording([x], sampling_frequency=FS)
    pos = np.array([[g * 300.0 + dx, dy] for g in groups for dx, dy in [(0, 0), (20, 0), (0, 20), (20, 20)]])
    pos = pos[: x.shape[1]]
    probe = Probe(ndim=2)
    probe.set_contacts(positions=pos, shapes="circle", shape_params={"radius": 5})
    probe.set_device_channel_indices(np.arange(x.shape[1]))
    rec.set_probe(probe)
    rec.set_property("group", np.repeat(np.unique(groups), 4)[: x.shape[1]])
    rec.set_channel_gains(0.195)
    rec.set_channel_offsets(0.0)
    rec.save(folder=path, format="binary", dtype="int16", n_jobs=1, progress_bar=False)
    return path


def test_chunk_noise_levels_leave_out_bad_intervals(tmp_path) -> None:
    rng = np.random.default_rng(0)
    n = int(30 * FS)
    x = rng.normal(0, 100, (n, 8))
    x[int(10 * FS):int(20 * FS), 4:] *= 5
    path = _save_recording(tmp_path / "rec", x.astype(np.int16), [0, 1])
    params = an.AnalyzerParams(noise_chunks=20)
    noise, table = an.chunk_noise_levels(path, {1: [(int(10 * FS), int(20 * FS))]}, params, workers=2)
    assert list(table["chunks"]) == [20] * 8
    assert np.all(np.abs(noise / 100 - 1) < 0.03)
    assert np.allclose(table["noise_uv"], noise * 0.195)


def test_multi_unit_admits_missed_spikes_and_not_claimed_ones(tmp_path) -> None:
    rng = np.random.default_rng(1)
    n = int(40 * FS)
    x = rng.normal(0, 10, (n, 4))
    t = np.arange(90) - 30
    shape = -np.exp(-(t / 4.0) ** 2) + 0.4 * np.exp(-((t - 12) / 8.0) ** 2)
    T = shape[:, None] * np.array([150.0, 110.0, 80.0, 60.0])
    spikes = np.sort(rng.choice(np.arange(200, n - 200, 400), 2000, replace=False))
    for f in spikes:
        x[f - 30:f + 60] += T
    path = _save_recording(tmp_path / "rec", x.astype(np.int16), [0])
    missed = spikes[::10]
    sorting = _sorting({0: np.setdiff1d(spikes, missed)}, [0])
    sorting.save(folder=tmp_path / "sorting")
    params = an.AnalyzerParams(mua_window_s=10.0, mua_template_min_spikes=20, mua_min_reference=100)
    noise = np.full(4, 10.0)
    an.run_multi_unit(path, tmp_path / "sorting", noise, {}, tmp_path / "mua", params, workers=2)
    events, table = an.multi_unit_events(tmp_path / "mua", params)
    got = events[0]
    near = np.abs(got[:, None] - missed[None, :]).min(axis=1) <= 2
    assert near.all()
    assert got.size >= 0.9 * missed.size
    assert table.loc[0, "reference"] >= 0.9 * (spikes.size - missed.size)
