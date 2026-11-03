import numpy as np
import pytest

pytest.importorskip("spikeinterface.preprocessing")
import spikeinterface as si
import spikeinterface.preprocessing as spre
from spikeinterface.preprocessing.phase_shift import apply_frequency_shift
from ttm.reference import (
    AlignedMedianReferenceRecording,
    estimate_wire_timing,
    estimate_wire_timing_windows,
)
from ttm.sorting import preprocess_for_sorting

FS = 30000.0


def _common_signal(n, seed=0):
    rng = np.random.default_rng(seed)
    s = np.zeros(n)
    t = np.arange(-15, 16)
    pulse = -np.exp(-(t / 3.0) ** 2)
    for i in rng.integers(100, n - 100, size=n // 500):
        s[i - 15:i + 16] += 200 * pulse * rng.uniform(0.5, 1.5)
    return s


def _skewed(n=60000, n_ch=12, seed=0):
    rng = np.random.default_rng(seed)
    s = _common_signal(n, seed)
    tau = rng.uniform(-0.6, 0.9, size=n_ch)
    g = rng.uniform(0.9, 1.1, size=n_ch)
    delayed = apply_frequency_shift(np.repeat(s[:, None], n_ch, axis=1), tau, axis=0)
    x = delayed * g[None, :] + rng.normal(scale=2.0, size=(n, n_ch))
    return x.astype("float32"), tau, g


def test_estimate_wire_timing_recovers_relative_offsets() -> None:
    x, tau, _ = _skewed()
    est, gain, n_event = estimate_wire_timing(x)
    assert n_event >= 3000
    np.testing.assert_allclose(est - est.mean(), tau - tau.mean(), atol=0.03)
    assert np.all((gain > 0.8) & (gain < 1.2))


def test_aligned_reference_removes_more_of_the_common_signal_than_the_median() -> None:
    x, _, _ = _skewed()
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    timing = estimate_wire_timing_windows(rec, window_frames=x.shape[0], block_frames=x.shape[0], n_blocks=1)
    aligned = AlignedMedianReferenceRecording(rec, timing["starts"], timing["shifts"], timing["gains"])
    y = aligned.get_traces()[5000:-5000]
    plain = (x - np.median(x, axis=1, keepdims=True))[5000:-5000]
    assert np.std(y) < 0.5 * np.std(plain)
    assert np.std(y) < 3.0


def test_window_boundary_and_chunking_do_not_change_the_traces() -> None:
    x, _, _ = _skewed(n=40000)
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    rec = si.NumpyRecording(spre.bandpass_filter(rec, freq_min=300.0, freq_max=6000.0).get_traces(), sampling_frequency=FS)
    rng = np.random.default_rng(1)
    starts = np.array([0, 20000])
    shifts = rng.uniform(-0.5, 0.5, size=(2, x.shape[1]))
    gains = rng.uniform(0.9, 1.1, size=(2, x.shape[1]))
    aligned = AlignedMedianReferenceRecording(rec, starts, shifts, gains)
    whole = aligned.get_traces()
    pieces = np.concatenate([aligned.get_traces(start_frame=a, end_frame=b)
                             for a, b in ((0, 7000), (7000, 19990), (19990, 20010), (20010, 40000))])
    np.testing.assert_allclose(pieces, whole, atol=1e-3)
    one = AlignedMedianReferenceRecording(rec, starts[:1], shifts[1:], gains[1:]).get_traces()
    np.testing.assert_allclose(whole[20000:], one[20000:], atol=1e-3)
    sub = aligned.get_traces(start_frame=100, end_frame=200, channel_ids=aligned.channel_ids[2:4])
    np.testing.assert_allclose(sub, whole[100:200, 2:4], atol=1e-3)


def test_global_aligned_needs_wire_timing() -> None:
    x, _, _ = _skewed(n=6000, n_ch=4)
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    with pytest.raises(ValueError, match="wire_timing"):
        preprocess_for_sorting(rec, cmr="global_aligned")


def test_windows_without_events_take_the_nearest_estimate() -> None:
    from ttm.reference import fill_missing_windows

    shifts = np.array([[np.nan, np.nan], [0.1, 0.2], [np.nan, np.nan], [np.nan, np.nan], [0.5, 0.6]])
    out = fill_missing_windows({"starts": np.arange(5) * 10, "shifts": shifts, "gains": np.where(np.isnan(shifts), np.nan, 1.0)})
    np.testing.assert_allclose(out["shifts"][:, 0], [0.1, 0.1, 0.1, 0.5, 0.5])
    assert out["filled"].tolist() == [True, False, True, True, False]


def test_pooling_over_windows_outvotes_a_single_outlier_block() -> None:
    from ttm.reference import pool_blocks

    good = (np.full((3, 2), -0.3), np.ones((3, 2)))
    bad = (np.full((1, 2), -3.0), np.full((1, 2), 0.3))
    empty = (np.zeros((0, 2)), np.zeros((0, 2)))
    shift, gain = pool_blocks([bad])
    np.testing.assert_allclose(shift, 3.0)
    shift, gain = pool_blocks([good, bad, empty, good])
    np.testing.assert_allclose(shift, 0.3)
    np.testing.assert_allclose(gain, 1.0)
    assert np.isnan(pool_blocks([bad, empty], min_blocks=3)[0]).all()


def test_no_window_with_events_falls_back_to_the_global_median_when_asked() -> None:
    from ttm.reference import fill_missing_windows

    timing = {"starts": np.arange(3) * 10, "shifts": np.full((3, 2), np.nan), "gains": np.full((3, 2), np.nan)}
    with pytest.raises(ValueError):
        fill_missing_windows(timing)
    out = fill_missing_windows(timing, if_none="identity")
    np.testing.assert_allclose(out["shifts"], 0.0)
    np.testing.assert_allclose(out["gains"], 1.0)
    assert out["filled"].all()


def test_pooled_windows_match_the_single_window_estimate_on_steady_offsets() -> None:
    x, tau, _ = _skewed(n=135000)
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    own = estimate_wire_timing_windows(rec, window_frames=45000, block_frames=45000, n_blocks=1)
    pooled = estimate_wire_timing_windows(rec, window_frames=45000, block_frames=45000, n_blocks=1,
                                          pool_windows=1)
    assert not pooled["filled"].any()
    for t in (own, pooled):
        est = -t["shifts"]
        np.testing.assert_allclose(est - est.mean(axis=1, keepdims=True), np.tile(tau - tau.mean(), (3, 1)),
                                   atol=0.05)
