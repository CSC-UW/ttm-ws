import numpy as np
import pytest

pytest.importorskip("spikeinterface.preprocessing")
import spikeinterface as si
import spikeinterface.preprocessing as spre
from ttm.sorting import preprocess_for_sorting

FS = 30000.0


def _grouped_recording(n_tetrodes=3, seconds=0.5, seed=0):
    rng = np.random.default_rng(seed)
    n_ch = 4 * n_tetrodes
    traces = rng.normal(size=(int(FS * seconds), n_ch)).astype("float32")
    traces += rng.normal(size=(traces.shape[0], 1)).astype("float32") * 3  # shared noise
    rec = si.NumpyRecording(traces, sampling_frequency=FS)
    rec.set_property("group", np.repeat(np.arange(n_tetrodes), 4))
    return rec


def test_cross_tetrode_cmr_is_median_of_other_tetrodes() -> None:
    rec = _grouped_recording()
    groups = rec.get_property("group")
    bp = spre.bandpass_filter(rec, freq_min=300.0, freq_max=6000.0, dtype="float32").get_traces()
    got = preprocess_for_sorting(rec, cmr="cross_tetrode").get_traces()
    expected = np.empty_like(bp)
    for g in np.unique(groups):
        ref = np.median(bp[:, groups != g], axis=1, keepdims=True)
        expected[:, groups == g] = bp[:, groups == g] - ref
    np.testing.assert_allclose(got, expected, rtol=0, atol=1e-5)
    assert np.array_equal(preprocess_for_sorting(rec, cmr="cross_tetrode").get_property("group"), groups)


def test_cross_tetrode_differs_from_global_and_unknown_cmr_raises() -> None:
    rec = _grouped_recording()
    cross = preprocess_for_sorting(rec, cmr="cross_tetrode").get_traces()
    glob = preprocess_for_sorting(rec, cmr="global").get_traces()
    assert not np.allclose(cross, glob)
    with pytest.raises(ValueError, match="Unknown cmr"):
        preprocess_for_sorting(rec, cmr="local")
