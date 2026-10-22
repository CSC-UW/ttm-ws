"""Group-sparse analyzers keep each unit's main channel on its own tetrode."""

import numpy as np
import spikeinterface as si

from ttm.sorting import create_group_sparse_analyzer, group_main_channel_indices

FS = 30000.0


def _recording_and_sorting():
    """Two tetrodes; unit 0 lives on tetrode 0 but leaks a larger spike onto tetrode 1's first wire."""
    rng = np.random.default_rng(0)
    n = int(20 * FS)
    x = rng.normal(0, 5, size=(n, 8)).astype("float32")
    t0 = np.arange(int(0.1 * FS), n - 100, int(0.05 * FS))
    t1 = t0 + int(0.02 * FS)
    for t in t0:
        x[t, 0] -= 100; x[t, 1] -= 60; x[t, 4] -= 150      # unit 0: extremum off its tetrode
    for t in t1:
        x[t, 5] -= 120; x[t, 6] -= 80                       # unit 1 on tetrode 1
    rec = si.NumpyRecording(x, sampling_frequency=FS)
    locations = np.column_stack([np.arange(8) * 10.0 + (np.arange(8) >= 4) * 300.0, np.zeros(8)])
    rec = rec.set_dummy_probe_from_locations(locations) or rec
    rec.set_property("group", np.array([0, 0, 0, 0, 1, 1, 1, 1]))
    sorting = si.NumpySorting.from_unit_dict([{0: t0, 1: t1}], sampling_frequency=FS)
    sorting.set_property("group", np.array([0, 1]))
    return rec, sorting


def test_main_channel_stays_within_group():
    rec, sorting = _recording_and_sorting()
    mask = np.array([[1, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 1, 1, 1]], dtype=bool)
    main = group_main_channel_indices(rec, sorting, mask, seed=0)
    assert main.tolist() == [0, 5]


def test_create_group_sparse_analyzer_accepts_leaking_unit():
    rec, sorting = _recording_and_sorting()
    az = create_group_sparse_analyzer(sorting, rec, format="memory", return_in_uV=False)
    assert az.sparsity.mask.sum(axis=1).tolist() == [4, 4]
    assert az.get_main_channels(outputs="index", with_dict=False).tolist() == [0, 5]
