"""Single-sample Open Ephys timestamp glitches are repaired; real gaps are not."""

import numpy as np

from ttm.spikeinterface import repair_timestamp_glitches

FS = 30000.0


def test_glitch_is_replaced_by_neighbour_mean():
    t = np.arange(1000) / FS
    t[400] += 85899.37926
    fixed = repair_timestamp_glitches(t, fs=FS)
    assert fixed == [400]
    assert np.allclose(np.diff(t), 1 / FS)


def test_real_gap_is_left_alone():
    t = np.arange(1000) / FS
    t[500:] += 2.0
    assert repair_timestamp_glitches(t, fs=FS) == []
    assert t[500] - t[499] > 1.0


def test_clean_vector_untouched():
    t = np.arange(1000) / FS
    before = t.copy()
    assert repair_timestamp_glitches(t, fs=FS) == []
    assert np.array_equal(t, before)
