import itertools

import numpy as np
import pandas as pd
import pytest
import spikeinterface as si

from ttm.bad_periods import (bad_mask, bin_features, bin_stats, detect_bad_periods, group_intervals,
                             mask_bad_periods, valid_seconds)

FS = 1000.0


def _recording(n_s=600, seed=0):
    rng = np.random.default_rng(seed)
    n = int(n_s * FS)
    shared = rng.standard_normal((n, 1)).astype(np.float32)
    x = (0.6 * shared + 0.8 * rng.standard_normal((n, 8))).astype(np.float32) * 10
    rec = si.NumpyRecording([x], sampling_frequency=FS)
    rec.set_property("group", np.repeat([0, 1], 4))
    return rec, x


def test_detects_one_loud_wire_and_a_shared_burst():
    rec, x = _recording()
    x = x.copy()
    x[int(200 * FS):int(230 * FS), 1] *= 6                     # one loud wire on group 0
    burst = np.random.default_rng(1).standard_normal(int(20 * FS)).astype(np.float32) * 20
    x[int(400 * FS):int(420 * FS), 4:8] += burst[:, None]      # a shared signal on group 1
    rec = si.NumpyRecording([x], sampling_frequency=FS)
    rec.set_property("group", np.repeat([0, 1], 4))
    st = bin_stats(rec, bin_s=0.1, chunk_s=10.0, n_jobs=1)
    f = bin_features(st, bin_s=1.0)
    df = detect_bad_periods(f, max_gap_s=2.0, pad_s=1.0)
    assert len(df) == 2
    a, b = df.iloc[0], df.iloc[1]
    assert a.group == 0 and a.kind == "loud" and a.channels == "1"
    assert a.start_s == pytest.approx(199.0) and a.end_s == pytest.approx(231.0)
    assert b.group == 1 and b.kind == "common" and b.channels == "4,5,6,7"
    assert b.start_s == pytest.approx(399.0) and b.end_s == pytest.approx(421.0)


def test_bin_stats_aggregates_exactly():
    rec, x = _recording(n_s=60)
    st = bin_stats(rec, bin_s=0.1, chunk_s=10.0, n_jobs=1)
    f = bin_features(st, bin_s=2.0)
    seg = x[int(10 * FS):int(12 * FS)].astype(np.float64)
    assert f["rms"][5] == pytest.approx(np.sqrt((seg ** 2).mean(axis=0)), rel=1e-5)
    c = np.corrcoef(seg[:, :4].T)
    num = (seg[:, 0] * seg[:, 1]).sum() / np.sqrt((seg[:, 0] ** 2).sum() * (seg[:, 1] ** 2).sum())
    assert f["corr"][5, 0] == pytest.approx(num, rel=1e-5)
    assert abs(num - c[0, 1]) < 0.05


def _table():
    return pd.DataFrame({"group": [1, 0, 1], "start_s": [10.0, 20.0, 12.0], "end_s": [15.0, 21.5, 18.0],
                         "kind": "loud", "channels": "", "max_loud": 5.0, "max_corr_excess": 0.0,
                         "flagged_s": 1.0})


def test_table_helpers():
    t = _table()
    iv = group_intervals(t, FS, int(100 * FS))
    assert iv[1] == [(10000, 18000)] and iv[0] == [(20000, 21500)]
    assert valid_seconds(t, 100.0, [0, 1, 2]) == {0: 98.5, 1: 92.0, 2: 100.0}
    m = bad_mask(t, 1, np.array([9999, 10000, 17999, 18000]), FS)
    assert m.tolist() == [False, True, True, False]


@pytest.mark.parametrize("mode", ["zeros", "noise"])
def test_mask_replaces_only_the_flagged_group(mode):
    rec, x = _recording(n_s=60)
    m = mask_bad_periods(rec, _table(), mode=mode, taper_frames=10, noise_levels=np.full(8, 3.0))
    y = m.get_traces(start_frame=0, end_frame=int(60 * FS))
    inside = slice(int(10 * FS), int(18 * FS))
    np.testing.assert_array_equal(y[inside, :4], x[inside, :4])
    if mode == "zeros":
        assert not y[inside, 4:].any()
    else:
        assert np.std(y[inside, 4:]) == pytest.approx(3.0, rel=0.05)
        assert abs(np.corrcoef(y[inside, 4], x[inside, 4])[0, 1]) < 0.05
    far = slice(0, int(9 * FS))
    np.testing.assert_array_equal(y[far], x[far])
    edge = y[int(18 * FS):int(18 * FS) + 10, 4] - x[int(18 * FS):int(18 * FS) + 10, 4]
    assert np.any(edge != 0) and np.all(y[int(18 * FS) + 10:int(19 * FS)] == x[int(18 * FS) + 10:int(19 * FS)])


def test_mask_is_chunk_independent_and_serializable(tmp_path):
    rec, _ = _recording(n_s=60)
    m = mask_bad_periods(rec, _table(), mode="noise", noise_levels=np.full(8, 3.0), seed=7)
    whole = m.get_traces(start_frame=int(9 * FS), end_frame=int(19 * FS))
    parts = np.vstack([m.get_traces(start_frame=a, end_frame=b)
                       for a, b in [(9000, 12345), (12345, 12346), (12346, 19000)]])
    np.testing.assert_array_equal(whole, parts)
    sub = m.get_traces(start_frame=int(9 * FS), end_frame=int(19 * FS), channel_ids=m.channel_ids[[5, 2]])
    np.testing.assert_array_equal(sub, whole[:, [5, 2]])
    m2 = si.load(m.to_dict())
    np.testing.assert_array_equal(m2.get_traces(start_frame=9000, end_frame=19000), whole)


def test_mask_with_overlapping_and_nearby_periods():
    rec, x = _recording(n_s=40)
    t = pd.DataFrame({"group": [0, 1, 0], "start_s": [10.0, 12.0, 16.01], "end_s": [16.0, 14.0, 17.0],
                      "kind": "loud", "channels": "", "max_loud": 5.0, "max_corr_excess": 0.0, "flagged_s": 1.0})
    m = mask_bad_periods(rec, t, mode="noise", taper_frames=30, noise_levels=np.full(8, 3.0), seed=3)
    whole = m.get_traces(start_frame=0, end_frame=int(30 * FS))
    cuts = [0, 9990, 12005, 13999, 16005, 16015, 30000]
    parts = np.vstack([m.get_traces(start_frame=a, end_frame=b) for a, b in itertools.pairwise(cuts)])
    np.testing.assert_array_equal(whole, parts)
    np.testing.assert_array_equal(whole[:9970], x[:9970])
    np.testing.assert_array_equal(whole[17030:], x[17030:30000])
    np.testing.assert_array_equal(whole[:11970, 4:], x[:11970, 4:])
    np.testing.assert_array_equal(whole[14030:, 4:], x[14030:30000, 4:])
    assert np.std(whole[10000:16000, :4]) == pytest.approx(3.0, rel=0.05)
    assert np.std(whole[12000:14000, 4:]) == pytest.approx(3.0, rel=0.05)
    assert np.all(whole[16000:16010, :4] != x[16000:16010, :4])


def test_cached_scan_is_reused_only_for_the_same_input(tmp_path):
    import json
    import types

    from ttm.cli.bad_periods import cached_key, scan_key

    exp = types.SimpleNamespace(store_path=lambda: tmp_path / "s.zarr")
    key = scan_key(exp, None, "global_aligned", ["CH22"])
    np.savez(tmp_path / "new.npz", ss=np.zeros(3), scan_key=json.dumps(key))
    np.savez(tmp_path / "old.npz", ss=np.zeros(3))
    assert cached_key(dict(np.load(tmp_path / "new.npz"))) == key
    assert cached_key(dict(np.load(tmp_path / "old.npz"))) is None
    assert scan_key(exp, None, "global", ["CH22"]) != key
    assert scan_key(exp, None, "global_aligned", []) != key
    assert scan_key(exp, tmp_path / "binary", "global_aligned", ["CH22"])["cmr"] == "binary"
