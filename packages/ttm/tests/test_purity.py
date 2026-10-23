"""Tests for ``ttm.purity``: per-spike scores, the purity table, templates and unit pairing.

A 20 s two-tetrode binary with three planted cells (two of different shape on
tetrode 0, one on tetrode 1) checks that every planted spike is scored with
its trough depth in noise units, the cosine a template plus white noise
should have and a residual energy near 1;
that a unit made by pooling the two tetrode-0 cells has a larger residual
energy than either cell; that the cache round-trips; and that pairing a
sorting with itself agrees completely.
"""

from __future__ import annotations

import numpy as np
import pytest
import spikeinterface as si
from ttm import capture as cp
from ttm import purity as pu

FS = 30000.0
SIGMA = 10.0


def _template(amp, profile, width=3.0, lag=12):
    t = np.arange(cp.N_SAMP) - cp.NBEFORE
    w = -amp * np.exp(-(t ** 2) / (2 * width ** 2)) + 0.3 * amp * np.exp(-((t - lag) ** 2) / 50.0)
    return (w[:, None] * np.asarray(profile)[None, :]).astype(np.float32)


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    rng = np.random.default_rng(0)
    n = int(20 * FS)
    x = rng.normal(0, SIGMA, size=(n, 8)).astype(np.float32)
    shapes = {0: (_template(200.0, [1.0, 0.5, 0.2, 0.1]), 0),
              1: (_template(150.0, [0.1, 0.2, 0.6, 1.0], width=2.0, lag=8), 0),
              2: (_template(180.0, [0.6, 1.0, 0.4, 0.2]), 1)}
    trains = {}
    for u, (t, g) in shapes.items():
        fr = np.sort(rng.choice(np.arange(300 + 37 * u, n - 300, 150), 250, replace=False))
        trains[u] = fr.astype(np.int64)
        for f in fr:
            x[f - cp.NBEFORE:f + cp.NAFTER, 4 * g:4 * g + 4] += t
    rec = si.NumpyRecording([x], sampling_frequency=FS)
    rec.set_property("group", np.repeat([0, 1], 4))
    rec.set_channel_gains(0.195)
    rec.set_channel_offsets(0.0)
    folder = tmp_path_factory.mktemp("purity") / "binary"
    rec.save(folder=folder, format="binary", n_jobs=1, progress_bar=False)
    sorting = si.NumpySorting.from_unit_dict([trains], sampling_frequency=FS)
    sorting.set_property("group", np.array([0, 0, 1]))
    pooled = si.NumpySorting.from_unit_dict([{0: np.sort(np.concatenate([trains[0], trains[1]])), 2: trains[2]}],
                                            sampling_frequency=FS)
    pooled.set_property("group", np.array([0, 1]))
    return str(folder), sorting, pooled, trains, shapes, n


def test_scores_of_planted_cells(planted, tmp_path):
    binary, sorting, _, trains, shapes, _ = planted
    cache = tmp_path / "scores.npz"
    frames, amps, rs, resid, ugroup = pu.score_sorting(binary, sorting, workers=1, cache=cache)
    assert ugroup == {0: 0, 1: 0, 2: 1}
    for u, (t, _) in shapes.items():
        assert np.array_equal(frames[u], trains[u])
        expected = -t.min() / SIGMA
        assert np.median(amps[u]) == pytest.approx(expected, rel=0.1)
        energy = float((t * t).sum())
        assert np.median(rs[u]) == pytest.approx(np.sqrt(energy / (energy + t.size * SIGMA ** 2)), abs=0.03)
        assert 0.8 < np.median(resid[u]) < 1.3
    again = pu.score_sorting(binary, sorting, workers=1, cache=cache)
    assert all(np.array_equal(again[3][u], resid[u]) for u in resid)


def test_cache_without_residuals_is_recomputed(planted, tmp_path):
    binary, sorting, _, _, _, _ = planted
    cache = tmp_path / "old.npz"
    np.savez(cache, **{f"frames_{u}": np.zeros(1, np.int64) for u in (0, 1, 2)},
             **{f"amps_{u}": np.zeros(1) for u in (0, 1, 2)}, **{f"rs_{u}": np.zeros(1) for u in (0, 1, 2)})
    frames, _, _, resid, _ = pu.score_sorting(binary, sorting, workers=1, cache=cache)
    assert frames[0].size == 250 and resid[0].size == 250
    assert "resid_0" in np.load(cache)


def test_frame_range_outside_reads_nothing(planted):
    binary, sorting, _, _, _, n = planted
    frames, _, _, _, _ = pu.score_sorting(binary, sorting, workers=1, frame_range=(10 * n, 11 * n))
    assert all(f.size == 0 for f in frames.values())


def test_pooled_cells_have_larger_spread(planted):
    binary, sorting, pooled, _, _, n = planted
    f1, a1, r1, e1, g1 = pu.score_sorting(binary, sorting, workers=1)
    _, _, _, e2, _ = pu.score_sorting(binary, pooled, workers=1)
    assert np.median(e2[0]) > 2 * max(np.median(e1[0]), np.median(e1[1]))
    tab = pu.purity_table(f1, a1, g1, n, 0.0, rs=r1, resid=e1)
    assert list(tab.unit) == [0, 1, 2]
    for col in ("rp_contamination", "ampsplit_ratio", "n_distinct", "cos_median", "cos_p10", "resid_median", "amp_cv"):
        assert col in tab
    assert (tab.rp_contamination.fillna(0) < 0.05).all()
    assert tab.set_index("unit").loc[0, "neighbours"].startswith("1:")


def test_templates_and_self_pairing(planted):
    binary, _, _, trains, shapes, n = planted
    T, _ = pu.unit_templates(binary, trains, {0: 0, 1: 0, 2: 1}, n, n_windows=4, win_s=2.0, workers=1)
    for u, (t, _) in shapes.items():
        assert T[u].shape == (cp.N_SAMP, 4)
        assert int(np.argmin(T[u].min(axis=1))) == cp.NBEFORE
        assert cp.template_cosine(T[u], t) > 0.98
    m = pu.match_units(("a", "b"), (trains, trains), ({0: 0, 1: 0, 2: 1},) * 2, (T, T), ([0, 1, 2], [0, 1, 2]))
    assert list(m.a) == [0, 1, 2] and list(m.b) == [0, 1, 2]
    assert (m.agreement == 1.0).all()
    assert pu.frac_within(np.array([100, 200]), np.array([110, 400])) == 0.5
