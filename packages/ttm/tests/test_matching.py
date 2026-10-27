"""Unit tests for ``ttm.matching`` (template matching per channel group).

Synthetic recording: three tetrodes of four channels, two units each, Gaussian noise. The per-group
matcher must give exactly what separate matchers on each tetrode's channels give, put spikes on the
right units, and leave a tetrode's spikes unchanged when another tetrode carries an artifact.
"""
import numpy as np
import pytest
import spikeinterface as si
from spikeinterface.core import Templates
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates
from ttm.matching import (
    PerGroupMatching,
    find_spikes_per_group,
    group_templates,
    template_groups,
)

FS = 30000.0
T, NBEFORE = 90, 30
N_TET = 3
AMP = {"amplitudes": [0.8, np.inf]}


def _wf(peak_ch, amp):
    t = np.arange(T) - NBEFORE
    base = amp * np.exp(-(t**2) / 8.0) - 0.4 * amp * np.exp(-((t - 6) ** 2) / 20.0)
    out = np.zeros((T, 4), dtype=np.float32)
    out[:, peak_ch] = base
    out[:, (peak_ch + 1) % 4] = 0.4 * base
    return out


def _bank():
    shapes = [_wf(u % 4, -120.0 - 30 * (u % 2)) for u in range(2 * N_TET)]
    dense = np.zeros((2 * N_TET, T, 4 * N_TET), dtype=np.float32)
    mask = np.zeros((2 * N_TET, 4 * N_TET), dtype=bool)
    for u, w in enumerate(shapes):
        g = u // 2
        dense[u, :, 4 * g: 4 * g + 4] = w
        mask[u, 4 * g: 4 * g + 4] = True
    packed = np.stack([dense[u][:, mask[u]] for u in range(len(mask))])
    return Templates(templates_array=packed, sampling_frequency=FS, nbefore=NBEFORE, is_in_uV=False,
                     sparsity_mask=mask, channel_ids=np.array([f"CH{c}" for c in range(4 * N_TET)]),
                     unit_ids=np.arange(10, 10 + 2 * N_TET)), dense


def _recording(dense, *, seconds=6.0, artifact_on=None, seed=0):
    rng = np.random.default_rng(seed)
    n = int(seconds * FS)
    x = rng.normal(0, 10.0, (n, 4 * N_TET)).astype(np.float32)
    truth = {}
    for u in range(len(dense)):
        times = np.sort(rng.choice(np.arange(200, n - 200, 120), size=int(15 * seconds), replace=False))
        truth[u] = times
        for t0 in times:
            x[t0 - NBEFORE: t0 - NBEFORE + T] += dense[u]
    if artifact_on is not None:
        a = np.arange(n // 3, 2 * n // 3)
        burst = 400.0 * np.sin(2 * np.pi * 900.0 * a / FS) * (rng.random(len(a)) < 0.6)
        for c in range(4 * artifact_on, 4 * artifact_on + 4):
            x[a, c] += burst.astype(np.float32)
    rec = si.NumpyRecording([x], FS, channel_ids=np.array([f"CH{c}" for c in range(4 * N_TET)]))
    rec.set_property("group", np.repeat(np.arange(N_TET), 4))
    return rec, truth


def _per_group_by_hand(rec, templates):
    groups = rec.get_property("group")
    ug = template_groups(templates, groups)
    out = []
    for g in range(N_TET):
        ch = np.flatnonzero(groups == g)
        units = np.flatnonzero(ug == g)
        s = find_spikes_from_templates(rec.select_channels(rec.channel_ids[ch]), group_templates(templates, units, ch),
                                       method="circus-omp", method_kwargs=AMP,
                                       job_kwargs={"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False})
        s["cluster_index"] = units[s["cluster_index"]]
        out.append(s)
    return np.concatenate(out)


def _key(spikes):
    return sorted(zip(spikes["sample_index"].tolist(), spikes["cluster_index"].tolist()))


JOB = {"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False}


def test_template_groups_and_restriction():
    templates, dense = _bank()
    groups = np.repeat(np.arange(N_TET), 4)
    assert template_groups(templates, groups).tolist() == [0, 0, 1, 1, 2, 2]
    sub = group_templates(templates, np.array([2, 3]), np.arange(4, 8))
    assert sub.unit_ids.tolist() == [12, 13]
    np.testing.assert_array_equal(sub.get_dense_templates(), dense[2:4, :, 4:8])


def test_template_across_groups_is_refused():
    templates, dense = _bank()
    mask = templates.sparsity_mask.copy()
    mask[0, 4] = True
    bad = Templates(templates_array=np.stack([dense[u][:, mask[u]] if mask[u].sum() == 5 else
                                              np.pad(dense[u][:, mask[u]], ((0, 0), (0, 1)))
                                              for u in range(len(mask))]),
                    sampling_frequency=FS, nbefore=NBEFORE, is_in_uV=False, sparsity_mask=mask,
                    channel_ids=templates.channel_ids, unit_ids=templates.unit_ids,
                    check_for_consistent_sparsity=False)
    with pytest.raises(ValueError, match="spans groups"):
        template_groups(bad, np.repeat(np.arange(N_TET), 4))


def test_equals_separate_matchers_per_tetrode():
    templates, dense = _bank()
    rec, _ = _recording(dense)
    got = find_spikes_per_group(rec, templates, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    assert _key(got) == _key(_per_group_by_hand(rec, templates))


def test_spikes_land_on_the_right_units():
    templates, dense = _bank()
    rec, truth = _recording(dense)
    got = find_spikes_per_group(rec, templates, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    for u, times in truth.items():
        mine = got["sample_index"][got["cluster_index"] == u]
        hit = np.isin(times, mine) | np.isin(times + 1, mine) | np.isin(times - 1, mine)
        assert hit.mean() > 0.95, (u, hit.mean())
    channel_group = rec.get_property("group")[got["channel_index"]]
    assert np.array_equal(channel_group, got["cluster_index"] // 2)


def test_artifact_on_one_tetrode_leaves_the_others_unchanged():
    templates, dense = _bank()
    clean, _ = _recording(dense, seed=1)
    noisy, _ = _recording(dense, seed=1, artifact_on=0)
    a = find_spikes_per_group(clean, templates, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    b = find_spikes_per_group(noisy, templates, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    others = lambda s: s[s["cluster_index"] >= 2]
    assert _key(others(a)) == _key(others(b))
    assert _key(a) != _key(b)


def test_parallel_equals_serial(tmp_path):
    templates, dense = _bank()
    rec, _ = _recording(dense)
    rec = rec.save(folder=tmp_path / "rec")
    one = find_spikes_per_group(rec, templates, method="circus-omp", method_kwargs=AMP, job_kwargs=JOB)
    two = find_spikes_per_group(rec, templates, method="circus-omp", method_kwargs=AMP,
                                job_kwargs={"n_jobs": 2, "chunk_duration": "1s", "progress_bar": False})
    assert _key(one) == _key(two)


def test_node_needs_matching_channel_order():
    templates, dense = _bank()
    rec, _ = _recording(dense)
    flipped = rec.select_channels(rec.channel_ids[::-1])
    with pytest.raises(ValueError, match="same channels"):
        PerGroupMatching(flipped, templates, method="circus-omp", **AMP)
