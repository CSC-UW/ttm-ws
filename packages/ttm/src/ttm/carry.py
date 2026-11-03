"""Per-tetrode carry-forward template matching over a long recording.

The sort of ``docs/sorting_pipeline.md``, steps 3 to 7, on a preprocessed
binary (``ttm.sorting.preprocess_for_sorting``):

- ``whitening_matrix``: one 4 x 4 ZCA block per tetrode from the first
  ``whiten_s`` seconds, samples near threshold crossings left out.
- ``seed_bank``: MountainSort5 on the first window, units with enough snr and
  spikes, near-identical units of a tetrode merged.
- ``run_carry``: per window, circus-omp on each tetrode's whitened channels
  with only its templates, admission of each fit in noise units, a relaxed
  rule inside laser pulses for units that fire through them, and each
  template blended toward the mean of its window's clean lead spikes.
  Every ``reseed_every`` windows a MountainSort5 sort of the window adds
  units unlike any carried on their tetrode. Each window writes its spikes
  and banks, so a run resumes at the first window without them.
- ``export_sorting``: the accepted spikes of every window as one sorting.

Parameters are those of ``CarryParams``. Templates are (n_samples, 4) on the
unit's tetrode, with the spike at sample ``nbefore``; inside the carry they are
held in whitened space (``T @ W``).
"""

from __future__ import annotations

import dataclasses
import json
import pathlib
import time

import numba
import numpy as np
import spikeinterface as si
from spikeinterface.core import Templates

from ttm import capture as cp

REASON_ACCEPTED, REASON_FLOOR, REASON_SHAPE, REASON_REFRACTORY, REASON_BAD_PERIOD = 0, 1, 2, 3, 5


@dataclasses.dataclass(frozen=True, kw_only=True)
class CarryParams:
    """Parameters of the carry-forward sort (the values adopted on 2026-11-03)."""

    window_s: float = 1800.0
    nbefore: int = 30
    nafter: int = 60
    whiten_s: float = 300.0
    whiten_detect_mad: float = 4.0
    whiten_pad: int = 30
    noise_block_s: float = 10.0
    noise_n_blocks: int = 6
    method: str = "circus-omp"
    match_amp_floor: float = 0.2
    match_snr_floor: float = 2.5
    snr_floor: float = 4.5
    cos_floor: float = 0.5
    rel_floor: float = 0.5
    rel_min_lead: int = 30
    lead_ms: float = 10.0
    refractory_ms: float = 1.5
    pulse_post_ms: float = 100.0
    pulse_snr: float = 2.5
    pulse_rel: float = 0.3
    pulse_cos: float = 0.5
    driven_start_ms: float = 15.0
    driven_ratio: float = 2.0
    driven_p: float = 1e-4
    driven_min_spikes: int = 20
    blend: float = 0.2
    rebuild_max_spikes: int = 500
    rebuild_min_accepted: int = 100
    rebuild_min_selected: int = 50
    rebuild_seed: int = 0
    reseed_every: int = 12
    detect_threshold: float = 5.5
    whitening_seed: int = 42
    seed_min_snr: float = 5.0
    seed_min_spikes: int = 100
    seed_dedup_cos: float = 0.95
    reseed_add_cos: float = 0.8
    max_shift_samples: int = 10
    template_seed: int = 0

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def in_intervals(frames: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Whether each frame lies in one of the sorted, disjoint intervals [starts, ends)."""
    i = np.searchsorted(starts, frames, side="right") - 1
    return (i >= 0) & (frames < ends[np.clip(i, 0, None)])


def read_pulses(path, fs: float) -> tuple[np.ndarray, np.ndarray]:
    """Onset and offset frames, sorted by onset, from a CSV with ``onset_s`` and ``offset_s`` (session seconds)."""
    import pandas as pd

    p = pd.read_csv(path).sort_values("onset_s", kind="stable")
    return (p.onset_s.to_numpy() * fs).astype(np.int64), (p.offset_s.to_numpy() * fs).astype(np.int64)


def read_bad_periods(path, fs: float) -> dict[int, list[tuple[int, int]]]:
    """``{group: [(start_frame, end_frame), ...]}`` from a ``ttm-bad-periods`` table."""
    from ttm.bad_periods import group_intervals, read_table

    return group_intervals(read_table(path), fs, 10 ** 12)


# --- bank ------------------------------------------------------------------------------------------

@dataclasses.dataclass
class Bank:
    """Templates of one carry: ``T[k]`` is unit ``ids[k]``'s (n_samples, 4) template on tetrode ``group[k]``."""

    ids: np.ndarray
    group: np.ndarray
    T: np.ndarray

    def copy(self) -> Bank:
        return Bank(self.ids.copy(), self.group.copy(), self.T.copy())

    def templates(self, recording, nbefore: int) -> Templates:
        groups = np.asarray(recording.get_property("group")).astype(int)
        mask = np.stack([groups == g for g in self.group])
        return Templates(templates_array=self.T.astype(np.float32), sampling_frequency=recording.get_sampling_frequency(),
                         nbefore=nbefore, is_in_uV=False, sparsity_mask=mask,
                         channel_ids=np.asarray(recording.channel_ids), unit_ids=self.ids.copy(), probe=None,
                         check_for_consistent_sparsity=True)

    def whitened(self, W: np.ndarray, groups: np.ndarray) -> Bank:
        """Templates moved to the space of ``x @ W`` (``W`` over all channels, ``groups`` per channel)."""
        out = self.copy()
        for k, g in enumerate(self.group):
            ch = np.flatnonzero(groups == g)
            out.T[k] = self.T[k] @ W[np.ix_(ch, ch)]
        return out

    def save(self, path, **extra) -> None:
        np.savez(path, ids=self.ids, group=self.group, T=self.T, **extra)

    @classmethod
    def load(cls, path) -> Bank:
        z = np.load(path)
        return cls(z["ids"], z["group"], z["T"])


# --- whitening -------------------------------------------------------------------------------------

def whitening_matrix(recording, params: CarryParams = CarryParams()) -> np.ndarray:
    """Block-diagonal ZCA whitening over all channels, one block per tetrode, as float32.

    Each block is C^(-1/2) of the covariance of the tetrode's samples in the first ``whiten_s``
    seconds that lie more than ``whiten_pad`` samples from any crossing of -``whiten_detect_mad``
    times the excised MAD noise. Applied as ``x @ W``.
    """
    groups = np.asarray(recording.get_property("group")).astype(int)
    end = int(params.whiten_s * recording.get_sampling_frequency())
    W = np.zeros((len(groups), len(groups)), np.float32)
    for g in np.unique(groups):
        x = cp.load_tetrode_traces(recording, int(g), 0, end).astype(np.float64)
        noise, _ = cp.noise_mad_excised(x)
        keep = cp.excision_mask(x, noise, detect_mad=params.whiten_detect_mad, pad_frames=params.whiten_pad)
        ev, V = np.linalg.eigh(np.cov(x[keep].T))
        ch = np.flatnonzero(groups == g)
        W[np.ix_(ch, ch)] = V @ np.diag(ev ** -0.5) @ V.T
    return W


def whitened_recording(recording, W: np.ndarray):
    from spikeinterface.preprocessing import whiten

    return whiten(recording, W=W, M=None, apply_mean=False, dtype="float32")


# --- MountainSort5 seeds -----------------------------------------------------------------------------

def unit_templates(sorting, recording, *, n_jobs: int = 16, with_snr: bool = True, seed: int = 0):
    """Dense (n_units, 90, n_channels) templates (1 ms before, 2 ms after, raw units) of each unit on
    its tetrode, and each unit's snr (None without ``with_snr``), from a group-sparse in-memory analyzer
    with seeded spike sampling."""
    from spikeinterface.core.template_tools import get_dense_templates_array

    from ttm.sorting import create_group_sparse_analyzer

    az = create_group_sparse_analyzer(sorting, recording, format="memory", return_in_uV=False, seed=seed)
    az.compute({"random_spikes": {"seed": seed}, "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
                "templates": {}, "noise_levels": {}}, n_jobs=n_jobs)
    snr = None
    if with_snr:
        az.compute({"quality_metrics": {"metric_names": ["snr", "firing_rate"]}})
        snr = az.get_extension("quality_metrics").get_data()["snr"].to_numpy()
    return get_dense_templates_array(az, return_in_uV=False), snr


def dedup_units(sorting, recording, *, cosine_min: float, max_shift_samples: int = 10, coincidence_ms: float = 0.3,
                seed: int = 0):
    """Units of one tetrode whose templates have shift-cosine ``cosine_min`` or more merged (union-find).

    Merged trains are unioned with spikes closer than ``coincidence_ms`` collapsed. Returns a
    ``NumpySorting`` with ids 0..n-1 and the ``group`` property.
    """
    from ttm.tracking import _UnionFind, cosine_from_templates

    dense, _ = unit_templates(sorting, recording, with_snr=False, seed=seed)
    rgroups = np.asarray(recording.get_property("group")).astype(int)
    uids = list(sorting.unit_ids)
    groups = np.asarray(sorting.get_property("group")).astype(int)
    union = _UnionFind()
    for u in uids:
        union.add(u)
    by_group: dict = {}
    for i, g in enumerate(groups):
        by_group.setdefault(int(g), []).append(i)
    for g, members in by_group.items():
        ch = rgroups == g
        for a in range(len(members)):
            ta = dense[members[a]][:, ch]
            for b in range(a + 1, len(members)):
                if cosine_from_templates(ta, dense[members[b]][:, ch], max_shift_samples=max_shift_samples) >= cosine_min:
                    union.union(uids[members[a]], uids[members[b]])
    comp: dict = {}
    for u in uids:
        comp.setdefault(union.find(u), []).append(u)
    fs = sorting.get_sampling_frequency()
    tol = int(coincidence_ms * 1e-3 * fs)
    trains, new_groups = {}, []
    for new_id, members in enumerate(comp.values()):
        allspk = np.sort(np.concatenate([sorting.get_unit_spike_train(u) for u in members]))
        if allspk.size and tol > 0:
            allspk = allspk[np.concatenate([[True], np.diff(allspk) > tol])]
        trains[new_id] = allspk
        new_groups.append(int(groups[uids.index(members[0])]))
    merged = si.NumpySorting.from_unit_dict([trains], sampling_frequency=fs)
    merged.set_property("group", np.asarray(new_groups))
    return merged


def sort_window(recording, sort_dir: pathlib.Path, params: CarryParams):
    """MountainSort5 (``ttm.tracking.sort_chunk``) of ``recording``, or the sort already in ``sort_dir``.

    Returns an in-memory sorting with int64 ids and the ``group`` property.
    """
    import spikeinterface.sorters as ss

    from ttm.tracking import sort_chunk, to_int_numpy_sorting

    sort_dir = pathlib.Path(sort_dir)
    if (sort_dir / "by_group").exists():
        parts = []
        for g in sorted(int(p.name) for p in (sort_dir / "by_group").iterdir()):
            sg = ss.read_sorter_folder(sort_dir / "by_group" / str(g), register_recording=False)
            sg.set_property("group", np.full(sg.get_num_units(), g))
            parts.append(sg)
        return to_int_numpy_sorting(si.aggregate_units(parts))
    return to_int_numpy_sorting(sort_chunk(recording, sort_dir, whitening_seed=params.whitening_seed,
                                           detect_threshold=params.detect_threshold))


def seed_bank(recording, sort_dir: pathlib.Path, params: CarryParams = CarryParams(), *, n_jobs: int = 32) -> Bank:
    """The starting bank (plain, not whitened) from MountainSort5 on the first window of ``recording``."""
    fs = recording.get_sampling_frequency()
    win0 = recording.frame_slice(0, min(int(params.window_s * fs), recording.get_num_frames()))
    win0.reset_times()
    ref0 = sort_window(win0, sort_dir, params)
    _, snr = unit_templates(ref0, win0, n_jobs=n_jobs, with_snr=True, seed=params.template_seed)
    nsp = np.array([len(ref0.get_unit_spike_train(u)) for u in ref0.unit_ids])
    conf = np.asarray(ref0.unit_ids)[(snr >= params.seed_min_snr) & (nsp >= params.seed_min_spikes)]
    seed = dedup_units(ref0.select_units(conf), win0, cosine_min=params.seed_dedup_cos,
                       max_shift_samples=params.max_shift_samples, seed=params.template_seed)
    dense, _ = unit_templates(seed, win0, n_jobs=n_jobs, with_snr=False, seed=params.template_seed)
    groups = np.asarray(seed.get_property("group")).astype(int)
    T = np.stack([dense[i][:, cp.tetrode_channel_indices(recording, g)] for i, g in enumerate(groups)])
    log(f"seed: {ref0.get_num_units()} units -> {len(conf)} with snr >= {params.seed_min_snr} and "
        f">= {params.seed_min_spikes} spikes -> {seed.get_num_units()} after merging at cosine {params.seed_dedup_cos}")
    return Bank(np.arange(len(groups), dtype=np.int64), groups, T.astype(np.float32))


def reseed(plain, whitened, bank: Bank, a: int, b: int, sort_dir: pathlib.Path, params: CarryParams, *,
           n_jobs: int) -> Bank:
    """``bank`` plus every cluster of a MountainSort5 sort of frames [a, b) of ``plain`` with snr
    ``seed_min_snr`` or more and ``seed_min_spikes`` or more whose whitened template has shift-cosine
    under ``reseed_add_cos`` to every template of its tetrode (ids from 2000 up)."""
    from ttm.tracking import cosine_from_templates

    win_p = plain.frame_slice(a, b)
    win_p.reset_times()
    win_w = whitened.frame_slice(a, b)
    win_w.reset_times()
    rs = sort_window(win_p, sort_dir, params)
    dense, snr = unit_templates(rs, win_w, n_jobs=n_jobs, with_snr=True, seed=params.template_seed)
    groups = np.asarray(whitened.get_property("group")).astype(int)
    rgrp = np.asarray(rs.get_property("group")).astype(int)
    out = bank.copy()
    next_id = max(2000, int(bank.ids.max()) + 1)
    added = []
    for j, cu in enumerate(rs.unit_ids):
        if snr[j] < params.seed_min_snr or len(rs.get_unit_spike_train(cu)) < params.seed_min_spikes:
            continue
        g = rgrp[j]
        t = dense[j][:, groups == g]
        best = max((cosine_from_templates(t, bt, max_shift_samples=params.max_shift_samples)
                    for bt in out.T[out.group == g]), default=-1.0)
        if best < params.reseed_add_cos:
            out = Bank(np.r_[out.ids, next_id], np.r_[out.group, g], np.concatenate([out.T, t[None].astype(np.float32)]))
            added.append((next_id, int(g), round(float(best), 3)))
            next_id += 1
    log(f"reseed after window {a // int(params.window_s * plain.get_sampling_frequency())}: "
        f"+{len(added)} units {added}")
    return out


# --- matching --------------------------------------------------------------------------------------

def amplitude_bounds(bank: Bank, sigma: np.ndarray, params: CarryParams) -> np.ndarray:
    """Per template (lower, upper) fit amplitude for circus-omp: max(``match_amp_floor``,
    ``match_snr_floor`` x noise / template norm), no upper bound."""
    tn = np.sqrt((bank.T.astype(np.float64) ** 2).sum(axis=(1, 2)))
    lo = np.maximum(params.match_amp_floor, params.match_snr_floor * sigma / tn)
    return np.stack([lo, np.full_like(lo, np.inf)], axis=1)


def match_window(recording, bank: Bank, amplitudes: np.ndarray, params: CarryParams, *, n_jobs: int) -> np.ndarray:
    """Spikes of every tetrode in ``recording`` (one window), each tetrode matched on its own."""
    from spikeinterface.core.node_pipeline import run_node_pipeline

    from ttm.matching import PerGroupMatching

    node = PerGroupMatching(recording, bank.templates(recording, params.nbefore), method=params.method,
                            amplitudes=amplitudes)
    spikes = run_node_pipeline(recording, [node], {"n_jobs": n_jobs, "chunk_duration": "1s", "progress_bar": False},
                               job_name="per-tetrode matching", squeeze_output=True, names=["spikes"])
    node.clean()
    return spikes


# --- per-tetrode numerics ----------------------------------------------------------------------------

@numba.njit(cache=True)
def add_spikes(out, frames, labels, amps, T, nb):
    """``out += sum_i amps[i] * T[labels[i]]`` placed with sample ``nb`` of the template at ``frames[i]``."""
    n = out.shape[0]
    ns = T.shape[1]
    nc = T.shape[2]
    for i in range(frames.size):
        s0 = frames[i] - nb
        for k in range(ns):
            j = s0 + k
            if j < 0 or j >= n:
                continue
            for c in range(nc):
                out[j, c] += amps[i] * T[labels[i], k, c]


def mf_noise(x: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Per template, the robust standard deviation (MAD / 0.6745) of the projection of ``x`` on the unit-norm template."""
    from scipy.signal import oaconvolve

    seg = x.astype(np.float64)
    out = np.zeros(len(T))
    for k, t in enumerate(T):
        th = t / max(np.linalg.norm(t), 1e-12)
        y = sum(oaconvolve(seg[:, c], th[::-1, c], mode="valid") for c in range(seg.shape[1]))
        out[k] = np.median(np.abs(y - np.median(y))) / 0.6745
    return out


def window_noise(recording, g: int, a: int, b: int, T: np.ndarray, params: CarryParams) -> np.ndarray:
    """``mf_noise`` of tetrode ``g``'s templates over ``noise_n_blocks`` blocks of ``noise_block_s`` spread over [a, b)."""
    nblk = int(params.noise_block_s * recording.get_sampling_frequency())
    starts = np.linspace(a, max(b - nblk, a), params.noise_n_blocks).astype(np.int64)
    x = np.concatenate([cp.load_tetrode_traces(recording, g, s, s + nblk) for s in starts])
    return mf_noise(x, T)


def spike_scores(x: np.ndarray, frames: np.ndarray, labels: np.ndarray, amps: np.ndarray, T: np.ndarray,
                 noise: np.ndarray, sigma: np.ndarray, nbefore: int, nafter: int) -> dict:
    """Per spike of one tetrode (``frames`` local to ``x``, ``labels`` index ``T``): snr (the
    neighbour-subtracted snippet's projection on the unit-norm template over ``sigma``), cosines of the
    raw and the neighbour-subtracted (peeled) snippet to the template, residual energy and the peeled
    trough in noise units. Also returns the residual traces under ``residual_traces``."""
    model = np.zeros_like(x)
    add_spikes(model, frames, labels, amps.astype(np.float32), T, nbefore)
    resid = x - model
    n = frames.size
    out = {k: np.full(n, np.nan) for k in ("snr", "cos_raw", "cos_pl", "resid", "trough")}
    ok = np.flatnonzero((frames >= nbefore) & (frames + nafter <= x.shape[0]))
    tn = np.sqrt((T.astype(np.float64) ** 2).sum(axis=(1, 2)))
    inv = 1.0 / noise[None, None, :]
    for b0 in range(0, ok.size, 100_000):
        sel = ok[b0:b0 + 100_000]
        lab = labels[sel]
        raw = cp.snippets(x, frames[sel], nbefore=nbefore, nafter=nafter)
        r = cp.snippets(resid, frames[sel], nbefore=nbefore, nafter=nafter)
        t = T[lab]
        pl = r + amps[sel, None, None].astype(np.float32) * t
        proj = np.einsum("ntc,ntc->n", pl, t) / tn[lab]
        out["snr"][sel] = proj / sigma[lab]
        out["cos_pl"][sel] = proj / np.maximum(np.sqrt(np.einsum("ntc,ntc->n", pl, pl)), 1e-12)
        out["cos_raw"][sel] = (np.einsum("ntc,ntc->n", raw, t) / tn[lab]
                               / np.maximum(np.sqrt(np.einsum("ntc,ntc->n", raw, raw)), 1e-12))
        out["resid"][sel] = ((r * inv) ** 2).mean(axis=(1, 2))
        out["trough"][sel] = -(pl * inv).min(axis=(1, 2))
    out["residual_traces"] = resid
    return out


def light_driven(fr, un, on, off, dur_s: float, fs: float, params: CarryParams) -> list[int]:
    """Labels that keep firing through the pulses: spikes from ``driven_start_ms`` after onset to offset
    at least ``driven_ratio`` times their rate outside pulses (onset to offset plus ``pulse_post_ms``),
    Poisson p under ``driven_p``, and ``driven_min_spikes`` or more."""
    from scipy.stats import poisson

    out = []
    post = int(params.pulse_post_ms * 1e-3 * fs)
    st = on + int(params.driven_start_ms * 1e-3 * fs)
    span = (off - st).sum() / fs
    out_s = dur_s - ((off + post) - on).sum() / fs
    for u in np.unique(un):
        f = fr[un == u]
        n_in = int(in_intervals(f, st, off).sum())
        rate = (~in_intervals(f, on, off + post)).sum() / out_s
        if (n_in >= params.driven_min_spikes and n_in >= params.driven_ratio * rate * span
                and poisson.sf(n_in - 1, rate * span) < params.driven_p):
            out.append(int(u))
    return out


def tetrode_window(x: np.ndarray, a: int, frames_abs: np.ndarray, labels: np.ndarray, amps: np.ndarray, T: np.ndarray,
                   sigma: np.ndarray, params: CarryParams, *, fs: float, pulses=None, bad=(), window_index: int = 0) -> dict:
    """Admission and template rebuild for one tetrode in one window.

    ``x`` is the tetrode's whitened traces of the window starting at frame ``a``; ``frames_abs`` the
    matcher's spikes (absolute frames, sorted), ``labels`` indices into ``T`` (the tetrode's whitened
    templates), ``sigma`` each template's matched-filter noise. ``pulses`` are (onset, offset) frames
    or None; ``bad`` the tetrode's bad periods as (start, end) frames. Returns per-spike ``scores`` and
    ``reason`` (0 accepted, 1 snr or relative floor, 2 shape, 3 refractory, 5 bad period), ``newT``
    (NaN where a template was not rebuilt), ``n_used``, ``noise`` and ``driven`` (indices into ``T``).
    """
    nb, na = params.nbefore, params.nafter
    noise, _ = cp.noise_mad_excised(x)
    fl = frames_abs - a
    sc = spike_scores(x, fl, labels, amps, T, noise, sigma, nb, na)
    resid = sc.pop("residual_traces")
    reason = np.zeros(fl.size, np.int8)
    reason[~(sc["snr"] >= params.snr_floor)] = REASON_FLOOR
    if len(bad):
        st, en = np.array(bad).T
        reason[(reason == 0) & in_intervals(frames_abs, st, en)] = REASON_BAD_PERIOD
    reason[(reason == 0) & ~(sc["cos_pl"] >= params.cos_floor)] = REASON_SHAPE
    post = int(params.pulse_post_ms * 1e-3 * fs)
    outside = np.ones(fl.size, bool)
    if pulses is not None:
        outside = ~in_intervals(frames_abs, pulses[0], pulses[1] + post)
    lead_frames = int(params.lead_ms * 1e-3 * fs)
    amp_ref = np.ones(len(T))
    for k_u in range(len(T)):
        idx = np.flatnonzero((labels == k_u) & (reason == 0))
        if idx.size < params.rel_min_lead:
            continue
        idx = idx[np.argsort(fl[idx], kind="stable")]
        lead = np.r_[True, np.diff(fl[idx]) > lead_frames] & outside[idx]
        if lead.sum() >= params.rel_min_lead:
            amp_ref[k_u] = float(np.median(amps[idx[lead]]))
    reason[(reason == 0) & ~(amps >= params.rel_floor * amp_ref[labels])] = REASON_FLOOR
    score = np.nan_to_num(sc["snr"], nan=-np.inf)
    rf = int(params.refractory_ms * 1e-3 * fs)

    def refractory(units):
        for k_u in units:
            idx = np.flatnonzero((labels == k_u) & (reason == 0))
            if idx.size > 1:
                o = idx[np.argsort(fl[idx], kind="stable")]
                keep = cp.resolve_refractory(fl[o], score[o], rf)
                reason[o[~keep]] = REASON_REFRACTORY

    refractory(range(len(T)))
    driven = []
    if pulses is not None:
        on, off = pulses
        m = (on >= a) & (on < a + x.shape[0])
        if m.any():
            on, off = on[m] - a, off[m] - a
            near_p = in_intervals(fl, on, off + post)
            ok = ((reason != REASON_BAD_PERIOD) & (sc["snr"] >= params.pulse_snr)
                  & (amps >= params.pulse_rel * amp_ref[labels]) & (sc["cos_pl"] >= params.pulse_cos))
            driven = light_driven(fl[ok], labels[ok], on, off, x.shape[0] / fs, fs, params)
            relax = np.isin(labels, driven) & near_p & (reason != 0) & ok
            reason[relax] = 0
            refractory(driven)
    acc = reason == 0
    newT = np.full_like(T, np.nan)
    n_used = np.zeros(len(T), np.int64)
    rng = np.random.default_rng(params.rebuild_seed + window_index)
    rej = ~acc
    if rej.any():
        add_spikes(resid, fl[rej], labels[rej], amps[rej].astype(np.float32), T, nb)
    for k_u in range(len(T)):
        if acc[labels == k_u].sum() < params.rebuild_min_accepted:
            continue
        idx = np.flatnonzero(acc & (labels == k_u))
        f = frames_abs[idx]
        keep = np.r_[np.inf, np.diff(f)] > lead_frames
        if pulses is not None:
            keep &= ~in_intervals(f, pulses[0], pulses[1] + post)
        if len(bad):
            st, en = np.array(bad).T
            keep &= ~in_intervals(f, st, en)
        sel = idx[keep]
        sel = sel[(fl[sel] >= nb) & (fl[sel] + na <= x.shape[0])]
        if sel.size < params.rebuild_min_selected:
            continue
        if sel.size > params.rebuild_max_spikes:
            sel = np.sort(rng.choice(sel, params.rebuild_max_spikes, replace=False))
        sn = cp.snippets(resid, fl[sel], nbefore=nb, nafter=na) + amps[sel, None, None].astype(np.float32) * T[k_u][None]
        est = sn.astype(np.float64).mean(axis=0)
        newT[k_u] = (1 - params.blend) * T[k_u] + params.blend * est
        n_used[k_u] = sel.size
    return {"scores": sc, "reason": reason, "newT": newT, "n_used": n_used, "noise": noise, "driven": driven}


# --- the carry ---------------------------------------------------------------------------------------

def _open(binary, W):
    rec = si.load(binary)
    return rec, whitened_recording(rec, W)


def _tetrode_job(binary, W, g, a, b, frames_abs, local, amps, T, sigma, params, pulses, bad, window_index):
    _, rec = _open(binary, W)
    x = cp.load_tetrode_traces(rec, g, a, b)
    o = np.argsort(frames_abs, kind="stable")
    inv = np.empty_like(o)
    inv[o] = np.arange(o.size)
    r = tetrode_window(x, a, frames_abs[o], local[o], amps[o], T, sigma, params, fs=rec.get_sampling_frequency(),
                       pulses=pulses, bad=bad, window_index=window_index)
    r["scores"] = {k: v[inv] for k, v in r["scores"].items()}
    r["reason"] = r["reason"][inv]
    return r


def run_carry(binary, W: np.ndarray, bank0: Bank, out_dir: pathlib.Path, *, params: CarryParams = CarryParams(),
              windows: range | None = None, pulses=None, bad_periods: dict | None = None,
              reseed_dir: pathlib.Path | None = None, n_jobs: int = 16, workers: int = 8) -> None:
    """Carry ``bank0`` (whitened templates) over the windows of the binary at ``binary`` whitened with ``W``.

    Per window ``w`` writes ``spikes_wNN.npz`` (frame, unit, a, snr, cos_raw, cos_pl, resid, trough,
    reason), ``bankin_wNN.npz`` (the bank matched, with each template's noise ``sigma``) and
    ``bank_wNN.npz`` (the bank after the rebuild and any re-seed). A window whose spikes and bank exist
    is skipped and its bank read, so a run resumes. Re-seed sorts go to ``reseed_dir/reseed_w<w>``
    (default ``out_dir/reseed_sorts``) and are read from there when present.
    """
    from joblib import Parallel, delayed

    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    reseed_dir = pathlib.Path(reseed_dir) if reseed_dir is not None else out_dir / "reseed_sorts"
    plain, rec = _open(binary, W)
    fs = rec.get_sampling_frequency()
    win_frames = int(params.window_s * fs)
    n_frames = rec.get_num_frames()
    if windows is None:
        windows = range(int(np.ceil(n_frames / win_frames)))
    bad_periods = bad_periods or {}
    bank = bank0.copy()
    for w in windows:
        f_sp, f_bk = out_dir / f"spikes_w{w:02d}.npz", out_dir / f"bank_w{w:02d}.npz"
        if f_sp.exists() and f_bk.exists():
            bank = Bank.load(f_bk)
            continue
        t0 = time.time()
        a, b = w * win_frames, min((w + 1) * win_frames, n_frames)
        groups = sorted(set(bank.group.tolist()))
        sig = Parallel(n_jobs=workers)(
            delayed(window_noise)(rec, g, a, b, bank.T[bank.group == g], params) for g in groups)
        sigma = np.zeros(len(bank.ids))
        for g, s in zip(groups, sig):
            sigma[bank.group == g] = s
        win = rec.frame_slice(a, b)
        sp = match_window(win, bank, amplitude_bounds(bank, sigma, params), params, n_jobs=n_jobs)
        t_match = time.time() - t0
        fr = sp["sample_index"].astype(np.int64) + a
        ci = sp["cluster_index"].astype(np.int64)
        amp = sp["amplitude"].astype(np.float64)
        jobs = []
        for g in groups:
            gi = np.flatnonzero(bank.group == g)
            m = bank.group[ci] == g
            jobs.append(delayed(_tetrode_job)(binary, W, g, a, b, fr[m], np.searchsorted(gi, ci[m]), amp[m],
                                              bank.T[gi], sigma[gi], params, pulses, bad_periods.get(int(g), []), w))
        res = Parallel(n_jobs=workers)(jobs)
        cols = {k: [] for k in ("frame", "unit", "a", "snr", "cos_raw", "cos_pl", "resid", "trough", "reason")}
        newbank = bank.copy()
        n_used = np.zeros(len(bank.ids), np.int64)
        noise, driven = {}, []
        for g, r in zip(groups, res):
            gi = np.flatnonzero(bank.group == g)
            m = bank.group[ci] == g
            cols["frame"].append(fr[m])
            cols["unit"].append(bank.ids[ci[m]])
            cols["a"].append(amp[m])
            for k in ("snr", "cos_raw", "cos_pl", "resid", "trough"):
                cols[k].append(r["scores"][k])
            cols["reason"].append(r["reason"])
            for j, k_u in enumerate(gi):
                if np.isfinite(r["newT"][j]).all():
                    newbank.T[k_u] = r["newT"][j].astype(np.float32)
                    n_used[k_u] = r["n_used"][j]
            noise[g] = r["noise"]
            driven += [int(bank.ids[gi[k]]) for k in r["driven"]]
        arrays = {k: np.concatenate(v) for k, v in cols.items()}
        o = np.argsort(arrays["frame"], kind="stable")
        np.savez(f_sp, **{k: v[o] for k, v in arrays.items()})
        bank.save(out_dir / f"bankin_w{w:02d}.npz", sigma=sigma)
        driven = np.array(sorted(driven), np.int64)
        extra = dict(noise=np.stack([noise[g] for g in groups]), noise_groups=np.array(groups), driven=driven)
        acc = arrays["reason"] == 0
        if params.reseed_every and w > 0 and w % params.reseed_every == 0:
            newbank = reseed(plain, rec, newbank, a, b, reseed_dir / f"reseed_w{w}", params, n_jobs=n_jobs)
            n_used = np.r_[n_used, np.zeros(len(newbank.ids) - len(n_used), np.int64)]
        newbank.save(f_bk, n_used=n_used, **extra)
        log(f"w{w:02d}: {acc.sum():,} accepted of {acc.size:,} fits; driven {driven.tolist()}; "
            f"{(n_used > 0).sum()}/{len(bank.ids)} templates re-estimated; match {t_match / 60:.1f} min, "
            f"total {(time.time() - t0) / 60:.1f} min")
        bank = newbank


def export_sorting(out_dir: pathlib.Path, folder: pathlib.Path, fs: float):
    """The accepted spikes of every ``spikes_wNN.npz`` in ``out_dir`` as one ``NumpySorting`` saved to
    ``folder``, with the last bank's unit ids and ``group`` property."""
    out_dir = pathlib.Path(out_dir)
    files = sorted(out_dir.glob("spikes_w*.npz"))
    fr, un = [], []
    for f in files:
        z = np.load(f)
        acc = z["reason"] == 0
        fr.append(z["frame"][acc])
        un.append(z["unit"][acc])
    fr, un = np.concatenate(fr), np.concatenate(un)
    last = Bank.load(sorted(out_dir.glob("bank_w*.npz"))[-1])
    o = np.argsort(fr, kind="stable")
    s = si.NumpySorting.from_samples_and_labels([fr[o]], [un[o]], sampling_frequency=fs, unit_ids=last.ids)
    s.set_property("group", last.group)
    s.save(folder=folder, overwrite=True)
    log(f"{len(files)} windows, {fr.size:,} spikes, {len(last.ids)} units -> {folder}")
    return s


def write_run_info(path: pathlib.Path, params: CarryParams, **inputs) -> None:
    """``carry_params.json``: the parameters, the inputs and the package versions; refuses to
    overwrite one with different parameters, so a resumed run cannot mix them."""
    import importlib.metadata as md

    path = pathlib.Path(path)
    if path.exists():
        old = json.loads(path.read_text())
        if old.get("params") != params.to_dict():
            raise ValueError(f"{path} holds different parameters; use another run directory")
        return
    versions = {}
    for name in ("ttm", "spikeinterface", "mountainsort5", "numpy"):
        try:
            versions[name] = md.version(name)
        except md.PackageNotFoundError:
            versions[name] = None
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"params": params.to_dict(), "inputs": {k: str(v) for k, v in inputs.items()},
                                "versions": versions, "written": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2) + "\n")
