"""Shared code of the carry-forward bake-off on TTM-003's 48 h recording.

A carry-forward runs a template matcher over 1800 s windows and rebuilds each
template from the spikes it took in the window before the next. This module
holds one carry loop with the parts the bake-off varies:

- the matcher: circus-omp on all channels at once (June) or on each tetrode's
  channels on its own, with a per-template amplitude floor;
- admission: a floor relative to the template (June: a >= 0.8), a floor on
  the matched-filter signal-to-noise (the event's projection on the
  unit-norm template over that projection's noise, measured per template and
  window), a shape gate and a refractory rule, the last two applied after
  matching;
- re-estimation: which accepted spikes a template is rebuilt from, from what
  snippets (raw, or with every other fitted spike subtracted), with which
  estimator, and how much of the old template is kept.

Per spike it records the fit amplitude, the signal-to-noise, cosines to the
assigned template of the raw and of the peeled snippet (the residual plus
the spike's own fit), the residual energy and the trough, so admission
rules can also be compared after the fact. Paths are literal tononi-2 paths.
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
from spikeinterface.core.node_pipeline import run_node_pipeline
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates
from spikeinterface.sortingcomponents.matching.base import BaseTemplateMatching, _base_matching_dtype
from spikeinterface.sortingcomponents.matching.method_list import matching_methods

from ttm import capture as cp

FS = 30000.0
WIN = int(1800 * FS)
NB, NA = 30, 60
N_SAMP = NB + NA
SESSION = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12")
SORTINGS = SESSION / "sortings_seed42_pcafix"
BINARY = SORTINGS / "mp_baseline" / "binary"
BAKE_ROOT = SORTINGS / "carry_bakeoff"
BAKE = BAKE_ROOT / "v2"
HYBRID = BAKE / "hybrid_binary"
BAD_PERIODS = SESSION / "bad_periods" / "bad_periods.csv"
RESULTS = pathlib.Path("/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/subjects/ttm-003/results")
PULSES = RESULTS / "stimulation_pulses_2026-10-20_09-02-12.csv"
BLOCKS = RESULTS / "stimulation_blocks_2026-10-20_09-02-12.csv"
JOB = {"chunk_duration": "1s", "progress_bar": False}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --- inputs ----------------------------------------------------------------------------------

def laser_pulses() -> tuple[np.ndarray, np.ndarray]:
    """Onset and offset frames of every delivered, non-glitch laser pulse."""
    import pandas as pd

    blocks = pd.read_csv(BLOCKS)
    delivered = set(blocks.loc[blocks.delivered, "ttl_train"])
    p = pd.read_csv(PULSES)
    p = p[~p.glitch & p.train.isin(delivered)]
    return (p.onset_s.to_numpy() * FS).astype(np.int64), (p.offset_s.to_numpy() * FS).astype(np.int64)


def in_intervals(frames: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Whether each frame lies in one of the sorted, disjoint intervals [starts, ends)."""
    i = np.searchsorted(starts, frames, side="right") - 1
    return (i >= 0) & (frames < ends[np.clip(i, 0, None)])


@dataclasses.dataclass
class Bank:
    """Templates of one carry: ``T[k]`` is unit ``ids[k]``'s (N_SAMP, 4) template on tetrode ``group[k]``."""
    ids: np.ndarray
    group: np.ndarray
    T: np.ndarray

    def copy(self) -> "Bank":
        return Bank(self.ids.copy(), self.group.copy(), self.T.copy())

    def templates(self, recording) -> Templates:
        groups = np.asarray(recording.get_property("group")).astype(int)
        mask = np.stack([groups == g for g in self.group])
        return Templates(templates_array=self.T.astype(np.float32), sampling_frequency=FS, nbefore=NB,
                         is_in_uV=False, sparsity_mask=mask, channel_ids=np.asarray(recording.channel_ids),
                         unit_ids=self.ids.copy(), probe=None, check_for_consistent_sparsity=True)

    def save(self, path, **extra) -> None:
        np.savez(path, ids=self.ids, group=self.group, T=self.T, **extra)

    @classmethod
    def load(cls, path) -> "Bank":
        z = np.load(path)
        return cls(z["ids"], z["group"], z["T"])


def bank_from_dense(dense: np.ndarray, mask: np.ndarray, ids: np.ndarray, rec_groups: np.ndarray) -> Bank:
    group = np.array([int(rec_groups[np.flatnonzero(m)[0]]) for m in mask])
    T = np.stack([dense[i][:, rec_groups == g] for i, g in enumerate(group)]).astype(np.float32)
    return Bank(np.asarray(ids).astype(np.int64), group, T)


# --- per-tetrode matching with per-template settings ------------------------------------------

class PerGroupMatchingArrays(BaseTemplateMatching):
    """``ttm.matching.PerGroupMatching`` that also slices per-template ``amplitudes`` (an (n, 2)
    array over the whole bank) to each group's templates."""

    def __init__(self, recording, templates, return_output=True, method="circus-omp", group_property="group",
                 **method_kwargs):
        from ttm.matching import group_templates, template_groups

        BaseTemplateMatching.__init__(self, recording, templates, return_output=return_output)
        groups = np.asarray(recording.get_property(group_property))
        unit_group = template_groups(templates, groups)
        amps = method_kwargs.pop("amplitudes", None)
        self.parts = []
        for g in np.unique(unit_group):
            channels = np.flatnonzero(groups == g)
            units = np.flatnonzero(unit_group == g)
            sub = recording.select_channels(np.asarray(recording.channel_ids)[channels])
            kw = dict(method_kwargs)
            if amps is not None:
                a = np.asarray(amps, dtype=np.float64)
                kw["amplitudes"] = a[units] if a.ndim == 2 else list(a)
            node = matching_methods[method](sub, templates=group_templates(templates, units, channels), **kw)
            self.parts.append((channels, units, node))
        self.margin = max(node.get_margin() for _, _, node in self.parts)

    def get_margin(self):
        return self.margin

    def compute_matching(self, traces, start_frame, end_frame, segment_index):
        found = []
        for channels, units, node in self.parts:
            s = node.compute_matching(np.ascontiguousarray(traces[:, channels]), start_frame, end_frame,
                                      segment_index)
            s["cluster_index"] = units[s["cluster_index"]]
            s["channel_index"] = channels[s["channel_index"]]
            found.append(s)
        spikes = np.concatenate(found)
        return spikes[np.argsort(spikes["sample_index"], kind="stable")]

    def clean(self):
        for _, _, node in self.parts:
            node.clean()


def match(recording, templates: Templates, *, per_tetrode: bool, amplitudes, n_jobs: int,
          method: str = "circus-omp", extra: dict | None = None) -> np.ndarray:
    kw = {"amplitudes": amplitudes, **(extra or {})}
    job = {"n_jobs": n_jobs, **JOB}
    if not per_tetrode:
        if isinstance(amplitudes, np.ndarray):
            kw["amplitudes"] = amplitudes
        return find_spikes_from_templates(recording, templates, method=method, method_kwargs=kw, job_kwargs=job)
    node = PerGroupMatchingArrays(recording, templates, method=method, **kw)
    spikes = run_node_pipeline(recording, [node], job, job_name="per-tetrode matching", squeeze_output=True,
                               names=["spikes"])
    node.clean()
    return spikes


# --- per-tetrode numerics ----------------------------------------------------------------------

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


def mf_noise(x: np.ndarray, T: np.ndarray, *, block_s: float = 10.0, n_blocks: int = 6) -> np.ndarray:
    """Per template, the noise of its matched filter: the robust standard deviation (MAD / 0.6745) of
    the projection of ``x`` on the unit-norm template, over ``n_blocks`` blocks spread over ``x``."""
    from scipy.signal import oaconvolve

    nblk = int(block_s * FS)
    starts = np.linspace(0, max(x.shape[0] - nblk, 0), n_blocks).astype(np.int64)
    seg = np.concatenate([x[s:s + nblk] for s in starts]).astype(np.float64)
    out = np.zeros(len(T))
    for k, t in enumerate(T):
        th = t / max(np.linalg.norm(t), 1e-12)
        y = sum(oaconvolve(seg[:, c], th[::-1, c], mode="valid") for c in range(seg.shape[1]))
        out[k] = np.median(np.abs(y - np.median(y))) / 0.6745
    return out


def spike_scores(x: np.ndarray, frames: np.ndarray, labels: np.ndarray, amps: np.ndarray, T: np.ndarray,
                 noise: np.ndarray, sigma: np.ndarray) -> dict:
    """Per spike of one tetrode (``frames`` local to ``x``, ``labels`` index ``T``): the peeled
    snippet's projection on the unit-norm template over ``sigma`` (snr), cosines of the raw and the
    peeled snippet to the template, residual energy (mean of (residual / noise)^2 over the window)
    and the peeled trough in noise units. Also returns the residual traces."""
    model = np.zeros_like(x)
    add_spikes(model, frames, labels, amps.astype(np.float32), T, NB)
    resid = x - model
    n = frames.size
    out = {k: np.full(n, np.nan) for k in ("snr", "cos_raw", "cos_pl", "resid", "trough")}
    ok = np.flatnonzero((frames >= NB) & (frames + NA <= x.shape[0]))
    tn = np.sqrt((T.astype(np.float64) ** 2).sum(axis=(1, 2)))
    inv = 1.0 / noise[None, None, :]
    for b0 in range(0, ok.size, 100_000):
        sel = ok[b0:b0 + 100_000]
        lab = labels[sel]
        raw = cp.snippets(x, frames[sel])
        r = cp.snippets(resid, frames[sel])
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


# --- re-estimation recipes --------------------------------------------------------------------

@dataclasses.dataclass
class Recipe:
    """How a template is rebuilt from its accepted spikes in a window.

    select: subset of the accepted spikes ("all"; "lead": no spike of the unit in the preceding
        ``lead_ms``; "nostim": outside every laser pulse and ``post_ms`` after it; "full": fit
        amplitude a >= 0.8, June's circus population, which truncates the amplitude distribution;
        "lead_nostim"; "isolated": no spike of another unit of the tetrode within ``iso_ms``).
    snippet: "raw" or "peeled" (every other fitted spike subtracted).
    estimator: "mean", "median" or "trimmed" (20 % trimmed mean per sample and wire).
    blend: weight of the new estimate; 1 replaces the template (June), 0 freezes it.
    """
    select: str = "all"
    snippet: str = "raw"
    estimator: str = "mean"
    blend: float = 1.0
    max_spikes: int = 500
    min_spikes: int = 100
    scale: str = "mean"
    lead_ms: float = 10.0
    post_ms: float = 100.0
    iso_ms: float = 1.5
    seed: int = 0


_BAD = {}


def bad_intervals() -> dict:
    """``{group: [(start, end), ...]}`` frames of the session's bad-period table (cached)."""
    if not _BAD:
        from ttm.bad_periods import group_intervals, read_table

        _BAD.update(group_intervals(read_table(BAD_PERIODS), FS, 10 ** 12))
        _BAD.setdefault(-1, [])
    return _BAD


def light_driven(fr, un, on, off, dur_s, *, start_ms=15.0, min_ratio=2.0, p_max=1e-4, min_spikes=20):
    """Labels that keep firing through the pulses: spikes from ``start_ms`` after onset to offset at
    least ``min_ratio`` times their rate outside pulses and onset to offset plus 100 ms (Poisson,
    p < ``p_max``). Frames and pulses in the same frame of reference; ``dur_s`` the span's length."""
    from scipy.stats import poisson

    out = []
    st = on + int(start_ms * 1e-3 * FS)
    span = (off - st).sum() / FS
    out_s = dur_s - ((off + int(0.1 * FS)) - on).sum() / FS
    for u in np.unique(un):
        f = fr[un == u]
        n_in = int(in_intervals(f, st, off).sum())
        rate = (~in_intervals(f, on, off + int(0.1 * FS))).sum() / out_s
        if n_in >= min_spikes and n_in >= min_ratio * rate * span and poisson.sf(n_in - 1, rate * span) < p_max:
            out.append(int(u))
    return out


def truncated_mean(x: np.ndarray, c: float) -> float:
    """Maximum-likelihood mean of a normal truncated below at ``c`` (sd at least 1, the noise), from
    the samples ``x >= c`` that survived the truncation."""
    from scipy.optimize import minimize
    from scipy.stats import norm

    x = x[np.isfinite(x) & (x >= c)]
    if x.size < 20:
        return float(np.mean(x)) if x.size else np.nan

    def nll(th):
        mu, ls = th
        sd = 1.0 + np.exp(ls)
        return -(norm.logpdf(x, mu, sd).sum() - x.size * norm.logsf(c, mu, sd))

    m0, s0 = float(np.mean(x)), max(float(np.std(x)) - 1.0, 0.1)
    r = minimize(nll, [m0, np.log(s0)], method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
    return float(r.x[0])


def select_spikes(recipe: Recipe, frames_abs: np.ndarray, labels: np.ndarray, amps: np.ndarray, u: int,
                  pulses: tuple[np.ndarray, np.ndarray] | None, group: int = -1) -> np.ndarray:
    """Indices (into the tetrode's spikes) of unit ``u``'s spikes a template is rebuilt from."""
    idx = np.flatnonzero(labels == u)
    if idx.size == 0 or recipe.select == "all":
        return idx
    keep = np.ones(idx.size, bool)
    f = frames_abs[idx]
    if "lead" in recipe.select:
        gap = np.r_[np.inf, np.diff(f)]
        keep &= gap > recipe.lead_ms * 1e-3 * FS
    if "nostim" in recipe.select and pulses is not None:
        on, off = pulses
        keep &= ~in_intervals(f, on, off + int(recipe.post_ms * 1e-3 * FS))
    if "nobad" in recipe.select:
        iv = bad_intervals().get(int(group), [])
        if iv:
            st, en = np.array(iv).T
            keep &= ~in_intervals(f, st, en)
    if recipe.select == "full":
        keep &= amps[idx] >= 0.8
    if recipe.select == "isolated":
        others = np.sort(frames_abs[labels != u])
        if others.size:
            j = np.clip(np.searchsorted(others, f), 1, others.size - 1)
            d = np.minimum(np.abs(others[j] - f), np.abs(others[j - 1] - f))
            keep &= d > recipe.iso_ms * 1e-3 * FS
    return idx[keep]


def estimate(snips: np.ndarray, estimator: str) -> np.ndarray:
    if estimator == "mean":
        return snips.mean(axis=0)
    if estimator == "median":
        return np.median(snips, axis=0)
    if estimator == "trimmed":
        from scipy.stats import trim_mean

        return trim_mean(snips, 0.2, axis=0)
    raise ValueError(estimator)


# --- one window on one tetrode, after matching ---------------------------------------------------

def tetrode_window(rec, g: int, a: int, b: int, frames_abs: np.ndarray, labels: np.ndarray, amps: np.ndarray,
                   T: np.ndarray, unit_ids: np.ndarray, admission: dict, recipe: Recipe,
                   pulses, sigma: np.ndarray | None) -> dict:
    """Scores, admission and re-estimation for tetrode ``g`` in frames [a, b).

    ``labels`` index ``T`` / ``unit_ids`` (this tetrode's units). Returns per-spike arrays, the
    accepted mask, the reason codes (0 accepted, 1 snr or relative floor, 2 shape, 3 refractory, 5 in a bad
    period of the tetrode), and the
    re-estimated templates (NaN where a unit was not re-estimated)."""
    x = cp.load_tetrode_traces(rec, g, a, b)
    noise, _ = cp.noise_mad_excised(x)
    if sigma is None:
        sigma = mf_noise(x, T)
    fl = frames_abs - a
    sc = spike_scores(x, fl, labels, amps, T, noise, sigma)
    resid = sc.pop("residual_traces")
    if admission.get("sigma_from") == "resid":
        # the projection noise of each template on what the window's fits leave
        nblk = int(10 * FS)
        starts = np.linspace(0, max(x.shape[0] - nblk, 0), 6).astype(np.int64)
        sig_r = mf_noise(np.concatenate([resid[s0:s0 + nblk] for s0 in starts]), T, block_s=60, n_blocks=1)
        sc["snr"] = sc["snr"] * sigma[labels] / sig_r[labels]
        sigma = sig_r
    reason = np.zeros(fl.size, np.int8)
    k = admission.get("snr_floor")
    if k is not None:
        reason[(reason == 0) & ~(sc["snr"] >= k)] = 1
    if admission.get("amp_rel") is not None:
        reason[(reason == 0) & ~(amps >= admission["amp_rel"])] = 1
    if admission.get("drop_bad"):
        iv = bad_intervals().get(int(g), [])
        if iv:
            st, en = np.array(iv).T
            reason[(reason == 0) & in_intervals(frames_abs, st, en)] = 5
    gate = admission.get("gate")
    if gate is not None:
        kind, thr = gate
        if kind == "resid":
            bad = ~(sc["resid"] <= thr)
        elif kind == "resid_rel":
            # residual energy over the median of the tetrode's fits at snr 4 or more (a correct fit's level)
            ref = sc["snr"] >= 4
            ref = np.nanmedian(sc["resid"][ref]) if ref.any() else np.nanmedian(sc["resid"])
            bad = ~(sc["resid"] <= thr * ref)
        else:
            bad = ~(sc[kind] >= thr)
        reason[(reason == 0) & bad] = 2
    refr = admission.get("refractory_ms")
    score = np.nan_to_num(sc["snr"], nan=-np.inf)

    def refractory(units):
        if refr is None:
            return
        rf = int(refr * 1e-3 * FS)
        for k_u in units:
            idx = np.flatnonzero((labels == k_u) & (reason == 0))
            if idx.size > 1:
                o = idx[np.argsort(fl[idx], kind="stable")]
                keep = cp.resolve_refractory(fl[o], score[o], rf)
                reason[o[~keep]] = 3

    refractory(range(len(unit_ids)))
    driven = []
    pr = admission.get("pulse_rule")
    if pr is not None and pulses is not None:
        on, off = pulses
        m = (on >= a) & (on < b)
        if m.any():
            on, off = on[m] - a, off[m] - a
            near_p = in_intervals(fl, on, off + int(0.1 * FS))
            driven = light_driven(fl[reason == 0], labels[reason == 0], on, off, (b - a) / FS)
            kp, ap, cp_ = pr
            relax = (np.isin(labels, driven) & near_p & (reason != 0) & (reason != 5) & (sc["snr"] >= kp) & (amps >= ap)
                     & (sc["cos_pl"] >= cp_))
            reason[relax] = 0
            refractory(driven)
    acc = reason == 0
    newT = np.full_like(T, np.nan)
    n_used = np.zeros(len(unit_ids), np.int64)
    rng = np.random.default_rng(recipe.seed + a // WIN)
    if recipe.snippet == "peeled":
        # the residual over accepted spikes only: rejected fits go back into the data
        rej = ~acc
        if rej.any():
            add_spikes(resid, fl[rej], labels[rej], amps[rej].astype(np.float32), T, NB)
    for k_u in range(len(unit_ids)):
        if acc[labels == k_u].sum() < recipe.min_spikes:
            continue
        sel = select_spikes(recipe, frames_abs, np.where(acc, labels, -1), amps, k_u, pulses, g)
        sel = sel[(fl[sel] >= NB) & (fl[sel] + NA <= x.shape[0])]
        if sel.size < min(recipe.min_spikes, 50):
            continue
        if sel.size > recipe.max_spikes:
            sel = np.sort(rng.choice(sel, recipe.max_spikes, replace=False))
        if recipe.snippet == "peeled":
            sn = cp.snippets(resid, fl[sel]) + amps[sel, None, None].astype(np.float32) * T[k_u][None]
        else:
            sn = cp.snippets(x, fl[sel])
        est = estimate(sn.astype(np.float64), recipe.estimator)
        if getattr(recipe, "scale", "mean") == "trunc":
            # the admission floors cut the unit's amplitude distribution below; rescale the estimate by
            # the truncated-normal mean of the selected spikes' snr over their plain mean
            tn_u = float(np.sqrt((T[k_u].astype(np.float64) ** 2).sum()))
            c = max(admission.get("snr_floor") or 0.0, admission.get("amp_rel", 0.0) * tn_u / sigma[k_u])
            x_sel = sc["snr"][sel]
            mu = truncated_mean(x_sel, c)
            if np.isfinite(mu) and np.nanmean(x_sel) > 0:
                est = est * float(np.clip(mu / np.nanmean(x_sel), 0.5, 1.2))
        newT[k_u] = (1 - recipe.blend) * T[k_u] + recipe.blend * est
        n_used[k_u] = sel.size
    return {"scores": sc, "reason": reason, "newT": newT, "n_used": n_used, "noise": noise, "sigma": sigma,
            "driven": [int(unit_ids[k_u]) for k_u in driven]}


# --- the carry ----------------------------------------------------------------------------------

@dataclasses.dataclass
class Arm:
    name: str
    per_tetrode: bool = True
    amp_rel: float = 0.8
    snr_floor: float | None = None
    gate: tuple | None = None
    refractory_ms: float | None = None
    recipe: Recipe = dataclasses.field(default_factory=Recipe)
    masked: bool = False
    hybrid: bool = True
    method: str = "circus-omp"
    whiten: bool = False
    reseed_every: int = 0
    reseed_sorts: str | None = None
    sigma_from: str = "raw"
    match_floor: tuple | None = None
    pulse_rule: tuple | None = None
    drop_bad: bool = False

    def to_json(self) -> dict:
        d = dataclasses.asdict(self)
        return d


def open_recording(arm: Arm):
    rec = si.load(HYBRID if arm.hybrid else BINARY)
    if arm.masked:
        from ttm.bad_periods import mask_bad_periods, read_table

        rec = mask_bad_periods(rec, read_table(BAD_PERIODS))
    if arm.whiten:
        from spikeinterface.preprocessing import whiten

        rec = whiten(rec, W=whitening_matrix(arm.hybrid), M=None, apply_mean=False, dtype="float32")
    return rec


def whitening_matrix(hybrid: bool = True) -> np.ndarray:
    """Block-diagonal (64, 64) ZCA whitening, one 4x4 block per tetrode, from the covariance of the
    spike-excised samples of the first 300 s (cached in ``whiten_W[_orig].npz``). Applied as ``x @ W``,
    it leaves each tetrode's noise decorrelated across wires with unit variance."""
    path = BAKE / ("whiten_W.npz" if hybrid else "whiten_W_orig.npz")
    if path.exists():
        return np.load(path)["W"]
    rec = si.load(HYBRID if hybrid else BINARY)
    groups = np.asarray(rec.get_property("group")).astype(int)
    W = np.zeros((len(groups), len(groups)), np.float32)
    for g in np.unique(groups):
        x = cp.load_tetrode_traces(rec, int(g), 0, int(300 * FS)).astype(np.float64)
        noise, _ = cp.noise_mad_excised(x)
        keep = cp.excision_mask(x, noise, detect_mad=4.0, pad_frames=30)
        ev, V = np.linalg.eigh(np.cov(x[keep].T))
        ch = np.flatnonzero(groups == g)
        W[np.ix_(ch, ch)] = V @ np.diag(ev ** -0.5) @ V.T
    np.savez(path, W=W)
    return W


def whiten_bank(bank: Bank, W: np.ndarray, groups: np.ndarray) -> Bank:
    out = bank.copy()
    for k, g in enumerate(bank.group):
        ch = np.flatnonzero(groups == g)
        out.T[k] = bank.T[k] @ W[np.ix_(ch, ch)]
    return out


def amplitude_bounds(bank: Bank, sigma: np.ndarray, arm: Arm) -> np.ndarray | list:
    a_rel, k = getattr(arm, "match_floor", None) or (arm.amp_rel, arm.snr_floor)
    if k is None:
        return [a_rel, float("inf")]
    tn = np.sqrt((bank.T.astype(np.float64) ** 2).sum(axis=(1, 2)))
    # with the noise taken from the residual the binding floor is applied after matching
    k = k * (0.7 if arm.sigma_from == "resid" else 1.0)
    lo = np.maximum(a_rel, k * sigma / tn)
    return np.stack([lo, np.full_like(lo, np.inf)], axis=1)


def run_carry(arm: Arm, bank0: Bank, out_dir: pathlib.Path, *, windows: range, n_jobs: int = 16,
              workers: int = 8) -> None:
    """Carry ``bank0`` over ``windows`` with ``arm``; one npz of spikes and one of templates per window
    in ``out_dir`` (a window whose spikes file exists is reused, so a run can resume)."""
    from joblib import Parallel, delayed

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "arm.json").write_text(json.dumps(arm.to_json(), indent=2))
    rec = open_recording(arm)
    pulses = laser_pulses()
    bank = bank0.copy()
    for w in windows:
        f_sp, f_bk = out_dir / f"spikes_w{w:02d}.npz", out_dir / f"bank_w{w:02d}.npz"
        if f_sp.exists() and f_bk.exists():
            bank = Bank.load(f_bk)
            continue
        t0 = time.time()
        a, b = w * WIN, min((w + 1) * WIN, rec.get_num_frames())
        win = rec.frame_slice(a, b)
        groups = sorted(set(bank.group.tolist()))
        # matched-filter noise per template, measured on this window before matching
        sig = Parallel(n_jobs=workers)(
            delayed(_sigma)(rec, g, a, b, bank.T[bank.group == g]) for g in groups)
        sigma = np.zeros(len(bank.ids))
        for g, s in zip(groups, sig):
            sigma[bank.group == g] = s
        sp = match(win, bank.templates(rec), per_tetrode=arm.per_tetrode,
                   amplitudes=amplitude_bounds(bank, sigma, arm), n_jobs=n_jobs, method=arm.method)
        t_match = time.time() - t0
        fr = sp["sample_index"].astype(np.int64) + a
        ci = sp["cluster_index"].astype(np.int64)
        amp = sp["amplitude"].astype(np.float64)
        res = Parallel(n_jobs=workers)(
            delayed(_tetrode_job)(arm, g, a, b, fr[bank.group[ci] == g], ci[bank.group[ci] == g], amp[bank.group[ci] == g],
                                  bank, pulses, sigma[bank.group == g])
            for g in groups)
        cols = {k: [] for k in ("frame", "unit", "a", "snr", "cos_raw", "cos_pl", "resid", "trough", "reason")}
        newbank = bank.copy()
        n_used = np.zeros(len(bank.ids), np.int64)
        noise = {}
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
        arrays = {k: np.concatenate(v) for k, v in cols.items()}
        o = np.argsort(arrays["frame"], kind="stable")
        np.savez(f_sp, **{k: v[o] for k, v in arrays.items()})
        bank.save(f_bk.with_name(f"bankin_w{w:02d}.npz"), sigma=sigma)
        driven = np.array(sorted(u for r in res for u in r.get("driven", [])), np.int64)
        newbank.save(f_bk, n_used=n_used, noise=np.stack([noise[g] for g in groups]), noise_groups=np.array(groups),
                     driven=driven)
        acc = arrays["reason"] == 0
        log(f"{arm.name} w{w:02d}: {acc.sum():,} accepted of {acc.size:,} fits; driven {driven.tolist()}; "
            f"{(n_used > 0).sum()}/{len(n_used)} templates re-estimated; match {t_match / 60:.1f} min, "
            f"total {(time.time() - t0) / 60:.1f} min")
        if arm.reseed_every and w > 0 and w % arm.reseed_every == 0:
            newbank = reseed(rec, newbank, a, b, pathlib.Path(arm.reseed_sorts) / f"reseed_w{w}", n_jobs=n_jobs)
            newbank.save(f_bk, n_used=np.r_[n_used, np.zeros(len(newbank.ids) - len(n_used), np.int64)],
                         noise=np.stack([noise[g] for g in groups]), noise_groups=np.array(groups))
        bank = newbank


def reseed(rec, bank: Bank, a: int, b: int, sort_dir: pathlib.Path, *, n_jobs: int, add_cos: float = 0.8,
           min_snr: float = 5.0, min_spikes: int = 100) -> Bank:
    """June's re-seed (``_mp_common.windowed_carry_forward_reseed``): MountainSort5 on window [a, b)
    (read from ``sort_dir`` when it holds a sort of that window), and every cluster with snr 5 or more
    and 100 spikes or more whose template has shift-cosine under ``add_cos`` to every template of its
    tetrode joins the bank (ids from 2000 up)."""
    import sys

    import spikeinterface.sorters as ss

    sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "spike-sorting"))
    from _mp_common import build_templates_object
    from ttm.tracking import cosine_from_templates, sort_chunk, to_int_numpy_sorting

    win = rec.frame_slice(a, b)
    win.reset_times()
    if (sort_dir / "by_group").exists():
        parts = []
        for g in sorted(int(p.name) for p in (sort_dir / "by_group").iterdir()):
            sg = ss.read_sorter_folder(sort_dir / "by_group" / str(g), register_recording=False)
            sg.set_property("group", np.full(sg.get_num_units(), g))
            parts.append(sg)
        rs = to_int_numpy_sorting(si.aggregate_units(parts))
    else:
        rs = to_int_numpy_sorting(sort_chunk(win, sort_dir))
    _, raz = build_templates_object(rs, win, with_snr=True, n_jobs=n_jobs)
    snr = raz.get_extension("quality_metrics").get_data()["snr"].to_numpy()
    from spikeinterface.core.template_tools import get_dense_templates_array

    dense = get_dense_templates_array(raz, return_in_uV=False)
    groups = np.asarray(rec.get_property("group")).astype(int)
    rgrp = np.asarray(rs.get_property("group")).astype(int)
    out = bank.copy()
    next_id = max(2000, int(bank.ids.max()) + 1)
    added = []
    for j, cu in enumerate(rs.unit_ids):
        if snr[j] < min_snr or len(rs.get_unit_spike_train(cu)) < min_spikes:
            continue
        g = rgrp[j]
        t = dense[j][:, groups == g]
        best = max((cosine_from_templates(t, bt, max_shift_samples=10) for bt in out.T[out.group == g]),
                   default=-1.0)
        if best < add_cos:
            out = Bank(np.r_[out.ids, next_id], np.r_[out.group, g], np.concatenate([out.T, t[None].astype(np.float32)]))
            added.append((next_id, int(g), round(float(best), 3)))
            next_id += 1
    log(f"reseed after window {a // WIN}: +{len(added)} units {added}")
    return out


def _sigma(rec, g, a, b, T, *, block_s: float = 10.0, n_blocks: int = 6):
    nblk = int(block_s * FS)
    starts = np.linspace(a, max(b - nblk, a), n_blocks).astype(np.int64)
    x = np.concatenate([cp.load_tetrode_traces(rec, g, s, s + nblk) for s in starts])
    return mf_noise(x, T, block_s=block_s * n_blocks, n_blocks=1)


def _tetrode_job(arm, g, a, b, frames_abs, ci, amps, bank, pulses, sigma):
    rec = open_recording(arm)
    gi = np.flatnonzero(bank.group == g)
    local = np.searchsorted(gi, ci)
    o = np.argsort(frames_abs, kind="stable")
    inv = np.empty_like(o)
    inv[o] = np.arange(o.size)
    adm = {"snr_floor": arm.snr_floor, "gate": arm.gate, "refractory_ms": arm.refractory_ms,
           "sigma_from": getattr(arm, "sigma_from", "raw"), "amp_rel": arm.amp_rel,
           "pulse_rule": getattr(arm, "pulse_rule", None), "drop_bad": getattr(arm, "drop_bad", False)}
    r = tetrode_window(rec, g, a, b, frames_abs[o], local[o], amps[o], bank.T[gi], bank.ids[gi], adm,
                       arm.recipe, pulses, sigma)
    r["scores"] = {k: v[inv] for k, v in r["scores"].items()}
    r["reason"] = r["reason"][inv]
    return r
