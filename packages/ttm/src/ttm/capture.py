"""Shape-based spike capture and purity report for tetrode sortings.

The sorter decides which events exist by amplitude (a detection threshold of
5.5 MAD on the whitened signal) and which cluster they join by PCA features
that weight amplitude heavily. Two things fall through that net on the
TTM-003 tetrodes: a cell whose spikes shrink to half size inside a burst
(TT6) and a cell whose spikes sit at about 4.3 MAD (TT10). This module
replaces the two decisions for one tetrode at a time:

1. ``noise_mad_excised``: a per-wire noise estimate in which the samples around
   every threshold crossing are left out of the median, so a busy wire does not
   inflate its own threshold (``OPEN_CONCERNS.md``, concern 2).
2. ``detect_troughs``: every negative trough below ``threshold_mad`` times that
   noise on any wire, one event per trough.
3. ``cosine_to_templates``: the cosine of each event's four-wire snippet over
   the full template window (-1.0 to +2.0 ms) to every template of the same
   tetrode, with the amplitude left free.
4. ``assign_events``: each event goes to the template it matches best when that
   cosine clears ``gate_r``; events that would break the host's refractory
   period lose to the better-matching event of the pair; the rest form a
   per-tetrode multi-unit pool.
5. ``match_wobble``: the alternative in which SpikeInterface's wobble matcher
   both detects and assigns (templates scaled to each spike), gated by the
   same cosine, so the two can be compared on the same traces.
6. ``unit_purity`` and ``plot_unit_purity``: refractory contamination, the
   amplitude-split cross-correlogram, cross-correlograms with every same-tetrode
   neighbour with the dip-or-filled verdict, within-cluster waveform overlays
   and the amplitude-by-shape distribution, so purity is judged by pictures
   as well as numbers.

Everything works on a ``(n_samples, 4)`` float32 array of one tetrode's
preprocessed traces (bandpass and common median reference applied, raw ADC
units) plus the recording's gain; nothing here reads a store or a manifest.
Templates are ``(n_units, n_samp, 4)`` arrays aligned so the four-wire trough
sits at ``nbefore``.
"""

from __future__ import annotations

import dataclasses
import math

import numpy as np

FS = 30000.0
NBEFORE = 30   # 1.0 ms
NAFTER = 60    # 2.0 ms
N_SAMP = NBEFORE + NAFTER


# --- traces ---------------------------------------------------------------------------------

def tetrode_channel_indices(recording, group: int) -> np.ndarray:
    """Channel indices of ``group`` (0-based tetrode index) in ``recording``."""
    groups = np.asarray(recording.get_property("group")).astype(int)
    return np.flatnonzero(groups == group)


def load_tetrode_traces(recording, group: int, start_frame: int | None = None,
                        end_frame: int | None = None) -> np.ndarray:
    """``(n_samples, 4)`` float32 traces of one tetrode (raw units, as stored)."""
    idx = tetrode_channel_indices(recording, group)
    ids = np.asarray(recording.channel_ids)[idx]
    x = recording.get_traces(start_frame=start_frame, end_frame=end_frame, channel_ids=list(ids))
    return np.ascontiguousarray(x, dtype=np.float32)


# --- noise ----------------------------------------------------------------------------------

def noise_mad(traces: np.ndarray, *, subsample: int = 10) -> np.ndarray:
    """Plain per-wire noise: median absolute deviation divided by 0.6745. The median is taken
    over every ``subsample``-th sample; on a 32 min TTM-003 tetrode the 1-in-10 estimate differs
    from the full one by under 0.15 %."""
    return np.median(np.abs(traces[::subsample]), axis=0) / 0.6745


def excision_mask(traces: np.ndarray, noise: np.ndarray, *, detect_mad: float = 4.0,
                  pad_frames: int = 30) -> np.ndarray:
    """Boolean mask of samples farther than ``pad_frames`` from every crossing below
    ``-detect_mad * noise`` on any wire."""
    below = (traces < -detect_mad * noise).any(axis=1)
    starts = np.flatnonzero(below)
    n = traces.shape[0]
    if starts.size == 0:
        return np.ones(n, dtype=bool)
    edges = np.bincount(np.clip(starts - pad_frames, 0, n), minlength=n + 1)
    edges -= np.bincount(np.clip(starts + pad_frames + 1, 0, n), minlength=n + 1)
    return np.cumsum(edges[:-1]) == 0


def noise_mad_excised(traces: np.ndarray, *, detect_mad: float = 4.0, pad_frames: int = 30,
                      n_iter: int = 2, subsample: int = 10) -> tuple[np.ndarray, float]:
    """Per-wire noise with spikes excised: samples within ``pad_frames`` of any crossing below
    ``-detect_mad`` times the running noise estimate, on any wire, are dropped from the median.
    Repeated ``n_iter`` times from the plain estimate. Returns (noise, fraction of samples excised).
    """
    noise = noise_mad(traces, subsample=subsample)
    keep = np.ones(traces.shape[0], dtype=bool)
    for _ in range(n_iter):
        keep = excision_mask(traces, noise, detect_mad=detect_mad, pad_frames=pad_frames)
        if keep.sum() < 1000:
            break
        noise = np.median(np.abs(traces[keep][::subsample]), axis=0) / 0.6745
    return noise, float(1.0 - keep.mean())


def noise_by_block(traces: np.ndarray, *, block_frames: int, detect_mad: float = 4.0,
                   threshold_mad: float = 3.5, high_mad: float = 7.0) -> dict[str, np.ndarray]:
    """Plain and excised noise per block plus two crossing rates per block, per wire, so the
    question "does the noise estimate track the firing rate" can be answered on one channel over
    time. Returns ``{"plain": (n_blocks, 4), "excised": (n_blocks, 4), "crossings_per_s": (n_blocks, 4),
    "high_crossings_per_s": (n_blocks, 4)}``. Crossings are counted per wire against the
    session-wide excised noise times ``threshold_mad`` (near the detection floor, where noise
    crossings contribute) and times ``high_mad`` (spikes only), both fixed across blocks.
    """
    n = traces.shape[0] // block_frames
    ref, _ = noise_mad_excised(traces, detect_mad=detect_mad)
    plain = np.zeros((n, traces.shape[1]))
    exc = np.zeros((n, traces.shape[1]))
    rate = np.zeros((n, traces.shape[1]))
    high = np.zeros((n, traces.shape[1]))
    for b in range(n):
        x = traces[b * block_frames:(b + 1) * block_frames]
        plain[b] = noise_mad(x)
        exc[b] = noise_mad_excised(x, detect_mad=detect_mad)[0]
        below = x < -threshold_mad * ref
        rate[b] = (below[1:] & ~below[:-1]).sum(axis=0) / (block_frames / FS)
        below = x < -high_mad * ref
        high[b] = (below[1:] & ~below[:-1]).sum(axis=0) / (block_frames / FS)
    return {"plain": plain, "excised": exc, "crossings_per_s": rate, "high_crossings_per_s": high}


# --- detection ------------------------------------------------------------------------------

def detect_troughs(traces: np.ndarray, noise: np.ndarray, *, threshold_mad: float = 3.5,
                   dead_frames: int = 9, nbefore: int = NBEFORE, nafter: int = NAFTER
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Negative troughs below ``threshold_mad`` times ``noise`` on any wire.

    An event is a sample at which the four-wire minimum (each wire in units of its own noise) is
    below the threshold and is the minimum within ``dead_frames`` on either side, so two troughs
    closer than ``dead_frames`` (0.3 ms by default) count once. Events whose template window would
    leave the array are dropped. Returns (frames, amplitude in MAD units of the deepest wire).
    """
    from scipy.ndimage import minimum_filter1d

    m = (traces / noise[None, :]).min(axis=1)
    local = minimum_filter1d(m, 2 * dead_frames + 1, mode="nearest")
    frames = np.flatnonzero((m < -threshold_mad) & (m <= local))
    if frames.size > 1:  # a flat minimum can fire twice inside one dead window
        frames = frames[np.r_[True, np.diff(frames) > dead_frames]]
    ok = (frames >= nbefore) & (frames + nafter <= traces.shape[0])
    frames = frames[ok]
    return frames, -m[frames]


def snippets(traces: np.ndarray, frames: np.ndarray, *, nbefore: int = NBEFORE,
             nafter: int = NAFTER) -> np.ndarray:
    """``(n, nbefore + nafter, 4)`` windows around ``frames`` (all must be in bounds)."""
    cols = np.arange(-nbefore, nafter)
    return traces[frames[:, None] + cols[None, :], :]


# --- templates ------------------------------------------------------------------------------

def templates_from_trains(traces: np.ndarray, trains: dict, *, nbefore: int = NBEFORE,
                          nafter: int = NAFTER, max_spikes: int = 1000, seed: int = 0,
                          align: bool = True) -> tuple[np.ndarray, list]:
    """Mean four-wire template per unit from up to ``max_spikes`` random spikes of each train
    (frames into ``traces``). With ``align`` the template is shifted so its four-wire trough sits at
    ``nbefore`` (sorters align to their own peak convention). Returns (templates, unit ids).
    """
    rng = np.random.default_rng(seed)
    ids = sorted(trains)
    out = np.zeros((len(ids), nbefore + nafter, traces.shape[1]), np.float32)
    pad = 10
    for i, u in enumerate(ids):
        fr = np.asarray(trains[u], np.int64)
        fr = fr[(fr >= nbefore + pad) & (fr + nafter + pad <= traces.shape[0])]
        if fr.size > max_spikes:
            fr = np.sort(rng.choice(fr, max_spikes, replace=False))
        if fr.size == 0:
            continue
        wide = snippets(traces, fr, nbefore=nbefore + pad, nafter=nafter + pad).mean(axis=0)
        shift = 0
        if align:
            k = int(np.argmin(wide[pad:pad + nbefore + nafter].min(axis=1)))
            shift = int(np.clip(k - nbefore, -pad, pad))
        out[i] = wide[pad + shift:pad + shift + nbefore + nafter]
    return out, ids


def template_cosine(a: np.ndarray, b: np.ndarray, *, max_shift: int = 3) -> float:
    """Cosine between two ``(n_samp, 4)`` templates, maximised over integer shifts."""
    best = -1.0
    n = a.shape[0]
    for s in range(-max_shift, max_shift + 1):
        if s >= 0:
            x, y = a[s:], b[:n - s]
        else:
            x, y = a[:n + s], b[-s:]
        d = np.linalg.norm(x) * np.linalg.norm(y)
        if d > 0:
            best = max(best, float((x * y).sum() / d))
    return best


# --- cosine and assignment ------------------------------------------------------------------

def cosine_to_templates(traces: np.ndarray, frames: np.ndarray, templates: np.ndarray, *,
                        nbefore: int = NBEFORE, s_max: int = 1, batch: int = 50_000) -> np.ndarray:
    """``(n_events, n_templates)`` cosine of each event's full-window snippet to every template,
    the maximum over shifts of ``-s_max..s_max`` samples. NaN where the shifted window leaves the array.
    """
    n_samp = templates.shape[1]
    nafter = n_samp - nbefore
    tnorm = np.sqrt((templates ** 2).sum(axis=(1, 2)))
    tflat = templates.reshape(templates.shape[0], -1)
    out = np.full((frames.size, templates.shape[0]), np.nan)
    valid = (frames - nbefore - s_max >= 0) & (frames + nafter + s_max <= traces.shape[0])
    idx = np.flatnonzero(valid)
    for b0 in range(0, idx.size, batch):
        sel = idx[b0:b0 + batch]
        best = np.full((sel.size, templates.shape[0]), -np.inf)
        for s in range(-s_max, s_max + 1):
            sn = snippets(traces, frames[sel] + s, nbefore=nbefore, nafter=nafter)
            flat = sn.reshape(sel.size, -1)
            snorm = np.sqrt((flat ** 2).sum(axis=1))
            r = (flat @ tflat.T) / (snorm[:, None] * tnorm[None, :])
            best = np.maximum(best, r)
        out[sel] = best
    return out


def noise_corrected_cosine(traces: np.ndarray, frames: np.ndarray, templates: np.ndarray,
                           noise: np.ndarray, *, nbefore: int = NBEFORE, s_max: int = 1,
                           floor: float = 0.25, batch: int = 50_000) -> np.ndarray:
    """Cosine with the noise energy removed from the snippet norm.

    A small spike of exactly the template's shape has a raw cosine well below one because the
    snippet's energy is mostly noise: with template energy ``E`` at the spike's size and noise
    energy ``N = n_samp * sum(noise ** 2)``, the raw cosine is about ``sqrt(E / (E + N))`` (0.75 for
    a TT6 spike at half size, 0.6 at a third). Dividing by ``sqrt(max(|s|^2 - N, floor * |s|^2))``
    instead of ``|s|`` estimates the cosine the snippet's signal part has to the template, so the
    gate can ask for shape agreement at any amplitude. Its noise is about ``1 / sqrt(E)`` (in MAD
    units), so it is only informative for events several MAD above the noise.
    """
    n_samp = templates.shape[1]
    nafter = n_samp - nbefore
    N = n_samp * float((np.asarray(noise, np.float64) ** 2).sum())
    tnorm = np.sqrt((templates ** 2).sum(axis=(1, 2)))
    tflat = templates.reshape(templates.shape[0], -1)
    out = np.full((frames.size, templates.shape[0]), np.nan)
    valid = (frames - nbefore - s_max >= 0) & (frames + nafter + s_max <= traces.shape[0])
    idx = np.flatnonzero(valid)
    for b0 in range(0, idx.size, batch):
        sel = idx[b0:b0 + batch]
        best = np.full((sel.size, templates.shape[0]), -np.inf)
        for s in range(-s_max, s_max + 1):
            sn = snippets(traces, frames[sel] + s, nbefore=nbefore, nafter=nafter)
            flat = sn.reshape(sel.size, -1)
            ssq = (flat.astype(np.float64) ** 2).sum(axis=1)
            signal = np.sqrt(np.maximum(ssq - N, floor * ssq))
            r = (flat @ tflat.T) / (signal[:, None] * tnorm[None, :])
            best = np.maximum(best, r)
        out[sel] = best
    return out


def resolve_refractory(frames: np.ndarray, score: np.ndarray, refr_frames: int) -> np.ndarray:
    """Boolean mask of events that survive the refractory rule: events closer than ``refr_frames``
    form a chain and only the chain's highest ``score`` survives. ``frames`` must be sorted.
    """
    keep = np.ones(frames.size, dtype=bool)
    if frames.size < 2:
        return keep
    close = np.diff(frames) < refr_frames
    chain_id = np.r_[0, np.cumsum(~close)]
    for c in np.flatnonzero(np.bincount(chain_id) > 1):
        members = np.flatnonzero(chain_id == c)
        keep[members] = False
        keep[members[np.argmax(score[members])]] = True
    return keep


@dataclasses.dataclass
class Assignment:
    """Per-event outcome of ``assign_events``. ``label`` is the template index or -1; ``reason`` is
    ``accepted``, ``below_gate`` (best cosine under the gate) or ``refractory`` (lost a chain)."""
    frames: np.ndarray
    amp_mad: np.ndarray
    label: np.ndarray
    best: np.ndarray
    r_best: np.ndarray
    reason: np.ndarray

    def train(self, template_index: int) -> np.ndarray:
        return self.frames[self.label == template_index]

    def trains(self, n_templates: int) -> dict[int, np.ndarray]:
        return {i: self.train(i) for i in range(n_templates)}

    @property
    def pool(self) -> np.ndarray:
        """Frames left for the multi-unit pool (every event not accepted by a template)."""
        return self.frames[self.label < 0]


def assign_events(frames: np.ndarray, amp_mad: np.ndarray, cosines: np.ndarray, *,
                  gate_r: float = 0.6, refr_frames: int = 45, gate: np.ndarray | None = None,
                  gate_value: float = 0.8) -> Assignment:
    """Competitive shape assignment: each event goes to its best-cosine template if that cosine
    reaches ``gate_r``; among the events a template accepts, those closer than ``refr_frames``
    (1.5 ms at 30 kHz by default) keep only the better match. Events must be sorted by frame.
    ``gate`` is an optional second ``(n_events, n_templates)`` score (the noise-corrected cosine)
    that the winning template must also reach (``gate_value``).
    """
    n = frames.size
    label = np.full(n, -1, np.int64)
    reason = np.full(n, "below_gate", dtype=object)
    if cosines.shape[1] == 0 or n == 0:
        return Assignment(frames, amp_mad, label, np.full(n, -1), np.full(n, np.nan), reason)
    r_best = np.nanmax(np.where(np.isnan(cosines), -np.inf, cosines), axis=1)
    r_best = np.where(np.isfinite(r_best), r_best, np.nan)
    best = np.where(np.isnan(r_best), -1, np.nanargmax(np.where(np.isnan(cosines), -np.inf, cosines), axis=1))
    cand = np.isfinite(r_best) & (r_best >= gate_r)
    if gate is not None:
        cand &= np.isfinite(gate[np.arange(n), np.clip(best, 0, None)]) & (
            gate[np.arange(n), np.clip(best, 0, None)] >= gate_value)
    for u in np.unique(best[cand]):
        idx = np.flatnonzero(cand & (best == u))
        keep = resolve_refractory(frames[idx], r_best[idx], refr_frames)
        label[idx[keep]] = u
        reason[idx[keep]] = "accepted"
        reason[idx[~keep]] = "refractory"
    return Assignment(frames, amp_mad, label, best, r_best, reason)


def score_events(traces: np.ndarray, frames: np.ndarray, templates: np.ndarray, noise: np.ndarray, *,
                 nbefore: int = NBEFORE, s_max: int = 1, floor: float = 0.25, batch: int = 50_000) -> dict:
    """Raw cosine of every event to every template (max over shifts), and for the best template the
    shift, the amplitude scale ``a = conv / |t|^2`` and the noise-corrected cosine. One pass over the
    snippets serves ``cosine_to_templates`` and ``noise_corrected_cosine`` at once."""
    n_samp = templates.shape[1]
    nafter = n_samp - nbefore
    N = n_samp * float((np.asarray(noise, np.float64) ** 2).sum())
    tsq = (templates.astype(np.float64) ** 2).sum(axis=(1, 2))
    tnorm = np.sqrt(tsq)
    tflat = templates.reshape(templates.shape[0], -1)
    n, k = frames.size, templates.shape[0]
    r = np.full((n, k), np.nan)
    shift = np.zeros(n, np.int64)
    a = np.full(n, np.nan)
    rc = np.full(n, np.nan)
    valid = (frames - nbefore - s_max >= 0) & (frames + nafter + s_max <= traces.shape[0])
    idx = np.flatnonzero(valid)
    for b0 in range(0, idx.size, batch):
        sel = idx[b0:b0 + batch]
        best_r = np.full((sel.size, k), -np.inf)
        best_conv = np.zeros((sel.size, k))
        best_shift = np.zeros((sel.size, k), np.int64)
        best_ssq = np.zeros((sel.size, k))
        for sh in range(-s_max, s_max + 1):
            sn = snippets(traces, frames[sel] + sh, nbefore=nbefore, nafter=nafter)
            flat = sn.reshape(sel.size, -1)
            ssq = (flat.astype(np.float64) ** 2).sum(axis=1)
            conv = (flat @ tflat.T).astype(np.float64)
            rr = conv / (np.sqrt(ssq)[:, None] * tnorm[None, :])
            better = rr > best_r
            best_r = np.where(better, rr, best_r)
            best_conv = np.where(better, conv, best_conv)
            best_shift = np.where(better, sh, best_shift)
            best_ssq = np.where(better, ssq[:, None], best_ssq)
        r[sel] = best_r
        j = np.argmax(best_r, axis=1)
        rows = np.arange(sel.size)
        shift[sel] = best_shift[rows, j]
        a[sel] = best_conv[rows, j] / tsq[j]
        ssq_j = best_ssq[rows, j]
        signal = np.sqrt(np.maximum(ssq_j - N, floor * ssq_j))
        rc[sel] = best_conv[rows, j] / (signal * tnorm[j])
    return {"r": r, "shift": shift, "a": a, "rc": rc}


def subtract_templates(residual: np.ndarray, frames: np.ndarray, labels: np.ndarray, a: np.ndarray,
                       templates: np.ndarray, *, nbefore: int = NBEFORE) -> None:
    """Subtract ``a[i] * templates[labels[i]]`` at ``frames[i]`` from ``residual`` in place."""
    n_samp = templates.shape[1]
    for f, u, s in zip(frames, labels, a):
        residual[f - nbefore:f - nbefore + n_samp] -= np.float32(s) * templates[u]


@dataclasses.dataclass
class Capture:
    """Result of ``capture_tetrode``: every detected event with its outcome.

    ``label`` is the template index or -1; ``reason`` is ``accepted``, ``below_gate`` or
    ``refractory``; ``r`` is the raw cosine to the assigned (or best) template, ``rc`` the
    noise-corrected one, ``a`` the amplitude scale relative to the template, ``pass_index`` the
    peeling pass that accepted the event (0 for the first) and ``amp_mad`` the trough depth in MAD.
    """
    frames: np.ndarray
    amp_mad: np.ndarray
    label: np.ndarray
    best: np.ndarray
    r: np.ndarray
    rc: np.ndarray
    a: np.ndarray
    reason: np.ndarray
    pass_index: np.ndarray
    noise: np.ndarray
    noise_plain: np.ndarray
    excised_fraction: float
    template_ids: list

    def train(self, i: int) -> np.ndarray:
        return self.frames[self.label == i]

    def trains(self) -> dict:
        return {i: self.train(i) for i in range(len(self.template_ids))}

    def amps(self) -> dict:
        return {i: self.amp_mad[self.label == i] for i in range(len(self.template_ids))}

    def rs(self) -> dict:
        return {i: self.r[self.label == i] for i in range(len(self.template_ids))}

    def rejected_by_refractory(self) -> dict:
        return {i: int(((self.best == i) & (self.reason == "refractory")).sum())
                for i in range(len(self.template_ids))}

    @property
    def pool(self) -> np.ndarray:
        return self.frames[self.label < 0]


def capture_tetrode(traces: np.ndarray, templates: np.ndarray, template_ids: list, *,
                    threshold_mad: float = 3.5, gate_r: float = 0.6, gate_rc: float | None = None,
                    refr_frames: int = 45, n_pass: int = 3, dead_frames: int = 9,
                    nbefore: int = NBEFORE, s_max: int = 1, noise: np.ndarray | None = None,
                    detect_mad: float = 4.0, margin: float = 0.0) -> Capture:
    """Detect every trough at ``threshold_mad`` times the excised noise and assign it by shape.

    Assignment runs in ``n_pass`` peeling passes: an event accepted in one pass has its scaled
    template subtracted from a working copy of the traces, so a smaller spike inside its 3 ms window
    is scored on a cleaner snippet in the next pass (a spike 1 to 2 ms from a bigger one otherwise
    fails the gate and the pair is lost, which hollows every cross-correlogram). ``gate_r`` is the
    raw cosine the best template must reach; ``gate_rc`` (optional) is the noise-corrected cosine it
    must also reach; ``margin`` is how much the best template's cosine must exceed the runner-up's
    (0 accepts every winner). After the passes each template keeps, among its accepted events closer
    than ``refr_frames``, only the best match.
    """
    noise_plain = noise_mad(traces)
    if noise is None:
        noise, frac = noise_mad_excised(traces, detect_mad=detect_mad)
    else:
        frac = float("nan")
    frames, amp = detect_troughs(traces, noise, threshold_mad=threshold_mad, dead_frames=dead_frames,
                                 nbefore=nbefore, nafter=templates.shape[1] - nbefore)
    n = frames.size
    label = np.full(n, -1, np.int64)
    best = np.full(n, -1, np.int64)
    r = np.full(n, np.nan)
    rc = np.full(n, np.nan)
    a = np.full(n, np.nan)
    reason = np.full(n, "below_gate", dtype=object)
    pass_index = np.full(n, -1, np.int64)
    if templates.shape[0] == 0 or n == 0:
        return Capture(frames, amp, label, best, r, rc, a, reason, pass_index, noise, noise_plain, frac, list(template_ids))
    residual = traces.copy()
    pending = np.arange(n)
    for p in range(n_pass):
        if pending.size == 0:
            break
        sc = score_events(residual, frames[pending], templates, noise, nbefore=nbefore, s_max=s_max)
        rr = np.where(np.isnan(sc["r"]), -np.inf, sc["r"])
        j = np.argmax(rr, axis=1)
        r_best = rr[np.arange(pending.size), j]
        ok = np.isfinite(r_best) & (r_best >= gate_r)
        if gate_rc is not None:
            ok &= np.isfinite(sc["rc"]) & (sc["rc"] >= gate_rc)
        if margin > 0 and rr.shape[1] > 1:
            second = np.partition(rr, -2, axis=1)[:, -2]
            ok &= (r_best - second) >= margin
        best[pending] = np.where(np.isfinite(r_best), j, -1)
        r[pending] = np.where(np.isfinite(r_best), r_best, np.nan)
        rc[pending] = sc["rc"]
        a[pending] = sc["a"]
        acc = pending[ok]
        label[acc] = j[ok]
        reason[acc] = "accepted"
        pass_index[acc] = p
        subtract_templates(residual, frames[acc] + sc["shift"][ok], j[ok], sc["a"][ok], templates, nbefore=nbefore)
        pending = pending[~ok]
    for u in np.unique(label[label >= 0]):
        idx = np.flatnonzero(label == u)
        keep = resolve_refractory(frames[idx], r[idx], refr_frames)
        label[idx[~keep]] = -1
        reason[idx[~keep]] = "refractory"
    return Capture(frames, amp, label, best, r, rc, a, reason, pass_index, noise, noise_plain, frac, list(template_ids))


# --- wobble arm ------------------------------------------------------------------------------

def match_wobble(recording_tetrode, templates: np.ndarray, *, nbefore: int = NBEFORE,
                 threshold_factor: float = 0.35, n_jobs: int = 8,
                 refractory_frames: int = 10) -> tuple[np.ndarray, np.ndarray]:
    """Run SpikeInterface's wobble matcher on a four-channel recording with the given templates
    (raw units, aligned at ``nbefore``). The objective threshold is ``threshold_factor`` times the
    median template energy, the permissive admit the June work settled on; the cosine gate is applied
    by the caller with ``cosine_to_templates``. Returns (frames, template index) sorted by frame.
    """
    import spikeinterface as si
    from spikeinterface.core import Templates
    from spikeinterface.sortingcomponents.matching import find_spikes_from_templates

    n_units = templates.shape[0]
    mask = np.ones((n_units, 4), dtype=bool)
    bank = Templates(templates_array=np.asarray(templates, np.float32), sampling_frequency=FS,
                     nbefore=nbefore, is_in_uV=False, sparsity_mask=mask,
                     channel_ids=np.asarray(recording_tetrode.channel_ids),
                     unit_ids=np.arange(n_units), probe=None, check_for_consistent_sparsity=True)
    tsq = np.median((templates ** 2).sum(axis=(1, 2)))
    kwargs = {"parameters": {"threshold": float(threshold_factor * tsq), "approx_rank": 4,
                             "amplitude_variance": 1.0, "jitter_factor": 8, "max_iter": 1000,
                             "refractory_period_frames": int(refractory_frames),
                             "visibility_threshold": 1.0, "scale_min": 0.0, "scale_max": float("inf")},
              "engine": "numpy", "torch_device": "cpu", "shared_memory": True}
    spikes = find_spikes_from_templates(
        recording_tetrode, bank, method="wobble", method_kwargs=kwargs,
        job_kwargs={"n_jobs": n_jobs, "chunk_duration": "1s", "progress_bar": False})
    uf = "cluster_index" if "cluster_index" in spikes.dtype.names else "unit_index"
    order = np.argsort(spikes["sample_index"], kind="stable")
    return spikes["sample_index"][order].astype(np.int64), spikes[uf][order].astype(np.int64)


# --- purity ---------------------------------------------------------------------------------

def rp_contamination(train: np.ndarray, n_frames: int, *, refr_frames: int = 30,
                     censored_frames: int = 0) -> tuple[float, int]:
    """Llobet et al. refractory-period contamination (the SpikeInterface ``rp_contamination``
    formula) for a sorted frame train. Returns (contamination, number of violations)."""
    n = train.size
    if n <= 1:
        return float("nan"), 0
    isi = np.diff(np.sort(train))
    nv = int(((isi > censored_frames) & (isi < refr_frames)).sum())
    denom = 1 - nv * (n_frames - 2 * n * censored_frames) / (n ** 2 * (refr_frames - censored_frames))
    return (1.0 if denom < 0 else 1 - math.sqrt(denom)), nv


def ccg_lags(a: np.ndarray, b: np.ndarray, max_lag: int) -> np.ndarray:
    """All lags ``b - a`` within ``max_lag`` frames for two sorted frame trains."""
    if a.size == 0 or b.size == 0:
        return np.zeros(0, np.int64)
    lo = np.searchsorted(b, a - max_lag)
    hi = np.searchsorted(b, a + max_lag)
    parts = [b[s:e] - t for t, s, e in zip(a, lo, hi) if e > s]
    return np.concatenate(parts) if parts else np.zeros(0, np.int64)


def ccg_verdict(a: np.ndarray, b: np.ndarray, *, refr_ms: float = 1.5, flank_ms=(5.0, 25.0),
                min_flank: int = 30, dip: float = 0.30, fill: float = 0.70,
                dead_ms: float = 0.0) -> dict:
    """Cross-correlogram arbiter between two frame trains: the ratio of coincidences per ms within
    ``refr_ms`` to those in the flank band. Below ``dip`` the pair shares a refractory period
    (``duplicate``, one cell split in two); above ``fill`` they fire independently (``distinct``);
    between is ``ambiguous``; too few flank coincidences is ``segregated``. The thresholds are the
    June arbiter's (``investigations/spike-sorting/_assignment_eval.py``). ``dead_ms`` is the
    detector's dead time: lags inside it cannot occur for any pair (one event per trough), so the
    central window is ``(dead_ms, refr_ms]`` and is normalised by its own width.
    """
    max_lag = int(flank_ms[1] * FS / 1000)
    lags = np.abs(ccg_lags(np.sort(a), np.sort(b), max_lag)) / FS * 1000
    central = int(((lags > dead_ms) & (lags <= refr_ms)).sum())
    flank = int(((lags >= flank_ms[0]) & (lags <= flank_ms[1])).sum())
    if flank < min_flank:
        return {"ratio": float("nan"), "central": central, "flank": flank, "verdict": "segregated"}
    ratio = (central / (2 * (refr_ms - dead_ms))) / (flank / (2 * (flank_ms[1] - flank_ms[0])))
    verdict = "duplicate" if ratio < dip else "distinct" if ratio > fill else "ambiguous"
    return {"ratio": float(ratio), "central": central, "flank": flank, "verdict": verdict}


def amplitude_split(train: np.ndarray, amp: np.ndarray, n_frames: int, *, lower_q: float = 0.35,
                    upper_q: float = 0.65, core_ms: float = 1.5, window_ms: float = 20.0,
                    bin_ms: float = 0.5, dead_ms: float = 0.0) -> dict:
    """The large-versus-small test: split a unit's spikes at the ``lower_q`` and ``upper_q``
    amplitude quantiles (``amp`` positive, larger is bigger), cross-correlate the two populations
    and compare the count within ``core_ms`` of zero lag with the independent expectation. One cell
    gives a ratio near zero (its small spikes cannot fall inside its large spikes' refractory
    period); two cells give about one. Returns the ratio, the counts, the lag histogram and the
    2 to 3 ms ratio (burst following).
    """
    if train.size < 50:
        return {"ratio": float("nan"), "n_large": 0, "n_small": 0}
    lo, hi = np.quantile(amp, [lower_q, upper_q])
    large = np.sort(train[amp >= hi])
    small = np.sort(train[amp <= lo])
    max_lag = int(window_ms * FS / 1000)
    lags = ccg_lags(large, small, max_lag) / FS * 1000
    edges = np.arange(-window_ms, window_ms + bin_ms, bin_ms) - bin_ms / 2
    counts, _ = np.histogram(lags, edges)
    centers = edges[:-1] + bin_ms / 2
    expected_per_bin = large.size * small.size * (bin_ms / 1000) / (n_frames / FS)
    core = (np.abs(centers) <= core_ms) & (np.abs(centers) > dead_ms)
    obs = int(counts[core].sum())
    exp = float(expected_per_bin * core.sum())
    burst = (centers > 1.75) & (centers < 3.25)
    return {"ratio": float(obs / exp) if exp > 0 else float("nan"), "observed": obs, "expected": exp,
            "n_large": int(large.size), "n_small": int(small.size), "large_min": float(hi),
            "small_max": float(lo), "lag_ms": centers, "counts": counts,
            "expected_per_bin": expected_per_bin,
            "burst_2_3ms_ratio": float(counts[burst].mean() / expected_per_bin) if expected_per_bin > 0 else float("nan")}


def refractory_conflict_ratio(n_accepted: int, n_rejected: int, n_frames: int, refr_frames: int) -> float:
    """Refractory conflicts a template saw (events it matched best that fell inside another accepted
    event's refractory window) relative to the number two independent trains of the accepted size
    would produce by chance. Near zero: the candidates respected one refractory period. Near one or
    above: the template was collecting from more than one cell. Replaces refractory contamination
    for a train whose refractory period was enforced during assignment."""
    expected = n_accepted ** 2 * 2 * refr_frames / max(n_frames, 1)
    return float(n_rejected / expected) if expected > 0 else float("nan")


def unit_purity(trains: dict, amps: dict, n_frames: int, *, refr_frames: int = 30,
                dead_ms: float = 0.0, rejected: dict | None = None,
                assign_refr_frames: int = 45) -> dict:
    """Per-unit refractory contamination, refractory-conflict ratio (when ``rejected`` gives the
    per-unit count of events the assignment's refractory rule discarded), amplitude-split ratio and
    the CCG verdict against every other unit of the same tetrode (``trains`` and ``amps`` are keyed
    alike; frames sorted)."""
    out = {}
    ids = sorted(trains)
    for u in ids:
        c, nv = rp_contamination(trains[u], n_frames, refr_frames=refr_frames)
        split = amplitude_split(trains[u], amps[u], n_frames, dead_ms=dead_ms)
        neigh = {v: ccg_verdict(trains[u], trains[v], dead_ms=dead_ms) for v in ids if v != u}
        conflict = (refractory_conflict_ratio(trains[u].size, rejected.get(u, 0), n_frames, assign_refr_frames)
                    if rejected is not None else float("nan"))
        out[u] = {"n_spikes": int(trains[u].size), "rp_contamination": c, "rp_violations": nv,
                  "refractory_conflict_ratio": conflict,
                  "ampsplit_ratio": split.get("ratio"), "ampsplit_burst_ratio": split.get("burst_2_3ms_ratio"),
                  "neighbours": {v: d["verdict"] for v, d in neigh.items()},
                  "neighbour_ratios": {v: d["ratio"] for v, d in neigh.items()}}
    return out


def plot_unit_purity(ax_grid, traces: np.ndarray, gain: float, u: int, trains: dict, amps: dict,
                     r: dict, templates: np.ndarray, template_ids: list, *, nbefore: int = NBEFORE,
                     n_frames: int | None = None, n_overlay: int = 120, seed: int = 0,
                     highlight: np.ndarray | None = None, highlight_label: str = "in pulse",
                     title: str = "", dead_ms: float = 0.0, n_rejected: int | None = None,
                     assign_refr_frames: int = 45) -> None:
    """Draw one unit's purity page into a ``(2, 3)`` axes grid: waveform overlays with the template,
    amplitude against cosine, amplitude over time, ISI histogram with the contamination and the
    refractory conflict ratio, the amplitude-split cross-correlogram and the neighbour
    cross-correlograms with verdicts. ``amps`` are in MAD units (as ``Capture.amps`` gives them).
    ``highlight`` is a boolean mask over the unit's spikes (for example those inside a light pulse)
    drawn in a second colour.
    """
    import matplotlib.pyplot as plt  # noqa: F401  (ensures a backend is loaded by the caller)

    n_frames = n_frames or traces.shape[0]
    rng = np.random.default_rng(seed)
    fr, am, rr = trains[u], amps[u], r[u]
    ti = template_ids.index(u)
    hl = np.zeros(fr.size, bool) if highlight is None else highlight
    tt = (np.arange(templates.shape[1]) - nbefore) / FS * 1000
    ax = ax_grid[0, 0]
    pick = rng.choice(fr.size, min(n_overlay, fr.size), replace=False) if fr.size else np.zeros(0, int)
    if pick.size:
        sn = snippets(traces, fr[pick], nbefore=nbefore, nafter=templates.shape[1] - nbefore) * gain
        for k, j in enumerate(pick):
            for c in range(4):
                ax.plot(tt + c * 3.2, sn[k, :, c], color="tab:orange" if hl[j] else "0.6", lw=0.3, alpha=0.5)
    for c in range(4):
        ax.plot(tt + c * 3.2, templates[ti, :, c] * gain, color="k", lw=1.2)
    ax.set_title(f"{title} unit {u}: {fr.size} spikes; {min(n_overlay, fr.size)} overlaid (orange = {highlight_label})", fontsize=8)
    ax.set_xlabel("ms, wires offset by 3.2 ms"); ax.set_ylabel("uV")
    ax = ax_grid[0, 1]
    ax.scatter(rr[~hl], am[~hl], s=2, c="0.4", alpha=0.4, label="outside")
    if hl.any():
        ax.scatter(rr[hl], am[hl], s=3, c="tab:orange", alpha=0.6, label=highlight_label)
    ax.set_xlabel("cosine to assigned template"); ax.set_ylabel("trough amplitude (MAD)"); ax.legend(fontsize=7)
    ax.set_title("amplitude by shape", fontsize=8)
    ax = ax_grid[0, 2]
    ax.scatter(fr / FS / 60, am, s=1, c=np.where(hl, "tab:orange", "0.4"), alpha=0.4)
    ax.set_xlabel("min"); ax.set_ylabel("trough amplitude (MAD)"); ax.set_title("amplitude over time", fontsize=8)
    ax = ax_grid[1, 0]
    isi = np.diff(fr) / FS * 1000
    c, nv = rp_contamination(fr, n_frames)
    ax.hist(isi[isi < 50], bins=np.arange(0, 50, 0.25), color="k")
    ax.axvline(1.0, color="r", lw=0.6)
    conflict = ("" if n_rejected is None else
                f"; refractory conflict ratio {refractory_conflict_ratio(fr.size, n_rejected, n_frames, assign_refr_frames):.2f} ({n_rejected} lost)")
    ax.set_title(f"ISI; rp contamination {c:.3f} ({nv} < 1 ms){conflict}", fontsize=8); ax.set_xlabel("ms")
    ax = ax_grid[1, 1]
    sp = amplitude_split(fr, am, n_frames, dead_ms=dead_ms)
    if "lag_ms" in sp:
        ax.bar(sp["lag_ms"], sp["counts"], width=0.5, color="k")
        ax.axhline(sp["expected_per_bin"], color="gray", lw=0.8)
        ax.set_title(f"large (>= {sp['large_min']:.1f} MAD, n={sp['n_large']}) vs small (<= {sp['small_max']:.1f} MAD, n={sp['n_small']}): "
                     f"zero-lag ratio {sp['ratio']:.2f}, 2-3 ms {sp['burst_2_3ms_ratio']:.2f}", fontsize=7)
        ax.set_xlabel("ms (small relative to large)")
    ax = ax_grid[1, 2]
    others = [v for v in sorted(trains) if v != u]
    max_lag = int(25 * FS / 1000)
    for k, v in enumerate(others):
        lags = ccg_lags(fr, trains[v], max_lag) / FS * 1000
        h, e = np.histogram(lags, np.arange(-25, 25.5, 0.5))
        verdict = ccg_verdict(fr, trains[v], dead_ms=dead_ms)
        base = h.mean() if h.size else 1
        ax.plot(e[:-1] + 0.25, h / max(base, 1e-9) + 1.5 * k, lw=0.7,
                label=f"unit {v}: {verdict['verdict']} ({verdict['ratio']:.2f})" if np.isfinite(verdict["ratio"]) else f"unit {v}: {verdict['verdict']}")
    ax.set_title("cross-correlograms with same-tetrode units (stacked, normalised)", fontsize=8)
    ax.set_xlabel("ms"); ax.legend(fontsize=6)
