"""Small-spike capture on what a full-size sorting leaves, for units that answer the light.

A full-size sorting (a matcher that accepts a spike at 0.8 of its template's amplitude or more,
then a residual pass for large unclaimed events) keeps a light-driven cell's first spikes after
onset and drops the attenuated spikes it fires through the pulse. This module recovers those one
tetrode and one window at a time, and only for units that answer the light, so a template whose
cell has gone does not keep collecting look-alike events:

1. ``subtract_sorted``: every spike the full-size sorting assigned is fitted with its unit's
   template (least-squares scale, best shift within one sample) and subtracted, deepest first.
   What is left is the residual of the full-size stages.
2. ``residual_events``: troughs of that residual below ``threshold_mad`` times the excised noise of
   the original traces (``ttm.capture.noise_mad_excised``), minus those within the detector's dead
   time of a sorted spike (the remains of an imperfect subtraction).
3. ``offer_events``: each event is scored against every template of the tetrode with the
   full-window cosine (amplitude free); the best template is offered the event when the
   noise-corrected cosine reaches ``gate_rc`` (``ttm.capture.score_events`` gives both cosines),
   and offered events are peeled off the residual over ``n_pass`` passes as in
   ``ttm.capture.capture_tetrode``. Every template competes, so an event that looks most like a
   unit that is not admitted is not handed to one that is.
4. ``light_response`` and ``admit``: whether a unit answers the light at one train, from per-pulse
   spike counts in the pulse (``in_start_s`` after onset to offset) against the same length of
   time before onset, on the unit's sorted spikes plus the events it was offered.
5. ``governing_train`` and ``accept_offered``: an offered event is kept only if its unit was
   admitted at the most recent train before it (before the first train, the first train decides),
   is not within the refractory period of one of the unit's sorted spikes, and wins its chain among
   the unit's other new events.

Nothing here reads a store, a manifest or a stimulation table; times and traces come from the
caller. Templates are ``(n_units, n_samp, 4)`` arrays aligned at ``nbefore`` as in
``ttm.capture``.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from ttm import capture as cp

FS = cp.FS


def subtract_sorted(traces: np.ndarray, frames: np.ndarray, labels: np.ndarray, templates: np.ndarray, *,
                    nbefore: int = cp.NBEFORE, s_max: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Residual of ``traces`` after each sorted spike (``frames``, template index ``labels``; -1 is
    left in place) is fitted with its template and subtracted. The fit is the scale
    ``a = <snippet, t> / |t|^2`` (at least zero) at the shift within ``s_max`` samples that maximises
    the inner product, computed on the residual as it stands, deepest spike first, so the second of
    two overlapping spikes is fitted after the first is removed. Returns (residual, a per spike;
    NaN where the spike was not subtracted)."""
    residual = np.array(traces, dtype=np.float32, copy=True)
    n_samp = templates.shape[1]
    tsq = (templates.astype(np.float64) ** 2).sum(axis=(1, 2))
    a_out = np.full(frames.size, np.nan)
    ok = (labels >= 0) & (frames - nbefore - s_max >= 0) & (frames - nbefore + n_samp + s_max <= traces.shape[0])
    idx = np.flatnonzero(ok)
    if idx.size == 0:
        return residual, a_out
    depth = traces[frames[idx]].min(axis=1)
    for i in idx[np.argsort(depth, kind="stable")]:
        f, t = int(frames[i]), templates[labels[i]]
        best_c, best_s = -np.inf, 0
        for s in range(-s_max, s_max + 1):
            c = float(np.vdot(residual[f + s - nbefore:f + s - nbefore + n_samp], t))
            if c > best_c:
                best_c, best_s = c, s
        a = max(best_c / tsq[labels[i]], 0.0)
        a_out[i] = a
        if a > 0:
            residual[f + best_s - nbefore:f + best_s - nbefore + n_samp] -= np.float32(a) * t
    return residual, a_out


def near(frames: np.ndarray, reference: np.ndarray, tol: int) -> np.ndarray:
    """Boolean per frame: a reference frame within ``tol`` frames (both arrays sorted)."""
    if reference.size == 0 or frames.size == 0:
        return np.zeros(frames.size, bool)
    j = np.searchsorted(reference, frames)
    lo = reference[np.clip(j - 1, 0, reference.size - 1)]
    hi = reference[np.clip(j, 0, reference.size - 1)]
    return np.minimum(np.abs(frames - lo), np.abs(hi - frames)) <= tol


def residual_events(residual: np.ndarray, noise: np.ndarray, sorted_frames: np.ndarray, *,
                    threshold_mad: float, dead_frames: int = 9, nbefore: int = cp.NBEFORE,
                    nafter: int = cp.NAFTER) -> tuple[np.ndarray, np.ndarray]:
    """Troughs of the residual below ``threshold_mad`` times ``noise`` (``ttm.capture.detect_troughs``),
    dropping those within ``dead_frames`` of a sorted spike. Returns (frames, depth in MAD units)."""
    frames, amp = cp.detect_troughs(residual, noise, threshold_mad=threshold_mad, dead_frames=dead_frames,
                                    nbefore=nbefore, nafter=nafter)
    keep = ~near(frames, np.sort(np.asarray(sorted_frames, np.int64)), dead_frames)
    return frames[keep], amp[keep]


@dataclasses.dataclass
class Offer:
    """Every residual event with the template it was scored best against (index, -1 if none), the raw
    and noise-corrected cosines to it, the fitted scale, whether it was offered (best template's
    noise-corrected cosine at the gate) and the peeling pass that offered it."""
    frames: np.ndarray
    amp_mad: np.ndarray
    best: np.ndarray
    r: np.ndarray
    rc: np.ndarray
    a: np.ndarray
    offered: np.ndarray
    pass_index: np.ndarray


def offer_events(residual: np.ndarray, frames: np.ndarray, amp_mad: np.ndarray, templates: np.ndarray,
                 noise: np.ndarray, *, gate_rc: float | None = 0.8, gate_r: float | None = None,
                 n_pass: int = 3, nbefore: int = cp.NBEFORE, s_max: int = 1) -> Offer:
    """Competitive shape assignment of residual events with peeling. The best template (by raw
    cosine over shifts) is offered an event when its noise-corrected cosine reaches ``gate_rc`` and,
    if given, its raw cosine reaches ``gate_r``; offered events are subtracted (scaled) from a working
    copy and the rest scored again, ``n_pass`` passes in all. ``residual`` is not modified."""
    n = frames.size
    best = np.full(n, -1, np.int64)
    r = np.full(n, np.nan)
    rc = np.full(n, np.nan)
    a = np.full(n, np.nan)
    offered = np.zeros(n, bool)
    pass_index = np.full(n, -1, np.int64)
    if templates.shape[0] == 0 or n == 0:
        return Offer(frames, amp_mad, best, r, rc, a, offered, pass_index)
    work = np.array(residual, dtype=np.float32, copy=True)
    pending = np.arange(n)
    for p in range(n_pass):
        if pending.size == 0:
            break
        sc = cp.score_events(work, frames[pending], templates, noise, nbefore=nbefore, s_max=s_max)
        rr = np.where(np.isnan(sc["r"]), -np.inf, sc["r"])
        j = np.argmax(rr, axis=1)
        r_best = rr[np.arange(pending.size), j]
        ok = np.isfinite(r_best) & np.isfinite(sc["rc"])
        if gate_rc is not None:
            ok &= sc["rc"] >= gate_rc
        if gate_r is not None:
            ok &= r_best >= gate_r
        best[pending] = np.where(np.isfinite(r_best), j, -1)
        r[pending] = np.where(np.isfinite(r_best), r_best, np.nan)
        rc[pending] = sc["rc"]
        a[pending] = sc["a"]
        acc = pending[ok]
        offered[acc] = True
        pass_index[acc] = p
        cp.subtract_templates(work, frames[acc] + sc["shift"][ok], j[ok], sc["a"][ok], templates, nbefore=nbefore)
        pending = pending[~ok]
    return Offer(frames, amp_mad, best, r, rc, a, offered, pass_index)


def light_response(spike_s: np.ndarray, onsets: np.ndarray, offsets: np.ndarray, *, in_start_s: float = 0.010,
                   onset_window_s: tuple[float, float] = (0.001, 0.010)) -> dict:
    """Per-pulse counts of a unit's spikes (sorted times, s) in the pulse (``in_start_s`` after onset
    to offset) and in a window of the same length ending at onset ("pre-pulse"; inside a train this
    is between pulses, not a baseline), with the rates, their ratio, a one-sided exact test that the
    in-pulse rate exceeds the pre-pulse rate given the total count, and the paired per-pulse sign
    test. Also the onset rate (``onset_window_s``) for reference."""
    from scipy.stats import binomtest

    onsets = np.asarray(onsets, float)
    offsets = np.asarray(offsets, float)
    if onsets.size == 0:
        return {"n_pulses": 0}
    length = offsets - onsets - in_start_s
    n_in = np.searchsorted(spike_s, offsets) - np.searchsorted(spike_s, onsets + in_start_s)
    n_pre = np.searchsorted(spike_s, onsets) - np.searchsorted(spike_s, onsets - length)
    n_on = np.searchsorted(spike_s, onsets + onset_window_s[1]) - np.searchsorted(spike_s, onsets + onset_window_s[0])
    t = float(length.sum())
    tot_in, tot_pre = int(n_in.sum()), int(n_pre.sum())
    p = float(binomtest(tot_in, tot_in + tot_pre, 0.5, alternative="greater").pvalue) if tot_in + tot_pre else 1.0
    up, down = int((n_in > n_pre).sum()), int((n_in < n_pre).sum())
    p_sign = float(binomtest(up, up + down, 0.5, alternative="greater").pvalue) if up + down else 1.0
    pre_hz, in_hz = tot_pre / t, tot_in / t
    return {"n_pulses": int(onsets.size), "pre_hz": pre_hz, "in_hz": in_hz,
            "onset_hz": float(n_on.sum() / (onsets.size * (onset_window_s[1] - onset_window_s[0]))),
            "ratio": in_hz / pre_hz if pre_hz > 0 else (np.inf if in_hz > 0 else np.nan),
            "n_in": tot_in, "n_pre": tot_pre, "p_count": p, "pulses_up": up, "pulses_down": down, "p_sign": p_sign}


def admit(resp: dict, *, min_ratio: float, max_p: float, min_in_hz: float) -> bool:
    """The admission rule: in-pulse rate at least ``min_ratio`` times the pre-pulse rate, the count
    test below ``max_p`` and the in-pulse rate at least ``min_in_hz``."""
    if not resp.get("n_pulses"):
        return False
    return bool(resp["ratio"] >= min_ratio and resp["p_count"] <= max_p and resp["in_hz"] >= min_in_hz)


def governing_train(times_s: np.ndarray, train_starts_s: np.ndarray) -> np.ndarray:
    """Index of the most recent train start at or before each time; times before the first train get
    the first train (index 0). ``train_starts_s`` sorted."""
    k = np.searchsorted(np.asarray(train_starts_s, float), np.asarray(times_s, float), side="right") - 1
    return np.clip(k, 0, None)


def accept_offered(offer: Offer, eligible: np.ndarray, sorted_trains: dict, template_ids: list, *,
                   refr_frames: int = 45) -> np.ndarray:
    """Boolean mask over ``offer``'s events: offered, ``eligible`` (per event, the event's unit admitted
    at its governing train), farther than ``refr_frames`` from every sorted spike of that unit
    (``sorted_trains[unit]``, sorted frames), and the best raw cosine in its chain of the unit's new
    events closer than ``refr_frames`` (``ttm.capture.resolve_refractory``)."""
    keep = offer.offered & np.asarray(eligible, bool)
    for i, u in enumerate(template_ids):
        idx = np.flatnonzero(keep & (offer.best == i))
        if idx.size == 0:
            continue
        host = sorted_trains.get(u, np.zeros(0, np.int64))
        clash = near(offer.frames[idx], host, refr_frames)
        keep[idx[clash]] = False
        idx = idx[~clash]
        ok = cp.resolve_refractory(offer.frames[idx], offer.r[idx], refr_frames)
        keep[idx[~ok]] = False
    return keep
