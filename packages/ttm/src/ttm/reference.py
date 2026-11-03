"""Median reference with each wire aligned in time to the median first.

The channels of a headstage are not sampled at the same instant, and each
electrode may delay its signal further. A signal common to every wire (the
reference electrode's own) then reaches each wire at its own offset from the
median of all wires, and subtracting the median leaves part of it. Per wire,
a least-squares fit of the bandpassed trace x to the median m and its slope m'
(central difference), x = g m + b m', gives the offset: a wire that records
the common signal ``tau`` samples after the median has x(t) = g m(t - tau),
close to g m - g tau m' while tau is short against the signal, so
tau = -b / g. The aligned reference shifts each wire by -tau (a
fractional-sample shift by FFT, as spikeinterface's ``phase_shift``), takes
the median of the shifted wires and subtracts each shifted wire's
least-squares gain on that median times the median.

Offsets drift over hours, so they are estimated per window and applied
piecewise: ``estimate_wire_timing_windows`` returns the window starts with
one row of shifts and gains per window (each the median over a few blocks
spread across the window, optionally pooled with the neighbouring windows'
blocks), and ``AlignedMedianReferenceRecording``
applies the row of the window each sample falls in.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
from spikeinterface.core import get_chunk_with_margin
from spikeinterface.preprocessing.basepreprocessor import BasePreprocessor, BasePreprocessorSegment
from spikeinterface.preprocessing.phase_shift import apply_frequency_shift

EDGE = 3000


def event_mask(m: np.ndarray, event_mad: float = 5.0, pad: int = 15) -> np.ndarray:
    """Samples within ``pad`` of an excursion of ``m`` beyond ``event_mad`` of its MAD (either sign)."""
    mad = np.median(np.abs(m - np.median(m))) / 0.6745
    hit = np.abs(m - np.median(m)) > event_mad * mad
    k = np.ones(2 * pad + 1)
    return np.convolve(hit.astype(float), k, mode="same") > 0


def estimate_wire_timing(
    traces: np.ndarray, edge: int = EDGE, n_iter: int = 3, event_mad: float = 5.0, pad: int = 15,
    min_event_samples: int = 3000,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Per-wire offset after the median (samples), gain on the aligned median, and event samples used.

    ``traces`` is (samples, channels) of bandpassed, unreferenced data. The fits use only samples
    within ``pad`` of the median's excursions beyond ``event_mad`` of its MAD, where the large common
    events dominate it (all samples if there are fewer than ``min_event_samples``, which the caller
    sees in the returned count); between events the median
    follows smaller common signals whose spread across wires differs. The first-order fit overstates
    large offsets, so each further iteration shifts the wires by the offsets so far and fits what is
    left against the median of the shifted wires. ``edge`` samples at each end, where the FFT shift
    wraps around, are left out.
    """
    x = np.asarray(traces, dtype=np.float64)
    keep = event_mask(np.median(x, axis=1), event_mad, pad)
    keep[:edge] = False
    keep[-edge:] = False
    n_event = int(keep.sum())
    if n_event < min_event_samples:
        keep[edge:-edge] = True
    tau = np.zeros(x.shape[1])
    for _ in range(n_iter):
        y = apply_frequency_shift(x, -tau, axis=0) if tau.any() else x
        m = np.median(y, axis=1)
        A = np.stack([m, np.gradient(m)], axis=1)[keep]
        coef, *_ = np.linalg.lstsq(A, y[keep], rcond=None)
        tau = tau - coef[1] / coef[0]
    y = apply_frequency_shift(x, -tau, axis=0)[keep]
    my = np.median(y, axis=1)
    gain = (y * my[:, None]).sum(axis=0) / float((my**2).sum())
    return tau, gain, n_event


def estimate_wire_timing_blocks(
    recording, start_frame: int, end_frame: int, block_frames: int, *, n_blocks: int = 6,
    min_event_samples: int = 3000, segment_index: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Offsets (``tau``, samples) and gains of the blocks of one window that have enough events.

    ``n_blocks`` blocks of ``block_frames`` centred at (i + 1) / (n_blocks + 1) of the window, each
    fitted with ``estimate_wire_timing``; returns (k, channels) arrays for the k blocks with at least
    ``min_event_samples`` event samples (k may be 0). A window shorter than ``n_blocks`` blocks gets
    as many as fit, at least one, so copies of one span do not count as several blocks.
    ``recording`` is bandpassed and unreferenced.
    """
    span = end_frame - start_frame
    n_blocks = min(n_blocks, max(1, span // block_frames))
    taus, gains = [], []
    for i in range(n_blocks):
        c = start_frame + span * (i + 1) // (n_blocks + 1)
        b0 = max(start_frame, c - block_frames // 2)
        b1 = min(end_frame, b0 + block_frames)
        x = recording.get_traces(segment_index=segment_index, start_frame=int(b0), end_frame=int(b1))
        tau, gain, n_event = estimate_wire_timing(x, edge=min(EDGE, (b1 - b0) // 10),
                                                  min_event_samples=min_event_samples)
        if n_event >= min_event_samples:
            taus.append(tau)
            gains.append(gain)
    n_ch = recording.get_num_channels()
    return np.array(taus).reshape(-1, n_ch), np.array(gains).reshape(-1, n_ch)


def estimate_wire_timing_window(
    recording, start_frame: int, end_frame: int, block_frames: int, *, n_blocks: int = 6,
    min_event_samples: int = 3000, segment_index: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Shifts to apply (``-tau``, samples) and gains for one window, or NaN rows.

    The per-wire median over the blocks of ``estimate_wire_timing_blocks``, so one block inside an
    artefact does not decide the window when others qualify. NaN when no block has enough events.
    """
    taus, gains = estimate_wire_timing_blocks(recording, start_frame, end_frame, block_frames, n_blocks=n_blocks,
                                              min_event_samples=min_event_samples, segment_index=segment_index)
    return pool_blocks([(taus, gains)])


def pool_blocks(blocks: list[tuple[np.ndarray, np.ndarray]], min_blocks: int = 1) -> tuple[np.ndarray, np.ndarray]:
    """Shifts (``-tau``) and gains as the per-wire median over every block in ``blocks``.

    ``blocks`` is a list of (taus, gains) pairs of (k, channels) arrays, as from
    ``estimate_wire_timing_blocks``; NaN rows when fewer than ``min_blocks`` blocks are given.
    """
    taus = np.concatenate([t for t, _ in blocks])
    gains = np.concatenate([g for _, g in blocks])
    if taus.shape[0] < max(min_blocks, 1):
        nan = np.full(taus.shape[1], np.nan)
        return nan, nan.copy()
    return -np.median(taus, axis=0), np.median(gains, axis=0)


def fill_missing_windows(timing: dict, if_none: str = "raise") -> dict:
    """Rows of NaN (windows without enough events) take the nearest window's row that has values.

    When no window has values, ``if_none="raise"`` raises ValueError and ``"identity"`` gives every
    window shift 0 and gain 1 (the global median).
    """
    shifts, gains = timing["shifts"].copy(), timing["gains"].copy()
    missing = np.isnan(shifts).any(axis=1)
    ok = np.flatnonzero(~missing)
    if ok.size == 0:
        if if_none != "identity":
            raise ValueError("no window has enough common events to estimate the wire timing")
        return {**timing, "shifts": np.zeros_like(shifts), "gains": np.ones_like(gains), "filled": missing}
    for w in np.flatnonzero(missing):
        src = ok[np.argmin(np.abs(ok - w))]
        shifts[w], gains[w] = shifts[src], gains[src]
    return {**timing, "shifts": shifts, "gains": gains, "filled": missing}


def estimate_wire_timing_windows(
    recording, window_frames: int, block_frames: int, *, n_blocks: int = 6, pool_windows: int = 2,
    min_blocks: int = 3, if_none: str = "identity", segment_index: int = 0, n_jobs: int = 1,
) -> dict:
    """Wire timing for each window of ``window_frames`` from frame 0, with ``fill_missing_windows``.

    Each window's row is the per-wire median over the qualifying blocks
    (``estimate_wire_timing_blocks``) of the window and of the ``pool_windows`` windows on each
    side, NaN (then filled from the nearest window) when fewer than ``min_blocks`` qualify, and
    shift 0 and gain 1 everywhere when no window has enough (``if_none``). Offsets change little
    over hours, and pooling keeps one block inside an artefact from deciding a window that has no
    other qualifying block. ``pool_windows=0, min_blocks=1`` uses the window's own blocks only.
    ``n_jobs`` windows are fitted at a time, on threads.

    Returns ``{"starts", "shifts", "gains", "filled"}``: the window start frames and, per window,
    each wire's shift to apply (samples), gain, and whether the row was taken from another window.
    """
    n = recording.get_num_frames(segment_index)
    starts = np.arange(0, n, window_frames, dtype=np.int64)

    def fit(s):
        return estimate_wire_timing_blocks(recording, int(s), int(min(s + window_frames, n)), block_frames,
                                           n_blocks=n_blocks, segment_index=segment_index)

    if n_jobs > 1:
        with ThreadPoolExecutor(n_jobs) as ex:
            blocks = list(ex.map(fit, starts))
    else:
        blocks = [fit(s) for s in starts]
    rows = [pool_blocks(blocks[max(0, w - pool_windows):w + pool_windows + 1], min_blocks)
            for w in range(starts.size)]
    return fill_missing_windows(
        {"starts": starts, "shifts": np.array([r[0] for r in rows]), "gains": np.array([r[1] for r in rows])},
        if_none=if_none)


def slice_wire_timing(timing: dict, start_frame: int, end_frame: int) -> dict:
    """The rows of ``timing`` that cover frames ``start_frame`` to ``end_frame``, with starts from the slice's frame 0.

    For referencing a crop of the recording with the timing estimated on all of it, so pooling over
    neighbouring windows still sees the windows outside the crop.
    """
    starts = np.asarray(timing["starts"])
    w0 = int(np.searchsorted(starts, start_frame, side="right") - 1)
    w1 = int(np.searchsorted(starts, end_frame, side="left"))
    rows = slice(max(w0, 0), max(w1, w0 + 1))
    out = {k: np.asarray(v)[rows].copy() for k, v in timing.items() if k in ("starts", "shifts", "gains", "filled")}
    out["starts"] = np.maximum(out["starts"] - start_frame, 0)
    return out


class AlignedMedianReferenceRecording(BasePreprocessor):
    """Shift each wire by its window's shift, then subtract its gain times the median of the shifted wires.

    ``starts`` (n_windows,) are window start frames, the first 0; ``shifts`` and ``gains`` are
    (n_windows, n_channels). Traces are read with ``margin_ms`` on each side (tapered, zeros past the
    ends) for the FFT shift; a chunk that crosses a window start is computed once per window.
    """

    def __init__(self, recording, starts, shifts, gains, margin_ms: float = 40.0, dtype="float32"):
        starts = np.asarray(starts, dtype=np.int64)
        shifts = np.asarray(shifts, dtype=np.float64)
        gains = np.asarray(gains, dtype=np.float64)
        n_ch = recording.get_num_channels()
        assert starts.ndim == 1 and starts[0] == 0 and np.all(np.diff(starts) > 0), "starts must rise from 0"
        assert shifts.shape == gains.shape == (starts.size, n_ch), "one row of shifts and gains per window"
        assert recording.get_num_segments() == 1, "one segment"
        margin = int(margin_ms * recording.get_sampling_frequency() / 1000.0)
        BasePreprocessor.__init__(self, recording, dtype=dtype)
        self.add_recording_segment(
            AlignedMedianReferenceSegment(recording.segments[0], starts, shifts, gains, margin, np.dtype(dtype))
        )
        self._kwargs = dict(recording=recording, starts=starts.tolist(), shifts=shifts.tolist(),
                            gains=gains.tolist(), margin_ms=float(margin_ms), dtype=str(np.dtype(dtype)))


class AlignedMedianReferenceSegment(BasePreprocessorSegment):
    def __init__(self, parent_segment, starts, shifts, gains, margin, dtype):
        BasePreprocessorSegment.__init__(self, parent_segment)
        self.starts, self.shifts, self.gains, self.margin, self.dtype = starts, shifts, gains, margin, dtype

    def get_traces(self, start_frame, end_frame, channel_indices):
        x, left, right = get_chunk_with_margin(
            self.parent_recording_segment, start_frame, end_frame, None, self.margin,
            dtype=np.float64, add_zeros=True, window_on_margin=True,
        )
        n = end_frame - start_frame
        out = np.empty((n, x.shape[1]), dtype=np.float64)
        k0 = int(np.searchsorted(self.starts, start_frame, side="right") - 1)
        k1 = int(np.searchsorted(self.starts, end_frame - 1, side="right") - 1)
        for k in range(k0, k1 + 1):
            a = max(start_frame, int(self.starts[k])) - start_frame
            b = (min(end_frame, int(self.starts[k + 1])) if k + 1 < self.starts.size else end_frame) - start_frame
            y = apply_frequency_shift(x, self.shifts[k], axis=0)[left + a:left + b]
            out[a:b] = y - self.gains[k][None, :] * np.median(y, axis=1, keepdims=True)
        if channel_indices is not None:
            out = out[:, channel_indices]
        return out.astype(self.dtype)


def aligned_median_reference(recording, wire_timing: dict, **kwargs) -> AlignedMedianReferenceRecording:
    """``AlignedMedianReferenceRecording`` from the dict of ``estimate_wire_timing_windows``."""
    return AlignedMedianReferenceRecording(
        recording, wire_timing["starts"], wire_timing["shifts"], wire_timing["gains"], **kwargs
    )
