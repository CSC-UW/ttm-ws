"""Transient bad periods on tetrodes: detection, an interval table, and masking for sorting.

A bad period is a stretch of time in which one tetrode's wires carry a signal
the sorter must not learn from: a burst shared by all four wires, or one or
two wires several times louder than usual. It is a fact about the recording,
scoped to a tetrode and found before sorting, not a property of any unit.

1. ``bin_stats``: one pass over a preprocessed recording that keeps, per short
   bin, each wire's sum of squares, the cross-products of every within-tetrode
   wire pair, each wire's peak and its count of negative threshold crossings.
   Sums aggregate exactly to any coarser bin, so the scan is run once.
2. ``bin_features``: loudness of each wire relative to its own running median,
   the median wire-to-wire correlation of each tetrode relative to its own
   running median, and the crossing rate, at a chosen bin size.
3. ``detect_bad_periods``: flags per tetrode and bin, joined into intervals,
   written as a table with one row per tetrode and interval.
4. ``mask_bad_periods``: a lazy recording in which each flagged tetrode's
   channels read as zero inside its intervals, with short cosine tapers at the
   edges, so every sorting stage that reads the recording skips them.
5. ``valid_seconds``: the time each tetrode was not masked, for firing rates
   and presence.

The table is the source of truth; it lives next to the recording it was
measured on, as CSV.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from spikeinterface.preprocessing.basepreprocessor import BasePreprocessor, BasePreprocessorSegment


# --- one pass over the recording ------------------------------------------------------------

def group_pairs(groups: np.ndarray) -> np.ndarray:
    """``(n_pairs, 3)`` int array of ``(group, channel_i, channel_j)`` for every within-group pair,
    groups in ascending order and pairs in index order."""
    groups = np.asarray(groups).astype(int)
    rows = []
    for g in np.unique(groups):
        idx = np.flatnonzero(groups == g)
        rows += [(g, i, j) for i, j in itertools.combinations(idx, 2)]
    return np.asarray(rows, dtype=np.int64).reshape(-1, 3)


def _chunk_stats(recording, start: int, end: int, bin_frames: int, pairs: np.ndarray,
                 cross_level: np.ndarray) -> dict:
    x = recording.get_traces(start_frame=start, end_frame=end, return_in_uV=False).astype(np.float32)
    nb = x.shape[0] // bin_frames
    x = x[:nb * bin_frames]
    xb = x.reshape(nb, bin_frames, -1)
    below = x < -cross_level
    onset = below[1:] & ~below[:-1]
    onset = np.vstack([below[:1], onset]).reshape(nb, bin_frames, -1)
    gram = np.matmul(xb.transpose(0, 2, 1), xb)
    return {
        "ss": np.diagonal(gram, axis1=1, axis2=2).copy(),
        "xp": gram[:, pairs[:, 1], pairs[:, 2]],
        "peak": np.abs(xb).max(axis=1),
        "cross": onset.sum(axis=1).astype(np.uint16),
    }


def bin_stats(recording, *, bin_s: float = 0.1, chunk_s: float = 60.0, cross_mad: float = 5.0,
              noise_levels: np.ndarray | None = None, start_frame: int = 0, end_frame: int | None = None,
              n_jobs: int = 32) -> dict:
    """Per-bin sufficient statistics of a single-segment preprocessed recording (raw units).

    Returns a dict of arrays with ``n_bins`` rows: ``ss`` (sum of squares per channel), ``xp``
    (sum of cross-products per within-group pair, columns as in ``pairs``), ``peak`` (largest
    absolute value per channel), ``cross`` (onsets of excursions below ``-cross_mad`` times the
    channel's noise level), plus ``pairs``, ``groups``, ``channel_ids``, ``noise_levels``,
    ``bin_frames``, ``start_frame`` and ``fs``. A trailing partial bin is dropped.
    """
    from joblib import Parallel, delayed
    from spikeinterface.core import get_noise_levels

    fs = recording.get_sampling_frequency()
    groups = np.asarray(recording.get_property("group")).astype(int)
    pairs = group_pairs(groups)
    if noise_levels is None:
        noise_levels = get_noise_levels(recording, return_in_uV=False)
    noise_levels = np.asarray(noise_levels, dtype=np.float32)
    end_frame = recording.get_num_frames() if end_frame is None else int(end_frame)
    bin_frames = int(round(bin_s * fs))
    chunk_frames = int(round(chunk_s / bin_s)) * bin_frames
    starts = range(int(start_frame), end_frame, chunk_frames)
    parts = Parallel(n_jobs=n_jobs)(
        delayed(_chunk_stats)(recording, s, min(s + chunk_frames, end_frame), bin_frames, pairs,
                              cross_mad * noise_levels)
        for s in starts)
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    out.update(pairs=pairs, groups=groups, channel_ids=np.asarray(recording.channel_ids).astype(str),
               noise_levels=noise_levels, bin_frames=bin_frames, start_frame=int(start_frame),
               fs=float(fs), cross_mad=float(cross_mad))
    return out


# --- features -------------------------------------------------------------------------------

def _rebin(a: np.ndarray, k: int, how: str) -> np.ndarray:
    n = a.shape[0] // k
    a = a[:n * k].reshape(n, k, *a.shape[1:])
    return a.max(axis=1) if how == "max" else a.sum(axis=1, dtype=np.float64)


def bin_features(stats: dict, *, bin_s: float = 1.0) -> dict:
    """Features at ``bin_s`` (a multiple of the stats' bin) from ``bin_stats`` output.

    ``rms`` per channel (raw units); ``corr`` per within-group pair; ``corr_med`` per group, the
    median over its pairs; ``cross_hz`` per channel; ``peak`` per channel (raw units); ``t``, the
    bin start in seconds of the recording; plus ``group_ids``, ``groups``, ``channel_ids``, ``bin_s``.
    """
    fs = stats["fs"]
    base_bin_s = int(stats["bin_frames"]) / fs
    k = int(round(bin_s / base_bin_s))
    if abs(k * base_bin_s - bin_s) > 1e-9:
        raise ValueError(f"bin_s={bin_s} is not a multiple of the stats' {base_bin_s} s bin")
    n_frames = k * int(stats["bin_frames"])
    ss = _rebin(stats["ss"], k, "sum")
    xp = _rebin(stats["xp"], k, "sum")
    pairs = stats["pairs"]
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = (xp / np.sqrt(ss[:, pairs[:, 1]] * ss[:, pairs[:, 2]])).astype(np.float32)
    group_ids = np.unique(pairs[:, 0])
    corr_med = np.stack([np.median(corr[:, pairs[:, 0] == g], axis=1) for g in group_ids], axis=1)
    return dict(
        rms=np.sqrt(ss / n_frames).astype(np.float32), corr=corr, corr_med=corr_med.astype(np.float32),
        cross_hz=(_rebin(stats["cross"].astype(np.float32), k, "sum") / bin_s).astype(np.float32),
        peak=_rebin(stats["peak"], k, "max"),
        t=int(stats["start_frame"]) / fs + np.arange(ss.shape[0]) * bin_s,
        group_ids=group_ids, groups=np.asarray(stats["groups"]).astype(int),
        channel_ids=np.asarray(stats["channel_ids"]).astype(str), bin_s=float(bin_s))


# --- detection ------------------------------------------------------------------------------

def _running_median(a: np.ndarray, window_bins: int) -> np.ndarray:
    """Centred running median along axis 0 over non-NaN values, at least a quarter of the window
    present; gaps left by too few values are filled from the nearest estimate."""
    df = pd.DataFrame(a).rolling(window_bins, center=True, min_periods=max(1, window_bins // 4)).median()
    return df.ffill().bfill().to_numpy()


def _dilate(mask: np.ndarray, n: int) -> np.ndarray:
    if n <= 0:
        return mask
    return pd.DataFrame(mask.astype(np.float32)).rolling(2 * n + 1, center=True,
                                                         min_periods=1).max().to_numpy() > 0


def _runs(mask: np.ndarray, max_gap: int) -> list[tuple[int, int]]:
    """``[start, end)`` index runs of True, joining runs separated by at most ``max_gap`` False."""
    k = np.flatnonzero(mask)
    if k.size == 0:
        return []
    cut = np.flatnonzero(np.diff(k) > max_gap + 1)
    return [(int(r[0]), int(r[-1]) + 1) for r in np.split(k, cut + 1)]


def baselines(features: dict, *, baseline_s: float = 1800.0, wire_loud: float = 2.0,
              exclude_pad_s: float = 30.0, n_iter: int = 2) -> dict:
    """Running medians of each wire's rms and each tetrode's median correlation over ``baseline_s``,
    re-estimated ``n_iter`` times without the bins in which any wire of the tetrode is more than
    ``wire_loud`` times its baseline (padded by ``exclude_pad_s``), so a long bad stretch does not
    raise its own baseline. Returns ``rms_base`` and ``corr_base``."""
    rms = features["rms"].astype(np.float64)
    corr_med = features["corr_med"].astype(np.float64)
    groups, gids = features["groups"], features["group_ids"]
    w = max(1, int(round(baseline_s / features["bin_s"])))
    pad = int(round(exclude_pad_s / features["bin_s"]))
    rms_base = _running_median(rms, w)
    corr_base = _running_median(corr_med, w)
    for _ in range(n_iter):
        excl = np.stack([(rms[:, groups == g] / rms_base[:, groups == g]).max(axis=1) > wire_loud
                         for g in gids], axis=1)
        excl = _dilate(excl, pad)
        r = rms.copy()
        for gi, g in enumerate(gids):
            r[excl[:, gi][:, None] & (groups == g)[None, :]] = np.nan
        c = corr_med.copy()
        c[excl] = np.nan
        rms_base = _running_median(r, w)
        corr_base = _running_median(c, w)
    return dict(rms_base=rms_base.astype(np.float32), corr_base=corr_base.astype(np.float32))


TABLE_COLUMNS = ["group", "start_s", "end_s", "kind", "channels", "max_loud", "max_corr_excess",
                 "flagged_s"]


def detect_bad_periods(features: dict, base: dict | None = None, *, wire_loud: float = 3.0,
                       wire_loud_extend: float = 2.0, common_corr_excess: float = 0.25,
                       common_loud: float = 2.0,
                       global_min_groups: int | None = None, max_gap_s: float = 5.0,
                       pad_s: float = 1.0, group_names: dict | None = None) -> pd.DataFrame:
    """Bad periods per tetrode from ``bin_features`` output.

    A bin of a tetrode is flagged ``loud`` when any wire's rms exceeds ``wire_loud`` times its
    baseline, ``common`` when the tetrode's median wire-to-wire correlation exceeds its baseline
    by ``common_corr_excess`` while its median wire is ``common_loud`` times its baseline (a
    shared signal that adds power), and ``global`` on every tetrode when at least
    ``global_min_groups`` tetrodes are flagged in that bin (default: three quarters of them,
    rounded up, so a fault shared by the tetrodes of one connector half stays local).
    Flagged bins of a tetrode join one interval when they are closer than ``max_gap_s`` or when
    some wire stays ``wire_loud_extend`` times its baseline in every bin between them; the
    interval runs from its first to its last flagged bin and is padded by ``pad_s``.

    One row per tetrode and interval: ``group``, ``start_s``, ``end_s``, ``kind`` (the kinds
    present, "+"-joined, most frequent first), ``channels`` (wires that were loud in it, or every
    wire for ``common`` and ``global``), ``max_loud``, ``max_corr_excess``, ``flagged_s``, and
    ``tetrode`` when ``group_names`` maps group to a name.
    """
    base = baselines(features) if base is None else base
    rms, groups, gids = features["rms"], features["groups"], features["group_ids"]
    ch_ids, bin_s, t = features["channel_ids"], features["bin_s"], features["t"]
    loud = rms / base["rms_base"]
    cex = features["corr_med"] - base["corr_base"]
    n_bins = rms.shape[0]
    kinds = np.zeros((n_bins, len(gids)), dtype=np.int8)  # 0 none, 1 loud, 2 common, 3 global
    for gi, g in enumerate(gids):
        lg = loud[:, groups == g]
        common = (cex[:, gi] > common_corr_excess) & (np.median(lg, axis=1) > common_loud)
        kinds[common, gi] = 2
        kinds[lg.max(axis=1) > wire_loud, gi] = 1
    if global_min_groups is None:
        global_min_groups = int(np.ceil(0.75 * len(gids)))
    is_global = (kinds > 0).sum(axis=1) >= global_min_groups
    kinds[is_global, :] = 3
    names = {1: "loud", 2: "common", 3: "global"}
    gap, pad = int(round(max_gap_s / bin_s)), int(round(pad_s / bin_s))
    rows = []
    for gi, g in enumerate(gids):
        chans = np.flatnonzero(groups == g)
        flagged = kinds[:, gi] > 0
        extend = flagged | (loud[:, chans].max(axis=1) > wire_loud_extend)
        for a, b in _runs(extend, gap):
            seeds = np.flatnonzero(flagged[a:b])
            if seeds.size == 0:
                continue
            a, b = a + int(seeds[0]), a + int(seeds[-1]) + 1
            k = kinds[a:b, gi]
            present = [names[v] for v in sorted(set(k[k > 0].tolist()), key=lambda v: -(k == v).sum())]
            hot = (loud[a:b][:, chans] > wire_loud) & (k == 1)[:, None]
            wires = chans[hot.any(axis=0)] if present == ["loud"] else chans
            rows.append(dict(
                group=int(g), start_s=float(t[max(a - pad, 0)]),
                end_s=float(t[min(b + pad, n_bins) - 1] + bin_s), kind="+".join(present),
                channels=",".join(ch_ids[wires]), max_loud=float(np.nanmax(loud[a:b][:, chans])),
                max_corr_excess=float(np.nanmax(cex[a:b, gi])), flagged_s=float((k > 0).sum() * bin_s)))
    df = pd.DataFrame(rows, columns=TABLE_COLUMNS)
    if group_names is not None:
        df.insert(1, "tetrode", df["group"].map(group_names))
    return df.sort_values(["start_s", "group"], ignore_index=True)


# --- the table ------------------------------------------------------------------------------

def read_table(path) -> pd.DataFrame:
    """A bad-period table written by ``detect_bad_periods(...).to_csv(path, index=False)``."""
    return pd.read_csv(path, dtype={"channels": str})


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out = []
    for a, b in sorted(intervals):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def group_intervals(table: pd.DataFrame, fs: float, n_frames: int) -> dict[int, list[tuple[int, int]]]:
    """``{group: [(start_frame, end_frame), ...]}``, clipped to the recording and merged."""
    out = {}
    for g, sub in table.groupby("group"):
        iv = [(max(0, int(np.floor(a * fs))), min(n_frames, int(np.ceil(b * fs))))
              for a, b in zip(sub["start_s"], sub["end_s"])]
        out[int(g)] = _merge([(a, b) for a, b in iv if b > a])
    return out


def valid_seconds(table: pd.DataFrame, duration_s: float, group_ids) -> dict[int, float]:
    """Seconds of each group not covered by its bad periods."""
    out = {}
    for g in group_ids:
        sub = table[table["group"] == g]
        iv = _merge([(max(0.0, a), min(duration_s, b)) for a, b in zip(sub["start_s"], sub["end_s"])])
        out[int(g)] = float(duration_s - sum(max(0.0, b - a) for a, b in iv))
    return out


def bad_mask(table: pd.DataFrame, group: int, frames: np.ndarray, fs: float) -> np.ndarray:
    """True for each of ``frames`` that falls inside a bad period of ``group``."""
    sub = table[table["group"] == group]
    iv = _merge([(a * fs, b * fs) for a, b in zip(sub["start_s"], sub["end_s"])])
    if not iv:
        return np.zeros(len(frames), dtype=bool)
    starts = np.array([a for a, _ in iv])
    ends = np.array([b for _, b in iv])
    k = np.searchsorted(starts, frames, side="right") - 1
    return (k >= 0) & (frames < ends[np.clip(k, 0, None)])


# --- masking --------------------------------------------------------------------------------

def clean_noise_levels(recording, periods: list, *, n_chunks: int = 60, chunk_frames: int = 10000,
                       seed: int = 0) -> np.ndarray:
    """Per-channel MAD noise (raw units) over random chunks of segment 0 that avoid each channel's
    periods (``[(start_frame, end_frame, channel_indices), ...]``)."""
    n = recording.get_num_frames(0)
    rng = np.random.default_rng(seed)
    starts = np.sort(rng.integers(0, max(1, n - chunk_frames), size=n_chunks))
    n_ch = recording.get_num_channels()
    keep = np.ones((n_chunks, n_ch), dtype=bool)
    for a, b, chans in periods:
        hit = (starts < b) & (starts + chunk_frames > a)
        keep[np.ix_(hit, np.asarray(chans, dtype=int))] = False
    vals = [[] for _ in range(n_ch)]
    for i, s in enumerate(starts):
        x = recording.get_traces(segment_index=0, start_frame=int(s), end_frame=int(s) + chunk_frames,
                                 return_in_uV=False)
        for c in np.flatnonzero(keep[i]):
            vals[c].append(x[:, c])
    return np.array([np.median(np.abs(np.concatenate(v))) / 0.6745 if v else np.nan for v in vals],
                    dtype=np.float32)


class MaskedPeriodsRecording(BasePreprocessor):
    """Replaces chosen channels of a recording during chosen periods (segment 0 only).

    ``periods`` is ``[(start_frame, end_frame, [channel_index, ...]), ...]``. In ``mode="noise"``
    the samples are replaced by white Gaussian noise at each channel's ``noise_levels`` (computed
    with ``clean_noise_levels`` when not given), so noise estimates made on the result are not
    pulled down; in ``mode="zeros"`` they read as zero. Either way the replacement fades in and
    out over ``taper_frames`` outside each period with a raised-cosine weight. The noise is the
    same whichever chunk boundaries a reader uses.
    """

    def __init__(self, recording, periods, mode: str = "noise", taper_frames: int = 30,
                 noise_levels=None, seed: int = 0):
        if recording.get_num_segments() != 1:
            raise ValueError("MaskedPeriodsRecording supports one segment")
        if mode not in ("noise", "zeros"):
            raise ValueError(f"mode must be 'noise' or 'zeros', not {mode!r}")
        periods = sorted((int(a), int(b), sorted(int(c) for c in ch)) for a, b, ch in periods)
        if mode == "noise" and noise_levels is None:
            noise_levels = clean_noise_levels(recording, periods, seed=seed)
        BasePreprocessor.__init__(self, recording)
        levels = None if noise_levels is None else np.asarray(noise_levels, dtype=np.float32)
        self.add_recording_segment(_MaskedPeriodsSegment(recording.segments[0], periods, mode,
                                                         int(taper_frames), levels, int(seed),
                                                         int(recording.sampling_frequency),
                                                         recording.get_num_channels()))
        self._kwargs = dict(recording=recording, periods=[[a, b, ch] for a, b, ch in periods], mode=mode,
                            taper_frames=int(taper_frames),
                            noise_levels=None if noise_levels is None else [float(v) for v in noise_levels],
                            seed=int(seed))


class _MaskedPeriodsSegment(BasePreprocessorSegment):
    def __init__(self, parent_segment, periods, mode, taper_frames, noise_levels, seed, block_frames,
                 n_channels):
        BasePreprocessorSegment.__init__(self, parent_segment)
        self.periods = periods
        self.starts = np.array([p[0] for p in periods], dtype=np.int64)
        self.ends = np.array([p[1] for p in periods], dtype=np.int64)
        self.max_len = int((self.ends - self.starts).max()) if periods else 0
        self.mode = mode
        self.taper = taper_frames
        self.noise_levels = noise_levels
        self.seed = seed
        self.block = block_frames
        self.n_ch = n_channels

    def _noise(self, start: int, end: int, chans: np.ndarray, need: np.ndarray) -> np.ndarray:
        """Noise for frames [start, end) of ``chans``, drawn only in blocks where ``need`` (frames x
        chans) has a True; zero elsewhere."""
        out = np.zeros((end - start, len(chans)), dtype=np.float32)
        for blk in range(start // self.block, (end - 1) // self.block + 1):
            b0 = blk * self.block
            a, b = max(start, b0), min(end, b0 + self.block)
            for k, c in enumerate(chans):
                if not need[a - start:b - start, k].any():
                    continue
                z = np.random.default_rng((self.seed, blk, int(c))).standard_normal(self.block, dtype=np.float32)
                out[a - start:b - start, k] = z[a - b0:b - b0]
        return out * self.noise_levels[chans]

    def _weight(self, a: int, b: int, start: int, end: int) -> np.ndarray:
        w = np.zeros(end - start, dtype=np.float32)
        lo, hi = max(start, a - self.taper), min(end, b + self.taper)
        t = np.arange(lo, hi)
        v = ((t >= a) & (t < b)).astype(np.float32)
        if self.taper > 0:
            ramp_in = t < a
            ramp_out = t >= b
            v[ramp_in] = 0.5 - 0.5 * np.cos(np.pi * (t[ramp_in] - (a - self.taper)) / self.taper)
            v[ramp_out] = 0.5 + 0.5 * np.cos(np.pi * (t[ramp_out] - b + 1) / self.taper)
        w[lo - start:hi - start] = v
        return w

    def get_traces(self, start_frame, end_frame, channel_indices):
        start_frame = 0 if start_frame is None else int(start_frame)
        end_frame = self.get_num_samples() if end_frame is None else int(end_frame)
        traces = self.parent_recording_segment.get_traces(start_frame, end_frame, channel_indices)
        if not self.periods:
            return traces
        lo, hi = start_frame - self.taper, end_frame + self.taper
        first = np.searchsorted(self.starts, lo - self.max_len, side="left")
        last = np.searchsorted(self.starts, hi, side="left")
        hits = [i for i in range(first, last) if self.ends[i] > lo]
        if not hits:
            return traces
        all_ch = np.arange(self.n_ch) if channel_indices is None else np.arange(self.n_ch)[channel_indices]
        pos = {int(c): k for k, c in enumerate(all_ch)}
        spans = []
        for i in hits:
            a, b, chans = self.periods[i]
            cols = [pos[c] for c in chans if c in pos]
            s0, s1 = max(start_frame, a - self.taper), min(end_frame, b + self.taper)
            if cols and s0 < s1:
                spans.append((s0, s1, a, b, cols))
        if not spans:
            return traces
        spans.sort(key=lambda s: s[0])
        out = np.array(traces, dtype=traces.dtype, copy=True)
        j = 0
        while j < len(spans):
            s0, s1 = spans[j][0], spans[j][1]
            k = j + 1
            while k < len(spans) and spans[k][0] < s1:
                s1 = max(s1, spans[k][1])
                k += 1
            col_w = {}
            for _, _, a, b, cols in spans[j:k]:
                w = self._weight(a, b, s0, s1)
                for col in cols:
                    col_w[col] = w if col not in col_w else np.maximum(col_w[col], w)
            cols = np.array(sorted(col_w))
            weight = np.stack([col_w[c] for c in cols], axis=1)
            x = traces[s0 - start_frame:s1 - start_frame]
            mixed = x[:, cols].astype(np.float32) * (1.0 - weight)
            if self.mode == "noise":
                mixed += weight * self._noise(s0, s1, all_ch[cols], weight > 0)
            out[s0 - start_frame:s1 - start_frame, cols] = mixed.astype(traces.dtype, copy=False)
            j = k
        return out


def mask_bad_periods(recording, table: pd.DataFrame, *, mode: str = "noise", taper_frames: int = 30,
                     noise_levels=None, seed: int = 0) -> MaskedPeriodsRecording:
    """Lazy recording in which each group's channels are replaced inside that group's bad periods."""
    groups = np.asarray(recording.get_property("group")).astype(int)
    iv = group_intervals(table, recording.sampling_frequency, recording.get_num_frames(0))
    periods = [(a, b, np.flatnonzero(groups == g).tolist()) for g, lst in iv.items() for a, b in lst]
    return MaskedPeriodsRecording(recording, periods, mode=mode, taper_frames=taper_frames,
                                  noise_levels=noise_levels, seed=seed)
