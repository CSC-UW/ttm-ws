"""The analyzer stage of ``ttm-sort`` (``docs/sorting_pipeline.md``, step 8).

On the binary a sorting was made from:

- ``chunk_noise_levels``: each channel's noise, the mean MAD of the 0.5 s chunks of a grid outside its
  tetrode's bad intervals, and the mean standard deviation of the same chunks (for sd_ratio).
- ``run_multi_unit`` and ``multi_unit_events``: per window, the detected events that no sorted spike of
  their tetrode claims, scored against the window's templates after the tetrode's sorted spikes around
  them are subtracted; per tetrode, those at or above the score that ``mua_loss`` % of the events one
  sorted spike claims fall below, with a trough of ``mua_floor`` or more after the subtraction.
- ``with_multi_unit``: the sorting plus one multi-unit unit per tetrode (property ``is_mua``).
- ``unit_periods``: each unit's tetrode's time outside its bad intervals, as SpikeInterface periods.
- ``all_spike_pc_metrics``: isolation distance, L-ratio, d' and silhouette on every spike's principal
  components, multi-unit spikes counted among the other spikes.
- ``set_sd_ratio_noise``: the noise SpikeInterface's sd_ratio divides by, fixed.
- ``build_analyzer``: the SortingAnalyzer with its extensions, these noise levels, the metrics over the
  periods (sd_ratio with the chunks' standard deviation) and the all-spike PC metrics.

Parameters are those of ``AnalyzerParams``. Noise levels and the multi-unit detection are in ADC counts of
the binary; the analyzer is in µV.
"""

from __future__ import annotations

import dataclasses
import json
import pathlib
import time

import numpy as np
import pandas as pd

from ttm.carry import in_intervals, log, open_sorting_recording

MAD_TO_SD = 0.6744897501960817

METRIC_NAMES = ("num_spikes", "firing_rate", "presence_ratio", "snr", "isi_violation", "rp_violation",
                "sliding_rp_violation", "synchrony", "firing_range", "amplitude_cv", "amplitude_cutoff",
                "noise_cutoff", "amplitude_median", "sd_ratio")
PC_METRIC_NAMES = ("isolation_distance", "l_ratio", "d_prime", "silhouette")


@dataclasses.dataclass(frozen=True)
class AnalyzerParams:
    seed: int = 0
    noise_chunks: int | None = None
    noise_chunk_s: float = 0.5
    mua_window_s: float = 1800.0
    mua_detect: float = 5.5
    mua_floor: float = 7.0
    mua_radius_um: float = 40.0
    mua_claim_ms: float = 0.5
    mua_ms_before: float = 1.0
    mua_ms_after: float = 2.0
    mua_template_min_spikes: int = 100
    mua_template_spikes: int = 500
    mua_extend_ms: float = 3.5
    mua_fit_ms: tuple[float, float] = (-1 / 3, 2 / 3)
    mua_trough_ms: float = 0.1
    mua_loss: float = 5.0
    mua_min_reference: int = 1000
    max_spikes_per_unit: int = 5000
    ms_before: float = 1.0
    ms_after: float = 2.0
    n_pcs: int = 5
    pc_mode: str = "by_channel_local"
    correlogram_window_ms: float = 100.0
    correlogram_bin_ms: float = 1.0
    metric_names: tuple[str, ...] = METRIC_NAMES

    def to_dict(self) -> dict:
        return {k: list(v) if isinstance(v, tuple) else v for k, v in dataclasses.asdict(self).items()}


def samples(ms: float, fs: float) -> int:
    return round(ms * fs / 1000)


def bad_arrays(bad_periods: dict) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """``{group: (starts, ends)}`` sorted, from ``{group: [(start, end), ...]}``."""
    out = {}
    for g, v in bad_periods.items():
        if len(v):
            s = np.array(sorted(v), np.int64)
            out[int(g)] = (s[:, 0], s[:, 1])
    return out


def in_bad(frames: np.ndarray, g: int, bad: dict) -> np.ndarray:
    if g not in bad:
        return np.zeros(frames.size, bool)
    return in_intervals(frames, *bad[g])


# --- noise levels ----------------------------------------------------------------------------------

def _chunk_noise(args):
    rec_path, drop, starts, length = args
    rec = open_sorting_recording(rec_path, drop)
    mad = np.empty((len(starts), rec.get_num_channels()), np.float32)
    sd = np.empty_like(mad)
    for i, s in enumerate(starts):
        x = rec.get_traces(start_frame=int(s), end_frame=int(s) + length, return_in_uV=False)
        med = np.median(x, axis=0)
        mad[i] = np.median(np.abs(x - med), axis=0) / MAD_TO_SD
        sd[i] = x.std(axis=0, dtype=np.float64)
    return mad, sd


def chunk_noise_levels(rec_path, bad_periods: dict, params: AnalyzerParams | None = None, *,
                       drop_channel_ids=(), workers: int = 32) -> tuple[np.ndarray, pd.DataFrame]:
    """Each channel's noise in ADC counts: the mean, over the chunks of ``noise_chunk_s`` on a grid
    that do not overlap its tetrode's bad intervals, of the chunk's MAD / 0.6745. With ``noise_chunks``
    None (the default) every such chunk counts. The table also holds the mean of the same chunks'
    standard deviations (``std_counts``, ``std_uv``), the noise SpikeInterface's sd_ratio is defined
    with.

    With ``noise_chunks`` a number, one seeded permutation of the grid serves every tetrode: each takes
    its first ``noise_chunks`` chunks clear of its bad intervals. Returns the levels and a table of them
    per channel with the chunks used."""
    from concurrent.futures import ProcessPoolExecutor

    params = params or AnalyzerParams()
    drop = [str(c) for c in drop_channel_ids]
    rec = open_sorting_recording(rec_path, drop)
    fs = rec.get_sampling_frequency()
    length = samples(params.noise_chunk_s * 1000, fs)
    n_grid = rec.get_num_frames() // length
    groups = np.asarray(rec.get_property("group")).astype(int)
    bad = bad_arrays(bad_periods)
    order = np.random.default_rng(params.seed).permutation(n_grid)
    starts = order * length
    pick = {}
    for g in np.unique(groups).tolist():
        clear = np.ones(n_grid, bool)
        if g in bad:
            s, e = bad[g]
            i0 = np.searchsorted(e, starts, "right")
            hit = (i0 < s.size) & (s[np.minimum(i0, s.size - 1)] < starts + length)
            clear = ~hit
        pick[g] = np.sort(order[clear]) if params.noise_chunks is None else order[clear][:params.noise_chunks]
    used = np.unique(np.concatenate(list(pick.values())))
    batches = np.array_split(used * length, max(1, min(workers * 4, used.size)))
    with ProcessPoolExecutor(min(workers, len(batches))) as ex:
        res = list(ex.map(_chunk_noise, [(str(rec_path), drop, b, length) for b in batches]))
    mads = np.concatenate([r[0] for r in res])
    sds = np.concatenate([r[1] for r in res])
    noise = np.empty(len(groups), np.float64)
    std = np.empty(len(groups), np.float64)
    n_used = np.empty(len(groups), np.int64)
    for ch, g in enumerate(groups.tolist()):
        idx = np.searchsorted(used, pick[g])
        noise[ch] = mads[idx, ch].astype(np.float64).mean()
        std[ch] = sds[idx, ch].astype(np.float64).mean()
        n_used[ch] = idx.size
    gain = rec.get_channel_gains()
    table = pd.DataFrame({"channel": np.asarray(rec.channel_ids).astype(str), "group": groups,
                          "noise_counts": noise, "noise_uv": noise * gain, "std_counts": std,
                          "std_uv": std * gain, "chunks": n_used})
    return noise, table


# --- multi-unit events -----------------------------------------------------------------------------

_W: dict = {}


def _mua_init(rec_path, drop, sorting_path, noise, bad, params):
    import spikeinterface as si

    rec = open_sorting_recording(rec_path, drop)
    s = si.load(sorting_path)
    groups = np.asarray(rec.get_property("group")).astype(int)
    uids = np.asarray(s.unit_ids)
    ugroup = np.asarray(s.get_property("group")).astype(int)
    trains = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in uids.tolist()}
    pool_f, pool_u = {}, {}
    for g in np.unique(ugroup).tolist():
        us = uids[ugroup == g].tolist()
        f = np.concatenate([trains[u] for u in us])
        u = np.concatenate([np.full(trains[u].size, i, np.int64) for i, u in zip(np.flatnonzero(ugroup == g), us)])
        o = np.argsort(f, kind="stable")
        pool_f[g], pool_u[g] = f[o], u[o]
    _W.update(rec=rec, groups=groups, uids=uids, ugroup=ugroup, trains=trains, pool_f=pool_f, pool_u=pool_u,
              noise=np.asarray(noise, np.float64), bad=bad, params=params)


def _residual_scorer(xt, ra, T, tn, E, ecol, pf, pu, off, fit, ext_n, before, after, claim):
    """Functions scoring windows of the tetrode traces ``xt`` (frame ``ra`` at row 0) raw and after
    subtracting the tetrode's sorted spikes (``pf``, unit indices ``pu``) as extended templates ``E``."""
    n_t = xt.shape[0]

    def take(frames, rel):
        return xt[np.clip((frames - ra)[:, None] + rel, 0, n_t - 1)].astype(np.float32)

    def cosines(s):
        dots = np.einsum("ntc,utc->nu", s, T)
        sn = np.sqrt(np.einsum("ntc,ntc->n", s, s))
        return dots, dots / np.maximum(sn[:, None] * tn[None], 1e-12)

    def residual(frames, s):
        i0 = np.searchsorted(pf, frames - before - ext_n, "left")
        cnt = np.searchsorted(pf, frames + after - 1 + ext_n, "right") - i0
        pe = np.repeat(np.arange(frames.size), cnt)
        pj = np.repeat(i0, cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        sp, su = pf[pj], pu[pj]
        keep = (np.abs(sp - frames[pe]) > claim) & (ecol[su] >= 0)
        pe, sp, ti = pe[keep], sp[keep], ecol[su[keep]]
        model = np.zeros_like(s)
        if sp.size:
            tf = E[ti][:, ext_n + fit[0]:ext_n + fit[-1] + 1]
            alpha = np.einsum("ntc,ntc->n", take(sp, fit), tf) / np.maximum(np.einsum("ntc,ntc->n", tf, tf), 1e-12)
            j = (frames[pe] - before - sp)[:, None] + np.arange(before + after) + ext_n
            inside = (j >= 0) & (j <= 2 * ext_n)
            contrib = E[ti[:, None], np.clip(j, 0, 2 * ext_n)] * (inside * alpha[:, None])[:, :, None]
            starts = np.flatnonzero(np.r_[True, pe[1:] != pe[:-1]])
            model[pe[starts]] = np.add.reduceat(contrib, starts, axis=0)
        return s - model

    def score_sub(frames):
        res = residual(frames, take(frames, off))
        return res, cosines(res)[1]

    return take, score_sub


def _mua_window(w: int, out_dir: str) -> tuple[int, int, int, float]:
    from spikeinterface.sortingcomponents.peak_detection import detect_peaks

    t0 = time.time()
    P: AnalyzerParams = _W["params"]
    rec, groups, noise, bad = _W["rec"], _W["groups"], _W["noise"], _W["bad"]
    fs = rec.get_sampling_frequency()
    nfr = rec.get_num_frames()
    before, after = samples(P.mua_ms_before, fs), samples(P.mua_ms_after, fs)
    claim, ext_n = samples(P.mua_claim_ms, fs), samples(P.mua_extend_ms, fs)
    half = samples(P.mua_trough_ms, fs)
    fit = np.arange(samples(P.mua_fit_ms[0], fs), samples(P.mua_fit_ms[1], fs) + 1)
    win = samples(P.mua_window_s * 1000, fs)
    a, b = w * win, min((w + 1) * win, nfr)
    rng = np.random.default_rng([P.seed, w])
    ra, rb = max(0, a - 2 * ext_n), min(nfr, b + 2 * ext_n)
    x = np.ascontiguousarray(rec.get_traces(start_frame=ra, end_frame=rb, return_in_uV=False))
    da, db = max(0, a - before - after), min(nfr, b + before + after)
    peaks = detect_peaks(rec.frame_slice(da, db), method="locally_exclusive",
                         method_kwargs={"peak_sign": "neg", "detect_threshold": P.mua_detect,
                                        "radius_um": P.mua_radius_um, "noise_levels": noise},
                         job_kwargs={"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False})
    f = peaks["sample_index"].astype(np.int64) + da
    ch = peaks["channel_index"].astype(np.int64)
    amp = -peaks["amplitude"].astype(np.float64) / noise[ch]
    keep = (f >= a) & (f < b) & (f - before >= 0) & (f + after <= nfr) & (amp >= P.mua_floor)
    f, ch, amp = f[keep], ch[keep], amp[keep]
    g = groups[ch]
    out = {k: [] for k in ("frame", "group", "amp", "score", "trough", "ref_group", "ref_score")}
    off = np.arange(-before, after)
    ext = np.arange(-ext_n, ext_n + 1)
    for k in np.unique(g).tolist():
        if k not in _W["pool_f"]:
            continue
        m = np.flatnonzero((g == k) & ~in_bad(f, k, bad))
        pf, pu = _W["pool_f"][k], _W["pool_u"][k]
        claims = np.searchsorted(pf, f[m] + claim, "right") - np.searchsorted(pf, f[m] - claim, "left")
        m = m[claims <= 1]
        claims = claims[claims <= 1]
        if m.size == 0:
            continue
        cols = np.flatnonzero(groups == k)
        xt = np.ascontiguousarray(x[:, cols])
        T, E, ids = [], [], []
        ecol = np.full(len(_W["uids"]), -1)
        for ui in np.flatnonzero(_W["ugroup"] == k).tolist():
            t = _W["trains"][_W["uids"][ui]]
            t = t[(t >= a) & (t < b) & (t - before >= 0) & (t + after <= nfr)]
            if t.size == 0:
                continue
            scorable = t.size >= P.mua_template_min_spikes
            if scorable and t.size > P.mua_template_spikes:
                t = np.sort(rng.choice(t, P.mua_template_spikes, replace=False))
            e = xt[np.clip((t - ra)[:, None] + ext, 0, xt.shape[0] - 1)].astype(np.float32).mean(axis=0)
            ecol[ui] = len(E)
            E.append(e)
            if scorable:
                T.append(e[ext_n - before:ext_n + after])
                ids.append(ui)
        if not T:
            continue
        T, E = np.stack(T), np.stack(E)
        tn = np.sqrt((T ** 2).sum(axis=(1, 2)))
        _, score_sub = _residual_scorer(xt, ra, T, tn, E, ecol, pf, pu, off, fit, ext_n, before, after, claim)
        nz = noise[cols].astype(np.float32)
        score = np.empty(m.size, np.float32)
        trough = np.empty(m.size, np.float32)
        for c0 in range(0, m.size, 50_000):
            fr = f[m[c0:c0 + 50_000]]
            res, c = score_sub(fr)
            score[c0:c0 + fr.size] = c.max(axis=1)
            trough[c0:c0 + fr.size] = (-res[:, before - half:before + half + 1] / nz).max(axis=(1, 2))
        un = claims == 0
        out["frame"].append(f[m][un])
        out["group"].append(np.full(int(un.sum()), k, np.int16))
        out["amp"].append(amp[m][un].astype(np.float32))
        out["score"].append(score[un])
        out["trough"].append(trough[un])
        out["ref_group"].append(np.full(int((~un).sum()), k, np.int16))
        out["ref_score"].append(score[~un])
    dtypes = {"frame": np.int64, "group": np.int16, "amp": np.float32, "score": np.float32, "trough": np.float32,
              "ref_group": np.int16, "ref_score": np.float32}
    arrs = {k: np.concatenate(v) if v else np.zeros(0, dtypes[k]) for k, v in out.items()}
    path = pathlib.Path(out_dir) / f"w{w:03d}.npz"
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, a=a, b=b, **arrs)
    tmp.replace(path)
    return w, arrs["frame"].size, arrs["ref_score"].size, time.time() - t0


def run_multi_unit(rec_path, sorting_path, noise: np.ndarray, bad_periods: dict, out_dir,
                   params: AnalyzerParams | None = None, *, drop_channel_ids=(), workers: int = 24) -> None:
    """Score every window's unclaimed and reference events into ``out_dir/wNNN.npz`` (skipping windows
    already there). Each worker holds one window of the binary in memory. ``out_dir/noise_counts.npy``
    keeps the noise the windows were scored with; other noise raises ValueError."""
    from concurrent.futures import ProcessPoolExecutor

    params = params or AnalyzerParams()
    drop = [str(c) for c in drop_channel_ids]
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    noise_path = out_dir / "noise_counts.npy"
    if noise_path.exists():
        if not np.allclose(np.load(noise_path), np.asarray(noise, np.float64), rtol=1e-9, atol=0):
            raise ValueError(f"{out_dir} was scored with other noise levels; remove it to rescore")
    elif any(out_dir.glob("w*.npz")):
        raise ValueError(f"{out_dir} holds windows without noise_counts.npy; remove it to rescore")
    else:
        np.save(noise_path, np.asarray(noise, np.float64))
    rec = open_sorting_recording(rec_path, drop)
    win = samples(params.mua_window_s * 1000, rec.get_sampling_frequency())
    n_win = int(np.ceil(rec.get_num_frames() / win))
    todo = [w for w in range(n_win) if not (out_dir / f"w{w:03d}.npz").exists()]
    log(f"multi-unit: {n_win} windows, {len(todo)} to score, {workers} workers")
    if not todo:
        return
    bad = bad_arrays(bad_periods)
    with ProcessPoolExecutor(min(workers, len(todo)), initializer=_mua_init,
                             initargs=(str(rec_path), drop, str(sorting_path), noise, bad, params)) as ex:
        for w, nu, nr, dt in ex.map(_mua_window, todo, [str(out_dir)] * len(todo)):
            log(f"multi-unit window {w}: {nu:,} unclaimed, {nr:,} reference, {dt:.0f} s")


def multi_unit_events(out_dir, params: AnalyzerParams | None = None) -> tuple[dict, pd.DataFrame]:
    """``{group: frames}`` of the admitted events, from the window files of ``run_multi_unit``, and a
    table per tetrode of the cutoff and the counts. A tetrode with fewer than ``mua_min_reference``
    reference events gets no cutoff and no events."""
    params = params or AnalyzerParams()
    files = sorted(pathlib.Path(out_dir).glob("w[0-9][0-9][0-9].npz"))
    d = {k: [] for k in ("frame", "group", "amp", "score", "trough", "ref_group", "ref_score")}
    for f in files:
        z = np.load(f)
        for k, v in d.items():
            v.append(z[k])
    d = {k: np.concatenate(v) for k, v in d.items()}
    rows, events = [], {}
    for k in sorted(set(np.unique(d["group"]).tolist()) | set(np.unique(d["ref_group"]).tolist())):
        ref = d["ref_score"][d["ref_group"] == k]
        un = d["group"] == k
        row = {"group": k, "tetrode": k + 1, "reference": int(ref.size), "unclaimed": int(un.sum())}
        if ref.size >= params.mua_min_reference:
            cut = float(np.percentile(ref, params.mua_loss))
            ok = un & (d["score"] >= cut)
            adm = ok & (d["trough"] >= params.mua_floor)
            row.update(cutoff=cut, at_or_above_cutoff=int(ok.sum()), trough_below_floor=int((ok & ~adm).sum()),
                       admitted=int(adm.sum()))
            events[k] = np.sort(d["frame"][adm])
        else:
            row.update(cutoff=np.nan, at_or_above_cutoff=0, trough_below_floor=0, admitted=0)
        rows.append(row)
    return events, pd.DataFrame(rows)


def with_multi_unit(sorting, events: dict, *, unit_id_base: int | None = None):
    """``sorting`` plus one unit per tetrode holding ``events[group]``, with properties ``group`` and
    ``is_mua``. A multi-unit unit's id is ``unit_id_base`` plus the tetrode number (group + 1); the
    default base is the first power of ten above every unit id."""
    import spikeinterface as si

    ids = np.asarray(sorting.unit_ids)
    if not np.issubdtype(ids.dtype, np.integer):
        raise ValueError("multi-unit ids need integer unit ids")
    base = unit_id_base if unit_id_base is not None else 10 ** len(str(int(ids.max()) + 100))
    mua_groups = [k for k in sorted(events) if events[k].size]
    mua_ids = [base + k + 1 for k in mua_groups]
    if set(mua_ids) & set(ids.tolist()):
        raise ValueError("multi-unit ids collide with unit ids")
    sv = sorting.to_spike_vector()
    frames = [sv["sample_index"]] + [events[k] for k in mua_groups]
    labels = [ids[sv["unit_index"]]] + [np.full(events[k].size, i, ids.dtype) for k, i in zip(mua_groups, mua_ids)]
    fr, lab = np.concatenate(frames), np.concatenate(labels)
    o = np.argsort(fr, kind="stable")
    all_ids = np.concatenate([ids, np.asarray(mua_ids, ids.dtype)])
    out = si.NumpySorting.from_samples_and_labels([fr[o]], [lab[o]], sampling_frequency=sorting.sampling_frequency,
                                                  unit_ids=all_ids)
    group = np.concatenate([np.asarray(sorting.get_property("group")).astype(int), np.asarray(mua_groups, int)])
    out.set_property("group", group)
    out.set_property("is_mua", np.r_[np.zeros(ids.size, bool), np.ones(len(mua_ids), bool)])
    return out


# --- periods ---------------------------------------------------------------------------------------

def unit_periods(sorting, n_frames: int, bad_periods: dict) -> np.ndarray:
    """SpikeInterface unit periods: for each unit, the stretches of [0, n_frames) outside its tetrode's
    bad intervals (single segment)."""
    from spikeinterface.core.base import unit_period_dtype

    bad = bad_arrays(bad_periods)
    group = np.asarray(sorting.get_property("group")).astype(int)
    rows = []
    for ui, g in enumerate(group.tolist()):
        edges = [0]
        if g in bad:
            for s, e in zip(*bad[g]):
                edges += [max(0, int(s)), min(n_frames, int(e))]
        edges.append(n_frames)
        for s, e in zip(edges[::2], edges[1::2]):
            if e > s:
                rows.append((0, s, e, ui))
    return np.array(rows, dtype=unit_period_dtype)


# --- all-spike PC metrics --------------------------------------------------------------------------

def _sq_mahalanobis(X: np.ndarray, mu: np.ndarray, VI: np.ndarray, step: int = 2_000_000) -> np.ndarray:
    out = np.empty(len(X))
    for i in range(0, len(X), step):
        d = X[i:i + step] - mu
        out[i:i + step] = ((d @ VI) * d).sum(axis=1)
    return out


def mahalanobis_all(own: np.ndarray, other: np.ndarray) -> tuple[float, float]:
    """Isolation distance and L-ratio of ``own`` against ``other`` (rows are spikes), as SpikeInterface's
    ``mahalanobis_metrics`` defines them, but with isolation distance missing when the unit has at least
    as many spikes as the others, where it is undefined."""
    import scipy.stats

    if len(own) < 2 or len(other) < 1:
        return np.nan, np.nan
    try:
        VI = np.linalg.inv(np.cov(own.T))
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    d2 = _sq_mahalanobis(other, own.mean(0), VI)
    l_ratio = float(scipy.stats.chi2.sf(d2, own.shape[1]).sum() / len(own))
    n = len(own)
    iso = float(np.partition(d2, n - 1)[n - 1]) if n < len(other) else np.nan
    return iso, l_ratio


def d_prime_all(X: np.ndarray, is_own: np.ndarray) -> float:
    """d' of SpikeInterface's ``d_prime_metric`` (a two-class linear discriminant of the unit against the
    rest), from the discriminant's closed form: the direction is the inverse of the class covariances,
    weighted by class size, times the difference of the class means. Oriented so that d' is positive."""
    n_a = int(is_own.sum())
    n_b = len(X) - n_a
    s1, s2 = X.sum(axis=0), X.T @ X
    a = X[is_own]
    s1a, s2a = a.sum(axis=0), a.T @ a
    mu_a, mu_b = s1a / n_a, (s1 - s1a) / n_b
    cov_a = s2a / n_a - np.outer(mu_a, mu_a)
    cov_b = (s2 - s2a) / n_b - np.outer(mu_b, mu_b)
    w = np.linalg.solve((n_a * cov_a + n_b * cov_b) / (n_a + n_b), mu_a - mu_b)
    p = X @ w
    pa, pb = p[is_own], p[~is_own]
    return float((pa.mean() - pb.mean()) / np.sqrt((pa.var() + pb.var()) / 2))


def silhouette_all(X: np.ndarray, labels: np.ndarray, unit: int) -> float:
    """SpikeInterface's ``simplified_silhouette_score``: against the other cluster whose centroid is
    nearest the unit's spikes on average."""
    own = X[labels == unit]
    a = np.linalg.norm(own - own.mean(0), axis=1)
    best, b = np.inf, None
    for v in np.unique(labels):
        if v == unit:
            continue
        dist = np.linalg.norm(own - X[labels == v].mean(0), axis=1)
        if dist.mean() < best:
            best, b = dist.mean(), dist
    if b is None:
        return np.nan
    return float(np.mean((b - a) / np.maximum(b, a)))


def all_spike_pc_metrics(sorting, recording, scratch_dir, params: AnalyzerParams | None = None, *,
                         n_jobs: int = 48, keep_projections: bool = False) -> pd.DataFrame:
    """PC metrics on every spike's projection. An analyzer in ADC counts (in memory) fits the PCA on up
    to ``max_spikes_per_unit`` random spikes per unit and projects every spike (``run_for_all_spikes``
    reads ADC counts whatever ``return_in_uV`` says). Per tetrode, each unit is measured against every
    other spike of the tetrode; units with ``is_mua`` count as other spikes and get no metrics, nor do
    units with fewer than two spikes or fewer than two other spikes on their tetrode."""
    params = params or AnalyzerParams()
    from ttm.sorting import create_group_sparse_analyzer

    scratch_dir = pathlib.Path(scratch_dir)
    scratch_dir.mkdir(parents=True, exist_ok=True)
    az = create_group_sparse_analyzer(sorting, recording, format="memory", return_in_uV=False, n_jobs=n_jobs,
                                      seed=params.seed)
    az.compute({"random_spikes": {"method": "uniform", "max_spikes_per_unit": params.max_spikes_per_unit,
                                  "seed": params.seed},
                "waveforms": {"ms_before": params.ms_before, "ms_after": params.ms_after},
                "principal_components": {"n_components": params.n_pcs, "mode": params.pc_mode}},
               n_jobs=n_jobs, chunk_duration="1s", progress_bar=False)
    path = scratch_dir / "all_pcs.npy"
    t0 = time.time()
    az.get_extension("principal_components").run_for_all_spikes(file_path=path, n_jobs=n_jobs,
                                                                chunk_duration="1s", progress_bar=False)
    log(f"all-spike projections in {(time.time() - t0) / 60:.1f} min -> {path}")
    pc = pc_metrics_from_projections(np.load(path, mmap_mode="r"), sorting, az.sparsity.mask.sum(axis=1))
    if not keep_projections:
        path.unlink()
        if not any(scratch_dir.iterdir()):
            scratch_dir.rmdir()
    return pc


def pc_metrics_from_projections(P: np.ndarray, sorting, n_channels: np.ndarray) -> pd.DataFrame:
    """The metrics of ``all_spike_pc_metrics`` from projections ``P`` (spike, component, channel) in the
    order of ``sorting``'s spike vector, a unit's first ``n_channels[unit]`` channels being its own."""
    labels = sorting.to_spike_vector()["unit_index"]
    group = np.asarray(sorting.get_property("group")).astype(int)
    keys = sorting.get_property_keys()
    is_mua = np.asarray(sorting.get_property("is_mua")).astype(bool) if "is_mua" in keys else np.zeros(group.size, bool)
    rows = []
    for g in np.unique(group).tolist():
        units = np.flatnonzero(group == g)
        nch = int(n_channels[units[0]])
        sel = np.flatnonzero(np.isin(labels, units))
        X = np.asarray(P[sel][:, :, :nch]).reshape(sel.size, nch * P.shape[1]).astype(np.float64)
        lab = labels[sel]
        for u in units.tolist():
            row = {"unit_id": sorting.unit_ids[u], **{m: np.nan for m in PC_METRIC_NAMES}}
            own = lab == u
            if not is_mua[u] and own.sum() >= 2 and (~own).sum() >= 2:
                row["isolation_distance"], row["l_ratio"] = mahalanobis_all(X[own], X[~own])
                row["d_prime"] = d_prime_all(X, own)
                row["silhouette"] = silhouette_all(X, lab, u)
            rows.append(row)
        log(f"PC metrics on group {g}: {len(units)} units, {sel.size:,} spikes")
    return pd.DataFrame(rows).set_index("unit_id")


# --- the analyzer ----------------------------------------------------------------------------------

def set_sd_ratio_noise(recording, noise_std_uv: np.ndarray) -> None:
    """Make SpikeInterface's sd_ratio divide by ``noise_std_uv`` (per channel, µV). Its
    ``compute_sd_ratio`` takes the noise from ``get_noise_levels(method="std")``, which returns the
    recording property ``noise_level_std_scaled`` when it is set and otherwise draws 20 random chunks
    with no fixed seed and keeps them for the process."""
    recording.set_property("noise_level_std_scaled", np.asarray(noise_std_uv, np.float64))


def build_analyzer(sorting, recording, folder, noise_uv: np.ndarray, periods: np.ndarray, pc_metrics: pd.DataFrame,
                   params: AnalyzerParams | None = None, *, noise_std_uv: np.ndarray | None = None,
                   n_jobs: int = 48, note: str = ""):
    """The SortingAnalyzer in ``folder`` (zarr, µV, every unit sparse to its tetrode) with random spikes,
    waveforms, templates, the given noise levels, spike amplitudes, correlograms, ISI histograms, unit
    locations, template similarity, template metrics, principal components (on the random spikes, for
    viewing) and quality metrics over ``periods``, the all-spike PC metrics added to them. sd_ratio, if
    listed, divides by ``noise_std_uv`` (``set_sd_ratio_noise``). Everything is
    computed in memory and written to zarr once: SpikeInterface 0.105.4 appends each chunk's waveforms
    and amplitudes to a zarr analyzer's arrays one at a time, rewriting the open zarr chunk each time."""
    params = params or AnalyzerParams()
    from ttm.sorting import create_group_sparse_analyzer

    if "sd_ratio" in params.metric_names and noise_std_uv is None:
        raise ValueError("sd_ratio needs noise_std_uv")
    job = {"n_jobs": n_jobs, "chunk_duration": "1s", "progress_bar": False}
    az = create_group_sparse_analyzer(sorting, recording, format="memory", return_in_uV=True, n_jobs=n_jobs,
                                      seed=params.seed)
    steps = [
        ("random_spikes", {"method": "uniform", "max_spikes_per_unit": params.max_spikes_per_unit,
                           "seed": params.seed}),
        ("waveforms", {"ms_before": params.ms_before, "ms_after": params.ms_after}),
        ("templates", {"operators": ["average", "std", "median"]}),
        ("noise_levels", {"random_slices_kwargs": {"num_chunks_per_segment": 1, "seed": params.seed}}),
    ]
    for name, kw in steps:
        t0 = time.time()
        az.compute(name, **kw, **job)
        log(f"analyzer {name} in {time.time() - t0:.0f} s")
    ext = az.get_extension("noise_levels")
    ext.set_data("noise_levels", np.asarray(noise_uv, np.float32))
    ext.params = {**ext.params, "ttm": note or "set by ttm.analyzer.chunk_noise_levels"}
    ext.save()
    for name, kw in [
        ("spike_amplitudes", {}),
        ("correlograms", {"window_ms": params.correlogram_window_ms, "bin_ms": params.correlogram_bin_ms}),
        ("isi_histograms", {"window_ms": 50.0, "bin_ms": 1.0}),
        ("unit_locations", {"method": "monopolar_triangulation"}),
        ("template_similarity", {"method": "cosine"}),
        ("template_metrics", {"include_multi_channel_metrics": True}),
        ("principal_components", {"n_components": params.n_pcs, "mode": params.pc_mode}),
    ]:
        t0 = time.time()
        az.compute(name, **kw, **job)
        log(f"analyzer {name} in {time.time() - t0:.0f} s")
    t0 = time.time()
    if noise_std_uv is not None:
        set_sd_ratio_noise(az.recording, noise_std_uv)
    az.compute("quality_metrics", metric_names=list(params.metric_names), periods=periods, **job)
    log(f"analyzer quality_metrics in {time.time() - t0:.0f} s")
    qm = az.get_extension("quality_metrics")
    df = qm.get_data().copy()
    for m in PC_METRIC_NAMES:
        df[m] = pc_metrics.reindex(df.index)[m].to_numpy()
    qm.set_data("metrics", df)
    qm.params = {**qm.params, "ttm_all_spike_pc_metrics": list(PC_METRIC_NAMES)}
    if noise_std_uv is not None:
        qm.params["ttm_sd_ratio_noise"] = "std_uv of ttm.analyzer.chunk_noise_levels"
    qm.save()
    t0 = time.time()
    az = az.save_as(format="zarr", folder=folder)
    log(f"analyzer written in {time.time() - t0:.0f} s")
    return az


def write_params(path, params: AnalyzerParams, **inputs) -> None:
    """``analyzer_params.json``: the parameters, inputs and package versions."""
    import importlib.metadata as md

    versions = {}
    for name in ("ttm", "spikeinterface", "numpy", "scikit-learn"):
        try:
            versions[name] = md.version(name)
        except md.PackageNotFoundError:
            versions[name] = None
    pathlib.Path(path).write_text(json.dumps({"params": params.to_dict(), "inputs": {k: str(v) for k, v in inputs.items()},
                                              "versions": versions, "written": time.strftime("%Y-%m-%dT%H:%M:%S")},
                                             indent=2) + "\n")
