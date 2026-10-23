"""Purity numbers and pages for every unit of a long tetrode sorting, and pairing with another.

Per 1800 s window and tetrode the preprocessed binary is read once and the
spike-excised noise estimated (``ttm.capture.noise_mad_excised``); every
spike of every unit of the tetrode gets its trough depth in noise units, its
cosine to the unit's window-local template and its residual energy: the
mean over the four-wire snippet of ``(snippet - a * template)^2 / noise^2``
with ``a`` the spike's least-squares amplitude scale, so a spike that is the
template plus noise scores about 1 and extra shape scores more. From the whole
trains ``ttm.capture.unit_purity`` gives per unit the refractory
contamination, the amplitude-split ratios and the cross-correlogram verdict
against every same-tetrode unit. ``unit_page`` draws one page per unit
(waveform overlays, amplitude by shape, amplitude over time, ISI, the
amplitude split and the neighbour cross-correlograms). ``match_units`` pairs
the units of two sortings of the same recording by template cosine on the
same tetrode and measures spike agreement.

Nothing here reads a manifest: callers pass the binary the sorting was
produced from and, for sortings that cover part of the recording, the frame
range and the duration the rates are computed over.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

from ttm import capture as cp

FS = 30000.0
WIN_S = 1800.0
MAX_POINTS = 20_000   # scatter points per page; a 48 h unit has millions of spikes


def residual_energy(traces: np.ndarray, frames: np.ndarray, template: np.ndarray, noise: np.ndarray, *,
                    s_max: int = 1) -> np.ndarray:
    """Per event: the mean over the four-wire window of ``((snippet - a * template) / noise)^2``,
    ``a`` the least-squares scale, minimised over shifts ``-s_max..s_max`` (the shifts
    ``cp.cosine_to_templates`` maximises over). About 1 for the template plus noise, more for any
    other shape or an overlapping spike; NaN where the shifted window leaves the array."""
    nbefore, nafter = cp.NBEFORE, template.shape[0] - cp.NBEFORE
    tsq = max(float((template * template).sum()), 1e-12)
    out = np.full(frames.size, np.nan)
    ok = (frames - nbefore - s_max >= 0) & (frames + nafter + s_max <= traces.shape[0])
    if not ok.any():
        return out
    best = np.full(int(ok.sum()), np.inf)
    for s in range(-s_max, s_max + 1):
        sn = cp.snippets(traces, frames[ok] + s, nbefore=nbefore, nafter=nafter)
        a = np.einsum("ntc,tc->n", sn, template) / tsq
        e = (((sn - a[:, None, None] * template[None]) / noise[None, None, :]) ** 2).mean(axis=(1, 2))
        best = np.minimum(best, e)
    out[ok] = best
    return out


def window_scores(binary: str, group: int, a: int, b: int, trains: dict) -> dict:
    """Per unit with spikes in frames [a, b) of tetrode ``group``: (absolute frames, trough
    amplitude in excised-MAD units, cosine to the window template, residual energy in noise units).
    Units with fewer than 20 spikes in the window get no template unless no unit has 20; their
    cosine and residual are NaN."""
    import spikeinterface as si

    rec = si.load(binary)
    x = cp.load_tetrode_traces(rec, group, a, b)
    noise, _ = cp.noise_mad_excised(x)
    local = {u: fr[(fr >= a + cp.NBEFORE) & (fr + cp.NAFTER <= b)] - a for u, fr in trains.items()}
    local = {u: fr for u, fr in local.items() if fr.size}
    out = {"group": group, "start": a, "noise": noise, "units": {}}
    if not local:
        return out
    T, ids = cp.templates_from_trains(x, {u: fr for u, fr in local.items() if fr.size >= 20} or local)
    for u, fr in local.items():
        sn = cp.snippets(x, fr)
        amp = -(sn / noise[None, None, :]).min(axis=(1, 2))
        if u in ids:
            t = T[ids.index(u)]
            r = cp.cosine_to_templates(x, fr, t[None], s_max=1)[:, 0]
            resid = residual_energy(x, fr, t, noise, s_max=1)
        else:
            r = np.full(fr.size, np.nan)
            resid = np.full(fr.size, np.nan)
        out["units"][u] = (fr + a, amp, r, resid)
    return out


def score_sorting(binary: str, sorting, workers: int, cache: pathlib.Path | None = None,
                  frame_range: tuple[int, int] | None = None) -> tuple[dict, dict, dict, dict, dict]:
    """Per unit: frames, trough amplitudes (MAD), cosines and residual energies of every spike, and
    the unit's tetrode. ``frame_range`` limits the windows read to those overlapping it (window
    edges stay on the 1800 s grid from frame 0). Cached in ``cache`` (npz); a cache without
    residuals is recomputed."""
    import spikeinterface as si
    from joblib import Parallel, delayed

    units = [int(u) for u in sorting.unit_ids]
    groups = np.asarray(sorting.get_property("group")).astype(int)
    ugroup = dict(zip(units, groups.tolist()))
    if cache is not None and cache.exists():
        z = np.load(cache)
        if all(f"resid_{u}" in z for u in units):
            return ({u: z[f"frames_{u}"] for u in units}, {u: z[f"amps_{u}"] for u in units},
                    {u: z[f"rs_{u}"] for u in units}, {u: z[f"resid_{u}"] for u in units}, ugroup)
        print(f"{cache.name} has no residual energies; scoring again", flush=True)
    nfr = si.load(binary).get_num_frames()
    trains = {u: np.sort(sorting.get_unit_spike_train(u).astype(np.int64)) for u in units}
    wf = int(WIN_S * FS)
    lo, hi = frame_range if frame_range is not None else (0, nfr)
    tasks = [(g, s, min(s + wf, nfr)) for s in range(0, nfr, wf) if s < hi and min(s + wf, nfr) > lo
             for g in sorted(set(groups.tolist()))]
    res = Parallel(n_jobs=workers, batch_size=1)(
        delayed(window_scores)(binary, g, s, e, {u: trains[u] for u in units if ugroup[u] == g})
        for g, s, e in tasks)
    acc = {k: {u: [] for u in units} for k in ("frames", "amps", "rs", "resid")}
    for r in res:
        for u, vals in r["units"].items():
            for k, v in zip(("frames", "amps", "rs", "resid"), vals):
                acc[k][u].append(v)

    def cat(k, u, dtype):
        return np.concatenate(acc[k][u]) if acc[k][u] else np.zeros(0, dtype)

    frames = {u: cat("frames", u, np.int64) for u in units}
    amps = {u: cat("amps", u, float) for u in units}
    rs = {u: cat("rs", u, float) for u in units}
    resid = {u: cat("resid", u, float) for u in units}
    if cache is not None:
        np.savez(cache, **{f"frames_{u}": frames[u] for u in units}, **{f"amps_{u}": amps[u] for u in units},
                 **{f"rs_{u}": rs[u] for u in units}, **{f"resid_{u}": resid[u] for u in units})
    return frames, amps, rs, resid, ugroup


def frac_within(ref: np.ndarray, other: np.ndarray, tol: int = 15) -> float:
    """Fraction of ``ref`` frames with a frame of ``other`` (sorted) within ``tol`` frames."""
    if ref.size == 0 or other.size == 0:
        return 0.0
    j = np.searchsorted(other, ref)
    lo = other[np.clip(j - 1, 0, other.size - 1)]
    hi = other[np.clip(j, 0, other.size - 1)]
    return float((np.minimum(np.abs(ref - lo), np.abs(ref - hi)) <= tol).mean())


def _r(x, nd):
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


def purity_table(frames: dict, amps: dict, ugroup: dict, n_frames: int, dead_ms: float,
                 rs: dict | None = None, resid: dict | None = None) -> pd.DataFrame:
    """One row per unit. ``n_frames`` is the duration the sorting covers (rates, refractory
    contamination and the amplitude-split expectation are computed over it). With ``rs`` and
    ``resid`` the within-cluster spread columns are added: median and 10th percentile cosine,
    median residual energy, and the coefficient of variation of the trough amplitude."""
    rows = []
    for g in sorted(set(ugroup.values())):
        members = [u for u in frames if ugroup[u] == g]
        pur = cp.unit_purity({u: frames[u] for u in members}, {u: amps[u] for u in members}, n_frames,
                             dead_ms=dead_ms)
        for u in members:
            d = pur[u]
            verdicts = list(d["neighbours"].values())
            row = {"unit": u, "tetrode": f"TT{g + 1}", "n_spikes": d["n_spikes"],
                   "rate_hz": round(d["n_spikes"] / (n_frames / FS), 3),
                   "amp_mad_median": _r(np.median(amps[u]), 2) if amps[u].size else None,
                   "rp_contamination": _r(d["rp_contamination"], 4),
                   "rp_violations": d["rp_violations"],
                   "ampsplit_ratio": _r(d["ampsplit_ratio"], 3),
                   "ampsplit_burst_ratio": _r(d.get("ampsplit_burst_ratio"), 3),
                   "n_duplicate": verdicts.count("duplicate"), "n_ambiguous": verdicts.count("ambiguous"),
                   "n_distinct": verdicts.count("distinct"), "n_segregated": verdicts.count("segregated"),
                   "neighbours": " ".join(f"{v}:{d['neighbours'][v][:3]}" for v in sorted(d["neighbours"]))}
            if rs is not None and resid is not None:
                ok = np.isfinite(rs[u])
                row.update({
                    "n_scored": int(ok.sum()),
                    "cos_median": _r(np.median(rs[u][ok]), 4) if ok.any() else None,
                    "cos_p10": _r(np.quantile(rs[u][ok], 0.1), 4) if ok.any() else None,
                    "resid_median": _r(np.median(resid[u][ok]), 3) if ok.any() else None,
                    "amp_cv": _r(np.std(amps[u]) / np.mean(amps[u]), 3) if amps[u].size > 1 else None})
            rows.append(row)
    return pd.DataFrame(rows)


def _window_snippets(binary: str, group: int, a: int, b: int, trains: dict, per_unit: int, pad: int) -> dict:
    import spikeinterface as si

    rec = si.load(binary)
    x = cp.load_tetrode_traces(rec, group, a, b)
    out = {}
    rng = np.random.default_rng(a)
    for u, fr in trains.items():
        loc = fr[(fr >= a + cp.NBEFORE + pad) & (fr + cp.NAFTER + pad < b)] - a
        if loc.size == 0:
            continue
        if loc.size > per_unit:
            loc = np.sort(rng.choice(loc, per_unit, replace=False))
        out[u] = cp.snippets(x, loc, nbefore=cp.NBEFORE + pad, nafter=cp.NAFTER + pad)
    return out


def unit_templates(binary: str, trains: dict, ugroup: dict, nfr: int, *, n_windows: int = 8, win_s: float = 60.0,
                   max_spikes: int = 1000, workers: int = 16,
                   frame_range: tuple[int, int] | None = None) -> tuple[dict, dict]:
    """Four-wire templates of every unit from ``n_windows`` evenly spaced windows of ``win_s`` inside
    ``frame_range`` (default the whole recording; single reads for sparse units), aligned at the
    trough. Returns ({unit: (90, 4)}, {unit: n spikes used})."""
    import spikeinterface as si
    from joblib import Parallel, delayed

    pad = 10
    wf = int(win_s * FS)
    lo, hi = frame_range if frame_range is not None else (0, nfr)
    starts = np.linspace(lo, hi - wf, n_windows + 2)[1:-1].astype(int)
    groups = sorted(set(ugroup.values()))
    tasks = [(g, int(s)) for s in starts for g in groups]
    per_unit = max(max_spikes // n_windows, 20)
    res = Parallel(n_jobs=workers)(delayed(_window_snippets)(binary, g, s, s + wf,
                                                             {u: trains[u] for u in trains if ugroup[u] == g},
                                                             per_unit, pad) for g, s in tasks)
    acc = {u: [] for u in trains}
    for r in res:
        for u, sn in r.items():
            acc[u].append(sn)
    rec = si.load(binary)
    templates, n_used = {}, {}
    rng = np.random.default_rng(0)
    for u, tr in trains.items():
        sn = np.concatenate(acc[u]) if acc[u] else np.zeros((0, cp.N_SAMP + 2 * pad, 4), np.float32)
        if sn.shape[0] < 50 and tr.size:
            fr = tr[(tr >= cp.NBEFORE + pad) & (tr + cp.NAFTER + pad < nfr)]
            pick = np.sort(rng.choice(fr, min(200, fr.size), replace=False)) if fr.size else fr
            idx = cp.tetrode_channel_indices(rec, ugroup[u])
            ids = list(np.asarray(rec.channel_ids)[idx])
            extra = [rec.get_traces(start_frame=f - cp.NBEFORE - pad, end_frame=f + cp.NAFTER + pad,
                                    channel_ids=ids).astype(np.float32) for f in pick]
            if extra:
                sn = np.concatenate([sn, np.stack(extra)])
        n_used[u] = int(sn.shape[0])
        if sn.shape[0] == 0:
            templates[u] = np.zeros((cp.N_SAMP, 4), np.float32)
            continue
        wide = sn.mean(axis=0)
        k = int(np.argmin(wide[pad:pad + cp.N_SAMP].min(axis=1)))
        shift = int(np.clip(k - cp.NBEFORE, -pad, pad))
        templates[u] = wide[pad + shift:pad + shift + cp.N_SAMP]
    return templates, n_used


def match_units(names: tuple[str, str], frames: tuple[dict, dict], ugroups: tuple[dict, dict],
                templates: tuple[dict, dict], units: tuple[list, list]) -> pd.DataFrame:
    """Pair each unit of the first sorting with the same-tetrode unit of the second whose template
    is most similar (``cp.template_cosine``); ``agreement`` is the fraction of the smaller train's
    spikes within 0.5 ms of a spike of the other, ``in_<name>`` the fraction of one train found in
    the other. Second-sorting units no first unit chose get a row of their own."""
    (n1, n2), (f1, f2), (g1, g2), (T1, T2), (s1, s2) = names, frames, ugroups, templates, units
    rows, used = [], set()
    for u in s1:
        cands = [v for v in s2 if g2[v] == g1[u]]
        if not cands:
            rows.append({n1: u, n2: None, "cosine": None, "agreement": 0.0, f"n_{n1}": int(f1[u].size)})
            continue
        cos = np.array([cp.template_cosine(T1[u], T2[v]) for v in cands])
        v = cands[int(np.argmax(cos))]
        used.add(v)
        small, big = (f1[u], f2[v]) if f1[u].size <= f2[v].size else (f2[v], f1[u])
        rows.append({n1: u, n2: v, "cosine": round(float(cos.max()), 4), "agreement": round(frac_within(small, big), 4),
                     f"n_{n1}": int(f1[u].size), f"n_{n2}": int(f2[v].size),
                     f"in_{n2}": round(frac_within(f1[u], f2[v]), 4), f"in_{n1}": round(frac_within(f2[v], f1[u]), 4)})
    for v in s2:
        if v not in used:
            rows.append({n1: None, n2: v, "cosine": None, "agreement": 0.0, f"n_{n2}": int(f2[v].size)})
    return pd.DataFrame(rows)


def unit_page(axes, binary: str, u: int, frames: dict, amps: dict, rs: dict, ugroup: dict, templates: dict,
              gain: float, nfr: int, dead_ms: float, title: str, duration_frames: int | None = None) -> None:
    """One purity page into a (2, 3) axes grid from whole-recording trains; the overlays come from
    the 1800 s window with most spikes. ``nfr`` is the recording's length (window bounds),
    ``duration_frames`` the duration the statistics use (default ``nfr``)."""
    import spikeinterface as si

    dur = duration_frames or nfr
    fr, am, rr = frames[u], amps[u], rs[u]
    g = ugroup[u]
    wf = int(WIN_S * FS)
    tt = (np.arange(cp.N_SAMP) - cp.NBEFORE) / FS * 1000
    ax = axes[0, 0]
    if fr.size:
        w = int(np.bincount(fr // wf).argmax())
        a_, b_ = w * wf, min((w + 1) * wf, nfr)
        rec = si.load(binary)
        x = cp.load_tetrode_traces(rec, g, a_, b_)
        loc = fr[(fr >= a_ + cp.NBEFORE) & (fr + cp.NAFTER <= b_)] - a_
        rng = np.random.default_rng(0)
        pick = np.sort(rng.choice(loc, min(120, loc.size), replace=False)) if loc.size else loc
        if pick.size:
            sn = cp.snippets(x, pick) * gain
            for k in range(pick.size):
                for c in range(4):
                    ax.plot(tt + c * 3.2, sn[k, :, c], color="0.6", lw=0.3, alpha=0.5)
        ax.set_title(f"{title} unit {u} (TT{g + 1}): {fr.size} spikes; {pick.size} overlaid from the "
                     f"{a_ / FS / 3600:.1f} h window", fontsize=8)
    for c in range(4):
        ax.plot(tt + c * 3.2, templates[u][:, c] * gain, color="k", lw=1.2)
    ax.set_xlabel("ms, wires offset by 3.2 ms"); ax.set_ylabel("uV")
    sub = (np.sort(np.random.default_rng(1).choice(fr.size, min(fr.size, MAX_POINTS), replace=False))
           if fr.size else np.zeros(0, int))
    ax = axes[0, 1]
    ax.scatter(rr[sub], am[sub], s=2, c="0.4", alpha=0.4, rasterized=True)
    ax.set_xlabel("cosine to the unit's window template"); ax.set_ylabel("trough amplitude (MAD)")
    ax.set_title(f"amplitude by shape ({sub.size:,} of {fr.size:,} spikes)", fontsize=8)
    ax = axes[0, 2]
    ax.scatter(fr[sub] / FS / 3600, am[sub], s=1, c="0.4", alpha=0.4, rasterized=True)
    ax.set_xlabel("session hour"); ax.set_ylabel("trough amplitude (MAD)"); ax.set_title("amplitude over time", fontsize=8)
    ax = axes[1, 0]
    isi = np.diff(fr) / FS * 1000
    c, nv = cp.rp_contamination(fr, dur)
    ax.hist(isi[isi < 50], bins=np.arange(0, 50, 0.25), color="k")
    ax.axvline(1.0, color="r", lw=0.6)
    ax.set_title(f"ISI; rp contamination {c:.3f} ({nv} < 1 ms)", fontsize=8); ax.set_xlabel("ms")
    ax = axes[1, 1]
    sp = cp.amplitude_split(fr, am, dur, dead_ms=dead_ms)
    if "lag_ms" in sp:
        ax.bar(sp["lag_ms"], sp["counts"], width=0.5, color="k")
        ax.axhline(sp["expected_per_bin"], color="gray", lw=0.8)
        ax.set_title(f"large (>= {sp['large_min']:.1f} MAD, n={sp['n_large']}) vs small (<= {sp['small_max']:.1f} MAD, "
                     f"n={sp['n_small']}): zero-lag ratio {sp['ratio']:.2f}, 2-3 ms {sp['burst_2_3ms_ratio']:.2f}",
                     fontsize=7)
        ax.set_xlabel("ms (small relative to large)")
    ax = axes[1, 2]
    others = [v for v in sorted(frames) if v != u and ugroup[v] == g]
    max_lag = int(25 * FS / 1000)
    for k, v in enumerate(others):
        lags = cp.ccg_lags(fr, frames[v], max_lag) / FS * 1000
        h, e = np.histogram(lags, np.arange(-25, 25.5, 0.5))
        verdict = cp.ccg_verdict(fr, frames[v], dead_ms=dead_ms)
        base = h.mean() if h.size else 1
        lab = f"unit {v}: {verdict['verdict']}" + (f" ({verdict['ratio']:.2f})" if np.isfinite(verdict["ratio"]) else "")
        ax.plot(e[:-1] + 0.25, h / max(base, 1e-9) + 1.5 * k, lw=0.7, label=lab)
    ax.set_title("cross-correlograms with same-tetrode units (stacked, normalised)", fontsize=8)
    ax.set_xlabel("ms")
    if others:
        ax.legend(fontsize=6)


def write_pages(path: pathlib.Path, units, binary: str, frames: dict, amps: dict, rs: dict, ugroup: dict,
                templates: dict, gain: float, nfr: int, dead_ms: float, title: str,
                duration_frames: int | None = None) -> int:
    """A PDF with one ``unit_page`` per unit (by tetrode, then id); units without spikes are
    skipped. Returns the number of pages."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    n = 0
    with PdfPages(path) as pdf:
        for u in sorted(units, key=lambda u: (ugroup[u], u)):
            if frames[u].size == 0:
                continue
            fig, axes = plt.subplots(2, 3, figsize=(20, 10))
            unit_page(axes, binary, u, frames, amps, rs, ugroup, templates, gain, nfr, dead_ms, title,
                      duration_frames=duration_frames)
            fig.tight_layout()
            pdf.savefig(fig, dpi=110)
            plt.close(fig)
            n += 1
    return n
