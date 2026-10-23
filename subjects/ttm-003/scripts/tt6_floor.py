"""TT6 cell: what a 3.5 against a 4.5 excised-MAD detection floor admits, window by window.

For one tetrode and one OPTO2 capture template (the TT6 cell, unit 22 by
default) the capture of ``opto_capture.py`` is run twice on each window of
the experiment's preprocessed binary, at the two detection floors with
every other setting equal (OPTO2 templates, excised noise, cosine gate,
peeling passes, refractory rule). The spikes the lower run assigns to the
unit are split at the higher floor into the class only the lower floor
admits (trough between the floors, "low") and the class both admit
("high"), and the two classes are compared: trough-amplitude histogram,
waveform overlays and normalised mean shapes, their cross-correlogram and
each class's autocorrelogram before the refractory rule (the rule hollows
the centre by construction, so the events it removed are put back for
these pictures), the laser PSTH of each class where the window holds a
delivered train, the firing rate of each class over the window, and the
capture's own purity pictures (neighbour cross-correlograms with verdicts,
refractory conflict ratio, amplitude-split) at each floor.

Windows: ``full`` (the whole recording), ``b<N>`` (``--window-s`` seconds
centred on TTL train N), ``h<X>`` (starting at session hour X) or
``s<A>-<B>`` (seconds). Writes ``results/tt6_floor_<E>.csv`` (one row per
window), ``figures/tt6_floor/<E>_<window>.png``,
``figures/tt6_floor/<E>_<window>_purity.pdf`` (one page per floor) and
``figures/tt6_floor/summary.png`` over every ``results/tt6_floor_*.csv``
present.

Usage: python subjects/ttm-003/scripts/tt6_floor.py --subject TTM-003 --experiment TTM-003-OPTO2 --windows full
       python subjects/ttm-003/scripts/tt6_floor.py --subject TTM-003 --experiment TTM-003-48H --windows b1 b4 b26 b33 b35 b41
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered
from ttm.subject import find_subject

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
DEAD_MS = 9 / FS * 1000
BIN = 0.005
COL = {"low": "tab:red", "high": "tab:blue", "run_hi": "0.3"}


def delivered_pulses(session: str, pattern: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pulses, trains) of the delivered trains: for a session with a block table its delivered
    blocks of 10 pulses or more, otherwise the trains of ``pattern`` at the highest power."""
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    pulses = pulses[~pulses.glitch]
    blocks = RESULTS / f"stimulation_blocks_{session}.csv"
    if blocks.exists():
        tr = pd.read_csv(blocks)
        tr = tr[tr.delivered & (tr.n_pulses >= 10)]
    else:
        tr = pd.read_csv(RESULTS / "stimulation_trains.csv")
        tr = tr[(tr.session == session) & (tr.ttl_pattern == pattern)]
        tr = tr[tr.power_mW == tr.power_mW.max()]
    tr = tr.reset_index(drop=True)
    return pulses[pulses.train.isin(tr.ttl_train)].reset_index(drop=True), tr


def all_trains(session: str) -> pd.DataFrame:
    blocks = RESULTS / f"stimulation_blocks_{session}.csv"
    if blocks.exists():
        return pd.read_csv(blocks)
    tr = pd.read_csv(RESULTS / "stimulation_trains.csv")
    return tr[tr.session == session].reset_index(drop=True)


def window_bounds(spec: str, trains: pd.DataFrame, nfr: int, win_s: float) -> tuple[int, int]:
    if spec == "full":
        return 0, nfr
    if spec.startswith("b"):
        row = trains[trains.ttl_train == int(spec[1:])].iloc[0]
        mid = (row.onset_s + row.offset_s) / 2
        start = mid - win_s / 2
    elif spec.startswith("h"):
        start = float(spec[1:]) * 3600
    elif spec.startswith("s"):
        lo, hi = spec[1:].split("-")
        return max(int(float(lo) * FS), 0), min(int(float(hi) * FS), nfr)
    else:
        raise ValueError(spec)
    a = max(int(start * FS), 0)
    return a, min(a + int(win_s * FS), nfr)


def templates_for(tt: str, d: pathlib.Path, wires_here: list[str]) -> tuple[np.ndarray, list]:
    z = np.load(d / f"templates_{tt}.npz")
    order = [[str(w) for w in z["wires"]].index(w) for w in wires_here]
    return z["templates"][:, :, order], [int(u) for u in z["unit_ids"]]


def frac_matched(ref: np.ndarray, other: np.ndarray, tol: int = 15) -> float:
    if ref.size == 0 or other.size == 0:
        return float("nan")
    j = np.searchsorted(other, ref)
    lo = other[np.clip(j - 1, 0, other.size - 1)]
    hi = other[np.clip(j, 0, other.size - 1)]
    return float((np.minimum(np.abs(ref - lo), np.abs(ref - hi)) <= tol).mean())


def rate_in(ts: np.ndarray, on: np.ndarray, lo: float, hi: float) -> float:
    if len(on) == 0:
        return float("nan")
    i0 = np.searchsorted(ts, on + lo); i1 = np.searchsorted(ts, on + hi)
    return float((i1 - i0).sum() / (len(on) * (hi - lo)))


def responses(ts: np.ndarray, on: np.ndarray, dur: float) -> dict:
    return {"base_hz": rate_in(ts, on, -0.05, 0.0), "onset_hz": rate_in(ts, on, 0.001, 0.010),
            "in_pulse_hz": rate_in(ts, on, 0.010, dur), "post_hz": rate_in(ts, on + dur, 0.0, 0.05)}


def psth(ts: np.ndarray, on: np.ndarray, dur: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(-0.1, dur + 0.3 + BIN / 2, BIN)
    i0 = np.searchsorted(ts, on + edges[0]); i1 = np.searchsorted(ts, on + edges[-1])
    rel = np.concatenate([ts[s:e] - t for s, e, t in zip(i0, i1, on)]) if len(on) else np.zeros(0)
    return edges[:-1] + BIN / 2, np.histogram(rel, edges)[0] / (max(len(on), 1) * BIN)


def ccg_hist(a: np.ndarray, b: np.ndarray, window_ms: float = 25.0, bin_ms: float = 0.5):
    lags = cp.ccg_lags(np.sort(a), np.sort(b), int(window_ms * FS / 1000)) / FS * 1000
    edges = np.arange(-window_ms, window_ms + bin_ms, bin_ms) - bin_ms / 2
    h, _ = np.histogram(lags, edges)
    return edges[:-1] + bin_ms / 2, h


def burst_ratio(centers: np.ndarray, h: np.ndarray) -> float:
    flank = h[(np.abs(centers) >= 5) & (np.abs(centers) <= 25)].mean()
    core = h[(centers > 1.75) & (centers < 3.25)].mean()
    return float(core / flank) if flank > 0 else float("nan")


def neighbour_panel(ax, C: cp.Capture, ui: int, ids: list, title: str) -> dict:
    trains = C.trains()
    fr = trains[ui]
    out = {}
    max_lag = int(25 * FS / 1000)
    for k, v in enumerate(j for j in range(len(ids)) if j != ui):
        lags = cp.ccg_lags(fr, trains[v], max_lag) / FS * 1000
        h, e = np.histogram(lags, np.arange(-25, 25.5, 0.5))
        verdict = cp.ccg_verdict(fr, trains[v], dead_ms=DEAD_MS)
        out[ids[v]] = verdict["verdict"]
        lab = f"unit {ids[v]}: {verdict['verdict']}" + (f" ({verdict['ratio']:.2f})" if np.isfinite(verdict["ratio"]) else "")
        ax.plot(e[:-1] + 0.25, h / max(h.mean(), 1e-9) + 1.5 * k, lw=0.7, label=lab)
    ax.set_title(title, fontsize=8); ax.set_xlabel("ms"); ax.legend(fontsize=5, ncol=2)
    return out


def run_window(binary: str, group: int, a: int, b: int, T: np.ndarray, ids: list, unit: int, floors: tuple,
               split: float, pulses: pd.DataFrame, trains: pd.DataFrame, gate_r: float, n_pass: int,
               refr_frames: int, figdir: pathlib.Path, exp_name: str, spec: str, make_figs: bool) -> dict:
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    rec = si.load(binary)
    gain = float(rec.get_property("gain_to_uV")[0])
    x = cp.load_tetrode_traces(rec, group, a, b)
    n = b - a
    ui = ids.index(unit)
    lo_f, hi_f = floors
    C, pur = {}, {}
    for f in floors:
        C[f] = cp.capture_tetrode(x, T, ids, threshold_mad=f, gate_r=gate_r, n_pass=n_pass, refr_frames=refr_frames)
        pur[f] = cp.unit_purity(C[f].trains(), C[f].amps(), n, dead_ms=DEAD_MS, rejected=C[f].rejected_by_refractory(),
                                assign_refr_frames=refr_frames)[ui]
    L = C[lo_f]
    acc = L.label == ui
    cand = (L.best == ui) & (L.reason != "below_gate")  # accepted plus those the refractory rule removed
    is_low = L.amp_mad < split
    cls = {"low": acc & is_low, "high": acc & ~is_low}
    cand_cls = {"low": cand & is_low, "high": cand & ~is_low}
    fr = {k: L.frames[m] for k, m in cls.items()}
    cfr = {k: L.frames[m] for k, m in cand_cls.items()}
    hi_train = C[hi_f].train(ui)
    t_off = a / FS
    row = {"experiment": exp_name, "window": spec, "start_h": round(a / FS / 3600, 3), "end_h": round(b / FS / 3600, 3),
           "window_s": round(n / FS, 1), "noise_excised_uV": round(float(L.noise.max() * gain), 2),
           "noise_plain_uV": round(float(L.noise_plain.max() * gain), 2), "excised_fraction": round(L.excised_fraction, 3),
           "split_mad": split}
    for f in floors:
        tag = f"{f:g}"
        row[f"n_events_{tag}"] = int(C[f].frames.size)
        row[f"n_unit_{tag}"] = int(C[f].train(ui).size)
        row[f"rate_unit_hz_{tag}"] = round(C[f].train(ui).size / (n / FS), 3)
        row[f"conflict_ratio_{tag}"] = round(pur[f]["refractory_conflict_ratio"], 3)
        row[f"refractory_lost_{tag}"] = int(C[f].rejected_by_refractory()[ui])
        row[f"ampsplit_ratio_{tag}"] = round(pur[f]["ampsplit_ratio"], 3) if pur[f]["ampsplit_ratio"] is not None else None
        row[f"ampsplit_burst_{tag}"] = round(pur[f]["ampsplit_burst_ratio"], 3) if pur[f]["ampsplit_burst_ratio"] is not None else None
        nb = pur[f]["neighbours"]
        row[f"neighbours_dup_amb_dist_{tag}"] = f"{sum(v == 'duplicate' for v in nb.values())}/{sum(v == 'ambiguous' for v in nb.values())}/{sum(v == 'distinct' for v in nb.values())}"
    row["n_low"] = int(fr["low"].size); row["n_high"] = int(fr["high"].size)
    row["frac_low"] = round(fr["low"].size / max(acc.sum(), 1), 3)
    row["high_matched_by_hi_run"] = round(frac_matched(fr["high"], hi_train), 3)
    row["hi_run_matched_by_high"] = round(frac_matched(hi_train, fr["high"]), 3)
    n_cand = int(cand.sum())
    for k in ("low", "high"):
        lost = int((cand_cls[k] & (L.reason == "refractory")).sum())
        row[f"{k}_candidates"] = int(cand_cls[k].sum())
        row[f"{k}_refractory_lost"] = lost
        row[f"{k}_refractory_lost_frac"] = round(lost / max(cand_cls[k].sum(), 1), 4)
        expected = cand_cls[k].sum() * n_cand * 2 * refr_frames / n
        row[f"{k}_conflict_ratio"] = round(lost / expected, 3) if expected > 0 else None
        isi = np.diff(cfr[k]) / FS * 1000
        row[f"{k}_isi_under_1.5ms_frac"] = round(float((isi < 1.5).mean()), 4) if isi.size else None
        v = cp.ccg_verdict(cfr[k], cfr[k], dead_ms=DEAD_MS)
        row[f"{k}_acg_ratio"] = round(v["ratio"], 3) if np.isfinite(v["ratio"]) else None
        row[f"{k}_r_median"] = round(float(np.median(L.r[cls[k]])), 3) if cls[k].any() else None
    v = cp.ccg_verdict(cfr["low"], cfr["high"], dead_ms=DEAD_MS)
    row["ccg_low_high_ratio"] = round(v["ratio"], 3) if np.isfinite(v["ratio"]) else None
    row["ccg_low_high_verdict"] = v["verdict"]
    cc, ch = ccg_hist(cfr["high"], cfr["low"])
    row["ccg_low_high_2_3ms"] = round(burst_ratio(cc, ch), 3)
    # pulses inside the window
    pp = pulses[(pulses.onset_s >= t_off + 0.1) & (pulses.offset_s + 0.3 <= b / FS)]
    on = pp.onset_s.values - t_off
    dur = float(np.median(pp.offset_s.values - pp.onset_s.values)) if len(pp) else float("nan")
    row["n_pulses"] = int(len(pp)); row["pulse_s"] = round(dur, 3) if len(pp) else None
    row["trains"] = " ".join(str(t) for t in sorted(pp.train.unique())) if len(pp) else ""
    tb = trains[(trains.offset_s > t_off) & (trains.onset_s < b / FS)]
    outside = np.ones(n, bool)
    for t0, t1 in zip(tb.onset_s.values - t_off - 1, tb.offset_s.values - t_off + 1):
        outside[max(int(t0 * FS), 0):min(int(t1 * FS), n)] = False
    s_out = outside.sum() / FS
    ts = {k: fr[k] / FS for k in fr}
    ts["run_hi"] = hi_train / FS
    for k in ts:
        row[f"{k}_rate_outside_trains_hz"] = round(float(outside[(ts[k] * FS).astype(int).clip(0, n - 1)].sum() / s_out), 3) if s_out > 0 else None
        if len(pp):
            r = responses(ts[k], on, dur)
            for kk, vv in r.items():
                row[f"{k}_{kk}"] = round(vv, 3)
            row[f"{k}_in_pulse_over_base"] = round(r["in_pulse_hz"] / r["base_hz"], 3) if r["base_hz"] > 0 else None
    if not make_figs:
        return row
    figdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    tt_ms = (np.arange(cp.N_SAMP) - cp.NBEFORE) / FS * 1000
    label = f"{exp_name} {spec} ({row['start_h']:.2f} to {row['end_h']:.2f} h), unit {unit}"
    fig, ax = plt.subplots(3, 4, figsize=(24, 15))
    # amplitude histogram
    A = ax[0, 0]
    edges = np.arange(3.0, 30.25, 0.25)
    A.hist(L.amp_mad[acc], edges, color="0.6", label=f"floor {lo_f:g}: {int(acc.sum())} spikes")
    A.hist(C[hi_f].amp_mad[C[hi_f].label == ui], edges, histtype="step", color="k", lw=1.0, label=f"floor {hi_f:g}: {hi_train.size} spikes")
    for f, c in ((lo_f, COL["low"]), (hi_f, COL["high"])):
        A.axvline(f, color=c, lw=1.0)
    A.set_yscale("log"); A.set_xlabel("trough amplitude (MAD of excised noise, deepest wire)"); A.legend(fontsize=7)
    A.set_title(f"{label}\nexcised noise {np.round(L.noise.astype(float) * gain, 1).tolist()} uV, plain {np.round(L.noise_plain.astype(float) * gain, 1).tolist()} uV", fontsize=8)
    # overlays
    for j, k in enumerate(("low", "high")):
        A = ax[0, 1 + j]
        pick = rng.choice(fr[k].size, min(80, fr[k].size), replace=False) if fr[k].size else np.zeros(0, int)
        if pick.size:
            sn = cp.snippets(x, fr[k][pick]) * gain
            for s in sn:
                for c in range(4):
                    A.plot(tt_ms + c * 3.2, s[:, c], color=COL[k], lw=0.3, alpha=0.4)
        for c in range(4):
            A.plot(tt_ms + c * 3.2, T[ui, :, c] * gain, color="k", lw=1.2)
        rng_txt = f"{lo_f:g} to {split:g} MAD" if k == "low" else f">= {split:g} MAD"
        A.set_title(f"{k}: {rng_txt}, n={fr[k].size}, {pick.size} overlaid, template in black", fontsize=8)
        A.set_xlabel("ms, wires offset by 3.2 ms"); A.set_ylabel("uV")
    A = ax[0, 3]
    for k in ("low", "high"):
        if fr[k].size:
            m = cp.snippets(x, fr[k][rng.choice(fr[k].size, min(2000, fr[k].size), replace=False)]).mean(axis=0)
            for c in range(4):
                A.plot(tt_ms + c * 3.2, m[:, c] / -m.min(), color=COL[k], lw=1.0, label=k if c == 0 else None)
    for c in range(4):
        A.plot(tt_ms + c * 3.2, T[ui, :, c] / -T[ui].min(), color="k", lw=1.0, ls="--", label="template" if c == 0 else None)
    A.set_title("mean waveform of each class, normalised to its trough", fontsize=8); A.legend(fontsize=7)
    # CCG and ACGs before the refractory rule
    A = ax[1, 0]
    A.bar(cc, ch, width=0.5, color="k")
    A.set_title(f"CCG low relative to high, before the refractory rule: ratio {row['ccg_low_high_ratio']} ({v['verdict']}), 2 to 3 ms {row['ccg_low_high_2_3ms']}", fontsize=8)
    A.set_xlabel("ms"); A.axvspan(-1.5, 1.5, color="y", alpha=0.2)
    for j, k in enumerate(("high", "low")):
        A = ax[1, 1 + j]
        c_, h_ = ccg_hist(cfr[k], cfr[k])
        h_[np.abs(c_) < 0.25] = 0
        A.bar(c_, h_, width=0.5, color=COL[k])
        A.set_title(f"ACG {k} ({cfr[k].size} candidates): ratio {row[f'{k}_acg_ratio']}, lost to rule {row[f'{k}_refractory_lost']} "
                    f"({100 * row[f'{k}_refractory_lost_frac']:.1f} %, {row[f'{k}_conflict_ratio']} of chance)", fontsize=8)
        A.set_xlabel("ms"); A.axvspan(-1.5, 1.5, color="y", alpha=0.2)
    A = ax[1, 3]
    if len(pp):
        for k, c in ((("high", COL["high"]), ("low", COL["low"]), ("run_hi", COL["run_hi"]))):
            tbins, r = psth(ts[k], on, dur)
            A.plot(tbins * 1000, r, color=c, lw=1.0, ls="--" if k == "run_hi" else "-",
                   label=f"{k if k != 'run_hi' else f'floor {hi_f:g} run'}: base {row[f'{k}_base_hz']:.1f}, onset {row[f'{k}_onset_hz']:.1f}, in {row[f'{k}_in_pulse_hz']:.1f}, post {row[f'{k}_post_hz']:.1f} /s")
        A.axvspan(0, dur * 1000, color="c", alpha=0.15)
        A.set_title(f"laser PSTH, {len(pp)} pulses of {dur * 1000:.0f} ms (trains {row['trains']}), 5 ms bins", fontsize=8)
        A.set_xlabel("ms from pulse onset"); A.set_ylabel("spikes/s"); A.legend(fontsize=6)
    else:
        A.set_title("no delivered pulses in this window", fontsize=8)
    # rates over the window
    A = ax[2, 0]
    e = np.arange(0, n / FS + 30, 30)
    for k in ("high", "low"):
        A.plot(e[:-1] / 60, np.histogram(ts[k], e)[0] / 30, color=COL[k], lw=0.8, label=f"{k}: {row[f'{k}_rate_outside_trains_hz']} /s outside trains")
    for t0, t1 in zip(tb.onset_s.values - t_off, tb.offset_s.values - t_off):
        A.axvspan(t0 / 60, t1 / 60, color="c", alpha=0.2)
    A.set_xlabel("min into window"); A.set_ylabel("spikes/s (30 s bins)"); A.legend(fontsize=7)
    A.set_title("firing rate of each class over the window (trains shaded)", fontsize=8)
    for j, f in enumerate(floors):
        neighbour_panel(ax[2, 1 + j], C[f], ui, ids,
                        f"floor {f:g}: neighbour CCGs (stacked, normalised); conflict ratio {row[f'conflict_ratio_{f:g}']}, "
                        f"ampsplit {row[f'ampsplit_ratio_{f:g}']}, dup/amb/dist {row[f'neighbours_dup_amb_dist_{f:g}']}")
    A = ax[2, 3]
    for k in ("high", "low"):
        isi = np.diff(cfr[k]) / FS * 1000
        A.hist(isi[isi < 20], np.arange(0, 20, 0.25), histtype="step", color=COL[k], lw=1.0,
               label=f"{k}: {100 * (row[f'{k}_isi_under_1.5ms_frac'] or 0):.2f} % under 1.5 ms")
    A.axvline(1.5, color="y", lw=0.8); A.set_xlabel("ISI (ms), candidates before the rule"); A.legend(fontsize=7)
    A.set_title("ISI of each class before the refractory rule", fontsize=8)
    fig.tight_layout(); fig.savefig(figdir / f"{exp_name}_{spec}.png", dpi=90); plt.close(fig)
    with PdfPages(figdir / f"{exp_name}_{spec}_purity.pdf") as pdf:
        for f in floors:
            fig, axes = plt.subplots(2, 3, figsize=(20, 10))
            tr_, am_, rs_ = C[f].trains(), C[f].amps(), C[f].rs()
            hl = np.zeros(tr_[ui].size, bool)
            for p0, p1 in zip(on, on + dur):
                hl |= (tr_[ui] / FS > p0 + 0.001) & (tr_[ui] / FS < p1)
            cp.plot_unit_purity(axes, x, gain, ui, tr_, am_, rs_, T, list(range(len(ids))), n_frames=n, highlight=hl,
                                title=f"{label} floor {f:g}:", dead_ms=DEAD_MS, n_rejected=C[f].rejected_by_refractory()[ui],
                                assign_refr_frames=refr_frames)
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    return row


def summary_figure(path: pathlib.Path, floors: tuple) -> None:
    files = sorted(RESULTS.glob("tt6_floor_*.csv"))
    df = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    df = df.sort_values(["experiment", "start_h"]).reset_index(drop=True)
    lab = [f"{e.removeprefix('TTM-003-')}\n{w}\n{h:.1f} h" for e, w, h in zip(df.experiment, df.window, df.start_h)]
    xs = np.arange(len(df))
    lo_f, hi_f = floors
    fig, ax = plt.subplots(2, 3, figsize=(20, 9))
    A = ax[0, 0]
    A.bar(xs, df.n_high, color=COL["high"], label=f"high (>= {hi_f:g} MAD)")
    A.bar(xs, df.n_low, bottom=df.n_high, color=COL["low"], label=f"low ({lo_f:g} to {hi_f:g} MAD): gained by {lo_f:g}")
    for i, (nl, nh, f) in enumerate(zip(df.n_low, df.n_high, df.frac_low)):
        A.text(i, nl + nh, f"{100 * f:.0f} %", ha="center", va="bottom", fontsize=7)
    A.set_ylabel("unit spikes in the window"); A.legend(fontsize=7); A.set_title(f"spikes of the unit at floor {lo_f:g}, by class", fontsize=9)
    A = ax[0, 1]
    A.plot(xs, df.low_refractory_lost_frac * 100, "o-", color=COL["low"], label="low: lost to the 1.5 ms rule (%)")
    A.plot(xs, df.high_refractory_lost_frac * 100, "o-", color=COL["high"], label="high: lost to the rule (%)")
    A.plot(xs, df["low_isi_under_1.5ms_frac"] * 100, "s--", color=COL["low"], label="low: own ISI under 1.5 ms (%)")
    A.plot(xs, df["high_isi_under_1.5ms_frac"] * 100, "s--", color=COL["high"], label="high: own ISI under 1.5 ms (%)")
    A.set_ylabel("%"); A.legend(fontsize=7); A.set_title("refractory violations of each class before the rule", fontsize=9)
    A = ax[0, 2]
    A.plot(xs, df.ccg_low_high_ratio, "o-", color="k", label="CCG low x high, 0.3 to 1.5 ms over 5 to 25 ms")
    A.plot(xs, df.low_acg_ratio, "s--", color=COL["low"], label="ACG low")
    A.plot(xs, df.high_acg_ratio, "s--", color=COL["high"], label="ACG high")
    A.plot(xs, df.ccg_low_high_2_3ms, "^:", color="0.5", label="CCG low x high, 2 to 3 ms")
    A.axhline(0.3, color="g", lw=0.6); A.axhline(0.7, color="r", lw=0.6)
    A.set_ylabel("ratio to flank"); A.legend(fontsize=7); A.set_title("the dip: below 0.3 one refractory period, above 0.7 none", fontsize=9)
    A = ax[1, 0]
    for k, c, m in (("high", COL["high"], "o"), ("low", COL["low"], "o"), ("run_hi", COL["run_hi"], "x")):
        col = f"{k}_in_pulse_over_base"
        if col in df:
            A.plot(xs, df[col], m + "-", color=c, label=f"{k if k != 'run_hi' else f'floor {hi_f:g} run'}")
    A.axhline(1, color="0.5", lw=0.6); A.set_yscale("log")
    A.set_ylabel("in-pulse rate / baseline rate"); A.legend(fontsize=7); A.set_title("laser response of each class (windows with a train)", fontsize=9)
    A = ax[1, 1]
    A.plot(xs, df[f"conflict_ratio_{lo_f:g}"], "o-", color=COL["low"], label=f"conflict ratio, floor {lo_f:g}")
    A.plot(xs, df[f"conflict_ratio_{hi_f:g}"], "o-", color=COL["high"], label=f"conflict ratio, floor {hi_f:g}")
    A.plot(xs, df[f"ampsplit_ratio_{lo_f:g}"], "s--", color=COL["low"], label=f"amplitude-split ratio, floor {lo_f:g}")
    A.plot(xs, df[f"ampsplit_ratio_{hi_f:g}"], "s--", color=COL["high"], label=f"amplitude-split ratio, floor {hi_f:g}")
    A.set_ylabel("ratio"); A.legend(fontsize=7); A.set_title("the capture's purity numbers at each floor", fontsize=9)
    A = ax[1, 2]
    A.plot(xs, df.noise_excised_uV, "o-", color="k", label="excised noise, max wire (uV)")
    A.plot(xs, df.noise_plain_uV, "o--", color="0.5", label="plain noise (uV)")
    A2 = A.twinx()
    A2.plot(xs, df.high_matched_by_hi_run, "s-", color=COL["high"], label=f"high spikes also in the floor-{hi_f:g} run")
    A2.set_ylim(0.8, 1.01); A2.set_ylabel("fraction")
    A.set_ylabel("uV"); A.legend(fontsize=7, loc="upper left"); A2.legend(fontsize=7, loc="lower right")
    A.set_title("noise of the window and the agreement of the two runs on the high class", fontsize=9)
    for A in ax.flat:
        A.set_xticks(xs); A.set_xticklabels(lab, fontsize=6)
    fig.suptitle(f"TT6 cell: detection floor {lo_f:g} against {hi_f:g} excised MAD, every window ({len(df)} windows, {len(files)} experiments)", fontsize=11)
    fig.tight_layout(); fig.savefig(path, dpi=100); plt.close(fig)


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--tetrode", default="TT6")
    p.add_argument("--unit", type=int, default=22, help="template id in the OPTO2 capture templates")
    p.add_argument("--templates-experiment", default="TTM-003-OPTO2")
    p.add_argument("--templates-out", default="capture", help="capture variant under sortings_opto whose templates_<TT>.npz are used")
    p.add_argument("--binary", default=None, help="preprocessed binary (default: sortings_opto/asis/preprocessed if present, else mp_baseline/binary)")
    p.add_argument("--windows", nargs="+", default=["full"])
    p.add_argument("--window-s", type=float, default=1800.0)
    p.add_argument("--floors", nargs=2, type=float, default=[3.5, 4.5])
    p.add_argument("--split", type=float, default=None, help="amplitude (MAD) that splits the classes; default the higher floor")
    p.add_argument("--gate-r", type=float, default=0.6)
    p.add_argument("--n-pass", type=int, default=3)
    p.add_argument("--refr-ms", type=float, default=1.5)
    p.add_argument("--pattern", default="A", help="pulse pattern used for the PSTH in a session without a block table")
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--no-figures", action="store_true")
    a = p.parse_args()
    exp = experiment_from_args(a)
    floors = (min(a.floors), max(a.floors))
    split = a.split if a.split is not None else floors[1]
    binary = a.binary
    if binary is None:
        pp = exp.sorting_dir("asis", subdir="sortings_opto") / "preprocessed"
        binary = str(pp if pp.exists() else exp.sorting_dir("mp_baseline") / "binary")
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    group = int(a.tetrode.removeprefix("TT")) - 1
    wires = [str(c) for c in np.asarray(rec.channel_ids)[cp.tetrode_channel_indices(rec, group)]]
    tdir = find_subject(exp.subject).experiment(a.templates_experiment).sorting_dir(a.templates_out, subdir="sortings_opto")
    T, ids = templates_for(a.tetrode, tdir, wires)
    pulses, trains = delivered_pulses(exp.session, a.pattern)
    every_train = all_trains(exp.session)
    bounds = [window_bounds(w, every_train, nfr, a.window_s) for w in a.windows]
    figdir = exp.figures_dir / "tt6_floor"
    refr_frames = int(round(a.refr_ms * FS / 1000))
    print(f"{a.tetrode} unit {a.unit} from {tdir}; {len(a.windows)} windows on {binary}", flush=True)
    rows = Parallel(n_jobs=min(a.workers, len(a.windows)))(
        delayed(run_window)(binary, group, b0, b1, T, ids, a.unit, floors, split, pulses, trains, a.gate_r, a.n_pass,
                            refr_frames, figdir, exp.experiment, w, not a.no_figures)
        for w, (b0, b1) in zip(a.windows, bounds))
    df = pd.DataFrame(rows)
    path = RESULTS / f"tt6_floor_{exp.experiment}.csv"
    df.to_csv(path, index=False)
    pd.set_option("display.width", 250)
    cols = ["window", "start_h", "noise_excised_uV", f"n_unit_{floors[0]:g}", f"n_unit_{floors[1]:g}", "n_low", "frac_low",
            "low_refractory_lost_frac", "high_refractory_lost_frac", "ccg_low_high_ratio", "low_acg_ratio", "high_acg_ratio",
            "low_in_pulse_over_base", "high_in_pulse_over_base", f"conflict_ratio_{floors[0]:g}", f"conflict_ratio_{floors[1]:g}"]
    print(df[[c for c in cols if c in df]].to_string(index=False))
    print(f"wrote {path}")
    if not a.no_figures:
        summary_figure(figdir / "summary.png", floors)
        print(f"figures in {figdir}")


if __name__ == "__main__":
    main()
