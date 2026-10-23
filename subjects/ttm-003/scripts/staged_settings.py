"""Stage 3 settings of the staged pipeline, measured window by window before any full run.

Stage 3 (``ttm.staged_capture``) offers the troughs that the June stages leave to the tetrode's
templates with the amplitude free, gated on the noise-corrected cosine, and keeps them only for
units that answer the light at the most recent delivered train. This script measures what each
setting admits. For each window (a short session whole, or the 1800 s pipeline window of the 48 h
recording that holds a named train) and tetrode it subtracts the June stage-2 sorting
(``mp_baseline/assembled_prod_sua<tag>``) from the preprocessed binary, detects the residual's
troughs at each ``--floors`` value of the spike-excised noise, and runs the offer under each gate
arm: the noise-corrected cosine at each ``--gates`` value and, for reference, the raw cosine 0.6 the
candidate used. For every unit and every delivered train in the window the light test of
``ttm.staged_capture.light_response`` is computed on the unit's June spikes, on the events it was
offered (as if admitted: refractory against its June spikes and among themselves) and on both.

The tagged roles (TT6 cell, TT10 driven cell, TT4 and TT11 candidates) are the June units whose
window template best matches the OPTO2 capture templates (mutual best, as
``tagged_cells_blocks.py``). The sorting-free control on a short session is the capture of the
same window with the OPTO2 templates themselves at 3.5 MAD and raw cosine 0.6 (what
``tagged_cells_xsession.py`` does around each 48 h block, whose table is reused for the 48 h).

Writes ``results/staged_settings_<E>.csv`` (one row per window, tetrode, floor, arm, train and
unit), ``results/staged_settings_<E>_control.csv``, the per-event tables of ``--keep-tetrodes``
under ``<session>/sortings_seed42_pcafix/staged_settings/`` and, with ``--figures``, per tagged
cell and floor a page of what each arm admits and leaves (PSTH around the pulses, waveform
overlays, amplitude distributions) under ``figures/staged_settings/``.

Usage: python subjects/ttm-003/scripts/staged_settings.py --subject TTM-003 --experiment TTM-003-OPTO1 --base-tag _smoke
       python subjects/ttm-003/scripts/staged_settings.py --subject TTM-003 --experiment TTM-003-48H --windows b1 b2 b3 b4 b25 b26 b27
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm import config
from ttm import staged_capture as sc
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from tagged_cells_blocks import TAGGED, reference_templates  # noqa: E402

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
WIN_S = 1800.0
MIN_SPIKES_TEMPLATE = 100
OPTO2_SESSION = "2026-10-16_14-44-22"
DEAD_FRAMES = 9
REFR_FRAMES = 45
BIN = 0.005


def delivered_trains(session: str) -> pd.DataFrame:
    """Delivered trains of ten pulses or more: the block table where one exists, else the hand-logged trains."""
    path = RESULTS / f"stimulation_blocks_{session}.csv"
    if path.exists():
        df = pd.read_csv(path)
        df = df[df.delivered & (df.n_pulses >= 10)].copy()
    else:
        df = pd.read_csv(RESULTS / "stimulation_trains.csv")
        df = df[df.session == session].copy()
    df["train"] = df.ttl_train
    return df.reset_index(drop=True)


def window_spans(spec: list[str], trains: pd.DataFrame, nfr: int) -> list[tuple[str, int, int]]:
    out = []
    wf = int(WIN_S * FS)
    for s in spec:
        if s == "full":
            out.append(("full", 0, nfr))
        elif s.startswith("b"):
            row = trains[trains.train == int(s[1:])].iloc[0]
            a = int(row.onset_s * FS) // wf * wf
            out.append((s, a, min(a + wf, nfr)))
        elif s.startswith("w"):
            a = int(s[1:]) * wf
            out.append((s, a, min(a + wf, nfr)))
        else:
            raise ValueError(s)
    seen, uniq = set(), []
    for name, a, b in out:
        if (a, b) not in seen:
            seen.add((a, b)); uniq.append((name, a, b))
    return uniq


def tagged_roles(T: np.ndarray, ids: list, tt: str, refs: dict, wires_here: list[str]) -> dict:
    """{unit: (role, opto2 unit, cosine)} for the units of one tetrode matched to the OPTO2 templates."""
    if tt not in refs or not ids:
        return {}
    Tref, rids, wires = refs[tt]
    order = [wires.index(w) for w in wires_here]
    cos = np.array([[cp.template_cosine(T[i], Tref[k][:, order]) for k in range(len(rids))] for i in range(len(ids))])
    best = {ids[i]: (rids[int(np.argmax(cos[i]))], float(cos[i].max())) for i in range(len(ids))}
    out = {u: ("", b[0], b[1]) for u, b in best.items()}
    for role, (rtt, rid) in TAGGED.items():
        if rtt != tt or rid not in rids:
            continue
        k = rids.index(rid)
        mutual = [i for i in range(len(ids)) if best[ids[i]][0] == rid]
        pool = mutual or list(range(len(ids)))
        i = max(pool, key=lambda i: cos[i, k])
        out[ids[i]] = (role, rid, float(cos[i, k]))
    return out


def responses_rows(base: dict, sets: dict, pulses_by_train: dict, t0_s: float) -> list[dict]:
    rows = []
    for train, (on, off) in pulses_by_train.items():
        row = dict(base, train=int(train))
        for name, fr in sets.items():
            r = sc.light_response(np.sort(fr) / FS + t0_s, on, off)
            for k in ("pre_hz", "onset_hz", "in_hz", "ratio", "p_count", "p_sign", "n_in", "n_pre"):
                row[f"{name}_{k}"] = r.get(k)
        rows.append(row)
    return rows


def window_task(binary: str, group: int, name: str, a: int, b: int, trains: dict, floors: list, arms: list,
                pulses_by_train: dict, refs: dict, keep_dir: str | None, gain: float) -> dict:
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    t_start = time.perf_counter()
    rec = si.load(binary)
    tt = f"TT{group + 1}"
    wires = [str(c) for c in np.asarray(rec.channel_ids)[cp.tetrode_channel_indices(rec, group)]]
    x = cp.load_tetrode_traces(rec, group, a, b)
    noise, _ = cp.noise_mad_excised(x)
    local = {u: fr[(fr >= a) & (fr < b)] - a for u, fr in trains.items()}
    all_sorted = np.sort(np.concatenate(list(local.values()))) if local else np.zeros(0, np.int64)
    bank = {u: fr for u, fr in local.items() if fr.size >= MIN_SPIKES_TEMPLATE}
    out = {"rows": [], "info": []}
    if not bank:
        return out
    T, ids = cp.templates_from_trains(x, bank)
    roles = tagged_roles(T, ids, tt, refs, wires)
    lab_of = {u: i for i, u in enumerate(ids)}
    frames = np.concatenate([bank[u] for u in ids])
    labels = np.concatenate([np.full(bank[u].size, lab_of[u]) for u in ids])
    order = np.argsort(frames, kind="stable")
    t0 = time.perf_counter()
    residual, amp_fit = sc.subtract_sorted(x, frames[order], labels[order], T)
    t_sub = time.perf_counter() - t0
    for floor in floors:
        t0 = time.perf_counter()
        ev, amp = sc.residual_events(residual, noise, all_sorted, threshold_mad=floor, dead_frames=DEAD_FRAMES)
        t_det = time.perf_counter() - t0
        for arm_name, kind, gate in arms:
            t0 = time.perf_counter()
            off = sc.offer_events(residual, ev, amp, T, noise, gate_rc=gate if kind == "rc" else None,
                                  gate_r=gate if kind == "raw" else None)
            trial = sc.accept_offered(off, np.ones(off.frames.size, bool), bank, ids, refr_frames=REFR_FRAMES)
            t_off = time.perf_counter() - t0
            for i, u in enumerate(ids):
                s3 = off.frames[trial & (off.best == i)]
                left = off.frames[~off.offered & (off.best == i)]
                offered_i = off.frames[off.offered & (off.best == i)]
                lost = int(offered_i.size - s3.size)
                ccg = cp.ccg_verdict(offered_i, bank[u], dead_ms=DEAD_FRAMES / FS * 1000)
                role, ou, oc = roles.get(u, ("", None, None))
                base = {"window": name, "start_s": a / FS, "tetrode": tt, "unit": int(u), "role": role, "opto2_unit": ou,
                        "opto2_cos": None if oc is None else round(oc, 4), "trough_uV": round(float(T[i].min() * gain), 1),
                        "noise_uV_max": round(float(noise.max() * gain), 2), "floor": floor, "arm": arm_name,
                        "gate_kind": kind, "gate": gate, "n_stage2": int(bank[u].size), "n_offered": int((off.offered & (off.best == i)).sum()),
                        "n_stage3": int(s3.size), "n_left": int(left.size),
                        "stage3_amp_median": round(float(np.median(off.amp_mad[trial & (off.best == i)])), 2) if s3.size else None,
                        "n_refractory_lost": lost,
                        "conflict_ratio": round(lost / max(offered_i.size * (s3.size + bank[u].size) * 2 * REFR_FRAMES / (b - a), 1e-9), 3),
                        "ccg_offered_vs_s2": None if not np.isfinite(ccg["ratio"]) else round(ccg["ratio"], 3),
                        "ccg_offered_vs_s2_verdict": ccg["verdict"]}
                sets = {"s2": bank[u], "s3": s3, "all": np.concatenate([bank[u], s3]), "left": left}
                out["rows"] += responses_rows(base, sets, pulses_by_train, a / FS)
            out["info"].append({"window": name, "tetrode": tt, "floor": floor, "arm": arm_name, "n_events": int(ev.size),
                                "n_offered": int(off.offered.sum()), "n_trial": int(trial.sum()), "n_units": len(ids),
                                "s_subtract": round(t_sub, 1), "s_detect": round(t_det, 1), "s_offer": round(t_off, 1),
                                "fit_a_median": round(float(np.nanmedian(amp_fit)), 3)})
            if keep_dir is not None:
                np.savez(pathlib.Path(keep_dir) / f"{name}_{tt}_f{floor:g}_{arm_name}.npz", frames=off.frames + a,
                         amp_mad=off.amp_mad, best=np.array([ids[j] if j >= 0 else -1 for j in off.best]), r=off.r, rc=off.rc,
                         a=off.a, offered=off.offered, trial=trial, pass_index=off.pass_index, templates=T,
                         template_ids=np.array(ids), noise=noise, roles=json.dumps({str(k): v for k, v in roles.items()}))
    out["seconds"] = time.perf_counter() - t_start
    return out


def control_task(binary: str, group: int, name: str, a: int, b: int, refs: dict, pulses_by_train: dict) -> list[dict]:
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    rec = si.load(binary)
    tt = f"TT{group + 1}"
    if tt not in refs:
        return []
    wires = [str(c) for c in np.asarray(rec.channel_ids)[cp.tetrode_channel_indices(rec, group)]]
    Tref, rids, rw = refs[tt]
    order = [rw.index(w) for w in wires]
    T = np.ascontiguousarray(Tref[:, :, order], np.float32)
    x = cp.load_tetrode_traces(rec, group, a, b)
    C = cp.capture_tetrode(x, T, rids, threshold_mad=3.5, gate_r=0.6)
    rows = []
    for role, (rtt, rid) in TAGGED.items():
        if rtt != tt or rid not in rids:
            continue
        fr = C.train(rids.index(rid))
        rows += responses_rows({"window": name, "tetrode": tt, "role": role, "opto2_unit": rid, "n_window": int(fr.size)},
                               {"control": fr}, pulses_by_train, a / FS)
    return rows


def psth(ts: np.ndarray, on: np.ndarray, dur: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(-0.1, dur + 0.2 + BIN / 2, BIN)
    i0 = np.searchsorted(ts, on + edges[0]); i1 = np.searchsorted(ts, on + edges[-1])
    rel = np.concatenate([ts[s:e] - t for s, e, t in zip(i0, i1, on)]) if len(on) else np.zeros(0)
    return edges[:-1] + BIN / 2, np.histogram(rel, edges)[0] / (max(len(on), 1) * BIN)


def cell_pages(exp, keep_dir: pathlib.Path, binary: str, base_tag: str, windows: list, floors: list, arms: list,
               pulses: pd.DataFrame, figdir: pathlib.Path, psth_trains: list[int]) -> None:
    """Per tagged cell, window and floor: one column per arm; PSTH of admitted, left and June spikes around the
    ``psth_trains`` pulses; overlays of admitted in-pulse, admitted outside and left events; amplitudes."""
    rec = si.load(binary)
    gain = float(rec.get_property("gain_to_uV")[0])
    sua = si.load(str(exp.sorting_dir("mp_baseline") / f"assembled_prod_sua{base_tag}"))
    s2 = {int(u): np.sort(sua.get_unit_spike_train(u).astype(np.int64)) for u in sua.unit_ids}
    pp = pulses[pulses.train.isin(psth_trains)]
    tt_ms = (np.arange(cp.N_SAMP) - cp.NBEFORE) / FS * 1000
    rng = np.random.default_rng(0)
    for name, a, b in windows:
        on_all = pp.onset_s.values
        on = on_all[(on_all * FS >= a + 0.2 * FS) & (on_all * FS < b - 0.5 * FS)]
        dur = float(np.median((pp.offset_s - pp.onset_s).values)) if len(pp) else 0.2
        for role, (tt, _) in TAGGED.items():
            if role not in ("TT6_cell", "TT10_driven"):
                continue
            g = int(tt.removeprefix("TT")) - 1
            x = None
            for floor in floors:
                files = [keep_dir / f"{name}_{tt}_f{floor:g}_{arm}.npz" for arm, _, _ in arms]
                if not all(f.exists() for f in files):
                    continue
                z0 = np.load(files[0])
                roles = {int(k): v for k, v in json.loads(str(z0["roles"])).items()}
                unit = next((u for u, v in roles.items() if v[0] == role), None)
                if unit is None:
                    continue
                if x is None:
                    x = cp.load_tetrode_traces(rec, g, a, b)
                ti = list(z0["template_ids"]).index(unit)
                tpl = z0["templates"][ti]
                fig, axes = plt.subplots(5, len(arms), figsize=(3.2 * len(arms), 15), squeeze=False)
                for c, ((arm, kind, gate), f) in enumerate(zip(arms, files)):
                    z = np.load(f)
                    m_acc = z["trial"] & (z["best"] == unit)
                    m_left = ~z["offered"] & (z["best"] == unit)
                    fr_acc, fr_left = z["frames"][m_acc], z["frames"][m_left]
                    ts_acc, ts_left = fr_acc / FS, fr_left / FS
                    ts2 = s2[unit][(s2[unit] >= a) & (s2[unit] < b)] / FS
                    ax = axes[0, c]
                    for ts, col, lab in ((ts2, "k", "June"), (ts_acc, "tab:red", "admitted"), (ts_left, "tab:blue", "left")):
                        tb, rr = psth(np.sort(ts), on, dur)
                        ax.plot(tb * 1000, rr, color=col, lw=0.9, label=lab)
                    ax.axvspan(0, dur * 1000, color="gold", alpha=0.25, lw=0)
                    ax.set_title(f"{arm}: admitted {fr_acc.size}, left {fr_left.size}", fontsize=8)
                    ax.set_xlabel("ms from onset"); ax.set_ylabel("spikes/s") if c == 0 else None
                    if c == 0:
                        ax.legend(fontsize=6)
                    in_p = np.zeros(fr_acc.size, bool)
                    if on.size:
                        k = np.searchsorted(on, ts_acc) - 1
                        ok = k >= 0
                        in_p[ok] = (ts_acc[ok] - on[k[ok]] >= 0.01) & (ts_acc[ok] - on[k[ok]] < dur)
                    for row, (sel, col, lab) in enumerate(((fr_acc[in_p], "tab:orange", "admitted in pulse"),
                                                           (fr_acc[~in_p], "tab:red", "admitted outside pulses"),
                                                           (fr_left, "tab:blue", "left (best = cell, under the gate)")), start=1):
                        ax = axes[row, c]
                        loc = sel - a
                        loc = loc[(loc >= cp.NBEFORE) & (loc + cp.NAFTER <= x.shape[0])]
                        pick = np.sort(rng.choice(loc, min(60, loc.size), replace=False)) if loc.size else loc
                        if pick.size:
                            sn = cp.snippets(x, pick) * gain
                            for k_ in range(pick.size):
                                for w in range(4):
                                    ax.plot(tt_ms + w * 3.2, sn[k_, :, w], color=col, lw=0.3, alpha=0.4)
                            mean = sn.mean(axis=0)
                            scale = tpl.min() * gain / mean.min() if mean.min() < 0 else 1.0
                            for w in range(4):
                                ax.plot(tt_ms + w * 3.2, tpl[:, w] * gain / scale, color="k", lw=1.0)
                        ax.set_title(f"{lab}: {sel.size}", fontsize=7)
                        ax.set_ylim(-180, 90)
                    ax = axes[4, c]
                    bins = np.arange(floor, 30, 0.25)
                    ax.hist(z["amp_mad"][m_acc], bins, color="tab:red", alpha=0.6, label="admitted")
                    ax.hist(z["amp_mad"][m_left], bins, color="tab:blue", alpha=0.5, label="left")
                    ax.set_xlabel("residual trough (MAD)"); ax.set_xlim(floor, 14)
                    if c == 0:
                        ax.legend(fontsize=6)
                fig.suptitle(f"{exp.experiment} {name} {role} (June unit {unit}, {tt}), floor {floor:g} MAD: what each gate admits and "
                             f"leaves; PSTH over {on.size} pulses of trains {psth_trains}; black = template scaled to the mean", fontsize=10)
                fig.tight_layout()
                fig.savefig(figdir / f"{exp.experiment}_{name}_{role}_f{floor:g}.png", dpi=90)
                plt.close(fig)


def outside_trains(frames: np.ndarray, trains: pd.DataFrame, margin_s: float = 1.0) -> np.ndarray:
    ts = frames / FS
    m = np.ones(ts.size, bool)
    for on, off in zip(trains.onset_s.values, trains.offset_s.values):
        m &= ~((ts >= on - margin_s) & (ts < off + margin_s))
    return m


def rule_table(experiments: list, floor: float, arm: str) -> pd.DataFrame:
    """The light test at the chosen floor and gate for every unit and train of the finished runs, with the
    count test recomputed from the counts and the admission of ``mp_staged.py``'s defaults. Writes
    ``results/staged_settings_rule.csv`` and ``figures/staged_settings/rule.png`` (under the last experiment)."""
    from scipy.stats import binomtest
    from ttm.subject import resolve_experiment
    from mp_staged import MAX_P, MIN_IN_HZ, MIN_RATIO

    parts = []
    for name in experiments:
        d = pd.read_csv(RESULTS / f"staged_settings_{name}.csv")
        d = d[(d.floor == floor) & (d.arm == arm)].copy()
        d["role"] = d.role.fillna("")
        d["experiment"] = name
        d["p_count"] = [binomtest(int(i), int(i + q), 0.5, alternative="greater").pvalue if i + q > 0 else 1.0
                        for i, q in zip(d.all_n_in, d.all_n_pre)]
        d["admitted"] = (d.all_ratio >= MIN_RATIO) & (d.p_count < MAX_P) & (d.all_in_hz >= MIN_IN_HZ)
        parts.append(d)
    t = pd.concat(parts)
    cols = ["experiment", "train", "tetrode", "unit", "role", "trough_uV", "noise_uV_max", "s2_pre_hz", "s2_onset_hz",
            "s2_in_hz", "s3_pre_hz", "s3_in_hz", "all_pre_hz", "all_in_hz", "all_ratio", "all_n_in", "all_n_pre", "p_count", "admitted"]
    t = t[cols].sort_values(["experiment", "train", "tetrode", "unit"])
    t.to_csv(RESULTS / "staged_settings_rule.csv", index=False, float_format="%.4g")
    fig, axes = plt.subplots(1, len(experiments), figsize=(6 * len(experiments), 6), squeeze=False)
    for ax, name in zip(axes[0], experiments):
        e = t[t.experiment == name]
        trains = sorted(e.train.unique())
        xpos = {tr: i for i, tr in enumerate(trains)}
        for (tt, u), g in e.groupby(["tetrode", "unit"]):
            role = ",".join(sorted(set(g.role) - {""}))
            x = [xpos[v] for v in g.train]
            y = g.all_ratio.clip(0.05, 60)
            if tt == "TT6" and role:
                ax.plot(x, y, "o-", color="tab:red", lw=1.6, label=f"{tt} u{u} {role}", zorder=5)
            elif role:
                ax.plot(x, y, "s-", lw=1.0, label=f"{tt} u{u} {role}", zorder=4)
            else:
                ax.plot(x, y, "-", color="0.6", lw=0.5, zorder=1)
                ax.scatter(np.array(x)[g.admitted.values], y[g.admitted.values], s=10, color="0.3", zorder=2)
        ax.axhline(MIN_RATIO, color="k", ls="--", lw=0.8)
        ax.set_yscale("log"); ax.set_xticks(range(len(trains))); ax.set_xticklabels([str(v) for v in trains])
        ax.set_xlabel("train"); ax.set_ylabel("in-pulse / pre-pulse rate (June plus trial stage-3 spikes)")
        ax.set_title(f"{name}: floor {floor:g}, {arm}; dots = admitted untagged units", fontsize=9)
        ax.legend(fontsize=6, loc="lower left")
    fig.tight_layout()
    exp = resolve_experiment(subject="TTM-003", experiment=experiments[-1])
    fig.savefig(exp.figures_dir / "staged_settings" / "rule.png", dpi=100); plt.close(fig)
    return t


def summarize(experiments: list, floors: list, arms: list) -> None:
    """Per short session, tagged cell, floor and arm: what the arm admits in and outside the pulses, what it
    leaves, and how its events outside the trains sit against the June unit's own spikes (cross-correlogram
    0.3 to 1.5 ms over 5 to 25 ms, next to the June unit's autocorrelogram ratio). Writes
    ``results/staged_settings_summary.csv`` and ``figures/staged_settings/summary_<E>.png``."""
    from ttm.subject import resolve_experiment

    rows = []
    for name in experiments:
        exp = resolve_experiment(subject="TTM-003", experiment=name)
        d = pd.read_csv(RESULTS / f"staged_settings_{name}.csv")
        ctl = pd.read_csv(RESULTS / f"staged_settings_{name}_control.csv")
        trains = delivered_trains(exp.session)
        tag = "_smoke" if name == "TTM-003-OPTO1" else ""
        sua = si.load(str(exp.sorting_dir("mp_baseline") / f"assembled_prod_sua{tag}"))
        keep_dir = exp.sorting_dir("staged_settings")
        nfr = si.load(str(exp.sorting_dir("mp_baseline") / "binary")).get_num_frames()
        t_out = float(outside_trains(np.arange(0, nfr, int(FS)), trains).sum())
        for role in ("TT6_cell", "TT10_driven", "TT4_candidate", "TT11_candidate"):
            sub = d[d.role == role]
            if not len(sub):
                continue
            unit = int(sub.unit.iloc[0])
            tt = sub.tetrode.iloc[0]
            host = np.sort(sua.get_unit_spike_train(unit).astype(np.int64))
            host_out = host[outside_trains(host, trains)]
            acg = cp.ccg_verdict(host_out, host_out, dead_ms=DEAD_FRAMES / FS * 1000)["ratio"]
            for floor in floors:
                for arm, kind, gate in arms:
                    f = keep_dir / f"full_{tt}_f{floor:g}_{arm}.npz"
                    s = sub[(sub.floor == floor) & (sub.arm == arm)]
                    if not f.exists() or not len(s):
                        continue
                    z = np.load(f)
                    off = z["frames"][z["offered"] & (z["best"] == unit)]
                    o_out = off[outside_trains(off, trains)]
                    acc = z["frames"][z["trial"] & (z["best"] == unit)]
                    ccg = cp.ccg_verdict(host_out, o_out, dead_ms=DEAD_FRAMES / FS * 1000)
                    row = {"experiment": name, "role": role, "unit": unit, "tetrode": tt, "floor": floor, "arm": arm,
                           "gate_kind": kind, "gate": gate, "host_acg_ratio_outside": round(acg, 3),
                           "offered_vs_host_ratio_outside": round(ccg["ratio"], 3) if np.isfinite(ccg["ratio"]) else None,
                           "admitted_outside_hz": round(acc[outside_trains(acc, trains)].size / t_out, 2),
                           "host_outside_hz": round(host_out.size / t_out, 2),
                           "conflict_ratio": float(s.conflict_ratio.iloc[0])}
                    for pat, grp in s.groupby("train_pattern"):
                        lab = f"{pat}{grp.train_power_mW.max():g}"
                        row[f"{lab}_s2_in_hz"] = round(grp.s2_in_hz.mean(), 1)
                        row[f"{lab}_admitted_in_hz"] = round(grp.s3_in_hz.mean(), 1)
                        row[f"{lab}_admitted_pre_hz"] = round(grp.s3_pre_hz.mean(), 1)
                        row[f"{lab}_left_in_hz"] = round(grp.left_in_hz.mean(), 1)
                        row[f"{lab}_left_pre_hz"] = round(grp.left_pre_hz.mean(), 1)
                        c = ctl[(ctl.role == role) & ctl.train.isin(grp.train.unique())]
                        row[f"{lab}_control_in_hz"] = round(c.control_in_hz.mean(), 1) if len(c) else None
                    rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "staged_settings_summary.csv", index=False)
    for name in experiments:
        exp = resolve_experiment(subject="TTM-003", experiment=name)
        e = out[out.experiment == name]
        roles = list(e.role.unique())
        fig, axes = plt.subplots(len(roles), 3, figsize=(17, 3.6 * len(roles)), squeeze=False)
        for r_, role in enumerate(roles):
            for floor, ls in zip(floors, ("-", "--")):
                q = e[(e.role == role) & (e.floor == floor) & (e.gate_kind == "rc")].sort_values("gate")
                raw = e[(e.role == role) & (e.floor == floor) & (e.gate_kind == "raw")]
                for col, lab, c in (("A4.2_admitted_in_hz", "admitted in pulse", "tab:red"), ("A4.2_left_in_hz", "left in pulse", "tab:blue"),
                                    ("A4.2_admitted_pre_hz", "admitted pre-pulse", "tab:orange")):
                    if col in q:
                        axes[r_, 0].plot(q.gate, q[col], ls, color=c, marker=".", label=f"{lab}, floor {floor:g}")
                        if len(raw):
                            axes[r_, 0].plot([0.6], raw[col], "x" if floor == floors[0] else "+", color=c, ms=8)
                if "A4.2_control_in_hz" in q and len(q):
                    axes[r_, 0].axhline(q["A4.2_control_in_hz"].iloc[0], color="k", lw=0.8, ls=":")
                axes[r_, 1].plot(q.gate, q.offered_vs_host_ratio_outside, ls, marker=".", color="tab:purple", label=f"floor {floor:g}")
                axes[r_, 1].axhline(q.host_acg_ratio_outside.iloc[0] if len(q) else np.nan, color="k", lw=0.8, ls=":")
                axes[r_, 2].plot(q.gate, q.admitted_outside_hz, ls, marker=".", color="tab:green", label=f"floor {floor:g}")
                axes[r_, 2].axhline(q.host_outside_hz.iloc[0] if len(q) else np.nan, color="k", lw=0.8, ls=":")
            axes[r_, 0].set_title(f"{role} (June unit {e[e.role == role].unit.iloc[0]}), 207 ms 4.2 mW pulses\n"
                                  f"dotted = sorting-free control; x / + = raw cosine 0.6 at floor 3.5 / 4.5", fontsize=8)
            axes[r_, 0].legend(fontsize=6); axes[r_, 0].set_xlabel("noise-corrected cosine gate"); axes[r_, 0].set_ylabel("spikes/s")
            axes[r_, 1].set_title("outside trains: offered events against the June unit's spikes,\n"
                                  "CCG 0.3-1.5 ms over 5-25 ms (dotted = June unit's own ACG)", fontsize=8)
            axes[r_, 1].legend(fontsize=6); axes[r_, 1].set_xlabel("noise-corrected cosine gate")
            axes[r_, 2].set_title("stage-3 spikes outside trains, per s\n(dotted = June unit's rate there)", fontsize=8)
            axes[r_, 2].legend(fontsize=6); axes[r_, 2].set_xlabel("noise-corrected cosine gate")
        fig.suptitle(f"{name}: stage-3 gate and floor for the tagged cells")
        fig.tight_layout(); fig.savefig(exp.figures_dir / "staged_settings" / f"summary_{name}.png", dpi=100); plt.close(fig)
    print(out.to_string())


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--base-tag", default="", help="suffix of the June outputs (assembled_prod_sua<tag>)")
    p.add_argument("--windows", nargs="+", default=["full"])
    p.add_argument("--floors", nargs="+", type=float, default=[3.5, 4.5])
    p.add_argument("--gates", nargs="+", type=float, default=[0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95])
    p.add_argument("--no-raw-arm", action="store_true", help="omit the raw-cosine 0.6 reference arm")
    p.add_argument("--tetrodes", nargs="+", default=None, help="default: every tetrode not marked bad")
    p.add_argument("--keep-tetrodes", nargs="+", default=["TT6", "TT10", "TT4", "TT11"])
    p.add_argument("--psth-trains", nargs="+", type=int, default=None, help="trains whose pulses the pages average (default: pattern A)")
    p.add_argument("--workers", type=int, default=48)
    p.add_argument("--figures", action="store_true")
    p.add_argument("--figures-only", action="store_true")
    p.add_argument("--out-tag", default="")
    p.add_argument("--summarize", nargs="+", default=None, metavar="EXPERIMENT",
                   help="only tabulate and plot finished runs: the tagged cells of the short sessions, and the light test")
    p.add_argument("--rule-floor", type=float, default=3.5)
    p.add_argument("--rule-arm", default="rc0.8")
    a = p.parse_args()
    arms = ([] if a.no_raw_arm else [("raw0.6", "raw", 0.6)]) + [(f"rc{g:g}", "rc", g) for g in a.gates]
    if a.summarize:
        short = [e for e in a.summarize if e != "TTM-003-48H"]
        if short:
            summarize(short, a.floors, arms)
        rule_table(a.summarize, a.rule_floor, a.rule_arm)
        return
    exp = experiment_from_args(a)
    data_root = config.get().data_root
    base_dir = exp.sorting_dir("mp_baseline")
    binary = str(base_dir / "binary")
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    gain = float(rec.get_property("gain_to_uV")[0])
    trains = delivered_trains(exp.session)
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    pulses = pulses[~pulses.glitch & pulses.train.isin(trains.train)]
    windows = window_spans(a.windows, trains, nfr)
    arms = ([] if a.no_raw_arm else [("raw0.6", "raw", 0.6)]) + [(f"rc{g:g}", "rc", g) for g in a.gates]
    keep_dir = exp.sorting_dir(f"staged_settings{a.out_tag}")
    keep_dir.mkdir(parents=True, exist_ok=True)
    figdir = exp.figures_dir / "staged_settings"
    figdir.mkdir(parents=True, exist_ok=True)
    pattern_col = "ttl_pattern" if "ttl_pattern" in trains else None
    psth_trains = a.psth_trains or (trains[trains[pattern_col] == "A"].train.tolist() if pattern_col else trains.train.tolist())
    if a.figures_only:
        cell_pages(exp, keep_dir, binary, a.base_tag, windows, a.floors, arms, pulses, figdir, psth_trains)
        return
    refs = reference_templates([data_root / OPTO2_SESSION / "sortings_opto" / "capture",
                                data_root / OPTO2_SESSION / "sortings_opto" / "capture_tt10seed"])
    sua = si.load(str(base_dir / f"assembled_prod_sua{a.base_tag}"))
    ugroup = {int(u): int(g) for u, g in zip(sua.unit_ids, sua.get_property("group"))}
    s2 = {int(u): np.sort(sua.get_unit_spike_train(u).astype(np.int64)) for u in sua.unit_ids}
    bad = {int(t.removeprefix("TT")) - 1 for t in exp.badTetrodes}
    groups = ([int(t.removeprefix("TT")) - 1 for t in a.tetrodes] if a.tetrodes else
              [g for g in sorted(set(np.asarray(rec.get_property("group")).astype(int).tolist())) if g not in bad])
    keep = {int(t.removeprefix("TT")) - 1 for t in a.keep_tetrodes}
    tasks, ctasks = [], []
    for name, wa, wb in windows:
        pbt = {}
        for t, ptr in pulses.groupby("train"):
            on, off = ptr.onset_s.values, ptr.offset_s.values
            m = (on * FS >= wa + 0.3 * FS) & (off * FS < wb - 0.3 * FS)
            if m.sum() >= 10:
                pbt[int(t)] = (on[m], off[m])
        for g in groups:
            tr = {u: fr for u, fr in s2.items() if ugroup[u] == g}
            tasks.append(delayed(window_task)(binary, g, name, wa, wb, tr, a.floors, arms, pbt, refs,
                                              str(keep_dir) if g in keep else None, gain))
        if exp.experiment != "TTM-003-48H":
            for g in keep:
                ctasks.append(delayed(control_task)(binary, g, name, wa, wb, refs, pbt))
    print(f"{len(tasks)} tasks ({len(windows)} windows x {len(groups)} tetrodes), floors {a.floors}, arms {[x[0] for x in arms]}; "
          f"{len(ctasks)} control tasks; {a.workers} workers", flush=True)
    t0 = time.perf_counter()
    res = Parallel(n_jobs=a.workers, batch_size=1)(tasks + ctasks)
    print(f"done in {(time.perf_counter() - t0) / 60:.1f} min", flush=True)
    rows = [r for out in res[:len(tasks)] for r in out["rows"]]
    info = [r for out in res[:len(tasks)] for r in out["info"]]
    df = pd.DataFrame(rows)
    meta = trains.set_index("train")
    for col in ("onset_s", "power_mW", "n_pulses"):
        if col in meta:
            df[f"train_{col}"] = df.train.map(meta[col])
    if pattern_col:
        df["train_pattern"] = df.train.map(meta[pattern_col])
    df.round(4).to_csv(RESULTS / f"staged_settings_{exp.experiment}{a.out_tag}.csv", index=False)
    pd.DataFrame(info).to_csv(RESULTS / f"staged_settings_{exp.experiment}{a.out_tag}_timing.csv", index=False)
    if ctasks:
        pd.DataFrame([r for rr in res[len(tasks):] for r in rr]).round(4).to_csv(
            RESULTS / f"staged_settings_{exp.experiment}{a.out_tag}_control.csv", index=False)
    print(f"wrote results/staged_settings_{exp.experiment}{a.out_tag}.csv ({len(df)} rows)", flush=True)
    if a.figures:
        cell_pages(exp, keep_dir, binary, a.base_tag, windows, a.floors, arms, pulses, figdir, psth_trains)


if __name__ == "__main__":
    main()
