"""Staged sorting: the June stages unchanged, then small-spike capture for units that answer the light.

Stages 1 and 2 are the June baseline as ``mp_baseline.py`` ran it (circus-omp
carry-forward at amplitudes of 0.8 of the template or more with per-window
re-estimation and re-seeding every 12 windows, the CCG-guarded merge, residual
capture of unclaimed events of 10 MAD or more at cosine 0.8): this driver reads
their output ``mp_baseline/assembled_prod_sua<base-tag>`` and the preprocessed
binary and does not re-run them. Stage 3 is ``ttm.staged_capture``: per 1800 s
window and tetrode the stage-2 spikes are fitted and subtracted, the residual's
troughs at ``--floor`` times the spike-excised noise are offered to the best
template of the tetrode (every template competes, amplitude free) when the
noise-corrected cosine reaches ``--gate-rc``, with three peeling passes; an
offered event is kept only for a unit that answered the light at the most
recent delivered train (``--min-ratio``, ``--max-p``, ``--min-in-hz``; before
the first train the first train decides) and that keeps its refractory period
against its stage-2 spikes and its other new spikes. Stage 4 is the June
multi-unit bucket on what remains, after the CCG-guarded merge of the single
units. Stages, cached under ``<session>/sortings_seed42_pcafix/mp_staged<tag>/``:

  offer     ``offer/TT<n>_w<NNN>.npz`` (every residual event with its best
            template, both cosines and whether it was offered), first for the
            windows that hold a delivered train (all tetrodes), then for the
            other windows on tetrodes with a unit admitted there
  admit     ``admission.csv`` (the light test per train and unit),
            ``assembled_stage3`` (stage-2 units with their stage-3 spikes added),
            ``stage3_counts.npz``, ``stage3_summary.json``
  merge     the June CCG-guarded merge of the single units (propose 0.90,
            fallback 0.95): ``assembled_staged_sua``, ``staged_merges.json``
  mua       the June per-tetrode multi-unit bucket (unclaimed 5.5 MAD peaks of
            7 MAD or more whose best template cosine is 0.55 to 0.8, script 97):
            ``assembled_staged``
  analyzer  ``analyzer_staged.zarr`` as ``mp_baseline.py`` builds it
  score     the June three-axis scoreboard: ``production_scoreboard.json``

``--pulses-only`` keeps stage-3 spikes only inside the delivered pulses, so
every unit is its June self outside them; with ``--offer-from`` it reuses
another run's offers and only admits, merges and buckets again.
``--between-floor`` and ``--between-gate`` (and ``--inpulse-*``) keep only
the offered events that also clear a stricter depth and noise-corrected
cosine between (inside) the pulses: a filter on the saved offers, which
agrees with re-running the offer at that setting on 94 to 99 % of events
(``staged_split.py --check``).

Timing per stage goes to ``staged_timing.json``.

Usage: python subjects/ttm-003/scripts/mp_staged.py --subject TTM-003 --experiment TTM-003-48H
           [--stages offer,admit,merge,mua,analyzer,score] [--max-windows N] [--tag T] [--workers 48]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import resource
import sys
import time

import numpy as np
import pandas as pd
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm import staged_capture as sc
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered
from ttm.sorting import create_group_sparse_analyzer
from ttm.tracking import _UnionFind

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
SORTING_DIR = pathlib.Path(__file__).resolve().parents[3] / "investigations" / "spike-sorting"
sys.path.insert(0, str(SORTING_DIR))

FS = 30000.0
WIN_S = 1800.0
MIN_SPIKES_TEMPLATE = 100
PROPOSE_COS, FALLBACK_COS = 0.90, 0.95
STAGES = ["offer", "admit", "merge", "mua", "analyzer", "score"]
DEAD_FRAMES = 9
TRAIN_MARGIN_S = 1.0
# stage-3 settings chosen on the October 16 sessions and 48 h blocks 1 to 4 and 25 to 27
# (staged_settings.py; docs/reports/2026-10-23_ttm003_48h_staged_pipeline.md)
FLOOR = 3.5        # excised MAD
GATE_RC = 0.8      # noise-corrected cosine
MIN_RATIO = 2.0    # light test: in-pulse over pre-pulse rate
MAX_P = 1e-3       # light test: one-sided count test
MIN_IN_HZ = 0.0    # light test: in-pulse rate (not used)


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _dir_bytes(path: pathlib.Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.exists() else 0


def _maxrss_gb() -> float:
    return max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1e6


def _p97(out: pathlib.Path, n_jobs: int):
    spec = importlib.util.spec_from_file_location("p97", SORTING_DIR / "97_production_pipeline.py")
    p97 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p97)
    p97.OUT = out
    p97.N_JOBS = n_jobs
    return p97


def offer_dir_of(exp, out: pathlib.Path, a) -> pathlib.Path:
    """Where the offer files are: this run's ``offer/`` or, with ``--offer-from``, another run's."""
    return exp.sorting_dir(f"mp_staged{a.offer_from}") / "offer" if a.offer_from is not None else out / "offer"


def in_pulses(frames: np.ndarray, pulses: pd.DataFrame) -> np.ndarray:
    """Boolean per frame: inside a delivered pulse (onset to offset)."""
    on, off = pulses.onset_s.values, pulses.offset_s.values
    t = frames / FS
    k = np.searchsorted(on, t, side="right") - 1
    ok = k >= 0
    m = np.zeros(t.size, bool)
    m[ok] = t[ok] < off[k[ok]]
    return m


def delivered_trains(session: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(trains, pulses) of the delivered trains of ten pulses or more, from the block table where one
    exists (trains the laser did not deliver are excluded there) or the hand-logged trains."""
    path = RESULTS / f"stimulation_blocks_{session}.csv"
    if path.exists():
        tr = pd.read_csv(path)
        tr = tr[tr.delivered & (tr.n_pulses >= 10)].copy()
    else:
        tr = pd.read_csv(RESULTS / "stimulation_trains.csv")
        tr = tr[tr.session == session].copy()
    tr = tr.sort_values("onset_s").reset_index(drop=True)
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    pulses = pulses[~pulses.glitch & pulses.train.isin(tr.ttl_train)].reset_index(drop=True)
    return tr, pulses


def offer_task(binary: str, group: int, wi: int, a: int, b: int, trains: dict, gain: float, out_npz: str, *,
               floor: float, gate_rc: float, n_pass: int, refr_frames: int) -> dict:
    """One window of one tetrode: templates from the stage-2 trains, subtraction, residual troughs, offer."""
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    t0 = time.perf_counter()
    rec = si.load(binary)
    x = cp.load_tetrode_traces(rec, group, a, b)
    noise, _ = cp.noise_mad_excised(x)
    local = {u: fr[(fr >= a) & (fr < b)] - a for u, fr in trains.items()}
    all_sorted = np.sort(np.concatenate(list(local.values()))) if local else np.zeros(0, np.int64)
    bank = {u: fr for u, fr in local.items() if fr.size >= MIN_SPIKES_TEMPLATE}
    info = {"group": group, "window": wi, "n_events": 0, "n_offered": 0, "n_units": len(bank), "noise_uV": (noise * gain).tolist()}
    if not bank:
        np.savez(out_npz, frames=np.zeros(0, np.int64), amp_mad=np.zeros(0, np.float32), best=np.zeros(0, np.int32),
                 r=np.zeros(0, np.float32), rc=np.zeros(0, np.float32), offered=np.zeros(0, bool), trial=np.zeros(0, bool),
                 pass_index=np.zeros(0, np.int8), template_ids=np.zeros(0, np.int64), trough_uV=np.zeros(0), noise=noise)
        return info
    T, ids = cp.templates_from_trains(x, bank)
    lab = {u: i for i, u in enumerate(ids)}
    frames = np.concatenate([bank[u] for u in ids])
    labels = np.concatenate([np.full(bank[u].size, lab[u]) for u in ids])
    order = np.argsort(frames, kind="stable")
    residual, fit_a = sc.subtract_sorted(x, frames[order], labels[order], T)
    del x
    ev, amp = sc.residual_events(residual, noise, all_sorted, threshold_mad=floor, dead_frames=DEAD_FRAMES)
    off = sc.offer_events(residual, ev, amp, T, noise, gate_rc=gate_rc, n_pass=n_pass)
    trial = sc.accept_offered(off, np.ones(off.frames.size, bool), bank, ids, refr_frames=refr_frames)
    uid = np.asarray(ids, np.int64)
    np.savez(out_npz, frames=off.frames + a, amp_mad=off.amp_mad.astype(np.float32),
             best=np.where(off.best >= 0, uid[np.clip(off.best, 0, None)], -1).astype(np.int32),
             r=off.r.astype(np.float32), rc=off.rc.astype(np.float32), offered=off.offered, trial=trial,
             pass_index=off.pass_index.astype(np.int8), template_ids=uid, trough_uV=T.min(axis=(1, 2)) * gain, noise=noise)
    info.update(n_events=int(off.frames.size), n_offered=int(off.offered.sum()), n_trial=int(trial.sum()),
                fit_a_median=float(np.nanmedian(fit_a)), seconds=round(time.perf_counter() - t0, 1))
    return info


def stage_offer(exp, base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    if a.offer_from is not None:
        _log(f"offer: reuse the offers of mp_staged{a.offer_from}")
        return
    offer_dir = out / "offer"
    offer_dir.mkdir(exist_ok=True)
    binary = str(base_dir / "binary")
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    gain = float(rec.get_property("gain_to_uV")[0])
    sua = si.load(str(base_dir / f"assembled_prod_sua{a.base_tag}"))
    ugroup = {int(u): int(g) for u, g in zip(sua.unit_ids, sua.get_property("group"))}
    s2 = {int(u): np.sort(sua.get_unit_spike_train(u).astype(np.int64)) for u in sua.unit_ids}
    groups = sorted(set(np.asarray(rec.get_property("group")).astype(int).tolist()) -
                    {int(t.removeprefix("TT")) - 1 for t in exp.badTetrodes})
    wf = int(WIN_S * FS)
    starts = list(range(0, nfr, wf))[:a.max_windows] if a.max_windows else list(range(0, nfr, wf))
    trains, _ = delivered_trains(exp.session)
    train_win = sorted({int(s * FS) // wf for s in trains.onset_s if int(s * FS) // wf < len(starts)})
    refr = int(round(a.refr_ms * FS / 1000))
    kw = dict(floor=a.floor, gate_rc=a.gate_rc, n_pass=a.n_pass, refr_frames=refr)

    def run(pairs, label):
        todo = [(g, wi) for g, wi in pairs if not (offer_dir / f"TT{g + 1}_w{wi:03d}.npz").exists()]
        _log(f"offer ({label}): {len(todo)} of {len(pairs)} window-tetrode tasks to run, {a.workers} workers; "
             f"floor {a.floor} MAD, noise-corrected gate {a.gate_rc}, {a.n_pass} passes")
        t0 = time.perf_counter()
        infos = Parallel(n_jobs=a.workers, batch_size=1)(
            delayed(offer_task)(binary, g, wi, starts[wi], min(starts[wi] + wf, nfr),
                                {u: s2[u] for u in s2 if ugroup[u] == g}, gain, str(offer_dir / f"TT{g + 1}_w{wi:03d}.npz"), **kw)
            for g, wi in todo)
        _log(f"offer ({label}): {len(todo)} tasks in {(time.perf_counter() - t0) / 60:.1f} min; "
             f"{sum(i['n_events'] for i in infos):,} events, {sum(i['n_offered'] for i in infos):,} offered")
        return infos

    infos = run([(g, wi) for wi in train_win for g in groups], "windows with a delivered train")
    adm = admission_table(exp, base_dir, out, a, starts, nfr)
    need = governed_windows(adm, trains, starts, nfr)
    infos += run([(g, wi) for (g, wi) in sorted(need) if wi not in train_win], "other windows with an admitted unit")
    log = out / "offer_tasks.json"
    old = json.loads(log.read_text()) if log.exists() else []
    log.write_text(json.dumps(old + infos, indent=0) + "\n")


def admission_table(exp, base_dir: pathlib.Path, out: pathlib.Path, a, starts: list, nfr: int) -> pd.DataFrame:
    """The light test at every delivered train for every unit with a template in the train's window: the unit's
    stage-2 spikes plus its trial stage-3 events (offered, refractory as if admitted)."""
    trains, pulses = delivered_trains(exp.session)
    sua = si.load(str(base_dir / f"assembled_prod_sua{a.base_tag}"))
    ugroup = {int(u): int(g) for u, g in zip(sua.unit_ids, sua.get_property("group"))}
    s2 = {int(u): np.sort(sua.get_unit_spike_train(u).astype(np.int64)) for u in sua.unit_ids}
    wf = int(WIN_S * FS)
    rows = []
    for tr in trains.itertuples(index=False):
        wi = int(tr.onset_s * FS) // wf
        if wi >= len(starts):
            continue
        pp = pulses[pulses.train == tr.ttl_train]
        on, off = pp.onset_s.values, pp.offset_s.values
        lo, hi = int((on[0] - 2) * FS), int((off[-1] + 2) * FS)
        for g in sorted(set(ugroup.values())):
            f = offer_dir_of(exp, out, a) / f"TT{g + 1}_w{wi:03d}.npz"
            if not f.exists():
                continue
            z = np.load(f)
            m = (z["frames"] >= lo) & (z["frames"] < hi) & z["trial"]
            for u in z["template_ids"]:
                u = int(u)
                new = z["frames"][m & (z["best"] == u)]
                old = s2[u][(s2[u] >= lo) & (s2[u] < hi)]
                r_all = sc.light_response(np.sort(np.concatenate([old, new])) / FS, on, off)
                r_old = sc.light_response(old / FS, on, off)
                ok = sc.admit(r_all, min_ratio=a.min_ratio, max_p=a.max_p, min_in_hz=a.min_in_hz)
                rows.append({"train": int(tr.ttl_train), "onset_s": float(on[0]), "onset_h": round(float(on[0]) / 3600, 4),
                             "window": wi, "tetrode": f"TT{g + 1}", "unit": u, "admitted": ok,
                             "pre_hz": r_all["pre_hz"], "in_hz": r_all["in_hz"], "onset_hz": r_all["onset_hz"],
                             "ratio": r_all["ratio"], "p_count": r_all["p_count"], "pulses_up": r_all["pulses_up"],
                             "pulses_down": r_all["pulses_down"], "s2_pre_hz": r_old["pre_hz"], "s2_in_hz": r_old["in_hz"],
                             "s2_onset_hz": r_old["onset_hz"], "n_stage3_trial": int(new.size)})
    adm = pd.DataFrame(rows)
    adm.round(5).to_csv(out / "admission.csv", index=False)
    n_adm = adm.groupby("train").admitted.sum() if len(adm) else pd.Series(dtype=int)
    _log(f"admit: light test at {adm.train.nunique() if len(adm) else 0} trains; admitted units per train "
         f"{n_adm.astype(int).to_dict()}")
    return adm


def governed_windows(adm: pd.DataFrame, trains: pd.DataFrame, starts: list, nfr: int) -> set:
    """(group, window) pairs in which some unit is admitted under a governing train."""
    wf = int(WIN_S * FS)
    bounds = trains.onset_s.values - TRAIN_MARGIN_S
    order = trains.ttl_train.values
    need = set()
    for wi, s in enumerate(starts):
        e = min(s + wf, nfr)
        ks = set(np.unique(sc.governing_train(np.array([s / FS, (e - 1) / FS]), bounds)).tolist())
        ks |= {k for k, t in enumerate(bounds) if s / FS <= t < e / FS}
        for k in ks:
            sub = adm[(adm.train == order[k]) & adm.admitted]
            for tt in sub.tetrode.unique():
                need.add((int(tt.removeprefix("TT")) - 1, wi))
    return need


def stage_admit(exp, base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    asm = out / "assembled_stage3"
    if asm.exists():
        _log(f"admit: reuse {asm.name}")
        return
    rec = si.load(str(base_dir / "binary"))
    nfr = rec.get_num_frames()
    wf = int(WIN_S * FS)
    starts = list(range(0, nfr, wf))[:a.max_windows] if a.max_windows else list(range(0, nfr, wf))
    adm = admission_table(exp, base_dir, out, a, starts, nfr)
    trains, _ = delivered_trains(exp.session)
    bounds = trains.onset_s.values - TRAIN_MARGIN_S
    order = list(trains.ttl_train.values)
    admitted = {(int(r.train), int(r.unit)) for r in adm[adm.admitted].itertuples(index=False)}
    sua = si.load(str(base_dir / f"assembled_prod_sua{a.base_tag}"))
    ugroup = {int(u): int(g) for u, g in zip(sua.unit_ids, sua.get_property("group"))}
    s2 = {int(u): np.sort(sua.get_unit_spike_train(u).astype(np.int64)) for u in sua.unit_ids}
    all_u = np.array(sorted(s2), np.int64)
    adm_mat = np.zeros((len(order), all_u.size), bool)
    for t, u in admitted:
        adm_mat[order.index(t), np.searchsorted(all_u, u)] = True
    refr = int(round(a.refr_ms * FS / 1000))
    new = {u: [] for u in s2}
    counts = np.zeros((len(s2), len(starts)), np.int64)
    uidx = {u: i for i, u in enumerate(sorted(s2))}
    tally = {"offered": 0, "offered_admitted_unit": 0, "accepted": 0, "refractory_lost": 0}
    _, pulses = delivered_trains(exp.session)
    for f in sorted(offer_dir_of(exp, out, a).glob("TT*_w*.npz")):
        g = int(f.stem.split("_")[0].removeprefix("TT")) - 1
        wi = int(f.stem.split("_w")[1])
        if wi >= len(starts):
            continue
        z = np.load(f)
        if z["frames"].size == 0:
            continue
        ids = [int(u) for u in z["template_ids"]]
        k = sc.governing_train(z["frames"] / FS, bounds)
        best = z["best"]
        col = np.searchsorted(all_u, np.clip(best, 0, None))
        elig = (best >= 0) & adm_mat[k, np.clip(col, 0, all_u.size - 1)]
        ids_arr = np.asarray(ids, np.int64)
        srt = np.argsort(ids_arr)
        pos = srt[np.clip(np.searchsorted(ids_arr[srt], best), 0, ids_arr.size - 1)]
        lab = np.where(best >= 0, pos, -1)
        off = sc.Offer(frames=z["frames"], amp_mad=z["amp_mad"], best=lab, r=z["r"].astype(float), rc=z["rc"].astype(float),
                       a=np.zeros(z["frames"].size), offered=z["offered"], pass_index=z["pass_index"])
        inside = in_pulses(z["frames"], pulses)
        offered = z["offered"].copy()
        for region, fl, gt in ((inside, a.inpulse_floor, a.inpulse_gate), (~inside, a.between_floor, a.between_gate)):
            if fl is not None:
                offered &= ~region | (z["amp_mad"] >= fl)
            if gt is not None:
                offered &= ~region | (z["rc"] >= gt)
        off.offered = offered
        keep = sc.accept_offered(off, elig, {u: s2[u] for u in ids}, ids, refr_frames=refr)
        lost = offered & elig & ~keep
        if a.pulses_only:
            keep &= inside
            lost &= inside
        tally["offered"] += int(offered.sum())
        tally["offered_admitted_unit"] += int((offered & elig).sum())
        tally["accepted"] += int(keep.sum())
        tally["refractory_lost"] += int(lost.sum())
        for u in ids:
            fr = z["frames"][keep & (best == u)]
            if fr.size:
                new[u].append(fr)
                counts[uidx[u], wi] += fr.size
    trains_out = {}
    for u in sorted(s2):
        add = np.concatenate(new[u]) if new[u] else np.zeros(0, np.int64)
        trains_out[u] = np.unique(np.concatenate([s2[u], add]))
    sorting = si.NumpySorting.from_unit_dict([trains_out], sampling_frequency=FS)
    sorting.set_property("group", np.asarray([ugroup[u] for u in sorted(s2)]))
    sorting.save(folder=str(asm))
    np.savez(out / "stage3_counts.npz", unit_ids=np.array(sorted(s2)), counts=counts, window_start_frames=np.array(starts))
    per_unit = {str(u): int(counts[uidx[u]].sum()) for u in sorted(s2) if counts[uidx[u]].sum()}
    summary = {"params": {"floor": a.floor, "gate_rc": a.gate_rc, "n_pass": a.n_pass, "refr_ms": a.refr_ms,
                          "min_ratio": a.min_ratio, "max_p": a.max_p, "min_in_hz": a.min_in_hz,
                          "train_margin_s": TRAIN_MARGIN_S, "min_spikes_template": MIN_SPIKES_TEMPLATE,
                          "pulses_only": a.pulses_only, "offer_from": a.offer_from,
                          "between_floor": a.between_floor, "between_gate": a.between_gate,
                          "inpulse_floor": a.inpulse_floor, "inpulse_gate": a.inpulse_gate},
               "windows": len(starts), "tally": tally, "stage2_spikes": int(sum(t.size for t in s2.values())),
               "stage3_spikes_per_unit": per_unit,
               "admitted_per_train": {str(t): sorted(int(u) for tt, u in admitted if tt == t) for t in order}}
    (out / "stage3_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    _log(f"admit: {tally['accepted']:,} stage-3 spikes kept for {len(per_unit)} units "
         f"({tally['offered']:,} offered, {tally['offered_admitted_unit']:,} to admitted units, {tally['refractory_lost']:,} lost to the refractory rule)")


def stage_merge(base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    from _assignment_eval import ccg_guarded_merge  # noqa: E402

    merged_dir = out / "assembled_staged_sua"
    if merged_dir.exists():
        _log(f"merge: reuse {merged_dir.name}")
        return
    rec = si.load(str(base_dir / "binary"))
    s3 = si.load(str(out / "assembled_stage3"))
    ids = [int(u) for u in s3.unit_ids if s3.get_unit_spike_train(u).size > 0]
    s3 = s3.select_units(ids)
    _log(f"merge: {len(ids)} single units -> CCG-guarded merge (propose >= {PROPOSE_COS}, fallback >= {FALLBACK_COS})")
    merged, merges = ccg_guarded_merge(s3, rec, propose_cos=PROPOSE_COS, fallback_cos=FALLBACK_COS, win_s=WIN_S, n_jobs=a.n_jobs)
    union = _UnionFind()
    for u in ids:
        union.add(u)
    for x, y, _, _ in merges:
        union.union(x, y)
    comp = {}
    for u in ids:
        comp.setdefault(union.find(u), []).append(u)
    id_map = {int(u): int(new) for new, mem in enumerate(comp.values()) for u in mem}
    merged.save(folder=str(merged_dir))
    (out / "staged_merges.json").write_text(json.dumps(
        {"propose_cos": PROPOSE_COS, "fallback_cos": FALLBACK_COS, "n_before": len(ids), "n_after": int(merged.get_num_units()),
         "merges": [{"a": int(x), "b": int(y), "verdict": v, "cosine": round(c, 4)} for x, y, v, c in merges],
         "id_map": id_map}, indent=1) + "\n")
    _log(f"merge: {len(ids)} -> {merged.get_num_units()} units ({len(merges)} pair merges)")


def stage_mua(base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    final = out / "assembled_staged"
    if final.exists():
        _log(f"mua: reuse {final.name}")
        return
    p97 = _p97(out, a.n_jobs)
    rec = si.load(str(base_dir / "binary"))
    z = np.load(base_dir / "spike_coverage.npz", mmap_mode="r")
    peaks = (np.asarray(z["peak_sample"]), np.asarray(z["peak_group"]), np.asarray(z["amp_mad"]))
    starts = p97.window_starts(rec.get_num_frames(), max_windows=a.max_windows)
    sua = si.load(str(out / "assembled_staged_sua"))
    p97.stage_mua(rec, sua, peaks, starts, p97.MUA_THETA_DEFAULT, "assembled_staged", force=False)


def stage_analyzer(base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    path = out / "analyzer_staged.zarr"
    if path.exists():
        _log(f"analyzer: reuse {path.name}")
        return
    rec = si.load(str(base_dir / "binary"))
    sorting = si.load(str(out / "assembled_staged"))
    _log(f"analyzer: {sorting.get_num_units()} units -> {path}")
    analyzer = create_group_sparse_analyzer(sorting, rec, format="zarr", folder=str(path), overwrite=False,
                                            return_in_uV=True, n_jobs=a.n_jobs)
    analyzer.compute({
        "random_spikes": {"method": "uniform", "max_spikes_per_unit": 500},
        "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
        "templates": {"operators": ["average", "std"]},
        "noise_levels": {},
        "correlograms": {},
        "isi_histograms": {},
        "template_similarity": {},
        "principal_components": {"n_components": 5, "mode": "by_channel_local"},
        "spike_amplitudes": {},
    }, n_jobs=a.n_jobs, progress_bar=False)
    analyzer.compute("quality_metrics", metric_names=[
        "firing_rate", "presence_ratio", "snr", "isi_violation", "rp_violation", "sliding_rp_violation",
        "d_prime", "mahalanobis", "nearest_neighbor", "silhouette", "amplitude_cutoff", "amplitude_median"])
    _log(f"analyzer: extensions {sorted(analyzer.get_loaded_extension_names())}")


def stage_score(base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    p97 = _p97(out, a.n_jobs)
    rec = si.load(str(base_dir / "binary"))
    z = np.load(base_dir / "spike_coverage.npz", mmap_mode="r")
    peaks = (np.asarray(z["peak_sample"]), np.asarray(z["peak_group"]), np.asarray(z["amp_mad"]))
    starts = p97.window_starts(rec.get_num_frames(), max_windows=a.max_windows)
    full = si.load(str(out / "assembled_staged"))
    is_mua = np.asarray(full.get_property("is_mua")).astype(bool)
    sua = full.select_units([u for u, m in zip(full.unit_ids, is_mua) if not m])
    _log(f"score: {sua.get_num_units()} single units and the whole sorting over {len(starts)} windows")
    p97.stage_score(rec, [("staged_sua", sua, True, True), ("staged_sua+mua", full, False, False)], peaks, starts, tag="")


def main() -> None:
    line_buffered()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--max-windows", type=int, default=None)
    ap.add_argument("--tag", default="", help="suffix of the output directory (smoke runs use _smoke)")
    ap.add_argument("--base-tag", default="", help="suffix of the June outputs read (assembled_prod_sua<base-tag>)")
    ap.add_argument("--workers", type=int, default=48, help="offer tasks in parallel (one window of one tetrode each)")
    ap.add_argument("--n-jobs", type=int, default=16, help="SpikeInterface jobs for the merge, the bucket and the analyzer")
    ap.add_argument("--floor", type=float, default=FLOOR, help="stage-3 detection floor, excised MAD")
    ap.add_argument("--gate-rc", type=float, default=GATE_RC, help="noise-corrected cosine an offered event must reach")
    ap.add_argument("--n-pass", type=int, default=3)
    ap.add_argument("--refr-ms", type=float, default=1.5)
    ap.add_argument("--min-ratio", type=float, default=MIN_RATIO, help="light test: in-pulse over pre-pulse rate")
    ap.add_argument("--max-p", type=float, default=MAX_P, help="light test: one-sided count test")
    ap.add_argument("--min-in-hz", type=float, default=MIN_IN_HZ, help="light test: in-pulse rate")
    ap.add_argument("--pulses-only", action="store_true",
                    help="keep stage-3 spikes only inside delivered pulses (the units are the June ones elsewhere)")
    ap.add_argument("--between-floor", type=float, default=None,
                    help="keep between-pulse offered events only at this depth (excised MAD) or more; a post hoc filter")
    ap.add_argument("--between-gate", type=float, default=None,
                    help="keep between-pulse offered events only at this noise-corrected cosine or more")
    ap.add_argument("--inpulse-floor", type=float, default=None, help="the same inside the delivered pulses")
    ap.add_argument("--inpulse-gate", type=float, default=None, help="the same inside the delivered pulses")
    ap.add_argument("--offer-from", default=None, metavar="TAG",
                    help="read the offer files of mp_staged<TAG> instead of computing them (same floor and gate)")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    unknown = set(stages) - set(STAGES)
    if unknown:
        ap.error(f"unknown stages {sorted(unknown)}; choose from {STAGES}")
    base_dir = exp.sorting_dir("mp_baseline")
    out = exp.sorting_dir(f"mp_staged{a.tag}")
    out.mkdir(parents=True, exist_ok=True)
    timing_path = out / "staged_timing.json"
    timing = json.loads(timing_path.read_text()) if timing_path.exists() else []
    for stage in STAGES:
        if stage not in stages:
            continue
        t0 = time.time()
        before = _dir_bytes(out)
        {"offer": lambda: stage_offer(exp, base_dir, out, a), "admit": lambda: stage_admit(exp, base_dir, out, a),
         "merge": lambda: stage_merge(base_dir, out, a), "mua": lambda: stage_mua(base_dir, out, a),
         "analyzer": lambda: stage_analyzer(base_dir, out, a), "score": lambda: stage_score(base_dir, out, a)}[stage]()
        row = {"stage": stage, "max_windows": a.max_windows, "workers": a.workers,
               "started": time.strftime("%F %T", time.localtime(t0)), "hours": round((time.time() - t0) / 3600, 3),
               "max_rss_gb": round(_maxrss_gb(), 2), "disk_added_gb": round((_dir_bytes(out) - before) / 1e9, 2)}
        timing.append(row)
        timing_path.write_text(json.dumps(timing, indent=2))
        _log(f"{stage}: {row['hours']:.3f} h, max rss {row['max_rss_gb']:.1f} GB, disk +{row['disk_added_gb']:.2f} GB")


if __name__ == "__main__":
    main()
