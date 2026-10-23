"""Candidate 48 h sorting: the June baseline's seeds with shape-based capture in place of the matcher.

Reads the baseline's preprocessed binary and its carry-forward sorting
``assembled_reseed_c12`` (``mp_baseline.py``), so the windows, the seed units
and their identities are the baseline's and only detection and assignment
differ. Per 1800 s window and tetrode, the units with at least
``MIN_SPIKES_TEMPLATE`` baseline spikes in the window give window-local
templates (``ttm.capture.templates_from_trains``), and ``capture_tetrode``
detects every trough at ``--thresh-mad`` times the spike-excised noise,
assigns it by full-window cosine (gate ``--gate-r``, three peeling passes)
with the amplitude free, and enforces each unit's refractory period. Bad
tetrodes of the manifest and templates whose trough is shallower than
``--min-trough-uV`` are not captured, and a tetrode window with no template
above that floor is not detected at all (so a quiet tetrode does not fill
its multi-unit unit with noise troughs); a tetrode's unassigned events form
one ``is_mua`` unit. Stages, cached under
``<session>/sortings_seed42_pcafix/mp_candidate/``:

  capture   ``assembled_capture`` (baseline unit ids plus one multi-unit
            unit per tetrode), ``capture_TT<n>.npz`` (every event with its
            outcome), ``capture_counts.npz`` (spikes per unit and window),
            ``capture_summary.json``
  merge     the June CCG-guarded merge of the single units (propose 0.90,
            fallback 0.95): ``assembled_candidate`` (merged units renumbered,
            multi-unit units appended), ``candidate_merges.json`` (pairs and
            the old-to-new id map)
  analyzer  ``analyzer_candidate.zarr`` as ``mp_baseline.py`` builds it
  score     the June three-axis scoreboard (script 97's stage_score) of the
            single units (coverage, assignment purity, identity) and of the
            whole sorting (coverage only), against the baseline's peaks
            (``spike_coverage.npz``): ``production_scoreboard.json``

Timing per stage goes to ``candidate_timing.json``.

Usage: python subjects/ttm-003/scripts/mp_candidate.py --subject TTM-003 --experiment TTM-003-48H
           [--stages capture,merge,analyzer] [--max-windows N] [--tag T] [--workers 48] [--n-jobs 16]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import resource
import shutil
import sys
import time

import numpy as np
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered
from ttm.sorting import create_group_sparse_analyzer
from ttm.tracking import _UnionFind

SORTING_DIR = pathlib.Path(__file__).resolve().parents[3] / "investigations" / "spike-sorting"
sys.path.insert(0, str(SORTING_DIR))

FS = 30000.0
WIN_S = 1800.0
BASE_NAME = "assembled_reseed_c12"
MIN_SPIKES_TEMPLATE = 100
PROPOSE_COS, FALLBACK_COS = 0.90, 0.95
STAGES = ["capture", "merge", "analyzer", "score"]


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _dir_bytes(path: pathlib.Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.exists() else 0


def _maxrss_gb() -> float:
    return max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss) / 1e6


def capture_task(binary: str, group: int, a: int, b: int, trains: dict, gain: float, *, thresh_mad: float,
                 gate_r: float, n_pass: int, refr_frames: int, min_trough_uV: float) -> dict:
    """One window of one tetrode: window-local templates from the baseline trains, then the capture.
    ``trains`` are whole-recording frames of this tetrode's units; returns absolute frames."""
    rec = si.load(binary)
    x = cp.load_tetrode_traces(rec, group, a, b)
    local = {u: fr[(fr >= a) & (fr < b)] - a for u, fr in trains.items()}
    local = {u: fr for u, fr in local.items() if fr.size >= MIN_SPIKES_TEMPLATE}
    out = {"group": group, "start": a, "end": b, "units": {}, "skipped_shallow": [], "n_events": 0,
           "pool": np.zeros(0, np.int64), "no_templates": False}
    if local:
        T, ids = cp.templates_from_trains(x, local)
        deep = np.array([T[i].min() * gain <= -min_trough_uV for i in range(len(ids))])
        out["skipped_shallow"] = [int(u) for u, d in zip(ids, deep) if not d]
        T, ids = T[deep], [u for u, d in zip(ids, deep) if d]
    if not local or not ids:
        # No template above the floor in this window: nothing is detected, so a tetrode without
        # units (or with only noise-sized ones) does not fill its multi-unit unit with noise troughs.
        out["no_templates"] = True
        return out
    C = cp.capture_tetrode(x, T, ids, threshold_mad=thresh_mad, gate_r=gate_r, refr_frames=refr_frames, n_pass=n_pass)
    rej = C.rejected_by_refractory()
    for i, u in enumerate(ids):
        m = C.label == i
        out["units"][int(u)] = {"frames": C.frames[m] + a, "amp_mad": C.amp_mad[m], "r": C.r[m],
                                "n_rejected": int(rej[i]), "n_seed": int(local[u].size),
                                "trough_uV": float(T[i].min() * gain)}
    out["n_events"] = int(C.frames.size)
    out["pool"] = C.pool + a
    out["noise_uV"] = (C.noise * gain).tolist()
    out["accepted_by_pass"] = np.bincount(C.pass_index[C.pass_index >= 0], minlength=n_pass).tolist()
    return out


def stage_capture(exp, base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    asm_dir = out / "assembled_capture"
    if asm_dir.exists():
        _log(f"capture: reuse {asm_dir.name}")
        return
    binary = str(base_dir / "binary")
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    gain = float(rec.get_property("gain_to_uV")[0])
    rec_groups = np.asarray(rec.get_property("group")).astype(int)
    base = si.load(str(base_dir / BASE_NAME))
    ugroup = {int(u): int(g) for u, g in zip(base.unit_ids, base.get_property("group"))}
    trains = {int(u): np.sort(base.get_unit_spike_train(u).astype(np.int64)) for u in base.unit_ids}
    bad = {int(t.removeprefix("TT")) - 1 for t in exp.badTetrodes}
    groups = [g for g in sorted(set(rec_groups.tolist())) if g not in bad]
    wf = int(WIN_S * FS)
    starts = list(range(0, nfr, wf))
    if a.max_windows:
        starts = starts[:a.max_windows]
    refr_frames = int(round(a.refr_ms * FS / 1000))
    tasks = [(g, s, min(s + wf, nfr)) for s in starts for g in groups]
    _log(f"capture: {len(starts)} windows x {len(groups)} tetrodes (bad {sorted(bad)}), {len(base.unit_ids)} baseline units, "
         f"{a.workers} workers; {a.thresh_mad} MAD, gate {a.gate_r}, {a.n_pass} passes, min trough {a.min_trough_uV} uV")
    t0 = time.perf_counter()
    results = Parallel(n_jobs=a.workers, batch_size=1)(
        delayed(capture_task)(binary, g, s, e, {u: trains[u] for u in trains if ugroup[u] == g}, gain,
                              thresh_mad=a.thresh_mad, gate_r=a.gate_r, n_pass=a.n_pass, refr_frames=refr_frames,
                              min_trough_uV=a.min_trough_uV)
        for g, s, e in tasks)
    _log(f"capture: {len(tasks)} tasks in {(time.perf_counter() - t0) / 60:.1f} min")

    unit_frames = {u: [] for u in trains}
    unit_amp = {u: [] for u in trains}
    unit_r = {u: [] for u in trains}
    unit_rej = {u: 0 for u in trains}
    counts = np.zeros((len(trains), len(starts)), np.int64)
    uidx = {u: i for i, u in enumerate(sorted(trains))}
    pool = {g: [] for g in groups}
    per_tt = {g: {"frames": [], "label": [], "amp_mad": [], "r": [], "window": []} for g in groups}
    summary = {"windows": len(starts), "tetrodes": [f"TT{g + 1}" for g in groups], "bad_tetrodes": exp.badTetrodes,
               "params": {"thresh_mad": a.thresh_mad, "gate_r": a.gate_r, "n_pass": a.n_pass, "refr_ms": a.refr_ms,
                          "min_trough_uV": a.min_trough_uV, "min_spikes_template": MIN_SPIKES_TEMPLATE},
               "n_events": 0, "n_accepted": 0, "n_pool": 0, "shallow_skips": {}, "windows_without_templates": {}, "per_tetrode": {}}
    for res in results:
        g, wi = res["group"], starts.index(res["start"])
        summary["n_events"] += res["n_events"]
        summary["n_pool"] += int(res["pool"].size)
        pool[g].append(res["pool"])
        tt = per_tt[g]
        tt["frames"].append(res["pool"]); tt["label"].append(np.full(res["pool"].size, -1)); tt["amp_mad"].append(np.full(res["pool"].size, np.nan))
        tt["r"].append(np.full(res["pool"].size, np.nan)); tt["window"].append(np.full(res["pool"].size, wi))
        if res["no_templates"]:
            summary["windows_without_templates"][f"TT{g + 1}"] = summary["windows_without_templates"].get(f"TT{g + 1}", 0) + 1
        for u in res["skipped_shallow"]:
            summary["shallow_skips"][str(u)] = summary["shallow_skips"].get(str(u), 0) + 1
        for u, d in res["units"].items():
            unit_frames[u].append(d["frames"]); unit_amp[u].append(d["amp_mad"]); unit_r[u].append(d["r"])
            unit_rej[u] += d["n_rejected"]
            counts[uidx[u], wi] = d["frames"].size
            summary["n_accepted"] += int(d["frames"].size)
            tt["frames"].append(d["frames"]); tt["label"].append(np.full(d["frames"].size, u)); tt["amp_mad"].append(d["amp_mad"])
            tt["r"].append(d["r"]); tt["window"].append(np.full(d["frames"].size, wi))
    for g in groups:
        tt = per_tt[g]
        fr = np.concatenate(tt["frames"]) if tt["frames"] else np.zeros(0, np.int64)
        order = np.argsort(fr, kind="stable")
        np.savez(out / f"capture_TT{g + 1}.npz", frames=fr[order],
                 label=np.concatenate(tt["label"])[order], amp_mad=np.concatenate(tt["amp_mad"])[order],
                 r=np.concatenate(tt["r"])[order], window=np.concatenate(tt["window"])[order])
        summary["per_tetrode"][f"TT{g + 1}"] = {"n_events": int(fr.size), "n_pool": int(sum(p.size for p in pool[g]))}
    new_trains, new_groups, is_mua, rejected = {}, [], [], {}
    for u in sorted(trains):
        fr = np.sort(np.concatenate(unit_frames[u])) if unit_frames[u] else np.zeros(0, np.int64)
        new_trains[u] = fr
        new_groups.append(ugroup[u]); is_mua.append(False); rejected[str(u)] = unit_rej[u]
    next_id = max(trains) + 1
    mua_ids = {}
    for g in groups:
        fr = np.sort(np.concatenate(pool[g])) if pool[g] else np.zeros(0, np.int64)
        new_trains[next_id] = fr
        new_groups.append(g); is_mua.append(True); mua_ids[f"TT{g + 1}"] = next_id
        next_id += 1
    sorting = si.NumpySorting.from_unit_dict([new_trains], sampling_frequency=FS)
    sorting.set_property("group", np.asarray(new_groups))
    sorting.set_property("is_mua", np.asarray(is_mua))
    sorting.save(folder=str(asm_dir))
    np.savez(out / "capture_counts.npz", unit_ids=np.array(sorted(trains)), counts=counts,
             window_start_frames=np.array(starts),
             base_counts=np.stack([np.histogram(trains[u], bins=np.r_[starts, nfr])[0] for u in sorted(trains)]))
    summary["n_units"] = len(trains)
    summary["mua_unit_ids"] = mua_ids
    summary["refractory_rejected"] = rejected
    summary["units_with_no_capture"] = [int(u) for u in sorted(trains) if new_trains[u].size == 0]
    summary["baseline_spikes"] = int(sum(t.size for t in trains.values()))
    (out / "capture_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    _log(f"capture: {summary['n_events']:,} events, {summary['n_accepted']:,} accepted by {len(trains)} units "
         f"(baseline had {summary['baseline_spikes']:,}), {summary['n_pool']:,} pooled into {len(groups)} multi-unit units; "
         f"{len(summary['units_with_no_capture'])} units never captured; shallow skips {len(summary['shallow_skips'])} units")


def stage_merge(base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    from _assignment_eval import ccg_guarded_merge  # noqa: E402

    merged_dir = out / "assembled_candidate"
    if merged_dir.exists():
        _log(f"merge: reuse {merged_dir.name}")
        return
    rec = si.load(str(base_dir / "binary"))
    cap = si.load(str(out / "assembled_capture"))
    is_mua = np.asarray(cap.get_property("is_mua")).astype(bool)
    groups = np.asarray(cap.get_property("group")).astype(int)
    sua_ids = [int(u) for u, m in zip(cap.unit_ids, is_mua) if not m]
    nonempty = [u for u in sua_ids if cap.get_unit_spike_train(u).size > 0]
    sua = cap.select_units(nonempty)
    _log(f"merge: {len(nonempty)} single units (of {len(sua_ids)}; empty ones dropped) -> CCG-guarded merge "
         f"(propose >= {PROPOSE_COS}, fallback >= {FALLBACK_COS})")
    merged, merges = ccg_guarded_merge(sua, rec, propose_cos=PROPOSE_COS, fallback_cos=FALLBACK_COS,
                                       win_s=WIN_S, n_jobs=a.n_jobs)
    union = _UnionFind()
    for u in nonempty:
        union.add(u)
    for x, y, _, _ in merges:
        union.union(x, y)
    comp = {}
    for u in nonempty:
        comp.setdefault(union.find(u), []).append(u)
    id_map = {}
    for new_id, mem in enumerate(comp.values()):
        for u in mem:
            id_map[int(u)] = int(new_id)
    assert all(merged.get_unit_spike_train(id_map[u]).size >= 1 for u in nonempty[:5])
    trains = {int(u): merged.get_unit_spike_train(u).astype(np.int64) for u in merged.unit_ids}
    new_groups = list(np.asarray(merged.get_property("group")).astype(int))
    flags = [False] * len(trains)
    next_id = len(trains)
    mua_map = {}
    for u, m, g in zip(cap.unit_ids, is_mua, groups):
        if m and cap.get_unit_spike_train(u).size > 0:
            trains[next_id] = cap.get_unit_spike_train(u).astype(np.int64)
            new_groups.append(int(g)); flags.append(True); mua_map[int(u)] = next_id
            next_id += 1
    final = si.NumpySorting.from_unit_dict([trains], sampling_frequency=FS)
    final.set_property("group", np.asarray(new_groups))
    final.set_property("is_mua", np.asarray(flags))
    final.save(folder=str(merged_dir))
    (out / "candidate_merges.json").write_text(json.dumps(
        {"propose_cos": PROPOSE_COS, "fallback_cos": FALLBACK_COS, "n_before": len(nonempty), "n_after": int(merged.get_num_units()),
         "merges": [{"a": int(x), "b": int(y), "verdict": v, "cosine": round(c, 4)} for x, y, v, c in merges],
         "id_map": id_map, "mua_id_map": mua_map, "dropped_empty": [u for u in sua_ids if u not in nonempty]}, indent=1) + "\n")
    _log(f"merge: {len(nonempty)} -> {merged.get_num_units()} units ({len(merges)} pair merges) plus {len(mua_map)} non-empty multi-unit units")


def stage_analyzer(base_dir: pathlib.Path, out: pathlib.Path, a) -> None:
    path = out / "analyzer_candidate.zarr"
    if path.exists():
        _log(f"analyzer: reuse {path.name}")
        return
    rec = si.load(str(base_dir / "binary"))
    sorting = si.load(str(out / "assembled_candidate"))
    _log(f"analyzer: {sorting.get_num_units()} units -> {path}")
    analyzer = create_group_sparse_analyzer(sorting, rec, format="zarr", folder=str(path), overwrite=True,
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
    import importlib.util

    spec = importlib.util.spec_from_file_location("p97", SORTING_DIR / "97_production_pipeline.py")
    p97 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(p97)
    p97.OUT = out
    p97.N_JOBS = a.n_jobs
    rec = si.load(str(base_dir / "binary"))
    z = np.load(base_dir / "spike_coverage.npz", mmap_mode="r")
    peaks = (np.asarray(z["peak_sample"]), np.asarray(z["peak_group"]), np.asarray(z["amp_mad"]))
    starts = p97.window_starts(rec.get_num_frames(), max_windows=a.max_windows)
    full = si.load(str(out / "assembled_candidate"))
    is_mua = np.asarray(full.get_property("is_mua")).astype(bool)
    sua = full.select_units([u for u, m in zip(full.unit_ids, is_mua) if not m])
    _log(f"score: {sua.get_num_units()} single units and the whole sorting over {len(starts)} windows")
    p97.stage_score(rec, [("candidate_sua", sua, True, True), ("candidate_sua+mua", full, False, False)], peaks, starts, tag="")


def main() -> None:
    line_buffered()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--max-windows", type=int, default=None)
    ap.add_argument("--tag", default="", help="suffix of the output directory (smoke runs use _smoke)")
    ap.add_argument("--workers", type=int, default=48, help="capture tasks in parallel (one window of one tetrode each)")
    ap.add_argument("--n-jobs", type=int, default=16, help="SpikeInterface jobs for the merge and the analyzer")
    ap.add_argument("--thresh-mad", type=float, default=3.5)
    ap.add_argument("--gate-r", type=float, default=0.6)
    ap.add_argument("--n-pass", type=int, default=3)
    ap.add_argument("--refr-ms", type=float, default=1.5)
    ap.add_argument("--min-trough-uV", type=float, default=40.0, help="templates shallower than this are not captured")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    unknown = set(stages) - set(STAGES)
    if unknown:
        ap.error(f"unknown stages {sorted(unknown)}; choose from {STAGES}")
    base_dir = exp.sorting_dir("mp_baseline")
    out = exp.sorting_dir(f"mp_candidate{a.tag}")
    out.mkdir(parents=True, exist_ok=True)
    timing_path = out / "candidate_timing.json"
    timing = json.loads(timing_path.read_text()) if timing_path.exists() else []
    for stage in STAGES:
        if stage not in stages:
            continue
        t0 = time.time()
        before = _dir_bytes(out)
        if stage == "capture":
            stage_capture(exp, base_dir, out, a)
        elif stage == "merge":
            stage_merge(base_dir, out, a)
        elif stage == "analyzer":
            stage_analyzer(base_dir, out, a)
        elif stage == "score":
            stage_score(base_dir, out, a)
        row = {"stage": stage, "max_windows": a.max_windows, "workers": a.workers,
               "started": time.strftime("%F %T", time.localtime(t0)), "hours": round((time.time() - t0) / 3600, 3),
               "max_rss_gb": round(_maxrss_gb(), 2), "disk_added_gb": round((_dir_bytes(out) - before) / 1e9, 2)}
        timing.append(row)
        timing_path.write_text(json.dumps(timing, indent=2))
        _log(f"{stage}: {row['hours']:.3f} h, max rss {row['max_rss_gb']:.1f} GB, disk +{row['disk_added_gb']:.2f} GB")


if __name__ == "__main__":
    main()
