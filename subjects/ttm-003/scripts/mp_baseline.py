"""Run the June 2026 TTM-001 production matching-pursuit pipeline on a TTM-003 session, unchanged.

The June scripts in ``investigations/spike-sorting/`` carry literal TTM-001
paths. This driver imports them as modules, points their output directory at
``<session>/sortings_seed42_pcafix/mp_baseline/`` and keeps every June
parameter, so the result is the documented baseline for the untagged
population. Stages, each cached and skipped when its output exists:

  materialize  bandpass 300 to 6000 Hz plus global common median reference
               over the whole session, float32 binary (``binary/``)
  carry        seed bank from the first 1800 s window's MountainSort5 sort
               (snr >= 5, >= 100 spikes, dedup cosine 0.95), circus-omp
               carry-forward (amplitudes >= 0.8) with per-window template
               re-estimation and re-seeding every 12 windows; the recipe of
               script 58 -> ``assembled_reseed_c12``, ``reseed_c12.npz``
  peaks        locally exclusive negative peaks at 5.5 MAD, script 64 ->
               ``spike_coverage.npz``
  production   CCG-guarded merge, residual single-unit capture, per-tetrode
               multi-unit bucket and the three-axis scoreboard, script 97 ->
               ``assembled_prod``, ``production_scoreboard.json``
  analyzer     curation analyzer, script 98 -> ``analyzer_prod.zarr``

Wall-clock, peak resident memory and disk per stage go to
``baseline_timing.json``.

Usage: python subjects/ttm-003/scripts/mp_baseline.py --subject TTM-003 --experiment TTM-003-48H
           [--stages materialize,carry,peaks,production,analyzer] [--max-windows N]
           [--n-jobs 16] [--materialize-jobs 96] [--tag ""]
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import resource
import shutil
import sys
import time

import numpy as np
import spikeinterface as si

from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered
from ttm.sorting import create_group_sparse_analyzer
from ttm.tracking import Chunk, materialize_chunk, sort_chunk, to_int_numpy_sorting

SORTING_DIR = pathlib.Path(__file__).resolve().parents[3] / "investigations" / "spike-sorting"
sys.path.insert(0, str(SORTING_DIR))
from _mp_common import build_templates_object, dedup_sorting, windowed_carry_forward_reseed  # noqa: E402

FS = 30000.0
STAGES = ["materialize", "carry", "peaks", "production", "analyzer"]
BASE_NAME = "assembled_reseed_c12"
# June parameters (scripts 58, 61, 64, 97, 98)
WINDOW_S = 1800.0
SEED_WINDOW_S = 1800.0
MIN_SPIKES = 100
SEED_MIN_SNR = 5.0
DEDUP_COSINE = 0.95
RESEED_EVERY = 12
RESEED_ADD_COS = 0.8
AMP = {"amplitudes": [0.8, float("inf")]}
DETECT_THRESH = 5.5


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, SORTING_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _dir_bytes(path: pathlib.Path) -> int:
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) if path.exists() else 0


def _maxrss_gb() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    kids = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return max(own, kids) / 1e6


def stage_materialize(exp, out: pathlib.Path, *, n_jobs: int):
    bin_dir = out / "binary"
    if bin_dir.exists():
        rec = si.load(bin_dir)
        _log(f"materialize: reuse {bin_dir} ({rec.get_num_frames() / FS / 3600:.2f} h)")
        return rec
    store = exp.store_path()
    nfr = si.read_zarr(str(store)).get_num_frames()
    _log(f"materialize: {store.name} whole session, {nfr / FS / 3600:.2f} h, {n_jobs} jobs -> {bin_dir}")
    chunk = Chunk(index=0, start_frame=0, end_frame=int(nfr), fs=FS)
    return materialize_chunk(store, chunk, bin_dir, cmr="global", n_jobs=n_jobs)


def stage_carry(rec, out: pathlib.Path, *, n_jobs: int, max_windows: int | None):
    asm_dir = out / BASE_NAME
    if asm_dir.exists():
        s = si.load(asm_dir)
        _log(f"carry: reuse {BASE_NAME} ({s.get_num_units()} units)")
        return s
    n_win = int(np.ceil(rec.get_num_frames() / FS / WINDOW_S))
    _log(f"carry: {n_win} windows of {WINDOW_S:.0f} s; reseed every {RESEED_EVERY} (add cos < {RESEED_ADD_COS}); "
         f"circus-omp amplitudes >= {AMP['amplitudes'][0]}")
    seed_frames = min(int(SEED_WINDOW_S * FS), rec.get_num_frames())
    win0 = rec.frame_slice(0, seed_frames)
    win0.reset_times()
    shutil.rmtree(out / "seed_sort", ignore_errors=True)
    ref0 = to_int_numpy_sorting(sort_chunk(win0, out / "seed_sort"))
    _, az = build_templates_object(ref0, win0, with_snr=True, n_jobs=n_jobs)
    snr = az.get_extension("quality_metrics").get_data()["snr"].to_numpy()
    nsp = np.array([len(ref0.get_unit_spike_train(u)) for u in ref0.unit_ids])
    conf_ids = np.asarray(ref0.unit_ids)[(snr >= SEED_MIN_SNR) & (nsp >= MIN_SPIKES)]
    seed_sorting = dedup_sorting(ref0.select_units(conf_ids), win0, cosine_min=DEDUP_COSINE)
    init_templates, _ = build_templates_object(seed_sorting, win0, with_snr=False, n_jobs=n_jobs)
    _log(f"carry: seed sort {ref0.get_num_units()} units -> {len(conf_ids)} confident -> "
         f"{seed_sorting.get_num_units()} after dedup {DEDUP_COSINE}")
    shutil.rmtree(out / "reseed_sorts", ignore_errors=True)
    asm, counts, births, births_cos = windowed_carry_forward_reseed(
        rec, init_templates, window_s=WINDOW_S, method="circus-omp", method_kwargs=AMP,
        shape_gate_r=None, wobble_factor=None, n_jobs=n_jobs, min_spikes_reestimate=MIN_SPIKES,
        reseed_min_spikes=MIN_SPIKES, reseed_every_windows=RESEED_EVERY, reseed_add_cos=RESEED_ADD_COS,
        reseed_dir=out / "reseed_sorts", max_windows=max_windows)
    all_ids = np.asarray(sorted(counts))
    cnt_mat = np.stack([counts[u] for u in all_ids])
    present = (cnt_mat >= 20).mean(axis=1)
    _log(f"carry: {seed_sorting.get_num_units()} seed -> {len(all_ids)} units ({len(births)} born by "
         f"re-seeding); median present fraction {np.median(present):.2f}")
    asm.save(folder=asm_dir, overwrite=True)
    bids = list(births)
    np.savez(out / "reseed_c12.npz", all_ids=all_ids, counts=cnt_mat, birth_ids=np.array(bids),
             birth_win=np.array([births[u] for u in bids]), birth_cos=np.array([births_cos[u] for u in bids]),
             seed_ids=np.asarray(seed_sorting.unit_ids))
    return asm


def stage_peaks(rec, out: pathlib.Path, *, n_jobs: int):
    npz = out / "spike_coverage.npz"
    if npz.exists():
        _log(f"peaks: reuse {npz.name}")
        return
    from spikeinterface.core import get_noise_levels
    from spikeinterface.sortingcomponents.peak_detection import detect_peaks

    rec_groups = np.asarray(rec.get_property("group"))
    noise = get_noise_levels(rec, return_in_uV=False)
    _log(f"peaks: locally_exclusive, neg, {DETECT_THRESH} MAD, radius 40 um, {n_jobs} jobs")
    peaks = detect_peaks(rec, method="locally_exclusive",
                         method_kwargs={"peak_sign": "neg", "detect_threshold": DETECT_THRESH, "radius_um": 40.0,
                                        "noise_levels": noise},
                         job_kwargs={"n_jobs": n_jobs, "chunk_duration": "1s", "progress_bar": False})
    peak_s = peaks["sample_index"].astype(np.int64)
    peak_ch = peaks["channel_index"].astype(np.int64)
    amp_mad = np.abs(peaks["amplitude"]) / noise[peak_ch]
    _log(f"peaks: {len(peak_s):,} events ({len(peak_s) / (rec.get_num_frames() / FS):.0f}/s); "
         f"amp MAD median {np.median(amp_mad):.1f}, p90 {np.percentile(amp_mad, 90):.1f}")
    np.savez(npz, peak_sample=peak_s, peak_group=rec_groups[peak_ch].astype(np.int64), amp_mad=amp_mad)


def stage_production(rec, out: pathlib.Path, *, n_jobs: int, max_windows: int | None, tag: str):
    p97 = _load_script("97_production_pipeline")
    p97.OUT = out
    p97.N_JOBS = n_jobs
    peaks = p97.load_peaks()
    starts = p97.window_starts(rec.get_num_frames(), max_windows=max_windows)
    _log(f"production: base {BASE_NAME}, theta {p97.MUA_THETA_DEFAULT}, {len(starts)} windows, tag {tag!r}")
    merged = p97.stage_merge(rec, BASE_NAME, f"assembled_prod_merge{tag}", force=False)
    sua = p97.stage_residual(rec, merged, peaks, starts, f"assembled_prod_sua{tag}", force=False)
    prod = p97.stage_mua(rec, sua, peaks, starts, p97.MUA_THETA_DEFAULT, f"assembled_prod{tag}", force=False)
    variants = [("base", si.load(out / BASE_NAME), True, True), ("prod_sua", sua, True, True),
                ("prod_sua+mua", prod, False, False)]
    p97.stage_score(rec, variants, peaks, starts, tag=tag)
    _log(f"production: assembled_prod{tag} has {prod.get_num_units()} units "
         f"({sua.get_num_units()} single plus per-tetrode multi-unit)")
    return prod


def stage_analyzer(rec, out: pathlib.Path, *, n_jobs: int, tag: str):
    analyzer_path = out / f"analyzer_prod{tag}.zarr"
    if analyzer_path.exists():
        _log(f"analyzer: reuse {analyzer_path.name}")
        return
    sorting = si.load(out / f"assembled_prod{tag}")
    _log(f"analyzer: {sorting.get_num_units()} units -> {analyzer_path}")
    analyzer = create_group_sparse_analyzer(sorting, rec, format="zarr", folder=str(analyzer_path),
                                            overwrite=True, return_in_uV=True, n_jobs=n_jobs)
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
    }, n_jobs=n_jobs, progress_bar=False)
    analyzer.compute("quality_metrics", metric_names=[
        "firing_rate", "presence_ratio", "snr", "isi_violation", "rp_violation", "sliding_rp_violation",
        "d_prime", "mahalanobis", "nearest_neighbor", "silhouette", "amplitude_cutoff", "amplitude_median"])
    _log(f"analyzer: extensions {sorted(analyzer.get_loaded_extension_names())}")


def main() -> None:
    line_buffered()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--max-windows", type=int, default=None, help="smoke test: first N windows only")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--materialize-jobs", type=int, default=96)
    ap.add_argument("--tag", default="", help="suffix for production outputs (smoke runs use _smoke)")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    stages = [s.strip() for s in a.stages.split(",") if s.strip()]
    unknown = set(stages) - set(STAGES)
    if unknown:
        ap.error(f"unknown stages {sorted(unknown)}; choose from {STAGES}")
    out = exp.sorting_dir("mp_baseline")
    out.mkdir(parents=True, exist_ok=True)
    timing_path = out / "baseline_timing.json"
    timing = json.loads(timing_path.read_text()) if timing_path.exists() else []

    rec = None
    for stage in STAGES:
        if stage not in stages:
            continue
        t0 = time.time()
        before = _dir_bytes(out)
        if stage == "materialize":
            rec = stage_materialize(exp, out, n_jobs=a.materialize_jobs)
        else:
            rec = rec if rec is not None else si.load(out / "binary")
            if stage == "carry":
                stage_carry(rec, out, n_jobs=a.n_jobs, max_windows=a.max_windows)
            elif stage == "peaks":
                stage_peaks(rec, out, n_jobs=a.n_jobs)
            elif stage == "production":
                stage_production(rec, out, n_jobs=a.n_jobs, max_windows=a.max_windows, tag=a.tag)
            elif stage == "analyzer":
                stage_analyzer(rec, out, n_jobs=a.n_jobs, tag=a.tag)
        row = {"stage": stage, "tag": a.tag, "max_windows": a.max_windows, "started": time.strftime("%F %T", time.localtime(t0)),
               "hours": round((time.time() - t0) / 3600, 3), "max_rss_gb": round(_maxrss_gb(), 2),
               "disk_added_gb": round((_dir_bytes(out) - before) / 1e9, 2)}
        timing.append(row)
        timing_path.write_text(json.dumps(timing, indent=2))
        _log(f"{stage}: {row['hours']:.3f} h, max rss {row['max_rss_gb']:.1f} GB, disk +{row['disk_added_gb']:.1f} GB")


if __name__ == "__main__":
    main()
