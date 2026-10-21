"""Spike sort one TTM-003 stimulation session with MountainSort5 scheme 2.

Preprocessing is the library's (bandpass 300-6000 Hz, global common median
reference) with an optional artifact step that replaces a short window
around every laser onset and offset by a linear interpolation before
filtering (``--blank-ms PRE POST``). The preprocessed recording is
materialized to ``<session>/sortings_opto/<variant>/preprocessed`` and
kept, the per-tetrode sort goes to ``by_group/``, and a SortingAnalyzer
with waveforms, templates, amplitudes, correlograms and quality metrics to
``analyzer``.

Usage:
  python subjects/ttm-003/scripts/opto_sort.py --subject TTM-003 --experiment TTM-003-OPTO1 --variant asis
  python subjects/ttm-003/scripts/opto_sort.py ... --variant blanked --blank-ms 0.5 2.0
"""

from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import time

import numpy as np
import pandas as pd
import spikeinterface as si
import spikeinterface.preprocessing as spre
import spikeinterface.sorters as ss

from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered
from ttm.sorting import preprocess_for_sorting

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
SORTINGS_SUBDIR = "sortings_opto"


def laser_edge_frames(rec, session: str) -> np.ndarray:
    df = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    df = df[~df.glitch]
    t = rec.get_times()
    edges = np.sort(np.r_[df.onset_s.values, df.offset_s.values])
    return np.searchsorted(t, t[0] + edges)


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis", help="output name under sortings_opto/")
    p.add_argument("--blank-ms", nargs=2, type=float, default=None, metavar=("PRE", "POST"),
                   help="linearly interpolate PRE ms before to POST ms after each laser edge")
    p.add_argument("--detect-threshold", type=float, default=5.5)
    p.add_argument("--n-jobs", type=int, default=32)
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()
    exp = experiment_from_args(a)

    out = exp.sorting_dir(a.variant, subdir=SORTINGS_SUBDIR)
    if out.exists() and a.overwrite:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    job_kwargs = dict(n_jobs=a.n_jobs, chunk_duration="10s", progress_bar=False)

    rec = si.read_zarr(str(exp.store_path()))
    if a.blank_ms is not None:
        frames = laser_edge_frames(rec, exp.session)
        rec = spre.remove_artifacts(rec, list_triggers=frames.tolist(), ms_before=a.blank_ms[0],
                                    ms_after=a.blank_ms[1], mode="linear")
        print(f"blanking {len(frames)} laser edges, -{a.blank_ms[0]} / +{a.blank_ms[1]} ms")
    pp_dir = out / "preprocessed"
    if pp_dir.exists():
        pp = si.load(str(pp_dir))
    else:
        t0 = time.perf_counter()
        pp = preprocess_for_sorting(rec, cmr="global").save(format="binary", folder=str(pp_dir), **job_kwargs)
        print(f"materialized preprocessed recording in {time.perf_counter() - t0:.0f}s")

    sorting_dir = out / "by_group"
    if (out / "sorting").exists():
        sorting = si.load(str(out / "sorting"))
        print(f"loaded existing sorting with {len(sorting.unit_ids)} units")
    else:
        t0 = time.perf_counter()
        sorting = ss.run_sorter_by_property(
            "mountainsort5", pp, grouping_property="group", folder=str(sorting_dir),
            engine="joblib", engine_kwargs={"n_jobs": 16}, verbose=False,
            scheme="2", filter=False, whiten=True, whitening_seed=42,
            detect_threshold=a.detect_threshold, detect_sign=-1,
        )
        print(f"sorted {len(sorting.unit_ids)} units in {time.perf_counter() - t0:.0f}s")
        sorting = sorting.save(folder=str(out / "sorting"), overwrite=True)

    t0 = time.perf_counter()
    analyzer = si.create_sorting_analyzer(sorting, pp, format="binary_folder", folder=str(out / "analyzer"),
                                          sparse=True, overwrite=True, **job_kwargs)
    analyzer.compute({
        "random_spikes": {"max_spikes_per_unit": 1000, "seed": 42},
        "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
        "templates": {"operators": ["average", "std"]},
        "noise_levels": {},
        "spike_amplitudes": {},
        "correlograms": {"window_ms": 100.0, "bin_ms": 1.0},
        "isi_histograms": {"window_ms": 100.0, "bin_ms": 1.0},
        "template_metrics": {},
        "quality_metrics": {"metric_names": ["num_spikes", "firing_rate", "presence_ratio", "snr",
                                             "isi_violation", "rp_violation", "amplitude_cutoff",
                                             "amplitude_median", "sliding_rp_violation"]},
    }, **job_kwargs)
    print(f"analyzer in {time.perf_counter() - t0:.0f}s -> {out / 'analyzer'}")

    qm = analyzer.get_extension("quality_metrics").get_data()
    qm["group"] = pd.Series(sorting.get_property("group"), index=sorting.unit_ids)
    qm.to_csv(RESULTS / f"units_{exp.experiment}_{a.variant}.csv")
    summary = {"experiment": exp.experiment, "variant": a.variant, "blank_ms": a.blank_ms,
               "n_units": int(len(sorting.unit_ids)),
               "units_per_tetrode": {int(g): int(n) for g, n in qm.groupby("group").size().items()}}
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
