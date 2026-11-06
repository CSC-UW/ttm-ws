"""``ttm-sort``: the per-tetrode carry-forward sort of one experiment (``docs/sorting_pipeline.md``).

Stages, each resumable, in ``<session>/sortings_carry/<run>/``:

  binary    bandpass and aligned median reference (``preprocess_for_sorting`` defaults), rounded to
            int16 ADC counts (``quantize_int16``; ``--binary-dtype float32`` keeps float32) and
            written to ``binary/``, its wire timing to ``wire_timing.npz``
  seed      MountainSort5 on window 0 (``seed_sort/``), the starting bank ``bank0.npz``, and the
            whitening ``whitening.npz``
  carry     the windows (``carry/spikes_wNN.npz``, ``bank_wNN.npz``), re-seed sorts in ``reseed_sorts/``
  export    the accepted spikes as a sorting in ``sorting/``
  analyzer  noise levels (``noise_levels.csv``), the multi-unit events (``multi_unit/``, cutoffs in
            ``multi_unit/cutoffs.csv``), the sorting with one multi-unit unit per tetrode
            (``sorting_mua/``), the all-spike PC metrics (``pc_metrics.csv``) and the SortingAnalyzer
            (``analyzer.zarr``); ``analyzer_params.json`` records its parameters (``ttm.analyzer``)
  all       the five in order

``carry_params.json`` records the parameters, inputs and package versions. Bad periods are read from
the experiment's ``bad_periods/bad_periods.csv`` when it exists (``ttm-bad-periods``). The manifest's
``badChannels`` and every channel of its ``badTetrodes`` are removed before the reference and are not
sorted (also from a ``--binary`` that holds them), unless ``--keep-bad-channels``. Laser pulses are
read from ``--pulses`` (CSV with ``onset_s`` and ``offset_s`` in session seconds); without it the rule
for pulses is not applied.

Example: ttm-sort all --subject TTM-003 --experiment TTM-003-48H --pulses PULSES.csv
         ttm-sort analyzer --subject TTM-003 --experiment TTM-003-48H --n-jobs 48
"""

from __future__ import annotations

import argparse
import pathlib
import time

import numpy as np

from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered


def _windows(spec: str | None):
    if spec is None:
        return None
    lo, hi = (int(v) for v in spec.split("-"))
    return range(lo, hi + 1)


def main(argv: list[str] | None = None) -> None:
    line_buffered()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["binary", "seed", "carry", "export", "analyzer", "all"])
    add_experiment_args(ap)
    ap.add_argument("--run", default="main", help="run directory name under <session>/sortings_carry/")
    ap.add_argument("--binary", default=None, help="an existing preprocessed binary to sort instead of binary/")
    ap.add_argument("--pulses", default=None, help="CSV of laser pulses (onset_s, offset_s)")
    ap.add_argument("--bad-periods", default=None, help="bad-period table (default: the experiment's)")
    ap.add_argument("--no-bad-periods", action="store_true", help="sort without a bad-period table")
    ap.add_argument("--bank", default=None,
                    help="starting bank (plain templates, made on the same channels) instead of bank0.npz")
    ap.add_argument("--whitening", default=None, help="whitening npz (key W) instead of whitening.npz")
    ap.add_argument("--reseed-sorts", default=None, help="directory of re-seed sorts instead of reseed_sorts/")
    ap.add_argument("--windows", default=None, help="carry only windows A-B (inclusive)")
    ap.add_argument("--binary-dtype", choices=["int16", "float32"], default="int16",
                    help="sample type of the binary stage's output")
    ap.add_argument("--keep-bad-channels", action="store_true",
                    help="keep the manifest's badChannels and badTetrodes in the reference and the sort")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--analyzer-workers", type=int, default=24,
                    help="processes for the multi-unit windows and noise levels (one window in memory each)")
    ap.add_argument("--keep-pc-projections", action="store_true",
                    help="keep the all-spike PC projections (analyzer_work/all_pcs.npy)")
    a = ap.parse_args(argv)

    import spikeinterface as si

    from ttm import carry as cr
    from ttm.sorting import bad_channel_ids

    exp = experiment_from_args(a)
    store = si.read_zarr(str(exp.store_path()))
    drop = [] if a.keep_bad_channels else bad_channel_ids(store, exp.badChannels, exp.badTetrodes)
    cr.log(f"bad channels left out: {[str(c) for c in drop] or 'none'}")
    run_dir = exp.sorting_dir(a.run, subdir="sortings_carry")
    run_dir.mkdir(parents=True, exist_ok=True)
    params = cr.CarryParams()
    binary = pathlib.Path(a.binary) if a.binary else run_dir / "binary"
    bad_path = None if a.no_bad_periods else (pathlib.Path(a.bad_periods) if a.bad_periods
                                               else exp.bad_periods_dir / "bad_periods.csv")
    if bad_path is not None and not bad_path.exists():
        if a.bad_periods:
            raise FileNotFoundError(bad_path)
        bad_path = None
    whitening = pathlib.Path(a.whitening) if a.whitening else run_dir / "whitening.npz"
    bank_path = pathlib.Path(a.bank) if a.bank else run_dir / "bank0.npz"
    cr.write_run_info(run_dir / "carry_params.json", params, store=exp.store_path(), binary=binary,
                      dropped_channels=[str(c) for c in drop],
                      binary_dtype=a.binary_dtype if not a.binary else "as given",
                      bad_periods=bad_path, pulses=a.pulses, bank=bank_path, whitening=whitening,
                      reseed_sorts=a.reseed_sorts or run_dir / "reseed_sorts")
    stages = ["binary", "seed", "carry", "export", "analyzer"] if a.stage == "all" else [a.stage]

    if "binary" in stages and not (binary / "binary.json").exists():
        from ttm.sorting import estimate_recording_wire_timing, preprocess_for_sorting, quantize_int16

        store.reset_times()
        raw = store.remove_channels(drop) if drop else store
        t0 = time.time()
        timing = estimate_recording_wire_timing(raw)
        np.savez(run_dir / "wire_timing.npz", **timing, channel_ids=np.asarray(raw.channel_ids).astype(str))
        cr.log(f"wire timing in {(time.time() - t0) / 60:.1f} min; windows filled from a neighbour: "
               f"{np.flatnonzero(timing['filled']).tolist()}")
        pp = preprocess_for_sorting(raw, wire_timing=timing)
        if a.binary_dtype == "int16":
            pp = quantize_int16(pp)
        t0 = time.time()
        pp.save(format="binary", folder=str(binary), dtype=a.binary_dtype, n_jobs=64, chunk_duration="1s",
                progress_bar=False, overwrite=True)
        cr.log(f"{a.binary_dtype} binary in {(time.time() - t0) / 60:.1f} min -> {binary}")

    if "seed" in stages:
        rec = cr.open_sorting_recording(binary, drop)
        if not whitening.exists():
            np.savez(whitening, W=cr.whitening_matrix(rec, params))
            cr.log(f"whitening -> {whitening}")
        if not bank_path.exists():
            cr.seed_bank(rec, run_dir / "seed_sort", params).save(bank_path)
            cr.log(f"seed bank -> {bank_path}")

    if "carry" in stages:
        rec = cr.open_sorting_recording(binary, drop)
        W = np.load(whitening)["W"]
        groups = np.asarray(rec.get_property("group")).astype(int)
        bank = cr.Bank.load(bank_path).whitened(W, groups)
        fs = rec.get_sampling_frequency()
        pulses = cr.read_pulses(a.pulses, fs) if a.pulses else None
        bad = cr.read_bad_periods(bad_path, fs) if bad_path is not None else {}
        cr.log(f"carry {len(bank.ids)} units; pulses {0 if pulses is None else pulses[0].size}; "
               f"bad periods {'none' if bad_path is None else bad_path}")
        cr.run_carry(binary, W, bank, run_dir / "carry", params=params, windows=_windows(a.windows), pulses=pulses,
                     bad_periods=bad, reseed_dir=a.reseed_sorts or run_dir / "reseed_sorts", drop_channel_ids=drop,
                     n_jobs=a.n_jobs,
                     workers=a.workers)

    if "export" in stages:
        rec = si.load(binary)
        cr.export_sorting(run_dir / "carry", run_dir / "sorting", rec.get_sampling_frequency())

    if "analyzer" in stages:
        _analyzer_stage(a, run_dir, binary, drop, bad_path)


def _analyzer_stage(a, run_dir: pathlib.Path, binary: pathlib.Path, drop: list, bad_path) -> None:
    import pandas as pd
    import spikeinterface as si

    from ttm import analyzer as an
    from ttm import carry as cr

    params = an.AnalyzerParams()
    folder = run_dir / "analyzer.zarr"
    if (folder / "extensions" / "quality_metrics").exists():
        cr.log(f"analyzer exists: {folder}")
        return
    if folder.exists():
        raise FileExistsError(f"{folder} exists without quality metrics (an interrupted build); remove it to rebuild")
    rec = cr.open_sorting_recording(binary, drop)
    fs = rec.get_sampling_frequency()
    bad = cr.read_bad_periods(bad_path, fs) if bad_path is not None else {}
    sorting_path = run_dir / "sorting"
    an.write_params(run_dir / "analyzer_params.json", params, binary=binary, sorting=sorting_path,
                    bad_periods=bad_path, dropped_channels=[str(c) for c in drop])

    noise_csv = run_dir / "noise_levels.csv"
    if noise_csv.exists():
        table = pd.read_csv(noise_csv)
        if list(table["channel"].astype(str)) != [str(c) for c in rec.channel_ids]:
            raise ValueError(f"{noise_csv} is for other channels")
    else:
        t0 = time.time()
        _, table = an.chunk_noise_levels(binary, bad, params, drop_channel_ids=drop, workers=a.analyzer_workers)
        table.to_csv(noise_csv, index=False)
        cr.log(f"noise levels in {time.time() - t0:.0f} s -> {noise_csv}")
    noise = table["noise_counts"].to_numpy()

    mua_dir = run_dir / "multi_unit"
    an.run_multi_unit(binary, sorting_path, noise, bad, mua_dir, params, drop_channel_ids=drop,
                      workers=a.analyzer_workers)
    events, cutoffs = an.multi_unit_events(mua_dir, params)
    cutoffs.to_csv(mua_dir / "cutoffs.csv", index=False)
    cr.log("multi-unit cutoffs:\n" + cutoffs.round(3).to_string(index=False))
    sorting = an.with_multi_unit(si.load(sorting_path), events)
    sorting.save(folder=run_dir / "sorting_mua", overwrite=True)
    sorting = si.load(run_dir / "sorting_mua")
    periods = an.unit_periods(sorting, rec.get_num_frames(), bad)

    pc_csv = run_dir / "pc_metrics.csv"
    counts = pd.Series(sorting.count_num_spikes_per_unit(outputs="dict"), name="num_spikes")
    if pc_csv.exists():
        pc = pd.read_csv(pc_csv).set_index("unit_id")
        if "num_spikes" not in pc or not pc["num_spikes"].reindex(counts.index).equals(counts.astype(pc["num_spikes"].dtype)):
            raise ValueError(f"{pc_csv} is for other spike counts than {run_dir / 'sorting_mua'}; remove it")
    else:
        t0 = time.time()
        pc = an.all_spike_pc_metrics(sorting, rec, run_dir / "analyzer_work", params, n_jobs=a.n_jobs,
                                     keep_projections=a.keep_pc_projections)
        pc["num_spikes"] = counts.reindex(pc.index).to_numpy()
        pc.to_csv(pc_csv)
        cr.log(f"all-spike PC metrics in {(time.time() - t0) / 60:.1f} min -> {pc_csv}")

    t0 = time.time()
    an.build_analyzer(sorting, rec, folder, table["noise_uv"].to_numpy(), periods, pc, params, n_jobs=a.n_jobs)
    cr.log(f"analyzer in {(time.time() - t0) / 60:.1f} min -> {folder}")


if __name__ == "__main__":
    main()
