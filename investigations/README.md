---
title: Investigations
updated: 2026-10-27
---

# Investigations

Closed method studies. Each directory is a record of a question asked on the
TTM-001 and TTM-003 data in 2026, the scripts and notebooks that answered it, and the
findings document that holds the answer. The decision each one produced lives
in the library or the workspace pins; the directory is not maintained as
runnable code. Paths inside the scripts and notebooks are literal `tononi-2`
paths as they were run (`/nvme/neuropixels/tetrode_data`,
`/Volumes/neuropixel_archive/tetrode_data`, `/nvme/neuropixels/tmp`; the raw
sessions moved to `/Volumes/tetrode/<SUBJECT>/` on 2026-10-19 and the old
location was emptied the next day), and the
READMEs' `cd gfys_workspace && uv run ...` instructions describe the
environment of the time; run anything from this workspace's environment
instead.

| Directory | Question | Outcome | Where the decision lives |
|---|---|---|---|
| [`acquisition-noise/`](acquisition-noise/) | Source of the ~400 Hz comb and ~9.87 Hz transient in the Open Ephys tetrode recordings | Board LEDs and a clock-locked divider ruled out; the comb stays at ~400 Hz across sample rates, pointing to an external source on the shared ground. Report status: needs verification. | Hardware; nothing in the library |
| [`zarr-compression/`](zarr-compression/README.md) | How to chunk and compress 30 kHz tetrode data for sorting by tetrode | One tetrode per channel chunk, 30 s time chunks, channels in tetrode order, the real sync clock kept; WavPack bps 2.25 (7.1x) or blosc-zstd (1.73x, lossless) for traces | `ttm.spikeinterface.convert_recording`, `ttm-convert` |
| [`spike-sorting/`](spike-sorting/README.md) | MountainSort5 and Kilosort4 on 48 h of tetrode data: reproducibility, lossless vs lossy, block and training duration, cross-chunk unit tracking, matching pursuit | The full 48 h MS5 sort is bit-reproducible once both RNGs are controlled; KS4 clustering does not fit a 48 h tetrode on the 32 GB GPU; chunk+match tracking and the matching-pursuit production sort (`97_production_pipeline.py`) with its open concerns | `ttm.sorting`, `ttm.tracking`, mountainsort5 0.5.9 (the deterministic PCA solver, merged upstream), the `kilosort` fork pin, the `mp_production_sorting` registry entry |
| [`unitrefine-advisory/`](unitrefine-advisory/) | Is UnitRefine's noise/neural classifier usable on tetrodes? | The stock model is near chance on tetrodes; retrained on 23 geometry-free, duration-robust features it ranks units usefully but its probabilities are uncalibrated | `ttm.unitrefine_advisory`; the `unitrefine_*` unit properties on the curation analyzer |
| [`si-frame-slice-memory/`](si-frame-slice-memory/frame_slice_memory_FINDINGS.md) | Why `run_sorter_by_property` worker memory scales with the whole parent recording | `BinaryFolderRecording` eagerly loads the full float64 time vector on every reconstruction; `np.load(..., mmap_mode="r")` fixes it, a genuine crop plus `reset_times()` works around it | SpikeInterface (upstream); the sorting drivers crop before sorting |
| [`wobble-fullsize-ttm001/`](wobble-fullsize-ttm001/README.md) | Track G: does wobble beat circus-omp for the full-size stages of the June pipeline, judged on TTM-001's 48 h over 12 windows across a re-seed? | Wobble as June configured it (threshold 0.35 or 0.25 of the median template energy, raw cosine gate 0.60) loses: lower assignment purity, its bank shrinks window by window, and the gate hollows every same-tetrode correlogram so the CCG-guarded merge joins independent cells. With the gate scored on the residual and re-estimation anchored to full-size spikes it is untested (`docs/reports/2026-10-23_ttm001_wobble_fullsize_span.md`) | Nothing changes in the pipeline; `ttm.purity` (purity pages for any sorting) came out of it |
| [`bad-periods-ttm003/`](bad-periods-ttm003/README.md) | Why masking TTM-003's bad periods changes the sorting on tetrodes that have no bad period, and what TT15's units are | circus-omp matches all tetrodes of a 1 s chunk together, so the 2287 to 2696 s event on seven tetrodes costs the clean ones most of their spikes unless masked; TT15's units follow a signal synchronous with the acquisition clock (`docs/reports/2026-10-26_ttm003_bad_periods.md`) | `ttm.bad_periods`, `ttm-bad-periods`; `clock_lock` and the event coverage in `subjects/ttm-003/scripts/bad_periods_compare.py` |
| [`per-tetrode-matching/`](per-tetrode-matching/README.md) | Does circus-omp run on each tetrode's channels on its own keep clean tetrodes' spikes during multi-tetrode artifacts, and what does it change in clean data? | It keeps them (42,730 TT5 to TT12 spikes in 60 s of the unmasked event against 4,943 joint), but on clean data it adds 25 % spikes, mostly poor fits and events shared with other tetrodes, which the joint matcher's competition across tetrodes kept out; amplitude and failure settings do not restore joint's precision. 14-window carry-forwards pending | `ttm.matching`, `mp_baseline.py --per-tetrode-matching` (off by default) |

`spike-sorting/` keeps the tracking (31-45) and matching-pursuit (51-98)
threads together because they share `_track_eval.py`, `_wobble_eval.py`,
`_assignment_eval.py` and `_scoreboard.py`; three tests under
`packages/ttm/tests/` import those helpers. The UnitRefine
scripts (46-50) import only the library and sit on their own.
