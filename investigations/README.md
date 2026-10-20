---
title: Investigations
updated: 2026-10-19
---

# Investigations

Closed method studies. Each directory is a record of a question asked on the
TTM-001 data in 2026, the scripts and notebooks that answered it, and the
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
| [`zarr-compression/`](zarr-compression/README.md) | How to chunk and compress 30 kHz tetrode data for sorting by tetrode | One tetrode per channel chunk, 30 s time chunks, channels in tetrode order, the real sync clock kept; WavPack bps 2.25 (7.1x) or blosc-zstd (1.73x, lossless) for traces | `tetrode_analyses.spikeinterface.convert_recording`, `ttm-convert` |
| [`spike-sorting/`](spike-sorting/README.md) | MountainSort5 and Kilosort4 on 48 h of tetrode data: reproducibility, lossless vs lossy, block and training duration, cross-chunk unit tracking, matching pursuit | The full 48 h MS5 sort is bit-reproducible once both RNGs are controlled; KS4 clustering does not fit a 48 h tetrode on the 32 GB GPU; chunk+match tracking and the matching-pursuit production sort (`97_production_pipeline.py`) with its open concerns | `tetrode_analyses.sorting`, `tetrode_analyses.tracking`, mountainsort5 0.5.9 (the deterministic PCA solver, merged upstream), the `kilosort` fork pin, the `mp_production_sorting` registry entry |
| [`unitrefine-advisory/`](unitrefine-advisory/) | Is UnitRefine's noise/neural classifier usable on tetrodes? | The stock model is near chance on tetrodes; retrained on 23 geometry-free, duration-robust features it ranks units usefully but its probabilities are uncalibrated | `tetrode_analyses.unitrefine_advisory`; the `unitrefine_*` unit properties on the curation analyzer |
| [`si-frame-slice-memory/`](si-frame-slice-memory/frame_slice_memory_FINDINGS.md) | Why `run_sorter_by_property` worker memory scales with the whole parent recording | `BinaryFolderRecording` eagerly loads the full float64 time vector on every reconstruction; `np.load(..., mmap_mode="r")` fixes it, a genuine crop plus `reset_times()` works around it | SpikeInterface (upstream); the sorting drivers crop before sorting |

`spike-sorting/` keeps the tracking (31-45) and matching-pursuit (51-98)
threads together because they share `_track_eval.py`, `_wobble_eval.py`,
`_assignment_eval.py` and `_scoreboard.py`; three tests under
`packages/tetrode-analyses/tests/` import those helpers. The UnitRefine
scripts (46-50) import only the library and sit on their own.
