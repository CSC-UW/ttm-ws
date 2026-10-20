---
title: Subject TTM-001
updated: 2026-10-19
---

# TTM-001

16-tetrode Open Ephys subject. Manifests in `experiments/`; the pipeline that
consumes them is `docs/pipeline.md`.

Raw acquisitions: `/Volumes/tetrode/TTM1/<session>/` (moved from
`/Volumes/neuropixel_archive/tetrode_data/<session>/` on 2026-10-19 and the old
location emptied on 2026-10-20; the investigations' scripts and notebooks still
name it).

| Experiment | Session | What exists on the server |
|---|---|---|
| `TTM-NOD` | `2026-05-27_09-07-52` (48.3 h, two Open Ephys experiments with a 58 s gap) | blosc-zstd and WavPack stores, 625 Hz and 125 Hz LFPs, `experiment_params.json`, SWA stores and figures, EMG store and figures, the MountainSort5 and Kilosort4 sortings and the matching-pursuit production sort under `sortings_seed42_pcafix/`, viz bundle, curation manifests |
| `TTM-Test` | `2026-05-19_18-09-51` (three 60 s recordings) | `slice_table.csv`, `experiment_params.json` |

Everything under `sortings_seed42_pcafix/` and the method findings are
documented in `investigations/spike-sorting/` and `investigations/zarr-compression/`.
`notebooks/subj01.ipynb` is the early exploration notebook.

Products written before 2026-10-19 were made by the per-script drivers that
the `ttm-*` entry points replaced; the entry points reproduce the same
outputs from the same library functions, and `experiment_params.json` was
rewritten by `ttm-experiment-params` (same values, relative paths).
