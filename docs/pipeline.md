---
title: Per-subject pipeline and command-line tools
updated: 2026-10-19
---

# Per-subject pipeline

Every tool takes the experiment it works on from a manifest and every path
from the per-machine config (`docs/SETUP.md`, section 5). Select the
experiment with `--manifest subjects/<id>/experiments/<NAME>.toml`, or
`--subject ID --experiment NAME`, or `TTM_SUBJECT`/`TTM_EXPERIMENT` in the
environment. There is no default.

| Step | Command | Reads | Writes |
|---|---|---|---|
| 1. Convert | `ttm-convert [--compressor blosc-zstd] [--compressor wavpack --bps 2.25]` | raw acquisition on the `archive` remote | `<session>/<session>.<compressor>.zarr`, `slice_table.csv`, `conversion_results.json` |
| 2. Experiment params | `ttm-experiment-params` | manifest + `slice_table.csv` (or the acquisition's first `sync_messages.txt`) | `<subject>/<experiment>/experiment_params.json` |
| 3. LFP | `ttm-make-lfp [--subsample 125] [--check-plot]` | blosc-zstd store | `<session>.lfp.zarr` (625 Hz), `<session>.lfp.125hz.zarr` (tetrode leads) |
| 4. SWA | `ttm-swa all` (or `inst`, `stft`, `plot`) | 625 Hz LFP, params | `delta.idelta.zarr`, `stft_bandpowers.zarr`, `figures/swa_*.png` |
| 5. EMG | `ttm-emg --full` (default: 30 min window, plot only) | blosc-zstd store, `slice_table.csv`, `delta.idelta.zarr` | `synthetic_emg_methods.zarr`, `figures/emg_*.png` |
| 6. Sort | `investigations/spike-sorting/97_production_pipeline.py` (the production matching-pursuit sort; not yet an entry point) | store | `<session>/sortings_seed42_pcafix/<sorting>/` |
| 7. Curate | `ttm-curation-manifest`, `ttm-fetch-bundle curation`, `ttm-curate` | analyzer + the recording it was sorted from | `curation_bundles/<sorting>.manifest.json` |
| 8. View | `ttm-viz-bundle`, `ttm-fetch-bundle viz`, `ttm-loupe` | EMG, 125 Hz LFP, sorting | `<session>/viz_bundle/` |

`<session>` is `<data_root>/<session>/`, `<subject>/<experiment>` is
`<data_root>/<subject>/<experiment>/`. Steps 1 to 6 run on the server; 7 and
8 are hand-offs to a laptop. The library functions behind each step are in
`packages/tetrode-analyses/README.md`.

## Conversion (step 1)

`get_recording` concatenates every Open Ephys experiment of the session on the
real per-sample sync clock, reorders channels to tetrode order, attaches the
tetrode `ProbeGroup`, and annotates provenance; `convert_recording` writes a
Zarr store chunked one tetrode per chunk. The chunking and compressor choice
is the result of `investigations/zarr-compression/`
(blosc-zstd lossless at 1.73x, WavPack bps 2.25 lossy at 7.1x; sorting was
validated on both). Pass `--compressor` twice to write both stores.

The 48 h TTM-001 conversion took about 2 h per store with 32 jobs. Check
`df -h` on the data volume first: the blosc store is ~352 GB.

## LFP (step 3)

Two anti-aliased decimation stages (30000 -> 3750 -> 625 Hz), float32
microvolts, the session time vector carried through. `--subsample 125` keeps
one lead channel per tetrode and decimates by 5; the viewer draws that store.
`--check-plot` writes trace and PSD overlays to `<experiment>/figures/`.

## SWA and EMG (steps 4 and 5)

`ttm-swa inst` is bandpass 0.5-4 Hz + Hilbert envelope, decimated to 125 Hz,
averaged per tetrode (the `wisc_ecephys_tools` rat recipe without bipolar
referencing). `ttm-swa stft` is the `findlay2025a` 4 s DPSS STFT recipe,
q=2, bands delta/vlad/theta/sigma/gamma, averaged per tetrode. `ttm-swa plot`
draws both delta estimates with the light/dark cycle and deprivation window
from `experiment_params.json`; outliers are replaced at plot time only.

`ttm-emg` computes `ecephys.xrsig.synthetic_emg` on one channel each of
TT1/TT8/TT16, decimated to 1500 Hz, both estimators, each Open Ephys
experiment separately so the inter-experiment gap survives on the session
timeline. Run it without `--full` first (a 30 min window 10 h in) to check
the channels, then `--full` to write the store.

## Curation (step 7)

On the server, `ttm-curate` opens the selected experiment's chunk-tracked
analyzer (`sortings_seed42_pcafix/tracked_48h/analyzer_clustered.zarr` by
default, or `--analyzer PATH`) in spikeinterface-gui. `--mode web` (default)
serves on `localhost:8000` without opening a browser; tunnel with
`ssh -N -L 8000:localhost:8000 <you>@tononi-2`. `--mode desktop` needs a
display. `--style grahams_curation` applies the layout, per-view settings and
unit-list columns shipped in `tetrode_analyses/assets/`.

For a laptop, curation data are too large to stage: the analyzer is ~2.5 GB
and the recording it was sorted from ~352 GB. `ttm-curation-manifest` writes
a JSON listing those two items with the relative path each must land at, so
the analyzer's relative link to its recording keeps resolving after the copy.
Then on the laptop:

```bash
ttm-fetch-bundle curation --subject TTM-001 --experiment TTM-NOD --dest /Volumes/Drive/2026-05-27_09-07-52
ttm-curate --bundle /Volumes/Drive/2026-05-27_09-07-52 --mode desktop --style grahams_curation
```

`--no-traces` at manifest, fetch or launch time skips the recording and
curates without the raw-traces view. Only the blosc-zstd store is ever
offered for traces: it is the exact data the sort was produced from.

The analyzer's persisted unit properties (`tier`, `tier_level`, `n_chunks`,
`track_hours` from `investigations/spike-sorting/45_add_curation_unit_properties.py`;
`unitrefine_neural_prob`, `unitrefine_label` from
`tetrode_analyses.unitrefine_advisory`) are sortable columns in the unit
list. The UnitRefine columns are an advisory ranking: sort by
`unitrefine_neural_prob`, do not threshold at the 0.5 label.

## Visualization (step 8)

`ttm-viz-bundle` assembles `<session>/viz_bundle/`: the EMG store copied
verbatim, the 125 Hz LFP re-exported as a plain xarray zarr (a
SpikeInterface-saved zarr is not directly `xr.open_zarr`-able), the sorting
flattened to `spikes.parquet` in session seconds (one `sample_index_to_time`
pass over the whole spike vector), the loupe `state_definitions.json`, and a
`manifest.json` with per-tetrode separator boundaries. The default sorting
is `blosc-43200s-train3600s`; override with `--sorting`.

```bash
ttm-fetch-bundle viz --subject TTM-001 --experiment TTM-NOD --dest /Volumes/Drive/ttm_nod_viz
ttm-loupe --data-dir /Volumes/Drive/ttm_nod_viz
```

`ttm-loupe` needs only the `viz` extra: no SpikeInterface, no source
recording, no server mounts. Panes: both EMG estimators, the 16 LFP leads
coloured by tetrode with the Open Ephys "Classic" palette, and a single-pane
raster with horizontal separators between tetrodes. Interval scoring uses
the bundled state definitions.

## Fetching

`ttm-fetch-bundle` uses rclone and the `server` remote of the per-machine
config, an rclone `remote:path` pointing at the server's data root (for
example `tononi-2:/nvme/neuropixels/tetrode_data`). `--no-check-dest` is the
default for a first one-way copy; `--check` checksums for resumable re-runs.

## Checking a re-run

`scripts/snapshot_derived.py snapshot` summarises an experiment's derived
products into a JSON file and `compare` diffs two summaries. Baselines for
TTM-001 and what they showed are in `docs/baselines/README.md`.
