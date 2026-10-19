---
title: Tetrode preprocessing + sorting study (session 2026-05-27_09-07-52)
updated: 2026-06-10
---

# Tetrode preprocessing + sorting

End-to-end study of one 48 h, 16-tetrode Open Ephys session
(`2026-05-27_09-07-52`, TTM-001 / TTM-NOD): convert + compress to a
sort-ready Zarr store, derive LFPs, spike sort by tetrode group with
MountainSort5, and characterize what perturbs the sort.

> Renamed 2026-06-10 from `chunking_and_compression` and split into the
> subfolders below. The original folder had grown from a compression benchmark
> into four distinct threads; scripts are standalone drivers (they import from
> the `tetrode_analyses` package, not each other), so the split moved files only.

## Subfolders

| folder | thread | canonical doc |
|---|---|---|
| [`compression/`](compression/README.md) | Chunking + compressor benchmarks; the chunking/compression decision (channel_chunk=4, WavPack bps=2.25 vs blosc-zstd) now implemented by `ttm-convert` | `compression/README.md` |
| [`sorting/`](sorting/README.md) | MountainSort5 sort by tetrode group; lossless-vs-lossy agreement, determinism, block/training-duration sweeps, the ms5 int32 ceiling, parameter comparison | `sorting/SORTING_COMPARISON_FINDINGS.md`, `sorting/ms5_parameter_comparison.md` |
| `lfp/`, `visualization/`, `curation/` | Replaced 2026-10-19 by the `ttm-make-lfp`, `ttm-viz-bundle`/`ttm-fetch-bundle`/`ttm-loupe` and `ttm-curation-manifest`/`ttm-curate` entry points | `../../docs/pipeline.md` |
| [`si_frame_slice_memory/`](si_frame_slice_memory/) | Upstream SpikeInterface bug: `frame_slice` / `BinaryFolderRecording` worker memory scales with the full parent recording, not the slice (root-caused while measuring `training_duration` footprint) | `si_frame_slice_memory/frame_slice_memory_FINDINGS.md` |

## Pipeline order

1. `compression/` — `01`–`07`, `17`: the benchmarks. The conversion itself
   (formerly `06_convert.py` / `08_convert_session.py`) is `ttm-convert`.
2. LFP — `ttm-make-lfp` (formerly `lfp/09`, `lfp/13`).
3. `sorting/` — `10`–`24`: sort each store, compare lossless vs lossy, then sweep
   sorter parameters (determinism, block/training duration) and curate.
4. Visualization — `ttm-viz-bundle` → `ttm-fetch-bundle viz` → `ttm-loupe`.
5. Curation — `ttm-curate` on the server, or `ttm-curation-manifest` →
   `ttm-fetch-bundle curation` → `ttm-curate --bundle` locally.

Numeric script prefixes (`01`–`24`) are the original chronological order across
all threads; they are preserved within each subfolder, so the sequence now reads
across folders rather than within one.

## Environment

Run through the workspace env, e.g.
`cd gfys_workspace && uv run python ../tetrode_analyses/analyses/tetrode_preprocessing_and_sorting/compression/02_bench_matrix.py`.
Scratch/experiment data live under `/nvme/neuropixels/`.
