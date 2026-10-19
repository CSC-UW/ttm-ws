---
title: Setting up the workspace
updated: 2026-10-19
---

# Setting up the workspace

The short way: clone, `cd` in, launch Claude Code, and say "set up the
workspace". The agent follows `.claude/skills/setup-workspace/SKILL.md`,
which is the same procedure as below.

## 1. Prerequisites

- git
- [uv](https://docs.astral.sh/uv/) 0.9 or newer: `curl -LsSf https://astral.sh/uv/install.sh | sh`
- Python 3.14 or newer. uv downloads one if none is found: `uv python install 3.14`

## 2. Environment

```bash
uv sync --all-packages --group dev
```

Always pass `--all-packages`. A plain `uv sync` installs only the root
project and uninstalls the workspace members.

The first sync downloads torch (about 1 GB on Linux). Pinned packages resolve
from GitHub over HTTPS; no sibling checkouts and no credentials are needed.

The environment deliberately runs zarr 2 (see "Pinned packages"). Do not
"upgrade" it to zarr 3: spikeinterface can read but not write Zarr with
zarr 3, and the conversion and LFP steps write Zarr.

## 3. Tests

```bash
uv run pytest -q
```

Two tests are marked `requires_nfs` and skip themselves when the production
paths are not mounted. On `tononi-2` they run.

## 4. GPU

```bash
uv run python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
```

On `tononi-2` (Tesla V100) this must report a CUDA 12 build and `True`.
On a laptop without CUDA it reports a CPU build and `False`; everything except
Kilosort4 works there.

## 5. Data and the per-machine config

Nothing in the library or the `ttm-*` tools names a data path. Each machine
says where its data are in `~/.config/ttm/config.toml` (or the file named by
`$TTM_CONFIG`):

```toml
data_root = "/nvme/neuropixels/tetrode_data"

[remotes]
archive = "/Volumes/neuropixel_archive/tetrode_data"
```

`data_root` holds `<session>/` (Zarr stores, LFPs, sortings, bundles) and
`<subject>/<experiment>/` (`experiment_params.json`, band powers, EMG).
`remotes.archive` is where raw Open Ephys acquisitions are mounted; it is
needed only to convert. The values above are the server, `tononi-2`.

A laptop that receives server-built bundles instead sets `data_root` to the
external drive and names the server as an rclone remote:

```toml
data_root = "/Volumes/Drive/ttm"

[remotes]
server = "tononi-2:/nvme/neuropixels/tetrode_data"
```

Environment variables override the file: `TTM_DATA_ROOT`, `TTM_CACHE_DIR`,
`TTM_REMOTE_ARCHIVE`, `TTM_REMOTE_SERVER`. There is no default; a tool run
without a `data_root` stops with the template above.

Check the config with:

```bash
uv run python -c "from tetrode_analyses import config; print(config.get())"
```

The subject to work on is never defaulted either: every tool takes
`--manifest`, or `--subject` and `--experiment`, or `TTM_SUBJECT` and
`TTM_EXPERIMENT`. Manifests live in `subjects/<id>/experiments/`.

## 6. Where to read next

`README.md`, then `CLAUDE.md`, then `docs/pipeline.md` (the per-subject
pipeline and every `ttm-*` command), then `packages/tetrode-analyses/README.md`,
then `analyses/tetrode_preprocessing_and_sorting/README.md`.

## Pinned packages

| Package | Pin | Why this and not PyPI |
|---|---|---|
| `ecephys` | `CSC-UW/ecephys` `develop` | Lab signal-processing library; not on PyPI. Supplies bandpass, Hilbert, decimation, STFT for `tetrode_analyses.power`. |
| `loupe` | `CSC-UW/loupe` `main` | The local viewer used by the visualization bundle. Not on PyPI. |
| `mountainsort5` | `grahamfindlay/mountainsort5` `feature/deterministic-pca` | Seeds the PCA solver so repeat sorts are bit-identical; upstream 0.5.8 is not deterministic. |
| `kilosort` | `grahamfindlay/Kilosort` `perf/clustering_chunk_size` | Adds `clustering_chunk_size`, which lets 48 h tetrode clustering fit in 32 GB of GPU memory. Stock 4.1.7 rejects the key. |
| `spikeinterface-gui` | PyPI `>=0.14.0` | The lab's fork is retired: its eleven commits (the `--recording` CLI fixes and the desktop performance work) were merged upstream in PRs 274 and 275, and 0.14.0 also carries the int32 scrollbar fix for recordings longer than 2^31 samples. |
| `torch` | `==2.9.1` (constraint) | The CUDA 12 build; CUDA 13 wheels are compiled for sm_75 and newer and fail on the V100. |
| `faiss-cpu` | `==1.12.0` (constraint) | The only release with a Python 3.14 wheel. |
| `zarr` | `>=2.18,<3` (override) | spikeinterface 0.105 requires zarr 3 on Python 3.14 but raises `NotImplementedError` on every Zarr write with it. zarr 2 installs on 3.14; the override removes spikeinterface's marker-conditional requirement. `gfys_workspace` runs zarr 3 and therefore cannot run the converters. |
| `numcodecs` | `<0.16` (override) | zarr 2 needs the pre-0.16 blosc API (`cbuffer_sizes`). |
| `xarray` | `<2026.8` (constraint) | 2026.9 dropped zarr 2 (its backend passes `zarr_format=`, which zarr 2 rejects); 2026.7.0 is the newest release that reads and writes with the zarr override. Every xarray store here (band powers, EMG, viz bundle) depends on it. |

To work on one of them for a session:

```bash
uv pip install -e ../NAME
```

The next `uv sync --all-packages` restores the pin. To change the pin for
everyone: edit `rev` in the root `pyproject.toml`, run
`uv lock --upgrade-package NAME`, then `uv sync --all-packages`, and commit
both files.
