---
title: Setting up the workspace
updated: 2026-10-22
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

After the checkout is moved (as `~/projects/ece/ttm-ws` became
`~/projects/ttm-ws` in October 2026) the console scripts in `.venv/bin` keep
the old interpreter path in their shebang and `uv run pytest` fails with
"Failed to spawn"; `uv run python -m pytest` works, and the next
`uv sync --all-packages --group dev` rewrites the scripts.

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
archive = "/Volumes/tetrode"
```

`data_root` holds `<session>/` (Zarr stores, LFPs, sortings, bundles) and
`<subject>/<experiment>/` (`experiment_params.json`, band powers, EMG).
`remotes.archive` is where raw Open Ephys acquisitions are mounted, laid out
as `<SUBJECT>/<session>/Record Node */` (`TTM1`, `TTM2`, `TTM3`, and
`Unknown` for sessions not yet attributed to a subject); each manifest's
`acquisition` is a path under it. Remotes are needed only to convert. The
values above are the server, `tononi-2`.

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
uv run python -c "from ttm import config; print(config.get())"
```

The subject to work on is never defaulted either: every tool takes
`--manifest`, or `--subject` and `--experiment`, or `TTM_SUBJECT` and
`TTM_EXPERIMENT`. Manifests live in `subjects/<id>/experiments/`.

### tononi-2

The server's facts, dated 2026-10-22, for anyone sharing it:

| Volume | Size | Free | Use |
|---|---|---|---|
| `/nvme` | 5.8 TB | 2.1 TB | `data_root`. TTM-001's session directory alone is 2.2 TB (store plus 1.3 TB sorting cache). Often near full; `df -h /nvme` before writing. |
| `/ssd-raid0` | 7.0 TB | 6.0 TB | Allowed for stores and sorting caches. Slower than `/nvme`, much faster than spinning disk. The TTM-003 48 h session directory `2026-10-20_09-02-12` lives at `/ssd-raid0/neuropixels/tetrode_data/` and is reached through a symlink of the same name under `data_root`, so no tool needs a second root. |
| `/Volumes/tetrode` | CIFS | | `remotes.archive`. Writing needs the `paxilline` group (`sg paxilline -c '...'` in a shell that predates the membership) and `rsync --inplace`. |

Until the 48 h session's `structure.oebin` on `/Volumes/tetrode` is repaired in
place (`docs/developer_notes/openephys_blank_channel_names.md`), its store was
converted from `<session>/oe_mirror/`, a symlink mirror of the acquisition with
a repaired copy of the oebin, selected with
`TTM_REMOTE_ARCHIVE=/ssd-raid0/neuropixels/tetrode_data/2026-10-20_09-02-12/oe_mirror`;
the store's provenance records the mirror's `settings.xml` path. The mirror can
go once the in-place repair is done.

Other machines use it at the same time: never kill processes by pattern
(`pkill -f` matches the shell that runs it), find PIDs with
`pgrep -f "pat[t]ern"`, walk the parent chain, and kill by PID only when
the owner is known. Use absolute paths and `git -C` in every command; a
shell's working directory does not persist between agent calls. The GPU is
one Tesla V100 with 32 GB (section 4).

## 6. Where to read next

`README.md`, then `CLAUDE.md`, then `docs/pipeline.md` (the per-subject
pipeline and every `ttm-*` command), then `packages/ttm/README.md`,
then `investigations/README.md`. Tool limits that shaped the pipeline:
`docs/developer_notes/sorting_tool_limits.md`.

## Pinned packages

| Package | Pin | Why this and not PyPI |
|---|---|---|
| `ecephys` | `CSC-UW/ecephys` `develop` | Lab signal-processing library; not on PyPI. Supplies bandpass, Hilbert, decimation, STFT for `ttm.power`. |
| `loupe` | `CSC-UW/loupe` `main` | The local viewer used by the visualization bundle. Not on PyPI. |
| `mountainsort5` | PyPI `>=0.5.9` | The lab's deterministic-PCA fork was merged upstream (flatironinstitute/mountainsort5 PR 65) and released in 0.5.9, so repeat sorts are bit-identical with the stock package. The `whitening_seed` half of that fix is spikeinterface PR 4600, in 0.105.0. |
| `kilosort` | `grahamfindlay/Kilosort` `perf/clustering_chunk_size` | Adds `clustering_chunk_size`, which lets 48 h tetrode clustering fit in 32 GB of GPU memory. Stock 4.1.7 rejects the key. Rebased onto upstream main (4.1.7 plus the swarmsplitter CCG fix from PR #1043) on 2026-10-19; upstream PR MouseLand/Kilosort#1037 (opened 2026-06-20) is still open. |
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
