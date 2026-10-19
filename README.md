# ttm-ws

Tetrode (TTM) electrophysiology workspace: multi-day Open Ephys tetrode
recordings converted to sort-ready Zarr, LFP and band-power products, spike
sorting by tetrode group, cross-chunk unit tracking, and curation/visualization
hand-offs between the server and a laptop.

## Getting started

```bash
git clone https://github.com/CSC-UW/ttm-ws.git
cd ttm-ws
claude
> set up the workspace
```

The agent follows [`docs/SETUP.md`](docs/SETUP.md); that page has the same
steps for doing it by hand.

## Layout

| Path | What it is |
|---|---|
| `packages/tetrode-analyses/` | The shared library (`import tetrode_analyses`). Subject-agnostic: Open Ephys loading and channel maps, SpikeInterface recording assembly and Zarr conversion, LFP, band power, MountainSort5/Kilosort4 drivers, unit tracking, the UnitRefine advisory, SWA plotting. Its own README documents the API. |
| `subjects/<id>/` | One subject: `subject.toml`, `experiments/<NAME>.toml` manifests, notebooks, a README of what exists on the server. `ttm-001` (TTM-NOD, TTM-Test) and the `ttm-002` skeleton. |
| `investigations/` | Closed method studies (acquisition noise, Zarr compression, spike sorting and tracking, the UnitRefine advisory, the SpikeInterface memory bug): the question, the scripts as run, the findings. Indexed in `investigations/README.md`. |
| `docs/` | `SETUP.md`, `pipeline.md` (every `ttm-*` command), developer notes, reports, plans. |

## Environment

```bash
uv sync --all-packages --group dev
uv run pytest -q
```

Always pass `--all-packages`: a plain `uv sync` installs only the root and
uninstalls workspace members. External packages (`ecephys`, `loupe`, the
`mountainsort5` and `kilosort` forks) are pinned by
commit in `pyproject.toml` and resolve without credentials or sibling
checkouts. Python 3.14 or newer is required.

## Running the pipeline

```bash
uv run ttm-convert --subject TTM-001 --experiment TTM-NOD
uv run ttm-experiment-params --subject TTM-001 --experiment TTM-NOD
uv run ttm-make-lfp --subject TTM-001 --experiment TTM-NOD
uv run ttm-swa all --subject TTM-001 --experiment TTM-NOD
```

Every tool takes its experiment from a manifest under `subjects/` and every
path from the per-machine config `~/.config/ttm/config.toml`. The full list,
including the curation and visualization hand-offs to a laptop, is
[`docs/pipeline.md`](docs/pipeline.md).

## Data

Production data live on the server `tononi-2`: raw Open Ephys acquisitions
under `/Volumes/neuropixel_archive/tetrode_data/<session>/`, derived products
under `/nvme/neuropixels/tetrode_data/<session>/` and
`/nvme/neuropixels/tetrode_data/<subject>/<experiment>/`. Those two roots
are the server's config; nothing tracked names them except the closed
studies under `investigations/`, which are records of runs on that machine.

`/nvme` is a shared volume. The materialised 48 h preprocessing cache that
sorting needs is about 1.3 TB; check free space before a sort.

### Repository size

No single file over 5 MB is committed, and the repository does not grow past
200 MB, without first discussing git LFS. `git count-objects -vH` reports the
current size.

## Design

[`docs/plans/2026-10-16_ttm-ws-restructure.md`](docs/plans/2026-10-16_ttm-ws-restructure.md)
holds the restructure plan and its phase status.
