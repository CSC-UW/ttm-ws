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
| `analyses/` | Transitional. The per-session studies and drivers as they stood before the restructure. They move to `subjects/` and `investigations/` in a later phase (see the plan below); until then their READMEs still say to run through `gfys_workspace`, which no longer applies: run them from this workspace's environment. |
| `subjects/<id>/` | One subject's analysis: manifests, notebooks, committed result summaries. Not yet populated. |
| `investigations/` | Method bake-offs whose product was a decision baked into the library. Not yet populated. |
| `docs/` | Setup, developer notes, reports, plans, acceptance baselines. |

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

## Data

Production data live on the server `tononi-2`: raw Open Ephys acquisitions
under `/Volumes/neuropixel_archive/tetrode_data/<session>/`, derived products
under `/nvme/neuropixels/tetrode_data/<session>/` and
`/nvme/neuropixels/tetrode_data/<subject>/<experiment>/`. Until the
per-machine config lands (plan, Phase 2), those paths are literal in the
scripts under `analyses/`.

`/nvme` is a shared volume. The materialised 48 h preprocessing cache that
sorting needs is about 1.3 TB; check free space before a sort.

### Repository size

No single file over 5 MB is committed, and the repository does not grow past
200 MB, without first discussing git LFS. `git count-objects -vH` reports the
current size.

## Design

[`docs/plans/2026-10-16_ttm-ws-restructure.md`](docs/plans/2026-10-16_ttm-ws-restructure.md)
holds the restructure plan and its phase status.
