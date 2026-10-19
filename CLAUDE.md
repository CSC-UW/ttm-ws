# Conventions

## Structure
- `packages/tetrode-analyses` is subject-agnostic. Subject- and
  session-specific facts belong in `subjects/<id>/` manifests, machine facts
  in `~/.config/ttm/config.toml`, never in library code. The library never
  imports a subject package, and `grep -rn '"/nvme\|"/Volumes' packages
  subjects` must stay empty.
- New per-subject processing is a `ttm-*` entry point in
  `tetrode_analyses.cli` that selects its experiment with
  `add_experiment_args` and takes paths from `Experiment` (`docs/pipeline.md`).
- `analyses/` is transitional: the closed method studies as run. Move them
  out per the plan; do not add new work there.
- Sync with `uv sync --all-packages --group dev`. Plain `uv sync` uninstalls
  the workspace members.
- This workspace has its own environment. Do not run its code through
  `gfys_workspace` (`uv run --project ../gfys_workspace`); that workspace
  only keeps an editable pointer to the library during the transition.
- When this directory is opened from inside `projects/ece/`, the
  `gfys_workspace/AGENTS.md` rules also load. Where they disagree on commands
  or environment, this file wins for work under this directory.

## Dependencies
- External packages are pinned by commit in the root `pyproject.toml`. To
  bump: edit `rev`, run `uv lock --upgrade-package NAME`, commit both files.
- Never edit a pinned package from this workspace. For a session of local
  hacking, `uv pip install -e ../NAME`; the next sync restores the pin.
- `torch==2.9.1` and `faiss-cpu==1.12.0` are constraints for the server GPU
  (Tesla V100) and Python 3.14 wheels; `zarr<3` and `numcodecs<0.16` are
  overrides that keep spikeinterface able to write Zarr. Do not loosen any
  of them without reading `docs/SETUP.md`, "Pinned packages".
- `spikeinterface` comes from PyPI (0.105.2), not the editable upstream clone
  `gfys_workspace` uses. `set_probegroup` is in place and returns `None` there.
- Every Zarr write goes through zarr 2 here: spikeinterface stores, and the
  xarray stores (`to_zarr`/`open_zarr`) that depend on the xarray pin.

## Git
- Commit directly to `main` in small, self-contained commits. No feature
  branches: several people and agents work this repo concurrently.
- Leave edits unstaged for the user to commit unless told otherwise.

## Data
- Production data are on `tononi-2` (`/nvme/neuropixels/tetrode_data`,
  `/Volumes/neuropixel_archive/tetrode_data`), which is what that machine's
  `~/.config/ttm/config.toml` says. `/nvme` is shared and often near full;
  the sorting cache is ~1.3 TB. Check `df -h /nvme` first.
- `experiment_params.json` on the data volume is written from the manifest by
  `ttm-experiment-params`; edit the manifest, then rewrite it.
- Bulk products never enter git. Results that are committed are small
  summaries (JSON, CSV, Markdown) next to the code that produced them.

## Docs
- Markdown notes take the front matter schema used across the lab
  workspaces: `title`, `updated`, and `status` only when not current
  (`historical | superseded | hypothesis | needs_verification`).
- Interpretive or scientific claims in tracked docs need the user's
  confirmation; flag unconfirmed text as conjecture.
- Avoid narrative comments in code; put the story in commit messages or
  `docs/developer_notes/`.

## Repository size
- No single file over 5 MB; repository under 200 MB without first discussing
  git LFS. `git count-objects -vH`.
