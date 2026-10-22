# Conventions

## Structure
- `packages/ttm` is subject-agnostic. Subject- and
  session-specific facts belong in `subjects/<id>/` manifests, machine facts
  in `~/.config/ttm/config.toml`, never in library code. The library never
  imports a subject package, and `grep -rn '"/nvme\|"/Volumes' packages
  subjects` must stay empty.
- New per-subject processing is a `ttm-*` entry point in
  `ttm.cli` that selects its experiment with
  `add_experiment_args` and takes paths from `Experiment` (`docs/pipeline.md`).
- `investigations/<question>/` holds closed method studies as run: literal
  server paths, no maintenance. New bake-offs get a new directory there and a
  row in `investigations/README.md`; their outcome goes into the library.
- Sync with `uv sync --all-packages --group dev`. Plain `uv sync` uninstalls
  the workspace members.
- This workspace has its own environment and is standalone. Do not run its
  code through the lab's `gfys_workspace` environment; nothing there depends
  on this workspace. The rules in this file are complete for work here.

## Dependencies
- External packages are pinned by commit in the root `pyproject.toml`. To
  bump: edit `rev`, run `uv lock --upgrade-package NAME`, commit both files.
- Never edit a pinned package from this workspace. For a session of local
  hacking, `uv pip install -e ../NAME`; the next sync restores the pin.
- `torch==2.9.1`, `faiss-cpu==1.12.0` and `xarray<2026.8` are constraints
  for the server GPU (Tesla V100), Python 3.14 wheels and zarr 2; `zarr<3`
  and `numcodecs<0.16` are overrides that keep spikeinterface able to write
  Zarr. Do not loosen any of them without reading `docs/SETUP.md`, "Pinned
  packages". Before adding a fork pin, check whether the change is already
  upstream: the spikeinterface-gui and mountainsort5 forks both were.
- `spikeinterface` comes from PyPI (0.105.2), not the editable upstream clone
  `gfys_workspace` uses. `set_probegroup` is in place and returns `None` there.
- Every Zarr write goes through zarr 2 here: spikeinterface stores, and the
  xarray stores (`to_zarr`/`open_zarr`) that depend on the xarray pin.

## Git
- Commit directly to `main` in small, self-contained commits. No feature
  branches: several people and agents work this repo concurrently.
- Leave edits unstaged for the user to commit unless told otherwise.
- Never push, and never open or edit a pull request, unless the user asks
  for it in that conversation.
- Review a commit message as a hostile reader before proposing it: tense
  (present for behaviour that persists, past for what is removed,
  imperative for the fix), "so" only for causation, no vague time words,
  no "every" where a partial case exists. Body lines under 72 characters.

## Data
- Production data are on `tononi-2` (`/nvme/neuropixels/tetrode_data`; raw
  acquisitions on `/Volumes/tetrode/<SUBJECT>/`), which is what that
  machine's `~/.config/ttm/config.toml` says. Writing to `/Volumes/tetrode`
  needs the `paxilline` group (`sg paxilline -c '...'` in a shell that
  predates the membership). `/nvme` is shared and often near full;
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
- Everything written into the repo must read as the user's own writing:
  no bold-label sentence openers or bold list labels, no rhetorical
  italics, no em-dashes, no "[CONFIRM]" blocks, no closing offers.
- Code behaviour and environment facts go in `docs/developer_notes/`,
  completed investigations in `investigations/` or `docs/reports/`,
  future work in `docs/plans/`. Machine facts go in `docs/SETUP.md` or
  the per-machine config, never in this file or the library.

## Analysis practice
- Compare, don't assert. When two methods compete, run them head to head
  on the data and show what each admits or produces, not only summary
  metrics. Before calling something an open problem, check whether an
  earlier script in `investigations/` already prototyped the answer.
- An approximation (a windowed preview, a subsample, a cached proxy) is
  adopted only after its agreement with the exact result is measured.
- Measure before diagnosing slowness and before claiming a speed-up;
  report the before and after numbers together.
- Anything tied to a sorting (curation, re-attached traces, metrics) uses
  only the recording store the sort was produced from, never a lossy or
  alternative variant of the same session.
- When a data source's validity is uncertain, build around it so it can be
  swapped, and label results derived from it as provisional.
- Compute added to an existing analysis gets an opt-out toggle, and
  repeated random draws are capped at the number of distinct outcomes.
- A window inside a stimulus is never a "baseline" or "quiet"; name it by
  its position relative to onset and offset and say what it may contain.
- Destructive actions (deleting data, killing processes, overwriting
  stores) need the user's explicit confirmation; look at the target first.

## Repository size
- No single file over 5 MB; repository under 200 MB without first discussing
  git LFS. `git count-objects -vH`.

## Current work
- The TTM-003 48 h spike-sorting plan, written 2026-10-22 for the agent
  that continues after the workspace moves to `~/projects/ttm-ws/`:
  `docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`. Read it after this
  file. It carries the TTM-003 short-session findings, the TTM-001
  pipeline's state, the data and machine facts, and the working rules
  that used to come from the parent workspace.
