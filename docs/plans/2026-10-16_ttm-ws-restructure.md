---
title: Restructure tetrode_analyses into the ttm-ws machine-portable workspace
updated: 2026-10-19
---

# Restructure `tetrode_analyses` into `ttm-ws`

Decisions 1, 2, 6 and 8 accepted by the user 2026-10-16, with the `ttm`
import rename made the first follow-up (Phase 6). Modelled on
`grahamfindlay/testi-ripples-ws` (private; read 2026-10-16 at its 2026-10-07
head). Phase status is tracked inline below.

## Verdict

Appropriate, and timely. Three facts measured today argue for doing it now
rather than later:

1. A second subject exists. `/nvme/neuropixels/tetrode_data/TTM-002/TTM-002-SD`
   (session `2026-09-24_09-03-22`, 48.23 h, 64 ch) was converted, QC'd and
   LFP'd on 2026-09-30..10-02 by `txiao37`, from a standalone clone at
   `/Users/txiao37@ad.wisc.edu/tetrode/tetrode_analyses/`. The repo's own
   layout still assumes one session, `2026-05-27_09-07-52` (TTM-001/TTM-NOD).
2. That clone carries uncommitted, unpushed changes (no fork, no branch on
   GitHub; 2 commits behind `origin/main`). Decision 2026-10-16: this work is
   not adopted. TTM-002 will be re-analysed from the archive with our own
   pipeline, and the TTM-002 products currently on `/nvme` are not treated
   as a baseline for anything. The restructure proceeds without waiting on
   that clone; if its owner later pushes, rebasing onto the new layout is
   their task.
3. The external-clone install mode is in real use and is under-specified. The
   standalone `uv.lock` resolves `mountainsort5`, `kilosort`, `torch`,
   `spikeinterface` from PyPI; the forks (`mountainsort5@feature/deterministic-pca`,
   `Kilosort@perf/clustering_chunk_size`), `torch==2.9.1` (CUDA 12 for the
   V100) and `faiss-cpu==1.12.0` live only in `gfys_workspace/pyproject.toml`.
   A collaborator's clone therefore runs a different sorter than the one the
   findings docs describe. A workspace root is the right place for those pins.

## What the reference workspace does

Distilled from `testi-ripples-ws` (README, CLAUDE.md, root `pyproject.toml`,
`docs/SETUP.md`, `.claude/skills/setup-workspace/SKILL.md`,
`docs/superpowers/specs/*`):

| Convention | Mechanism |
|---|---|
| One git repo, one `uv.lock` | Root `[tool.uv.workspace] members = ["packages/*", "subjects/*", "investigations/*"]`; always `uv sync --all-packages` |
| Shared library is subject-agnostic | `packages/<lib>/`; the grep that must stay empty is an import of a subject package |
| One directory per subject | `subjects/<id>/{pyproject.toml, subject.toml, src/<id>/ingest.py, notebooks/, results/, README.md}`; dependency arrow `subjects/* -> lib -> external` only |
| Method bake-offs separated from subject results | `investigations/<question>/` with a README table: question, outcome, decision baked into the library |
| Machine facts live outside the repo | `~/.config/<name>/config.toml` (`data_root`, `[remotes]`), env overrides, no default, no hostname sniffing; `lib.config` raises if unset |
| No default subject | `--manifest` > `--subject` > `$ENV` > error |
| External tools pinned, not vendored | Root `[tool.uv.sources]` with `git = ..., rev = <sha>` over HTTPS; temporary local hacking via `uv pip install -e ../X` |
| Agent-first onboarding | `README` says `claude` > "set up the workspace"; a committed project skill walks prerequisites, sync, tests, config, data; `docs/SETUP.md` is the human copy; `.claude/settings.json` allowlists `uv`/`git` read commands |
| Size rule | No file > 5 MB, repo < 200 MB without discussing LFS |
| Acceptance baseline | `scripts/snapshot_derived.py --compare` hashes derived tables before/after, with a README stating exactly what was and was not re-run |

## Current state of `tetrode_analyses` (measured)

- 197 tracked files; 125 of them under
  `analyses/tetrode_preprocessing_and_sorting/sorting/` (scripts 10–98,
  helpers, findings docs, JSON/CSV results). Library: 11 modules in
  `src/tetrode_analyses/`. Tests: 6 files, 2 marked `requires_nfs`.
- Loose object store 11.9 MiB; largest blobs are four ~1.3 MB notebooks.
  `sorting/refs/` (74 MB of PDFs) and all `*.log`/`*.png` are gitignored.
  Well inside the size rule.
- Two acquisition sessions appear in code: `2026-05-19_18-09-51` (artifact
  notebooks, `test_channel_map`) and `2026-05-27_09-07-52` (everything else).
  The archive now holds nine session directories.
- Machine coupling, by count of literal paths in tracked `.py`/`.ipynb`/`.md`:
  97 × `"/nvme/neuropixels/tetrode_data`, 27 × `"/nvme/neuropixels/tmp`,
  6 × `"/Volumes/neuropixel_archive/tetrode_data`, plus `tononi-2` and an
  rclone remote named `tononi-2` in `curation/` and `visualization/`. In the
  library itself only two: `experiment.DEFAULT_DATA_ROOT` and
  `unitrefine_advisory`'s `TETRODE_MODELS_DIR`/`~/.cache` fallback.
- Subject coupling in the library: none by name. `experiment.py` is already
  keyed by `(subject, experiment)`; the session id is only linked through
  `experiment_params.json`. The three `experiment_params.json` files on disk
  hold absolute `/nvme/...` paths (`openephys_session`, `lfp_zarr`).
- Dependency wiring: `gfys_workspace` consumes it as an editable path override
  with extras `[analysis,sorting,curation,advisory]`; nothing else in the
  workspace imports it. `power.py` imports `ecephys.npsig`/`ecephys.xrsig.core`;
  one script (`67_np_spike_coverage.py`) imports `wisc_ecephys_tools`.
- Workspace docs that name its paths and must follow a move: `AGENTS.md` tier
  list; 4 registry entries (`subsampled_lfp`, `loupe_viz_bundle`,
  `mp_production_sorting`, `experiment_params`); 5 developer notes/reports.
- All upstream repos needed for HTTPS pins are public (`CSC-UW/{tetrode_analyses,
  ecephys, loupe, wisc_ecephys_tools, emg_from_lfp}`, `grahamfindlay/{mountainsort5,
  Kilosort, spikeinterface-gui}`). Neither `grahamfindlay/ttm-ws` nor
  `CSC-UW/ttm-ws` exists yet.

## Target layout

```
ttm-ws/
  pyproject.toml            # workspace members; git pins; torch/faiss constraints
  uv.lock
  README.md  CLAUDE.md  .gitignore
  .claude/settings.json  .claude/skills/setup-workspace/SKILL.md
  docs/SETUP.md  docs/developer_notes/  docs/reports/  docs/plans/  docs/baselines/
  packages/tetrode-analyses/            # the library; import name unchanged
    pyproject.toml  src/tetrode_analyses/  tests/
  subjects/ttm-001/                     # TTM-NOD
    pyproject.toml  subject.toml  experiments/TTM-NOD.toml
    src/ttm001/  notebooks/  results/  README.md
  subjects/ttm-002/                     # TTM-002-SD + the three test experiments
    ...
  investigations/
    acquisition-noise/        # 5-19 + spectral_and_artifact + subj01-artifact notebooks + the 05-30 report
    zarr-compression/         # compression/ (01–08, 17, results, README)
    ms5-ks4-sorting/          # sorting/ 10–45, findings docs, helpers
    unit-tracking/            # or kept inside ms5-ks4-sorting; TRACKING_FINDINGS.md
    unitrefine-advisory/      # 46–50
    matching-pursuit/         # 51–98, MATCHING_PURSUIT_FINDINGS.md, OPEN_CONCERNS.md
    si-frame-slice-memory/    # unchanged
  scripts/snapshot_derived.py
```

What moves where, by current directory:

| Now | Goes to | Why |
|---|---|---|
| `src/`, `tests/`, `pyproject.toml` | `packages/tetrode-analyses/` | the subject-agnostic library |
| `analyses/swa/`, `analyses/emg/` | library CLI entry points (`ttm-swa`, `ttm-emg`) driven by a manifest; TTM-001 outputs recorded under `subjects/ttm-001/results/` | they are already generic over `(subject, experiment)` |
| `analyses/tetrode_preprocessing_and_sorting/{curation,visualization}/` | library CLI entry points (`ttm-curate`, `ttm-build-viz-bundle`, `ttm-fetch-bundle`, `ttm-loupe`) | server/laptop hand-off tools, not a study |
| `analyses/tetrode_preprocessing_and_sorting/lfp/` | library (`make_lfp` already is); scripts become `ttm-make-lfp` | same |
| `compression/`, `sorting/`, `si_frame_slice_memory/` | `investigations/*` as listed above, moved with `git mv` only | closed studies whose products are decisions and findings docs |
| `97_production_pipeline.py`, `98_build_analyzer_prod.py` | stay with `matching-pursuit/`; the registry entry and `subjects/ttm-001/README.md` point at them | the run is a TTM-001 result, the code is a method |
| `5-19/`, `analyses/investigate_spectrum_and_artifact/`, `analyses/subj01/*artifact*` | `investigations/acquisition-noise/` | one question, one report |
| `analyses/subj01/subj01.ipynb` | `subjects/ttm-001/notebooks/` | subject exploration |
| (nothing yet) | `subjects/ttm-002/` with manifests written from the archive session and our own `experiment_params.json` | TTM-002 is re-processed by our pipeline; existing on-disk products are not inputs |

## Design decisions (recommendation first)

1. Repo identity. Rename `CSC-UW/tetrode_analyses` to `CSC-UW/ttm-ws` in
   place and restructure on `main` with `git mv`. History is preserved,
   GitHub redirects the old URL, and the collaborator's remote keeps working.
   The alternative (new repo + subtree import, as testi did) buys nothing here
   because there is one repo to fold in, not several.
2. Library name. Keep import name `tetrode_analyses` and dist name
   `tetrode-analyses`; the workspace is `ttm-ws`. The reference kept its
   library name unchanged for the same reason: 125 scripts, four registry
   entries and five docs reference it, and a rename adds churn without
   changing a boundary. Renaming to `ttm` is Phase 6, the first follow-up
   once the layout has settled.
3. Per-machine config. `~/.config/ttm/config.toml` (or `$TTM_CONFIG`) with
   `data_root` and `[remotes] archive = "/Volumes/tetrode"` (the raw sessions moved
   there from `/Volumes/neuropixel_archive/tetrode_data` on 2026-10-19)
   (plus `scratch`, `tmp` as needed); env overrides `TTM_DATA_ROOT`,
   `TTM_REMOTE_<NAME>`; no default. `experiment.DEFAULT_DATA_ROOT` is deleted.
   On tononi-2 `data_root = /nvme/neuropixels/tetrode_data`.
4. Manifest. Keep the WNE hierarchy: `subjects/<subject>/subject.toml` for
   subject-level facts, `subjects/<subject>/experiments/<EXPERIMENT>.toml` per
   experiment, seeded from the on-disk TTM-001 `experiment_params.json` and
   the `ExperimentParams` dataclass. Paths inside manifests are relative to
   `data_root`/a remote; `experiment_params.json` keeps being written to the
   data volume (it is a data product other code reads) but with relative paths.
   Resolution: `--manifest` > `--subject/--experiment` > `$TTM_SUBJECT`/
   `$TTM_EXPERIMENT` > error.
5. Pins. Root `[tool.uv.sources]`: `ecephys` (git rev on `develop`), `loupe`,
   `mountainsort5` fork, `kilosort` fork, `spikeinterface-gui` fork,
   `wisc_ecephys_tools` (for script 67 only; consider dropping that import).
   Root `[tool.uv] constraint-dependencies`: `torch==2.9.1`, `faiss-cpu==1.12.0`,
   moved out of `gfys_workspace`. `spikeinterface` needs a decision: the
   workspace uses the local editable `ece/spikeinterface` (upstream main plus
   local fixes); the standalone lock uses PyPI. Pin a PyPI release ≥ 0.106 and
   run the suite; if a local fix is still required, pin the upstream git rev
   that has it.
6. `gfys_workspace` relationship. During the transition keep the editable
   override, pointing at `../ttm-ws/packages/tetrode-analyses`; once
   `ttm-ws`'s own environment reproduces the production validation commands,
   drop `tetrode_analyses` and the `tetrodes` extra from `gfys_workspace`
   entirely. Its tier-3 line in `AGENTS.md` then says the project lives in its
   own workspace.
7. De-literalising paths. Library, subject packages and CLI tools: zero
   literal data paths (gated by `grep -rn '"/nvme\|"/Volumes' packages subjects`
   being empty). Closed investigations: moved as-is, with each README gaining
   a line that paths are as-run on tononi-2 and the directory is a record.
   Rewriting ~100 literals in closed studies would be churn on code nobody
   will run again.
8. Data fetch. Do not port `fetch.toml`/`fetch-data`. The existing
   `download_bundle.py` and `download_curation_bundle.py` are the real need
   (server-built bundle, rclone to a laptop); fold them into one
   `ttm-fetch-bundle` reading the remote name from config. Revisit if raw
   Open Ephys data ever needs to flow to laptops.

## Phases

Each phase is one or a few commits on `main` (the package's existing rule:
no feature branches, small commits, because several people and agents work
it concurrently). The user commits; edits are left in the working tree.

### Phase 1: workspace skeleton

Status 2026-10-16: done on disk, uncommitted. Library moved with `git mv`;
root and package `pyproject.toml` written; `uv.lock` resolves 278 packages
with the forks, `spikeinterface` 0.105.2 from PyPI, `zarr` 3.4.1, the
CUDA 12 `torch` 2.9.1; `requires-python` raised to 3.14 because
`spikeinterface` 0.105.2 pins `zarr<3` on 3.13. README, CLAUDE.md,
`docs/SETUP.md`, `.claude/settings.json` and the setup skill written. The
GitHub rename and the local directory rename (`tetrode_analyses` ->
`ttm-ws`) are left to the user: both affect other sessions and the remote.
Three tests located the sorting helpers via `parents[1]`; now `parents[3]`,
to be revisited when Phase 3 moves those helpers.

Acceptance: a copy of the working tree in a directory with no sibling
checkouts syncs in 56 s (warm cache) and passes 50/50 tests, forks and the
CUDA 12 torch resolved. `gfys_workspace` relocked (path-only diff plus the
pins uv now reads from this root) and still imports the library; its venv
needed no reinstall.

Two pre-existing failures surfaced and were fixed here because the
acceptance test needs a green suite. (a) `set_probegroup` is in place in
spikeinterface 0.105 and returns `None`; `io.attach_tetrode_probegroup`
assigned its result. (b) spikeinterface requires zarr 3 on Python 3.14 but
raises `NotImplementedError` on every Zarr write with it, so the converters
and LFP writer were dead in any zarr 3 environment, `gfys_workspace`
included. Root `override-dependencies` now force `zarr<3` and
`numcodecs<0.16`, which zarr 2 needs; 5 previously failing tests pass.
`ruff check packages/` reports 42 findings, identical before and after the
move; not addressed.

2026-10-19: the `spikeinterface-gui` fork pin is dropped for PyPI 0.14.0.
All eleven fork commits were merged upstream verbatim (PRs 274 and 275), and
0.14.0 includes the int32 scrollbar fix (PR 263). Suite still 50/50; the
curation launchers' `run_mainwindow` call signature is unchanged, but the
`grahams_curation` layout/settings presets have not been exercised against
0.14.0's new merge-focus and curation-panel views. The zarr override stays:
upstream main (0.106.0, unreleased) carries the same marker-conditional pins
and the `requires_zarr_write` test marker, and zarr v3 write support is the
open PR SpikeInterface/spikeinterface#4260 (opened 2025-12-11, still under
review 2026-10-19). `gfys_workspace` still pins the sigui fork for `sigui2`.

2026-10-19: GitHub repo renamed to `CSC-UW/ttm-ws` by the user; local
directory renamed `tetrode_analyses` -> `ttm-ws`; `gfys_workspace` paths
follow (`../ttm-ws/packages/tetrode-analyses`, commit d2c1c6e). A uv venv
does not survive a directory move (absolute shebangs and `pyvenv.cfg`), so
`.venv` was recreated. Phase 1 complete.

- Rename the GitHub repo; update the local checkout's directory name and
  `gfys_workspace` source/override paths.
- Root `pyproject.toml` (members, sources, constraints, `testpaths`),
  `git mv src tests pyproject.toml README.md -> packages/tetrode-analyses/`,
  package `pyproject.toml` loses its direct-git sibling URLs (root sources
  now apply to every member).
- `CLAUDE.md`, `README.md`, `docs/SETUP.md`, `.claude/settings.json`,
  `.claude/skills/setup-workspace/SKILL.md`, `.gitignore`.
- Acceptance: in a scratch directory with no sibling checkouts,
  `git clone` + `uv sync --all-packages` + `uv run pytest -q` passes, and the
  resolved `mountainsort5`/`kilosort`/`torch` are the fork and the CUDA-12
  build. `gfys_workspace` still syncs.

### Phase 2: config and manifest

Status 2026-10-19: built, uncommitted. `tetrode_analyses.config` and
`tetrode_analyses.subject` as specified; `ExperimentParams` is now the
manifest `Experiment` plus `t0_unix`, with `session`, `acquisition` and
`oe_stream_index` in place of the absolute `openephys_session`/`lfp_zarr`
(old files are rejected on load). Ten entry points under
`tetrode_analyses.cli` (`ttm-convert`, `ttm-experiment-params`,
`ttm-make-lfp`, `ttm-swa`, `ttm-emg`, `ttm-curation-manifest`, `ttm-curate`,
`ttm-viz-bundle`, `ttm-loupe`, `ttm-fetch-bundle`) replace `swa/`, `emg/`,
`lfp/`, `compression/06+08`, `curation/*.py` and `visualization/*.py`, which
are deleted; the curation presets and loupe state definitions ship in
`tetrode_analyses/assets/`. `subjects/ttm-001` (TTM-NOD, TTM-Test) and the
`subjects/ttm-002` skeleton are workspace members. `docs/pipeline.md` is the
operator document. The TTM-001 `slice_table.csv` had the pre-June column
names (`experiment_name`, `software_time_s`); renamed in place to the
current `oe_*` schema, which `read_session_t0_unix` now reads.

Acceptance: `ttm-experiment-params` reproduces the TTM-001 light/dark
periods and anchor; `ttm-curation-manifest` reproduces the June manifest
except `created`; `ttm-viz-bundle` into a scratch directory reproduces the
June bundle manifest (81,276,400 spikes, 232 units, same separator
boundaries); `ttm-swa plot` and `ttm-emg --plot-only` regenerate the TTM-001
figures. Seven `gfys_workspace` registry entries updated (producers,
consumers, the `experiment_params` schema).

Found on the way: xarray 2026.9 (locked in Phase 1) cannot open or write
Zarr with zarr 2, so every xarray store was unreachable from this
environment and no test covered it. `xarray<2026.8` is now a root
constraint (2026.7.0 verified). TTM-002's archived `structure.oebin` names
all 64 acquisition-board channels with a blank string and neo refuses to
parse the session; `ttm-convert` needs a channel-name repair before TTM-002
can be processed (noted in `subjects/ttm-002/README.md`).

- `tetrode_analyses.config` (data_root, remotes, env overrides, error on
  missing). Delete `DEFAULT_DATA_ROOT`; `unitrefine_advisory` cache dir reads
  config with the `~/.cache` fallback kept (it is a cache, not data).
- `tetrode_analyses.subject` (`Subject`, `Experiment`, discovery of
  `subjects/*/`, resolution order above).
- `subjects/ttm-001/` with manifests; regenerate its on-disk
  `experiment_params.json` with relative paths and update the
  `experiment_params` registry entry's schema. `subjects/ttm-002/` gets a
  manifest skeleton pointing at the archive session, filled in when we
  re-process it.
- Convert `swa/`, `emg/`, `lfp/`, `compression/06_convert.py`,
  `08_convert_session.py`, `curation/` and `visualization/` into
  manifest-driven CLI entry points. Together these are the end-to-end
  pipeline a new subject needs (convert, LFP, experiment params, SWA, EMG,
  bundles).

### Phase 3: investigations

Status 2026-10-19: done. `investigations/{acquisition-noise, zarr-compression,
spike-sorting, unitrefine-advisory, si-frame-slice-memory}`, moved with
`git mv`; `investigations/README.md` carries the question/outcome/decision
table. The tracking and matching-pursuit threads stay together in
`spike-sorting/` because four matching-pursuit scripts import `_track_eval`
alongside the tracking scripts; only the self-contained UnitRefine scripts
(46-50) were split out. The 2026-05-30 noise report moved from
`docs/reports/` into `acquisition-noise/`. `analyses/` and `5-19/` are gone.

`git mv` the closed studies into `investigations/<question>/`, write
`investigations/README.md` (question / outcome / decision table; the material
is already in the findings docs), add the as-run note to each README.

### Phase 4: workspace follow-through

`gfys_workspace/AGENTS.md`, the four registry entries, the five
developer notes/reports, and the memory topic files that cite
`tetrode_analyses/analyses/...` paths. Then decide whether `gfys_workspace`
keeps a dependency at all (recommendation: no).

Status 2026-10-20: done. `gfys_workspace` dropped the `tetrodes` extra, the
`tetrode_analyses` editable source and the mountainsort5, kilosort and
spikeinterface-gui fork overrides (gfys commit 6b1cf92); torch, faiss and
the CUDA wheels left its environment, and sigui2 resolves upstream
spikeinterface-gui 0.14.0. AGENTS.md names `ttm-ws` as a separate
workspace. Registry entries and developer notes were repointed earlier
(gfys a764031 through 33a184b).

### Phase 5: acceptance baseline

A small `scripts/snapshot_derived.py` that hashes, for TTM-001: the
production sorting's unit count and `is_mua` split (the registry's
validation command), `experiment_params.json` content, and summary stats
of `delta.idelta.zarr` and `stft_bandpowers.zarr`; snapshot before Phase 2,
compare after Phase 2 and after re-running the SWA entry points on TTM-001.
`docs/baselines/README.md` states what was re-run and what was only moved.

### Phase 6 (first follow-up): rename the import to `ttm`

`tetrode_analyses` -> `ttm` (dist `ttm`, directory `packages/ttm/`). One
mechanical commit: `git mv`, `sed` over imports in the library, subject
packages, investigations and docs, registry and workspace-doc references,
then the Phase 5 snapshot compare. Done after Phases 2 to 5 so that scripts
are rewritten once, not twice.

## Risks and open items

- A standalone clone with unpushed edits to `experiment.py`/`io.py` exists
  (see Verdict). It is deliberately ignored; the only exposure is a later
  push to `main` that conflicts, which review would catch.
- TTM-002's current `/nvme` products (blosc store, LFP, idelta,
  `experiment_params.json`) were made outside this repo's tracked code and
  are superseded when we re-process. Disk is tight (see below), so decide
  when to delete them.
- The zarr 2 environment also pins xarray (`<2026.8`); the two move together
  until spikeinterface writes zarr 3.
- `spikeinterface` pin: resolved in Phase 1. PyPI 0.105.2 passes the full
  suite, including the two `requires_nfs` tests on tononi-2, after the
  `set_probegroup` fix. Upstream main pins zarr the same way, so a git pin
  would not have avoided the zarr override.
- zarr split: this workspace runs zarr 2, `gfys_workspace` zarr 3. Zarr
  stores written here are v2 and readable by both. Code in this workspace
  must not use zarr 3 APIs (`zarr.storage.LocalStore`, `create_array`); the
  compression benchmark scripts that call `zarr.DirectoryStore` work again
  here.
- `mountainsort5` reads `recording._recording_segments`; current SI keeps that
  as a backward-compatibility property, so the old venv patch is no longer
  needed. Worth a one-line developer note.
- `torch==2.9.1` is pinned for the V100. On a GPU-less laptop it resolves the
  CPU wheel; untested. The setup skill should run a `torch.cuda.is_available()`
  report and say which is expected.
- `/nvme` is a shared 5.8 TB volume at 90 % full (measured today). Not a
  restructure problem, but `docs/SETUP.md` must say the materialised 48 h
  preprocessing cache is ~1.3 TB.
- Opening `ttm-ws` from inside `projects/ece/` applies both the workspace
  `AGENTS.md` and `ttm-ws/CLAUDE.md`; the two must not contradict on
  command rules (`uv run --project gfys_workspace` versus `uv sync --all-packages`).
