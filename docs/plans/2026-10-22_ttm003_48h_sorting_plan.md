---
title: TTM-003 48 h spike sorting plan (SOM+ capture, single-unit purity, 48 h tracking)
updated: 2026-10-22
---

# TTM-003 48 h spike sorting plan

Written 2026-10-22 for the agent that continues after the workspace moves
to `~/projects/ttm-ws/`. Read `CLAUDE.md` first (conventions, git rules,
analysis practice), then `docs/SETUP.md` (environment, per-machine config,
the tononi-2 notes), then this file, then the documents it points at.
Nothing on the 48 h recording has been started; the user said "don't
start yet" on 2026-10-22, and the fresh agent starts.

## 1. Goal, in the user's words

Graham Findlay (gfindlay@wisc.edu), 2026-10-22:

> Eventually we will analyze 48h continuous recordings of TTM3, which
> will have pulses for optotagging delivered every 2h, and we will be
> especially interested in tracking the activity of the SOM+ cells, and
> seeing what they do over the 48h. So we both need to capture their
> spikes, and also have good single-unit cluster purity, and we need to
> do it over 48h.

> I do also care very much about the SOM- units for TTM3. I just care
> _more_ about the SOM+ cells.

The deliverable is a 48 h sorting of the whole population in which the
tagged SOM+ cells are captured completely, including the spikes the
current sorter drops; every unit offered as single is credibly single by
criteria the user accepts; and identity is tracked across the 48 h with
the tagging blocks as waypoints.

## 2. What the TTM-003 short sessions taught

Two optogenetic sessions of 2026-10-16 (32 and 21 min; 207 ms pulses at
about 1 Hz and 64 ms pulses at 2.86 Hz; 0.5 to 9.8 mW; optogenetic
stimulation of SOM+ cells per the user, opsin and line not recorded in
the repo). Everything is in `subjects/ttm-003/` (README, scripts,
`results/`, the stimulation log) and the findings doc
`subjects/ttm-003/optogenetic_sessions_findings.md`, status
needs_verification, interpretive sentences marked as conjecture. The
points the pipeline must answer to:

1. No spike-band light artifact at any power (findings section 1).
   Sorting ran unblanked and should stay unblanked.
2. The population response saturates at the lowest power: 40 of 56 and 44
   of 56 units silenced within 5 ms, to 1 to 5 % of baseline through the
   pulse, recovering 45 to 105 ms after offset with a rebound, the same at
   0.5 mW as at 9.8 (section 4). Every unit therefore has a per-pulse
   signature, not only the tagged ones.
3. The TT6 SOM+ candidate (unit 21 in session 1, 22 in session 2) fires
   one or two full-size spikes 3.5 to 8 ms after onset and, in the 5.5
   MAD sort, nothing more. It in fact keeps firing through the pulse at
   130 to 160 per second with every spike shrunk to half amplitude. A
   re-sort of TT6 alone at threshold 3.5 folds those spikes into the
   cell's cluster, and the large-versus-small cross-correlogram shows one
   cell, burst-attenuated (sections 5, 5a, 5b).
4. TT10 carries a second light-driven cell at about 4.3 MAD whose
   waveform differs from the reference unit's. Lower-threshold
   MountainSort5 re-sorts absorb its spikes into a multi-unit cluster
   with a silenced neighbour; one session at one threshold split it out
   cleanly, and the split did not survive other settings (section 5c).
   Lowering the detection threshold improves capture and only by chance
   improves isolation.
5. What the sort must therefore do, for every unit: detect at about
   3.5 MAD of a noise estimate that busy channels do not inflate, assign
   by waveform shape with the amplitude left free, reject by the host
   unit's refractory period, arbitrate between same-tetrode neighbours
   with cross-correlograms, and judge purity by pictures of
   within-cluster waveform spread and by the amplitude-split test, not
   by a cosine scalar.

Tool limits that bear on this (MountainSort5's int32 ceiling and 5.5
default, the time-vector memory trap, Kilosort4's GPU wall, fictional
geometry, the amplitude gate of the production matcher, stream and
tetrode naming): `docs/developer_notes/sorting_tool_limits.md`.
Literature on sorting with optogenetic pulses:
`docs/reports/2026-10-21_optotagging_literature_survey.md`.

## 3. The TTM-001 48 h pipeline and where it stopped

TTM-001 session `2026-05-27_09-07-52` (TTM-NOD, 48.3 h) is the only 48 h
sort done. The record is `investigations/spike-sorting/` (README,
`MATCHING_PURSUIT_FINDINGS.md`, `RESIDUAL_CAPTURE_PLAN.md`,
`OPEN_CONCERNS.md`); the deliverable `assembled_prod` (71 single plus 16
per-tetrode `is_mua` units) and its analyzers live under that session's
`sortings_seed42_pcafix/track_eval/mp_long_s2000_d170000/` on the server.

The chain as run: bandpass 300 to 6000 Hz plus common median reference,
materialised once (`_mp_common.materialize_span`,
`ttm.sorting.preprocess_for_sorting`); seed templates from a MountainSort5
scheme 2 sort at threshold 5.5, deduplicated (scripts 51, 52); windowed
circus-omp carry-forward over 1800 s windows with per-window template
re-estimation, re-seeding every 12 windows, bank dedup at cosine 0.95
(script 58) giving `assembled_reseed_c12`; then
`97_production_pipeline.py`: CCG-guarded merge, residual single-unit
capture (unclaimed events of 10 MAD or more, cosine 0.8, refractory
disposal), per-tetrode multi-unit bucket (7 MAD or more, cosine 0.55 to
0.8), three-axis scoring (coverage, assignment purity, identity); 3.2 h
for those stages over 95 windows, checkpointed. `98_build_analyzer_prod.py`
builds the curation analyzers.

Why it fails on the TTM-003 cells as it stands: circus-omp accepts a
spike only at 0.8 of the template amplitude or more, and the TT6 in-pulse
spikes are at 0.5; detection is at 5.5 MAD with floors of 10 and 7, and
the TT10 cell is at 4.3. The June prior that near-threshold events are
unclaimable multi-unit activity is what these cells violate.

Two things to carry from June. First, the wobble matcher with a fixed
per-spike cosine gate (shape, amplitude free) was built and wired
(`_mp_common.run_matching(shape_gate_r=)`,
`tests/test_wobble_primary_wiring.py`) and beat circus-omp on the
precision and coverage frontier, but was never adopted: the full 48 h
go/no-go was not run, and the June 19 review found those metrics do not
measure per-unit assignment purity. Script 85 showed the full-window
cosine assigns 54 % of disputed low-amplitude spikes to a same-tetrode
neighbour (72 to 80 % on crowded tetrodes), which refractory contamination
cannot see; the gate as implemented deletes rather than reassigns.
Second, `OPEN_CONCERNS.md` (user-confirmed): the user does not believe
the production clusters are clean single units (within-cluster waveform
variability too high, purity bar too loose), and MAD noise estimates are
inflated on busy channels. Both unresolved, both apply here.

Reusable, tested helpers: `_assignment_eval.py`
(`best_template_for_events`, `ccg_verdict_pair`, `ccg_guarded_merge`,
`window_bank`), `_mp_common.py` (`run_matching`, `per_spike_cosine`,
`all_template_cosines`, `competitive_reassign`), `_scoreboard.py`,
`97_production_pipeline.py` (`refractory_violation`, `stage_residual`).
Three tests under `packages/ttm/tests/` import them.

## 4. Data

- The two short sessions are converted, sorted and analysed;
  `subjects/ttm-003/README.md` lists what exists and where.
- The 48 h recording, `2026-10-20_09-02-12`, is described in that README:
  47.98 h, one recording, blank channel names that need the documented
  repair before conversion, an IMU stream, TTL events present, no
  manifest, no store. The pulse pattern, power and block times are not in
  the repo.
- Storage: the blosc-zstd store will be about 350 GB and the preprocessed
  sorting cache about 1.3 TB (TTM-001's numbers). Both fit on
  `/ssd-raid0`; only the store fits comfortably on `/nvme`
  (`docs/SETUP.md`, tononi-2). Do not delete anything under TTM-001's
  directory to make room.
- `scripts/stim_trains.py` iterates over `stimulation_log.toml` and
  raises when the hand log's train count differs from the TTL's, so the
  48 h session needs either a hand-log section or a code path for
  sessions without one. Expect the sub-millisecond glitch before each
  train and the line 2 crosstalk seen in October.

## 5. The plan

Five tracks. A and B are independent and start first. C needs B. D
needs A and C. E needs D.

### Track A. Capture and purity module, built on the 20 min sessions

A library function that, given a preprocessed recording, seed templates
and their spike trains, assigns low-amplitude and burst-attenuated spikes
to the right unit by shape, rejects by refractory period, arbitrates
between same-tetrode neighbours, and reports purity in terms the user
accepts. Needed for every unit; proven on the tagged ones.

1. Read `MATCHING_PURSUIT_FINDINGS.md` ("Per-unit normalized gate",
   "Wobble as primary", "Phase B ladder" including the script 85 result)
   and `OPEN_CONCERNS.md`, then the helpers listed in section 3. Reuse
   them.
2. Detection with a noise estimate that busy channels do not inflate
   (spikes excised or threshold crossings masked before the MAD), at
   about 3.5 MAD. Measure on TT10 (busy) and TT6 (quiet) whether
   per-channel MAD tracks firing rate before and after; this answers
   open concern 2.
3. Assignment: cosine of each event to every same-tetrode template over
   the full window (-1.0 to +2.0 ms; June showed tight windows cost real
   spikes), amplitude free, offered competitively to the best template,
   accepted if it clears the gate and does not violate the host's
   refractory period, otherwise left for a multi-unit bucket. Compare
   the wobble matcher as detector-plus-assigner against this two-stage
   version head to head on the same windows, with plots of what each
   admits.
4. Purity, per unit and as pictures: refractory contamination; the
   large-versus-small amplitude-split cross-correlogram
   (`opto_tt_amplitude_split.py`); cross-correlograms against every
   same-tetrode neighbour with the CCG arbiter's verdicts; within-cluster
   waveform overlays and amplitude-by-shape distributions. Ask the user
   what the June "cosine purity bar" meant before encoding one.
5. Validation:
   - TT6: the in-pulse half-amplitude spikes are claimed (126 to 141 per
     second in the 9.8 mW pulse; zero-lag large-versus-small CCG at a
     tenth of independence) and the silenced CH49 neighbour (re-sort
     cluster 9, -180 uV, 3.4 Hz) gains no in-pulse spikes.
   - TT10: the driven cell is separated from the silenced CH19 neighbour
     in both sessions; session 2's threshold-4.0 clusters 3 and 4 show
     what separated looks like (flat CCG outside trains, opposite pulse
     responses).
   - All 56 units: pulse responses (`opto_population.py`, `opto_tag.py`)
     before and after; a silenced unit that gains in-pulse spikes has
     absorbed a neighbour.
   - Cross-session: the sessions are 40 min apart and share cells (21 to
     22, 34 to 35 matched by waveform). Match units across them and check
     that pulse signatures agree; the smallest identity test available.
6. Deliver as a module in `packages/ttm/src/ttm/` with tests that run on
   the two sessions (marked as needing the data volume), a developer
   note, and a section in the TTM-003 findings doc. The subject scripts
   stay the drivers; no new numbered investigation scripts.

Done when the validations pass with figures the user has seen, and the
module runs on an 1800 s window of the 48 h store in minutes.

### Track B. 48 h ingestion

1. Repair the oebin as `openephys_blank_channel_names.md` describes
   (under `sg paxilline`; keep `structure.oebin.orig`), and check with
   `read_openephys`: `CH1` to `CH64`, one segment, 47.98 h, gain 0.195.
2. Settle the storage location with the user (section 7), then
   `ttm-convert` with blosc-zstd (about 2 h with 32 jobs). WavPack is
   optional.
3. Manifest (`subjects/ttm-003/experiments/`, name and fields to confirm;
   `subjects/ttm-001/experiments/` shows the 48 h schema in use) and
   `ttm-experiment-params`; then `ttm-make-lfp` with `--subsample 125`,
   `ttm-swa all`, `ttm-emg` (30 min check, then `--full`).
4. Pulse tables from TTL line 1 (gap over 2 s between trains, drop
   sub-millisecond glitches): number of blocks, pulses per block, width,
   rate, block times in session seconds and local time. Add the hand log
   if the user has one.
5. Per-block checks with the existing scripts generalised to many
   blocks: artifact check (`opto_onset_zoom.py`), suppression per block
   (`opto_explore.py`, `opto_population.py`). This shows before any sort
   whether the signature is present at every block and whether it
   changes with time of day or state.
6. Channel QC: dead or saturated wires, the 400 Hz comb
   (`investigations/acquisition-noise/`), whether CH2 and TT14 are still
   bad. Write the lists into the manifest.

Done when the store, LFP, SWA, EMG, params, pulse tables and per-block
figures exist and the subject README documents the 48 h experiment the
way it documents the short ones.

### Track C. Baseline: the TTM-001 production pipeline on the 48 h

Run the June pipeline unchanged as the baseline for the untagged
population and as the thing the candidate must beat.

1. Materialise the preprocessed binary for the whole 48 h where the store
   is; it is shared by tracks C and D.
2. Seed bank, carry-forward with reseed, `97_production_pipeline.py`,
   `98_build_analyzer_prod.py`. The scripts carry literal TTM-001 paths;
   parametrise by `Experiment` rather than editing paths in place, and
   keep the June parameters so the baseline is the documented one.
3. Score with the three-axis scoreboard and with the tagged-cell tests
   from track A at every pulse block: does the unit matched to the TT6
   cell show the onset spike and the in-pulse train, and does the TT10
   cell exist at all. Expected and worth recording: the baseline keeps
   the onset spikes and drops the trains.
4. Record wall-clock, peak memory and disk per stage.

Done when a baseline `assembled_prod` for the 48 h exists with analyzers,
a scoreboard, and the tagged-cell checks tabulated per block.

### Track D. Candidate pipeline and the comparison

The baseline framework with track A's module in place of the
amplitude-gated matcher and the 5.5/10/7 MAD floors:

1. Track A's detection and noise estimate in every window.
2. Matching with the amplitude left free (whichever of wobble-with-gate
   or detect-and-assign won in track A), competitive among same-tetrode
   templates, refractory disposal, inside the windowed carry-forward
   with template re-estimation and re-seeding.
3. The CCG-guarded merge, with the user's purity criteria from track A
   applied to every merge and every unit and the waveform spread
   reported, so open concern 1 is answered by evidence.
4. Identity across the 48 h with the pulse blocks as waypoints every
   2 h. For tagged cells the tagging signature must hold at every block
   the unit is present. For other units, test whether the suppression
   profile (silence latency, in-pulse fraction, recovery, rebound) is
   stable across blocks; conjecture: it can serve as a weaker identity
   check. Measure first on the two short sessions (track A, step 5).
5. Compare baseline and candidate on the same 48 h: scoreboard;
   tagged-cell capture and continuity per block; purity pictures for a
   random sample of units and for every unit the two disagree on;
   identity continuity (presence fraction, swaps caught by waypoints).
   Present as a decision dashboard, not pass/fail.

Done when the user has chosen a pipeline on the evidence, recorded with
the numbers under `investigations/` or `docs/reports/`.

### Track E. Deliverables

The chosen sorting goes through the existing hand-offs: curation analyzer
(group sparsity, `is_mua`, `unitrefine_*` columns), `ttm-curation-manifest`
and `ttm-curate`, `ttm-viz-bundle` and `ttm-loupe`. Register it where the
June work's `mp_production_sorting` entry lives (there is no
`docs/artifacts/registry.yaml` in this workspace yet; create one or point
at the registry's new home). Then the science, which is the user's
question and not this plan's.

### Parallelism

A and B can run as two agents at once; they touch different files
(`packages/ttm` and `subjects/ttm-003/results` for A; the data volume,
the manifest and `stim_trains.py` for B). C is one long job after B. D
after A and C. Several agents at once: small commits straight to `main`,
`git -C` with absolute paths, one agent at a time owning write-heavy jobs
on the data volume.

## 6. Context specific to this hand-off

- In the October TTM-003 sessions the user let the agent commit as it
  went, against the `CLAUDE.md` default of leaving edits unstaged. Ask
  once which the user wants now.
- The collaborator `txiao37`'s standalone clone and the TTM-002 products
  on the server are not adopted; TTM-002 is re-analysed with this
  pipeline when its turn comes.
- `ttm-ws` main was pushed on 2026-10-22 after this plan; earlier October
  commits were local until then.

## 7. Questions for the user before starting

1. Storage: `data_root` stays on `/nvme` with the 48 h session directory
   symlinked to `/ssd-raid0`, or a second root in the config? The symlink
   needs no code.
2. The 48 h experiment's manifest name and fields (light/dark schedule,
   any deprivation, state scoring plans).
3. A hand log for the 48 h tagging blocks (pattern, power, times), if one
   exists.
4. What the June "cosine purity bar" meant, so track A encodes the right
   criterion (recorded as unresolved in `OPEN_CONCERNS.md`).
5. Whether the agent commits as it goes or leaves edits unstaged.
