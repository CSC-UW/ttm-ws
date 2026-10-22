---
title: TTM-003 48 h spike sorting plan (SOM+ capture, single-unit purity, 48 h tracking)
updated: 2026-10-22
---

# TTM-003 48 h spike sorting plan

Written 2026-10-22 for the agent that picks this work up after the
workspace moves from `~/projects/ece/ttm-ws/` to `~/projects/ttm-ws/`.
That agent will not have the `gfys_workspace` rules or the previous
agent's memory, so this document carries everything it needs: the goal,
what the TTM-003 short sessions taught, what the TTM-001 pipeline is and
where it stopped, the data and machine facts, the plan, and the working
rules. Read `CLAUDE.md` first, then this file, then the two findings
documents it points at.

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

So the deliverable is a 48 h sorting of the whole population in which
the tagged SOM+ cells are captured completely (including the spikes the
current sorter drops, section 2), every unit offered as single is
credibly single by criteria the user accepts, and identity is tracked
across the 48 h with the tagging blocks as waypoints. The 48 h recording
exists (section 4.3). Nothing has been started on it. The user said
"don't start yet" on 2026-10-22; the fresh agent starts.

## 2. What the TTM-003 short sessions taught (2026-10-21 and 22)

Two optogenetic sessions on 2026-10-16, each 20 to 32 min, with laser
trains of 207 ms pulses at about 1 Hz (pattern A) and 64 ms pulses at
2.86 Hz (pattern B), at 0.5 to 9.8 mW, optogenetic stimulation of SOM+
cells (per the user; the opsin and the animal line are not recorded in
the repo; see `subjects/ttm-003/optogenetic_stimulation.md` and
`stimulation_log.toml`). Full write-up with numbers:
`subjects/ttm-003/optogenetic_sessions_findings.md` (status
needs_verification; interpretive sentences there are marked as
conjecture). The points that bear on the pipeline:

1. There is no spike-band light artifact. Edge-locked high-pass averages
   at pulse onset and offset are 3 to 10 uV against 2 uV noise at both
   0.5 and 9.8 mW. There is a slow evoked potential of hundreds of uV
   that is weakly power dependent (conjecture: physiological). Sorting
   ran unblanked and should stay unblanked.
2. The population response saturates at the lowest power used: 40 of 56
   units (session 1) and 44 of 56 (session 2) are silenced within 5 ms
   of onset, to 1 to 5 % of baseline through the pulse, recover 45 to
   105 ms after offset and rebound 1.2 to 1.9 times, the same at 0.5 mW
   as at 9.8 mW. Suppression is spatial (TT4 to TT13 yes, TT1 to 3, 15,
   16 no). Every unit therefore has a per-pulse signature, not only the
   tagged ones.
3. The TT6 SOM+ candidate (session 1 unit 21, session 2 unit 22; 13 Hz,
   -110 uV on CH49, SNR 9.5) fires one or two full-size spikes 3.5 to
   8 ms after onset and then, in the 5.5 MAD sort, nothing. The user
   suspected the sorter was dropping the rest of a train, and was right:
   the cell keeps firing through the pulse at 130 to 160 per second with
   every spike shrunk to about -55 uV (half amplitude). A re-sort of TT6
   alone at threshold 3.5 folds those spikes into the cell's cluster
   (13 Hz becomes 36 and 46 Hz; 141 and 126 per second through the 9.8 mW
   pulse, 340 per second in the first bins, silent for 100 ms after
   offset). The large and small spikes avoid each other at zero lag in
   their cross-correlogram (35 observed against 316 expected in session
   1, 31 against 303 in session 2) and the small ones follow the large
   ones at 2 to 3 ms: one cell, burst-attenuated. Scripts
   `opto_tt6_inpulse.py`, `opto_tt_resort.py`, `opto_tt_amplitude_split.py`;
   findings doc sections 5, 5a, 5b.
4. TT10 carries a second light-driven cell that the 5.5 sort never
   isolates. Its spikes are about -50 uV on a 12 uV noise MAD (4.3 MAD),
   broader than the reference unit's (trough to peak 0.83 vs 0.63 ms).
   MountainSort5 re-sorts at 3.0 to 4.0 give one multi-unit cluster per
   wire (refractory contamination 0.2 to 0.5); the CH19 cluster absorbs
   the driven spikes but mixes them with a silenced neighbour. Only
   session 2 at threshold 4.0 split out a clean driven unit (contamination
   0.05; flat cross-correlogram against the neighbour outside the trains),
   and that split did not survive other thresholds or a longer snippet
   window. Findings doc section 5c. Lesson: lowering the detection
   threshold improves capture, and only by chance improves isolation.
5. "Hash" in these notes and in the TTM-001 documents means multi-unit
   activity: real spikes of several unresolved neurons, not noise. The
   June documents sometimes write "hash/noise"; keep the two apart.
6. Display and detection thresholds hide the problem. The galleries' grey
   unassigned ticks were at 4.5 MAD and showed nothing for the TT10
   spikes; `--thresh-mad 3.5` shows them.
7. The 5.5 MAD threshold is MountainSort5's default and the SpikeInterface
   wrapper's; MountainSort5 applies it to the whitened signal. The lab's
   `tetrode` property gives names TT1 to TT16 (1-based); the SpikeInterface
   `group` property is 0-based. Early figures were mislabeled by one
   tetrode and regenerated; use `TT{group + 1}` everywhere.
8. Bad channels in both sessions: CH2 (TT13) and CH10, CH12, CH14, CH16
   (all of TT14). They are in the manifests as `badChannels` and
   `badTetrodes`. Check the 48 h recording separately; the same wires are
   likely but not certain.
9. A SpikeInterface bug: with a time vector set in memory, `to_dict`
   carries `times_kwargs` (the whole vector) at every nesting level and
   the Zarr writer dumped it into `.zattrs` (4.9 GB). Fixed in
   `ttm.spikeinterface.convert_recording` (commit 9b138c4, test
   `packages/ttm/tests/test_convert_provenance.py`); upstream #4885 fixes
   only the top level.

Literature on preprocessing and sorting with optogenetic pulses (SALT,
latency and jitter, waveform correlation in versus out of pulse,
photoelectric artifact signatures, blanking practices):
`docs/reports/2026-10-21_optotagging_literature_survey.md`
(agent-compiled, needs_verification).

## 3. What the TTM-001 48 h pipeline is and where it stopped

TTM-001 session `2026-05-27_09-07-52` (experiment TTM-NOD, 48.3 h) was
sorted in May and June 2026. The record is `investigations/spike-sorting/`
(README, `SORTING_COMPARISON_FINDINGS.md`, `TRACKING_FINDINGS.md`,
`MATCHING_PURSUIT_FINDINGS.md`, `RESIDUAL_CAPTURE_PLAN.md`,
`OPEN_CONCERNS.md`). Paths in those scripts are literal `tononi-2` paths
as run. The production deliverable is registered as
`ttm.mp_production_sorting`: `assembled_prod` under
`/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52/sortings_seed42_pcafix/track_eval/mp_long_s2000_d170000/`
(87 units: 71 single plus 16 per-tetrode `is_mua` pseudo-units), with
curation analyzers `analyzer_prod*.zarr` beside it.

The chain, as run:

1. Preprocess: bandpass 300 to 6000 Hz plus common median reference,
   materialised once as a binary (`_mp_common.materialize_span`,
   `ttm.sorting.preprocess_for_sorting`). The TTM-003 opto sorts used
   `cmr="global"`; the library also supports the per-group
   "out_of_group" reference from SpikeInterface PR #4601.
2. Seed templates from a MountainSort5 scheme 2 sort at threshold 5.5
   on a drift-stable span, deduplicated (scripts 51 and 52).
3. Windowed carry-forward template matching with circus-omp over 1800 s
   windows, templates re-estimated per window, re-seeding every 12
   windows, seed-bank dedup at cosine 0.95 (`58_match_pursuit_reseed.py`,
   `_mp_common.windowed_carry_forward_reseed`). Output
   `assembled_reseed_c12`.
4. `97_production_pipeline.py`: CCG-guarded within-tetrode merge (cosine
   proposes at 0.90, CCG duplicate verdict accepts, cosine 0.95 fallback),
   then residual single-unit capture (unclaimed events of 10 MAD or more
   proposed to their best same-tetrode template at cosine 0.8 or more,
   folded in unless they break the host's 1.5 ms refractory period), then
   a per-tetrode multi-unit bucket (unclaimed events of 7 MAD or more with
   best cosine in [0.55, 0.8)), then scoring on three axes (coverage,
   per-unit assignment purity, identity stability). 3.2 h for these
   stages over 95 windows. Checkpointed; stages skip if their output
   exists.
5. `98_build_analyzer_prod.py` builds the curation analyzers;
   `ttm-curate` opens them.

Why it would fail on the TTM-003 cells as it stands:

- circus-omp accepts a spike only at 0.8 of the template amplitude or
  more (`amplitudes=[0.8, inf]`). The TT6 in-pulse spikes are at 0.5.
- Event detection for coverage, residual capture and the multi-unit bucket
  is the locally exclusive peak detector at 5.5 MAD, with floors of 10 and
  7 MAD. The TT10 driven cell is at 4.3 MAD and most of its spikes are
  under 5.5.
- The June work's stated prior was that near-threshold (5 to 7 MAD)
  events are multi-unit activity no isolation sorter claims, and that
  coverage above 12 MAD was "the meaningful number". That prior is what
  the TTM-003 cells violate.

Two further things from June the new work must carry:

- The wobble matcher with a fixed per-spike cosine gate (shape, amplitude
  free) was built and wired (`_mp_common.run_matching(shape_gate_r=)`,
  `tests/test_wobble_primary_wiring.py`) and beat circus-omp on the
  precision and coverage frontier (zero median refractory contamination at
  r >= 0.57 in every window, 92 to 97 % coverage above 10 MAD, 45 to 69 %
  retention of 5.5 to 10 MAD spikes). It was never adopted: the full 48 h
  go/no-go (Phase B, C, D) was not run, and the June 19 review found the
  metrics it won on do not measure per-unit assignment purity. Script 85
  showed the full-window cosine assigns 54 % of disputed low-amplitude
  spikes to a same-tetrode neighbour (72 to 80 % on crowded tetrodes), and
  refractory contamination cannot see that. The gate as implemented
  deletes, it does not reassign. Findings doc sections "Wobble as
  primary" and "Phase B ladder".
- `OPEN_CONCERNS.md` (user-confirmed, 2026-06-20): (1) the user does not
  believe the production clusters are clean single units, within-cluster
  waveform variability is too high, and the purity bar (a cosine
  threshold) is too loose; (2) MAD-based noise estimates are inflated on
  busy channels, so detection thresholds there are too high. Both are
  unresolved and both apply directly here.

Other recorded limits of the tools:

- MountainSort5 stores spike times as int32 in `detect_spikes`, so any
  single block or training recording over 2^31 samples (19.9 h at 30 kHz)
  crashes. This, not RAM, caps block and training duration.
- `BinaryFolderRecording` loads its full float64 time vector on every
  reconstruction (no mmap), so `run_sorter_by_property` on a frame slice
  of a long recording blows up memory. Crop to a genuine recording and
  `reset_times()` before sorting; the sorting drivers do this.
- Kilosort4's final clustering needs about 42 GiB for one tetrode's 48 h
  of spikes and does not fit the 32 GB V100; a fork pin with
  `clustering_chunk_size` exists. It is not a route for 48 h as is.
- MountainSort5 is bit-reproducible only with both the whitening seed and
  the PCA solver controlled (the fork pin, now upstream in 0.5.9).
- DARTsort and other geometric drift models were rejected: tetrode
  geometry is fictional (the ProbeGroup places tetrodes 300 um apart with
  10 um wires). Drift must be handled by template re-estimation, not
  localisation.
- UnitRefine's stock noise/neural classifier is near chance on tetrodes;
  the retrained advisory ranking (`ttm.unitrefine_advisory`) is usable
  for ordering, not thresholding.

## 4. Data and machine facts

### 4.1 Machine

tononi-2 (Linux, 1.5 TiB RAM, Tesla V100 32 GB). Per-machine config
`~/.config/ttm/config.toml`:

```toml
data_root = "/nvme/neuropixels/tetrode_data"

[remotes]
archive = "/Volumes/tetrode"
```

Volumes on 2026-10-22:

| Volume | Size | Free | Notes |
|---|---|---|---|
| `/nvme` | 5.8 TB | 2.1 TB | `data_root`. TTM-001's session dir alone is 2.2 TB (store plus 1.3 TB sorting cache). Shared, often near full. |
| `/ssd-raid0` | 7.0 TB | 6.0 TB | ext4, "fair game" per the user: slower than `/nvme`, much faster than spinning disk. Has `analysis/` and `neuropixels/`. |
| `/Volumes/tetrode` | CIFS | | Raw acquisitions, `<SUBJECT>/<session>/Record Node */`. Writing needs the `paxilline` group (`sg paxilline -c '...'` in a shell that predates the membership) and `rsync --inplace`. |

The machine is shared. Never kill processes by pattern; `pkill -f`
matches the agent's own shell (it has happened twice). Find PIDs with
`pgrep -f "pat[t]ern"`, walk the parent chain, kill by PID in a separate
call. Use absolute paths and `git -C` everywhere; the shell's working
directory drifts between calls.

### 4.2 TTM-003 short sessions (done)

| Experiment | Session | Store | Sortings | Figures |
|---|---|---|---|---|
| TTM-003-OPTO1 | `2026-10-16_14-06-44` (32 min, 7 trains) | `<data_root>/<session>/<session>.blosc-zstd.zarr` | `<session>/sortings_opto/asis/{preprocessed,by_group,sorting,analyzer,resort_TT6_thr{4,3.5},resort_TT10_thr{4,3.5,3,4_t2_45,3.5_t2_45}}` | `<data_root>/TTM-003/TTM-003-OPTO1/figures/` incl. `galleries/` |
| TTM-003-OPTO2 | `2026-10-16_14-44-22` (about 21 min, 5 trains, the 9.8 mW pattern B ones among them) | same layout | same layout | `<data_root>/TTM-003/TTM-003-OPTO2/figures/` |

Manifests: `subjects/ttm-003/experiments/TTM-003-OPTO{1,2}.toml`;
`experiment_params.json` written on the data volume. Small results in
`subjects/ttm-003/results/` (stimulation train and pulse tables, unit
metrics, optotag metrics per unit per train set, population PSTH
summaries, re-sort tables, in-pulse fingerprint JSONs). Scripts and the
order to run them: `subjects/ttm-003/README.md`. The scripts are
subject-specific and take `--subject TTM-003 --experiment <E>`.

### 4.3 TTM-003 48 h recording (not started)

- Raw: `/Volumes/tetrode/TTM3/2026-10-20_09-02-12/Record Node 102/experiment1/recording1/`,
  one experiment, one recording, 671 GB. `continuous.dat` is
  663,292,502,656 bytes = 47.98 h at 64 channels, 30 kHz, int16. Software
  Time 1792504932633 ms = 2026-10-20 14:02:12 UTC = 09:02:12 Chicago.
  Streams: `acquisition_board` (64 ch, gain 0.195 uV), `memory_usage`
  (1 ch, 100 Hz), `IMU_port_A` (18 ch, 100 Hz). Events: `acquisition_board/TTL`
  and `MessageCenter`. The recording finished the morning of 2026-10-22.
- The `structure.oebin` has blank channel names (all 64 are `" "`), the
  Acquisition Board plugin bug described in
  `docs/developer_notes/openephys_blank_channel_names.md`. neo refuses the
  stream until it is repaired. The note gives the deterministic repair
  (replace the N blank names by `CH1` to `CHN` in file order, keep
  `structure.oebin.orig`), done for eleven earlier recordings on
  2026-10-21 but not for this one. The repair script lived only in a
  scratchpad; rewrite it from the note, run it under `sg paxilline`, and
  check with `spikeinterface.extractors.read_openephys` before converting.
  The rig itself still needs the dropdown fixed (note, "Prevention").
- Streams are selected by name in `ttm.spikeinterface.get_recording`
  (`oe_stream`), not by index, because the IMU stream changes the order.
- Optotagging pulses every 2 h (user). The pulse pattern, power and
  number of blocks are not logged in the repo; derive the trains from
  TTL line 1 as `stim_trains.py` does. That script iterates over
  `stimulation_log.toml` sessions and raises if the hand log's train
  count differs from the TTL's, so either add a hand-log section for
  this session (ask the user for the log) or let the script run without
  one for sessions that have no entry. The sub-millisecond TTL glitch
  toggles before each train and the line 2 crosstalk seen in October are
  likely here too.
- Where to put it: the blosc-zstd store will be about 350 GB (TTM-001's
  was 352 GB for 48.3 h) and the preprocessed sorting cache about 1.3 TB.
  Both fit on `/ssd-raid0`; only the store fits comfortably on `/nvme`.
  Decision for the user (section 8): keep `data_root` on `/nvme` and
  symlink this session's directory to `/ssd-raid0`, or add a second root.
  Do not delete anything under TTM-001's directory to make room; the
  user has not asked for that.
- A manifest is needed: `subjects/ttm-003/experiments/TTM-003-48H.toml`
  (name to confirm with the user) with `session`, `acquisition =
  "TTM3/2026-10-20_09-02-12"`, bad channel lists after a QC pass, and
  whatever light/dark or deprivation fields the 48 h experiments carry
  (see `subjects/ttm-001/experiments/` for the schema in use and
  `ttm.subject.Experiment.manifest_fields()`).

## 5. The plan

Five tracks. A and B are independent and start first. C needs B. D
needs A and C. E needs D. Each track names its outputs and the test that
says it is done.

### Track A. Capture and purity module, built on the 20 min sessions

Purpose: a library function that, given a preprocessed recording, a set
of seed templates and their spike trains, assigns low-amplitude and
burst-attenuated spikes to the right unit by shape with the amplitude
left free, rejects by refractory period, arbitrates between same-tetrode
neighbours, and reports purity in terms the user accepts. It is needed
for every unit, not only the tagged ones; the tagged cells are where it
can be proven.

Steps:

1. Read `MATCHING_PURSUIT_FINDINGS.md` sections "Per-unit normalized
   gate", "Wobble as primary", "Phase B ladder" (including the script 85
   reassignment result) and `OPEN_CONCERNS.md`. Read
   `_assignment_eval.py` (`best_template_for_events`, `ccg_verdict_pair`,
   `ccg_guarded_merge`, `window_bank`), `_mp_common.py`
   (`run_matching`, `per_spike_cosine`, `all_template_cosines`,
   `competitive_reassign`) and `97_production_pipeline.py`
   (`refractory_violation`, `stage_residual`). These are tested helpers;
   reuse them rather than rewriting. Three tests in `packages/ttm/tests/`
   import them.
2. Build a detection stage that does not depend on an inflated MAD:
   estimate noise with spikes excised (or from the median absolute
   deviation of a trace with threshold crossings masked), and detect at
   about 3.5 MAD of that. Measure on TT10 (busy) and TT6 (quiet) whether
   the per-channel MAD correlates with firing rate before and after
   excision; this answers open concern 2.
3. Assignment: for each detected event, cosine to every same-tetrode
   template over the full window (-1.0 to +2.0 ms; the June tests showed
   tight windows cost real low-amplitude spikes), amplitude free; offer it
   to the best template (competitive, not delete-only); accept if cosine
   clears the gate and the event does not violate the host's refractory
   period; otherwise leave it for a multi-unit bucket. Compare the wobble
   matcher as the detector-plus-assigner against this two-stage version
   head to head on the same windows, with plots of what each admits (the
   user's rule: compare, don't assert; show pictures).
4. Purity criteria, each reported per unit and shown as pictures, not
   only as scalars: refractory contamination from the autocorrelogram;
   the large-versus-small amplitude split cross-correlogram
   (`opto_tt_amplitude_split.py` logic: same cell if the two halves share
   a refractory period); cross-correlograms against every same-tetrode
   neighbour with the CCG arbiter's verdicts; within-cluster waveform
   spread (overlays and amplitude-by-shape distributions, the thing the
   user said the cosine scalar hides). Ask the user what their "cosine
   purity" bar from June meant before encoding one; it is recorded as
   unresolved in `OPEN_CONCERNS.md`.
5. Validation on the opto sessions:
   - TT6 cell: the module must claim the in-pulse half-amplitude spikes
     (target: in-pulse rate within the 126 to 141 per second the 3.5
     re-sort found, zero-lag large-versus-small CCG at a tenth of
     independence) without admitting the silenced CH49 neighbour
     (re-sort cluster 9, -180 uV, 3.4 Hz) during pulses.
   - TT10: the module must separate the driven cell from the silenced
     CH19 neighbour in both sessions, which MountainSort5 managed once by
     luck; the session 2 threshold-4.0 clusters 3 and 4 are the reference
     for what "separated" looks like (flat CCG outside trains, opposite
     pulse responses).
   - All 56 units: each unit's pulse response (silence latency, in-pulse
     fraction, recovery, rebound; `opto_population.py` and `opto_tag.py`)
     before and after the module. A silenced unit that gains in-pulse
     spikes has absorbed a neighbour.
   - Cross-session: the two sessions are 40 min apart and share cells
     (unit 21 to 22, 34 to 35 matched by waveform). Match units across
     them by template and check that the pulse signatures agree. This is
     the smallest identity test available before the 48 h.
6. Deliver as `packages/ttm/src/ttm/<module>.py` with tests that run on
   the two opto sessions (mark them as requiring the data volume), a
   developer note in `docs/developer_notes/`, and a short findings
   section appended to the TTM-003 findings doc. Keep the subject scripts
   as the drivers; do not add more numbered investigation scripts.

Done when: the three validations pass with figures the user has looked
at, and the module runs on a 1800 s window of the 48 h store in minutes.

### Track B. 48 h ingestion

1. Repair the oebin (section 4.3). Check with `read_openephys` that the
   stream reads as `CH1` to `CH64`, one segment, 47.98 h, gain 0.195.
2. Settle the storage location with the user (section 8), then
   `ttm-convert` with blosc-zstd (about 2 h with 32 jobs; the WavPack
   store is optional and was only ever a sorting-validation variant).
   The converter writes `<session>/<session>.blosc-zstd.zarr`,
   `slice_table.csv`, `conversion_results.json` and the recording's
   provenance in the store attrs (without the time vector).
3. Manifest and `ttm-experiment-params`. Then `ttm-make-lfp`
   (`--subsample 125` too), `ttm-swa all`, `ttm-emg` (30 min check first,
   then `--full`). These are the standard steps and give the sleep and
   deprivation context the 48 h analysis will need.
4. Pulse tables: extend `stim_trains.py` to this session (TTL line 1,
   gap > 2 s between trains, drop sub-millisecond glitches). Report the
   number of blocks, pulses per block, pulse width and rate, and block
   times in session seconds and local time. If the user has a hand log
   for the 48 h session, add it to `stimulation_log.toml`.
5. Per-block checks with the existing scripts, generalised to many
   blocks: artifact check (`opto_onset_zoom.py`), population suppression
   per block (`opto_explore.py` MUA PSTH and `opto_population.py`).
   This establishes before any sorting whether the suppression signature
   is present at every block over 48 h, and whether anything changes
   with time of day or state.
6. Channel QC on the 48 h: dead or saturated wires, the 400 Hz comb noted
   in `investigations/acquisition-noise/`, and whether CH2 and TT14 are
   still bad. Write the lists into the manifest.

Done when: the store, LFP, SWA, EMG, params, pulse tables and per-block
figures exist and `subjects/ttm-003/README.md` documents the 48 h
experiment the way it documents the two short ones.

### Track C. Baseline: the TTM-001 production pipeline on the 48 h

Run the June pipeline as it is, as the baseline for the untagged
population and as the thing the candidate must beat. Do not improve it
in this track.

1. Materialise the preprocessed binary for the whole 48 h
   (`materialize_span` equivalent; the opto sessions' `opto_sort.py`
   shows the library call). This cache is the big disk item and is shared
   by tracks C and D; put it where the store is.
2. Seed bank, carry-forward with reseed, then `97_production_pipeline.py`
   and `98_build_analyzer_prod.py`. The scripts have literal TTM-001
   paths; parametrise them by `Experiment` rather than editing paths in
   place, and put the driver under `subjects/ttm-003/scripts/` or a
   `ttm-sort` entry point if it becomes general. Keep the June parameters
   (threshold 5.5, circus-omp amplitudes 0.8, residual floor 10 MAD,
   multi-unit floor 7 MAD, 1800 s windows, reseed every 12) so the
   baseline is the documented one.
3. Score it with the three-axis scoreboard (`_scoreboard.py`), and
   additionally with the tagged-cell tests from track A: at every pulse
   block, does the unit matched to the TT6 cell show the onset spike and
   the in-pulse train, and does the TT10 driven cell exist at all?
   Expected, and worth recording: the baseline captures the onset spikes
   and drops the in-pulse trains, and may not have the TT10 cell.
4. Record wall-clock, peak memory and disk for each stage in the findings
   doc. The June run took 3.2 h for the merge, residual and multi-unit
   stages; the carry-forward time is in `MATCHING_PURSUIT_FINDINGS.md`.

Done when: a baseline `assembled_prod` for the 48 h exists with analyzers
and a scoreboard, and the tagged-cell checks are tabulated per block.

### Track D. Candidate pipeline and the comparison

The candidate is the baseline framework with track A's module in place of
the amplitude-gated matcher and the 5.5/10/7 MAD floors:

1. Detection and noise estimate from track A step 2, in every window.
2. Matching with the amplitude left free (wobble with the cosine gate, or
   the two-stage detect-and-assign, whichever won in track A), with
   competitive assignment among same-tetrode templates and refractory
   disposal, inside the windowed carry-forward with template
   re-estimation and re-seeding.
3. The CCG-guarded merge as in June, but with the user's purity criteria
   from track A applied to every merge and every unit, and the
   within-cluster waveform spread reported, so open concern 1 is answered
   by evidence rather than by a threshold.
4. Identity across the 48 h: the pulse blocks are waypoints every 2 h.
   For the tagged cells, the tagging signature must hold at every block
   a unit is present. For every other unit, test whether its suppression
   profile (silence latency, in-pulse fraction, recovery time, rebound)
   is stable across blocks; conjecture: it can serve as a weaker identity
   check for untagged units. Measure on the two opto sessions first
   (track A step 5) whether matched units have stable signatures 40 min
   apart.
5. Compare baseline and candidate on the same 48 h: three-axis scoreboard;
   tagged-cell capture and continuity per block; purity pictures for a
   random sample of units and for every unit the two pipelines disagree
   on; identity continuity (presence fraction, swaps caught by the
   waypoints). Present as a decision dashboard, not a pass/fail.

Done when: the user has chosen a pipeline on the evidence, and the
choice and the numbers are in a findings doc under
`investigations/` or `docs/reports/`.

### Track E. Deliverables

The chosen sorting goes through the existing hand-offs: curation
analyzer (`98_build_analyzer_prod.py` equivalent, group sparsity,
`is_mua`, `unitrefine_*` advisory columns), `ttm-curation-manifest` and
`ttm-curate` for the user's curation pass, `ttm-viz-bundle` and
`ttm-loupe` for viewing with EMG and LFP. Register the sorting in the
registry the June work referenced (`docs/artifacts/registry.yaml` does
not yet exist in this workspace; create it or point at where the registry
lives now). Then the science: what the SOM+ cells do over 48 h, which is
the user's question and not this plan's.

### Parallelism

A and B can run as two agents at once; they touch different files
(`packages/ttm` plus `subjects/ttm-003/results` for A; the data volume,
the manifest and `stim_trains.py` for B). C is one long job after B. D
after A and C. If several agents work the repo at once: small commits
straight to `main`, no feature branches, `git -C` with absolute paths,
and one agent at a time owns write-heavy jobs on the data volume.

## 6. Working rules carried over

From `CLAUDE.md` (authoritative) plus the user's standing feedback that
lived outside this repo:

- Never push, never open a pull request, unless the user says so in
  that conversation. The user says "Push" explicitly.
- Commit to `main` directly in small, self-contained commits. `CLAUDE.md`
  says to leave edits unstaged for the user; in the October TTM-003
  sessions the user let the agent commit. Ask once.
- Commit messages: review them as a hostile reader would before
  proposing (tense, causation words, vagueness, "every" where a partial
  case exists); body lines under 72 characters.
- Destructive actions (deleting data, killing processes, overwriting
  stores) only with explicit confirmation, and look at the target first.
- No scientific or interpretive claim in tracked docs without the user's
  confirmation; mark unconfirmed text as conjecture. Operational facts
  need no confirmation.
- Prose in the repo must read as human-written: no bold-label sentence
  openers, no rhetorical italics, no em-dashes, no "[CONFIRM]" blocks.
- Compare, don't assert: when two methods compete, run them head to head
  on the data and show what each admits or produces, not only summary
  metrics. Before calling something an open problem, check whether an
  earlier script already prototyped the answer.
- Validate approximations: a windowed preview, a subsample or a cached
  proxy gets an explicit agreement measurement against the exact result
  before it is adopted.
- Do not call a window inside a stimulus "baseline" or "quiet"; name it
  by its position relative to onset and offset and say what it may
  contain.
- Flag provisional data sources and build around them so they can be
  swapped.
- Compute added to existing analyses gets an opt-out toggle and
  iteration caps.
- Use only the exact recording a sorting was produced from for anything
  tied to that sorting.
- Measure before diagnosing slowness or claiming a speed-up.
- Bulk products never enter git; committed results are small summaries
  next to the code that made them. No file over 5 MB.
- The collaborator `txiao37`'s standalone clone and TTM-002 products are
  not adopted; TTM-002 is re-analysed with this pipeline when its turn
  comes.

## 7. Pointers

| What | Where |
|---|---|
| Workspace rules | `CLAUDE.md`, `docs/SETUP.md`, `docs/pipeline.md`, `packages/ttm/README.md` |
| TTM-003 findings and logs | `subjects/ttm-003/optogenetic_sessions_findings.md`, `optogenetic_stimulation.md`, `stimulation_log.toml`, `README.md`, `results/`, `scripts/` |
| Optotagging literature | `docs/reports/2026-10-21_optotagging_literature_survey.md` |
| TTM-001 sorting record | `investigations/spike-sorting/README.md`, `MATCHING_PURSUIT_FINDINGS.md`, `RESIDUAL_CAPTURE_PLAN.md`, `OPEN_CONCERNS.md`, `SORTING_COMPARISON_FINDINGS.md`, `TRACKING_FINDINGS.md`, `KS4_FEASIBILITY_FINDINGS.md`, `ms5_parameter_comparison.md` |
| Reusable helpers | `investigations/spike-sorting/_mp_common.py`, `_assignment_eval.py`, `_scoreboard.py`, `_wobble_eval.py`, `_track_eval.py`; `packages/ttm/src/ttm/{sorting,tracking,spikeinterface,unitrefine_advisory}.py` |
| Blank channel names | `docs/developer_notes/openephys_blank_channel_names.md` |
| Analyzer group sparsity | `docs/developer_notes/analyzer_group_sparsity.md` |
| Frame-slice memory | `investigations/si-frame-slice-memory/` |
| Acquisition noise | `investigations/acquisition-noise/` |
| Previous restructure plan | `docs/plans/2026-10-16_ttm-ws-restructure.md` |

## 8. Questions for the user before starting

1. Storage: `data_root` stays on `/nvme` with the 48 h session directory
   symlinked to `/ssd-raid0`, or a second root in the config? The symlink
   needs no code.
2. The 48 h experiment's manifest name and fields (light/dark schedule,
   any deprivation, state scoring plans).
3. A hand log for the 48 h tagging blocks (pattern, power, times), if one
   exists.
4. What the June "cosine purity bar" meant, so track A encodes the right
   criterion (recorded as unresolved in `OPEN_CONCERNS.md`).
5. Whether the agent commits as it goes, as in October, or leaves edits
   unstaged as `CLAUDE.md` says.
