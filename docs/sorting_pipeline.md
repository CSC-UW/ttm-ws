---
title: Spike sorting pipeline for 48 h tetrode recordings
updated: 2026-11-03
status: needs_verification
---

# Spike sorting pipeline

What is done to a 48 h tetrode recording to get spike trains, in order,
with the parameters used. Each step has a collapsed "Why" note, and each
note ends with links to the studies behind the choice and the
alternatives they tested. The pipeline is the one that produced
`carry_A48_pt/assembled_reseed_c12` on TTM-003-48H (2026-11-03).

Where it stands:

- Steps 1 and 2 are library code (`ttm-bad-periods`,
  `ttm.sorting.preprocess_for_sorting`). Steps 3 to 7 still run from
  `investigations/carry-bakeoff/` (`run_aligned48.sh`) and are being
  moved into the library as one `ttm-*` command.
- It has run end to end on TTM-003-48H only.
- Steps 8 and 9 (analyzer, metrics, curation) are not yet defined for
  these sortings.

Data volumes and run times below are for TTM-003-48H (64 channels, 30 kHz,
48 h) on tononi-2.

| # | Step | Output |
|---|---|---|
| 0 | [Before sorting](#0-before-sorting) | Zarr store, `experiment_params.json` |
| 1 | [Bad periods](#1-bad-periods) | `bad_periods/bad_periods.csv` |
| 2 | [Filter and reference](#2-filter-and-reference) | preprocessed binary (1.3 TB) |
| 3 | [Whitening](#3-whitening) | one 4 x 4 matrix per tetrode |
| 4 | [Seed templates](#4-seed-templates) | template bank |
| 5 | [Carry-forward matching](#5-carry-forward-matching) | accepted spikes and updated templates per window |
| 6 | [Re-seeding](#6-re-seeding) | units added to the bank |
| 7 | [Export](#7-export) | SpikeInterface sorting |
| 8 | [Analyzer and metrics](#8-analyzer-and-metrics) | not yet defined |
| 9 | [Curation](#9-curation) | not yet defined |

## 0. Before sorting

- `ttm-convert`: Open Ephys session to a blosc-zstd Zarr store, one
  tetrode per chunk, channels in tetrode order, the tetrode `ProbeGroup`
  attached.
- `ttm-experiment-params`: writes `experiment_params.json` from the
  manifest, including the experimenter's `[[observations]]`
  (`docs/pipeline.md`, "Observations"). Read them before curating.

## 1. Bad periods

`ttm-bad-periods`, per tetrode:

- Input: bandpass 300 to 6000 Hz plus the global median reference (not
  the aligned reference of step 2).
- Statistics in 0.1 s bins, features in 1 s bins: each wire's rms and
  the median wire-to-wire correlation of the tetrode.
- Baselines: 30 min centred running median, estimated twice more without
  bins where a wire exceeds twice its baseline (padded 30 s).
- A 1 s bin is flagged when a wire's rms exceeds 3 times its baseline
  (`loud`), when the median correlation is 0.25 above baseline and the
  median wire is twice its baseline (`common`), or when 12 of 16 tetrodes
  are flagged (`global`, then every tetrode is flagged).
- Flagged bins within 5 s of each other form one interval, padded 1 s.

<details>
<summary>Why</summary>

The seed and re-seed sorts, the per-window noise estimates and the
template rebuild all learn from the data, so stretches in which a
tetrode's wires carry an artefact have to be known before any unit
exists. The table is per tetrode because on TTM-003 groups of tetrodes
fail for minutes at a time while the others stay clean. Steps 5 and 6 use
it to drop fits and to keep those stretches out of the rebuild.

Details and the checks of each threshold:
[`docs/developer_notes/bad_periods.md`](developer_notes/bad_periods.md),
[`docs/reports/2026-10-26_ttm003_bad_periods.md`](reports/2026-10-26_ttm003_bad_periods.md).

</details>

## 2. Filter and reference

`ttm.sorting.preprocess_for_sorting` (defaults), then written once to a
float32 binary of all 64 channels:

- Bandpass 300 to 6000 Hz, float32.
- Aligned median reference (`cmr="global_aligned"`). Each wire is
  shifted by its time offset from the median (FFT, fractional samples),
  and its gain times the median of the shifted wires is subtracted.
- Offsets and gains per 1800 s window. Each window has six blocks of 20 s,
  each fitted only on samples within 15 samples of a median excursion
  beyond 5 MAD, and only if there are 3000 such samples. The estimate is
  the per-wire median over the qualifying blocks of the window and the two
  windows on each side, needing at least 3. Without them it falls back
  to the global median (shift 0, gain 1).
- No channel is excluded. CH22 is in the manifest's `badChannels` but is
  referenced and sorted like the others.
- Cost: the timing estimate 13 min (16 threads); the binary about 40 min
  (64 jobs).

<details>
<summary>Why</summary>

Reference and ground are tied to a cerebellar screw, so its signal is on
every wire. The channels are not sampled at the same instant (about 1 us
apart within each bank of 32, plus each electrode's own delay), so
subtracting the plain median leaves 7 to 30 uV of each large common event
where the noise is 4 to 5 uV. On the tetrodes with few units, the matcher
fitted units to that remainder. Aligning each wire before subtracting
removes 70 to 90 % of it. On the hybrid ground truth it raised
small-unit precision from 0.866 to 0.909, and over the 48 h it removed
most spikes on TT3, TT13, TT14 and TT16 without changing the
light-driven cells. Pooling over neighbouring windows stops one block
inside an artefact from setting a window's offsets.

Studies: [`investigations/common-reference/`](../investigations/common-reference/README.md)
(the remainder, the estimator, held-out validation in three animals,
pooling); the aligned reference against the global median on the hybrid
and over the 48 h:
[`investigations/carry-bakeoff/`](../investigations/carry-bakeoff/README.md),
sections "The aligned reference on the hybrid" and "The aligned
reference over the 48 h". The cross-tetrode median was also tested
(`ttm.sorting`, `cmr="cross_tetrode"`); it is slower and does not remove
the screw's remainder.

</details>

## 3. Whitening

- One 4 x 4 ZCA matrix per tetrode, W = C^(-1/2), from the covariance of
  the samples of the first 300 s that lie more than 30 samples from any
  crossing of -4 MAD.
- Used for all 48 h. Matching reads the traces as x W and holds the
  templates as T W.

<details>
<summary>Why</summary>

Whitening was the largest single gain in the bake-off. It removes fits
whose shape is shared by the four wires, such as far-field spikes and
artefacts. A matrix per window gained nothing, and the first 300 s
matrix keeps cross-wire correlations at 0.15 or less on most tetrodes
over the 48 h.

Studies: [`docs/reports/2026-10-29_ttm003_carry_bakeoff.md`](reports/2026-10-29_ttm003_carry_bakeoff.md),
sections 6 and 10 ("Whitening per window").

</details>

## 4. Seed templates

- MountainSort5 scheme 2 on window 0 (the first 1800 s) of the binary,
  each tetrode separately (`ttm.tracking.sort_chunk`). Detection threshold
  5.5, negative peaks, MountainSort5's own whitening with seed 42.
- Units kept with snr 5 or more and 100 spikes or more.
- Duplicates on a tetrode merged at template cosine 0.95 or more
  (shifts up to 10 samples).
- Templates: 90 samples (30 before the spike, 60 after), on the unit's
  four wires.

<details>
<summary>Why</summary>

The carry-forward follows units it is given and finds no new ones on its
own, so it needs a starting set from a clustering sorter. MountainSort5
was the sorter of the June pipeline and of every comparison since. The
0.95 merge removes MountainSort5's oversplit twins, which otherwise
compete for the same spikes in the matcher. A 0.9 merge joined distinct
cells.

Studies: [`docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`](plans/2026-10-22_ttm003_48h_sorting_plan.md)
(the June pipeline and its findings); `investigations/spike-sorting/`
(sorter, scheme and block-length comparisons). One whitening for
clustering and matching is an open follow-on.

</details>

## 5. Carry-forward matching

For each 1800 s window (96 windows), and in each window for each tetrode
separately:

1. Noise per template: the robust standard deviation of the matched
   filter's output (the trace projected on the unit-norm template) over
   six blocks of 10 s of the window.
2. Matching: circus-omp (spikeinterface 0.105.2) on the tetrode's four
   whitened channels with only its own templates. Lower amplitude bound
   max(0.2, 2.5 x noise / template norm), so matching is permissive and
   the rules below decide.
3. Admission, per fit:
   - matched-filter snr 4.5 or more;
   - cosine 0.5 or more between the template and the neighbour-subtracted
     snippet (every other fitted spike removed);
   - not inside a bad period of the tetrode (step 1);
   - fit amplitude 0.5 or more of the unit's median amplitude over its
     lead spikes outside laser pulses in this window (when it has 30);
   - refractory 1.5 ms per unit: of two spikes closer than that, the one
     with the higher snr stays.
4. Laser pulses (recordings with stimulation only): a unit whose rate
   from 15 ms after onset to offset is at least twice its rate outside
   pulses (Poisson p < 1e-4, 20 spikes or more) gets a relaxed rule
   inside pulses (onset to offset plus 100 ms): snr 2.5, relative floor
   0.3, cosine 0.5, then the refractory rule again.
5. Template rebuild, per unit: from its lead spikes (no spike of the unit
   in the 10 ms before) outside laser pulses (plus 100 ms) and outside the
   tetrode's bad periods, neighbour-subtracted snippets, the mean of at
   most 500. Only with 100 or more such spikes. New template = 0.8 x old +
   0.2 x that mean.

Cost: 8.4 to 9.8 min per window, about 15 h for 48 h.

<details>
<summary>Why</summary>

A clustering sorter run on each window separately has to match its units
across windows, and a unit's waveform drifts over 48 h. A carry-forward
keeps each unit's identity by matching its template window to window and
updating it slowly.

- Per tetrode rather than all 64 channels at once: in circus-omp,
  competition across tetrodes and an early stop shared by the whole probe
  hide spikes on quiet tetrodes. Per tetrode has the higher recall on
  every ground-truth unit. With the aligned reference it has the higher
  small-unit F1 (0.811 against 0.794), and over the 48 h it keeps each
  light-driven cell on one unit, where joint matching finds 20 to 35 %
  fewer of their onset spikes.
- Admission in noise units (snr), not as a fraction of the template: a
  floor tied to the template's own scale admits background fits of the
  right size, and with a rebuild it ratchets near-threshold templates
  upward until they lose their units.
- The relative floor and the refractory rule remove second fits of one
  spike and fits far smaller than the unit's typical spike.
- The relaxed rule inside pulses recovers the driven cells' smaller
  spikes during high-rate firing (on the hybrid, half-size in-pulse
  spikes 0.13 to 0.29 recall under the outside rule, 0.60 to 0.63
  relaxed). It is limited to units that keep firing through the pulses
  because applying it to every unit admitted false spikes.
- Blend 0.2 from lead spikes outside pulses and bad periods: every
  full-replacement rebuild let a small unit be captured by background.
  Frozen templates did as well over 6 h but cannot follow 48 h of drift.
  The blend keeps spike times on the waveform: over the 48 h the
  per-tetrode carry's units move by 3 samples at most, where June's
  rebuild moved 9 units by 5 samples or more.

Studies: [`docs/reports/2026-10-29_ttm003_carry_bakeoff.md`](reports/2026-10-29_ttm003_carry_bakeoff.md)
(sections 4 to 10: matcher, admission, rebuild recipes, laser pulses, the
48 h, follow-up arms) and
[`investigations/carry-bakeoff/README.md`](../investigations/carry-bakeoff/README.md)
(every arm and the per-tetrode precision analysis, including the overlap
rule, which was declined). Template super-resolution:
[`investigations/template-superres/`](../investigations/template-superres/README.md)
(no gain).

</details>

## 6. Re-seeding

- After windows 12, 24, ..., 84: MountainSort5 as in step 4 on that
  window.
- A cluster with snr 5 or more and 100 spikes or more joins the bank if
  its template's cosine (shifts up to 10 samples) to every template on its
  tetrode is under 0.8. New units get ids from 2000 and are matched from
  the next window on.

<details>
<summary>Why</summary>

Cells can appear or become large enough to sort after window 0. The
re-seed adds them without disturbing the units already carried. The 0.8
bar keeps a drifted copy of a carried unit out.

Studies: the June pipeline's re-seed
(`investigations/spike-sorting/_mp_common.py`,
`windowed_carry_forward_reseed`); the bake-off kept it unchanged.

</details>

## 7. Export

- The accepted spikes of all windows are combined into one
  `NumpySorting` with a `group` property per unit (the tetrode), saved as
  `<session>/sortings_seed42_pcafix/<run>/assembled_reseed_c12`.
- Units are never merged, split or removed by the carry. A unit keeps its
  id for the 48 h.

## 8. Analyzer and metrics

Not yet defined. To be decided:

- The SortingAnalyzer is built on the binary of step 2 (the recording the
  sort came from), never on another preprocessing of the session.
- Quality metrics (presence ratio, firing rate, isolation) computed per
  vigilance state, which needs sleep scoring of the recording.
- Whether the June production post-processing (merges guarded by
  cross-correlograms, single- and multi-unit labels,
  `investigations/spike-sorting/97_production_pipeline.py`) applies to
  these sortings.

## 9. Curation

Not yet defined beyond the existing tools (`docs/pipeline.md`, step 7).
Known items for TTM-003-48H:

- TT15 unit 42 is suspected to follow a bursting 4.7 kHz line on CH20 and
  CH22 (conjecture, from its mean waveform).
- The TT6 light-driven cell may be split between units 16 and 17.
- TT1, TT2, TT3, TT13, TT14 and TT16 showed no spikes during the live
  recording (manifest observations).
