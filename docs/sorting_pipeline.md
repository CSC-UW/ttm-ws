---
title: Spike sorting pipeline for 48 h tetrode recordings
updated: 2026-11-03
status: needs_verification
---

# Spike sorting pipeline

What is done to a 48 h tetrode recording to get spike trains, in order,
with the parameters used. Each step has a collapsed "Why" note with the
measurements behind it, and each note ends with links to the studies that
made them, for readers who want the alternatives tested and the details.

## Setup and terms

The recordings are from rats with 16 tetrodes (64 channels) sampled at
30 kHz for about 48 h. Some recordings include optogenetic stimulation:
in TTM-003-48H, 207 ms light pulses at 1 Hz, in trains of 58 to 60
pulses, every 2 h. A light-driven cell is a unit whose firing rises
during the pulses.

- Tetrodes are TT1 to TT16. Channels keep their acquisition names (CH1
  to CH64; TT1 is CH34, CH36, CH38, CH40). A sorting's `group` property
  is the tetrode index from 0 (group 0 is TT1).
- Manifest: the per-recording TOML file in `subjects/<id>/experiments/`
  that names the session and lists its bad channels and tetrodes, the
  light schedule, and the experimenter's observations.
- Window: the recording is cut into consecutive 1800 s windows from its
  first sample, numbered from 0 (96 windows for 48 h).
- MAD noise: the median absolute deviation divided by 0.6745, a robust
  estimate of a trace's standard deviation.
- Template: a unit's mean waveform, 90 samples (30 before the spike time,
  60 after) on each wire of its tetrode. The bank is the set of templates
  carried from window to window, one per unit.
- Matched-filter snr: a spike's projection on its unit's template scaled
  to unit length, both whitened (step 3), divided by the MAD noise of
  that projection over six 10 s blocks of the window (step 5.1). Template snr is SpikeInterface's quality metric: the
  template's peak on its largest wire over that wire's MAD noise.
- Neighbour-subtracted snippet: the 90 samples around a spike after
  subtracting every other fitted spike, so overlapping spikes do not
  distort it. Cosines between a snippet and a template are taken over all
  samples and wires together.
- Lead spike: a spike with no spike of the same unit in the 10 ms
  before it, so it is not the smaller second spike of a burst.
- Hybrid ground truth: a copy of the first 12 windows (6 h) of
  TTM-003-48H's preprocessed data with 21 known units added: real
  templates moved to other tetrodes, with troughs of 5 to 14 MAD noise on
  their largest wire, bursts, and slow changes of size and shape. Three
  of them answer the real light pulses with one full-size spike followed
  by a train of half-size spikes, as a driven cell recorded in
  TTM-003 does. Small units have troughs of 6 MAD noise or less (11 of
  the 21), large units the rest. A unit's spike counts as found when it
  lies within 10 samples (0.33 ms) of a ground-truth spike; recall,
  precision and F1 (their harmonic mean) follow. With every choice
  below, the pipeline finds small units at recall 0.73 and precision
  0.91 (F1 0.81) and large units at 0.93 and 0.91 (F1 0.92). Comparisons
  in the "Why" notes say which other choices were in place when they
  were run.
- The earlier pipeline: the lab's previous sort of these recordings
  (2026-06), which matched templates on all 64 channels at once, admitted
  fits of at least 0.8 of the template's size, used the plain median
  reference and no whitening, and replaced each template every window
  by the mean of its accepted spikes.
- `<session>` is the session's directory on the data volume.

Run times are for TTM-003-48H on one server (Intel Xeon Platinum 8276M,
224 threads, 1.5 TB RAM; the GPU is not used).

## Steps

| # | Step | Command | Output |
|---|---|---|---|
| 0 | [Before sorting](#0-before-sorting) | `ttm-convert`, `ttm-experiment-params` | Zarr store, `experiment_params.json` |
| 1 | [Bad periods](#1-bad-periods) | `ttm-bad-periods` | table of bad intervals per tetrode |
| 2 | [Filter and reference](#2-filter-and-reference) | `ttm-sort binary` | preprocessed binary, int16 (about 660 GB) |
| 3 | [Whitening](#3-whitening) | `ttm-sort seed` | one matrix per tetrode |
| 4 | [Seed templates](#4-seed-templates) | `ttm-sort seed` | the starting bank |
| 5 | [Carry-forward matching](#5-carry-forward-matching) | `ttm-sort carry` | accepted spikes and updated bank per window |
| 6 | [Re-seeding](#6-re-seeding) | `ttm-sort carry` | units added to the bank |
| 7 | [Export](#7-export) | `ttm-sort export` | SpikeInterface sorting |
| 8 | [Analyzer and metrics](#8-analyzer-and-metrics) | | not yet defined |
| 9 | [Curation](#9-curation) | | not yet defined |

To run it on a recording that has been through step 0 (each command
resumes where it stopped):

```bash
ttm-bad-periods --subject TTM-003 --experiment TTM-003-48H
python subjects/ttm-003/scripts/laser_pulses.py --subject TTM-003 --experiment TTM-003-48H
ttm-sort all --subject TTM-003 --experiment TTM-003-48H \
  --pulses subjects/ttm-003/results/laser_pulses_2026-10-20_09-02-12.csv
```

- The pulse table is a CSV with one row per delivered pulse, `onset_s`
  and `offset_s` in seconds from the recording's first sample; the
  TTM-003 script writes it to `results/laser_pulses_<session name>.csv`.
  Each subject with stimulation needs a script that makes it from its
  TTL events. Without `--pulses`, every spike counts as outside pulses:
  step 5.4 is skipped, and in-pulse spikes enter the relative floor
  (5.3.4) and the template updates (5.5). A stimulated recording then
  loses most of its driven cells' smaller in-pulse spikes.
- `ttm-sort` writes to `<session>/sortings_carry/<run>/` (`--run`,
  default `main`). Its stages (`binary`, `seed`, `carry`, `export`) can
  be run one at a time. `--binary` sorts an existing preprocessed binary
  instead of writing one in step 2.
- End to end on TTM-003-48H: about 19 to 20 h (estimated from the
  steps' measured times): the offset fit 13 min, the binary 40 min, the
  seed sort 30 min, the 96 windows about 14.5 h, and the seven re-seed
  sorts about 3.5 h, run one after another.

## 0. Before sorting

- `ttm-convert`: the Open Ephys session to a Zarr store (blosc-zstd,
  lossless), raw int16 samples (0.195 uV per bit), unfiltered, channels
  in tetrode order with the tetrode layout attached. Steps 1 and 2 both
  start from it.
- `ttm-experiment-params`: writes `experiment_params.json` from the
  manifest, including the experimenter's observations
  (`docs/pipeline.md`, "Observations"). Read them before curating.

## 1. Bad periods

`ttm-bad-periods`, per tetrode, writes `bad_periods/bad_periods.csv`
(one row per tetrode and interval: tetrode, start and end in seconds,
`kind`, the wires involved):

- Input: the store bandpassed and referenced as in step 2 (aligned
  median, offsets fitted the same way), without the manifest's bad
  channels and tetrodes, computed from the store as it is read. When
  step 2's binary already exists, `--binary
  <session>/sortings_carry/<run>/binary` scans it instead.
- Per second: each wire's rms and the tetrode's median wire-to-wire
  correlation.
- Baselines: a centred 30 min running median of each wire's rms and of
  the tetrode's median correlation, computed three times; after the
  first and second pass, seconds in which any wire of the tetrode
  exceeds twice its baseline (padded 30 s on each side) are left out of
  all of the tetrode's baselines.
- A second is flagged `loud` when any wire's rms exceeds 3 times its
  baseline, `common` when the median correlation exceeds its baseline by
  0.25 and the median wire's rms is twice its baseline, and `global` on
  every tetrode when three quarters of the tetrodes (12 or more of 16)
  are flagged in it.
- Flagged seconds of a tetrode join into one interval when less than
  5 s apart, or when every second between them has a wire at twice its
  baseline; intervals are padded 1 s on each side.
- `bad_periods_params.json` records the input (store or binary,
  reference, channels left out). The scan's sums are cached in
  `bin_stats_0p1s.npz` and reused for the same input, so thresholds can
  be changed without rescanning.
- Cost: from the store, 57 min (the offset fit 13 min, the scan 43 min
  with 48 workers); from step 2's binary, 11 min.

<details>
<summary>Why</summary>

During movement or handling, a tetrode's wires can carry large signals
that are not spikes. Step 5 drops fits inside a tetrode's bad intervals
and leaves those intervals out of the template updates, so templates do
not learn the artefact. The table is per tetrode because groups of
tetrodes fail together for minutes while the others stay clean: on
TTM-003-48H the flagged time per tetrode is 0.01 to 1.3 % of the
recording, and 18 s are flagged on every tetrode.

The scan reads the aligned reference because that is what the sorter
reads. The detector was first run on the plain median of all channels.
On TTM-003-48H the aligned reference cut TT16's flagged time from
12,399 s (7.2 %) to 1,146 s, removed the shared-signal flags of TT13
(97 s) and nearly all of TT1's (55 to 3 s), and changed the other
tetrodes by 6 s or less; TT3's shared-signal flags (about 1,300 s)
remain. Leaving out CH22 (a bad channel of TT15) changes only TT15,
whose shared-signal flags on its three remaining wires rise from 6 to
1,067 s.

Details and checks of the thresholds:
[`docs/developer_notes/bad_periods.md`](developer_notes/bad_periods.md),
[`docs/reports/2026-10-26_ttm003_bad_periods.md`](reports/2026-10-26_ttm003_bad_periods.md),
[`investigations/carry-bakeoff/`](../investigations/carry-bakeoff/README.md)
("Bad periods on the aligned binary").

</details>

## 2. Filter and reference

Written once to a binary of int16 ADC counts:

- The manifest's bad channels, and every channel of its bad tetrodes,
  are removed first: they are neither in the reference nor sorted. A
  tetrode that loses a wire is sorted on the wires it has left; a
  tetrode with none is not sorted. `--keep-bad-channels` keeps them.
- Bandpass 300 to 6000 Hz: Butterworth of order 5 in each direction,
  applied forward and backward (zero phase), float32.
- Aligned median reference. Per window, each wire i has an offset
  tau_i (in samples, fractional) and a gain g_i. The wire is shifted
  earlier by its offset, y_i(t) = x_i(t + tau_i) (by FFT), and its gain
  times the median of all the recording's shifted wires (both banks)
  is subtracted:
  out_i(t) = y_i(t) - g_i median_j y_j(t). The output stays in the
  shifted time, which moves spike times by the offset: median 7
  microseconds, at most 40 (1.2 samples), on TTM-003-48H.
- Offsets and gains are fitted per window on six blocks of 20 s,
  centred at 1/7, 2/7, ..., 6/7 of the window. In each block, only
  samples within 15 samples of a large common event are used: those
  where the median of the wires departs from its own median by more than
  5 MAD noise, either sign. A block counts only if it has 3000 such
  samples. On those samples each wire is fitted by least squares as
  x_i = g m + b m' (m the median, m' its slope), which gives
  tau_i = -b / g; this is repeated three times on the shifted wires, and
  g_i is the least-squares gain on the shifted median.
- Each window's offsets and gains are the per-wire median over the
  counting blocks of that window and the two windows on each side (up to
  30 blocks; fewer at the ends of the recording), and need at least 3. A
  window without 3 takes the estimate of the nearest window that has
  them; if no window has them, the plain median is used (shift 0, gain
  1). Every window of TTM-003-48H had its own estimate.
- The result is clipped to the int16 range and rounded to the nearest
  ADC count (0.195 uV). The largest value on TTM-003-48H is about
  29,000 counts, inside the range of plus or minus 32,767.
  `--binary-dtype float32` keeps float32 instead.
- Cost: the offset and gain fit 13 min (16 threads); the binary about
  40 min (64 jobs).

<details>
<summary>Why</summary>

Reference and ground are tied to one skull screw over the cerebellum, so
whatever the screw picks up appears on every wire, in large common events
(tens to hundreds of uV). The headstage samples its channels one after
another, about 1 microsecond apart within each bank of 32, and each
electrode adds its own small delay, so on each wire the common event
arrives at its own offset from the median. Subtracting the plain median
then leaves 7 to 30 uV of each common event on the wires, where the
spike-band noise is 4 to 5 uV, and on tetrodes with few or no real cells
the matcher fitted units to these remainders. Shifting and scaling each
wire before subtracting removes 70 to 90 % of the remainder.

Measured against the plain median, with the matching of step 5:

- On the hybrid ground truth, small-unit precision rose from 0.866 to
  0.909 and small-unit F1 from 0.786 to 0.811; large units were
  unchanged (F1 0.911 and 0.917).
- Over the 48 h, on four tetrodes where the experimenter saw no spikes
  during the recording (TT3, TT13, TT14 and TT16), the spikes assigned
  fell from 7.8 M to 0.7 M.
  The light-driven cells' responses were unchanged: for example, the
  TT4 driven unit's rate in the first 10 ms after pulse onset ranged over
  the 24 two-hour blocks from 187 to 425 Hz with the plain median and
  from 202 to 429 Hz with the aligned one.

Pooling blocks over neighbouring windows keeps a single block inside an
artefact from setting a window's offsets. The estimator was tested on
TTM-003-48H and two other 48 h recordings (TTM-001-48H, TTM-002-48H),
estimating from some blocks of each window and measuring on the others.
Without pooling, a few windows whose estimate rested on one bad block
left their worst wire with more remainder than the plain median would
have (1.5 times in TTM-003, 2.7 times in TTM-002); with pooling, no
window did worse than 1.03 and 1.11 times.

The binary is int16 because float32 doubles its size (1.33 TB against
663 GB for TTM-003-48H) for precision the sort does not use. Rounding
moves a sample by at most 0.5 counts (rms 0.29), where the spike-band
noise is about 23 counts. The first four windows were carried on an
int16 copy and compared with the float32 run, next to a float32 rerun
that changed only the random seed used to pick spikes for the template
updates. By the fourth window, 94.5 % (int16) and 94.3 % (other seed)
of the float32 run's accepted spikes were identical, templates were
within 4.0 and 4.8 % of their peak, and the largest change of a unit
with over 10,000 spikes was 2.8 and 3.4 %. The spikes that differ are
borderline admissions: in the first window, those accepted by only one
of the runs have a median matched-filter snr of 5.3, against 6.4 for
spikes both accept.

Studies:
[`investigations/common-reference/`](../investigations/common-reference/README.md)
(the remainder, the estimator, held-out validation, pooling);
[`investigations/carry-bakeoff/`](../investigations/carry-bakeoff/README.md),
sections "The aligned reference on the hybrid", "The aligned
reference over the 48 h" and "An int16 binary".

</details>

## 3. Whitening

- One ZCA whitening matrix per tetrode, W = C^(-1/2), where C is the
  covariance between the tetrode's wires (4 x 4, or 3 x 3 with a wire
  removed) over the samples of the first 300 s that lie more than 30
  samples from any crossing of -4 MAD noise on any of its wires. The
  mean is not subtracted; no regularisation.
- One matrix for the whole recording. Matching reads the traces x
  (samples x wires) as x W and holds each template T (90 x wires) as
  T W.

<details>
<summary>Why</summary>

Whitening makes the noise uncorrelated across a tetrode's wires with
equal variance, so signals shared by all four wires (far-field activity,
artefacts) weigh less in the fit. On the hybrid ground truth, over 12
windows (plain median reference, templates held fixed, rules like
step 5's with the cosine bar at 0.6 unwhitened and 0.5 whitened),
whitening raised small-unit F1 from 0.60 to 0.79 and large-unit F1 from
0.79 to 0.92. A new matrix per
window changed no measure. With the first 300 s matrix, the largest
correlation between two whitened wires stays at 0.15 or less over the
48 h on most tetrodes; TT14 reaches 0.29 from hour 24, and single
windows of TT2, TT5 and TT10 reach 0.21 to 0.23.

Studies: [`docs/reports/2026-10-29_ttm003_carry_bakeoff.md`](reports/2026-10-29_ttm003_carry_bakeoff.md),
sections 6 and 10.

</details>

## 4. Seed templates

- MountainSort5, scheme 2, on window 0, each tetrode separately
  (`ttm.tracking.sort_chunk`): detection threshold 5.5 times the noise of
  MountainSort5's own whitened data, negative peaks, the random sample
  for MountainSort5's whitening drawn with seed 42; other parameters at
  the defaults of MountainSort5 0.5.9. It reads the preprocessed binary,
  bad periods included.
- Templates: the mean of up to 500 randomly chosen spikes (seeded) of
  each unit, from the preprocessed binary (not whitened), on its
  tetrode's wires. Step 5 whitens them (T W).
- Units with template snr 5 or more and 100 spikes or more are kept
  (TTM-003-48H: 54 of 59).
- Units of one tetrode whose templates have cosine 0.95 or more (at the
  best shift up to 10 samples) are merged: their spike trains are pooled,
  spikes within 0.3 ms of each other counted once, merging is
  transitive, and the merged unit's template is recomputed (TTM-003-48H:
  46 units after merging).

<details>
<summary>Why</summary>

Step 5 follows the units it is given and finds no new ones, so it needs
a starting set from a clustering sorter. MountainSort5 is the lab's
clustering sorter for tetrodes; the comparisons that chose it and its
settings are in `investigations/spike-sorting/`. MountainSort5 splits
some cells into near-identical clusters, which would then compete for
the same spikes in step 5; merging at cosine 0.95 removes them, while at
0.9, 19 % of merges joined two cells (the merged unit broke the
refractory period). The false-merge rate at 0.95 was not measured.

Studies: [`investigations/spike-sorting/`](../investigations/spike-sorting/README.md).

</details>

## 5. Carry-forward matching

For each window in turn, and in each window for each tetrode
separately, with the bank as it stands (window 0 uses the seed bank):

1. Noise: for each template, the MAD noise of the whitened traces
   projected on the unit-length template, over six 10 s blocks evenly
   spaced through the window.
2. Matching: circus-omp (SpikeInterface 0.105.4) on the tetrode's
   whitened wires with only its own templates, other parameters at
   their defaults (stop after 5 consecutive failed passes, `omp_min_sps`
   0.1, `relative_error` 5e-5, `rank` 5, `vicinity` 2), in 1 s chunks. A fit's amplitude is its scale relative to the
   template (1 is template size). The lower bound is
   max(0.2, 2.5 x noise / the template's L2 norm), with no upper bound,
   so matching is permissive and the rules below decide.
3. Admission, in this order:
   1. matched-filter snr 4.5 or more;
   2. the spike time not inside a bad interval of the tetrode (step 1);
   3. cosine 0.5 or more between the template and the
      neighbour-subtracted snippet;
   4. relative floor: among the unit's spikes that passed 1 to 3, take
      those that are lead spikes (counting only those spikes) and lie
      outside laser pulses (onset to offset plus 100 ms). With 30 or more,
      a spike needs an amplitude of at least 0.5 times their median
      amplitude; with fewer, at least 0.5 (half the template's size);
   5. refractory period 1.5 ms per unit: a run of spikes each closer than
      1.5 ms to the previous one keeps only its highest-snr spike, however
      long the run.
4. Laser pulses (only with a pulse table, only in windows with pulses).
   Candidates are fits outside bad intervals with matched-filter snr 2.5
   or more, amplitude at least 0.3 times the median of rule 3.4 (0.3
   with fewer than 30 lead spikes), and cosine 0.5 or more. A unit is
   driven in this window when its candidates from 15 ms after each onset
   (past the onset response, so the test asks whether it keeps firing
   through the pulse) to the offset number 20 or more, at least twice the count expected from its candidate rate
   outside the pulses (onset to offset plus 100 ms excluded), with a
   one-sided Poisson p under 1e-4. A driven unit's candidates from onset
   to offset plus 100 ms are accepted, and rule 3.5 is applied to its
   spikes again.
5. Template update, per unit with 100 or more accepted spikes in the
   window: from its accepted lead spikes outside laser pulses (plus
   100 ms), if there are 50 or more, up to 500 chosen at random (seeded by the window number); their
   neighbour-subtracted snippets (whitened, not re-aligned) are
   averaged, and the new template is 0.8 x the old + 0.2 x that mean.
   It is used from the next window.

Windows do not overlap, and matching does not read across a window's
edges. Fits in the first 30 or the last 60 samples of a window have no
complete snippet and are never accepted, so about 3 ms around each window
boundary is lost. Cost: 8.4 to 9.8 min per window (16 jobs), about 15 h for 48 h.

<details>
<summary>Why</summary>

A clustering sorter run on each window separately has to match its units
across windows while their waveforms drift. Carrying templates forward
keeps each unit's identity, and updating them slowly follows drift.

- One tetrode at a time, against matching all 64 channels in one
  circus-omp run: in a joint run, a spike on one tetrode can be
  suppressed by a larger coincident spike on another, and the search
  stops for the whole probe at once, which hides spikes on quiet
  tetrodes. On the hybrid, per tetrode has higher recall on every
  ground-truth unit (small units 0.732 against 0.681, large 0.929
  against 0.903) at lower precision (0.909 against 0.953, and 0.905
  against 0.941); F1 is 0.811 against 0.794 on small units and 0.917
  against 0.921 on large ones (both with every other choice of this
  document). Over the 48 h, the joint run finds 20 to
  35 % fewer of the driven cells' spikes in the first 10 ms after pulse
  onset, and divides one driven cell between two units.
- The admission bar is set in noise units (snr) rather than as a
  fraction of each template's size. A fraction of the template admits
  background fits of the right size and, combined with template updates,
  makes templates of units near the bar grow until they lose their
  spikes: on the hybrid under the earlier pipeline, three units' templates
  grew to 1.6 to 1.8 times their injected size and lost their spikes.
  With the bar in noise units, the 0.2 blend does as well as fixed
  templates (below).
- The relative floor compares a fit with the unit's own spikes in the
  same window, not with the template, and removes fits much smaller than
  the unit's typical spike, such as second fits of one spike. Its effect
  alone was not measured. Over TTM-003-48H the matcher made 841 M fits:
  17.5 % were accepted, 77.5 % failed the snr bar or the relative floor
  (recorded together), 3.5 % the cosine, 0.9 % the refractory period, and
  0.5 % fell in bad intervals.
- In-pulse rule: a driven cell firing at 100 to 200 Hz through a pulse
  fires smaller spikes than outside it. On the hybrid (window 0, plain
  median, unwhitened, fixed templates), the injected driven units'
  half-size in-pulse spikes were found at 0.13 to 0.29 recall under the
  normal rules and 0.60 to 0.63 under the in-pulse rule. Relaxing every
  unit raised the false spikes during pulses on the ground-truth units
  that fall silent in pulses from 2.9 to 11.3 per second in total;
  relaxing only driven units, to 3.6.
- Template update: on the hybrid over 12 windows (plain median,
  whitened, rules of step 5), the 0.2 blend gave
  small-unit F1 0.79 (precision 0.87), replacing the template each window
  0.72 (precision 0.74), and fixed templates 0.79 (0.89). Replacing let a
  small shrinking unit drift onto background spikes (precision 0.02 by
  window 6). Fixed templates were not run over the 48 h; there, the
  templates' cosine between window 0 and the last window has a median of
  0.92, and 22 of 46 starting units fall below 0.9. Over the 48 h the
  spike times of this pipeline's units stay within 3 samples of the same
  point on their waveforms, where the earlier pipeline, which replaced
  templates each window, moved 9 of its units by 5 samples or more (up
  to 1 ms).

Studies: [`docs/reports/2026-10-29_ttm003_carry_bakeoff.md`](reports/2026-10-29_ttm003_carry_bakeoff.md)
(sections 4 to 10) and
[`investigations/carry-bakeoff/README.md`](../investigations/carry-bakeoff/README.md).

</details>

## 6. Re-seeding

- After windows 12, 24, ..., 84: MountainSort5 as in step 4 on that
  window (bad periods included), keeping clusters with template snr 5 or
  more and 100 spikes or more; step 4's merge is not applied.
- In MountainSort5's order, a cluster joins the bank if its whitened
  template's cosine (best shift up to 10 samples) to every template on
  its tetrode, including clusters just added, is under 0.8. New units are numbered from 2000 (seed units
  from 0) and matched from the next window, so their spikes in the
  re-seed window itself are not recovered.

<details>
<summary>Why</summary>

Cells can appear, or grow large enough to sort, after window 0. Re-seeding
adds them without disturbing the units already carried, and the 0.8 bar
keeps a drifted copy of a carried unit out. It is stricter than the
seed's 0.95 merge, so a new cell whose template is within cosine 0.8 of
a carried unit's is not added. The interval and the bar were not
tuned. On TTM-003-48H the seed gave 46 units and the seven
re-seeds added 19.

</details>

## 7. Export

- The accepted spikes of all windows, as one SpikeInterface sorting in
  `<run>/sorting/`, with `group` (tetrode index) per unit. Spike times
  are sample indices from the recording's first sample.
- After seeding, no unit is merged, split or removed, and every unit
  keeps its id for the whole recording, including units that stop
  firing.
- The sorting does not contain new units' spikes in the window that
  found them (step 6), or spikes within about 3 ms of a window boundary
  (step 5).
- Bad intervals are not attached to the sorting. A tetrode is
  unobserved during its bad intervals, so firing rates and presence must
  be computed with `bad_periods/bad_periods.csv`.

## 8. Analyzer and metrics

Fixed:

- The SortingAnalyzer is built on the binary of step 2, the data the
  sort came from, never on another preprocessing of the session.

Open:

- Quality metrics (presence ratio, firing rate, isolation metrics such
  as isolation distance and L-ratio) computed per vigilance state, which
  needs the recording to be sleep scored, and excluding bad intervals.
- Whether an earlier post-processing script applies: it merges units
  whose cross-correlograms show they are one cell and labels units as
  single or multi-unit
  (`investigations/spike-sorting/97_production_pipeline.py`).

## 9. Curation

Not yet defined beyond the existing viewer (`ttm-curate`,
spikeinterface-gui; `docs/pipeline.md`). Notes for curating a particular
recording live with its subject (for TTM-003-48H,
`subjects/ttm-003/README.md`, "Curation notes").
