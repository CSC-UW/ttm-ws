---
title: Spike amplitudes of 0 or less and the fits in overlaps
updated: 2026-11-07
---

# Spike amplitudes of 0 or less and the fits in overlaps

On 2026-11-07 the analyzer stage stopped computing amplitude_cutoff
(`docs/sorting_pipeline.md`, step 8). This asks two questions behind that
decision. Why do SpikeInterface's spike amplitudes of the sorted units
include values of 0 or less? And since most of those spikes sit next to
another unit's spike: are the carry's accepted fits that lie within 1.5 ms
of another unit's accepted fit spikes of their unit, or the matcher fitting
the unit's template to what the other unit's subtraction left behind?

The sortings are the final ones of TTM-003's two 48 h recordings:
TTM-003-48H (session `2026-10-20_09-02-12`, 61 sorted units) and
TTM-003-48H-2 (session `2026-11-03_09-01-38`, 56 sorted units), each in
`sortings_carry/main` of its session directory, with the analyzers rebuilt
on 2026-11-07. The recordings are sampled at 30 kHz, so a sample is 0.033
ms.

Code, run from the repository root in this workspace's environment, each
script with an output folder as its last argument:

- `amp_le0.py` (folder A) and then `amp_le0_sub.py` (the same folder A,
  whose `amp_le0_snippets.csv` it reads);
- `overlap_fits.py` and `overlap_timing.py` (folder B);
- `amp_cutoff_check.py` and `positive_peak_units.py` (`results/`);
- `summarize.py A B results/`, which writes the other tables in `results/`.

The bulk outputs are not kept: `amp_le0_snippets.csv` (9.7 MB),
`amp_le0_sub.csv` (2.7 MB), the trace test's events (16 files, 10.3 MB), the
lag histograms (1.6 MB) and the counts by V's size (0.35 MB).

## Terms

- The carry, its fits, accepted fits, the fitted amplitude (1 is the
  template's size), the matched-filter snr (accepted fits need 4.5 or more
  outside the pulse rule) and the neighbour-subtracted snippet are those of
  `docs/sorting_pipeline.md` step 5; the matcher is circus-omp. Bad
  intervals are step 1's, the whitened traces step 3's, the multi-unit units
  step 8's. A multi-unit unit holds detected events and the
  refractory-rejected fits of its tetrode, which can include the tested
  unit's own.
- A frame is a spike's sample index.
- A spike amplitude is what SpikeInterface 0.105.4's `spike_amplitudes`
  extension stores: the trace in µV (step 2's filtered and referenced
  binary, not whitened) on the unit's main channel at one sample, the
  spike's frame plus the unit's peak offset, with nothing subtracted.
- A unit's main channel is chosen when the analyzer is built, as the channel
  of its tetrode where the mean of 100 random spikes is most negative.
- The peak offset is where the main channel's mean template (up to 5,000
  random spikes) has its largest absolute value, relative to the frame; the
  analyzers use SpikeInterface's default peak sign, "both". For seven units
  that value is a positive peak (`results/positive_peak_units.csv`): 40, 43,
  2008 and 2010 of TTM-003-48H and 37, 2003 and 2013 of TTM-003-48H-2, on
  TT15, TT3, TT16 and TT2.
- As in amplitude_cutoff, a unit's amplitudes are negated when their median
  is negative. "0 or less" is after that: at or above 0 µV for a unit read
  at a trough, at or below 0 µV for the seven read at a positive peak.
- The tetrodes fall into three groups: TT4 to TT12; TT1, TT2 and TT16, whose
  units are smaller (median amplitudes of 10 to 45 µV, against 42 to 167 µV
  on TT4 to TT12); and TT15, whose units follow a signal synchronous with
  the acquisition clock (`docs/reports/2026-10-26_ttm003_bad_periods.md`).
  TT3 and TT13 hold one unit of 34 spikes and one of 429 on TTM-003-48H and
  are left out.
- In an overlap, U is the unit whose fit is tested and V the other unit. The
  lag is U's frame minus V's. Pairs are ordered: an overlap can be tested
  twice, once with each unit as U.
- A fit's size is its fitted amplitude over the median of its unit's
  accepted fits in the same 1800 s window: small under 0.7, mid 0.7 to 1,
  large 1 or more.
- A percentile places a fit among its unit's accepted fits in the same
  window, on a scale of 0 to 1.

## Amplitudes of 0 or less

Counted per unit over all spikes by `amp_le0.py`
(`results/amp_le0_units.csv`; by group `results/le0_by_group.csv`).

| | TTM-003-48H | TTM-003-48H-2 |
|---|---|---|
| TT4 to TT12: units with any | 42 of 42 | 37 of 37 |
| TT4 to TT12: share of a unit's spikes, median (max) | 0.024 % (0.54 %) | 0.040 % (0.36 %) |
| TT4 to TT12: share of all their spikes | 0.066 % | 0.050 % |
| TT1, TT2, TT16: share of a unit's spikes, median | 1.2 % | 0.37 % |
| TT15: share of a unit's spikes, median | 1.3 % | 1.2 % |

- On TT4 to TT12 they fall in a median of 75 (TTM-003-48H) and 87
  (TTM-003-48H-2) of a unit's 96 windows, and a unit's busiest window holds
  a median 5 % of them.
- Four units have a peak offset of 15 to 24 samples (0.5 to 0.8 ms), all of
  them read at a positive peak: 43 (TT15) and 2010 (TT16) of TTM-003-48H and
  2003 (TT2) and 2013 (TT15) of TTM-003-48H-2, with 7.2, 7.1, 0.75 and 7.6 %
  of their spikes at 0 or less (`results/le0_far_peak_shift_units.csv`,
  which also lists TT3's unit 2008, offset 51 samples, 34 spikes).
- For 9 units, 6 of TTM-003-48H (TT3's unit 2008 among them) and 3 of
  TTM-003-48H-2, the main channel is not the channel where the mean
  template's peak is largest, taken on the side the amplitudes are read
  (columns `main_ch` and `template_extremum_ch` of
  `results/amp_le0_units.csv`).

`amp_le0.py` then takes up to 100 spikes at 0 or less per unit and up to 100
of the unit's other spikes as controls, and `amp_le0_sub.py` refits them
(`results/le0_context.csv`, `results/le0_subtraction.csv`). Units with a
peak offset of 10 samples or more are kept apart from their tetrode group
there. On TT4 to TT12 that is 7,404 spikes at 0 or less and 7,900 controls
from the two recordings together.

| TT4 to TT12 | at 0 or less | controls |
|---|---|---|
| another unit's sorted spike within 1 ms | 82 % | 21 % |
| an event of the tetrode's multi-unit unit within 1 ms | 24 % | 3.6 % |
| inside a laser pulse or 10 ms after it | 0.11 % | 0.15 % |
| best template match 0.1 ms or more off the amplitude's sample | 78 % | 5.3 % |
| best template match at the search limit of 10 samples | 39 % | 1.6 % |
| best template match, cosine (median) | 0.51 | 0.83 |
| cosine of the carry's whitened snippet with nothing subtracted under 0.307 | 87 % | 4.9 % |
| carry percentile: that cosine | 0.002 | 0.50 |
| carry percentile: cosine of the neighbour-subtracted snippet | 0.12 | 0.50 |
| carry percentile: snr | 0.26 | 0.50 |
| carry percentile: fitted amplitude | 0.20 | 0.49 |

The template match is the cosine of the tetrode's µV trace, 0.5 ms either
side of the amplitude's sample, with the unit's mean template, best over
shifts of up to 10 samples either way. For the spikes at 0 or less the
median absolute best shift is 8 samples (0.27 ms). The 82 % here pools up to
100 sampled spikes per unit; over every spike of TT4 to TT12 it is 91 % on
TTM-003-48H and 85 % on TTM-003-48H-2 (next section).

`amp_le0_sub.py` fits each spike's trace, 3 ms before to 5 ms after the
frame on the tetrode's channels, by least squares with the analyzer's mean
template (1 ms before to 2 ms after the frame) of every sorted spike of the
tetrode whose template reaches that span, so frames from 5 ms before to 6 ms
after, the spike's own included. It then takes the other spikes' fitted
share off the amplitude's sample. Multi-unit events are not subtracted.

- On TT4 to TT12, 92 % of the spikes at 0 or less have another unit's spike
  in that span, against 58 % of the controls. With the other spikes' share
  taken off, 78 % become positive (85 % of those with another unit's spike,
  1.3 % of the rest), 17 % reach half the unit's median amplitude, and the
  median is 0.24 of it. The spike's own fitted scale is a median 0.49,
  against 0.95 for the controls, all of which stay positive.
- On TT1, TT2 and TT16 (without 2003 and 2010), 43 % have another unit's
  spike in the span and 27 % become positive; on TT15 (without 43 and 2013),
  97 % and 65 %.

## Are the fits in overlaps spikes?

If the carry fitted U's template to what was left of V after V's
subtraction, a V that may also rescale, shift, widen or narrow should absorb
most of what U's fit explains. If U's fit is a spike of U, U should explain
as much beyond that V as a real spike of U added at the same lag does.

### Trace test

`overlap_fits.py` runs on the whitened traces the carry matched, with the
carry's templates for the window (the bank it matched with), in windows 8,
20, 32, 44, 56, 68, 80 and 92 of each recording, thirty 10 s blocks per
window, leaving out a tetrode's blocks that overlap its bad intervals. Every
fit of the tetrode, whatever its admission outcome, is subtracted except the
event's own.

- An event is an accepted fit of U with exactly one accepted fit of another
  unit V of the tetrode at a lag of 1.5 ms (45 samples) or less, and no
  other fit with snr 3 or more reaching its snippet, which runs from 80
  samples before V's frame to 190 after. Up to 60 events per pair and block
  are kept, and only in blocks where both U and V have lone fits: accepted
  fits with no other fit of snr 3 or more within 9 ms.
- Each event gets two controls, built from lone fits of the same block drawn
  at random with replacement. The twin is a lone fit of V plus a lone fit of
  U, scaled to the event's U amplitude and placed at the event's lag; both
  are traces, so the twin holds real spikes of U and V, with the U snippet's
  noise scaled along. V alone is the same lone fit of V by itself.
- On each snippet the test fits, by least squares, V's template with its
  first and second time derivatives at V's frame, with and without U's
  template at U's frame, each frame refined by up to 2 samples. The gain is
  the energy U's template adds, in whitened noise units (sum of squares). On
  V alone it is what U's template finds where U has no fit.
- The null level is the 99th percentile of the gain on V alone over all
  22,819 events of both recordings and every tetrode: 56.9.

`results/overlap_trace_summary.csv`, `results/overlap_trace_tails.csv`, per
pair `results/overlap_trace_per_pair.csv`:

| TT4 to TT12 | TTM-003-48H | TTM-003-48H-2 |
|---|---|---|
| events | 9,517 | 12,920 |
| gain on the event's snippet (median) | 271 | 364 |
| gain on the twin | 275 | 368 |
| gain on V alone | 4.7 | 6.3 |
| event over twin, median per event | 1.01 | 1.01 |
| events below the null level | 0.15 % | 0.22 % |
| twins below the null level | 1.2 % | 0.77 % |
| U's least-squares scale over the carry's amplitude | 1.00 | 1.00 |

- Over the events of TT4 to TT12 of both recordings, U's template finds 100
  or more in 0.12 % of the V-alone snippets and 150 or more in none of them;
  98.8 % of the events reach 100.
- Over the 157 ordered pairs of both recordings with 30 or more events, the
  median event-over-twin ratio per pair is 0.89 to 1.26.
- The ratio is 0.99 to 1.02 at every lag, 1.00 to 1.01 at every size of V's
  fit, and 0.98 and 1.02 for the 5 % of events with the lowest cosine of the
  whitened snippet with nothing subtracted (under 0.307), the range where 87
  % of the sampled spikes at 0 or less lie.
- Small fits of U (312 and 532 events) explain more than their twins (ratio
  1.28 and 1.24), and their least-squares scale is 1.09 and 1.10 times the
  carry's amplitude. The twin is scaled to the carry's amplitude, and a
  scale 1.09 to 1.10 times too small alone makes an energy ratio of about
  1.19 to 1.21.
- 1.6 and 1.9 % of the small fits are below the null level, against 7.7 and
  6.0 % of their twins. Those below are all within 0.3 ms of V: 8.2 % of 61
  and 8.5 % of 117 events there, against 9.8 and 17.1 % of their twins.
- On TT15 (40 and 233 events) the ratio is 1.06 and 1.05; on TT1, TT2 and
  TT16 (74 and 35) it is 0.97 and 0.87, with 15 % and 0 % of events below
  the null level against 22 % and 2.9 % of twins.

### Timing and size

`overlap_timing.py` counts, over all 96 windows of both recordings, U's
accepted fits at each lag from V's accepted fits within 5 ms, by the size of
U's fit. Per ordered pair on TT4 to TT12 with 20 or more U fits expected
within 1.5 ms, each count is divided by U's mean count per sample at 2 to 5
ms either side (`results/overlap_lag_profile.csv`). The table gives the
median over pairs: 172 (168 for small fits) on TTM-003-48H, 148 on
TTM-003-48H-2.

| lag of U from V (ms) | -1.5 to -0.73 | -0.70 to -0.33 | -0.30 to 0.30 | 0.33 to 0.70 | 0.73 to 1.5 |
|---|---|---|---|---|---|
| TTM-003-48H, small | 0.92 | 1.10 | 1.06 to 1.26 | 1.67 | 1.32 |
| TTM-003-48H, large | 1.09 | 1.43 | 1.13 to 1.26 | 1.39 | 1.02 |
| TTM-003-48H-2, small | 0.95 | 1.18 | 0.96 to 0.99 | 1.65 | 1.38 |
| TTM-003-48H-2, large | 0.96 | 1.33 | 0.72 to 0.96 | 1.28 | 0.93 |

The -0.30 to 0.30 ms column gives the range over its three bins (-0.30 to
-0.13, -0.10 to 0.10 and 0.13 to 0.30 ms).

- Small fits of U are more common 0.33 to 1.5 ms after V's fit than large
  ones, and not before it. In the trace test, small fits at those lags
  explain 1.18 to 1.38 times their twins' energy, none is below the null
  level, and their least-squares scale is 1.12 to 1.17 times the carry's
  amplitude.
- Sharp peaks at single lags are not particular to small fits
  (`results/overlap_lag_peaks.csv`). The test takes the 297 pairs with 500
  or more small fits within 1.5 ms and fewer small than large fits within 5
  ms, and thins the large fits binomially to the small fits' count. Within
  1.5 ms one lag stands above a 21-sample running mean (Poisson, Bonferroni
  over 91 lags, p < 0.01) in 68 and 57 % of the pairs for small fits and in
  59 and 56 % for thinned large fits. The two lag distributions differ
  (5-sample bins, chi-square p < 0.001) in 97 and 96 % of the pairs.
- How common U's fits within 1.5 ms are depends on the size of V's fit
  (`results/overlap_by_v_size_summary.csv`). The excess over U's rate at 2
  to 5 ms, per 1,000 pairings of a V fit with another unit of its tetrode,
  is given below for V's fit under 0.8, 0.8 to 0.95, 0.95 to 1.05, 1.05 to
  1.2 and over 1.2 times V's median. On TTM-003-48H it is 18.7, 17.0, 18.6,
  20.9 and 28.2; on TTM-003-48H-2, 15.8, 6.8, 5.4, 6.0 and 11.5. It is
  raised at both ends. Small fits of U account for 0.4 to 1.4 of it. In the
  trace test the event-over-twin ratio is 1.00 to 1.01 at every size of V's
  fit.

## amplitude_cutoff

SpikeInterface 0.105.4 computes amplitude_cutoff from the
`amplitude_scalings` extension when the analyzer has it and otherwise, with
a warning, from `spike_amplitudes`; the analyzer stage computes only
`spike_amplitudes`. In `amplitude_cutoff`
(`spikeinterface/metrics/quality/misc_metrics.py`) the metric:

- puts a unit's amplitudes in 100 bins between their minimum and maximum and
  smooths the counts with a Gaussian with a standard deviation of 3 bins;
- finds the last bin whose smoothed count is at least that of the lowest
  bin;
- sums the smoothed counts above that bin as the spikes missing below the
  lowest bin, so the upper tail stands in for a lower tail cut off by
  detection;
- returns that sum over itself plus the unit's spike count, at most 0.5, and
  NaN for units of fewer than 500 spikes (TT3's and TT13's).

`amp_cutoff_check.py` recomputes it over all spikes of the sorted units
(`results/amp_cutoff_check.csv` per unit; the table below is
`results/amp_cutoff_summary.csv` and the next one
`results/amp_cutoff_bands.csv`; the medians printed in
`results/amp_cutoff_check.txt` also count TT3's and TT13's units). The
stored metric covers the periods outside bad intervals; recomputed over all
spikes it differs from the stored value by at most 0.000024. Amplitudes are
given over the unit's median amplitude.

| 59 and 56 sorted units | TTM-003-48H | TTM-003-48H-2 |
|---|---|---|
| stored metric, median (max) | 0.00007 (0.13) | 0 (0.0004) |
| spikes in the lowest bin, median | 2 | 1 |
| smallest amplitude, median | -0.76 | -0.96 |
| largest amplitude, median (max) | 10.5 (64) | 5.1 (17) |
| 99.9th percentile, median | 2.7 | 2.6 |
| share of the bins above the 99.9th percentile, median | 65 % | 37 % |
| on amplitudes above 0: metric, median (max) | 0.031 (0.42) | 0.005 (0.15) |
| on amplitudes above 0: spikes in the lowest bin, median | 1,403 | 592 |
| on amplitudes above 0: width of the lowest bin, median | 0.105 | 0.051 |
| on amplitudes above 0: share of the bins above the 99.9th percentile, median | 74 % | 46 % |

- With every amplitude, the most negative ones lie far below the rest, the
  lowest bin is nearly empty, and the metric is about 0.
- Above 0 the lowest bin, from 0 to about 0.1 (TTM-003-48H) and 0.05
  (TTM-003-48H-2) of the median, is well filled, and most spikes from 0 to
  0.1 of the median have another unit's spike within 1 ms. On TT4 to TT12,
  over every spike, that share falls with amplitude:

| amplitude over the unit's median | 0 or less | 0 to 0.1 | 0.1 to 0.25 | 0.25 to 0.5 | over 0.5 |
|---|---|---|---|---|---|
| TTM-003-48H | 91 % | 79 % | 64 % | 40 % | 22 % |
| TTM-003-48H-2 | 85 % | 70 % | 53 % | 33 % | 18 % |

## Reading

Conjecture, not confirmed.

- On TT4 to TT12 the amplitudes of 0 or less are mostly spikes of their unit
  that overlap another unit's spike. SpikeInterface reads one sample of the
  trace with nothing subtracted, and the other spike's waveform at that
  sample lifts it to 0 or above; taking the other spike off brings most of
  them back above 0.
- The carry's fits in clean overlaps are spikes of their unit, not its
  template fitted to what V's subtraction left: on the whitened traces U
  explains as much beyond a flexible V as a real spike of U does, at every
  lag, size and size of V. No group of fits falls below the null level more
  often than its twins.
- The carry's amplitude for a small fit 0.33 to 1.5 ms after V's fit is 11
  to 14 % below its least-squares size, and small fits gather there. Why the
  carry fits them small was not measured.
- The refit in µV, a joint least-squares fit with the analyzer's 48 h mean
  templates at the sorted frames, puts the spikes at 0 or less at about half
  their unit's usual scale. The carry, with its per-window templates on
  whitened traces, places them at a median percentile of 0.20 of its fitted
  amplitudes. Any amplitude made the way the refit is may place spikes in
  overlaps lower than the carry does.
- Neither amplitude_cutoff on every amplitude nor on those above 0 estimates
  the spikes the carry missed. The first is set by the few most negative
  amplitudes, the second by the height of a lowest bin filled mostly with
  overlapping spikes, and the carry admits fits on their snr over the whole
  template, not on one sample.
- The trace test covers clean overlaps only, one other fit and nothing else
  strong nearby. Overlaps of three or more spikes were not tested.