---
title: Noise levels and duplicate units in the final TTM-003-48H sorting
updated: 2026-11-04
---

# Noise levels and duplicate units in the final TTM-003-48H sorting

Two questions left open for the analyzer stage (`docs/sorting_pipeline.md`,
step 8), asked of the final sorting (`sortings_carry/main`, 61 units,
138.7 M spikes) and its int16 binary: how SpikeInterface's noise levels,
the denominator of SNR, should be sampled; and whether the merge step of
`investigations/spike-sorting/97_production_pipeline.py` has anything to
merge.

Code: `noise_levels.py`, `ccg_pairs.py`, `zero_lag.py`, `double_fits.py`,
`coincident_shape.py`, `third_units.py`, `consistent_waveform.py`,
`event_refractory.py`, `leftover_waveform.py`, `lag_profiles.py`,
`unit13_partners.py`. Results in `results/`.

## Noise levels

`get_noise_levels` (method "mad") takes the MAD of each of 20 random
0.5 s chunks and averages the 20 values; the chunks can fall in bad
intervals. `noise_levels.py` computed the MAD of every 0.5 s chunk of a
grid on every channel in one pass over the binary (15.2 min with 48
workers). The reference is the mean over every chunk that does not overlap
a bad interval of the channel's tetrode. Draws of 20 to 10,000 chunks,
100 each, are subsets of that table, either of any chunk (as
SpikeInterface draws) or of clean chunks only. Five calls of
`get_noise_levels` itself (20 chunks, seeds 0 to 4) spread by 2.3 % on the
median channel, against 3.5 % in the simulated draws of any chunk, with
their largest errors on channels of TT1, TT15 and TT16 as in the simulation.

Relative error against the reference, per channel over the 100 draws: the
95th percentile of the absolute error, for the median channel and the worst
channel, and the largest bias of any channel:

| Chunks | Seconds | Draw | 95th percentile, median channel | 95th percentile, worst channel | Largest bias |
|---|---|---|---|---|---|
| 20 | 10 | any chunk | 6.7 % | 431 % | 40 % |
| 20 | 10 | clean chunks | 4.6 % | 9.4 % | 1.1 % |
| 100 | 50 | any chunk | 3.1 % | 104 % | 18 % |
| 100 | 50 | clean chunks | 2.2 % | 4.8 % | 0.3 % |
| 500 | 250 | any chunk | 1.3 % | 59 % | 20 % |
| 500 | 250 | clean chunks | 1.0 % | 2.0 % | 0.1 % |
| 2,000 | 1,000 | any chunk | 0.8 % | 35 % | 19 % |
| 2,000 | 1,000 | clean chunks | 0.5 % | 1.0 % | 0.0 % |
| 10,000 | 5,000 | any chunk | 0.4 % | 26 % | 19 % |
| 10,000 | 5,000 | clean chunks | 0.2 % | 0.4 % | 0.0 % |

The five real calls with 20 chunks overestimated CH36 (TT1) by up to
166 %, CH34 (TT1) by 147 %, CH28 (TT16) by 147 % and CH24 (TT15) by 80 %.
Averaging MADs is not robust to one chunk inside an artifact, and more
chunks do not remove the bias: the mean over every chunk, bad ones
included, is 1.19 times the clean mean on CH18 (TT15), 1.18 on CH36,
1.17 on CH24, 1.14 on CH28 and 1.12 on CH34, and at most 1.06 on any
other channel. Channels of TT5 to TT12, with almost no bad intervals,
are unaffected.

Two other definitions were looked at and not pursued. The MAD of all the
clean samples pooled is 0.95 to 1.02 times the mean of chunk MADs, since
the noise level varies over the 48 h; on int16 data it is also quantized to
whole counts, 6 to 9 % of the noise level on the quietest channels (MADs
of 11 to 17 counts). The median of chunk MADs is robust to artifact chunks
but is 0.89 to 1.04 times the clean mean, a different quantity.

## Cross-correlogram pairs

`ccg_pairs.py` judges all 126 same-tetrode pairs as `97_production_pipeline.py`
did (`_assignment_eval.ccg_verdict_pair`: coincidences within 1.5 ms per ms
over those 5 to 25 ms away per ms, in the 1800 s windows where both units
have at least 5 spikes; below 0.3 "duplicate", above 0.7 "distinct"), with
the template cosine (`ttm.tracking.cosine_from_templates`, 10-sample shift,
templates of up to 5,000 spikes per unit in µV). 200 random pairs of units on
different tetrodes, which the per-tetrode matcher fits separately, are the
control.

| Template cosine | duplicate | ambiguous | distinct |
|---|---|---|---|
| below 0.5 | 1 | 1 | 24 |
| 0.5 to 0.9 | 1 | 1 | 85 |
| 0.9 to 0.95 | 0 | 0 | 12 |
| 0.95 or more | 0 | 0 | 1 |

The two "duplicate" pairs are on TT15 (units 40 and 41, cosine 0.07; 41
and 42, cosine 0.51), whose units follow a signal locked to the acquisition
clock (`investigations/bad-periods-ttm003/`); their templates differ too much
for one cell. The control pairs gave 188 "distinct" and 12 without enough
co-activity, ratio median 1.05.

The merge rule of `97_production_pipeline.py` joins a pair with cosine 0.90
or more when the verdict is "duplicate", and any pair with cosine 0.95 or
more whatever the verdict (its docstring says only where the
cross-correlogram cannot decide, the code does not check). On this sorting
it would join one pair, units 21 and 22 on TT8 (cosine 0.967), whose
cross-correlogram ratio is 1.10, "distinct". No pair of the 13 with cosine
0.90 or more has a refractory dip; their ratios are 1.02 to 1.59.

## Coincidences within 0.1 ms

Ratios above 1 for similar units mean more coincidences within 1.5 ms than
5 to 25 ms away. At one-sample resolution, units 13 and 14 (TT6) coincide
at lag 0 at 4.9 times and at lags +1 and +2 samples at 15.3 and 14.5
times the rate 5 to 25 ms away, falling to 1.0 to 1.2 times within 0.2 ms either side; 13 with
15 and 13 with 16 look the same. Units 8 and 10 (TT5), 21 and 22 (TT8) and
9 and 2005 (TT5) show the opposite, a dip to 0.06 to 0.44 times within
0.3 ms; for the first two the rate from 0.3 to 1.5 ms is 1.2 to 1.4
times.

`zero_lag.py` counts, for every pair, the coincidences within 3 samples
(0.1 ms) beyond those expected from the rate 5 to 25 ms away, as a fraction
of the smaller unit's spikes:

| Pairs | Tetrodes | Pairs with an excess above 1 % | Largest excess | Median observed / expected |
|---|---|---|---|---|
| same tetrode | TT5 to TT12 | 18 of 76 | 3.8 % (13 and 14) | 1.35 |
| same tetrode | TT1 to TT4, TT13 to TT16 | 5 of 50 | 3.1 % | 0.74 |
| different tetrodes | both in TT5 to TT12 | 0 of 70 | 0.4 % | 0.99 |
| different tetrodes | at least one in TT1 to TT4, TT13 to TT16 | 12 of 130 | 40 % (unit 0 on TT1 with unit 44 on TT16) | 0.86 |

Across tetrodes, zero-lag coincidences occur only where the even-channel
tetrodes are involved: 40 % of unit 44's spikes (TT16) fall within 0.1 ms of
a spike of unit 0 (TT1), 26 % of unit 2010's (TT16) of unit 43 (TT15), 14 %
of unit 2014's (TT15) of unit 2002 (TT2). One cell never appears on two
tetrodes, so those units follow a shared signal; this fits the absence of
live spikes on TT1, TT2, TT3, TT13 and TT16 during the recording. Among
TT5 to TT12 no pair across tetrodes shares spikes at zero lag, so the
zero-lag excess within those tetrodes is not a signal common to the probe.

## Double-counted spikes?

The zero-lag excess within TT5 to TT12 looked like a known artifact of
template matching: Kilosort sometimes fits a second template to the
residual left after subtracting a large spike, which gives zero-lag peaks
between units (MouseLand/Kilosort2 issue 29, "Double-counted spikes",
2019). The Allen Institute's `ecephys_spike_sorting`
(`kilosort_postprocessing`) removes the later spike of any pair within 5
samples (0.166 ms) in one unit or in units whose peak channels are within 5
channels; SpikeInterface's `remove_duplicated_spikes` does the same within
a unit.

That is not what happens here. `double_fits.py` takes up to 1,000
coincident events (a spike of unit a with one of unit b within 3 samples)
from each of the 18 pairs with an excess above 1 %, and from 3 pairs on the
same tetrodes without an excess, and compares each event's waveform with a's
mean spike, b's (shifted by the lag) and their sum, the mean spikes taken
over spikes with no partner within 1.5 ms.

- The sum is the closest of the three for 69 to 96 % of the coincident
  events in 17 of the 18 pairs, and for 50 % in the other (units 2003 and
  2013 on TT5). Spikes of a with no partner are closest to a's mean spike
  69 to 94 % of the time.
- At a's peak sample and channel, coincident events are 1.28 to 1.94 times
  a's mean spike (2003 and 2013: 0.98); spikes without a partner are 0.89
  to 0.99 times.
- Amplitude alone cannot tell two spikes from one large spike of a whose
  residual was fitted with b. `coincident_shape.py` averages each pair's
  coincident events at their most common lag and fits that mean waveform
  with free amplitudes. With both units' mean spikes the fit leaves 0.5 to
  4.4 % of its energy, at amplitudes 0.68 to 1.68 for a and 0.79 to 1.36
  for b. The better of the two single-unit fits leaves 2.9 to 45 %, 2.3 to
  48 times as much. The control pairs' chance coincidences give the same
  picture (0.4 to 0.8 % with both, 4.5 to 16 % with one). One large spike
  of a would be fitted as well by a alone.
- Laser pulses are not involved: at most 0.3 % of a pair's coincident
  events fall within a pulse or the 50 ms after it, which cover 0.21 % of
  the recording.

So the coincident events are not one spike of a counted twice: they carry
two waveforms at about their usual amplitudes, as chance coincidences do.
Removing the later spike of each pair, as the Allen module does, would not
undo a double count here, and where the events are two cells firing
together it would remove real spikes. Why the two waveforms coincide
within 0.1 ms is taken up in the next section.

Within units nothing is duplicated: no unit has two spikes within 3
samples. 173,958 intervals are exactly 45 samples (1.5 ms), the shortest
the carry's refractory rule keeps; none is shorter and none crosses an
1800 s window edge.

## One waveform of their own?

Among the explanations left for the excess are two cells firing within
0.1 ms of each other, and a waveform the template bank lacks that the
matcher rebuilds from two templates. Three scripts look at the same 18
pairs. The controls are the 8 same-tetrode pairs on TT5 to TT12 whose
coincidences within 3 samples are 0.85 to 1.15 times the expected rate. In
one of them, 37 and 38 on TT12, that total hides an uneven lag profile:
0.02 and 0.04 times the expected rate at -3 and -2 samples, 1.4 to 1.9
times at +1 to +3.

`third_units.py` counts how often another unit of the tetrode fires within
1.5 ms of a's spike. In the 18 pairs that happens for 22 to 48 % of the
coincident events (median 33 %), for 19 to 48 % of a's spikes with b 5 to
25 ms away and none within 1.5 ms (median 35 %), and for a median 31 % of
a's spikes with no spike of b within 1.5 ms. Within 0.1 ms a third unit is
rarer at coincident events (median 2.2 %) than at the other two (5.0 and
4.7 %), except in 4 pairs: 11 and 2013 (13.7 % against 8.1 % with b 5 to
25 ms away), 30 and 32 (7.4 against 3.5 %), 2003 and 2013 (7.5 against
6.6 %) and 29 and 32 (6.1 against 5.8 %), so the coincident events are
rarely three spikes at the same instant. Within 1.5 ms the medians are
close, but per pair a third unit is at coincident events -3.4 to +16.3
percentage points more often than at a's spikes without b (11 and 2003:
+16.3, 11 and 2013: +15.0), and -9.7 to +8.3 points more often than at
a's spikes with b 5 to 25 ms away.
`coincident_shape.csv`'s column `third_unit_within_1p5ms` is a different
measure, the mean number of other units with a spike within 1.5 ms per
event at the most common lag.

`consistent_waveform.py` takes up to 3,000 coincident events per pair,
snippets read from the int16 binary, and fits each with free amplitudes
alpha and beta on a's and b's mean spikes (spikes with no partner within
1.5 ms).

- The lag profiles fall in two groups, by how far the best lag within 3
  samples stands above the lowest. In 8 pairs, called sharp below, the
  best lag is 3.7 to 14 times the lowest and 5.1 to 17 times the expected
  rate per sample: four on TT6 (unit 13 with 12, 14, 15 and 16), three on
  TT11 (32 with 29, 30 and 33), and 2003 and 2013 on TT5. Their peaks are
  2 to 4 samples wide at half height; 2003 and 2013 peak at +2 with +3
  still at 8.1 times; `lag_profiles.py` below finds 4.5 times at +4 and
  3.0 at +5 and +6. In 8
  others, called broad below, on TT5 and TT10, the best lag is 1.13 to
  1.93 times the lowest and every lag from -3 to +3 samples is at 1.7 to
  4.2 times the expected rate, so their excess is wider than 0.1 ms. 29
  and 33 and 11 and 2005 are in between (2.9 and 2.7 times).
- The mean waveform of the events at the most common lag, fitted with both
  mean spikes, leaves 1.7 to 4.8 % of its energy in 7 of the 8 sharp pairs
  (13 and 15: 0.9 %), 0.6 to 1.7 % in the broad pairs and 0.4 to 1.4 % in
  the controls. Adding the mean spikes' time derivatives, which absorb
  sub-sample misalignment, lowers these by at most 0.4 percentage points.
  The residual is not noise: the mean residuals of odd and of even events,
  fitted separately, have a cosine of 0.85 to 0.99 in the sharp pairs
  (0.96 to 0.99 in five). Controls reach 0.92 with smaller residuals, so
  the cosine shows the residual is real, not that it differs in kind from
  theirs.
- The correlation of alpha and beta within 1800 s windows, minus the same
  for a's spikes without a partner (fitted with b's mean spike at lag 0,
  the correlation noise gives through the overlap of the two mean spikes),
  does not separate excess pairs from controls: -0.34 to +0.18 against
  -0.19 to +0.25.
- The mean residual energy of coincident events is 0.65 to 0.90 of what
  two independent spikes would leave, taken as the mean for a's spikes
  without a partner plus that for b's, minus that for noise snippets; the
  controls give 0.84 to 0.96. 15 of the 18 excess pairs fall below every
  control, but sharp (0.65 to 0.90) and broad pairs (0.73 to 0.87) do not
  differ. The reference counts twice what a's and b's spikes carry besides
  their own waveforms, such as other units' spikes nearby, which may be
  why every pair falls below 1 (conjecture).

`event_refractory.py` asks whether the spikes in coincident events are
spikes of a and of b. A real spike of a unit sits in that unit's
refractory gap; a spike of another cell rebuilt from the two templates
does not. The measure is the rate of the unit's other spikes 1.5 to 3 ms
either side over the rate 20 to 50 ms away. Unit 11 (TT5) has a clear gap,
0.48 to 0.50 for its spikes with no spike of the partner within 1.5 ms,
and keeps it in its coincident events with 9, 2003, 2005 and 2013 (0.43
to 0.53), although 61 to 75 % of those coincidences are excess; those are
presumably unit 11's own spikes (conjecture). Unit 2013's gap is partly
filled in its coincident events (0.61 to 0.63 without the partner, 0.77
with 11, 0.89 with 2003). That is not decisive: in all 8 control pairs,
both units' coincident spikes have more of their own spikes at 1.5 to
3 ms than their spikes without the partner do (unit 17 with 18: 1.53
against 1.14), as in periods of higher firing. The units of the sharp
pairs on TT6 and TT11 have no gap beyond 1.5 ms (0.89 to 1.02 on TT6,
1.16 to 1.61 on TT11), so the test cannot tell there.

The broad pairs on TT5 and TT10 look like two cells firing together over a
range wider than 0.1 ms: their mean coincident waveform is fitted by the
two mean spikes about as well as chance coincidences are, and unit 11's
spikes in its coincident events keep its refractory gap (conjecture). For
the sharp pairs (on TT6, on TT11, and 2003 and 2013 on TT5) these tests do
not decide between two cells and a waveform the template bank lacks. Their
residual is small, the two mean spikes accounting for 95 % or more of the
mean waveform's energy, and could also come from changes in a's or b's
spike shape when both fire together. The last section weighs the
explanations.

## The leftover waveform

`leftover_waveform.py` looks at R, what the two mean spikes leave of the
mean coincident waveform W at the most common lag; "left" below is the
fraction of W's energy in R. It builds the mean spikes from 2,000 spikes
without a partner (`consistent_waveform.py` used 1,000) and averages up to
5,000 events at the most common lag (`consistent_waveform.py` used those
among its 3,000 events), so its fractions differ slightly: 0.8 to 4.5 % in
the sharp pairs. Cosines below are template cosines at the best shift
within 10 samples. With a figure path as argument, the script also draws
the waveforms; the figures are not committed.

- R's largest sample is 15 to 49 % of the peak of a's mean spike in the
  sharp pairs. R spreads over the four channels like a's and b's mean
  spikes, as it does in the controls (cosine of the per-channel RMS 0.86
  to 1.00 in the sharp pairs, 0.83 to 0.99 in the controls). Its cosine
  is at most 0.20 with the mean spike of any other unit of the tetrode,
  and at most 0.29 with the change in a's or b's spike within bursts (the
  mean of spikes without a partner less than 10 ms after the unit's
  previous spike, minus the mean of those more than 50 ms after). Adding
  that burst change to a fit that already has the first time derivatives
  lowers left by 1.7 percentage points in 29 and 32 and by at most 0.2 in
  the other sharp pairs.
- In 7 of the 8 sharp pairs R, on its main channel, is largest 4 to 5
  samples (0.13 to 0.17 ms) before a's spike and smallest 0 to 1 sample
  before it: W has a larger positive deflection just before the trough,
  and a trough 1.08 to 1.26 times as deep as the fit's. 2003 and 2013
  differ (largest 5 samples after a's spike, smallest 1 after). Adding the
  mean spikes' first and second time derivatives to the fit, which absorbs
  a difference in width or in the blur from spike-time jitter, leaves 0.2
  to 0.7 % in the four pairs on TT6, down from 0.8 to 2.9 %. The controls
  leave 0.2 to 1.1 % with these terms and 0.3 to 1.2 % without. In the
  three pairs on TT11 it leaves 1.1 to 2.2 %, down from 1.8 to 4.5 %; 32
  and 33 (2.2 %) and 29 and 32 (1.9 %) stay above every control. 2003 and
  2013 leave 1.1 %, down from 1.5 %.
- One of the three controls on TT6 shows a similar R: 14 and 16, at 0.94
  times the expected rate, largest 6 samples before a's spike and smallest
  2 before, a trough 1.15 times as deep, left 1.1 % falling to 0.3 % with
  the derivatives. The other two, 12 and 15 and 12 and 16, have troughs
  1.02 times as deep and leave 0.5 and 0.3 % without the derivatives and
  0.3 % each with them.
- The mean of the 10 % largest of 20,000 spikes of a without a partner,
  fitted with a's mean spike alone (amplitude 1.76 to 1.93 in the sharp
  pairs), leaves 1.1 to 6.8 % of its energy. On TT6 adding the first and
  second derivatives removes at most 0.1 percentage points of that, so the
  largest spikes of those units are not narrower. In 30 and 32 it removes
  0.7 points, about what it removes from that pair's coincident waveform.
  The residual of the large spikes has a cosine of at most 0.20 with R.
  Selecting on amplitude also selects noise and overlapping spikes, so this
  is a rough calibration.

So in the four sharp pairs on TT6, and in one of the three control pairs
there, the coincident waveform has a narrower and deeper trough than the
sum of the two mean spikes, which the largest single spikes of those units
do not show. What narrows it is open: spike times of two summed waveforms
may be more precise than those of single spikes, leaving their mean less
blurred, or the events may be a narrower waveform the template bank lacks
(both conjecture). On TT11, 32 and 33 and 29 and 32 keep 1.9 to 2.2 % that
the derivative terms do not absorb.

## Lag profiles out to 1.5 ms

`lag_profiles.py` counts b's spikes at each lag from -45 to +45 samples
(1.5 ms) around a's, over the mean count per lag 150 to 750 samples (5 to
25 ms) away on either side, for all 76 same-tetrode pairs on TT5 to TT12
and for the 70 pairs of units on two different tetrodes among TT5 to TT12
that `zero_lag.py` also counted (cross-tetrode pairs below). Means over
ranges of the absolute lag (`results/lag_profiles.csv`; each lag in
`results/lag_profiles_by_lag.csv`), for the 26 pairs compared above:

| Pairs | 0 to 3 samples (0 to 0.1 ms) | 4 to 8 (0.13 to 0.27 ms) | 9 to 15 (0.3 to 0.5 ms) | 16 to 30 (0.53 to 1 ms) | 31 to 45 (1.03 to 1.5 ms) |
|---|---|---|---|---|---|
| 8 controls | 0.88 to 1.09 | 0.73 to 1.73 | 1.17 to 2.12 | 1.21 to 1.67 | 0.86 to 1.29 |
| 8 broad | 1.97 to 3.55 | 1.68 to 3.55 | 1.69 to 2.70 | 1.35 to 1.80 | 0.71 to 1.53 |
| 8 sharp | 2.70 to 7.37 | 1.05 to 2.48 | 1.32 to 1.82 | 1.21 to 1.52 | 0.99 to 1.35 |
| 2 in between | 4.07 to 6.79 | 1.85 to 2.01 | 1.64 to 1.73 | 1.21 to 1.41 | 0.93 to 1.14 |

All 76 pairs fire together more often 0.3 to 0.5 ms apart than 5 to 25 ms
apart (1.08 to 2.70 times, median 1.52), and 75 do 0.53 to 1 ms apart.
Cross-tetrode pairs do so much less, 0.94 to 1.38 times (median 1.08) at
0.3 to 0.5 ms. Within 0.1 ms most same-tetrode pairs fall below their own
rate at 0.3 to 0.5 ms: 44 of the 76 do, and the 58 pairs without an excess
above 1 % sit at a median 0.68 times it (quartiles 0.29 and 0.99, lowest
0.08). That fits a matcher that misses spikes of two units that nearly
coincide (conjecture). The 8 controls show this dip partly because of how
they were chosen: their rate within 0.1 ms was fixed near the expected
one while their rate at 0.3 to 0.5 ms is 1.17 to 2.12 times the expected
one, which puts them at 0.51 to 0.76 times their own rate at 0.3 to
0.5 ms.

Measured against each pair's own rate at 0.3 to 0.5 ms, the 18 pairs with
an excess form a continuum: 1.01 to 1.75 times within 0.1 ms in the broad
pairs, 2.35 and 4.14 in the two in between and 1.87 to 4.73 in the sharp
pairs. In the sharp pairs the part above that rate holds 74 to 93 % of
their excess within 0.1 ms but 22 to 65 % of their excess within 1.5 ms.
Coincidences moved into the peak from nearby lags would leave a shortfall
next to it. The rate at 4 to 8 samples is below the rate at 0.3 to 0.5 ms
in the four sharp pairs on TT6, by at most 31 % of the peak's excess, and
not at all in the others. On TT6 the pairs without a peak fall below
their 0.3 to 0.5 ms rate at 4 to 8 samples by similar amounts (1.4 to 7.8
expected counts per lag summed over those lags, against 1.7 to 3.3 in the
sharp pairs), so the shortfall looks like the general dip near zero lag;
at most 31 % of a TT6 peak could be coincidences moved in from there.

On TT11 unit 32 is in three sharp pairs (with 29, 30 and 33), and 29 and
33 are one of the two pairs in between.

## Unit 13 and its four partners

Unit 13 is in all four sharp pairs on TT6, a hub of them, while the other
pairs there coincide within 0.1 ms at 0.35 to 1.43 times the expected
rate, none with an excess above 1 % (`results/zero_lag_pairs.csv`). If
one source, a cell
or other signal whose waveform the template bank lacks, were built by the
matcher from 13 and whichever of 12, 14, 15 and 16 fits best, 13's
coincident events with each partner would be the same waveform. If they
are 13's spike plus the partner's, they differ as the partners do.
`unit13_partners.py` averages, for each partner, up to 5,000 of 13's
spikes with a spike of the partner at the most common lag and no other
TT6 unit within 1.5 ms, aligned on 13's spike. It compares the averages
two by two, as the energy of their difference at the best whole-sample
shift over their mean energy, with the same for the sums of 13's and each
partner's mean spike (`results/unit13_partners.csv`). The ratio is the
averages' difference, less what the events' spread contributes to it, over
the sums' difference.

| Partners | Averages differ by | Sums differ by | Ratio, after the events' spread |
|---|---|---|---|
| 12 and 14 | 2.2 % | 2.7 % | 0.81 |
| 12 and 15 | 10.9 % | 9.5 % | 1.14 |
| 12 and 16 | 22.6 % | 14.8 % | 1.53 |
| 14 and 15 | 6.9 % | 3.4 % | 2.02 |
| 14 and 16 | 19.9 % | 10.8 % | 1.85 |
| 15 and 16 | 7.5 % | 5.6 % | 1.34 |

The spread of the events contributes at most 0.01 %. Two hubs whose
partners coincide with them at about the expected rate, where the events
are two spikes, give ratios of 1.70 (16 with 12 and 14) and 1.03 (12 with
15 and 16). Allowing a scale and a sub-sample shift, one average leaves
0.7 to 17.6 % of the other unfitted, against 1.8 to 13.8 % for the sums.
Each average is fitted by 13's mean spike at 0.69 to 1.04 times and the
partner's at 0.79 to 1.24 times, leaving 0.9 to 2.8 %
(`results/unit13_partners_pairs.csv`).

A fixed shared source would let the averages differ only through the
coincidences of two cells among the events, 9 to 28 % of them going by
the rate at the most common lag (5.1 to 15.3 times the expected rate)
against the co-firing at 0.3 to 0.5 ms (1.32 to 1.45 times), so by well
under the sums' difference (conjecture). They differ as much as two spikes
would: each average looks like 13's spike plus that partner's. A shared
source could still give this if its waveform varied from event to event
by as much as the partners differ, since the matcher picks for each event
the partner that fits it, and a separate missing waveform for each partner
would give it too. So the test makes one source split into 13 and a
partner unlikely on TT6 (conjecture) without ruling it out. Unit 32 on
TT11 was not tested this way.

## What the excess is

The coincident events are not one spike counted twice: they carry two
waveforms, mostly at amplitudes near the units' own (0.68 to 1.68 times
for a, the ends being 2003 and 2013 and 11 and 2005, and 0.79 to 1.36
for b). That leaves two explanations (conjecture): two neighbouring cells
firing together, or one source, a cell or other signal whose waveform the
template bank lacks, which the matcher builds from two templates.

What tells them apart:

- On TT6, 13's coincident events with each of its four partners look
  like 13's spike plus that partner's, as two spikes would, which makes
  one source shared by the four pairs unlikely (conjecture).
- Unit 11's spikes in its coincident events with 9, 2003, 2005 and 2013
  (broad and in-between pairs) keep its refractory gap, as its own spikes
  would (conjecture). Unit 2013's gap is partly filled.
- Same-tetrode pairs on TT5 to TT12 fire together 0.3 to 0.5 ms apart at a
  median 1.52 times the rate 5 to 25 ms apart, pairs on two different
  tetrodes there at 1.08 times. Local co-firing of this kind is what
  neighbouring cells with correlated input would show, though matching
  within a tetrode could also add some, for example by fitting a
  neighbour's template to the tail of a spike 10 to 15 samples later
  (both conjecture). In the broad pairs the rate within 0.1 ms is at most
  1.75 times each pair's own rate at 0.3 to 0.5 ms.
- Within 0.1 ms 44 of the 76 same-tetrode pairs fall below their own rate
  at 0.3 to 0.5 ms, as if the matcher missed coincidences there rather than
  added them (conjecture), and at most 31 % of a sharp peak could be
  coincidences moved in from nearby lags.

What remains open is the precision of the sharp pairs. Their mean rate
within 0.1 ms is 1.87 to 4.73 times their rate at 0.3 to 0.5 ms, and the
best lag 5.1 to 17 times the expected rate, in a peak 2 to 4 samples
(0.07 to 0.13 ms) wide. Correlated synaptic input would be expected to
give the broader co-firing rather than to lock two cells within 0.1 ms
(conjecture; whether electrically coupled cells fire this precisely is
not checked here). On TT11 two sharp pairs keep 1.9 to 2.2 % of their mean
coincident waveform that the derivative terms do not absorb, more than any
control, and unit 32, in three of them, was not tested as 13 was. 2003
and 2013 look more like one waveform the template bank lacks than two
full spikes (conjecture), though not like one spike counted twice either
(each mean spike alone leaves 12 to 13 % of their mean coincident
waveform, both together 1.5 %): the sum of the two
mean spikes is closest for only 50 % of their coincident events, which are
no larger than a's mean spike at a's peak (0.98 times), with a's fitted
amplitude 0.67 to 0.69.
