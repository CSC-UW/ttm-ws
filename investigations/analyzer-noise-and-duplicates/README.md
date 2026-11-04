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
`event_refractory.py`. Results in `results/`.

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
at lag 0 at 4.9 times and at lag +2 samples at 14.5 times the rate 5 to
25 ms away, falling to 1.0 to 1.2 times within 0.2 ms either side; 13 with
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
6.6 %) and 29 and 32 (6.1 against 5.8 %). Across the 18 pairs a third unit
is no more common at coincident events than near a's other spikes.
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
  still at 8.1 times and may extend beyond the 3 samples measured. In 8
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
spikes in its coincident events keep its refractory gap (conjecture). The
sharp pairs on TT6 and TT11 and 2003 and 2013 remain open. A peak 2 to 4
samples wide with a reproducible residual fits a waveform the template
bank lacks (conjecture), but the residual is small, the two mean spikes
accounting for 95 % or more of the mean waveform's energy, and could also
come from changes in a's or b's spike shape when both fire together.
