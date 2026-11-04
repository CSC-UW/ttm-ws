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

Code: `noise_levels.py`, `ccg_pairs.py`, `zero_lag.py`. Results in
`results/`.

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
Conjecture, not tested here: it is the matcher fitting one event with two
similar templates, which would count up to 3.8 % of a unit's spikes twice.
Looking at the residual of the coincident events would tell.
