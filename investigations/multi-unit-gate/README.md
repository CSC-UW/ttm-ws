---
title: Multi-unit cosine gate on the final TTM-003-48H sorting
updated: 2026-11-05
---

# Multi-unit cosine gate on the final TTM-003-48H sorting

The multi-unit step of `investigations/spike-sorting/97_production_pipeline.py`
(97 below) pools detected events that no sorted spike claims into one extra
unit per tetrode, when an event's waveform matches a template of its tetrode
well enough. 97 set "well enough" at a cosine of 0.55. That value was chosen
on TTM-001 in `investigations/spike-sorting/91_neural_noise_calibration.py`,
against random stretches of background trace: at 0.55, 3.7 % of those
stretches passed and 1.1 % of the spikes of 35 well-isolated units failed.

This asks where to set the gate on the final sorting of TTM-003's 48 h
recording (session `2026-10-20_09-02-12`, TTM-003-48H below;
`sortings_carry/main` in the session directory, 61 units, 138.7 M spikes),
per tetrode, from the share of real spikes it may lose. Events that one
sorted spike claims stand in for real spikes. The question was asked for the
analyzer stage (`docs/sorting_pipeline.md`, step 8). The replay of 97's
steps on the same sorting is in `investigations/unclaimed-events/`.

Code: `gate.py`, `summarize.py`, `refractory.py`, `lags.py`, `figure.py`,
`overlaps.py`, `neighbours.py`.
Results in `results/`.

## Method

`gate.py` works on the sort's int16 binary, one 1800 s window at a time.

- Events are negative peaks at 5.5 times each channel's noise
  (SpikeInterface's "locally_exclusive" detection, radius 40 µm, so a
  trough gives one event per tetrode), kept at 7 times the noise or more.
  Events inside a bad interval of their tetrode (the session's
  `bad_periods/bad_periods.csv`), or in the first 30 or last 59 samples of
  a window, are left out.
- Each channel's noise is fixed for the whole recording, as the analyzer
  stage is to have it: the mean MAD (as a standard deviation) of every
  0.5 s chunk outside the tetrode's bad intervals
  (`investigations/analyzer-noise-and-duplicates/results/noise_levels_channels.csv`).
  The stage is to draw 2,000 of those chunks; on every channel, 95 % of
  such draws come within 1.0 % of this mean. The replay in
  `investigations/unclaimed-events/` used each window's own noise instead.
- An event is unclaimed when no sorted spike of its tetrode is within
  0.5 ms of it, and a reference event when exactly one is. Events with two
  or more are left out: 3.72 M, 13.3 % of the events on TT5 to TT12.
- Both kinds get the score the gate uses: the best cosine between the
  event's trace (1 ms before to 2 ms after the trough, the tetrode's wires,
  no shift) and the templates of its tetrode. A template is the mean of up
  to 500 random spikes of a unit with 100 or more spikes in the window, as
  in 97.
- Background stretches are 1,000 random times per tetrode and window, at
  least 1 ms from every detection and sorted spike of the tetrode and
  outside its bad intervals, scored the same way.

The cutoff losing x % on a tetrode is the x-th percentile of its reference
events' scores (`summarize.py`). Counts and rates of kept events would drop
an event within 0.3 ms of the previous one, as 97 does; at every cutoff here
that removes none.

## Reference events and cutoffs

On TT5 to TT12:

| Tetrode | Reference events | Unclaimed events | Reference below 0.55 | Cutoff losing 1 % | Losing 5 % | Losing 10 % | Losing 5 %, per 6 h block |
|---|---|---|---|---|---|---|---|
| TT5 | 3.82 M | 169,454 | 19.0 % | 0.31 | 0.43 | 0.49 | 0.37 to 0.48 |
| TT6 | 4.91 M | 96,891 | 3.8 % | 0.44 | 0.57 | 0.64 | 0.56 to 0.60 |
| TT7 | 793,386 | 41,857 | 4.6 % | 0.44 | 0.56 | 0.62 | 0.52 to 0.58 |
| TT8 | 3.76 M | 97,132 | 2.9 % | 0.47 | 0.59 | 0.65 | 0.58 to 0.62 |
| TT9 | 31,082 | 4,477 | 9.2 % | 0.44 | 0.52 | 0.55 | 0.51 to 0.53 |
| TT10 | 3.74 M | 638,257 | 22.2 % | 0.34 | 0.43 | 0.48 | 0.42 to 0.45 |
| TT11 | 1.76 M | 106,733 | 11.4 % | 0.36 | 0.48 | 0.54 | 0.47 to 0.49 |
| TT12 | 4.26 M | 151,521 | 7.3 % | 0.43 | 0.52 | 0.57 | 0.51 to 0.55 |

97's 0.55 would lose 2.9 to 22.2 % of the reference events. The blocks are
eight runs of 12 windows.

Whether low-scoring reference events are real spikes matters for this
reading of the cutoffs. `refractory.py` counts the claiming unit's spikes
0.5 to 1.5 ms from each reference event, per ms, over its spikes 20 to 50 ms
away, per ms (`results/refractory_by_cosine.csv`). The ratio is 0.004 to
0.023 at a score of 0.7 or more and rises as the score falls: at most 0.13
from 0.2 to 0.4 (TT10 0.23), and 0.00 to 0.41 below 0.2 (TT10 0.41, TT7
0.30, TT8 0.25; TT9's zeros there rest on 16 to 119 events). The claiming
spike sits within 2 samples of its usual place relative to the event's
trough (the median over the unit's reference events) for 77 to 99.8 % of
reference events in every band from 0.2 up and 68 to 98 % below 0.2. So the
low-scoring reference events are mostly the claiming units' own spikes, at
their own troughs, scoring low, with more doubt below 0.2 (conjecture).

The cutoffs rise with amplitude on six of the eight tetrodes
(`results/by_amplitude.csv`). On TT6, TT8, TT10, TT11 and TT12 the cutoff
losing 5 % rises from 0.41 to 0.56 at 7 to 8 times the noise to 0.60 to 0.79
at 15 or more. On TT7 it rises from 0.53 to 0.71 at 12 to 15, then is 0.42
at 15 or more (394 reference events). On TT5 it stays at 0.42 to 0.47 up to
15 and is 0.55 above. TT9 has 185, 20 and 45 reference events at 10 to 12,
12 to 15 and 15 or more. 63 to 87 % of the unclaimed events are 7 to 8 times
the noise, and 86 to 95 % are 7 to 9. Weighting the reference events to the
unclaimed events' amplitudes lowers the cutoffs by 0.003 to 0.035.

A cutoff per amplitude band (7 to 8, 8 to 9, 9 to 10, 10 to 12, 12 to 15 and
15 or more times the noise, each the percentile of that band's reference
events, or the tetrode's cutoff where a band has fewer than 100) keeps more unclaimed events than one cutoff per tetrode
(`results/banded.csv`). Over TT5 to TT12, at the share losing 10 %, one
cutoff keeps 947,625 and the banded cutoffs 987,139; 17,090 events pass only
the one cutoff and 56,604 only the banded ones. At the share losing 5 % it is
1,088,839 against 1,114,664, with 13,508 and 39,333. With one cutoff per
tetrode losing 10 %, the reference events of 7 to 8 times the noise lose
10.5 to 16.6 % and those of 12 or more 0.5 to 5.8 % (38.5 % on TT9).

## What each cutoff keeps

The rate of the multi-unit unit on each tetrode, with no upper limit on the
score (97 stopped at 0.8). The sorted rate is all of the tetrode's units
together, over the time outside its bad intervals:

| Tetrode | Sorted rate | At 0.55 | At 1 % | At 5 % | At 10 % |
|---|---|---|---|---|---|
| TT5 | 84.3 Hz | 0.73 Hz | 0.96 Hz | 0.91 Hz | 0.84 Hz |
| TT6 | 108.5 Hz | 0.47 Hz | 0.53 Hz | 0.45 Hz | 0.39 Hz |
| TT7 | 32.9 Hz | 0.19 Hz | 0.22 Hz | 0.18 Hz | 0.15 Hz |
| TT8 | 71.3 Hz | 0.39 Hz | 0.48 Hz | 0.33 Hz | 0.24 Hz |
| TT9 | 8.0 Hz | 0.017 Hz | 0.023 Hz | 0.020 Hz | 0.017 Hz |
| TT10 | 95.1 Hz | 1.90 Hz | 3.52 Hz | 3.07 Hz | 2.65 Hz |
| TT11 | 68.4 Hz | 0.48 Hz | 0.60 Hz | 0.55 Hz | 0.49 Hz |
| TT12 | 76.1 Hz | 0.75 Hz | 0.86 Hz | 0.79 Hz | 0.71 Hz |

Over TT5 to TT12 that is 0.85 M events at 0.55, 1.24 M at 1 %, 1.09 M at
5 % and 0.95 M at 10 %, or 0.90, 1.32, 1.16 and 1.01 % of the sorted spikes
there (`results/summary.json`).

Background stretches pass 0.55 in 0.1 to 1.5 % of cases, the cutoff losing
1 % in 0.8 to 12.8 %, 5 % in 0.2 to 4.9 % and 10 % in 0.1 to 2.9 %
(`results/tetrodes.csv`). Their centres are more than 1 ms from any
detection, though a trace can reach one 1 to 2 ms later, so this says how
background scores, not how many non-neural events of 7 or more pass. Events
with a detection (5.5 times the noise or more, outside bad intervals) on 3
or more other tetrodes within 0.5 ms make up 2.1 to 5.9 % of the unclaimed
events below the cutoff losing 5 %, 1.8 to 2.8 % of those at or above it and
1.3 to 2.3 % of the reference events (on TT9, 16.2, 4.5 and 5.0 %).

## Sorted spikes around the unclaimed events

`lags.py` counts the sorted spikes of the tetrode at each lag from -3 to
3 ms of an event, over their rate 20 to 50 ms away (`results/lags.csv`,
`results/lag_hist.csv`). Within 0.5 ms a reference event has only its
claiming spike and an unclaimed event none, so the table starts at 0.5 ms.
The unclaimed events are grouped by where their score falls among the
tetrode's cutoffs. Ranges are over TT5 to TT12:

| Events | Number | With a sorted spike 0.5 to 1.5 ms away | Spikes 0.5 to 1.5 ms before, over the rate 20 to 50 ms away | After |
|---|---|---|---|---|
| Reference | 23.07 M | 1.8 to 23.2 % | 0.39 to 0.95 | 0.48 to 1.28 |
| Unclaimed, below the cutoff losing 1 % | 63,091 | 8.7 to 77.0 % | 1.18 to 2.34 | 3.82 to 7.17 |
| Unclaimed, 1 to 5 % | 154,392 | 8.8 to 68.8 % | 1.56 to 2.50 | 2.97 to 4.95 |
| Unclaimed, 5 to 10 % | 141,214 | 4.9 to 58.9 % | 1.17 to 2.15 | 1.96 to 3.50 |
| Unclaimed, at or above 10 % | 947,625 | 2.0 to 26.6 % | 0.76 to 2.35 | 0.60 to 1.14 |

Below the cutoff losing 10 %, a sorted spike follows the event 0.5 to 1.5 ms
later 1.96 to 7.17 times as often as 20 to 50 ms later. Below the cutoff
losing 1 %, it does so most often 1.07 to 1.30 ms later (0.60 ms on TT8).
The excess falls off after 1.5 ms: from 2 to 3 ms after the event, every
group of unclaimed events is at 0.98 to 2.94 times the rate 20 to 50 ms
away, and the reference events at 0.99 to 2.79.

At or above the cutoff losing 10 %, a sorted spike follows 0.60 to 1.14
times as often, within the reference events' 0.48 to 1.28; per tetrode that
is 0.77 to 1.41 times the reference events' ratio. A sorted spike precedes
these events more often than it precedes reference events on every
tetrode, 1.15 to 1.65 times as often, and 6.0 times on TT8, where it does
so most often by 1.43 ms (6.6 times the rate 20 to 50 ms away).

Against the unit each matches best, the unit's spikes are 1.33 to 6.80 times
as common 0.5 to 1.5 ms from unclaimed events of score 0.2 to 0.6 as 20 to
50 ms from them, and 0.38 to 1.97 times for those of 0.7 or more
(`results/refractory_by_cosine.csv`).

A sorted spike 1 ms after the event lies inside the event's 3 ms trace, off
centre, and lowers its score.

## Subtracting the neighbouring spikes

`overlaps.py` takes every unclaimed event with a sorted spike of its tetrode
0.5 to 1.5 ms before or after it, and as a control up to 2,000 such reference
events per tetrode and window. From the event's 3 ms trace it subtracts
every sorted spike of the tetrode within 3.5 ms of the trace, as its unit's
mean waveform in that 1800 s window scaled to the spike (a reference event's
claiming spike stays). It then measures the event's trough again and scores
it, before and after the subtraction, against the mean waveforms of the
tetrode's units with 100 or more spikes in the window
(`results/overlaps.csv`). Those mean waveforms come from a new random draw
of spikes, so the scores before the subtraction differ from `gate.py`'s,
by 0.03 or less for 98 % of the events in the window checked. 171,316
of the 358,697 unclaimed events below the cutoff losing 10 % have such a
neighbour.

The neighbouring spikes' mean waveforms do not account for the events. After
the subtraction, 91 to 100 % of the unclaimed events in each group still have
a trough of 7 or more times the noise (reference events 96 to 100 %), and at
most 0.27 % fall below 5.5. At the event's trough the subtracted waveforms
point upwards, by a median of 0.13 to 1.14 times the noise. A feature that
only some of a unit's spikes carry would not show in its mean waveform, so
this does not rule that out.

Once the neighbours are gone, the events score higher. For events with a
sorted spike after them only, per group the number of events (summed over
TT5 to TT12), the median score per tetrode before and after the
subtraction, and the share whose score after it reaches the tetrode's
cutoff losing 10 % (both ranges over TT5 to TT12):

| Unclaimed events | Number | Median score, before and after | Reaching the cutoff losing 10 % |
|---|---|---|---|
| Below the cutoff losing 1 % | 24,499 | 0.27 to 0.40, 0.39 to 0.66 | 11.5 to 52.1 % |
| 1 to 5 % | 46,486 | 0.39 to 0.53, 0.49 to 0.76 | 54.5 to 91.9 % |
| 5 to 10 % | 32,467 | 0.46 to 0.62, 0.55 to 0.80 | 91.4 to 99.0 % |
| Reference events, control | 609,825 | 0.56 to 0.66, 0.65 to 0.81 | 87.7 to 94.0 % |

Of the unclaimed events with a sorted spike before them only, 5.5 to 32.8 %
(0 of 6 on TT9), 37.3 to 69.3 % and 77.5 to 90.2 % (5 of 7 on TT9) of the
three groups reach the cutoff.

Scoring the unclaimed events after the subtraction, and requiring a trough
of 7 or more after it, would admit 112,821 events now below the cutoff
losing 10 % and drop 5,278 now at or above it, 11.3 % more than the 947,625
it keeps. At the share losing 5 % it would admit 82,024 and drop 6,520, or
6.9 % more (`results/overlaps_rescored.csv`). Those counts use the cutoffs
from unsubtracted reference events and leave out events whose nearest
sorted spike is 1.5 to 3.5 ms away, so they are approximate. Reference
events with a neighbour also score higher after the subtraction, so cutoffs
from subtracted reference events would sit a little higher and the counts
admitted are probably a little too high (conjecture).

The unit that best matches the event after the subtraction is the
neighbouring spike's unit more often than if the two were unrelated: for
events below the cutoff losing 10 % with a spike after them only, 0.23 to
0.72 of the time against 0.18 to 0.57 expected, with the largest
differences on TT6 (0.56 against 0.28) and TT8 (0.72 against 0.57). The expected share is the sum
over units of the unit's share of best matches times its share of
neighbouring spikes.

On TT8, the 7,109 unclaimed events at or above the cutoff losing 10 % with a
sorted spike before them only keep a trough of 7 or more in 99.6 % of cases
after the subtraction and score a median 0.80. Their best unit is the
preceding spike's unit 0.71 of the time, against 0.58 expected. These are
the events behind TT8's excess of sorted spikes before them in the lags
table above (6.0 times, most often 1.43 ms before).

So most unclaimed events with a sorted spike about 1 ms after them look like
separate spikes overlapping a sorted one, scored low because the sorted
spike sits in their trace, and the excess of sorted spikes after
low-scoring events reflects which events the overlap pushes below the
cutoff. It does not show the events to be parts of the following spike's
waveform, which admitting them would count twice (conjecture). Where the
event matches the neighbour's own unit, one cell would have fired twice
within 0.5 to 1.5 ms, inside its refractory period, if both came from it; a
second cell with a similar waveform would also fit (conjecture).

## Scoring every event after the subtraction

`gate.py` also scores every event and background stretch after the
subtraction: each unit's template is extended to 3.5 ms either side from
the same spikes, units with fewer than 100 spikes in the window get one
from all their spikes for the subtraction only, and a spike within 0.5 ms of
the event is not subtracted. The scores before the subtraction are
identical to those of the earlier run. `summarize.py`, `refractory.py`,
`lags.py` and `figure.py` take `--sub` to use the scores after the
subtraction and write files ending in `_sub`. `gate.py` took 20 min over
the 96 windows with 32 workers.

The reference events' scores rise and narrow too, so every cutoff rises
(`results/tetrodes_sub.csv`). The residual amplitude is the trough after
the subtraction, as in `overlaps.py`:

| Tetrode | Losing 1 % | Losing 5 % | Losing 10 % | Losing 5 %, per 6 h block | Reference below 0.55 | Multi-unit rate at 5 % | Same, residual of 7 or more |
|---|---|---|---|---|---|---|---|
| TT5 | 0.35 | 0.46 | 0.52 | 0.40 to 0.52 | 14.8 % | 0.92 Hz | 0.85 Hz |
| TT6 | 0.54 | 0.65 | 0.70 | 0.63 to 0.67 | 1.2 % | 0.48 Hz | 0.46 Hz |
| TT7 | 0.47 | 0.58 | 0.63 | 0.54 to 0.61 | 3.2 % | 0.18 Hz | 0.17 Hz |
| TT8 | 0.52 | 0.63 | 0.68 | 0.61 to 0.64 | 1.6 % | 0.39 Hz | 0.38 Hz |
| TT9 | 0.44 | 0.52 | 0.56 | 0.51 to 0.54 | 8.8 % | 0.02 Hz | 0.02 Hz |
| TT10 | 0.37 | 0.46 | 0.51 | 0.45 to 0.47 | 17.2 % | 3.12 Hz | 2.99 Hz |
| TT11 | 0.42 | 0.53 | 0.58 | 0.52 to 0.54 | 6.9 % | 0.54 Hz | 0.50 Hz |
| TT12 | 0.48 | 0.56 | 0.61 | 0.54 to 0.58 | 4.0 % | 0.80 Hz | 0.75 Hz |

Over TT5 to TT12 the cutoff losing 5 % keeps 1,112,413 unclaimed events,
2.2 % more than the same share keeps without the subtraction (1,088,839).
Per tetrode the rate changes by -1.5 % (TT7) to +17.8 % (TT8), and +5.2 % on
TT6. Requiring a residual of 7 or more as well, as the estimate above did,
keeps 1,056,089, 1.12 % of the sorted spikes and 3.0 % fewer than without
the subtraction, where that estimate, with the cutoffs held fixed, gave
6.9 % more. Background stretches pass the cutoff losing 5 % in 0.2 to
3.8 % of cases.

After the subtraction, the reference events still keep their units'
refractory period from a score of 0.2 up: the ratio of the unit's spikes
0.5 to 1.5 ms away to those 20 to 50 ms away is at most 0.11 there, and at
most 0.46 below 0.2 (TT7, 60 events). Against the unit that an unclaimed
event matches best, that ratio now rises with the score. It is 1.04 to 2.78
at scores of 0.5 to 0.6 and 1.79 to 4.71 at 0.7 or more; without the
subtraction it was 0.38 to 1.97 at 0.7 or more
(`results/refractory_by_cosine_sub.csv`).

The timing changes too (`results/lags_sub.csv`). Sorted spikes follow
unclaimed events of residual 7 or more at or above the cutoff losing 10 %
1.23 to 3.85 times as often 0.5 to 1.5 ms later as 20 to 50 ms later,
against the reference events' 0.48 to 1.28, and precede them 1.00 to 3.02
times as often. After, the highest are TT8 (3.85, peaking at 1.40 ms) and
TT6 (2.14); before, TT8 (3.02, at 1.43 ms) and TT10 (1.65). Between the
cutoffs losing 5 and 10 % the ratio after is 1.12 to 2.52, and below the
cutoff losing 1 % 1.11 to 2.79. At or above the cutoff losing 10 %, 3.9 to
52.3 % of the events have a sorted spike 0.5 to 1.5 ms away, against 2.0 to
26.6 % without the subtraction: with the overlap no longer lowering the
score, events next to a sorted spike now score high (conjecture).

`neighbours.py` checks, in every 8th window, whether the sorted spikes 0.5
to 1.5 ms after the admitted events (losing 5 %, residual 7 or more) are
spikes at their own time rather than the event's own spike placed late
(`results/neighbours.csv`). Per tetrode, and separately for events that do
and do not match the neighbour's unit best, their trough is a median 0.90 to
1.20 of their unit's template trough, against 0.85 to 1.06 for the sorted
spikes after reference events (leaving out TT9, with 11 and 33 events). The
admitted events match the following spike's unit best 0.25 to 0.74 of the
time, against 0.18 to 0.60 expected if the two were unrelated (the sum over
units of the unit's share of best matches times its share of following
spikes), most on TT6 (0.53 against 0.29) and TT8 (0.74 against 0.60).
Reference events do so 0.04 to 0.27 of the time, against 0.15 to 0.43
expected (TT9 left out throughout). Two cells with similar waveforms firing
about 1 ms apart would fit, as would one cell firing inside its refractory
period (conjecture).

Unclaimed events whose residual is below 7 are 5.9 % of the unclaimed
events (2.3 to 9.2 % per tetrode; reference events 1.6 to 5.1 %). Sorted
spikes 1.5 to 3 ms before them are 4.26 to 9.00 times as common as 20 to
50 ms away (31.1 on TT9), and 1.5 to 3 ms after them 2.70 to 6.44 times
(17.9 on TT9), peaking 1.83 to 2.53 ms away; for reference events the same
ratios are 1.00 to 1.44 and 0.93 to 1.41 (2.92 and 2.31 on TT9)
(`results/lag_hist_sub.csv`). The bandpass of step 2 (Butterworth of order
5, 300 to 6000 Hz, forward and backward) turns a negative pulse 1 to 6
samples wide into a trough with negative lobes 2.33 ms before and after it,
0.9 to 2.0 % of the trough, so a sorted spike's echo would deepen the trace
about 2 ms away. These
events may be crossings that a sorted spike's filter echo lifted over 7
(conjecture).

## Figures

`figure.py` draws, from the per-window files and `results/`, figures not
kept here:

- `cosine_distributions.png`: per tetrode, the score distributions of
  reference events, unclaimed events and background stretches, with the
  cutoffs losing 1, 5 and 10 % and 0.55.
- `lags.png`: the lag curves behind the table above, from -3 to 3 ms.
- With `--sub`, `figure.py` draws the same figures from the scores after
  the subtraction (`cosine_distributions_sub.png`, `lags_sub.png`,
  `gallery_sub_TTn.png`); the galleries then show the residual traces, and
  `lags_sub.png` adds the unclaimed events whose residual is below 7.
- `overlaps.png`, drawn by `overlaps.py` into its output directory: four
  random unclaimed events per tetrode below the cutoff losing 10 % with a
  sorted spike 0.5 to 1.5 ms after them, with the subtracted spikes and the
  residual.
- `gallery_TTn.png`, one per tetrode of TT5 to TT12: events by band of
  score (below 0.2, then 0.1 wide to 0.7, then 0.7 or more). Each column has
  four unclaimed events of 7 to 9 times the noise, two unclaimed of 9 or
  more and three reference events of 7 to 9, drawn at random from a sample
  of up to 3 per window, 0.05 band of score (of the score after the
  subtraction, since `gate.py` computes it) and amplitude class. Each trace
  is drawn with the best template scaled to it by least squares and ticks
  at the tetrode's other sorted spikes.

In the galleries, below a score of about 0.3 the unclaimed and reference
events alike sit in fast fluctuations of the trace and rarely resemble a
template, and many of the unclaimed ones have a large sorted spike later in
the trace. From about 0.5 most unclaimed events show a trough of the
template's shape on the template's wires. Several unclaimed events below
0.3, of 9 or more and on TT12 also of 7 to 9, are sharp transients of up to
several hundred µV, some seen on 10 or more other tetrodes (conjecture from
the galleries, not yet confirmed).

## Fixed noise against per-window noise

Over TT5 to TT12, fixing each channel's noise for the whole recording leaves
the number of unclaimed events of 7 or more almost unchanged (1,306,322
against 1,308,865 in the replay). Per tetrode it moves them by -14.5 % (TT9)
to +30.3 % (TT6), and +18.8 % on TT8. The events scoring 0.55 to 0.8 change
by +3.0 % in total and by -11.2 to +34.3 % per tetrode
(`results/tetrodes.csv`, against
`investigations/unclaimed-events/results/tetrodes.csv`); those counts also
differ by the random spikes drawn for the templates.

## Outcome

The gate scores each event after subtracting the tetrode's sorted spikes
around it, with one cutoff per tetrode: the score that 5 % of the reference
events, scored the same way, fall below (`docs/sorting_pipeline.md`,
step 8).

## Open

- Whether the 7 MAD floor also applies to the residual trough.
- What the events on TT6 and TT8 that match the unit of a sorted spike 0.5
  to 1.5 ms away are.
