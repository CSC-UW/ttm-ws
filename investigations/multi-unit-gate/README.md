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

Code: `gate.py`, `summarize.py`, `refractory.py`, `lags.py`, `figure.py`.
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

A sorted spike 1 ms after the event would lie inside the event's 3 ms trace,
off centre, and lower its score. Whether such an event is a separate spike,
from another cell or missed by the sort in an overlap, or part of the
following spike's own waveform, is not settled here. If it is part of the
waveform, admitting it to the multi-unit unit counts that spike twice
(conjecture).

## Figures

`figure.py` draws, from the per-window files and `results/`, figures not
kept here:

- `cosine_distributions.png`: per tetrode, the score distributions of
  reference events, unclaimed events and background stretches, with the
  cutoffs losing 1 and 5 % and 0.55.
- `lags.png`: the lag curves behind the table above, from -3 to 3 ms.
- `gallery_TTn.png`, one per tetrode of TT5 to TT12: events by band of
  score (below 0.2, then 0.1 wide to 0.7, then 0.7 or more). Each column has
  four unclaimed events of 7 to 9 times the noise, two unclaimed of 9 or
  more and three reference events of 7 to 9, drawn at random from a sample
  of up to 3 per window, 0.05 band of score and amplitude class. Each trace
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

## Open

- The share of reference events the cutoff may lose. At or above the
  cutoff losing 10 %, the rate of sorted spikes 0.5 to 1.5 ms after an
  unclaimed event is within the reference events' range, though sorted
  spikes precede those events more often than they precede reference
  events. Between the 5 and 10 % cutoffs, the rate after is 1.96 to 3.50
  times the rate 20 to 50 ms away.
- Whether the cutoff should rise with amplitude, as it does for the
  reference events on six of the eight tetrodes.
- What the events about 1 ms before a sorted spike are.
