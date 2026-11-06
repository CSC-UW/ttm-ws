---
title: Stimulation protocol for the first and last hour of a sleep deprivation
updated: 2026-11-06
status: needs_verification
---

# Stimulation protocol for the first and last hour of a sleep deprivation

This protocol asks whether laser pulses become more or less effective at
producing OFF periods and slow LFP waves (defined under Analysis)
between the first and the last hour of a 5 h sleep deprivation by novel
objects. Light is delivered only in those two hours,
so the three hours between them are not stimulated. The comparison is
each outcome's mean over the last hour minus its mean over the first,
with no assumption about the time course in between.

The precision figures come from TTM-003-48H-2
(`docs/reports/2026-11-05_ttm003_48h2_pulse_response.md`,
`subjects/ttm-003/results/sd_protocol_power_TTM-003-48H-2.json`) and
assume the same subject, light and implant. They are estimates, not
measurements of this protocol.

## Light

Each pulse repeats the TTM-003-48H-2 pulse. The power was about 0.3 mW.
The laser driver's level setting, which the hand log of 2026-10-16
records for those sessions, was not recorded for TTM-003-48H-2. Measure
the power at the fibre tip before and after the session. The width was
nominally 50 ms, and the acquisition board recorded 63.9 ms (median of
2,900 pulses, 59.9 to 78.1 ms), as it recorded about 64 ms for the
nominal 50 ms pulses of 2026-10-16 (pattern B in
`subjects/ttm-003/optogenetic_stimulation.md`). Which device setting
produced those pulses is not written down; confirm it before the
session. If the device that plays the schedule also drives the laser
TTL, a width that records about 64 ms on line 1 reproduces them.

## Windows

The sleep deprivation runs from 09:00 to 14:00 (lights on at 08:53).
Window 1 is 09:00:00.000 to 10:00:00.000 and window 2 is 13:00:00.000 to
14:00:00.000, local time. Each window's time zero is the stimulation
device's clock reading of its start; set that clock from the same time
source as the acquisition computer before the session. An offset of a
second between the two clocks does not matter, because the analysis uses
the event times the acquisition board records. No event is earlier than
5,000 ms or later than 3,595,000 ms after time zero (both bounds
inclusive).

The device plays the schedule's `offset_ms` from each window's time zero;
the schedule's `clock` column is for reference. If handling at 09:00
disturbs the first minutes, window 1 can start a few minutes late and
keep its 1 h length, and window 2 stays at 13:00 to 14:00. Generating the
schedule again with the same seed and the shifted window
(`--windows 09:05-10:05 13:00-14:00`) changes only the `clock` column.

## Events

Each window holds E events: R light pulses and S shams. A sham is a
scheduled event at which no light is delivered. It is marked on its own
TTL line and analysed like a pulse, and gives the rate of spontaneous OFF
periods under the same timing. The pulse rate is open; I suggest R = 480
with S = 480.

| Pulses per hour (R) | Shams per hour (S) | Events (E) | Median interval between events | 1st to 99th percentile | Median interval between pulses |
|---|---|---|---|---|---|
| 240 | 240 | 480 | 5.45 s | 1.07 to 31.63 s | 10.34 s |
| 480 | 480 | 960 | 2.90 s | 1.03 to 13.86 s | 5.35 s |
| 960 | 960 | 1,920 | 1.60 s | 1.01 to 4.99 s | 2.82 s |
| 960 | 240 | 1,200 | 2.37 s | 1.02 to 10.04 s | 2.89 s |

The intervals are from 20 draws of the rule below with a 1 s minimum
interval.

## Timing rule

Times are in milliseconds from the window's time zero. An event's time is
the rising edge of its TTL, on the laser line for a pulse and on the sham
line for a sham. The window length is W = 3,600,000, the margin
m = 5,000 and the minimum interval, onset to onset, d = 1,000.

1. The free time is F = W - 2m - (E - 1)d.
2. Draw E independent values uniformly on [0, F) and sort them, so that
   u(0) <= u(1) <= ... <= u(E - 1).
3. Event j, for j = 0 to E - 1, is at t(j) = m + u(j) + j d, rounded to
   the nearest millisecond.
4. Choose exactly S of the E events uniformly at random, without
   replacement. These are the shams and the rest are pulses.
5. Use the same times and the same labels in both windows.
6. Use a new random seed for every sleep deprivation.

What the rule gives:

- Every interval between events is at least 1 s. In TTM-003-48H-2,
  pulses 1 s apart showed no decline of the OFF share (defined under
  Expected precision) over 20 pulses: 0.363 for pulses 1 to 5 and 0.345
  for pulses 16 to 20 over all trains, a change in log-odds per pulse
  position of -0.004 (95 % interval -0.017 to 0.009). Three quarters of
  OFF periods ended by 126 ms after onset, and the mean LFP response was
  back at baseline by 338 ms, so the mean response to a pulse is over
  before the 500 ms that precede the next event, the longest stretch
  before onset that the analysis uses. The 1 s spacing there lasted for
  20 pulses; here intervals near 1 s are scattered through the hour
  (31 % of them are under 2 s at R = S = 480).
- The part of each interval beyond 1 s is close to exponential, with mean
  about F / E (2.74 s at R = S = 480). The chance of the next event does not
  grow while waiting. My expectation, untested, is that this leaves no
  rhythm to entrain to or anticipate.
- The counts are exact and the events are spread evenly across the hour
  on average.
- Because shams are drawn at random among all events, a pulse and a sham
  are almost equally likely to follow a sham (480 and 479 in 959 at
  R = S = 480), so the time since the previous light pulse has nearly the
  same distribution for both.
- With the same schedule in both windows, the windows differ only in when
  they occur.

`subjects/ttm-003/scripts/sd_stim_schedule.py` implements the rule. Its
`--sham-every k` makes one event in k a sham, so S = R / (k - 1), and R
must be a multiple of k - 1: k = 2 gives S = R and k = 5 gives the row
with 960 pulses and 240 shams. Its `--min-interval-s` sets d in seconds
(default 1). It writes a CSV with the columns `window` (from 1), `index`
(from 0 in each window), `offset_ms` (the time from time zero), `type`
(`pulse` or `sham`), `clock` (local clock time to the millisecond) and
`interval_ms` (from the previous event in the window, blank on each
window's first row), and a JSON file with the parameters, the seed and its checks: the
counts, the shortest, mean and longest interval, the margins, and the
median time since the previous light pulse for pulses and for shams
(5.15 and 5.42 s in the example below; over 300 seeds these medians
average 5.44 and 5.46 s, and their difference varies by 0.41 s, one
standard deviation, from seed to seed).

```
python subjects/ttm-003/scripts/sd_stim_schedule.py --pulses-per-hour 480 \
    --sham-every 2 --date 2026-12-01 --seed 20261201 --out sched_2026-12-01
```

The CSV is the schedule. Another implementation of the rule will not
reproduce a seed's draw, so the stimulation device should play the CSV's
`offset_ms` and `type` rather than draw its own times, and the CSV is
fixed and saved before the session.

## Outputs and recording

At each pulse, the device raises the laser line for the pulse width. The
laser line goes to the laser driver and to acquisition-board TTL line 1,
the input the laser TTL used in both 48 h recordings of TTM-003
(TTM-003-48H, at higher power, and TTM-003-48H-2). Lines are numbered as
in `subjects/ttm-003/optogenetic_stimulation.md`.

At each sham, the device raises a separate sham line for 10 ms. The sham
line goes only to acquisition-board TTL line 3 and never to the laser.
Its 10 ms width differs from the pulse width, so the two cannot be
confused. Line 2 stays unused, because it has carried sub-millisecond
crosstalk from line 1 (`subjects/ttm-003/optogenetic_stimulation.md`).

For the whole of each window, the device holds a window line high on
acquisition-board TTL line 4.

Save with the recording: the schedule's CSV and JSON (the JSON holds the
seed), the device's own record of what it played if it keeps one, the
device clock time at which each window started, the laser level, the
power at the fibre tip before and after, and any interruption.

## Checks

Before the session, with the laser disconnected, play the first 5 min of
the schedule into the acquisition board and confirm that the counts,
widths and order on lines 1 and 3 match the CSV, that line 4 is high
throughout, and that line 1 stays low at sham times. Fix any mismatch
before the sleep deprivation starts.

After the session, compare the events recorded on lines 1 and 3 with the
CSV: the counts, the order, and every interval between consecutive events
of lines 1 and 3 taken together, within 1 ms plus 50 parts per million of
its length (the two clocks may run at slightly different
rates). Log any missing or extra event. The analysis uses the recorded
times, so a timing deviation is tolerable; a missing or extra light pulse
must be known.

## Expected precision

The outcome with the best-understood noise is the OFF share, the
fraction of pulses followed by an OFF period. Threshold crossings are
negative crossings of 5 times the median absolute deviation on each wire
(300 to 6,000 Hz, with the median of the channels not marked bad in the
experiment's manifest subtracted), merged per
tetrode within 1 ms. A tetrode is suppressed when its crossing rate from
10 ms after onset to the offset is under half its rate in the 500 ms
before onset (in TTM-003-48H-2, TT5 to TT9 and TT11). An OFF period is a
gap of 50 ms or more in the suppressed tetrodes' pooled crossings that
spans onset + 20 ms.

In TTM-003-48H-2 the pulses came in trains of 20 at 1 Hz, one train every
20 min for 48 h. In the 16 trains of the sleep deprivation the OFF share
was 0.33 and the per-pulse variance 0.221. The OFF probability moves with
brain state: the variance of the train OFF share beyond its binomial part
was 0.010 there. The covariance of the OFF outcome between pulses of one
train, with the means of each light phase, dark phase and the sleep
deprivation, a linear trend in each and the mean of each pulse position
removed, was:

| Pulses apart | All 145 trains | 95 % interval | 16 sleep-deprivation trains | 95 % interval |
|---|---|---|---|---|
| 1 to 3 s | 0.017 | 0.010 to 0.025 | 0.012 | -0.003 to 0.029 |
| 4 to 7 s | 0.009 | 0.003 to 0.015 | 0.008 | -0.010 to 0.027 |
| 8 to 12 s | 0.008 | 0.001 to 0.015 | 0.015 | -0.002 to 0.035 |
| 13 to 19 s | 0.012 | 0.003 to 0.021 | -0.012 | -0.040 to 0.009 |

Between trains 20, 40, 60 and 80 min and 2, 3, 4, 5 and 6 h apart (each
within 10 min), the covariance of the train OFF shares (same means and
trends removed) lies between -0.0050 and +0.0007,
mean -0.0018, with standard errors of 0.0023 to 0.0030; only the 80 min
lag (-0.0050, standard error 0.0028) is more than one and a half
standard errors from zero, and it is negative. Fitted against these
covariances, a state component of variance 0.005 lasting about an hour
(exponential time constant 3,600 s) gives a chi-square of 14.7 over the 9
lags, against 7.4 with no such component; one of 0.010 lasting about
20 min (1,200 s) gives 11.6; and one lasting about 7 min (400 s) gives at
most 7.6, because its covariance has vanished by 20 min. The lags share
trains, so these figures are a guide rather than a test: the data
disfavour an hour-long component, weakly disfavour a 20 min one, and
cannot see a 7 min one.

How much of the state variance survives averaging over an hour depends on
how slowly state moves, which pulses spaced at 1 s or 20 min do not show
between 20 s and 20 min. The spontaneous activity of the suppressed
tetrodes fills that gap. I took the multi-unit rate (the square root of
the crossings in each 1 s bin, with crossings found as above except that
the median of the suppressed tetrodes' channels is subtracted and the
median absolute deviation is taken per 30 s) and the delta (0.5 to 4 Hz)
and sigma (10 to 15 Hz) shares of the 0.5 to 30 Hz LFP power in 4 s
windows stepped by 1 s. Of their variance at periods of 5 to 1,024 s,
through the 5 h sleep deprivation with each train and the 2 s before and
12 s after it excluded:

| Measure | Above 120 s | 10 to 120 s | 5 to 10 s |
|---|---|---|---|
| Multi-unit rate | 4.9 % | 67.6 % | 27.4 % |
| Delta share | 16.1 % | 61.1 % | 22.7 % |
| Sigma share | 10.7 % | 60.8 % | 28.5 % |

Over session hours 2 to 7 (the first light phase, which includes sleep),
the delta share has 46.7 % of that variance above 120 s. My conjecture is
that this comes from the alternation of sleep and wake, which the sleep
deprivation suppresses.

As I read the literature, it agrees on the fast part and is thin on the
slow part. In
awake mice, pupil-linked cortical state switches within 1 to 2 s (Reimer
et al. 2014, doi:10.1016/j.neuron.2014.09.033), and single-unit rate
spectra in prefrontal cortex follow a power law from 0.01 to 1 Hz with no
peak (Okun et al. 2019, doi:10.1093/cercor/bhz023). The 0.02 to 0.03 Hz
rhythm of sigma power and noradrenaline is described in NREM sleep (Lecci
et al. 2017, doi:10.1126/sciadv.1602026; Osorio-Forero et al. 2021,
doi:10.1016/j.cub.2021.09.041; Kjaerby et al. 2022,
doi:10.1038/s41593-022-01102-9), not in wake. A rhythm near 0.01 Hz
completes about 36 cycles in an hour and barely moves an hourly mean.
What matters is power at periods of 10 min or more, and I found no rodent
study that measures it in wake. Ultradian rhythms of 3 to 5 h in mouse
and rat activity (Dowse et al. 2010, doi:10.1242/jeb.038877; Stephenson
et al. 2012, doi:10.1177/0748730412461247) are weak in the light phase in
mice (Blum et al. 2014, doi:10.7554/eLife.05105).

Four models of the state covariance of the OFF outcome span these data,
the last only marginally consistent with them. Each is a sum of
exponentials in the lag between pulses, with a fast part of variance
0.007 and time constant 3 s, plus:

- A: 0.010 with time constant 16 s, a spectrum that halves at 0.01 Hz and
  has nothing slower;
- B: 0.0085 at 16 s and 0.0015 at 400 s, the slow part about the share of
  the delta share's variance above 120 s;
- C: 0.0075 at 16 s and 0.0025 at 1,200 s, a quarter of the slower
  variance lasting about 20 min, more than the spontaneous spectra
  suggest;
- D: 0.010 at 1,200 s, all of the slower variance lasting about 20 min,
  the worst case the covariances between trains only weakly disfavour.

Their covariance between pulses 1 to 3, 4 to 7, 8 to 12 and 13 to 19 s
apart is 0.013 to 0.014, 0.008 to 0.011, 0.006 to 0.010 and 0.004 to
0.010. All lie inside the 95 % intervals of both columns of the table
above, except D at 13 to 19 s in the sleep deprivation (0.010 against an
upper bound of 0.009).

With a binomial variance of 0.205 per pulse (the sleep deprivation's
per-pulse variance, 0.221, less the models' state variance at zero lag,
0.017), the standard error of the last-minus-first difference follows
from a model and a schedule. I computed it with pulse times in 1 s bins,
averaged over five schedules, with the two hours independent. The table
gives the smallest difference in OFF share that one sleep deprivation
detects (80 % power, two-sided p < 0.05, 2.8 standard errors) and, in
brackets, the sleep deprivations needed to detect 0.10 and 0.05 if the
true change is the same in each (from the unrounded differences):

| Pulses + shams per hour | A | B | C | D |
|---|---|---|---|---|
| 240 + 240 | 0.126 (2, 7) | 0.143 (3, 9) | 0.183 (4, 14) | 0.294 (9, 35) |
| 480 + 480 | 0.093 (1, 4) | 0.115 (2, 6) | 0.162 (3, 11) | 0.281 (8, 32) |
| 960 + 960 | 0.071 (1, 3) | 0.098 (1, 4) | 0.151 (3, 10) | 0.275 (8, 31) |
| 960 + 240 | 0.071 (1, 3) | 0.098 (1, 4) | 0.150 (3, 10) | 0.274 (8, 31) |

These figures are conjecture: they carry TTM-003-48H-2's OFF share and
state variance, and the shape of its spontaneous spectra, over to a new
session.

In TTM-003-48H-2, sham onsets were 40 times at 1 s spacing before each
train, and 0.5 % of those in the sleep deprivation were followed by an OFF
period (9 % over the whole recording). Treating shams as independent
binomial events, with no state covariance among themselves or with the
pulses: at 0.5 %, subtracting the shams widens the detectable
pulse-minus-sham difference by at most 3.1 %. At 9 %
it would widen it by 2 to 16 % at R = S = 240 or 480, depending on the
model, by 1 to 12 % at R = S = 960, and by 4 to 43 % at 960 pulses with
240 shams.

For scale, in rats sleep deprived with novel objects from light onset,
OFF periods in waking cortex became about 57 % more frequent from the
first to the fourth hour (57.7 % in frontal and 56.6 % in parietal
cortex; Vyazovskiy et al. 2011, doi:10.1038/nature10009). If
pulse-evoked OFF periods rose in the same proportion, the OFF share would
go from 0.33 to about 0.52, a difference of 0.19. This is conjecture,
since spontaneous and evoked OFF periods need not change alike.

## Analysis

The primary outcome is the OFF share of pulses, last hour minus first
hour, with an interval that allows for the correlation between nearby
pulses (a moving-block bootstrap within each hour with 5 min blocks,
checked with 2 and 10 min blocks). The same comparison applies to the
excess over shams, to the end of the OFF period, and to the two slow LFP
waves of the report: the mean LFP (0.5 to 30 Hz) of the suppressed
tetrodes from 10 to 100 ms (early wave) and from 100 to 300 ms (late wave)
after onset, each minus its mean in the 200 ms before onset. The
covariance of pulse outcomes by lag within each hour measures how long
state persists, and shows which of the four models held.

## Not decided

- The pulse rate. I suggest 480 per hour: under A it detects 0.093
  against 0.071 at 960 per hour, and under B 0.115 against 0.098, with
  half the stimulation; under C and D the rate matters little.
- Window 2 follows the 480 pulses of window 1, so any lasting effect of
  the first hour's stimulation on the last hour's response is mixed with
  time in the sleep deprivation. A sleep deprivation stimulated only in
  the last hour would measure it.
- A baseline-day control at the same clock hours, which would separate
  time in the sleep deprivation from time of day and from the novelty of
  the objects. It doubles the sleep deprivations needed and stimulates on
  the baseline day.
