---
title: Unclaimed events in the final TTM-003-48H sorting
updated: 2026-11-05
---

# Unclaimed events in the final TTM-003-48H sorting

What the two steps of `investigations/spike-sorting/97_production_pipeline.py`
(97 below) that follow its merge would do to the final TTM-003-48H sorting
(`sortings_carry/main` in the session directory, 61 units, 138.7 M
spikes). The question was asked for the analyzer stage
(`docs/sorting_pipeline.md`, step 8).

Both steps work on events, detected negative peaks, with amplitudes in
multiples of the channel's MAD noise. An event is claimed when a sorted
unit of its tetrode fires within 0.5 ms of it.

- The fold takes unclaimed events of 10 MAD or more. It adds each to the
  unit whose template it matches best, when the cosine is 0.8 or more and
  that unit has no spike within 1.5 ms of it.
- The multi-unit step takes the events of 7 MAD or more still unclaimed
  after the fold. It pools those with a weaker match (cosine 0.55 to 0.8)
  into one extra unit per tetrode.

Code: `replay.py`, `summarize.py`, `figure.py`. Results in `results/`.

## Replay

`replay.py` repeats both steps on the sort's int16 binary, one 1800 s
window at a time:

- Events are negative peaks at 5.5 times each channel's MAD noise
  (SpikeInterface's "locally_exclusive" detection, radius 40 µm), so a
  trough gives one event per tetrode.
- Templates are the mean of up to 500 random spikes of each unit with 100
  or more spikes in the window, 1 ms before to 2 ms after the spike.
- The cosine is taken over that 3 ms on the tetrode's wires, with no shift.
- An event within 0.3 ms of the previous one of the same train is
  dropped, as in 97. The train is the unit's for the fold and the
  tetrode's multi-unit train for the multi-unit step.

It departs from 97 in these ways:

- Detection. 97 took its events from one detection over the whole
  recording, with one SpikeInterface noise estimate (20 random chunks from
  anywhere). The replay detects in each window. Its noise there is the
  mean MAD of 40 random 0.5 s chunks of the window outside the tetrode's
  bad intervals. The 7 and 10 MAD floors use the same noise.
- Bad intervals. Events inside a bad interval of their tetrode are left
  out, since the sort has no spikes there.
- Templates. Both steps use the same templates. 97 rebuilt them after
  the fold, with the folded spikes added.
- Smaller differences. 97's refractory check sees only the unit's spikes
  in the same window, and it drops close multi-unit events over the whole
  recording; the replay does the first over the whole recording and the
  second within each window.

Events in the first 30 or last 60 samples of a window are also left out
of the counts. The sort has no spikes there, and 97 scores no event there
either.

The run took 70 to 93 s per window, 12 minutes with 12 workers.

For each unclaimed event of 7 MAD or more, `replay.py` also records the
following:

- the matcher's highest-snr fit on the tetrode within 0.5 ms, and the
  highest-snr fit of the best-matching unit there, each with the reason
  carry-forward matching rejected it (`docs/sorting_pipeline.md`, step 5);
- how many other tetrodes have an event within 0.5 ms, counting only
  events not left out above (for all events, this is kept as counts);
- whether it falls inside a laser pulse (onset to offset plus 100 ms).

The spikes each template is built from (up to 500 per unit and window,
"the units' own spikes" below) get the same cosines, amplitude and scale.
They are scored at the unit's spike times, and their amplitude is the
lowest sample within 5 samples of the spike time. Events are scored at
their own trough, with the amplitude that detection measured.

## How many events no unit claims

On TT5 to TT12, outside bad intervals:

| Tetrode | Events of 7 MAD or more | Unclaimed, 7 MAD or more | Unclaimed, 10 MAD or more | Folded | Multi-unit events | Multi-unit rate | Sorted rate |
|---|---|---|---|---|---|---|---|
| TT5 | 5.03 M | 3.4 % | 0.4 % | 732 | 115,079 | 0.67 Hz | 84.3 Hz |
| TT6 | 5.82 M | 1.3 % | 0.2 % | 700 | 39,876 | 0.23 Hz | 108.5 Hz |
| TT7 | 0.92 M | 4.6 % | 1.1 % | 202 | 27,282 | 0.16 Hz | 32.9 Hz |
| TT8 | 3.94 M | 2.1 % | 0.5 % | 426 | 47,065 | 0.27 Hz | 71.3 Hz |
| TT9 | 41,296 | 12.7 % | 29.4 % (149 of 507) | 0 | 3,380 | 0.02 Hz | 8.0 Hz |
| TT10 | 5.34 M | 12.6 % | 3.8 % | 639 | 332,694 | 1.93 Hz | 95.1 Hz |
| TT11 | 2.36 M | 4.2 % | 0.8 % | 610 | 67,933 | 0.39 Hz | 68.4 Hz |
| TT12 | 5.00 M | 3.3 % | 0.3 % | 316 | 129,626 | 0.75 Hz | 76.1 Hz |

Rates are over each tetrode's time outside its bad intervals; the sorted
rate is all its units together. On TT9 no unclaimed event of 10 MAD or
more reaches cosine 0.8, so nothing is folded there.

Five other tetrodes have units with templates: TT1, TT2, TT4, TT15 and
TT16. On these, 1.3 % (TT4) to 16.1 % (TT16) of events of 7 MAD or more
are unclaimed. TT3 and TT13 each have one unit with too few spikes for a
template in any window, so 98.7 % and 93.2 % of their few events (2,048
and 2,008) are unclaimed and none can be folded or pooled. TT14 has no
units.

Laser pulses cover 0.25 % of the recording. On each tetrode 0 to 0.33 %
of the folded and of the multi-unit events fall inside them.

## The fold

Over the 48 h the fold proposes 4,399 events and rejects 374 because the
unit has a spike within 1.5 ms. That leaves 4,025 folded (3,625 on TT5 to
TT12). No unit gains more than 0.019 % of its spikes; the most is unit 10
on TT5, with 437 events.

The folded events are larger than the units' typical spikes. Their scale
relative to the template (projection on the template over its squared
norm) has a median of 1.84 (10th to 90th percentile 1.42 to 2.50). Take
the 25 units on TT5 to TT12 with 20 or more folded events. A unit's
median folded event is 1.5 to 2.1 times its template. Only 0.6 to 8.2 %
of its own spikes (median 2.5 %) are at least that large.

`figure.py` draws the mean waveforms from the replay's per-window files;
the figure is not kept here. For each of the four units with the most
folded events, the folded events' mean keeps the template's spread over
the four wires. Its trough is about twice as deep and narrower. Some of
the narrowing may come from alignment: the events are aligned on their
own trough, the template on the unit's spike times.

The best template barely beats the next. Over the folded events, the best
template's cosine exceeds the second best by a median of 0.041. With the
best of shifts up to 3 samples, 15 % of the folded events would go to a
different unit, and the fold would propose 6,404 events instead of 4,399.

Whether a spike reaches cosine 0.8 depends on its amplitude. Among the
units' own spikes on TT5 to TT12, these shares have cosine 0.8 or more
with their own template:

| Amplitude | Cosine 0.8 or more |
|---|---|
| under 5.5 MAD | 2.8 % |
| 5.5 to 7 MAD | 11.7 % |
| 7 to 10 MAD | 25.5 % |
| 10 to 15 MAD | 45.5 % |
| 15 MAD or more | 60.9 % |

Applied to the sort's own spikes, the fold's rule (cosine 0.8 or more
with the unit's own template, and that template the best match) accepts
0 to 33 % (median 8.6 %) of each TT5 to TT12 unit's spikes.

Carry-forward matching fitted the unit within 0.5 ms of 91 % of the
folded events and rejected the fit. Of all folded events, 51 % had a fit
whose neighbour-subtracted snippet had snr under 4.5, 34 % one whose
snippet had cosine under 0.5 with the template, and 6 % one that lost to
the refractory rule.

The folded events' timing around the unit's spikes was compared with the
unit's own firing. If they were the unit's spikes, and the unit's largest
spikes fired like the rest, their cross-correlogram with the unit beyond
1.5 ms would look like its autocorrelogram. That assumption may not hold:
spike size can depend on recent firing. The 0.5 to 1.5 ms range is empty
by the fold's rule. The table gives the mean count per lag 1.5 to 5 ms
away, over the mean 20 to 50 ms away, for the 6 units with 200 or more
folded events:

| Unit | Folded | Folded with the unit | The unit's autocorrelogram | Folded with a same-tetrode neighbour, largest |
|---|---|---|---|---|
| 10 (TT5) | 437 | 1.98 | 1.54 | 2.45 |
| 13 (TT6) | 301 | 1.17 | 0.85 | 1.26 |
| 21 (TT8) | 210 | 1.17 | 0.98 | 1.12 |
| 22 (TT8) | 216 | 1.65 | 0.95 | 1.41 |
| 26 (TT10) | 431 | 1.36 | 1.19 | 1.40 |
| 29 (TT11) | 220 | 1.63 | 1.33 | 1.54 |

For three of the six, the folded events are more common 1.5 to 5 ms from
a neighbour's spikes than from the unit's; for the other three, the
reverse. In all six they are more common near the unit's spikes than the
unit's own spikes are. The last column is the largest over the unit's
neighbours (one on TT8, up to six elsewhere), so where there are several
it leans high.

On TT1, TT2, TT15 and TT16, 68 to 100 % of the folded events coincide
with events on 3 or more other tetrodes, against 0.5 to 3.8 % on TT5 to
TT12. On those four tetrodes the claimed events do so too, 31 to 72 %.

## The multi-unit step

The multi-unit units would hold 832,508 events (762,935 on TT5 to TT12),
0.60 % as many as the sorted units. Per tetrode that is 0.21 % (TT6) to
2.0 % (TT10) of the sorted rate on TT5 to TT12, and under 1.4 % on the
others.

Over all tetrodes the events are small, near the 7 MAD floor:

- Their median amplitude is 7.63 MAD (90th percentile 9.2).
- Their scale relative to the best template has a median of 1.16 (10th
  to 90th percentile 0.86 to 1.56).
- Their mean waveform on each of TT5 to TT12 looks, by eye, like the
  tetrode's templates.

The cosine band does not separate them from sorted spikes. On TT5 to
TT12, 37 to 91 % of the units' own spikes of 7 MAD or more have their best
cosine in the band, 0.55 to 0.8. The rest split into 4 to 58 % at 0.8 or
more and 3 to 24 % under 0.55. Of the unclaimed events of 7 MAD or more
there, 50 to 79 % fall in the band. The two are scored at different
alignments (above), so the comparison is approximate.

Over all tetrodes, matching fitted something within 0.5 ms of 99.7 % of
the multi-unit events and rejected the strongest fit:

- 50 % for shape;
- 37 % for snr under 4.5;
- 12 % for the refractory rule.

On TT5 to TT12, 1.8 to 4.9 % of the multi-unit events coincide with
events on 3 or more other tetrodes. Claimed events of the same tetrodes
do so in 1.4 to 5.6 %. Tetrode by tetrode, the multi-unit share is higher
on 7 of the 8, by 1.04 to 1.42 times, and lower on TT9 (0.88). On TT1,
TT2, TT15 and TT16 the multi-unit share is 47 to 68 %. An event on 4 or
more tetrodes at once may be a signal shared across the probe rather than
a spike near one tetrode (conjecture).

There are 36 units on TT5 to TT12 that are the best match of 200 or more
multi-unit events. For each, take the events' rate 0.5 to 1.5 ms from the
unit's spikes, over their rate 20 to 50 ms away. The units' own
autocorrelograms are at 0.01 to 0.05 there. Pairs of sorted units on one
tetrode of TT5 to TT12 are at 1.09 to 1.91 (median 1.33, 76 pairs) in the
same terms.

- 17 units are below 1 (0.31 to 0.99), under every sorted pair.
- 10 are at 1.12 to 1.89, within the pairs' range.
- 9 are at 2.06 to 5.69, above every pair: four of the five TT6 units
  (12, 13, 15, 16), unit 18 on TT7, unit 22 on TT8, and units 30, 32 and
  33 on TT11.

## What the numbers say

The fold would change no unit by more than 0.019 %. The events it would
add are 1.5 to 2.1 times the unit's template (median 1.84), larger than
all but 0.6 to 8.2 % of the unit's own spikes. Their best template
barely beats the next, and their timing ties them to the unit no more
consistently than to a neighbour. Nothing measured here shows that they
belong to the unit they would join (conjecture).

The multi-unit units would hold 0.60 % as many events as the sorted
units. Their cosine band holds 37 to 91 % of the sorted units' own
spikes of 7 MAD or more. Around the unit they match best, the events of
some units are rarer 0.5 to 1.5 ms away than 20 to 50 ms away, though
far more common there than the unit's own spikes. The events of other
units are more common there than for any sorted pair. So the bucket
looks like a mixture rather than one kind of event (conjecture).

The fold is not adopted (`docs/sorting_pipeline.md`, step 8).
Whether the multi-unit step is, is open.
