---
title: TTM-003-48H-2, response to the 0.3 mW pulses and the pulse rule
updated: 2026-11-06
status: needs_verification
---

# TTM-003-48H-2: response to the 0.3 mW pulses, and whether sorting needs a new pulse rule

TTM-003-48H-2 (session `2026-11-03_09-01-38`) carries 145 trains of 20
light pulses at 1 Hz, one train every 20 min for 48 h. The pulses are
nominally 50 ms (64 ms on the TTL) at about 0.3 mW, against
TTM-003-48H's 207 ms at 9.8 mW every 2 h
(`subjects/ttm-003/optogenetic_stimulation.md`). Before sorting it, this
asks whether the pulses evoke a neural response, in multi-unit activity,
sorted units and the LFP; whether the response calls for a change to the
sort's pulse rule (`docs/sorting_pipeline.md`, step 5.4); and whether
the pulses evoked OFF periods or slow waves more or less readily from
the start to the end of the sleep deprivation (09:00 to 14:00 on
2026-11-04, session hours 23.97 to 28.97, 16 trains).

In short: the pulses evoke a clear neural response, and no light artifact
larger than about 1 µV. Six tetrodes fire under half their rate through
the pulse and recover 50 to 100 ms after it, an OFF period follows 36 %
of pulses (9 % of sham onsets), and the LFP answers with a positive wave
during the pulse and a negative one after it. One TT4 unit fires
through the pulse in some trains, and one TT12 unit in two windows, with
full-size spikes. The OFF periods and LFP waves show no detectable change
from the start to the end of the sleep deprivation, within wide
intervals. No new pulse rule is needed. The existing rule adds fits to
the driven TT4 unit that are, by my reading, mostly false, so whether to
sort with the pulse table is an open choice (last section).

## The light edge

Averaged over all 2,900 pulses, the spike band (300 to 6000 Hz, no
reference) from 0.2 ms before to 0.7 ms after onset, where a
photoelectric artifact would sit, reaches at most 2.83 µV on any channel,
against 2.12 µV in the same average at sham times 500 ms before each
onset, and 1.4 µV at offset. The onset average exceeds the sham one on 61
of 63 channels, by a median of 0.5 µV and at most 1.1 µV; that may be a
small artifact, far below the spikes (troughs of 40 µV and more). Per
train, the largest absolute value of that window over the noise of the
train's average, on the largest channel, is 3.7 (median over trains),
against 3.6 in TTM-003-48H, whose pulses carry about 30 times the power.

## Multi-unit activity

Multi-unit crossings (5 MAD on any wire of a tetrode, merged within
1 ms), over every pulse, as a share of each tetrode's rate in the 500 ms
before onset:

| Tetrode | Rate in the 500 ms before onset (/s) | 1 to 10 ms after onset | 10 ms to the offset | 0 to 50 ms after the offset | 50 to 100 ms after the offset |
|---|---|---|---|---|---|
| TT1 | 7.9 | 1.05 | 0.99 | 1.05 | 1.11 |
| TT2 | 23.4 | 1.02 | 0.93 | 1.01 | 1.01 |
| TT4 | 61.3 | 1.77 | 0.79 | 0.37 | 0.91 |
| TT5 | 68.1 | 0.81 | 0.17 | 0.33 | 1.04 |
| TT6 | 17.3 | 0.97 | 0.25 | 0.63 | 1.33 |
| TT7 | 33.8 | 1.08 | 0.41 | 0.62 | 1.23 |
| TT8 | 46.7 | 0.91 | 0.47 | 0.90 | 1.32 |
| TT9 | 36.7 | 0.95 | 0.48 | 0.78 | 1.25 |
| TT10 | 53.9 | 0.90 | 0.72 | 0.85 | 1.08 |
| TT11 | 73.9 | 0.89 | 0.28 | 0.54 | 1.03 |
| TT12 | 59.4 | 1.00 | 0.82 | 0.82 | 1.02 |
| TT15 | 77.6 | 0.99 | 0.99 | 1.00 | 0.99 |
| TT16 | 11.1 | 0.88 | 0.93 | 0.99 | 1.11 |

TT3, TT13 and TT14 are left out: their median rate in the 60 s before
each train is under 5/s.

Pooled over all pulses, no tetrode's crossings rise from 10 ms after
onset to the offset; a light-driven unit on TT4 is outweighed there by
the tetrode's suppressed units (see Sorted units). Six tetrodes (TT5 to
TT9 and TT11) fall under half their rate from 10 ms after onset to the
offset. On TT5, TT6, TT7 and TT11 the low rate continues through the
50 ms after the offset (0.33 to 0.63), while TT8 and TT9 have mostly
recovered by then (0.90 and 0.78). All six recover 50 to 100 ms after the
offset, with a rebound of up to 1.33 times the rate before onset there
and at most 1.35 times in any 5 ms bin after that. TT4 alone answers the
onset: 1.77 times its rate from 1 to 10 ms after onset (0.42 extra
crossings per pulse), 0.79 from 10 ms to the offset and 0.37 in the 50 ms
after it.

## OFF periods and slow waves

An OFF period here is a gap of 50 ms or more between the pooled
crossings of the six suppressed tetrodes that spans onset plus 20 ms; it
ends where that gap ends, which can be after the pulse. The same is
measured at sham onsets, with no light: up to 40 onsets at the pulse
period, ending 2 s before each train's first pulse.

- An OFF period follows 36 % of pulses (1,031 of 2,900) and 9 % of sham
  onsets. It ends at a median of 104 ms after onset (quartiles 80 and
  126 ms), 40 ms after the pulse.
- The LFP (0.5 to 30 Hz, mean of the six tetrodes, relative to the
  200 ms before onset) rises during the pulse to 217 µV at 50 ms,
  crosses zero at 118 ms, reaches −144 µV at 172 ms and returns to zero
  at 338 ms. From 10 to 100 ms after onset (the early wave) it averages
  161 µV (sham −1 µV); from 100 to 300 ms (the late wave), −84 µV (sham
  0 µV). Both are larger on pulses with an OFF period (205 against
  137 µV, and −107 against −71 µV).
- My reading, a conjecture: the positive wave is the depth-positive
  side of a down state, as the silence it coincides with suggests, and
  the negative wave the return of firing after it.

The response varies with time and state across the 145 trains. Light and
dark follow the manifest (lights on at 08:53, off at 20:53; day 2 starts
at the second lights-on), and the state measure is the delta share
(0.5 to 4 Hz) of the six tetrodes' LFP power from 0.5 to 30 Hz in the
60 s before the train:

| Trains | n | OFF after pulses | OFF after sham onsets | Early wave (µV) | Late wave (µV) | Delta share before |
|---|---|---|---|---|---|---|
| Light, day 1 | 36 | 0.28 | 0.10 | 156 | −80 | 0.43 |
| Dark, day 1 | 36 | 0.39 | 0.08 | 201 | −105 | 0.40 |
| Sleep deprivation (light, day 2) | 16 | 0.33 | 0.005 | 133 | −51 | 0.32 |
| Light, day 2, after the sleep deprivation | 20 | 0.29 | 0.17 | 109 | −65 | 0.52 |
| Dark, day 2 | 36 | 0.43 | 0.09 | 168 | −91 | 0.46 |
| Light, day 3 | 1 | 0.55 | 0.00 | 161 | −74 | 0.28 |

Fitted together (least squares over the trains on session hour, delta
share, dark and sleep deprivation; 95 % intervals), the early wave falls
by 1.3 µV per hour (0.7 to 2.0) over the 48 h, is 60 µV larger in the
dark (40 to 80) and 121 µV smaller per unit of delta share (43 to 199),
and the sleep deprivation adds nothing beyond these (−9 µV, −40 to 22).
The share of pulses with an OFF period falls by 0.60 per unit of delta
share (0.38 to 0.82), rises by 0.003 per hour (0.001 to 0.005, about 0.13
over the 48 h), is 0.06 higher in the dark (0.01 to 0.12), and the sleep
deprivation adds −0.08 (−0.16 to 0.01). Sham OFF periods are rare in the
sleep deprivation (0.5 %, against 9 % overall); conjecture: the
spontaneous gaps are those of sleep, which the deprivation prevents.

## Across the sleep deprivation

The 16 trains of the sleep deprivation (73 to 88) show no detectable
change in the pulses' effect from start to end. The first five trains
against the last five, the rank correlation with time over all 16, and
the linear change over the 5 h with its 95 % bootstrap interval over
trains:

| Measure | First 5 | Last 5 | Spearman ρ (p) | Change over 5 h (95 % interval) |
|---|---|---|---|---|
| OFF after pulses | 0.38 (38 of 100) | 0.33 (33 of 100) | 0.01 (0.98) | −0.03 (−0.23 to 0.15) |
| OFF after sham onsets | 0.005 | 0.005 | −0.02 (0.95) | 0.00 (−0.02 to 0.02) |
| OFF end (ms from onset) | 95 | 89 | −0.09 (0.75) | −15 (−57 to 25) |
| Early wave (µV) | 138 | 131 | −0.03 (0.91) | −7 (−90 to 60) |
| Late wave (µV) | −51 | −57 | −0.03 (0.91) | −2 (−32 to 23) |
| Delta share before the train | 0.31 | 0.29 | −0.15 (0.57) | −0.03 (−0.13 to 0.09) |
| Six tetrodes' crossings before the train (/s) | 262 | 239 | −0.56 (0.02) | −26 (−47 to −6) |

The binomial 95 % intervals of the OFF shares are 0.28 to 0.48 and 0.24
to 0.43. With 16 trains of 20 pulses, the data cannot rule out a fall of
up to 0.23 or a rise of up to 0.15 in the share of pulses followed by an
OFF period, against 0.38 at the start. The six tetrodes' firing before
each train fell over the 5 h (by 9 %), with no detectable change in the
delta share.

## Sorted units

MountainSort5 with the seed sort's settings (`docs/sorting_pipeline.md`,
step 4: scheme 2, threshold 5.5), trained on 30 min sampled evenly over
each stretch, sorted two 5 h stretches of the preprocessed binary (step
2): the sleep deprivation (session hours 23.97 to 28.97, 320 pulses) and
the first 5 h of the recording (300 pulses). Units with template snr 5
or more and 100 spikes or more are kept. The counts below are of units
firing 0.5/s or more in the 500 ms before onsets, since a ratio or a
silence cannot be measured on the rest. The driven test is the pulse
rule's, applied to each unit's spikes per 1800 s window: from 15 ms
after each onset to the offset, 20 spikes or more, at least twice the
count expected from the unit's rate outside the pulses (onset to 100 ms
after the offset left out), and a one-sided Poisson p under 1e-4. The
rule itself applies it to looser candidates (matched-filter snr 2.5 or
more), which this does not reproduce. Unit numbers belong to each sort.

In the sleep deprivation, 53 units were kept, 39 of them at 0.5/s or
more.

- No unit passes the driven test in any of its 10 windows. The closest
  is TT4 unit 7 in the first window (24.0 to 24.5 h, trains 73 and 74):
  44 spikes from 15 ms after onset to the offset against 23.4 expected,
  1.88 times, with p 9e-5, so it fails only the bar of twice the
  expected count. In its other nine windows it has 0.2 to 1.5 times the
  expected count. Of the other units' 520 windows, two reach 20 spikes,
  at 1.60 and 1.45 times the expected count; higher ratios occur only on
  counts of 16 or fewer.
- Unit 7 answers the onset: 5.2 times its rate from 1 to 10 ms after
  onset and 1.33 times from 10 ms to the offset, as that response
  decays, then 0.37 and 0.62 times in the two 50 ms after the offset. No
  other unit exceeds 1.9 times from 1 to 10 ms. Whether it is the same
  cell as unit 6 of the first 5 h (below) was not checked.
- 19 of the 39 are suppressed (under half their rate from 10 ms after
  onset to the offset; median 0.21): 3 of 5 on TT4, 3 of 3 on TT5, 1 of
  1 on TT6, 3 of 4 on TT7, 2 of 2 on TT8, 2 of 3 on TT9 and 5 of 7 on
  TT11, none of the 14 on TT1, TT10 and TT12. Their firing regains half
  its rate by 67.5 to 102.5 ms after onset (median 72.5 ms; the first
  5 ms bin after the offset is centred on 67.5 ms), inside the rule's
  exclusion, which runs to 164 ms.
- For the 34 units with 10 or more spikes inside pulses, the median
  trough on the main wire inside pulses is 0.88 to 1.14 times the one
  outside (median 1.00). Unit 7's is 0.92 over 186 spikes, 74 of them
  from 1 to 10 ms after onset; by position in the pulse, 0.99, 0.87 and
  0.88 times (0 to 15, 15 to 30 and 30 ms to the offset).
- No unit's response changes detectably across the sleep deprivation:
  the rank correlation of a unit's per-train ratio (10 ms to the offset
  over the 500 ms before onset) with time has p under 0.05 for 1 of the
  39, about 2 expected by chance.

In the first 5 h, 60 units were kept, 41 of them at 0.5/s or more.

- TT4 unit 6 is light-driven: 10.8 times its rate from 1 to 10 ms after
  onset and 4.9 times from 10 ms to the offset. It passes the driven
  test in 6 of its 10 windows (0 to 2.5 h and 4.0 to 4.5 h, with 35 to
  57 spikes against 1.9 to 17 expected). It fails from 2.5 to 4.0 h,
  with 0 to 5 spikes, and from 4.5 to 5.0 h, where 18 spikes against 2.5
  expected fall under the minimum of 20. Its rate through the pulse
  varies from train to train between 0 and 52/s.
- Unit 6's spikes inside pulses are full size. On its main wire their
  median trough is 0.97, 0.96 and 0.90 times its median outside pulses
  from 0 to 15, 15 to 30 and 30 ms to the offset, and none of the 485 is
  under half. It fires up to 5 spikes in a pulse (2.6 on average in the
  pulses where it fires), at a median interval of 15 ms.
- TT12 unit 50 passes the driven test in 2 windows (1.5 to 2.5 h; 2.1
  times its rate through the pulse). Its spikes inside pulses are full
  size too (0.97 to 1.00 times by position in the pulse).
- 17 of the 41 are suppressed (median 0.21): 2 of 5 on TT4, 5 of 5 on
  TT5, 2 of 2 on TT6, 2 of 3 on TT7, 1 of 2 on TT8, 1 of 3 on TT9 and 4 of
  7 on TT11, none of the 14 on TT1, TT10 and TT12. Their firing regains
  half its rate by 67.5 to 107.5 ms after onset (median 77.5 ms).
- For the 35 units with 10 or more spikes inside pulses, the median
  trough inside pulses is 0.79 to 1.18 times the one outside (median
  1.00). The two lowest are TT7 unit 20 (0.79, 49 spikes; 2.8 times its
  rate from 1 to 10 ms after onset) and TT4 unit 8 (0.81, 26 spikes);
  none of their spikes inside pulses is under half.

A cell firing fast through a pulse can fire spikes too small for the
sort. To look for them, each unit's template was slid over the pulses,
from 15 ms after onset to the offset, and over the same window at sham
times 500 ms before each onset; the table gives the template-shaped
events (cosine 0.5 or more) of 0.2 to 0.4 times the template's size, per
second:

| Unit | Inside pulses | At sham times | Difference |
|---|---|---|---|
| TT4 unit 6, first 5 h (driven) | 38.2 | 22.7 | 15.5 |
| TT4 unit 7, sleep deprivation (not driven) | 43.3 | 23.3 | 20.0 |
| TT4 unit 9, first 5 h (suppressed) | 1.8 | 0.5 | 1.3 |
| TT5 unit 12, first 5 h (suppressed) | 16.8 | 9.3 | 7.6 |
| TT12 unit 50, first 5 h (driven) | 20.0 | 22.1 | −2.2 |

The driven unit 6 shows no larger excess of small events than unit 7 on
the same tetrode in the sleep deprivation, which passes no driven test
there (though it answers the onset and fires 1.33 times its rate through
the pulse), and the suppressed TT5 unit 12 shows one too. My reading, a
conjecture: the excess comes from the other cells' silence, not from
small spikes of a driven cell; with fewer large events in the window,
small template-shaped peaks are less often masked by a larger peak
within 1 ms.

Over the whole recording, TT4's multi-unit crossings from 15 ms after
onset to the offset exceed 1.2 times its rate before onset in 14 of the
145 trains (median train 0.71), against 3 to 8 trains on the tetrodes of
similar rate that are not suppressed (TT10, TT12 and TT15). The
tetrode's suppressed units outweigh the driven one in its pooled
crossings, so this does not show how often the cell is driven over the
48 h. Conjecture: its drive depends on state, as the train-to-train
variation of unit 6 suggests.

## TTM-003-48H for comparison

The same analysis of TTM-003-48H's 24 trains (1,432 pulses of 207 ms at
9.8 mW):

- Light-driven firing: TT4 fires 4.6 times its rate before onset through
  the pulse (15.9 times from 1 to 10 ms) and TT11 1.7 times (7.0 from 1
  to 10 ms). Seven other tetrodes with a rate of 5/s or more answer the
  onset at 1.6 to 5.4 times their rate.
- The suppressed tetrodes are TT5 to TT9 and TT12 (0.03 to 0.28 through
  the pulse). In TTM-003-48H-2, TT12 is barely suppressed (0.82) and TT11,
  driven in TTM-003-48H, is suppressed (0.28).
- An OFF period follows 68 % of pulses (sham 5 %) and ends at a median of
  287 ms after onset (quartiles 169 and 333 ms), 80 ms after the pulse.
- The early wave averages 134 µV (sham 7 µV). The late wave's window
  lies mostly inside that recording's pulse, so the two recordings' late
  waves are not comparable.

At about a thirtieth of the power, the TTM-003-48H-2 pulses drive a TT4
unit through the pulse in some windows (and a TT12 unit in two), at tens
of spikes per second rather than hundreds, with full-size spikes, and
evoke OFF periods
after half as many pulses that end sooner. Whether TT11's change from
driven to suppressed reflects the weaker light or a different position
of the tetrode is not known; the manifest's observation says TT6, TT7
and TT9 probably moved between the recordings.

## The pulse rule on the carry

The sort matches templates per tetrode in consecutive 1800 s windows,
carrying the templates from one window to the next (the carry,
`docs/sorting_pipeline.md`, step 5), starting from templates of a
MountainSort5 sort of the first window (the seed bank, step 4) in a
whitened space (step 3). Its normal rules accept a fit with
matched-filter snr 4.5 or more, an amplitude of at least half the unit's
median and a cosine of 0.5 or more to the template. Inside pulses
(onset to 100 ms after the offset), the pulse rule accepts a driven
unit's fits down to snr 2.5 and 0.3 of its median, with the same
cosine.

To see what the rule does here, the carry ran windows 0 to 4 (the first
2.5 h, where TT4 unit 6 is driven) twice from one seed bank and
whitening: once with the pulse table and once without. Its unit 4 holds
91 % of unit 6's spikes in those windows, and unit 6's spikes are 37 %
of its own (it fires at 11.7/s over the five windows); it is called
driven in windows 0, 1, 3 and 4, not in window 2, and no other unit is
called driven.

In window 0 the two runs make the same fits (the matching itself does
not use the pulses, and nothing earlier differs), so the difference is
the rule alone: it accepts 32 more of unit 4's fits and drops 1
(refractory), out of its 18,877. Its spikes that both runs accept from
15 ms after onset to the offset are full size (median amplitude 0.99 of
the template, snr 8.8, cosine 0.82). The 32 added fits are small and
poorly shaped: in each part of the pulse their median amplitude is 0.40
to 0.44 of the template (its accepted spikes outside pulses have a
median of 0.80), their snr 3.6 to 4.0 and their cosine 0.52. Of them, 5
fall 1 to 10 ms after onset, 8 from 15 ms to the offset, 18 in the 100
ms after the offset and 1 elsewhere in the first 15 ms. In the 100 ms
after the offset, where the unit fires at 6.5/s against 10.5/s over
window 0, they raise its count from 26 to 44. In windows 1, 3 and 4 the
fits differ slightly between the runs, and the counts differ the same
way: 4 to 8 more spikes from 15 ms to the offset and 7 or 8 more in the
100 ms after it. My reading, a conjecture: the added fits are mostly not
the driven cell's spikes, since its spikes inside pulses are full size
and the matched filter finds no more small events for it than for a unit
on the same tetrode that is not driven.

Without the pulse table, the relative floor and template updates also
take in the spikes from onset to 100 ms after offset, 3 to 7 s of each
1800 s window. Over the five windows the other units' accepted counts
differ between the runs by −4.6 to +3.0 % (by at most 1.9 % for units
with 10,000 spikes or more in a window). Two carries of TTM-003-48H that
differ only in the binary's sample type (int16 against float32) changed
a unit with over 10,000 spikes by up to 2.8 %
(`investigations/carry-bakeoff/`).

## Whether sorting needs a new pulse rule

No. The pulse rule (`docs/sorting_pipeline.md`, step 5.4) rescues the
smaller spikes of cells that fire at 100 to 200 Hz through a pulse. In
TTM-003-48H-2 the driven cells fire at tens of spikes per second, with
full-size spikes that the normal rules accept. Single units regain half
their rate by 108 ms after onset, and three quarters of the pooled OFF
periods end by 126 ms (median 104 ms), inside the 164 ms from onset that
the relative floor and template updates leave out. Nothing in the
binary, the bad periods or the seed depends on the pulses.

The existing rule still has an effect: in each window where the driven
unit passes the test, it adds fits to it that are, by my reading, mostly
false: 32 in window 0, where the two runs make the same fits, 18 of them
in the 100 ms after offset; in windows 1, 3 and 4, 4 to 8 more spikes
from 15 ms to the offset and 7 or 8 more in the 100 ms after it. That is
0.2 % of the unit's spikes, but it raises its count from 15 ms to the
offset by 8 to 14 % and in the 100 ms after the offset by 35 to 69 %, in
some windows and not in others. Two ways to sort TTM-003-48H-2:

- With the pulse table, as TTM-003-48H was sorted: the rule as it is,
  with those added fits.
- Without it: no added fits; pulse times enter the relative floor and
  template updates, which on windows 0 to 4 changed the other units by
  no more than the binary's sample type changed TTM-003-48H's.

Whether the 100 ms after offset adds false fits to TTM-003-48H's driven
units too, inside their longer OFF periods, was not measured.

Decided on 2026-11-06: neither way. The carry now gives the pulse table's
two uses separate spans and the pulse rule its own switch
(`docs/sorting_pipeline.md`, step 5.4), and TTM-003-48H-2 is to be sorted with the pulse table and
`--no-pulse-relax`: pulse times stay out of the relative floor and the
template updates, and the pulse rule does not run. In TTM-003-48H the
100 ms after the offset did add fits, 112 of them more than 5 ms after
the offset over the 48 h; the rule's span now ends 5 ms after the offset
(`investigations/pulse-rule-span/`).

## Scripts and data

- `subjects/ttm-003/scripts/opto_off_periods.py`: per train, from the
  session's recording store (the lossless Zarr copy of the raw data),
  multi-unit crossings per tetrode (300 to 6000 Hz, median of the
  good channels subtracted, 5 MAD), the LFP (0.5 to 30 Hz, mean per
  tetrode), the light-edge check, and the same at sham onsets; the
  summary gives the tetrode table, OFF periods, LFP waves, the light and
  dark groups, the fits and the sleep-deprivation trend. Results:
  `subjects/ttm-003/results/opto_off_periods_<E>.csv` and `.json`;
  figures in `<experiment>/figures/opto_off_periods/`.
- `subjects/ttm-003/scripts/opto_edges.py`: the spike-band average over
  every delivered pulse at onset, offset and sham times
  (`results/opto_edges_<E>.json`, `<experiment>/figures/opto_edges.png`).
- `subjects/ttm-003/scripts/opto_units.py`: MountainSort5 on a stretch of
  the preprocessed binary, and each unit's response, spike sizes inside
  pulses and the pulse rule's driven test per window
  (`results/opto_units_<E>_<name>*.csv`, `.json`; figures in
  `<experiment>/figures/opto_units_<name>/`). The stretches are `sd`
  (`--start-h 23.9727 --end-h 28.9727`) and `day1` (`--start-h 0
  --end-h 5`).
- `subjects/ttm-003/scripts/opto_small_spikes.py`: the matched filter of
  a unit's template inside pulses and at sham times
  (`results/opto_small_spikes_<E>.csv`).
- `subjects/ttm-003/scripts/opto_carry_rule.py`: the two carry runs with
  and without the pulse table compared, and the driven carry unit
  matched to the `day1` sort (`results/opto_carry_rule_<E>.csv`,
  `_added.csv`, `_match.csv`). The runs, in
  `<session>/sortings_carry/pulsetest/` and `nopulsetest/`:
  `ttm-sort seed --run pulsetest --binary <session>/sortings_carry/main/binary`,
  then `ttm-sort carry` with `--windows 0-4`, once with `--run pulsetest
  --pulses subjects/ttm-003/results/laser_pulses_2026-11-03_09-01-38.csv`
  and once with `--run nopulsetest` and that run's `--bank` and
  `--whitening`. Both runs were deleted on 2026-11-06. A rerun with the
  current code needs `--pulse-relax-post-ms 100` for the span they used;
  the current pulse rule also takes no fit inside a bad period, so the
  counts can differ slightly.
