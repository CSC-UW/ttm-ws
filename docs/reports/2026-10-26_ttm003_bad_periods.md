---
title: TTM-003 48 h, transient bad periods per tetrode and what masking them does to the sorting
updated: 2026-10-27
status: needs_verification
---

# TTM-003 48 h: transient bad periods, detected per tetrode and masked before sorting

Run on 2026-10-26 on `TTM-003-48H` (47.98 h). The user found TT4 at
2316 s full of noise on all four wires, with the sorter producing false
spikes there, and asked for a way to detect transiently bad channels and
keep them out of sorting without losing the rest of the recording. The
code is `ttm.bad_periods` and `ttm-bad-periods`
(`docs/developer_notes/bad_periods.md`); the sorting with the mask is
`subjects/ttm-003/scripts/mp_baseline.py --bad-periods`; the comparison is
`subjects/ttm-003/scripts/bad_periods_compare.py`. Measurements are stated
as such; sentences marked *conjecture* are the analyst's reading and have
not been confirmed by the user.

## 1. Summary

- The detector flags, per tetrode and second, a wire three times louder
  than its own 30 min running median or a signal shared by the four wires
  that adds power. On TTM-003 it flags 0.01 to 0.03 % of the 48 h on the
  eight tetrodes with most units (TT5 to TT12; under a minute each, all
  single electrical pops), 0.3 % on TT4, and 0.3 to 7 % on the tetrodes with
  few units. On TTM-001's 47 h it flags 56 s in all.
- TT4's bad stretch comes out as one interval, 2287 to 2696 s, covering
  both parts the user labelled: 2296 to 2359 s with all four wires
  correlated at 0.93 (normal 0.59), and from 2373 s with CH58 at about 10
  and CH64 at about 3.5 times their baseline while CH62 and CH60 stay near 1.2.
- In the June baseline the first two windows put 832,069 single-unit spikes
  into flagged time, 632,048 of them on TT4. Re-estimated on them, TT4's
  seed units 4 and 5 take on the artifact's waveform and never fire again;
  units 6 and 7 fall from 40,472 and 22,906 clean spikes in window 1 to 779
  and 171 in window 2 and need 3 to 4 h to come back to the level they have
  in the masked run.
- The same event costs the eight odd-channel tetrodes, which have no
  flagged time, most of their spikes for its 400 s: circus-omp fits all
  tetrodes of a 1 s chunk together, and the artifact fits on the
  even-channel tetrodes crowd theirs out. Matching only their own templates
  on the same seconds gives back the masked run's count.
- Masked, the same first 14 windows (7 h) put no spikes in flagged time,
  keep all four TT4 seed units present in every window with stable
  waveforms, give TT4 2.2 times as many clean-time spikes (1,332,262
  against 608,061), and claim 92.9 % of the clean-time events at 10 MAD or
  more against 82.3 % in the baseline and 80.6 % in a second unmasked run
  with another template seed.
- Over the full 48 h the baseline, and a second unmasked run with another
  template seed, put 1.45 M and 1.59 M single-unit spikes into flagged time
  and have 19 and 17 units with most of their spikes there; the masked run
  has none. At the second even-channel event (26.9 h) the unmasked runs
  lose ten and eight units on TT1 and TT13 to TT16 for good and claim almost
  none of TT1's, TT13's, TT15's and TT16's events for 3 h, until the next
  re-seed; masked, no unit stops. Pooled clean-time coverage at 10 MAD or more is 90.3 % masked,
  88.3 % baseline and 88.0 % control, and the masked run is at or above the
  control on every tetrode but TT5 (83.3 against 85.6 %), where the cause is
  not settled.
- The one baseline unit the masked run lacks, on TT15, fires like every
  TT15 unit of both runs at a fixed phase of a 300-sample cycle, the period
  of a signal on CH22 that is locked to the acquisition clock. In the 48 h
  baseline eight single units with 14.7 M spikes are locked this way; the
  mask does not address it.

## 2. Terms

- Bad period: an interval of the table, scoped to one tetrode. Flagged time
  of a tetrode is the union of its intervals; clean time is the rest.
- Loudness of a wire: its rms in a 1 s bin over its running median, the
  median taken over 30 min and re-estimated without the bins in which some
  wire of the tetrode was more than twice its baseline.
- Masked sorting: the June pipeline run on the preprocessed binary read
  through `mask_bad_periods`, in which each flagged tetrode's four channels
  are replaced inside its intervals by white noise at their clean noise
  level.
- Window: one of the 96 windows of 1800 s of the June carry-forward.
- Coverage: the share of detected events (negative peaks beyond 5.5 MAD,
  the baseline's `spike_coverage.npz`) that a single unit of the same
  tetrode claims within 0.5 ms; the scoreboard's axis A, restricted here to
  clean time and scored on the same events for both runs.

## 3. What the detector finds on TTM-003

Seconds flagged per tetrode over the 48 h (`<session>/bad_periods/bad_periods_summary.csv`):

| Tetrodes | flagged | longest interval |
|---|---|---|
| TT5 to TT12 | 21 to 47 s each (0.01 to 0.03 %) | 3 to 5 s |
| TT4 | 472 s (0.27 %) | 409 s |
| TT2, TT13 | 471 and 488 s | 404 and 333 s |
| TT1, TT15 | 1161 and 1280 s | 555 and 645 s |
| TT3, TT14 | 2166 and 2124 s | 408 and 653 s |
| TT16 | 12,399 s (7.2 %) | 555 s |

On TT5 to TT12 every flagged interval looked at was a single pop of a
millisecond or two, on one wire or on all four, many times larger than any
spike, sometimes with ringing; a few hit several tetrodes in the same
second. The 1 s bin and 1 s padding make each pop cost 3 s.

The tetrodes on even-numbered channels (TT1 to TT4, TT13 to TT16; every
wire of these eight tetrodes has an even channel number, every wire of TT5
to TT12 an odd one) fail together. From 2287 to about 2700 s seven of them
become 3 to 40 times louder at once, and from 96,900 to 97,560 s (26.9 h,
during the sleep deprivation) six of them; the odd-channel tetrodes show
nothing in either. *Conjecture: a fault on the half of the connector or
headstage that carries the even channels.* The same tetrodes, except TT4,
are the ones with few units in every sorting so far.

On those few-unit tetrodes most intervals are short: bursts of a few
seconds of high-frequency signal seen nearly identically on all four
wires, often in a rhythm of a few per second. Whether these are muscle
pickup or neural activity is not established; they are flagged because they
are shared and loud, not because they were shown to be artifact.

14 of the 24 delivered laser trains overlap a flagged interval, almost all
for a few seconds on TT3, TT14, TT15 or TT16; on the tetrodes with tagged
candidates only train 28 is touched (TT4 and TT5, 3 s at 50,976 s, one pop).

## 4. TTM-001 as a check

The same detector on TTM-001's 47.2 h preprocessed binary
(`mp_long_s2000_d170000`, which starts at session time 2000 s) flags 17
intervals, 56 s in all: 16 on TT6 and TT7, each a burst of about 50 ms
oscillating at a few hundred Hz on all four wires, often on both tetrodes at
once, and one pop on CH58 of TT4. Nothing on the other thirteen tetrodes.
The scan and table are in
`/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52/bad_periods_crop_s2000_d170000/`
(times relative to the crop).

## 5. How the mask enters the sorting

The June pipeline learns from the data at every stage: the seed sort and
the re-seeds (MountainSort5 per tetrode), the per-window template
re-estimation of the carry-forward (each template becomes the mean of the
spikes it took in the window, when it took 100 or more), the noise levels,
and the residual and multi-unit passes. `mask_bad_periods` sits between the
binary and all of them, so none of them sees a flagged tetrode inside its
intervals. The fill is noise rather than zeros so that MAD estimates over a
window that contains an interval are not pulled down (`docs/developer_notes/bad_periods.md`).

## 6. Head to head over the first 14 windows (7 h)

`mp_baseline_masked_w14` is the June carry-forward (seed sort, circus-omp
matching with per-window re-estimation, the re-seed at window 12) on the
masked binary, stopped after 14 windows; the baseline is the first 14
windows of `mp_baseline/assembled_reseed_c12`. Both are scored on the
unmasked binary (`results/bad_periods_*_assembled_reseed_c12_w14.*`).

Single-unit spikes in flagged and in clean time:

| Tetrode | flagged s | baseline, flagged | baseline, clean | masked, flagged | masked, clean |
|---|---|---|---|---|---|
| TT1 | 425 | 97,037 | 1,209 | 0 | 2,788 |
| TT2 | 434 | 6,451 | 3,535 | 0 | 4,407 |
| TT3 | 706 | 50,200 | 1,041 | 0 | 3,844 |
| TT4 | 412 | 632,048 | 608,061 | 0 | 1,332,262 |
| TT14 | 448 | 26,706 | 7,161 | 0 | 11,922 |
| TT15 | 502 | 8,488 | 1,769,182 | 0 | 939,911 |
| TT16 | 1221 | 10,214 | 8,677 | 0 | 21,635 |
| TT5 to TT13 | 0 to 6 each | 925 in all | 7,378,037 | 0 | 7,649,786 |

What happens on TT4 (`figures/bad_periods/counts_assembled_reseed_c12_w14_TT4.png`
and `templates_..._TT4.png` under the experiment's figures):

- Baseline. In window 1 seed units 4 and 5 take 513,657 and 31,321 spikes in
  flagged time, against 42,832 and 29,905 clean spikes over the 14 windows,
  nearly all of them in windows 0 and 1.
  The mean waveform of their window-1 spikes is the artifact, a ringing
  oscillation on all four wires, and that is the template they carry into
  window 2, where they match nothing: 0 and 1 spikes, and no more than 6
  per window to window 13. Their templates are never re-estimated again,
  because re-estimation needs 100 spikes. Units 6 and 7 take 45,009 and
  42,061 spikes in flagged time, their window-1 waveforms are distorted
  toward the artifact, and they take 779 and 171 spikes in window 2,
  against 40,472 and 22,906 clean spikes in window 1. Unit 6 is back to the
  masked run's count by window 6 (14,529 against 18,305), unit 7 by window 8
  or 9.
- Masked. Units 4, 5, 6 and 7 have 10,000 to 50,000 clean spikes in every
  window, and their mean waveforms in windows 0 to 3 are the same shape.
  Over the 14 windows unit 4 has 335,833 spikes (baseline 42,832; 98 % of
  the baseline's spikes are among them), unit 5 349,856.
- Both runs add one TT4 unit at the window-12 re-seed (baseline 49, masked
  47, 86 % of their spikes in common); it is a unit the seed sort did not
  find, not a reborn unit 4 or 5.

The tetrodes with few units show the same pattern on a smaller scale: their
seed units (0 on TT1, 1 on TT2, 3 on TT3) take 6,000 to 97,000 spikes in
flagged time in window 1, then go silent in the baseline; masked, they stay
present in 9 to 12 of the 14 windows. The median cosine of their spikes to
their own window template is 0.51 to 0.54 in the baseline and 0.74 to 0.77
masked.

The tetrodes with no flagged time lose spikes too. From 2287 to 2700 s,
while seven even-channel tetrodes carry the event, TT5 to TT12 have no
flagged time and the same data in both runs. In the baseline their units
take 2,000 to 9,000 spikes a minute from 2340 s on, against 13,600 to
29,500 a minute earlier in the window, and the share of their events at
10 MAD or more that a unit claims falls from about 95 % to 42 to 76 %.
Masked, it stays at 95.5 to 96.9 % and the spikes at 28,900 to 33,600 a
minute; 171,166 more TT5 to TT12 spikes than the baseline between 2280 and
2700 s (`investigations/bad-periods-ttm003/window1_coverage_output.txt`).
It is this event only: over the 14 windows, TT5 to TT12 events in seconds
when one to four other tetrodes are flagged are claimed at 94 to 95 % in
both runs, as in seconds with none flagged (93.6 and 93.8 %), and in the
408 s with five to eight flagged at 59.6 % in the baseline and 96.0 %
masked (`good_tetrode_coverage_by_flags_output.txt`).

The cause is the matcher. circus-omp fits the templates of all 64 channels
to each 1 s chunk together. With window-1 templates rebuilt from the
baseline's window-0 spikes, 2400 to 2460 s gives 4,943 TT5 to TT12 spikes
unmasked (and 133,863 artifact fits on the other tetrodes), 33,372 masked,
and 33,432 unmasked when only the TT5 to TT12 templates are matched
(`omp_coupling_output.txt`). Unmasked, the TT5 to TT12 spikes that remain
fall from 938 in the first tenth of a chunk to 316 in the last; masked they
are even (3,151 to 3,760). *Conjecture: the chunk's pursuit stops early. Each
pass walks its selected peaks in time order and abandons the rest of the
pass when a selected atom is dependent on those already chosen, and the
chunk ends after five passes that add no valid spike (`max_failures=5`).
With hundreds of artifact fits in a chunk, both happen before the clean
tetrodes' spikes are reached.*

Pooled over all tetrodes, single units claim 82.3 % of the clean-time
events at 10 MAD or more in the baseline and 92.9 % masked (all clean
events: 52.5 and 56.8 %;
`results/bad_periods_coverage_TTM-003-48H_assembled_reseed_c12_w14.csv`).
Per 60 s (`figures/bad_periods/coverage_assembled_reseed_c12_w14.png`), the
baseline falls to 42 % during the event, stays between 59 and 92 % (median
75 %) from 1 to 5 h while TT4's units come back, and is within 3 points of
the masked run from 5.3 h on; the masked run stays between 85 and 98 %.
Among the tetrodes with units, TT4 goes from 63.5 to 94.5 %, TT13 from 73.9
to 92.8 %, and TT5 to TT12 up by 0.8 to 1.9 points each (TT9, with 27 such
events, by 7.4).

To see how far two unmasked runs differ, the first 14 windows were sorted
again with only the random spikes of the templates changed
(`mp_baseline.py --template-seed 1`;
`results/bad_periods_*_assembled_reseed_c12_w14_seed1.*`). The control
repeats the baseline's damage: 853,837 single-unit spikes in flagged time,
558,349 clean spikes on TT4, 80.6 % of the clean events at 10 MAD or more
claimed, and the same collapse on TT5 to TT12 from 2287 to 2700 s. On TT5
to TT13 its units agree with the baseline's as closely as the masked run's
do: median share of spikes in common 0.833 against 0.830, median absolute
log2 ratio of clean spike counts 0.066 against 0.068. The masked run's
counts there are shifted up (median ratio 1.034 against 0.996; 9 of 32
units differ by more than 10 %, against 4 of 33), and 171,166 of its
271,749 extra clean spikes on TT5 to TT13 are the ones recovered between
2280 and 2700 s. The median cosine of their spikes to the window template
does not change: 0.718 masked against 0.721 for the same baseline units,
0.719 in the control against 0.717.

TT15 looked like the cost. Its three baseline seed units have 1.77 M clean
spikes; masked, the seed sort of TT15 (6 s masked, at 157 and 396 s, CH24
loud in time with bursts on TT1, TT3 and TT16) found two clusters instead
of three. Baseline units 42 and 43 are masked units 41 and 42 (99 % and
98 % of spikes in common), and baseline unit 41 (860,592 clean spikes in
7 h) has no counterpart. Section 7 shows that all of these units follow a
signal locked to the acquisition clock. The three-cluster split is also
fragile without the mask: re-sorted with whitening seeds 1 to 4, or with
the seed window moved by 5, 10, 30 or 60 s, the baseline's smallest cluster
(442 spikes in the window, the seed of unit 41) comes out on its own in 2
of these 8 sorts, and in none of 5 masked ones
(`tt15_seed_stability_output.txt`). Without unit 41, TT15's clean-time
coverage goes from 59.8 to 58.6 % (45.7 to 42.3 % at 10 MAD or more); 95 %
of its spikes fall on no detected event. The seed sort of every tetrode
without masked time in window 0 (TT2, TT4 to TT14) is identical in both
runs.

## 7. TT15's units follow the acquisition clock

Every TT15 unit of both runs fires at a fixed phase of a 300-sample (10 ms)
cycle counted from the first sample. Baseline unit 43 has 99.2 % of its
spikes in one tenth of the cycle; units 41 and 42 have none in four and
five of the ten tenths; 66, 14 and 98 % of their intervals are within 0.1 ms
of a multiple of 10 ms. In these 14 windows, units with 100 spikes or
more on the other tetrodes have at most 20 % of their spikes in one tenth
(median 10.5 %, the uniform value being 10 %).

The traces show the source. The mean over 60 s of the 300-sample cycle is
10 times the MAD peak to peak on CH22, 6 to 7 on CH20 and 4.5 on CH18 (2 to
3 on CH40 of TT1 and on TT14), at 0.3, 8.3, 27.8 and 44.4 h, against a
5.5 MAD detection threshold. The CH22 cycle keeps its trough at sample 12
or 13, with zero lag against its shape at 0.3 h, from 0.3 to 47.2 h
(`investigations/bad-periods-ttm003/tt15_clock_lock_output.txt`). A
source that holds its phase to the sample over 48 h is synchronous with the
acquisition clock. That differs from the ~400 Hz comb studied on TTM-001's
rig, which kept its frequency when the sample rate was changed
(`investigations/acquisition-noise/`); CH22 is the channel `channel_qc.py`
found carrying a 400 Hz comb on this recording. *Conjecture: a 100 Hz signal
scheduled from the acquisition board's clock, such as the headstage IMU this
recording samples at 100 Hz.*

In the 48 h baseline `assembled_prod`, eight of TT15's single units (units
41, 42, 43, 58, 63, 73, 74 and 82, 14.7 M spikes) and its multi-unit unit
are locked this way; three units on TT1 and TT14 reach 21 to 30 %.
Subtracting the cycle mean does not bring CH22's 160 crossings a second
below -5.5 MAD down, so the channel carries more than the locked part. CH22 is the manifest's bad channel, but the June
pipeline sorts all 64 channels. The bad-period mask does not touch this,
since the signal never stops; `bad_periods_compare.py` reports it per unit
as `clock_lock`.

## 8. The full 48 h

`mp_baseline_masked` is the same carry-forward over all 96 windows; the
baseline is all of `mp_baseline/assembled_reseed_c12`
(`results/bad_periods_*_assembled_reseed_c12.*`). Its first 14 windows are
the 14-window run of section 6, spike for spike (51 of 51 units). The carry
took 7.7 h against the baseline's 6.8 h, both on a shared machine. The
baseline's seed bank of 45 units grew to 88 by re-seeding; the masked
run's 44 grew to 73.

Single-unit spikes in flagged and in clean time:

| Tetrode | flagged s | baseline, flagged | baseline, clean | masked, flagged | masked, clean |
|---|---|---|---|---|---|
| TT1 | 1161 | 220,836 | 31,177 | 0 | 33,747 |
| TT2 | 471 | 6,540 | 27,251 | 0 | 38,243 |
| TT3 | 2166 | 51,226 | 2,001 | 0 | 11,391 |
| TT4 | 472 | 638,685 | 9,220,349 | 0 | 9,782,111 |
| TT13 | 488 | 6,369 | 184,760 | 0 | 246,950 |
| TT14 | 2124 | 37,681 | 10,864 | 0 | 33,531 |
| TT15 | 1280 | 329,182 | 14,526,989 | 0 | 11,744,787 |
| TT16 | 12,399 | 138,065 | 85,689 | 0 | 134,280 |
| TT5 to TT12 | 21 to 47 each | 17,434 in all | 59,932,842 | 0 | 60,777,903 |

In the baseline 19 single units have more spikes in flagged time than in
clean time: seven seed units (0, 1, 3, 4 and 5 of section 6, 40 on TT14 and
44 on TT16) and twelve units born at re-seeds on TT1, TT14, TT15 and TT16.
The masked run has none.

The second even-channel event (96,900 to 97,560 s, windows 53 and 54) does
in the baseline what the first did to TT4, to more tetrodes. Ten units stop
for good there: five on TT15 (41, 42, 43 and 58, together 9.0 M clean
spikes, and 68), 44 and 52 on TT16, 39 on TT13, 53 on TT1 and 67 on TT14.
Until window 61, when units from the re-seed at window 60 begin on TT1,
TT13, TT15 and TT16, the baseline claims 0.0 to 0.4 % of the clean events
at 10 MAD or more on these four tetrodes in each of windows 55 to 60 (3 h);
masked, 94 to 100 % on TT13, 60 to 93 % on TT15, 28 to 38 % on TT16 and 6
to 17 % on TT1
(`investigations/bad-periods-ttm003/window_coverage_output.txt`). No masked
unit present before either event stops at it. The odd-channel tetrodes lose their
spikes during the event as in window 1, more completely: from 97,260 to
97,500 s TT5 to TT12 have 1 to 4 spikes a minute in the baseline against
23,600 to 30,400 masked, 261,201 fewer over the 660 s
(`event_coverage_output.txt`).

Pooled over all tetrodes, single units claim 88.3 % of the clean-time
events at 10 MAD or more in the baseline and 90.3 % masked (all clean
events: 53.6 and 55.8 %). In seconds when five to eight other tetrodes are
flagged, 40.4 against 92.4 %; three or four, 59.4 against 76.9 %; none,
89.0 against 90.7 %. Per tetrode:

| Tetrode | events | baseline % | masked % |
|---|---|---|---|
| TT1 | 74,332 | 22.7 | 26.9 |
| TT2 | 23,976 | 45.9 | 52.3 |
| TT3 | 5,565 | 11.0 | 44.2 |
| TT4 | 2,994,387 | 85.1 | 90.6 |
| TT5 | 1,431,183 | 86.2 | 83.3 |
| TT6 | 1,484,036 | 93.6 | 94.1 |
| TT7 | 108,502 | 95.9 | 96.0 |
| TT8 | 718,508 | 98.4 | 99.4 |
| TT9 | 530 | 77.9 | 81.9 |
| TT10 | 1,163,435 | 98.0 | 98.4 |
| TT11 | 492,943 | 82.4 | 84.8 |
| TT12 | 822,464 | 98.3 | 98.9 |
| TT13 | 29,248 | 86.7 | 95.8 |
| TT14 | 27,101 | 18.9 | 36.7 |
| TT15 | 93,564 | 67.7 | 73.4 |
| TT16 | 165,312 | 25.9 | 39.7 |

A second unmasked run over the 48 h with only the template seed changed
(`mp_baseline_seed1`; `results/bad_periods_*_assembled_reseed_c12_seed1.*`)
repeats the baseline: 1,593,598 single-unit spikes in flagged time, 17 units
with most of their spikes there, 88.0 % of the clean events at 10 MAD or
more claimed, and the same collapse at the second event (0.0 to 3.7 % on
TT1, TT13, TT15 and TT16 in each of windows 54 to 60, after eight units
stop there). Its coverage differs from the baseline's by at most 2.2
points on any tetrode (TT11). On TT5 to
TT13 the masked run's units agree with the baseline's as closely as the
control's do: median share of spikes in common 0.674 against 0.663, median
absolute log2 ratio of clean spike counts 0.087 against 0.105, 19 against
18 of 45 units off by more than 10 %.

The masked run claims at least as much as the control on every tetrode but
TT5: 83.3 % against 85.6 % (baseline 86.2 %). The masked run is below the
control in every window from 56 on, by 0.1 to 5.2 points, and 4.2 to 6.5
points below the baseline from window 89. TT5's clean events at 10 MAD or
more go from 2,500 to 12,000 a window between windows 56 and 69 to 23,000
to 49,000 from window 72, and its coverage then follows the unit each run
adds at the window-72 re-seed: in the even windows from 88 to 94 baseline
unit 79 takes 28,000 to 46,000 spikes a window, the control's 11,000 to
16,000 and the masked run's 14,000 to 19,000; the control too is 1.3 to 4.7
points below the baseline from window 89. TT5 has 39 s of flagged time,
so the mask reaches its sorting only through the other tetrodes; whether its
lower TT5 coverage is that or the spread between runs, one control does not
decide (`investigations/bad-periods-ttm003/window_coverage_output.txt`).

Of the single units, 11 in the baseline and 7 masked are clock-locked
(`clock_lock` above 0.2): in the baseline nine on TT15 and TT1's units 46
and 53, which have most of their spikes in flagged time; masked, all seven
on TT15. The masked run's TT15 has fewer clean spikes (11.7 M against
14.5 M) and claims more of its events (73.4 against 67.7 %).

## 9. What this does not do, and what is next

- The mask drops all four wires of a tetrode. On TT4 from 2373 s only CH58
  and CH64 are loud; the table lists the loud wires, so a per-wire mask is
  possible, but templates and cosines use all four wires and would have to
  be computed on the remaining ones.
- The staged pipeline (`mp_staged.py`) reads the binary by path in its
  workers and does not take the table yet; the pulses-only run would need it
  to sit on the masked baseline.
- The analyzer of a masked run reads the masked binary, so its trace view
  shows noise inside flagged intervals. The table can be drawn over the
  traces once the GUI has interval shading.
- The matcher couples tetrodes: an artifact on some tetrodes costs the
  others their spikes in the same 1 s chunks unless it is masked. A bad
  period the detector misses therefore costs more than its own tetrode.
  Matching each tetrode on its own channels (`ttm.matching`,
  `mp_baseline.py --per-tetrode-matching`) removes the coupling: during the
  2287 s event it keeps TT5 to TT12's spikes unmasked. On clean data it
  also adds about a quarter more spikes, mostly poorly shaped or shared with
  another tetrode, which the joint matcher's competition across tetrodes
  kept out (`investigations/per-tetrode-matching/`). It is off by default.
- TT15's clock-locked units remain in both runs. Dropping CH22, the
  manifest's bad channel, from TT15 before sorting, or subtracting the
  cycle mean where it is the main contaminant, has not been tried.
- Thresholds were set by looking at TTM-003 (TT4's labelled stretch, the
  pops on the good tetrodes, the bursts on TT16) and checked on TTM-001;
  they have not been validated against labels beyond TT4's.
- The detector misses the tail of the second event on TT1 and TT16: from
  97,462 to about 97,540 s, after their intervals end, their wires stay at
  1.3 to 2.5 times baseline with 68 to 270 crossings of -5 times the noise
  a second (clean-time median 0 or 1, 99.9th percentile 35 to 395). These
  are the lowest minutes of the masked run's coverage (61 and 65 %),
  because no unit claims the tail's events; its TT16 units take about 190
  spikes there.
