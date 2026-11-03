---
title: The median reference and per-wire timing, TTM-001 to TTM-003
updated: 2026-11-02
status: needs_verification
---

# The median reference and per-wire timing

The sorting binaries are bandpassed 300 to 6000 Hz and have the median of
all 64 wires subtracted at each sample (`cmr="global"` in `ttm.sorting`).
That default came in with commit 10fb71e (2026-06-03) on a measured cost
and was not reviewed then; on 2026-11-02 the user chose to keep it. In the
TTM-003 carry-forward bake-off (`investigations/carry-bakeoff/`) the units
of the tetrodes with few units turned out to follow an event on all 64
wires that the median leaves partly in place. This directory asks why, and
what it costs to reference differently.

Reference and ground are shorted together and tied to a cerebellar skull
screw; the headstage is an Open Ephys OEPS-6571 (64 channels), on a
2nd-generation acquisition board for TTM-001 and a 3rd-generation board
for TTM-002 and TTM-003 (from the user; the IMU stream in `settings.xml`
is present only for TTM-002 and TTM-003). In all three `settings.xml` the
Record Node comes before the Channel Map, so CH1 to CH64 are in headstage
order.

## Code

- `cmr_cost.py`: the time of the lazy chain (Zarr read, bandpass,
  reference) on 60 s of the TTM-003 store, with no reference, the global
  median and the cross-tetrode median (spikeinterface's `out_of_group`).
- `sampling_offsets.py`: each wire's gain and time offset against the
  median of all wires, in TTM-003's 48 h session, TTM-001's NOD session and
  TTM-002's SD session, at three to five time points
  (`results/sampling_offsets.csv`).
- `aligned_binary.py`, `aligned_check.py`, `gain_choice.py`,
  `common_events.py`, `event_state.py`, `timing_validation.py`,
  `timing_pooling.py`: the aligned reference, described below.
- In `investigations/carry-bakeoff/`: `reference_check.py` (the event
  before referencing) and `reference_skew.py` (what each reference leaves
  at the few-unit tetrodes' units' spikes, `results/reference_skew.txt`
  there).

Paths are literal tononi-2 paths.

## The event

At the spike times of TTM-003's few-unit tetrodes' units in the
per-tetrode carry (TT1 unit 0, TT2 unit 2, TT3 unit 3, TT13 unit 39, TT14
unit 40, TT15 unit 43, TT16 unit 44), the unreferenced mean has a trough of
-160 to -320 uV on all 64 wires, and the median carries it at the same
size. The median trace has 61 troughs per second below 5 of its MAD (108
uV); 76 to 96 % of these units' spikes fall within 1 ms of one, against 11
to 16 % for well-isolated units. Subtracting the median leaves -7 to -28 uV
and +5 to +29 uV on the unit's own tetrode, where the noise is 4 to 5 uV.
Those remainders are what the carry fits as units. Conjecture: the event
is the screw's own signal, which every wire records through the shared
reference.

## Why the median leaves part of it

Per wire, over 60 s of bandpassed, unreferenced data, a least-squares fit
x = g m + b m', where m is the median of all 64 wires and m' its slope
from sample to sample. g is the wire's gain on the common signal. A wire
that records the signal tau later than the median has x(t) = g m(t - tau),
close to g m - g tau m' while tau is short against the event (here under
35 us against about 0.2 ms), so tau = -b / g. Local spikes are uncorrelated
with m and average out over 1.8 M samples; the median trace's noise is
21.5 uV against 4 to 16 uV on single wires. The median has no sampling time
of its own, so only differences between wires mean anything.

| | TTM-001 NOD | TTM-002 SD | TTM-003 48H |
|---|---|---|---|
| CH1 to CH32, slope (correlation) | -0.97 us per channel (-0.94) | -0.92 (-0.87) | -1.23 (-0.95) |
| CH33 to CH64, slope (correlation) | -0.97 us per channel (-0.97) | -1.04 (-0.94) | -0.90 (-0.73) |
| offset range | -15 to +19 us | -16 to +21 us | -19 to +35 us |
| two blocks 2 min apart agree to (median) | 0.29 us | 0.27 us | 0.31 us |
| largest change of a wire over 48 h | 11 us | 13 us | 13 us |
| gain range | 0.94 to 1.04 | 0.93 to 1.30 | 0.93 to 1.14 |

The offset profiles correlate 0.74 to 0.93 between animals. Within each
bank of 32 channels the offset falls by about 1 us per channel in all
three animals, close to 1/(35 x 30 kHz) = 0.95 us. Conjecture: the
headstage's amplifier chip digitises its channels one after another,
about 1 us apart. On top of that ramp, wires depart from it by animal
(TTM-003 CH33 +12, CH39 +24, CH64 +24 us; TTM-001 CH17 to CH23 odd +5 to
+8; TTM-002 CH11, CH17 and CH19 +11 to +13, CH40 +12) and drift over hours
one tetrode at a time: in TTM-003 TT14's four wires (CH10, 12, 14, 16)
dropped 10 to 13 us between 12 and 24 h, and TT13's (CH2, 4, 6, 8) rose 8
to 12 us over the 48 h. Conjecture: each electrode's impedance with the
headstage's input capacitance delays its signal by a few us, and changes as
the impedance does.

A 300 uV event whose slope reaches tens of uV per us leaves a few to tens
of uV after subtracting a median that is 10 to 30 us out of step.

## What each reference leaves

`reference_skew.py`, five blocks of 60 s of TTM-003's window 48, the mean
on the unit's own tetrode at its accepted spikes (minimum and maximum, uV):

| unit | median of 64 | median of the other 60 | g x median | shifted, then g x median |
|---|---|---|---|---|
| TT1 unit 0 | -21.5 / +6.5 | -22.9 / +7.0 | -19.6 / +4.3 | -8.8 / +2.8 |
| TT3 unit 3 | -14.2 / +22.2 | -15.2 / +23.8 | -9.5 / +22.3 | -3.1 / +4.4 |
| TT13 unit 39 | -28.3 / +20.2 | -29.9 / +21.0 | -28.2 / +9.8 | -7.0 / +1.7 |
| TT14 unit 40 | -14.2 / +5.3 | -15.8 / +5.6 | -12.7 / +4.3 | -3.6 / +1.9 |
| TT15 unit 43 | -10.6 / +29.0 | -11.6 / +30.9 | -9.8 / +29.1 | -2.8 / +5.9 |
| TT16 unit 44 | -10.6 / +28.5 | -11.4 / +30.0 | -7.9 / +27.6 | -1.7 / +3.8 |
| TT2 unit 2 | -7.2 / +20.0 | -8.0 / +21.9 | -7.0 / +19.8 | -4.5 / +14.3 |
| TT6 unit 16 (isolated) | -101.9 / +29.3 | -102.6 / +29.7 | -101.9 / +29.3 | -101.7 / +29.3 |

"Shifted" moves each wire by -tau with spikeinterface's `phase_shift`
(a fractional-sample shift by FFT, the step IBL runs before destriping) and
then subtracts each wire's least-squares gain on the shifted median times
that median. It removes 70 to 90 % of the remainder on six of the seven
units and leaves well-isolated units unchanged; shifting the wrong way
makes the remainders larger. TT2 unit 2 keeps most of its positive lobe.
The cross-tetrode median leaves slightly more than the global one, since
it is computed from the same skewed wires.

The remainders are measured at times chosen because they looked like
spikes, so they overstate what is left.

## Cost

`cmr_cost.py`, one process, 60 s of 64 channels, with the machine under
load: bandpass alone 7.4 s, with the global median 8.3 s, with the
cross-tetrode median 47.4 s. The 48 h global binary took 0.765 h with 96
jobs, most of it writing 1.3 TB. The 2026-06-03 docstring's "~20x" was
against earlier cross-tetrode code that took 64 separate medians, replaced
on 2026-10-19 by spikeinterface's `out_of_group`. The cost of `phase_shift`
in the chain is not measured.

## The aligned reference

`ttm.reference` and `preprocess_for_sorting(cmr="global_aligned")`
estimate each wire's offset and gain per 1800 s window and apply them:
the wire is shifted by its offset, then its gain times the median of the
shifted wires is subtracted. The fit uses only samples within 0.5 ms of
the median's excursions beyond 5 MAD; on all samples, quiet blocks gave
gains up to 1.9, since between the large events the median follows
smaller common signals with another spread across wires. A window's row
is the per-wire median over six blocks of 20 s with at least 3000 such
samples; a single block at the middle of window 1 fell in the 2287 to
2696 s event on seven tetrodes and gave offsets up to 119 us off. Window
4 had no such block and takes window 5's row.

`aligned_binary.py` estimated TTM-003's 96 windows in 12 min
(`results/wire_timing_ttm003.csv`: shifts -40 to +20 us, gains 0.76 to
1.11, median change between windows 0.42 us; over the 48 h the wires of
TT13 and TT14 move up to 11 and 13.5 us, most others under 8) and wrote
the first 12 windows with it in 4.8 min with 64 jobs. `aligned_check.py`
(`results/aligned_check.txt`): at the same spike times the few-unit
tetrodes' units keep 3 to 27 % of their mean's energy, TT4 unit 5 44 %,
well-isolated units 99 to 100 % (cosine 0.997 or more). The carry
bake-off on that binary is in `investigations/carry-bakeoff/README.md`,
"The aligned reference on the hybrid".

`gain_choice.py` (`results/gain_choice.txt`) fits windows 0 and 5 twice,
on all samples of the block with the fewest event samples and on the
event samples of the others, and applies each to the other blocks. Within
0.5 ms of events the all-sample fit leaves more than the global median on
the wires where the two gains differ most (window 5: CH64 52.7 against
24.5 uV, CH22 38.7 against 19.2, CH39 34.4 against 21.2), the event fit
less (22.3, 19.0, 18.7). Over all wires the median RMS at events changes
little (14.2 to 14.1 uV), since local activity and noise dominate it.

## How often the events occur

`common_events.py` (`results/common_events.csv`) counts excursions of the
median beyond 5 of its MAD in the six blocks of each window:

| | TTM-003 48H | TTM-001 NOD | TTM-002 SD |
|---|---|---|---|
| windows | 96 | 96 | 24 |
| events per s, median over windows | 40 | 25 | 27 |
| range | 1.7 to 87 | 0.7 to 53 | 2.6 to 65 |
| events beyond 100 uV per s, median | 26 | 7 | 16 |
| 99th percentile of size | 438 uV | 225 uV | 332 uV |
| windows under 1 event per s | none | 8, 9, 17 | none |

The rate changes by tens per second between neighbouring windows. In
TTM-003 and TTM-001 it is highest through windows 48 to 57, the sleep
deprivation (09:00 to 14:00 on the second day, forced wake with novel
objects and the experimenter's hands in and out of the Faraday cage).

`event_state.py` (`results/event_state.csv`) compares the rate per window
with the synthetic EMG and slow-wave activity over the same blocks:

| | TTM-003 48H | TTM-001 NOD |
|---|---|---|
| Spearman with EMG | +0.88 | +0.88 |
| Spearman with slow-wave activity | -0.83 | -0.84 |
| the same outside the deprivation | +0.86 / -0.84 | +0.88 / -0.82 |
| events per s, lowest to highest quarter of slow-wave activity | 65, 46, 33, 18 | 41, 29, 20, 8 |

The EMG is not an independent measure: it is the correlation of 300 to
600 Hz between single unreferenced wires, which a signal common to all
wires raises. Slow-wave activity is in another band. The rate follows
state outside the deprivation as closely as across it, so the deprivation
has more events because it has more wake. Conjecture: the events are
muscle activity of the neck that the cerebellar screw records, which is
also why the synthetic EMG works on these recordings. Units fitted to the
remainder of the events would fire more in wake than in sleep.

The hybrid's windows 0 to 11 range from 1.7 to 80.7 events per s (median
30); two of the 96 windows have more.

## The estimate on held-out data, three animals

`timing_validation.py` (`results/timing_validation.csv`, per-block fits
in `mp_aligned/timing_validation_blocks.npz` under TTM-003's sortings) fits each of the six blocks of
every window on its own, then applies the per-wire median of the
qualifying blocks among 0, 2 and 4 to blocks 1, 3 and 5. The measure is
each wire's RMS within the event mask, aligned over global; it is diluted
by local activity and noise, as in `gain_choice.py`, so it changes by a
few per cent where the unit-triggered remainders above change several
fold.

| | TTM-003 48H | TTM-001 NOD | TTM-002 SD |
|---|---|---|---|
| windows with no qualifying block | 1 (window 4) | 3 (8, 9, 17) | 0 of 24 |
| a wire's offset between blocks, median spread | 0.29 us | 0.46 us | 0.38 us |
| median wire, own window | 0.948 | 0.984 | 0.958 |
| median wire, estimate of the window before | 0.948 | 0.985 | 0.961 |
| median wire, estimate of 8 windows (4 h) before | 0.948 | 0.985 | 0.955 |

An estimate from 4 h earlier does as well as the window's own. Where one
block decides the estimate it can be far off: in TTM-002 window 84 only
block 0 qualifies, and it gives CH24, CH28 and CH58 gains of 0.15 to 0.39
and offsets of -90 to -180 us where the other blocks give 0.93 to 1.0 and
-10 us; those wires end 2.2 to 2.7 times worse than with the global
median. TTM-003 window 18 leaves CH22, the wire with the 400 Hz comb,
1.5 times worse (block gains 1.0 to 1.7).

`timing_pooling.py` (`results/timing_pooling.csv`) takes the median over
the qualifying blocks of the window and of the two windows on each side
instead:

| | own window: worst wire, largest over windows | pooled over 5 windows | windows with a wire 5 % worse, own / pooled |
|---|---|---|---|
| TTM-003 48H | 1.504 | 1.028 | 6 / 0 |
| TTM-001 NOD | 1.007 | 1.004 | 0 / 0 |
| TTM-002 SD | 2.653 | 1.107 | 3 / 1 |

The median wire is unchanged (0.945 against 0.948, 0.984, 0.959 against
0.958), and every window has an estimate. The TTM-003 rows the 48 h
aligned binary was written with depart from a rolling median over five
windows by at most 5 us and 0.1 of gain (windows 0, 2, 18, 29 and 63).
`estimate_wire_timing_windows(pool_windows=2, min_blocks=3,
if_none="identity")` does this (not the default).

This is not adaptive gain control: the gain is a fixed number per wire
and window, and spikeinterface's `highpass_spatial_filter`, which
assumes a dense line of sites, does not fit tetrodes. Not done: the 48 h
on the aligned binary.
