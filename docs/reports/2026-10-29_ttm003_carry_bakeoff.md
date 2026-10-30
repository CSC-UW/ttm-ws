---
title: TTM-003 48 h, carry-forward matching per tetrode and the rebuild of templates
updated: 2026-10-30
status: needs_verification
---

# TTM-003 48 h: matching per tetrode, and how a carry-forward should rebuild its templates

Run from 2026-10-27 to 2026-10-29 on `TTM-003-48H`. Two questions from the
user: can circus-omp run on each tetrode on its own beat circus-omp on all
channels at once (the June pipeline's matcher), and what is the right way
for the carry-forward to rebuild its templates from one 1800 s window to the
next. The code and the full tables are in `investigations/carry-bakeoff/`
(README there); the 48 h comparisons use
`subjects/ttm-003/scripts/bad_periods_compare.py` and
`tagged_cells_blocks.py`. Measurements are stated as such; sentences marked
conjecture are the analyst's reading and have not been confirmed by the
user.

## 1. Summary

- Three changes, not the matcher, carry most of the gain over June's carry:
  whitening each tetrode's traces across its four wires, admission in noise
  units that does not depend on the template's own scale, and a rebuild that
  blends the template 0.2 toward the window's clean lead spikes. Over 12
  windows (6 h) of hybrid ground truth, per-tetrode matching with them
  reaches F1 0.79 on the units of 6 MAD or less and 0.91 on the larger ones;
  June's carry reaches 0.38 and 0.53.
- With those changes applied to both, joint and per-tetrode matching tie on
  F1 (0.78 and 0.92 joint) and trade recall for precision: per tetrode has
  the higher recall on all 21 ground-truth units and on the driven units'
  half-size in-pulse spikes (0.61 against 0.50), joint the higher precision
  (0.93 against 0.87 on the small units). Per tetrode needs no mask of the
  bad periods and costs about four times joint's matching time.
- Whitening is the largest single gain and helps the joint matcher too. It
  raises the matched-filter signal-to-noise of the seed templates only 1.1
  to 1.5 times; most of the gain is in which false fits survive.
  Conjecture: the false fits it removes are signals common to the four wires.
- Every rebuild recipe did worse than frozen templates while admission was
  measured against the template's scale: the rebuild averages spikes the
  template's own floor has already cut, so templates near the floor grow
  and lose their units (a ratchet). With admission made scale-free, a
  blend of 0.2 toward the mean of the window's lead spikes (outside laser
  pulses and the tetrode's bad periods, neighbour-subtracted snippets)
  matches frozen templates on ground truth and can follow drift; full
  replacement lets a small unit be captured.
- On the real 48 h the per-tetrode carry claims 84.0 % of the detected
  events in clean time (98.9 % of those of 10 MAD or more) against 55.0 %
  (90.4 %) in production, puts 64 spikes in the tetrodes' own bad periods,
  keeps all 68 units present, and keeps each light-driven cell's response on
  one unit for 48 h, where production hands it from unit to unit.
  Production's rebuild also moves its units' spike times along their
  waveforms (by window 94, 9 units by 5 samples or more, the largest by 30
  samples, 1.0 ms); the per-tetrode carry's move by 3 samples at most.
- On the tetrodes with few units the per-tetrode carry has 10 to 100 times
  production's spikes. These units follow events seen on many tetrodes at
  once: at their spike times the busy tetrodes TT4, TT8 and TT12 show
  troughs two to five times their own, and 61 to 94 % of their spikes
  coincide with a spike on another tetrode (chance about 40 %). Joint
  matching's competition across tetrodes kept these out; per tetrode needs
  a test across tetrodes that does the same.
- Open: whether joint matching with the same rules also keeps the driven
  cells on one unit and the spike times on the trough over 48 h (running);
  the relaxed rule inside laser pulses admits some false fits. A whitening
  matrix per window gains nothing over the first 300 s' matrix.

## 2. Terms

- Window: one of the 96 windows of 1800 s of the carry-forward. Templates
  are matched over a window, then rebuilt from what was accepted in it.
- Joint matching: circus-omp over all 64 channels at once, as the June
  pipeline runs it. Per-tetrode matching: circus-omp over one tetrode's four
  channels with only that tetrode's templates, sixteen times per window.
- Matched-filter signal-to-noise (snr): the projection of a spike's
  neighbour-subtracted snippet on the unit-norm template, over the robust
  standard deviation of that projection across 60 s of the window (six
  blocks of 10 s). The neighbour-subtracted snippet is the data around the
  spike with every other fitted spike removed; the code calls it peeled.
- Relative floor: the fit amplitude over the median fit amplitude of the
  unit's lead spikes outside pulses in the window.
- Whitening: per tetrode, a 4x4 matrix W = C^(-1/2) (ZCA) from the
  covariance C of the samples of the first 300 s that lie more than 30
  samples from any crossing of -4 MAD; traces are read as x W and templates
  carried as T W.
- Hybrid: a copy of the first 12 windows of the preprocessed binary with 21
  known units added (section 3).

## 3. Ground truth

`make_hybrid.py` adds 21 units to the first 12 windows (155 GB on the
session's sortings volume). Each is a seed-bank template moved to another
tetrode with its wires permuted, at 5 to 14 MAD; two are mixes with a unit
already on the tetrode (cosine about 0.87). They have rate swings, bursts
whose second to fourth spikes are 0.85, 0.7 and 0.6 of full size, a slow
change of size (down to 0.6) and of wire ratios (up to 20 %), size jitter
and sub-sample timing; two stop firing after window 5. Three answer the real
laser pulses with a full-size spike and then a train at half size, the
pattern of the TT6 cell; the others fall silent through each pulse. Eleven
units are 6 MAD or less ("small") and ten larger.

A found spike is one the unit seeded from the ground-truth unit accepted
within 10 samples (0.33 ms). Precision is the share of that unit's accepted
spikes that have a ground-truth spike within 10 samples.

## 4. What the joint matcher does that per tetrode does not

Read from spikeinterface 0.105.2's circus-omp:
- Candidate times are chosen across the probe: within about 3 ms, only the
  sample with the best template score over all channels is a candidate, and
  a fit is kept if its amplitude is above 0.25 of the template's.
- The rule that stops a chunk after a number of failed fits is shared by the
  whole probe.
- Fits outside the amplitude bounds are still subtracted.

The joint matcher's precision on clean data (`investigations/per-tetrode-matching/`)
therefore comes partly from an early stop and from competition across
tetrodes, which also hide spikes on tetrodes with few or small units. On
TT1 and TT16 joint matching found 0 and 20 %, and 0 and 53 %, of two
ground-truth units' spikes in clean time, where per tetrode found 96 to
98 % of the larger ones.

## 5. Admission

Admission rules were applied to permissive per-tetrode fits on frozen
templates (stage A, window 0 clean apart from laser train 1, window 1 with
the 2287 to 2696 s event on seven tetrodes). A floor as a fraction of the
template (June's a >= 0.8) admits every fit of the right scale, including
background fits of the right footprint; a floor in noise units (snr) with a
cosine gate on the neighbour-subtracted snippet and a 1.5 ms refractory rule
does better on both. On plain traces in window 0, snr >= 4, a >= 0.5 and
cosine >= 0.6 give small-unit recall 0.63 at precision 0.80 against joint's
0.66 at 0.70, and large-unit 0.87 at 0.83 against 0.73 at 0.86.

Per tetrode, a tetrode's own bad periods are dropped from its own fits, and
no mask is needed: counted outside each unit's bad periods, window 1 gives
large-unit recall 0.84 to 0.88 per tetrode against 0.71 joint masked.

## 6. Whitening

| window 0, frozen templates | small recall | small precision | large recall | large precision |
|---|---|---|---|---|
| joint, a >= 0.8 (June) | 0.66 | 0.70 | 0.73 | 0.86 |
| per tetrode, plain, snr >= 4, a >= 0.5, cosine >= 0.6 | 0.63 | 0.80 | 0.87 | 0.83 |
| per tetrode, whitened, snr >= 4.5, a >= 0.5, cosine >= 0.5 | 0.85 | 0.93 | 0.96 | 0.93 |

Whitened, the same rule set at cosine 0.6 instead of 0.5 gives precision
0.98 to 0.99 at lower recall. The seed templates' matched-filter
signal-to-noise rose only 1.1 to 1.5 times. Conjecture: the fits whitening
removes are far-field spikes and artifacts common to the four wires, whose
cosine with a template falls once the common part is scaled down.

The whitening of the first 300 s holds over the 48 h on most tetrodes: on
traces whitened with it, the largest correlation between two wires stays at
0.15 or less, except TT14 (0.28 to 0.29 from window 48 on) and single
windows of TT2, TT5 and TT10 (0.21 to 0.23) (`whitening_drift.py`). Section
10 tests a matrix per window.

## 7. Rebuilding templates

June's carry (joint, masked, mean of 500 spikes accepted at a >= 0.8, full
replacement) fails on the hybrid in three ways:
- Templates of units near the floor grow (GT20 to 1.8 times its injected
  size, GT0 1.6, GT8 1.7): the mean is taken over spikes the floor has
  already cut from below, the larger template raises the floor, and the
  unit's spikes drop out.
- Small units on busy tetrodes are captured by background spikes (GT16 on
  TT10 at precision 0.01 by window 2, with 37,000 spikes).
- Unmasked, the window-1 event rebuilds the event tetrodes' units from
  artifact, and they never fire again.

Every recipe tried with admission tied to the template's scale did worse
than frozen templates: all spikes or lead spikes only, raw or
neighbour-subtracted snippets, mean or median, replacement or blend, laser
pulses and bad periods left out. Two corrections of the scale itself failed:
a truncated-normal fit shrank real units' templates to its clip, since their
signal-to-noise distributions are mixtures, and a band symmetric about the
mode left most real units with none.

With admission made scale-free (snr floor, cosine gate, relative floor
against the median fit), means over the 12 windows:

| carry | small R / P | large R / P | small F1 | large F1 | half-size laser | event, large |
|---|---|---|---|---|---|---|
| June: joint, masked | 0.38 / 0.38 | 0.46 / 0.64 | 0.38 | 0.53 | 0.02 | 0.56 |
| joint, masked, whitened, June's rebuild | 0.62 / 0.79 | 0.78 / 0.95 | 0.69 | 0.86 | 0.02 | 0.77 |
| per tetrode, plain, frozen | 0.54 / 0.67 | 0.80 / 0.78 | 0.60 | 0.79 | 0.53 | 0.70 |
| per tetrode, whitened, frozen | 0.72 / 0.89 | 0.93 / 0.91 | 0.79 | 0.92 | 0.60 | 0.90 |
| per tetrode, whitened, blend 0.2 | 0.72 / 0.87 | 0.93 / 0.90 | 0.79 | 0.91 | 0.61 | 0.90 |
| per tetrode, whitened, blend 0.5 | 0.71 / 0.81 | 0.93 / 0.89 | 0.76 | 0.91 | 0.61 | 0.90 |
| per tetrode, whitened, full replacement | 0.70 / 0.74 | 0.93 / 0.89 | 0.72 | 0.91 | 0.60 | 0.89 |

The rebuild that holds is a blend of 0.2 toward the mean of the
neighbour-subtracted snippets of the window's lead spikes (no spike of the
unit in the 10 ms before), outside laser pulses (onset to offset plus 100
ms) and the tetrode's bad periods, at most 500 of them and only when there
are 100. Most units' templates stay at cosine 0.98 or more to the injected
ones under every recipe; the differences are in a few small units. GT20 (6
MAD, shrinking to 0.6 of its size) is captured under full replacement
(cosine 0.62, precision 0.02 by window 6) and not under blend 0.2. The
hybrid's drift did not need close tracking, so frozen templates do as well
as the blend over 6 h; over 48 h the real units drift more, and frozen
templates were not run there.

## 8. Laser pulses

Inside each pulse (onset to offset plus 100 ms), units that keep firing
through the pulses at twice their rate outside get a relaxed rule (snr >=
2.5, relative floor 0.3, cosine >= 0.5). On window 0, plain traces:

| in-pulse rule | GT6 half-size | GT13 half-size | false in-pulse spikes on silent ground-truth units |
|---|---|---|---|
| as outside | 0.29 | 0.13 | 2.9 /s |
| relaxed, every unit | 0.63 | 0.60 | 11.3 /s |
| relaxed, units that fire through the pulses | 0.63 | 0.60 | 3.6 /s |

The test for "fires through the pulses" first counted both rates under the
outside rule and on the real 48 h missed the TT6 cell, whose half-size
train is below that rule. It now counts both under the relaxed rule. That
version flags one silent ground-truth unit (GT1) in some windows of the
hybrid; rerun under it, the per-tetrode carry's ground-truth measures do not
change (section 10).

## 9. The 48 h

Two carries ran over the 96 windows of the preprocessed binary from
mp_baseline's seed bank, with June's re-seed every 12 windows from the
MountainSort5 sorts June's runs had made: per tetrode, whitened, blend 0.2,
unmasked (68 units after assembly, 161.6 M spikes); and joint, masked,
whitened, June's rebuild (100 units, 90.5 M). Production
(`mp_baseline_masked/assembled_reseed_c12`) has 73 units and 82.8 M spikes.

| | production | joint, whitened | per tetrode, whitened |
|---|---|---|---|
| single units | 73 | 100 | 68 |
| units present after window 1 | 71 | 88 | 68 |
| spikes in the tetrode's own bad periods | 0 | 77,897 | 64 |
| units with most spikes in bad periods | 0 | 16 | 0 |
| detected events in clean time claimed by a unit | 55.0 % | | 84.0 % |
| the same, events of 10 MAD or more | 90.4 % | | 98.9 % |

Detected events are negative peaks beyond 5.5 MAD on the tetrode
(production's `spike_coverage.npz`); one is claimed when a single unit of
that tetrode has a spike within 0.5 ms. The mask puts noise independent
across wires in place of a tetrode's bad periods; whitened with a matrix
built for correlated wires, that noise has 1.6 to 5.0 times the variance of
real whitened noise on TT16 and 2.1 to 5.8 on TT14 (`shared_events.py`).
The joint whitened carry with June's rules fits it there and re-seeds units
in TT16's 12,399 s of bad time. Per tetrode reads the real traces and drops
its fits in the tetrode's own bad periods.

Coverage at 10 MAD or more per tetrode, production against per tetrode:
TT1 26.7 and 90.3 %, TT2 52.7 and 95.0 %, TT3 40.9 and 88.4 %, TT14 43.8
and 81.0 %, TT16 36.4 and 86.0 %; on the eight busy tetrodes 83 to 99 %
against 96 to 99.9 %, and TT9 80.3 against 60.8 %.

### Spike times

Against each unit's mean waveform in its first window, by window 94
production has 9 units whose spike times have moved 5 samples or more along
the waveform (90th percentile 6 samples, largest 30), joint whitened 8 (3,
40), per tetrode none (1, 3). Production's TT4 unit 5 sits 11 samples (0.4
ms) off from window 24 on. Both carries that drift rebuild templates June's
way, and nothing in that rebuild keeps the waveform in place. An earlier
version of these numbers measured the trough's offset from sample 30, which
misreads positive-going units and overstated every carry's drift (production
11 units, joint whitened 17, per tetrode 4). Latencies read
from production's spike times carry this drift, and comparisons between
sortings at 0.33 ms undercount their shared spikes.

### Light-driven cells

The units that answer the pulse onsets (1 to 10 ms), per block:

| tetrode | production | joint, whitened | per tetrode, whitened |
|---|---|---|---|
| TT4 | 6 and 7 in blocks 1 and 2 (310, 205 Hz), then 5 (124 to 203) | 6 to block 27, then re-seeded 2010 | 7 in all 24 blocks (187 to 425); re-seeded 2007 also from block 28 |
| TT10 | 24 in blocks 1 to 4, then 26, then re-seeded 52 | 25 in all blocks (57 to 172) | 25 in all blocks (100 to 224) |
| TT11 | 28, from 159 down to under 30 by block 37 | 28 to block 3, then 29 | 28 in all blocks (119 to 222) |
| TT6 | 16 at 19 and 13 Hz in blocks 1 and 2, then 0 | 16 at 31 in block 1 | 16 at 176 to 238 in blocks 1 to 4, then 7 to 76 |

In window 0 the three sortings' units agree by id, and on TT4 the mean of
the spikes after the onsets matches unit 7's spontaneous waveform in every
sorting (cosine 0.75 to 0.89). Production's unit 5 has 9 such spikes in
window 0 and carries the response from block 3, while its spike times move
0.4 ms off its trough. Conjecture: June's rebuild let unit 5's template
move onto the driven cell.

The per-tetrode in-pulse rates use the relaxed rule and can include false
fits. The template lookup of `tagged_cells_blocks.py` picks other units than
these in the whitened sortings, so its `tagged` summaries understate them.

### The tetrodes with few units

Per tetrode has 10 to 100 times production's spikes on them (TT1 33,747
against 2,840,246; TT3 11,391 against 1,300,739; TT14 33,531 against
610,502; TT16 134,280 against 3,557,807). TT1 unit 0 goes from 0.11 Hz at a
median 10.1 MAD to 12.3 Hz at 5.3 MAD, median cosine 0.82 against 0.55.
The procedure does not make their refractory pairs: on the hybrid in
window 0 the 1.5 ms rule removes 0.00 to 0.04 of the fits of the
ground-truth units put on TT1, TT3, TT13, TT14 and TT16, and 0.25 to 0.38 of
the real units' fits on the same tetrodes (most units of the busy tetrodes:
0.01 to 0.08). On the real data the removed fits lie 0.3 to 1.5 ms from the
kept spike, not under 0.3 ms as double fits would; the units' small spikes
have the waveform and wire ratios of their large ones (cosine 0.96 to
0.99), and 0.03 to 0.34 of even their spikes of snr 12 or more are removed
(`refractory_pairs.py`).

They are not local cells (`shared_events.py`, window 48). At their spike
times the mean trace has troughs on other tetrodes as large as or larger
than on their own: TT14 unit 40 is -14 uV on TT14 and -69, -71 and -52 uV on
TT4, TT8 and TT12; TT1 unit 0 is -21 on TT1, -44 on TT4 and -39 on TT8.
Well-isolated units are about -1 uV off their own tetrode. 61 to 94 % of
these units' spikes have an accepted spike on another tetrode within 0.2 ms
(chance 39 to 46 %), most often on the other few-unit tetrodes; the
well-isolated units are at chance. These units follow events seen on many
tetrodes at once, and the joint matcher's competition across tetrodes kept
most of them out. Production's units on these tetrodes show the same
footprint across tetrodes, with far fewer spikes.

### Cost

Per-tetrode matching took 4.6 to 5.3 min per window with the machine to
itself, against 1.2 min joint; 16.6 h for the 48 h as run, alongside other
jobs.

## 10. Follow-up arms

### Matcher against rules

Joint matching, masked, on whitened traces with every rule of the
per-tetrode carry (matching floor, admission, in-pulse rule, the tetrode's
own bad periods dropped, blend 0.2), against the per-tetrode carry, means
over the 12 windows of the hybrid:

| carry | small R / P | large R / P | small F1 | large F1 | half-size laser | burst tails, large | event, large |
|---|---|---|---|---|---|---|---|
| joint, June's rules and rebuild | 0.62 / 0.79 | 0.78 / 0.95 | 0.69 | 0.86 | 0.02 | 0.39 | 0.77 |
| joint, per-tetrode rules | 0.67 / 0.93 | 0.90 / 0.94 | 0.78 | 0.92 | 0.50 | 0.89 | 0.88 |
| per tetrode | 0.72 / 0.87 | 0.93 / 0.90 | 0.79 | 0.91 | 0.61 | 0.93 | 0.90 |
| per tetrode, cosine >= 0.55, full replacement | 0.65 / 0.81 | 0.90 / 0.94 | 0.72 | 0.92 | 0.61 | | 0.88 |

The stricter gate was run only with full replacement; it lands on joint's
point for the large units, and for the small units the rebuild confounds
the comparison.

Per tetrode has the higher recall on every ground-truth unit (by up to 0.11,
most on GT5 on TT9, GT18 on TT12 and GT15 on TT7). Joint keeps the higher precision, most on GT1 on TT1 (0.88 against
0.75), GT4 on TT9 (0.94 against 0.77) and the two units mixed with a
neighbour. Conjecture: the competition across tetrodes still removes false
fits that the rules let through (section 9 shows such fits on the real
data), and with them true small spikes.

### Driven test

Rerun with the driven test that counts both rates under the relaxed rule,
the per-tetrode carry's ground-truth measures are unchanged to the third
decimal. The test flags GT1 in some windows.

### Whitening per window

A matrix per window (30 blocks of 10 s spread over the window, clear of
the tetrode's bad periods and the pulses) instead of the first 300 s':
- On the hybrid, small-unit precision 0.84 against 0.87 and F1 0.77
  against 0.79; large units and the event unchanged; half-size laser 0.64
  against 0.61.
- On the real windows 48 to 59, from the 48 h carry's bank after window 47:
  2 to 5 % more spikes on the busy tetrodes at a lower median cosine (TT8
  0.684 against 0.717), 16 % fewer on TT14, whose residual correlation is
  highest, with fewer refractory removals (0.17 against 0.21). The two
  runs share 55 to 95 % of each other's spikes per tetrode and their
  templates stay at cosine 0.98 or more.

Neither asks for re-whitening every window. Re-whitening a tetrode when its
residual correlation crosses a threshold is not tested.

### Still running

Joint matching with the per-tetrode rules over the 48 h
(`R48_W_joint_masked_s45_c50_lnb_b2`), to see whether the driven cells stay
on one unit, the spike times stay on the trough and the few-unit tetrodes'
units stay small under the joint matcher.

## 11. What this suggests

These are the analyst's reading, for the user to confirm.
- Whitening per tetrode, scale-free admission, the refractory rule and the
  blend 0.2 rebuild belong in the carry whichever matcher runs it.
- Between the matchers the choice is recall and no mask against precision
  and speed. The per-tetrode carry's real-data advantages so far (driven
  cells, spike times) may come from the rebuild rather than the matcher;
  the 48 h joint run decides it.
- Per tetrode needs a test across tetrodes in place of the joint
  matcher's competition: drop a fit when the same event is larger on
  another tetrode, or when many tetrodes fit it at once. On the hybrid this
  should recover joint's precision on GT1 and GT4 without its loss of
  recall; not yet tested.
- One whitening for clustering and matching (one matrix for the
  MountainSort5 seed and re-seed sorts and for the carry) is a follow-on.
