---
title: Carry-forward bake-off on hybrid ground truth, TTM-003 48 h
updated: 2026-11-02
status: needs_verification
---

# Carry-forward bake-off

Two questions on the first 12 windows (6 h) of TTM-003's 48 h recording:
can circus-omp run on each tetrode on its own match or beat circus-omp on all
channels at once (June), and how should a carry-forward rebuild its
templates from one window to the next? The per-tetrode carries of
`investigations/per-tetrode-matching/` ran away on tetrodes with few or small
units. That README names two things not tried: a floor in noise units
instead of a fraction of the template, and re-estimation from full-size
spikes only.

## Ground truth

`make_hybrid.py` adds 21 known units to a copy of the first 12 windows of
`mp_baseline/binary` (`carry_bakeoff/v2/hybrid_binary` on the session's
sortings, 155 GB).
- The units are seed-bank templates moved to other tetrodes with their wires
  permuted, at 5 to 14 MAD; two are mixes with a destination unit, at cosine
  about 0.87.
- Each has rate swings, bursts at 0.85, 0.7 and 0.6 of full size, a slow
  change of size and wire ratios, size jitter and sub-sample timing.
- Two stop firing after window 5.
- Three answer the real laser pulses with a full-size spike, then a train at
  half size (the TT6 cell's pattern). The other units fall silent through each
  pulse.

`gt.npz` and `gt_plan.json` record every spike, its class and size factor,
and the injected template per window. A first version (v1) had driven units
of 5 and 6 MAD, whose half-size spikes no rule could find; it was replaced.

## Code

- `_bake.py`: the carry loop.
  - Matching: circus-omp, joint or per tetrode, with a per-template amplitude
    floor.
  - Per spike: fit amplitude; matched-filter signal-to-noise (projection on
    the unit-norm template over that projection's noise, measured per
    template and window); cosine of the raw snippet and of the peeled
    snippet (residual plus the spike's own fit); residual energy; trough.
  - Admission rules and re-estimation recipes, with optional spatial
    whitening per tetrode.
- `run_arm.py`: the arms. Every arm starts from the same bank: mp_baseline's
  seed bank plus one template per ground-truth unit, from its window-0
  spikes.
- `admission.py`, `admission_units.py`: stage A. Rules applied to permissive
  per-tetrode fits on frozen templates, against the joint runs.
- `score.py`, `compare_carries.py`: stage B. Carries over 12 windows, scored
  per window on the ground truth and on proxies for the real units.

Paths are literal tononi-2 paths.

## Results on the 12 windows

Numbers are measurements on the hybrid binary. Recall is the share of a
ground-truth unit's spikes that the unit seeded from it found within
0.33 ms; precision is the share of that unit's spikes that are its own.
"Small" pools the units with a trough of 6 MAD or less (11 units), "large"
the others (10). Bulk products are under
`sortings_seed42_pcafix/carry_bakeoff/v2/` on the session's volume, figures in
its `figures/` (`ladder.png`, `template_scale.png`, `laser.png`).

### Stage A: what each matcher admits on frozen templates

`admission.py`, window 0 (clean apart from laser train 1) and window 1 (the
2287 to 2696 s event on seven tetrodes). Per-tetrode rules are applied to the
permissive per-tetrode fits; "snr" is the matched-filter signal-to-noise, `a`
the fit amplitude, `cos_pl` the cosine of the peeled snippet, and every
per-tetrode rule includes the refractory rule (1.5 ms).

| window 0 | small recall | small precision | large recall | large precision | burst tails, large | half-size laser, large |
|---|---|---|---|---|---|---|
| joint, a >= 0.8 (June) | 0.66 | 0.70 | 0.73 | 0.86 | 0.40 | 0.05 |
| joint masked | 0.66 | 0.70 | 0.73 | 0.87 | 0.40 | 0.05 |
| per tetrode, a >= 0.8 | 0.82 | 0.52 | 0.84 | 0.82 | 0.52 | 0.06 |
| per tetrode, snr >= 4, a >= 0.5, cos_pl >= 0.6 | 0.63 | 0.80 | 0.87 | 0.83 | 0.72 | 0.21 |
| per tetrode, whitened, snr >= 4.5, a >= 0.5, cos_pl >= 0.5 | 0.85 | 0.93 | 0.96 | 0.93 | 0.89 | 0.62 |
| per tetrode, whitened, the same with cos_pl >= 0.6 | 0.68 | 0.99 | 0.89 | 0.98 | 0.74 | 0.24 |

| window 1 | small recall | small precision | large recall | large precision | recall on TT5 to TT12 inside the event |
|---|---|---|---|---|---|
| joint | 0.43 | 0.17 | 0.54 | 0.61 | 0.08 |
| joint masked | 0.46 | 0.52 | 0.63 | 0.81 | 0.64 |
| per tetrode, a >= 0.8 | 0.72 | 0.09 | 0.82 | 0.25 | 0.75 |
| per tetrode, snr >= 4, a >= 0.5, cos_pl >= 0.6 | 0.41 | 0.43 | 0.82 | 0.76 | 0.53 |
| per tetrode, whitened, snr >= 4.5, a >= 0.5, cos_pl >= 0.5 | 0.65 | 0.72 | 0.93 | 0.85 | 0.77 |

On tetrodes with one or two units the joint matcher finds few of the
ground-truth spikes even in clean time (TT1: 0 and 20 %; TT16: 0 and 53 %),
where per tetrode finds 96 to 98 % of the larger ones
(`admission_units.py`). Matching per tetrode with its own bad periods dropped
needs no mask: counted outside each unit's own bad periods, window 1 gives
large-unit recall 0.84 to 0.88 per tetrode against 0.71 joint masked
(`admission_badperiods.py`).

Spatial whitening per tetrode (a 4x4 ZCA matrix from the spike-excised noise
of the first 300 s) raised precision more than any rule did. The
matched-filter signal-to-noise of the seed templates rose only 1.1 to 1.5
times on Gaussian-noise grounds (`whitening_drift.py` and the shrinkage
estimates), so the gain is mostly in what the false fits were. Conjecture:
they were signals common to the four wires (far-field spikes, artifacts),
which whitening removes.

### June's carry on the ground truth

Over 12 windows joint masked matching with June's re-estimation (mean of 500
random spikes taken at a >= 0.8) lost ground-truth recall from 0.66 to 0.30
(small) and 0.73 to 0.24 (large), and precision with it. Three mechanisms,
per unit (`score.py` tables):

- Templates of units near the floor creep up: the rebuild averages spikes
  that the template's own floor already cut below, the larger template
  raises the floor, and the unit's spikes drop out (GT20 to 1.8 times its
  injected size, GT0 1.6, GT8 1.7).
- Small units on busy tetrodes are captured by background spikes (GT16 on
  TT10 at precision 0.01 by window 2 with 37,000 spikes).
- Unmasked, the window-1 event rebuilt the templates of the event tetrodes'
  units from artifact (2 to 13 times their size); they took nothing again
  and were never rebuilt (recall 0 from window 2 on).

### Re-estimation

Every rebuild recipe tried (all spikes or lead spikes only, raw or peeled
snippets, mean or median, replace or blend, laser pulses and bad periods
left out) did worse than freezing the templates when admission depended on
the template's scale.
- A truncated-normal correction of the scale shrank real units' templates to
  its clip, because their signal-to-noise distributions are mixtures.
- A band symmetric about the mode of the signal-to-noise left most real
  units with no band.

Admission was then made independent of the template's scale. The relative
floor is measured against the median fit amplitude of the unit's lead
spikes in the window, alongside the matched-filter signal-to-noise and the
peeled cosine. Means over the 12 windows (`compare_carries.py`,
`results/carries_summary_final.csv`):

| arm | small R / P | large R / P | small F1 | large F1 | half-size laser | event recall, large |
|---|---|---|---|---|---|---|
| joint masked, June | 0.38 / 0.38 | 0.46 / 0.64 | 0.38 | 0.53 | 0.02 | 0.56 |
| joint masked, whitened, June | 0.62 / 0.79 | 0.78 / 0.95 | 0.69 | 0.86 | 0.02 | 0.77 |
| per tetrode, plain, frozen | 0.54 / 0.67 | 0.80 / 0.78 | 0.60 | 0.79 | 0.53 | 0.70 |
| per tetrode, plain, rebuilt (lead, no laser, no bad periods, peeled) | 0.46 / 0.42 | 0.72 / 0.68 | 0.43 | 0.69 | 0.53 | 0.63 |
| per tetrode, whitened, frozen | 0.72 / 0.89 | 0.93 / 0.91 | 0.79 | 0.92 | 0.60 | 0.90 |
| per tetrode, whitened, blend 0.2 | 0.72 / 0.87 | 0.93 / 0.90 | 0.79 | 0.91 | 0.61 | 0.90 |
| per tetrode, whitened, blend 0.5 | 0.71 / 0.81 | 0.93 / 0.89 | 0.76 | 0.91 | 0.61 | 0.90 |
| per tetrode, whitened, full replace | 0.70 / 0.74 | 0.93 / 0.89 | 0.72 | 0.91 | 0.60 | 0.89 |
| per tetrode, whitened, confident spikes, blend 0.5 | 0.70 / 0.77 | 0.93 / 0.91 | 0.73 | 0.92 | 0.60 | 0.90 |

The per-tetrode rows use snr >= 4 (plain) or 4.5 (whitened), cos_pl >= 0.6
(plain) or 0.5 (whitened), a >= 0.5 of the median fit, refractory 1.5 ms,
fits in the tetrode's own bad periods dropped, and the in-pulse rule below.

On the whitened carries most units' templates stay at cosine 0.98 or more to
the injected ones under every recipe; the difference between recipes is in a
few units. GT20 (6 MAD, shrinking to 0.6 of its size) is captured under full
replacement (cosine 0.62, precision 0.02 by window 6) and not under blend 0.2
or frozen. The ground-truth drift (size down to 0.6, a step to 0.75, wire
gains by up to 20 % over 6 h) did not need the templates to track it
closely; real units over 48 h may.

### Laser pulses

Inside each delivered pulse (onset to offset plus 100 ms), units that fire
through the pulses at twice their rate outside them get a relaxed rule (snr >=
2.5, a >= 0.3 of the median fit, cos_pl >= 0.5). On the frozen-template fits
of window 0 (`pulse_rules.py`), plain space:

| in-pulse rule | GT6 half-size | GT13 half-size | false in-pulse spikes on silent ground-truth units | real units' in/out rate, median |
|---|---|---|---|---|
| as outside | 0.29 | 0.13 | 2.9 /s | 0.31 |
| relaxed, every unit | 0.63 | 0.60 | 11.3 /s | 2.46 |
| relaxed, units that fire through the pulses | 0.63 | 0.60 | 3.6 /s | 0.31 |

The test first counted both rates under the standard rule. On the real 48 h
it then missed the TT6 cell (29 Hz outside the pulses, 1.6 times that inside
under the standard rule, its half-size train below the standard floor). The
test now counts both rates under the relaxed rule. The 12-window carries ran
with the first version.

### Real units in the same carries

The real units cannot be scored against ground truth. Refractory conflicts
are the one direct sign. Per accepted spike, pairs of a unit's candidate
spikes 0.5 to 1.5 ms apart, before the per-tetrode carries' refractory rule
removes one of each pair, over the 12 windows (unit-windows with 200 spikes
or more):

| arm | median | 90th percentile |
|---|---|---|
| joint masked, June | 0.007 | 0.033 |
| joint masked, whitened | 0.010 | 0.075 |
| per tetrode, whitened, frozen | 0.025 | 0.191 |
| per tetrode, whitened, blend 0.2 | 0.027 | 0.210 |

The high ratios sit on the single units of the few-unit tetrodes (TT16
unit 44, TT15 units 42 and 43, TT13 unit 39, TT1 unit 0, TT14 unit 40, TT3
unit 3, TT2 unit 1). These already have 8 to 18 % under whitened joint
matching; per tetrode doubles both their spike counts and their ratios
(TT15 unit 42: 151,510 spikes joint, 1,615,428 per tetrode). The
ground-truth units on the same tetrodes keep precision 0.9 to 1.0 per
tetrode. Conjecture: those real units are multi-unit, and per-tetrode
matching gives them more of the activity they already collect.

### Matcher against recipe

`W_joint_masked_s45_c50_lnb_b2` runs joint matching, masked, on whitened
traces with every rule of the per-tetrode blend 0.2 carry (matching floor,
admission, in-pulse rule, the tetrode's own bad periods dropped, blend 0.2).
`W_s45_c50_lnb_b2_dt2` is the per-tetrode carry again under the driven test
that counts both rates under the relaxed rule, and `W_s45_c50_lnb_b2_rw` the
same with a whitening matrix per window. Means over the 12 windows
(`results/carries_summary_followup.csv`):

| arm | small R / P | large R / P | small F1 | large F1 | half-size laser | burst tails, large | event recall, large |
|---|---|---|---|---|---|---|---|
| joint masked, whitened, June's rebuild | 0.62 / 0.79 | 0.78 / 0.95 | 0.69 | 0.86 | 0.02 | 0.39 | 0.77 |
| joint masked, whitened, per-tetrode rules | 0.67 / 0.93 | 0.90 / 0.94 | 0.78 | 0.92 | 0.50 | 0.89 | 0.88 |
| per tetrode, whitened, blend 0.2 | 0.72 / 0.87 | 0.93 / 0.90 | 0.79 | 0.91 | 0.61 | 0.93 | 0.90 |
| the same, current driven test | 0.72 / 0.87 | 0.93 / 0.90 | 0.79 | 0.91 | 0.61 | 0.93 | 0.90 |
| the same, whitening per window | 0.71 / 0.84 | 0.93 / 0.89 | 0.77 | 0.91 | 0.64 | 0.93 | 0.90 |
| per tetrode, cos_pl >= 0.55, full replacement | 0.65 / 0.81 | 0.90 / 0.94 | 0.72 | 0.92 | 0.61 | | 0.88 |

With the same rules, joint matching gives up recall for precision against
per tetrode, and F1 is the same. Per tetrode has the higher recall on all
21 ground-truth units (by 0.0001 to 0.11) and on the half-size in-pulse
spikes of all three driven units; joint has the higher precision, most on
GT1 on TT1 (0.88 against 0.75), GT4 on TT9 (0.94 against 0.77) and the two
mixed units. A stricter gate per tetrode (cos_pl >= 0.55, run only with
full replacement) lands on joint's point for the large units; for the small
units the comparison is confounded by the rebuild. Most of the gain over June's carry comes from
whitening, scale-free admission and the rebuild, which the joint matcher
takes as well as per tetrode does. The driven test that counts both rates
under the relaxed rule leaves the ground-truth measures unchanged; it flags
GT1, which is silent in pulses, in some windows.

## The 48 h

Two carries ran over the 96 windows of `mp_baseline/binary`
(`run_arm.py --bank <BAKE_ROOT>/bank0.npz --windows 0-95 --root
<BAKE_ROOT>/real48`), from mp_baseline's seed bank, with June's re-seed
every 12 windows read from the MountainSort5 sorts June's runs had already
made:
- `R48_W_s45_c50_lnb_b2`: per tetrode, whitened, blend 0.2, the rules of the
  12-window carries, unmasked, the driven test counting both rates under the
  relaxed rule;
- `R48_W_joint_masked_june`: joint, masked, whitened, June's re-estimation.

`export_sorting.py` assembled each into a sorting under the session's
`sortings_seed42_pcafix/` (`carry_R48_pt/assembled_reseed_c12`, 68 units and
161.6 M spikes; `carry_R48_joint/assembled_reseed_c12`, 100 units and 90.5 M
spikes). Production is `mp_baseline_masked/assembled_reseed_c12` (73 units,
82.8 M spikes). Re-seeding added 23 units per tetrode and 55 joint. Matching
took 4.6 to 5.3 min per 1800 s window per tetrode when the run had the
machine to itself and 1.2 min joint.

`subjects/ttm-003/scripts/bad_periods_compare.py` (tags `carry48_pt` and
`carry48_jw_pt`, results in `subjects/ttm-003/results/`):

| | production | joint, whitened | per tetrode, whitened |
|---|---|---|---|
| single units | 73 | 100 | 68 |
| units present after window 1 | 71 | 88 | 68 |
| spikes in the tetrode's own bad periods | 0 | 77,897 | 64 |
| units with most spikes in bad periods | 0 | 16 | 0 |
| clean events over threshold claimed by a unit | 55.0 % | | 84.0 % |
| the same, events of 10 MAD or more | 90.4 % | | 98.9 % |

The joint whitened carry's bad-period spikes are on TT16 (76,830) and TT14
(1,067). The mask replaces a tetrode's four channels inside its bad periods
with noise independent across wires at each wire's clean level; whitened
with a matrix built for correlated wires, that noise has 1.6 to 5.0 times
the variance of real whitened noise on TT16 and 2.1 to 5.8 on TT14, whose
wires are the most correlated (`shared_events.py`), and the joint whitened
carry with June's rules fits it and re-seeds units in it. The per-tetrode
sorting is not masked. Its
refractory measures in that script read 0 by construction, since its
refractory rule removes one spike of each pair.

On the tetrodes with few units the per-tetrode sorting has 10 to 100 times
production's spikes (TT1: 33,747 against 2,840,246; TT3: 11,391 against
1,300,739; TT14: 33,531 against 610,502; TT16: 134,280 against
3,557,807; TT13: 246,950 against 2,313,880). TT1 unit 0 goes from 0.11 Hz
at a median 10.1 MAD to 12.3 Hz at 5.3 MAD, with median cosine 0.82 against
0.55.

`refractory_pairs.py` asks whether the per-tetrode procedure makes
refractory pairs or the units have them. On the hybrid in window 0 the rule
removes 0.00 to 0.04 of the fits of the ground-truth units put on TT1, TT3,
TT13, TT14 and TT16, and 0.25 to 0.38 of the fits of the real units on the
same tetrodes (TT1 unit 0, TT3 unit 3, TT13 unit 39, TT14 unit 40, TT16
unit 44; on the eight busy tetrodes most real units lose 0.01 to 0.08). On
the real 48 h the removed fits lie 0.3 to 1.5 ms from the kept spike, spread
over that range rather than under 0.3 ms as double fits of one event would
be. With `--shapes`, these units' spikes of snr 4.5 to 6 have the same
waveform as their spikes of snr 12 or more (cosine 0.96 to 0.99 on plain
traces, the same ratios between wires), and even among the latter 0.03 to
0.34 are removed.

`shared_events.py` looks across tetrodes in window 48. At these units'
spike times the mean trace has a trough on other tetrodes as large as or
larger than on their own: TT14 unit 40 is -14 uV on TT14 and -69, -71 and
-52 uV on TT4, TT8 and TT12; TT1 unit 0 -21 on TT1, -44 on TT4, -39 on TT8.
Well-isolated units (TT4 unit 4, TT6 units 12 and 16, TT10 unit 26) are
about -1 uV off their own tetrode. 61 to 94 % of the few-unit tetrodes'
units' spikes have an accepted spike on another tetrode within 0.2 ms
(chance 39 to 46 %), mostly on the other few-unit tetrodes (TT15 unit 43
with TT16 at 0.86, TT14 unit 40 with TT13 at 0.71); the well-isolated
units are at chance. These units follow events seen on many tetrodes at
once, larger on the busy ones, and not local cells; the same footprint at
every size and the pairs under 1.5 ms come with that. The earlier
conjecture that each is several local cells with one footprint is
replaced by this. Production's units there are small; joint matching with
the per-tetrode rules fits these events too (TT1 1.60 M spikes against
2.84 M per tetrode and 33,747 in production), so June's rules, not the
joint matcher, kept most of them out.

Before the reference is subtracted these events are on every wire
(`reference_check.py`: five blocks of 60 s of window 48, bandpassed 300 to
6000 Hz from the Zarr store without the reference). At the spike times of
TT1 unit 0, TT2 unit 2, TT3 unit 3, TT13 unit 39, TT14 unit 40, TT15 unit
43 and TT16 unit 44 the mean has a trough of -160 to -320 uV on all 64
wires, the median wire within 7 % of the deepest, and the median of all
wires, which the sorted binary subtracts, carries it at the same size.
What the subtraction leaves on the unit's own tetrode is -7 to -28 uV and
+5 to +29 uV; TT16 unit 44's positive-going waveform is that remainder.
The median trace has 61 troughs per second below 5 of its MAD (108 uV),
and 76 to 96 % of these units' spikes (68 % of TT4 unit 5's) fall within
1 ms of one, against 11 to 16 % (chance) for well-isolated units.
Subtracting each wire's least-squares gain on the median (0.94 to 1.03)
times the median leaves the remainders nearly as they are. These
remainders are measured at times chosen because they looked like spikes,
so they overstate what is left. Reference and ground are shorted together
on a cerebellar skull screw (from the user); conjecture: the event is the
screw's own signal.

The median leaves part of it because each wire records it at its own
time offset from the median (`reference_skew.py`,
`results/reference_skew.txt`; across animals,
`investigations/common-reference/`): about 1 us later per channel number
within each bank of 32, plus electrode-specific offsets, -19 to +35 us in
all. Shifting each wire by its offset with spikeinterface's `phase_shift`
and subtracting its gain times the shifted median removes 70 to 90 % of
the remainder on TT1, TT3, TT13, TT14, TT15 and TT16 (TT16 unit 44 from
-10.6 and +28.5 uV to -1.7 and +3.8) and leaves well-isolated units
unchanged; TT2 unit 2 keeps most of its positive lobe. The median of the
other 60 wires leaves slightly more than the global median. Every result
here on the few-unit tetrodes is on the median-referenced binary.

`joint_competition.py` measures what joint matching's selection costs on
the hybrid. In each pass over a 1 s chunk circus-omp takes, within 3 ms
either way, only the best template score over the whole probe as a
candidate; a template on another tetrode is unchanged by that fit, so a
coincident spike there can be taken in a later pass, and is lost only if
the chunk stops first (five passes in a row with no new fit inside the
amplitude bounds). 98 % of ground-truth spikes have an accepted spike of
the per-tetrode arm on another tetrode within 3 ms. With the same rules,
joint matching's recall is within 0.01 of per tetrode's when none of those
is larger, and lower by 0.05 when one is (large units 0.84 against 0.89,
small 0.60 against 0.65); GT5 on TT9 drops from 0.59 to 0.48 and GT18 on
TT12 from 0.76 to 0.65.

### Spike times against the waveform

`stamp_drift.py` takes the mean of up to 300 of each unit's spikes on the
plain binary, with the spike time at sample 30, in windows 0, 24, 48, 72 and
94, and finds the lag that best aligns it (all four wires) with the unit's
mean in its first window (`results/stamp_drift.csv`):

| units shifted 5 samples or more | window 24 | window 48 | window 72 | window 94 |
|---|---|---|---|---|
| production | 1 | 2 | 7 | 9 |
| joint, whitened | 1 | 1 | 3 | 8 |
| joint, per-tetrode rules | 0 | 0 | 0 | 0 |
| per tetrode, whitened | 0 | 0 | 0 | 0 |
| per tetrode, re-centred | 1 | 4 | 6 | 6 |

In window 94 the 90th percentile of the shift is 6 samples in production, 3
joint whitened and 1 per tetrode, and the largest 31, 41 and 3; joint with
the per-tetrode rules has 1 and 3. Production's
largest are TT16 unit 43 (31 samples), TT11 unit 32 (21), TT5 unit 11 (19,
a 790 uV unit), TT7 unit 19 (16) and TT4 unit 5 (11 from window 24 on, 0.4
ms). Both carries that drift rebuild their templates June's way (mean of
500 accepted spikes, full replacement), and nothing in that rebuild keeps
the waveform in place; joint matching with the blend 0.2 rebuild does not
drift, so the drift comes from the rebuild and not the matcher. A unit's
spike times then carry its drift, and a comparison between sortings at
0.33 ms finds fewer shared spikes than there are.

The first version of this table measured the offset of the trough on the
deepest wire from sample 30, which misreads positive-going units (TT16 unit
44, TT15 unit 43, TT3 unit 3 and TT2 unit 2 in the per-tetrode carry, whose
maximum is more than twice their trough's depth): their trough is a small,
variable lobe. That measure gave production 11 units, joint whitened 17 and
per tetrode 4 in window 94, and a largest offset of 54 samples, for a unit
that appears only in window 94. Seeds and re-seeds come from MountainSort5
with `detect_sign=-1`, so positive-going units are stamped on their trough.
On TT16 unit 44 in window 48 the trough is -10.6 uV at sample 30 and the
peak +29.5 uV at sample 42; across spikes the trough's position has an
interquartile range of 11 samples and the peak's of 2 (`review_units.py`).
`Recipe.recenter` shifts each rebuilt template back so that its reference
extremum stays at sample 30, and moves a positive-going unit to its peak
once. It gains nothing (`results/carries_summary_rc.csv`). On the hybrid
it moves 1 to 3 templates per window per tetrode and changes no pooled
measure by more than 0.003. Joint with June's rebuild it moves 7 to 9 per
window, and ground-truth units 0 and 20 are lost with or without it
(template cosine 0.3 to 0.6 by window 11), with small-unit precision 0.02
to 0.05 lower from window 8. Over the real 48 h (`R48_W_s45_c50_lnb_b2_rc`,
`carry_R48_pt_rc`) it moves 0 to 3 templates per window, the light-driven
responses are those of the run without it, and it adds lags that the run
without it does not have: four units moved from trough to peak (TT2 unit
2, TT15 units 43 and 2005, TT16 unit 44; 8 to 30 samples) and two TT15
units whose trough stays at sample 30 while the rest of the waveform
changes (unit 42 by 13 samples, 2009 by 12).

### Light-driven cells

`tagged_units.py` reads the per-block tables of
`subjects/ttm-003/scripts/tagged_cells_blocks.py` and lists, per tetrode,
the units that answer the pulse onsets (1 to 10 ms) in some block at 50 Hz
or more (`results/tagged_units_onset.csv`). In window 0 all three sortings
start from the same seed bank and their units agree by id (TT4 units 4 to
7 share 61 to 100 % of production's spikes within 10 samples). Onset rates
per block, Hz, for the units that carry the response:

| tetrode | production | joint, whitened | per tetrode, whitened |
|---|---|---|---|
| TT4 | 6 and 7 in blocks 1 and 2 (310, 205), then 5 (124 to 203) | 6 to block 27 (138 to 324), then re-seeded 2010 | 7 in every block (187 to 425); re-seeded 2007 from block 28 |
| TT10 | 24 in blocks 1 to 4 (103 to 134), then 26, then re-seeded 52 | 25 in every block (57 to 172) | 25 in every block (100 to 224) |
| TT11 | 28, falling from 159 to under 30 by block 37 | 28 to block 3, then 29 | 28 in every block (119 to 222) |
| TT6 | 16 at 19 and 13 in blocks 1 and 2, then 0 | 16 at 31 in block 1, then under 15 | 16 at 176 to 238 in blocks 1 to 4, then 7 to 76 |

On TT4 in window 0 the mean of the spikes after the onsets matches unit 7's
spontaneous waveform in every sorting (cosine 0.75 to 0.89; `tagged_units.py
--waveforms`); production's unit 5 has 9 such spikes there. Production's
unit 5 takes the response from block 3, while its spike times move 0.4 ms
off its trough. Conjecture: June's rebuild let unit 5's template move onto
the driven cell. The per-tetrode carry keeps each driven response on one
unit through the 48 h. Its in-pulse rates use the relaxed in-pulse rule, so
part of them can be false fits.

The template lookup of `tagged_cells_blocks.py` (cosine to the OPTO2
session's templates) picks other units in the whitened sortings (TT4 unit
4 and TT10 unit 24 per tetrode), so its `tagged` summaries understate them.

### Whitening per window

`R48_W_s45_c50_lnb_b2_rw` reran windows 48 to 59 from the per-tetrode
carry's bank after window 47 with a whitening matrix per window (30 blocks
of 10 s spread over the window, clear of the tetrode's bad periods and the
pulses). `rewhiten_real.py` compares it with the same windows of the 48 h
run (`results/rewhiten_real.csv`). On the busy tetrodes it accepts 2 to 5 %
more spikes, at a lower median cos_pl (TT8 0.684 against 0.717, TT12 0.626
against 0.647); on TT14, whose residual correlation under the first 300 s'
matrix is highest, it accepts 16 % fewer (91,755 against 109,454) and
removes fewer by the refractory rule (0.17 against 0.21). The two runs
share 55 to 95 % of each other's spikes per tetrode, and their templates
stay at cosine 0.98 or more. On the hybrid it lowers small-unit precision
from 0.87 to 0.84. Nothing here asks for re-whitening every window; TT14
may be a case for re-whitening when the residual correlation crosses a
threshold, not tested.

### Joint matching with the per-tetrode rules over the 48 h

`R48_W_joint_masked_s45_c50_lnb_b2` (run 2026-10-30) asks whether the
per-tetrode carry's real-data behaviour (driven cells on one unit, spike
times on the trough, the few-unit tetrodes' units) comes from the matcher
or from the rules and the rebuild. It is the hybrid arm of the same name
over the 96 windows: joint, masked, whitened, with every rule of the
per-tetrode carry and the blend 0.2 rebuild, from the same seed bank and
re-seeds. `export_sorting.py` assembled it as
`carry_R48_joint_rules/assembled_reseed_c12` (77 units, 133.5 M spikes).
Matching took 1.9 to 2.0 min per window, against 9 to 10 min per tetrode
for `R48_W_s45_c50_lnb_b2_rc`, which ran beside it.

`bad_periods_compare.py` against the per-tetrode carry (tag
`carry48_jr_pt`): 71 single units against 68, no spike in the tetrodes'
own bad periods against 64, no unit mostly in them, refractory
contamination median 0 in both, 7 clock-locked units in both. Joint has
fewer spikes on every tetrode but TT15 (36,976,353 against 36,881,610),
most on the few-unit tetrodes and TT9: TT1 1.60 M against 2.84 M, TT2
1.36 M against 2.48 M, TT3 0.73 M against 1.30 M, TT9 0.49 M against 1.26
M, TT14 0.20 M against 0.61 M, TT16 2.93 M against 3.56 M; on the busy
tetrodes 9 to 33 % fewer (TT10 10.9 M against 16.3 M, TT6 17.0 M against
18.6 M). Its spike times do not drift (previous section).

Onset rates (1 to 10 ms, Hz) of the units that carry the light-driven
responses (`tagged_units.py`), per tetrode against joint with the same
rules:

| tetrode | per tetrode | joint, per-tetrode rules |
|---|---|---|
| TT4 | unit 7, 187 to 425 in every block | unit 7, 117 to 352 in every block |
| TT6 | unit 16, 176 to 238 in blocks 1 to 4 | unit 16, 90 to 180 in blocks 1 to 4 |
| TT10 | unit 25, 100 to 224 in every block | unit 25, 61 to 181 in every block |
| TT11 | unit 28, 119 to 222 in every block | unit 28, 96 to 117 in blocks 1 to 4, 22 to 50 from block 37; unit 29 at 83 to 106 in blocks 41 to 44 |

`driven_spikes.py` (`results/driven_spikes.txt`) follows TT11 unit 28's
25,519 in-pulse spikes of the per-tetrode carry into the joint carry: 31 %
accepted on unit 28, 11 % fitted by unit 28 and rejected (mostly by the
snr or relative floor), 28 % accepted on another TT11 unit (unit 29 for 97
% of them), 30 % not fitted. Accepted on unit 28 fall from 63 % in windows
0 to 23 to 5 % in windows 72 to 95, while unit 29's share rises from 3 to
47 % and the unfitted from 18 to 38 %. Unit 28's and 29's templates stay
at cosine 0.90 to 0.97 with their window-0 templates in both carries. A
larger accepted spike on another tetrode within 3 ms is there for 97 % of
the unfitted spikes and 91 % of those joint accepts, so competition across
tetrodes does not separate them.

In windows 60, 72 and 84 the mean of the spikes joint leaves unfitted has
cosine 0.959 with the OPTO2 template of the TT11 cell (unit 39) and 0.887
with unit 29's OPTO2 match (unit 41), and troughs -41, -38, -12 and -13 uV
against unit 28's -75, -58, -35 and -46 uV outside the pulses. The spikes
joint gives to unit 29 have cosine 0.951 and 0.910, troughs -51, -52, -23
and -18 uV, between unit 28 (-75, -58) and unit 29 (-36, -65) on the first
two wires. Per spike in window 84 (`review_driven.py`'s analyzer), a
spike's shape is closer to unit 28's than to unit 29's for 0.98 of unit
28's spikes outside the pulses and 0.10 of unit 29's, and for 0.49 of the
spikes joint gives to unit 29 and 0.61 of those it leaves unfitted. Their
scale on unit 28's template is 0.61 and 0.42, and the former falls from
0.69 in the first tenth of the 207 ms pulse to 0.57 in its second half.
The user reviewed the analyzer: the unfitted spikes are not obviously
multi-unit and could all be one unit, they shrink across the pulse, and
the spikes joint gives to unit 29 look like a mixture, on the whole closer
to unit 28 by the CH9 waveform, while in the principal components all
three groups mostly overlap.

On the hybrid, where the driven units' half-size in-pulse spikes are
known, joint with the per-tetrode rules finds 0.50 of them against 0.61
per tetrode (means over the 12 windows; medians over the windows with
pulses 0.40 and 0.52), and the first spike of each pulse 0.85 against
0.98 (medians).

With the same rules and rebuild the two matchers are equally free of bad
periods and drift, and both fit the few-unit tetrodes' common-mode events.
What separates them on real data is that joint finds fewer spikes,
including small in-pulse spikes of the driven cells. Conjecture: these
are the driven cells' spikes shrunk by adaptation at 100 to 200 Hz, as the
hybrid's half-size spikes are; why joint misses them is not established.
