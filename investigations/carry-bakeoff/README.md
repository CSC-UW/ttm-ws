---
title: Carry-forward bake-off on hybrid ground truth, TTM-003 48 h
updated: 2026-10-28
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
mixed units. A stricter gate per tetrode (cos_pl >= 0.55) moves it to
about joint's point. Most of the gain over June's carry comes from
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
(1,067): the probe-wide mask leaves out each tetrode's own bad periods, and
re-seeding there adds units that live in them. The per-tetrode sorting's
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
0.34 are removed. Conjecture: each of these units is several cells with one
footprint on its tetrode, which joint matching mostly suppressed and
per-tetrode matching admits.

### Spike times against the waveform

`stamp_drift.py` takes the mean of up to 300 of each unit's spikes on the
plain binary, with the spike time at sample 30, and finds the trough on the
deepest wire. In window 0 every sorting puts its spikes on the trough. Later
they move off it (`results/stamp_drift.csv`):

| units 5 samples or more off the trough | window 24 | window 48 | window 72 | window 94 |
|---|---|---|---|---|
| production | 2 | 4 | 5 | 11 |
| joint, whitened | 1 | 6 | 9 | 17 |
| per tetrode, whitened | 1 | 2 | 2 | 4 |

The 90th percentile of the offset in window 94 is 15 samples in production,
11 joint whitened and 3 per tetrode. Production's TT4 unit 5 stamps its
spikes 11 to 13 samples (0.4 ms) before its trough from window 24 on; TT11
unit 28 drifts by 9. Both carries that drift rebuild their templates June's
way (mean of 500 accepted spikes, full replacement), and nothing in that
rebuild keeps the trough at sample 30. A unit's spike times then carry its
drift, and a comparison between sortings at 0.33 ms finds fewer shared spikes
than there are.

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
or from the rules and the rebuild.
