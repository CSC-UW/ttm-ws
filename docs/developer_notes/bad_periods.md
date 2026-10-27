---
title: Transient bad periods per tetrode (ttm.bad_periods, ttm-bad-periods)
updated: 2026-10-26
status: needs_verification
---

What `ttm.bad_periods` measures, how it decides, how the table is used by
the sorting, and what it does not do. The driver is `ttm-bad-periods`
(`packages/ttm/src/ttm/cli/bad_periods.py`); the sorting that uses the table
is `subjects/ttm-003/scripts/mp_baseline.py --bad-periods`; the comparison
with the unmasked sorting is `subjects/ttm-003/scripts/bad_periods_compare.py`.

## Why a recording-level table and not valid periods

SpikeInterface's `valid_unit_periods` says where one unit is well isolated
and is consumed after sorting (quality metrics, the GUI's scatter views). A
bad period here is a fact about the recording: a stretch in which one
tetrode's wires carry a signal the sorter must not learn from. It has to be
known before any unit exists, because the seed and re-seed sorts, the
per-window template re-estimation, the noise estimates and the residual and
multi-unit stages all learn from the data. The table is therefore scoped to
a tetrode, written next to the recording, and applied to the recording the
sorter reads. Per-unit valid time follows from it (a unit on a masked
tetrode has no data there), which is what `valid_seconds` returns.

## The scan

`bin_stats` reads the preprocessed recording (bandpass 300 to 6000 Hz plus
global common median reference, the sorter's input) in 60 s chunks and keeps
per 0.1 s bin: each wire's sum of squares, the cross-product of every
within-tetrode wire pair (from one batched 64 x 64 Gram matrix per bin),
each wire's largest absolute value, and its count of excursions below -5
times its noise level. Sums aggregate exactly to any coarser bin, so the
scan runs once. On the 48 h TTM-003 binary (1.3 TB, 64 channels) it takes
3 minutes with 48 workers when the binary is in the page cache.

`bin_features` makes 1 s features: rms per wire, correlation per pair and
the median correlation per tetrode. `baselines` takes a centred 30 min
running median of each wire's rms and each tetrode's median correlation,
then re-estimates both twice without the bins in which any wire of the
tetrode was more than twice its baseline (padded by 30 s). Without the
second pass a 400 s bad stretch raises its own baseline: on TT4 the clean
seconds just before 2296 s read 0.8 times the single-pass baseline.

## The decision

Per tetrode and 1 s bin:

- `loud`: some wire's rms is more than 3 times its baseline;
- `common`: the median wire-to-wire correlation is 0.25 above its baseline
  and the median wire is twice its baseline (a shared signal that adds power);
- `global`: three quarters of the tetrodes (12 of 16) are flagged in that
  bin; every tetrode is then flagged.

Flagged bins of a tetrode within 5 s of each other, or separated only by
bins in which some wire stays twice its baseline, form one interval from the
first to the last flagged bin, padded by 1 s. The table has one row per
tetrode and interval: `group`, `tetrode`, `start_s`, `end_s`, `kind`,
`channels` (the loud wires, or all four), `max_loud`, `max_corr_excess`,
`flagged_s`.

The global threshold is three quarters, not half, because on TTM-003 the
tetrodes on even-numbered channels (TT1 to TT4, TT13 to TT16) fail together
for minutes at a time, six or seven of the eight at once (2287 to 2700 s,
96900 to 97560 s), while the eight odd-channel tetrodes stay clean; at
half, those minutes would be cut from the clean tetrodes too.

## Masking

`mask_bad_periods(recording, table)` returns a lazy `MaskedPeriodsRecording`
in which each listed tetrode's channels are replaced inside its intervals,
with a 1 ms raised-cosine fade at each edge. The default fill is white
Gaussian noise at each channel's MAD noise level measured on chunks that
avoid that channel's intervals (`clean_noise_levels`); zeros would pull
down every MAD taken over a window containing the interval (a window 22 %
masked underestimates the MAD by about a third), which lowers the detection
thresholds in the clean rest of that window. White noise at the MAD level
crosses 5.5 MAD about once per half hour per wire. The noise of each
1 s block and channel comes from its own seed, so a reader gets the same
samples whatever chunk boundaries it uses, and the recording serialises to
JSON and reloads (`si.load(m.to_dict())`), which the parallel matching and
the sorter wrappers need.

A masked 60 s read of the TTM-003 binary in which eight tetrodes are masked
takes about 4 s, most of it reading the binary; chunks without an interval
cost nothing extra.

## What it does not do

- It works at the tetrode level: a flag drops all four wires, including in
  `loud` intervals where only one or two are loud (TT4 from 2373 s, CH58 and
  CH64). The table lists the loud wires so a per-wire mask can be added.
- It does not judge whether a loud or shared stretch is neural. On the
  tetrodes with few units (TT1 to TT3, TT13 to TT16) many short `common`
  intervals are bursts of high-frequency signal seen equally on all four
  wires, a few seconds long; whether they are muscle pickup or neural is not
  established.
- One-millisecond electrical pops on the good tetrodes are flagged as 3 s
  intervals (1 s bin plus padding); on TTM-003 that is under a minute per
  good tetrode over 48 h.
- Times are seconds of the recording it scanned. `ttm-bad-periods --binary`
  on a binary cropped from the session (TTM-001's `mp_long_s2000_d170000`)
  writes times relative to the crop.
