---
title: spikeinterface issues and PRs for the PC metrics, run_for_all_spikes and noise levels
updated: 2026-11-04
---

# spikeinterface: PC metrics on capped draws, run_for_all_spikes in µV, and noise levels

Not started. Open these once the analyzer stage of `ttm-sort` is in use,
since it carries the workarounds, and check spikeinterface's current
`main` for each first. The measurements are in
`investigations/analyzer-random-spikes/`.

## 1. run_for_all_spikes ignores return_in_uV

A bug. Open an issue together with a PR that fixes it.

`_all_pc_extractor_chunk` in `postprocessing/principal_component.py`
reads `recording.get_traces(start_frame=..., end_frame=...,
segment_index=...)` without `return_in_uV`, while the PCA model it applies
was fitted on waveforms the waveforms extension extracted with the
analyzer's `return_in_uV`. With an analyzer in µV on a recording whose gain
is not 1 (an int16 recording, or a preprocessed one that keeps its gain),
it projects ADC counts on a model fitted in µV. On TTM-003's TT9 its
projections of the analyzer's own random spikes differ from the stored
`pca_projection` of the same spikes by up to 54.6 (stored values up to
11.4); with the analyzer in counts they are identical. `export_to_phy`
writes `pc_features.npy` with it. On 2026-11-04 `main` had the same call,
and a search of the issues found none about it.

The fix passes the analyzer's `return_in_uV` through the worker context to
`get_traces`. The test: on a recording with a gain other than 1 and an
analyzer in µV, `run_for_all_spikes` at the random spikes equals the stored
`pca_projection`.

Until it is fixed, all-spike projections need an analyzer in ADC counts.

## 2. PC metrics on a capped draw

A question of method. Raise it in an issue before writing code.

The PC metrics (isolation distance, L-ratio, d′, nearest-neighbour,
silhouette) use only the `random_spikes` draw, by default up to 500 spikes
per unit. Every unit with more spikes than the cap contributes the same
count, so L-ratio and isolation distance describe equal-sized samples of
each unit, not the recording's mix of spikes, and their values change with
the cap. On the final TTM-003-48H sorting (61 units, 138.7 M spikes) the
median unit's capped L-ratio and isolation distance are 17 to 42 % off the
values on every spike at every cap from 500 to 20,000 spikes per unit
(rank correlation 0.68 to 0.84). One unit's L-ratio is 100.4 on every spike and
3.1 capped at 5,000. Weighting each neighbour's drawn spikes by its total
over its drawn count, with the unit's total in the L-ratio denominator,
brings the median error to 0.5 % (isolation distance) and 1.7 % (L-ratio)
at 5,000 per unit (rank correlation 0.999); the largest is 6 % for
isolation distance, and the largest L-ratio errors are on values near
zero.

What a PR could do: weight the neighbours' spikes by the sorting's spike
counts in `mahalanobis_metrics`, or offer the metrics on
`run_for_all_spikes` projections. d′ and the nearest-neighbour metrics also
depend on the counts (one unit's d′ is 3.8 at 500 and 1.5 at 20,000) and
would need their own treatment. The nearest-neighbour metrics subsample to
a fixed 10,000 spikes in order, so with true proportions a small unit is
represented by a handful of spikes.

## 3. Isolation distance where it is undefined

A methodological error with a small fix, separate from item 2 so it can
land on its own: an issue together with a PR.

Isolation distance is the squared Mahalanobis distance of the n-th closest
spike of other units, n being the unit's spike count, and is undefined
when a unit has at least as many spikes as the other units together.
`mahalanobis_metrics` takes `n = min(own, other)` instead and returns the
distance of the farthest other spike, a number that reads as excellent
isolation but measures one outlying spike. On every spike of four TTM-003
units it gave 113 to 102,604. The default capped draw makes it more
common: on a tetrode with two units that both exceed the cap, both units
are drawn to the same count, so both get the farthest other spike (one
TTM-003 unit: 20.8 on every spike, 1,145 capped at 5,000).

The fix returns NaN when the unit's count is not smaller than the other
units' count, with a test on two clusters of equal size. Before opening
it, check the definition in Harris et al. (2001) and Schmitzer-Torbert et
al. (2005), and whether the farthest-spike fallback was chosen on purpose
(the git history of `pca_metrics.py`).

## 4. Noise levels average chunk MADs

A question of method. Raise it in an issue before writing code.

`get_noise_levels` (method "mad") takes the MAD of each random chunk and
averages the chunks' values. The MAD is robust within a chunk, but the
average is not: one chunk inside an artifact moves a channel's noise
level, and with it every SNR on that channel. On the final TTM-003-48H
binary, five calls with the default 20 chunks of 0.5 s overestimated
channels of tetrodes with bad intervals by up to 166 %. More chunks do not
remove it: the mean over every chunk, bad ones included, is up to 19 %
above the mean over clean chunks. Neither `get_noise_levels` nor the
`noise_levels` extension can leave periods out.

What a PR could do: accept periods to leave out (per channel or per
channel group), so the chunks are drawn elsewhere. A median of chunk MADs
is robust to artifact chunks but measures something else: on TTM-003-48H
it was 0.89 to 1.04 times the mean over clean chunks, because the noise
level varies over 48 h. Measurements: `investigations/analyzer-noise-and-duplicates/`.

## Meanwhile, in the analyzer stage

The analyzer stage here computes L-ratio and isolation distance on every
spike's projection (decided 2026-11-04), so it does not need the
weighting. It still needs the ADC-counts workaround of item 1, it
reports isolation distance as missing where item 3 applies, and it
computes noise levels itself outside bad intervals for item 4
(`docs/sorting_pipeline.md`, step 8).
