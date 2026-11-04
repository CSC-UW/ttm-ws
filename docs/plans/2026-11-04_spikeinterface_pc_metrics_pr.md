---
title: spikeinterface issues and PRs for the PC metrics and run_for_all_spikes
updated: 2026-11-04
---

# spikeinterface: PC metrics on capped draws, and run_for_all_spikes in µV

Not started. Open both once the analyzer stage of `ttm-sort` is in use,
since it carries the workarounds, and check spikeinterface's current
`main` for either first. The measurements are in
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

Isolation distance is undefined when a unit has at least as many spikes as
the other units on its tetrode together, but `mahalanobis_metrics` takes
`n = min(own, other)` and returns the distance of the farthest other spike
(113 to 102,604 on four TTM-003 units). Under a cap it also happens to
both units of a tetrode with two units that both exceed the cap: each is
drawn to the cap, so the value is again the farthest other spike (one
TTM-003 unit: 20.8 on every spike, 1,145 capped at 5,000).

What a PR could do: weight the neighbours' spikes by the sorting's spike
counts in `mahalanobis_metrics`; offer the metrics on `run_for_all_spikes`
projections; return NaN where isolation distance is undefined. d′ and the
nearest-neighbour metrics also depend on the counts (one unit's d′ is 3.8
at 500 and 1.5 at 20,000) and would need their own treatment.

The analyzer stage here computes L-ratio and isolation distance on every
spike's projection (decided 2026-11-04), so it does not need the
weighting. It still needs the ADC-counts workaround of item 1, and what it
reports where isolation distance is undefined is open
(`docs/sorting_pipeline.md`, step 8).
