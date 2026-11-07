---
title: spikeinterface issues and PRs for the PC metrics, run_for_all_spikes, noise levels, sd_ratio, zarr gathering and periods
updated: 2026-11-07
---

# spikeinterface: PC metrics on capped draws, run_for_all_spikes in µV, noise levels, sd_ratio, zarr gathering, and periods

Item 1 is open as PR #4925; the others are not started. Open these once the analyzer stage of `ttm-sort` is in use,
since it carries the workarounds, and check spikeinterface's current
`main` for each first. The measurements are in
`investigations/analyzer-random-spikes/`.

## 1. run_for_all_spikes ignores return_in_uV

A bug. Opened on 2026-11-07 as
[PR #4925](https://github.com/SpikeInterface/spikeinterface/pull/4925)
(branch `fix-pc-all-spikes-uv-segments` on the fork), with the context in
the PR body and no separate issue. The PR also fixes a second bug in the
same functions, below.

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

The second bug: `_init_work_all_pc_extractor` keeps only segment 0's
spikes, but the chunk executor runs over every segment. Chunks of a later
segment project that segment's traces at segment 0's spike times and
overwrite segment 0's rows, and the later segments' spikes stay zero.
`export_to_phy` asserts one segment, so only direct calls are affected,
and not TTM-003-48H, whose binary has one segment.

The fix passes the analyzer's `return_in_uV` through the worker context to
`get_traces`, and keeps the spike vector of every segment in the worker
with each segment's slice of it, so a chunk searches only its own
segment's spikes. The test extends `test_compute_for_all_spikes`: on a
copy of the two-segment test recording with gain 0.195 and offset -100,
`run_for_all_spikes` at the random spikes equals the stored
`pca_projection`. On `main` it fails (errors up to 179 on values up to 6.7,
and the second segment's spikes all zero); with the fix the difference is
0.

Until #4925 is in a release, all-spike projections need an analyzer in
ADC counts.

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

## 5. sd_ratio draws its own noise levels, unseeded

A question of method with a reproducibility bug in it. Raise it in an
issue; the seed alone could be a PR.

`compute_sd_ratio` (`metrics/quality/misc_metrics.py`) divides each unit's
amplitude spread by a noise level it gets from
`get_noise_levels(recording, method="std")` itself, not from the
analyzer's `noise_levels` extension. That call draws 20 random chunks of
0.5 s with `seed=None` (`get_random_sample_slices`) and caches the result
as a recording property (`noise_level_std_scaled`), so the value holds
within one process and changes between processes. On a 10 min slice of
TTM-003-48H, three builds in one process gave identical `sd_ratio`, and
builds in three processes differed for 58 of 74 units, by up to 0.26.
On the whole of TTM-003-48H, five seeds gave standard deviations up to
2.8 times apart on one channel (no unit's main channel) and up to 1.10
times apart on the units' main channels. The draw ignores `periods`, so
chunks inside bad intervals count. A noise level set on the
`noise_levels` extension (here the stage's, outside bad intervals) has no
effect on `sd_ratio`.

What a PR could do: pass the analyzer's seed and the metric's `periods` to
the draw, or take the noise as an argument (the analyzer stage here sets
it through the cache property instead; see the last section). Taking the
analyzer's `noise_levels` extension, as SNR does, would change the
metric. That extension is MAD-based (MAD / 0.6745, which equals the std
for Gaussian noise). `sd_ratio` is near 1 when a unit's amplitudes vary
only as much as additive noise makes them vary, and that noise includes
the other units' spikes that overlap the unit's own. A chunk's std counts
those spikes, and `compute_sd_ratio` subtracts from it only the estimated
variance of the unit's own spikes; the MAD is barely moved by sparse
spikes and leaves most of them out. On TTM-003-48H the std over clean
chunks is 1.01 to 1.32 times the MAD noise, and `sd_ratio` over the MAD
noise, without that subtraction, would be 1.02 to 1.32 times its value
over the std.

## 6. Zarr gathering appends one chunk at a time

A performance bug. Open an issue together with a PR.

For a zarr analyzer, the waveforms extension
(`extract_waveforms_to_single_buffer`, mode "zarr") and the node-pipeline
extensions such as spike amplitudes gather through `GatherToZarr`
(`core/node_pipeline.py`), whose `__call__` appends each worker's result
to the zarr array as it arrives. With the default 1 s chunks that is one
append per second of recording, and each append rewrites the zarr chunk
being filled (about 10 MB before compression), in the main process,
while the workers wait. On TTM-003-48H the waveforms of 354,478 random
spikes (172,733 appends) took 2.1 min per 7,281 rows; stopped after 28
of 49 zarr chunks, it was on course for about 100 min, against 89 s for
the same extension on an in-memory analyzer.
In a benchmark one append cost 36 to 85 ms for amplitude rows and 34 to
283 ms for waveform rows, rising as the open zarr chunk filled.

The fix buffers rows in `GatherToZarr` and appends whole zarr chunks (the
last one at `finalize_buffers`). The test: the same arrays as now, and an
append count equal to the number of zarr chunks.

## 7. Each quality metric selects the periods again

A performance issue. Open an issue, with a PR if the fix stays small.

With `periods`, ten of the fourteen metric functions the analyzer stage
uses call `sorting.select_periods(periods)` themselves (for example the
first lines of `compute_sd_ratio`). On the TTM-003-48H sorting with its
multi-unit units (139.9 M spikes, 3,691 periods, each unit's intervals
outside its tetrode's bad intervals) one selection takes 88 s, so the
repeated selections account for about 15 of the 47 min the 14 metrics
took. Where the rest goes was not measured: num_spikes alone took 208 s
and firing_rate 187 s.

What a PR could do: select the periods once in the `quality_metrics`
extension and pass the selected sorting to every metric.

## Meanwhile, in the analyzer stage

The analyzer stage here computes L-ratio and isolation distance on every
spike's projection (decided 2026-11-04), so it does not need the
weighting. It needs the ADC-counts workaround of item 1 until #4925 is
in a release, it
reports isolation distance as missing where item 3 applies, and it
computes noise levels itself outside bad intervals for item 4
(`docs/sorting_pipeline.md`, step 8). For item 5 it measures the mean
std of the same chunks and makes `sd_ratio` divide by it by setting the
recording property in which `get_noise_levels` caches its result
(`ttm.analyzer.set_sd_ratio_noise`). For item 6 it computes every
extension on an in-memory analyzer and writes the zarr store once with
`save_as`. It has no workaround for item 7.
