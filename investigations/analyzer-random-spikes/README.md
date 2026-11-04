---
title: Random spikes per unit for the SortingAnalyzer
updated: 2026-11-04
---

# Random spikes per unit for the SortingAnalyzer

SpikeInterface's analyzer draws up to 500 spikes per unit (`random_spikes`,
method `uniform`, its default) and computes waveforms, templates, principal
components and the PC-based quality metrics from those spikes only. The
question is whether 500 is enough for the final TTM-003-48H sorting, or
whether a larger draw changes the results.

## What was run

The final sorting (`sortings_carry/main`, 61 units, 138,680,873 spikes over
48 h) on the int16 binary it was sorted from, through
`create_group_sparse_analyzer` (each unit sparse to its tetrode) in µV:
waveforms from 1 ms before to 2 ms after, templates (average and std), noise
levels with a fixed seed, 5 PCs per channel (`by_channel_local`), and the
metrics SNR, isolation distance and L-ratio (`mahalanobis`), d′,
nearest-neighbour hit and miss rates, and silhouette. Caps of 500, 2,000,
5,000 and 20,000 spikes per unit, three seeds each. The reference is the
mean of the three 20,000 runs.

The comparison across caps uses the 56 units with more than 20,000 spikes,
which are capped at every size. The other five are reported separately
below, because a unit drawn in full keeps its own count while its
neighbours' counts grow, and that changes the PC metrics through their
definition rather than through sampling. The units alone on their tetrode
(TT3, 34 spikes; TT13, 429 spikes) have no isolation distance, L-ratio, d′
or silhouette, and their nearest-neighbour rates are set to 1 and 0.

Code: `random_spikes.py`. Results: `results/random_spikes_units.csv` (one
row per cap, seed and unit), `results/random_spikes_summary.csv`,
`results/random_spikes_small_units.csv`, `results/random_spikes_timing.csv`.

## Results

Median shift of each metric from the reference, relative to the reference
value, over the 56 units (in brackets, the share of units that moved up),
and the lowest Spearman correlation of the unit ranking with the
reference over the three seeds:

| Metric | Shift at 500 | 2,000 | 5,000 | Rank correlation at 500 | 2,000 | 5,000 |
|---|---|---|---|---|---|---|
| SNR | -0.2 % (43 %) | +0.1 % | -0.1 % | 0.995 | 0.995 | 0.998 |
| isolation distance | +3.2 % (79 %) | +0.2 % | 0.0 % | 0.970 | 0.979 | 0.990 |
| L-ratio | -8.4 % (13 %) | -1.2 % | -0.9 % | 0.970 | 0.983 | 0.993 |
| d′ | +1.0 % (70 %) | +0.2 % | 0.0 % | 0.983 | 0.988 | 0.993 |
| nearest-neighbour hit rate | -7.6 % (2 %) | -0.9 % | -0.2 % | 0.968 | 0.985 | 0.994 |
| nearest-neighbour miss rate | +14.5 % (88 %) | +0.5 % | -0.2 % | 0.956 | 0.967 | 0.988 |
| silhouette | +1.5 % (66 %) | -0.1 % | +0.3 % | 0.955 | 0.976 | 0.974 |

At 500 spikes, L-ratio, isolation distance and the nearest-neighbour rates
are shifted the same way for most units, beyond the scatter between seeds.
The largest median shift is 1.2 % at 2,000 (L-ratio) and 0.9 % at 5,000
(L-ratio). The scatter between
seeds, as a fraction of the spread between units, is at most 0.09 at 500,
0.06 at 2,000 and 0.045 at 5,000 for every metric.

The nearest-neighbour metrics subsample to 10,000 spikes per comparison,
taking every k-th spike in order (`max_spikes`, a fixed default). On a
tetrode of five units that limit is reached at 2,000 spikes per unit, so
above that the cap no longer changes how many spikes they see.

Templates (average, every sample of every channel of the tetrode):

| Cap | 1 - cosine to the reference, median (90th percentile) | Largest sample difference / template peak, median (90th percentile) |
|---|---|---|
| 500 | 0.0018 (0.0034) | 3.8 % (5.8 %) |
| 2,000 | 0.0005 (0.0008) | 1.9 % (3.0 %) |
| 5,000 | 0.0002 (0.0003) | 1.2 % (2.0 %) |
| 20,000, one seed against the mean of three | 0.00003 (0.00006) | 0.5 % (0.7 %) |

The units with fewer spikes than some cap:

| Unit | Tetrode | Spikes | L-ratio at 500 / 2,000 / 5,000 / 20,000 | Isolation distance | d′ | Nearest-neighbour hit rate |
|---|---|---|---|---|---|---|
| 2015 | TT2 | 1,813 | 1.04 / 1.25 / 3.12 / 12.31 | 20.7 / 19.6 / 15.4 / 11.5 | 1.57 / 1.59 / 1.59 / 1.56 | 0.32 / 0.31 / 0.16 / 0.05 |
| 2006 | TT9 | 3,679 | 0.11 / 0.21 / 0.99 / 4.20 | 44.5 / 36.4 / 20.2 / 13.6 | 3.79 / 3.27 / 1.67 / 1.50 | 0.83 / 0.85 / 0.84 / 0.44 |
| 2011 | TT2 | 16,931 | 0.61 / 0.70 / 0.67 / 0.77 | 25.1 / 24.1 / 24.5 / 23.2 | 2.57 / 2.51 / 2.47 / 2.39 | 0.39 / 0.40 / 0.42 / 0.40 |

Unit 2015's L-ratio rises twelvefold from 500 to 20,000 although its own
spikes are the same 1,813 from 2,000 up: its neighbours on TT2 contribute
more spikes at each cap. Its SNR stays at 4.0, since SNR depends only on
the template.

Cost per run, with 48 jobs and the binary in the page cache:

| Cap | Waveforms (float32) | Waveforms, templates, PCs | Metrics |
|---|---|---|---|
| 500 | 43 MB | 73 s | 8 s |
| 2,000 | 170 MB | 82 s | 28 s |
| 5,000 | 419 MB | 98 s | 32 s |
| 20,000 | 1.6 GB | 151 s | 46 s |

Every run also spent about 80 s finding each unit's main channel within its
tetrode (`group_main_channel_indices`, 500 spikes, seed 0), which does not
depend on the cap.

## Weighted draws against every spike

The PC metrics compare a unit's drawn spikes with the drawn spikes of the
other units on its tetrode. Under a cap every large unit contributes the
same count, so the metrics describe equal-sized samples, not the tetrode's
actual mix of spikes. L-ratio and isolation distance as defined use every
spike. A draw can still estimate those values if each drawn spike of a
neighbour stands for its unit's total over its drawn count. This part
measures that weighted estimate, and the plain capped one SpikeInterface
computes, against the values on every spike.

`weighted_vs_all.py` fits one PCA model (5 components per channel, up to
5,000 spikes per unit, seed 0, in ADC counts) and projects all 138.7 M
spikes on it with `run_for_all_spikes` (2.1 min with 48 jobs, 11 GB). The
exact values use every spike of the unit (mean and covariance) against
every spike of the other units on its tetrode. Each draw is a subset of the
same projections, at caps of 500, 2,000, 5,000 and 20,000 with five seeds,
so the comparison measures the estimator and not the PCA fit. Both
estimates take the unit's mean and covariance from its drawn spikes. The
plain one counts drawn spikes as SpikeInterface does. The weighted one
gives each drawn spike of neighbour v the weight (spikes of v) / (drawn
spikes of v), divides the L-ratio by the unit's total spike count, and
takes the isolation distance where the cumulative weight of the other
spikes, in order of distance, reaches the unit's total.

Isolation distance is undefined when a unit has at least as many spikes as
all the other units on its tetrode together; SpikeInterface then returns
the distance of the farthest other spike (`n = min(own, other)` in
`mahalanobis_metrics`), a maximum no draw can estimate. That excludes four
units from the isolation distance comparison: 0 (TT1), 22 (TT8), 23 (TT9)
and 44 (TT16), whose all-spike values are 113, 5,498, 58,193 and 102,604.

Median over units of the relative error (the median over seeds), and the
lowest Spearman correlation of the unit ranking with the all-spike values
over the five seeds; 55 units for isolation distance, 59 for L-ratio:

| Estimate | Cap | Isolation distance error | Rank correlation | L-ratio error | Rank correlation |
|---|---|---|---|---|---|
| plain | 500 | 18.6 % | 0.68 | 39.7 % | 0.71 |
| plain | 2,000 | 18.3 % | 0.71 | 42.1 % | 0.74 |
| plain | 5,000 | 18.0 % | 0.79 | 38.1 % | 0.82 |
| plain | 20,000 | 17.5 % | 0.81 | 32.3 % | 0.84 |
| weighted | 500 | 5.5 % | 0.988 | 13.1 % | 0.991 |
| weighted | 2,000 | 1.4 % | 0.995 | 3.5 % | 0.997 |
| weighted | 5,000 | 0.5 % | 0.999 | 1.7 % | 0.999 |
| weighted | 20,000 | 0.3 % | 0.9995 | 0.8 % | 0.999 |

The plain estimate does not approach the all-spike values as the cap
grows. Unit 2015 (TT2, 1,813 spikes) has an L-ratio of 100.4 on every
spike, 3.1 plain at 5,000 and 100.5 weighted. Unit 2010 (TT16, 230,475
spikes, beside a unit with 523,242) has an isolation distance of 20.8 on
every spike, 1,145 plain at 5,000 (both units capped to 5,000, so the
plain value is again the farthest other spike) and 21.0 weighted.

The weighted errors have one sign for most units: isolation distance
+5.5 % at 500, +1.4 % at 2,000, +0.5 % at 5,000; L-ratio -12.7 %, -3.4 %,
-1.4 %. The weighted sums over the neighbours' spikes are unbiased, so
these shifts come from the unit's own mean and covariance, estimated on its
drawn spikes. The largest weighted error is unit 2006 (TT9, 3,679 spikes): from
500 and 2,000 of its spikes its L-ratio is 15 and 33 against 135 on every
spike, and its isolation distance 14.1 and 12.7 against 7.9. At 5,000 it
is drawn in full and matches. At 5,000 the largest remaining relative
errors are on L-ratios near zero (unit 5: 0.0014 against 0.0020) and on
isolation distance 6 % (unit 5: 68.3 against 64.4).

## `run_for_all_spikes` ignores `return_in_uV`

SpikeInterface's `run_for_all_spikes` (`_all_pc_extractor_chunk`) reads
traces with `recording.get_traces(...)` and no `return_in_uV`, while the
PCA model it applies was fitted on waveforms extracted with the analyzer's
`return_in_uV`. On an int16 recording with an analyzer in µV it projects
ADC counts on a model fitted in µV. On TT9, its projections of the
analyzer's own random spikes differ from the analyzer's stored projections
of the same spikes by up to 54.6 (stored values up to 11.4, its values up
to 66.0); with the analyzer in counts they are identical
(`results/run_for_all_spikes_uv.csv`). `export_to_phy` uses it for
`pc_features.npy`. SpikeInterface's `main` branch had the same call on
2026-11-04, and no issue about it was found.
