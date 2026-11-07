---
title: Refractory-rejected fits in the final TTM-003-48H sorting
updated: 2026-11-07
---

# Refractory-rejected fits in the final TTM-003-48H sorting

The carry's admission (`docs/sorting_pipeline.md`, step 5, rule 3.5) keeps
the highest-snr fit of each run of a unit's fits closer than 1.5 ms to one
another and rejects the rest. Until 2026-11-07 the rejected fits were
discarded. This asks what they are: second fits of one spike, noise, or
spikes of the unit's waveform. The outcome is that the analyzer stage puts
them in their tetrode's multi-unit unit (step 8, `mua_refractory` in
`ttm.analyzer.AnalyzerParams`, decided 2026-11-07).

The sorting is the final one of TTM-003's 48 h recording (session
`2026-10-20_09-02-12`, TTM-003-48H below; `sortings_carry/main` in the
session directory, 61 units, 138.7 M spikes, 30 kHz). The scores are those
the carry saved for each fit in `carry/spikes_wNN.npz`, read in every 8th
1800 s window (windows 0, 8, ..., 88), unless stated otherwise.

Code: `refr_rejects.py`, `refr_rejects2.py`, `refr_rejects_mua.py`,
`refr_acg.py`, `refr_pairs.py`, `trace_pairs.py`, `trace_pairs_report.py`,
`summarize.py`, `mua_refr_effect.py`. Results in `results/`. The scripts
were run from this workspace's environment, each with an output folder as
its argument; the bulk outputs (the pairs table, 19 MB, and the trace
check's tables and traces, 0.5 GB) are not kept.

## Terms

- A rejected fit is one whose saved reason is 3: rule 3.5 rejected it.
  It passed rules 3.1 to 3.4, or, for a sustained responder near a laser
  pulse, the pulse rule's relaxed criteria (step 5, item 4), which
  TTM-003-48H used.
- A chain is a run of one unit's accepted and rejected fits, each closer
  than 1.5 ms (45 samples) to the previous one. Its survivor is the fit
  rule 3.5 kept, and a rejected fit's lag is its distance to the survivor.
- A unit's isolated spikes are its accepted fits in a chain of one.
- A percentile places a fit among the isolated spikes of its unit in the
  same window, on a scale of 0 to 1, unless stated otherwise.
- Amplitude is the scale the matcher fitted to the template, 1 being the
  template's size.

## The saved scores

In the 12 windows rule 3.5 rejected 823,315 fits (`refr_rejects.py`,
`summarize.py`; `results/lag_hist.csv`, whose shares are of all of them).

- The median lag is 0.73 ms. 8.4 % of lags are under 0.3 ms, 0.12 % under
  0.1 ms, and 4.7 % over 1.5 ms (in chains of three or more). Per sample,
  the share rises from 0.09 % at 0.07 ms to 0.57 % at 0.10 ms and 2.8 % at
  0.30 ms, is 2.8 to 3.6 % from 0.33 to 0.63 ms, and falls to 1.5 to 1.6 %
  between the bumps at 1.13 to 1.43 ms. On top of that it has bumps every 6 to 7 samples
  (0.6, 0.83, 1.03 to 1.07, 1.27 and 1.47 ms), which come from unit 42
  (below).
- Rejected fits sit at a median percentile of 0.32 in amplitude, 0.32 in
  snr and 0.21 in cosine, their survivors at 0.71, 0.73 and 0.65
  (`refr_rejects2.py`, `results/refr_rejects_controls.csv`; medians over
  unit-windows). Rule 3.5 keeps the higher snr of two fits; for the lower
  and higher of two spikes drawn at random from one unit these medians
  would be 0.29 and 0.71, close to what the rejected fits and survivors
  show. Isolated spikes with another unit's accepted spike within 1.5 ms
  sit at 0.52, 0.54 and 0.43: an overlapping spike lowers the cosine
  percentile by about 0.07, and rejected fits, which always overlap their
  survivor, sit 0.08 below 0.29.
- 14.0 % of rejected fits have another unit's accepted spike within
  0.5 ms, against 10.8 % of isolated spikes (`results/refr_rejects_units.csv`),
  so they coincide with other units only slightly more often.
- For a source of events firing independently of the unit to give as many
  rejected fits by chance coincidence within 1.5 ms, it would need, for
  the median unit-window, 0.84 times the unit's accepted count. 37 % of
  rejected fits are in unit-windows where it would need more than the
  unit's own count.
- 2.6 % of the rejected fits had an event of the tetrode's multi-unit unit
  within 10 samples (0.33 ms): 0.1 % at lags under 0.5 ms, 2.2 % at 0.5 to
  1.0 ms and 5.1 % at 1.0 to 1.5 ms (`refr_rejects_mua.py`,
  `results/refr_rejects_mua.txt`). The fits themselves were discarded; for
  those 2.6 % a detected event in the multi-unit unit stood in for them.
  The multi-unit step leaves out detected events within 0.5 ms of a sorted
  spike.

## Autocorrelograms

Autocorrelograms per unit and window, each divided by its own level at
lags of 50 to 100 ms, over the 576 unit-windows with 2,000 or more
accepted spikes (55 units; `refr_acg.py`, `refr_pairs.py`;
`results/refr_acg.csv`, `results/refr_pairs_acg.csv`). Values are medians
over unit-windows.

- With the rejected fits put back (every fit that passed rules 3.1 to
  3.4), the value is 0.06 at 0 to 0.2 ms, then 1.15, 1.70, 1.36 and 1.22
  in 0.2 ms bins to 1.0 ms, and 0.90 and 0.79 at 1.0 to 1.4 ms.
- The accepted spikes alone do not fall below their 50 to 100 ms level
  after the 1.5 ms that rule 3.5 imposes: 0.97 at 1.5 to 2.0 ms
  (quartiles 0.70 to 1.35), 1.57 at 2.0 to 2.6 ms, 1.20 at 2.6 to 3.4 ms,
  1.01 at 3.4 to 4.2 ms and 1.22 to 1.40 from 4.2 to 8 ms. The 1.5 to
  2.0 ms bin is the lowest of these.
- Accepted pairs 2.0 to 2.6 ms apart look like ordinary spikes of the
  unit. Placed among the unit's accepted spikes with no accepted spike of
  the unit within 10 ms, both spikes sit at median percentiles of 0.53 to
  0.55 in amplitude and snr and 0.47 to 0.48 in cosine
  (`results/refr_pairs_features.csv`).

## Trace check

`trace_pairs.py` asks of the traces whether a rejected fit and its
survivor are two events of the template's shape or one spike fitted
twice. It uses windows 8, 20, 32, 44, 56, 68, 80 and 92, thirty 10 s
blocks of each, leaving out a tetrode's blocks that overlap its bad
intervals, on the whitened traces with every other unit's fits
subtracted, whatever their reason. It takes three kinds of event, each
with no other fit of the unit at snr 3 or more in its window:

- rejected: a rejected fit and its survivor in a chain of two, lag 0.3 to
  1.5 ms (20,200 events, 53 units);
- far: two accepted spikes of the unit 2.0 to 2.6 ms apart, which rule
  3.5 accepts as two spikes (13,570);
- lone: an accepted spike with no fit of any unit on the tetrode within
  3 ms (156,196).

On a 270-sample window around each event it fits, by least squares, one
copy of the unit's template that may shift, widen or narrow (the template
and its first and second time derivatives), the same with the template of
any unit of the tetrode, and two copies at the fitted frames, each allowed
to move by up to 2 samples. Each event also gets a twin: the matcher's own
fits of the event (templates at their frames and amplitudes) added to a
stretch of the same traces without any fit of the unit at snr 3 or more.
A twin holds two events of the template's shape by construction, so a
measure is read against its twin.

`trace_pairs_report.py` (`results/trace_pairs_summary.csv`,
`results/trace_pairs_per_unit.csv`):

| | rejected 0.3-0.5 ms | 0.5-0.7 | 0.7-1.0 | 1.0-1.5 | far 2.0-2.6 |
|---|---|---|---|---|---|
| events | 4,960 | 4,011 | 4,973 | 6,256 | 13,570 |
| gain of two copies over one flexible copy, data | 121 | 147 | 149 | 160 | 183 |
| the same, twin | 118 | 151 | 155 | 165 | 193 |
| data over twin, median per event | 1.021 | 0.975 | 0.975 | 0.979 | 0.969 |
| one flexible copy of any unit fits better than two copies, data | 0.6 % | 0.3 % | 0.3 % | 0.1 % | 0.0 % |
| the same, twin | 0.2 % | 0.1 % | 0.1 % | 0.0 % | 0.0 % |
| misfit of two copies, data minus twin | 0.102 | 0.127 | 0.137 | 0.137 | 0.109 |

The gain is the energy (sum of squares, in whitened noise units, median
per event) that two copies explain beyond one flexible copy. The misfit is
the residual energy over the fitted spikes' span beyond what the rest of
the window predicts, over the energy of the fitted spikes; for lone spikes
the data exceed their twins by 0.080. The background is 1.01 to 1.04 per
sample around rejected events and 0.91 to 0.95 in their twins' stretches,
so each comparison is made within an event.

- In the data a second copy explains as much beyond one flexible copy as
  it does in the twins. What one spike fitted twice would give was not
  simulated. The misfit of the rejected pairs at 0.5 to 1.5 ms (0.13 to
  0.14 beyond their twins) is larger than that of the far pairs (0.11),
  so the matcher's two fits explain these events a little less well.
- Per unit, over the 39 units with 30 or more rejected events at 0.3 to
  0.7 ms, the median data-over-twin ratio is 1.00 (10 to 90 %: 0.95 to
  1.04, range 0.88 to 1.16). One flexible copy wins in more than
  (2 × the twins' share) + 2 percentage points in one unit, unit 26 of
  TT10 (3.2 % of 188 events against 0.0 %).
- `results/trace_pairs.png` shows the main wire's average in µV for six
  units at their most common lag, read off the figure. In five of them
  the data show two troughs, the first 60 to 88 µV deep and the second 55
  to 70 µV, against 70 to 92 µV for the unit's lone spike, and the twins
  follow the data closely. The sixth, unit 40 of TT15, has a narrow
  positive peak of about 35 µV as its lone spike, and its average has
  high-frequency ringing before the event.
- The lags of one unit are spread out (`results/lag_mode_per_unit.csv`).
  Over the 50 units with 300 or more chains of two at 0.3 to 1.5 ms, a
  median 12.9 % of a unit's lags lie within one sample of its most common
  lag (10 to 90 %: 10.6 to 18.0 %), against 8.1 % for lags spread evenly.
- Unit 42 of TT15 is the exception: 28.5 % of its 57,395 chains of two lie
  within one sample of 1.07 ms, and its lags come in a comb, with peaks at
  0.10, 0.40, 0.60 to 0.63, 0.83, 1.03 to 1.07, 1.27 and 1.47 ms, from
  0.40 ms on about 6.5 samples (0.22 ms) apart. TT15's units follow a signal synchronous
  with the acquisition clock (`docs/reports/2026-10-26_ttm003_bad_periods.md`).
  Unit 42's rejected fits may be its template matched one cycle apart on
  an oscillation (conjecture); the trace check cannot tell that from two
  events, since two shifted copies fit an oscillation too (its
  data-over-twin ratio is 0.999).
- The earlier fit is the larger in a median 53.8 % of a unit's chains of
  two (range 36.7 to 62.8 %).

## Lags under 0.3 ms

The trace check leaves out lags under 0.3 ms, where one flexible copy can
stand in for two. Chains of two at those lags look like those above it
(`results/pairs_by_lag.csv`): at 0.2 to 0.3 ms the rejected fit's median
amplitude is 0.85 and its cosine 0.55, against 0.87 and 0.55 at 0.3 to
0.4 ms. At 0.1 to 0.2 ms they are 0.78 and 0.55. 816 of the 668,403
chains of two (0.12 %) are under 0.1 ms. The multi-unit units take every
rejected fit, whatever its lag.

## Rejected fits in the multi-unit units

`mua_refr_effect.py` repeats the analyzer stage's multi-unit selection on
all 96 windows of TTM-003-48H and of TTM-003-48H-2 (session
`2026-11-03_09-01-38`, sorted 2026-11-06 with the same pipeline, without
the pulse rule), with and without the rejected fits
(`results/multi_unit_counts_*.csv`). A detected multi-unit event within
0.5 ms of a rejected fit of its tetrode gives way to it.

| | TTM-003-48H | TTM-003-48H-2 |
|---|---|---|
| sorted spikes | 138.68 M | 131.57 M |
| rejected fits | 6.86 M (4.9 %) | 5.62 M (4.3 %) |
| detected multi-unit events | 1.17 M | 0.61 M |
| detected events giving way | 172,456 | 74,761 |
| multi-unit events with the rejected fits | 7.85 M | 6.16 M |

On the other tetrodes with 1,000 or more sorted spikes the rejected fits
are 2.1 to 8.8 % of the sorted spikes, except TT2 of TTM-003-48H-2 at
0.3 %. On TTM-003-48H they are 15.6 % on TT1, 12.4 % on TT2 and 17.9 % on
TT16; on TTM-003-48H-2, 14.6 % on TT8 and 33.4 % on TT16. TT15 has 1.19 M
on TTM-003-48H and 1.49 M on TTM-003-48H-2. TT3 and TT13 of TTM-003-48H,
whose units are too small for the multi-unit step's templates, get
multi-unit units of 7 and 5 rejected fits.

## Reading

Conjecture, not confirmed. Outside TT15, the rejected pairs at 0.3 to
1.5 ms behave on the traces like two spike-sized events with the unit's
waveform. One neuron should not fire twice within that time, so a unit
with many of them may hold two or more neurons whose waveforms the
tetrode does not tell apart, firing within about 1 ms of each other more
often than chance. That the accepted autocorrelograms do not fall below
their 50 to 100 ms level after 1.5 ms fits that reading. If so, the
waveform cannot say which neuron fired either spike, which is why the
rejected fits go to the multi-unit unit rather than back to the unit.
