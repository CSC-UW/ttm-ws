---
title: TTM-001 48 h, wobble in place of circus-omp for the full-size stages over a 12-window span
updated: 2026-10-23
---

# Wobble for the full-size stages, TTM-001 span test

Track G, items 1 to 3, of `docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`.
Scripts and small results: `investigations/wobble-fullsize-ttm001/`. Bulk
products: `/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52/sortings_seed42_pcafix/track_eval/wobble_span/`.
Figures: `/nvme/neuropixels/tetrode_data/TTM-001/TTM-NOD/figures/wobble_span/`.
Numbers are measurements. Sentences marked "Conjecture" are the
analyst's reading and have not been confirmed by the user.

## Verdict

On this span wobble, in the configuration June settled on, loses to
circus-omp. It puts more of the large events into single units (95.4 % of
detected events at 10 MAD or more against 80.8 %), but assignment purity
falls (tight window 0.855 against 0.884, full window 0.943 against 0.969;
a re-run of June gives 0.886 and 0.970), within-cluster spread widens on
the large units (10th-percentile cosine 0.615 against 0.715), and 14 % of
June's spikes, which have their June unit's waveform at full size, go to a
neighbour (2.9 % between June and a re-run of June). Two mechanisms
are measured: inside the carry-forward the per-window re-estimation
shrinks wobble's templates (median trough 103 to 72 µV in six windows) and
the threshold, tied to the bank's median energy, shrinks with them; and
the cosine gate, scored on the raw 3 ms snippet, deletes spikes that have
a neighbour within that window, which hollows every same-tetrode
correlogram and leads the CCG-guarded merge to join cells that fire
independently.

A diagnostic arm that is not June's loop, wobble with templates
re-estimated only from spikes at amplitude scale 0.8 or more, holds its
bank and gives higher tight purity than June (0.909) with fewer units (59
against 74), and loses 48 % of June's spikes on the eleven templates its
energy threshold penalises. Whether that variant wins is undecided; what
would decide it is listed at the end. The full 95-window pass was not run.

## What was run

The span is windows 55 to 66 of the June carry-forward (1800 s each, 28.1
to 34.1 h of the recording, 6.0 h). June re-seeded after matching windows
12, 24, ..., 84 (`reseed_sorts_c12/reseed_w*` and the births in
`reseed_cadence_sweep.npz`); the re-seed after window 60 added ten units,
the most of any boundary, so the span has six windows on each side of it.
Window 55 and most of window 56 fall in the manifest's novel-object
period (23.95 to 28.95 h after the session start); the whole span is in
the light period.

The bank entering window 55 is rebuilt from June's `assembled_reseed_c12`:
each unit's template from the last window up to 54 in which it had at
least 100 spikes, computed as June's loop computed it (80 units). The
re-seed units are June's rule (snr 5 or more, 100 spikes or more,
shift-cosine under 0.8 to every live same-tetrode template) applied to
June's MountainSort5 re-sort of window 60, against the bank after window
60's re-estimation rebuilt the same way. The rule returns June's ten units,
ids 80 to 89, with June's at-birth cosines to four decimals. Every arm
starts from this bank and receives these ten units after window 60.

Stage-1 arms (`_span.ARMS`), each with 16 jobs:

- `circus`: circus-omp, amplitudes 0.8 and up, June's matcher, re-run as a
  check. On window 55 it reproduces 1,006,950 of June's 1,007,150 spikes
  at the same sample and unit; the exact agreement falls to 79 % by window
  66 as small differences pass through re-estimation.
- `w35`: wobble with the objective threshold at 0.35 times the window bank's
  median template energy, then the per-spike cosine to the assigned
  template at 0.60 or more (`run_matching(shape_gate_r=0.60)`), the rest of
  `wobble_method_kwargs` at June's values (approx_rank 4, jitter 8,
  amplitude variance 1). This is June's documented operating point:
  `MATCHING_PURSUIT_FINDINGS.md` gives the threshold 0.35 as the measured
  saturation of 10 MAD coverage (script 81) and the gate 0.60 as the user's
  choice over the measured knee of 0.57 (script 79).
- `w25`: the same with the threshold at 0.25, the lowest admit June
  mapped (script 81: the same pooled coverage as 0.35). June did not leave
  the gate open, but the admit is the knob that decides how much wobble's
  energy threshold asks of a small template (a template with energy `q`
  times the median must reach amplitude scale `(1 + f / q) / 2`), which is
  the weakness the plan names, so the two lowest mapped admits were run.
- `w35a`: `w35` with templates re-estimated only from spikes at amplitude
  scale 0.8 or more of the current template (the population circus-omp
  keeps), added once `w35` showed its bank shrinking. Not June's loop; a
  diagnostic.

Downstream, unchanged: `97_production_pipeline.py` imported with its output
directory redirected, run on the span's windows only, with June's detected
events (`spike_coverage.npz`) cut to the span. CCG-guarded merge (cosine
0.90 proposes, CCG duplicate or cosine 0.95 accepts), residual single-unit
pass (unclaimed events of 10 MAD or more offered to the best same-tetrode
template at cosine 0.8, kept if they do not break its 1.5 ms refractory
period), per-tetrode multi-unit bucket (7 MAD or more, cosine 0.55 to 0.8)
and the three-axis scoreboard. June's residual pass contains no matcher,
so "wobble in the residual pass" is that pass run on wobble's stage-1
output; no wobble-on-the-residual variant was built. The June arm of the
comparison is June's own `assembled_reseed_c12` cut to the span and taken
through the same downstream; the circus re-run taken through it gives the
noise floor.

Amplitude classes: a stage-1 template is large when its trough on its
deepest wire, at the span start, is at least 10 excised MAD (noise of the
first 5 min of window 55 with spikes left out, `ttm.capture.noise_mad_excised`),
small otherwise: 37 large, 53 small. The cut is the residual pass's event
floor. Eleven templates have energy low enough that wobble's threshold at
0.35 asks more than circus-omp's 0.8 of them (`bar35 > 0.8`); they are
reported separately. A merged unit takes the class of its member with the
most spikes.

## Stage 1: what each matcher admits

Spikes per window (`results/stage1_windows.csv`) and the median trough of
the bank entering each window:

| window | June | circus re-run | w35 | w25 | w35a | bank trough, circus / w35 / w25 / w35a (µV) |
|---|---|---|---|---|---|---|
| 55 | 1,007,150 | 1,007,016 | 1,525,906 | 1,768,470 | 1,525,906 | 103 / 103 / 103 / 103 |
| 56 | 1,019,237 | 1,017,346 | 1,848,414 | 2,261,864 | 1,458,810 | 109 / 100 / 97 / 110 |
| 58 | 576,472 | 567,487 | 1,329,591 | 1,682,988 | 798,314 | 97 / 80 / 75 / 101 |
| 60 | 524,872 | 524,860 | 1,450,535 | 1,868,047 | 708,141 | 94 / 72 / 66 / 106 |
| 62 | 652,692 | 648,765 | 1,589,916 | 2,060,431 | 858,534 | 108 / 75 / 71 / 122 |
| 64 | 782,238 | 776,798 | 1,821,746 | 2,329,137 | 1,061,145 | 107 / 76 / 73 / 105 |
| 66 | 766,903 | 753,147 | 1,807,558 | 2,312,599 | 971,054 | 104 / 73 / 67 / 105 |
| span | 8,672,048 | 8,610,896 | 19,407,226 | 24,434,936 | 11,958,696 | |

Unit by unit against June (same unit within 0.33 ms; `04_stage1_compare.py`,
`results/stage1_units_<arm>.csv`, `results/stage1_summary.json`). June
dropped and June reassigned are fractions of June's spikes; new and from a
neighbour are fractions of the arm's:

| arm | class | spikes / June | June dropped | June reassigned | new | from a neighbour |
|---|---|---|---|---|---|---|
| circus re-run | large | 0.995 | 3.9 % | 3.3 % | 3.3 % | 3.5 % |
| circus re-run | small | 0.991 | 7.6 % | 2.5 % | 6.9 % | 2.4 % |
| w35 | large | 1.85 | 3.2 % | 12.2 % | 42.7 % | 11.6 % |
| w35 | small | 2.56 | 8.2 % | 15.9 % | 64.9 % | 5.4 % |
| w25 | large | 2.09 | 3.2 % | 12.7 % | 48.1 % | 11.6 % |
| w25 | small | 3.42 | 7.5 % | 18.3 % | 73.1 % | 5.2 % |
| w35a | large | 1.28 | 5.9 % | 16.4 % | 26.0 % | 13.4 % |
| w35a | small | 1.46 | 19.9 % | 11.1 % | 44.8 % | 7.9 % |

On the eleven templates with `bar35 > 0.8`, w35a drops 47.6 % of June's
spikes (circus re-run 11.9 %, w35 17.1 %). Conjecture: the energy
threshold's weakness, which the plan expected only for small-spike
capture, also reaches the low-energy templates of a full-size stage. The median trough of wobble's new
spikes is 14.3 MAD on large units and 7.0 on small for w35 (shared spikes
18.4 and 8.3). Large units inflate as well as small: unit 21 (22 MAD
template) goes from 1,557 spikes to 124,249 in w35, unit 25 from 191 to
25,052, unit 37 from 7,219 to 113,865, unit 20 from 22,583 to 264,547.

The galleries (`stage1_gallery_<arm>_<class>_u<unit>.png`, four large and
four small units per arm, those where the two matchers differ most; rows:
spikes both keep, wobble's new spikes, wobble's spikes that June gave a
neighbour, June's spikes wobble dropped, June's spikes wobble gave a
neighbour; 40 snippets from window 61 each over June's template, with the
trough histogram) show:

- Unit 20 (TT6, w35): 248,100 new spikes at a median 7.2 MAD whose mean,
  in window 61, is 0.30 of the template's trough with a different wire
  footprint (cosine of per-wire trough depths 0.90; the template is -160,
  -55, -40 and -125 µV, the mean about -40 on every wire); and 11,372 of
  June's spikes, template-shaped, assigned to a neighbour.
- Unit 31 (TT7, w35): 169,456 new spikes, at 0.53 of the template's trough
  with the template's footprint (0.998) and waveform cosine 0.96; 12,050 of
  June's template-shaped spikes assigned to a neighbour.
- Unit 48 (TT11, w35): 1,539,590 new spikes, at 0.39 of the trough,
  footprint cosine 0.96; 100,001 of June's spikes assigned to a neighbour.
- Unit 57 (TT13, w35a): 369,917 new spikes at a median 9.4 MAD whose mean
  has the template's footprint (0.9997) and 0.93 of its trough, but
  waveform cosine 0.86.

Over every unit, not only the gallery's (`09_new_spike_shape.py`, window
61, `results/new_spike_shape.json`; mean of up to 2000 spikes per group
against the mean of June's spikes of the unit; medians over units):

| arm | group | units | waveform cosine | footprint cosine | trough ratio |
|---|---|---|---|---|---|
| circus re-run | new, large / small | 20 / 45 | 0.988 / 0.986 | 1.000 / 1.000 | 0.90 / 0.86 |
| w35 | new, large / small | 27 / 51 | 0.981 / 0.985 | 0.999 / 0.999 | 0.77 / 0.74 |
| w35a | new, large / small | 23 / 42 | 0.968 / 0.978 | 0.999 / 0.998 | 0.81 / 0.85 |
| w35 | June's, given to a neighbour, large / small | 24 / 47 | 0.980 / 0.979 | 0.999 / 0.999 | 1.01 / 1.10 |

On average, then, wobble's new spikes are smaller copies of the template;
the gate admits only events at cosine 0.60 or more to it, so their mean is
pulled toward it by selection, and the test cannot say whose spikes they
are. The spikes June gave the unit and wobble gives a neighbour have the
unit's waveform at its full size.

The overview (`stage1_overview.png`) has the count ratio per unit against
template amplitude, the bank's median trough per window and each wobble
bank's troughs relative to circus's at window 66.

## The full stages: single units, multi-unit bucket, scoring

Single units are the merge plus the residual pass (`assembled_prod_sua`),
the deliverable adds the per-tetrode multi-unit bucket (`assembled_prod`).
Coverage is the share of June's detected events (locally exclusive troughs
at 5.5 MAD, `spike_coverage.npz`) in the span that some unit on the same
tetrode claims within 0.5 ms. Purity is June's assignment purity: the share
of a unit's spikes whose best same-tetrode template (window-local
templates) is its own, over the full template window and over the tight
trough window, spike-weighted. CCG pairs are June's arbiter over all
same-tetrode pairs (central 1.5 ms against the 5 to 25 ms flank; under 0.3
duplicate, over 0.7 distinct). From `results/scoreboard.csv`:

| arm | single units | 10 MAD covered, single units | with bucket | 7 to 9 MAD, with bucket | purity, full | purity, tight | CCG duplicate pairs | CCG distinct pairs |
|---|---|---|---|---|---|---|---|---|
| June | 74 | 80.8 % | 95.6 % | 75.1 % | 0.969 | 0.884 | 30 | 69 of 163 |
| circus re-run | 76 | 80.8 % | 95.6 % | 75.1 % | 0.970 | 0.886 | 32 | 77 of 178 |
| w35 | 68 | 95.4 % | 96.8 % | 84.2 % | 0.943 | 0.855 | 55 | 22 of 136 |
| w25 | 62 | 95.8 % | 96.8 % | 85.3 % | 0.938 | 0.854 | 40 | 14 of 107 |
| w35a | 59 | 91.3 % | 95.5 % | 76.2 % | 0.962 | 0.909 | 45 | 22 of 110 |

Along the way (`results/summary.json`): the merge takes June's 90 stage-1
units to 74 with 17 merges, the circus re-run's to 76 with 15, w35's to 68
with 22, w25's to 62 with 29, w35a's to 59 with 38. June leaves 1.36 M
events of 10 MAD or more unclaimed and the residual pass recovers 20.0 % of
them; w35 leaves 265 k and recovers 0.6 %, w35a leaves 517 k and recovers
4.4 %. The multi-unit stage sees 4.55 M events of 7 MAD or more in June
(65.6 % to the bucket, 23.5 % noise), 1.66 M in w35 (41.3 %, 57.4 %) and
3.00 M in w35a (51.1 %, 45.2 %).

Per unit, by amplitude class (`results/units_<arm>.csv` from the purity
scoring, `results/units_by_class.csv`). Refractory contamination is
Llobet's at 1 ms over the span's 6 h; the 10th-percentile cosine and the
residual energy (the snippet minus the scaled window template, in noise
units, 1 for template plus noise) measure within-cluster waveform spread;
verdicts count each unit's same-tetrode neighbours (`ttm.capture.ccg_verdict`,
June's thresholds):

| class | arm | units | spikes | contamination, median | units >= 0.1 | purity tight, spike-weighted | units under 0.8 | cosine p10, median | residual energy, median | duplicate / ambiguous / distinct verdicts |
|---|---|---|---|---|---|---|---|---|---|---|
| large | June | 27 | 4.26 M | 0.072 | 13 | 0.939 | 3 | 0.715 | 2.65 | 4 / 41 / 71 |
| large | w35 | 21 | 6.28 M | 0.001 | 0 | 0.886 | 5 | 0.615 | 3.02 | 28 / 36 / 17 |
| large | w25 | 19 | 7.12 M | 0.001 | 0 | 0.877 | 4 | 0.610 | 2.97 | 16 / 33 / 13 |
| large | w35a | 19 | 4.39 M | 0.000 | 0 | 0.938 | 1 | 0.649 | 3.49 | 22 / 27 / 21 |
| small | June | 47 | 4.68 M | 0.201 | 34 | 0.835 | 14 | 0.642 | 2.47 | 56 / 88 / 66 |
| small | w35 | 47 | 13.07 M | 0.004 | 1 | 0.841 | 14 | 0.585 | 2.44 | 90 / 78 / 23 |
| small | w25 | 43 | 17.21 M | 0.003 | 0 | 0.844 | 14 | 0.576 | 2.19 | 70 / 71 / 11 |
| small | w35a | 40 | 7.56 M | 0.003 | 1 | 0.892 | 9 | 0.624 | 2.67 | 70 / 59 / 19 |

The units differ between arms (the merges differ), so a class row compares
populations, not the same cells. `results/pairs_june_<arm>.csv` pairs units
by their stage-1 members: 58 of June's 74 single units have a w35 unit
with the same members, and the median paired unit holds 82 % of the June
unit's spikes while only 36 % of its own spikes are in the June unit (w35a:
46 same, 72 % and 51 %). `results/match_june_<arm>.csv` pairs by template
cosine and spike agreement, as the TTM-003 pages pair them (median
agreement 0.81 for w35, 0.80 for w35a, 0.86 for w25).

The purity pages (`purity_<arm>.pdf`, one per single unit, the same kind for
every arm) show the gate's window in the wobble arms. On w35's unit 0
(TT1, stage-1 units 0, 2, 3 and 80, 1.31 M spikes) the interval histogram
is empty below about 1.5 ms and peaks at 2.5 ms, and the neighbour
correlograms have a peak at zero, dips at plus and minus 1 ms and peaks at
plus and minus 2.5 ms. June's unit 0 (stage-1 units 0, 1 and 2, 397 k
spikes) has its own central structure (its large-against-small
correlogram is empty from -2 to +1 ms as w35's is), so the pages have to
be read in pairs; the amplitude-split ratio does not separate the arms
(medians 0.091 June, 0.110 w35).

### The gate hollows every correlogram, and the merge follows it

Wobble's near-zero refractory contamination and its CCG verdicts come from
the cosine gate, not from the matcher (`07_center_hollowing.py`,
`08_gate_or_matcher.py`, `center_hollowing.png`):

- Within units at stage 1, intervals under 1 ms: 2.58 per 1000 spikes in
  June, 0.12 in w35, 0.06 in w35a. Between same-tetrode units, lags from
  0.25 to 2 ms number 0.58 of those from 2 to 3.75 ms (the same width) in
  June and 0.28 in w35 (0.25 in w35a), with a step at 2 ms, the end of the
  -1 to +2 ms window the cosine is computed over: the raw snippet of a
  spike with a neighbour inside that window scores low and the spike is
  deleted. Within 0.25 ms the wobble arms have an excess (13.6 per 1000
  spikes in w35, 4.7 in w35a, 1.1 in June): one event assigned to two
  units.
- The 110 same-tetrode pairs that June's stage-1 trains show firing
  independently (central/flank over 0.7, median 0.996) have a median ratio
  of 0.445 in w35 and 0.408 in w35a; 27 and 33 of them fall below 0.3, the
  arbiter's duplicate line. The circus re-run: 1.001, none.
- Before the gate, on the same windows and banks, wobble's pairs have a
  median ratio of 0.928 with 13 below 0.3 (circus-omp: 0.612, 55). After
  the gate, 0.353 and 99. The same gate applied after the fact to circus's
  output takes it to 0.145 and 189.
- The merge accepts a pair when its CCG reads duplicate or its templates
  reach cosine 0.95. Of w35a's 38 merges, 18 join pairs that w35a's own
  ungated output shows firing independently and 10 join pairs that June's
  trains show so (8 both); w35 has 12 and 5, w25 17 and 9, June and the
  circus re-run none (`results/merges_vs_ccg.csv`). In w35a, unit 80 (a
  window-60 re-seed on TT1) absorbs units 0, 1, 2 and 3.

The refractory and CCG columns above are therefore not comparable between
June and the wobble arms: the gate removes the coincidences both measures
count.

## Wall-clock

Stage 1, each arm alone on windows 55 to 57 (16 jobs; the other agent's
work kept the load average between 7 and 19; `span_<arm>_timing/timing.json`):

| arm | matching | gate (per-spike cosine) | re-estimation | per window |
|---|---|---|---|---|
| circus-omp | 380 s | | 8.5 s | 130 s |
| wobble, w35 | 170 s | 77.5 s | 10.1 s | 86 s |

Over the 12 windows with all four arms running at once, plus the other
agent's jobs (load average 50 to 108): circus 45.4 min (matching 41.5,
re-estimation 0.6), w35 32.5 min (20.5, gate 10.6, 0.9), w25 39.3 min
(25.7, 12.5, 1.0), w35a 29.0 min (22.3, 5.6, 0.8). The circus arm also
computed per-spike cosines for this study (3.2 min), not counted. June's
own circus-omp run left only directory timestamps: the 12-window segments
of `reseed_sorts_c12`, each 12 windows of matching and re-estimation plus
one MountainSort5 re-sort, took 24 to 35 min, 29 min for windows 49 to 60
and 33 min for 61 to 72 (`results/timing.json`), about 2.0 to 2.9 min per
window on a machine whose load then is unknown. Downstream, at the
concurrent load: merge 2.4 to 3.3 min, residual pass 0.9 to 1.6 min,
multi-unit bucket 1.5 to 3.1 min per arm; the scoreboard 11 to 40 min,
longest for the arms with the most spikes.

The span's products take 11 GB on `/nvme` for all arms (the w35 stage-1
sorting 445 MB, its downstream 2.3 GB).

## Readings (conjecture)

- Conjecture: the shrinking bank in w35 and w25 is the cosine gate and
  the re-estimation together. The gate admits a spike of any size if its
  shape fits; re-estimating from those spikes makes the template the mean of
  a smaller population; the threshold, a fraction of the median template
  energy, falls with it; and the next window admits smaller spikes again.
  w35a, which re-estimates from the spikes circus-omp would have kept,
  holds its bank, which fits this reading.
- Conjecture: part of what w35 adds to a unit is low-amplitude activity
  of cells without a template of their own, which June leaves to the
  multi-unit bucket (that stage receives 64 % fewer events in w35). On
  the units where the matchers differ most the new spikes' footprint
  differs from the template's (unit 20, 0.90), but over all units their
  mean is a smaller copy of the template, which selection by the gate
  would produce either way; it is not shown whose spikes they are.
- Conjecture: the spikes wobble moves from a unit to a neighbour (14 % of
  June's, at full size and with the unit's waveform) follow from the
  shrinking bank: a neighbour whose template has shrunk and drifted
  toward the unit's shape explains a full-size spike of the unit at a
  lower energy cost.
- Conjecture: the higher tight purity of w35a is at least partly the
  merge. Fewer, larger units leave fewer neighbours to out-match a spike,
  and 18 of its 38 merges join units that fire independently in its own
  ungated output. On large units (19 in w35a against June's 27) its
  purity equals June's while its within-cluster spread is wider.
- Conjecture: June's high refractory contamination (median 0.175, 47 of
  74 single units at 0.1 or more) is itself partly an artifact of
  circus-omp placing two spikes of one unit within 1 ms; it is not evidence
  for the wobble arms, whose zero comes from the gate.

## What would decide it

The span says wobble, as June configured it (threshold 0.35 or 0.25 of the
median energy, raw cosine gate 0.60, per-window re-estimation), loses:
lower assignment purity on large and on all units, wider within-cluster
spread on large units, 14 % of June's spikes moved to a neighbour, a bank
that shrinks window by window, and merges of cells that fire
independently. Its two measured advantages, coverage of
the 10 MAD events by single units and the near-zero refractory
contamination, are the same events moved out of the multi-unit bucket and
an artifact of the gate.

Whether wobble itself could win is undecided. Ungated, its correlograms
are fuller than circus-omp's (median ratio 0.93 against 0.61), so the
matcher does not hollow them; the gate does. What would decide it, on the
same span and the same bank, all scored as here:

1. The gate scored on the residual after wobble's own subtraction (or with
   the peel `ttm.capture` uses), so a spike with a neighbour in its window
   is judged on its own waveform, and with the CCG arbiter then reading
   non-hollowed trains. A run that keeps the gate as it is cannot be
   judged by the merge or by refractory contamination.
2. Re-estimation anchored to full-size spikes, as in w35a, or the
   threshold fixed at the span start, so the bank cannot shrink.
3. A per-template threshold (or the cosine gate alone with a low fixed
   admit) for the small templates that the median-energy threshold
   penalises.
4. Then the same tables, and the purity pages read by the user, with the
   merge's acceptances checked against ungated correlograms.

The full 95-window pass was not run and is not proposed on this result.
For scale: stage 1 alone ran at 86 s a window for wobble against 130 s for
circus-omp, about 2.3 against 3.4 h over 95 windows, plus seven
MountainSort5 re-sorts in either case.

## Files

In the repository (`investigations/wobble-fullsize-ttm001/`): `_span.py`,
scripts `01` to `09` (bank and re-seed reconstruction, stage 1, downstream,
stage-1 comparison and galleries, purity tables and pages, summary,
correlogram centres, gate against matcher, shape of what wobble adds) and
`results/`. The purity
pages use `ttm.purity`, lifted in this work from
`subjects/ttm-003/scripts/sorting_purity.py`, which now imports it.

On the data volume, under `.../track_eval/wobble_span/`: `bank_w55.pkl`,
`reseed_w60.pkl`, `template_amp_w55.npy`, `span_<arm>/` (sorting, raw
per-window matcher output, timing), `span_<arm>_timing/`, `down_<arm>/`
(stage-1 input, merge, single units, deliverable, scoreboard, merges,
purity scores) and the logs. Figures, under
`/nvme/neuropixels/tetrode_data/TTM-001/TTM-NOD/figures/wobble_span/`:
`stage1_overview.png`, `stage1_gallery_<arm>_<class>_u<unit>.png`,
`center_hollowing.png`, `purity_june.pdf`, `purity_w35.pdf`,
`purity_w25.pdf`, `purity_w35a.pdf` (one page per single unit).
