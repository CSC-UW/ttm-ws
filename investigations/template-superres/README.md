---
title: Sub-sample template fitting in the carry, TTM-003 hybrid
updated: 2026-11-02
---

# Sub-sample template fitting

The carry fits each template at integer sample shifts and rebuilds it from
the mean of integer-aligned snippets. This asks how much fitting at the
best fractional shift would change the scores admission uses, the
residual and the rebuilt templates, before building it into a matcher.

## Code

- `ceiling.py`: for arm `W_s45_c50_lnb_b2` on the aligned-reference hybrid
  with the re-sorted seed bank (`BAKE_REFERENCE=aligned_seed`, see
  `investigations/carry-bakeoff/`), windows 3 and 8, every tetrode:
  rebuilds the residual of all fits, takes each fit's peeled snippet and
  evaluates the template at delays of -0.5 to +0.5 samples (FFT shift,
  step 1/16, parabolic refinement). Also the mean of the snippets aligned
  by their own best delay and, for ground-truth units, by the injected
  delay. Up to 20,000 fits per unit and window, reasons 0 to 2.
- `summarize.py`: the tables in `results/summary_w03_w08.txt`.
- `fp_near_large.py`: ground-truth false positives against hits near a
  large unit's accepted spike.

Per-fit tables went to the session scratchpad and are not kept;
`results/units_w03_w08.csv` is per unit.

## Results

2.59 M fits scored, 585 k accepted, 67 units (21 ground truth). The
integer-shift snr and cos_pl recomputed this way match the stored values to
1e-6.

| accepted fits, median [p90] | small units (6 MAD or less) | large units |
|---|---|---|
| snr, relative gain | 0.3 % [1.8 %] | 0.25 % [1.2 to 1.4 %] |
| cos_pl gain | 0.002 [0.011] | 0.002 [0.010] |
| residual energy drop, ground-truth hits | 0.6 % | 1.1 % [7.7 %] |

Fits that pass snr >= 4.5 and cos_pl >= 0.5 at the best delay but not at
the integer one: small ground-truth units gain 208 hits on 20,774 (1.0 %)
and 69 false positives on 1,102 (6.3 %), so 75 % of what would be added
are hits against 95 % of what is admitted now; large ones gain 185 hits
(0.1 %) and 431 false positives (4.9 %).

The best delay correlates with the injected quarter-sample offset at 0.72
for large units and 0.50 for small (error 0.29 and 0.52 samples against a
spread of 0.39 and 0.60). Aligning snippets by their best delay raises the
mean's energy by 0.8 to 0.9 %, no more than the same procedure gives on a
template plus noise at the exact time (0.7 to 1.5 %); aligning by the
injected delay gives 1.1 % for large units and 0.4 % for small, with
cosine to the truth changed by 0.002.

On large units (trough 10 MAD or more, 181,679 fits) the share of fits
whose residual peaks at 5 noise units or more falls from 15.3 % to
12.3 %. Ground-truth false positives are 2 to 4 times more often than hits
within 60 samples of a large unit's spike; how much of that is jitter left
in the residual is not established.

Sub-sample refinement of accepted spikes for gating and the rebuild costs
about 8 s per tetrode and window and would add mostly false positives at
the present floors. Only fitting at sub-sample shifts inside the matcher
(a jittered bank in circus-omp, or wobble's jitter) could use the residual
effect on large units; not costed.
