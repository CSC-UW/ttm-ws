---
title: Bad periods on TTM-003 48 h, the diagnostics behind the report
updated: 2026-10-27
---

# Bad periods on TTM-003 48 h

One-off checks run on 2026-10-26 and 27 for
`docs/reports/2026-10-26_ttm003_bad_periods.md`, which compares the June
pipeline on TTM-003's 48 h binary with and without `ttm.bad_periods`'s mask.
The comparison itself is `subjects/ttm-003/scripts/bad_periods_compare.py`;
these scripts answer the questions it raised. Paths are literal `tononi-2`
paths as run. Each script's printed output is the `_output.txt` beside it.

| Script | Question | Answer |
|---|---|---|
| `window1_coverage.py` | Do TT5 to TT12, which have no flagged time in window 1, lose spikes in the baseline while the even-channel tetrodes carry the 2287 to 2696 s event? | Yes. Their spikes fall from about 25,000 to 2,000 to 5,000 a minute and the share of their events at 10 MAD or more that a unit claims from 95 % to 42 to 60 %; masked it stays at 95 to 97 %. The template-seed control behaves like the baseline. |
| `good_tetrode_coverage_by_flags.py` | Does any bad period elsewhere cost TT5 to TT12 spikes, or only that event? | Only when five to eight other tetrodes are flagged (408 s, all in that event): 60 % claimed against 96 % masked. With one to four others flagged the three runs claim 94 to 95 %, as with none. |
| `event_coverage.py` | Does the second even-channel event (96,900 to 97,560 s) do the same in the full 48 h runs? | Worse. From 97,260 to 97,500 s TT5 to TT12 have 1 to 4 spikes a minute and 0.0 to 0.2 % of their events claimed in the baseline, against 23,600 to 30,400 spikes and 96.9 to 97.7 % masked. Before and after the event the two runs claim within 1.2 points of each other. |
| `window_coverage.py` | Per tetrode and window over 48 h, how much of the clean-time events at 10 MAD or more do the baseline, the masked run and an unmasked control with another template seed claim? | At the second even-channel event both unmasked runs claim 0.0 to 3.7 % on TT1, TT13, TT15 and TT16 for windows 54 to 60, masked 5 to 100 %. Over 48 h the masked run claims at least as much as the control on every tetrode but TT5 (83.3 against 85.6 %), lower there in every window from 56 on. |
| `omp_coupling.py` | Is the loss caused by the matcher rather than the data? | Yes. On 2400 to 2460 s with window-1-like templates, circus-omp finds 4,943 TT5 to TT12 spikes unmasked, 33,372 masked, and 33,432 unmasked when only the TT5 to TT12 templates are matched. Unmasked, the spikes that remain are three times as frequent in the first tenth of each 1 s chunk as in the last. |
| `tt15_seed_stability.py` | Is the TT15 unit missing from the masked run a cost of the mask? | TT15's three-cluster seed sort is fragile: its smallest cluster comes out on its own in 2 of 8 unmasked re-sorts with other whitening seeds or a shifted window, and in none of 5 masked ones. |
| `tt15_clock_lock.py` | What are TT15's units? | All of them, in both runs, fire at a fixed phase of a 300-sample cycle. The cycle mean of CH22 is 10 MAD peak to peak and keeps its trough at sample 12 or 13 from 0.3 to 47.2 h, so its source is synchronous with the acquisition clock. |

`bad_periods_compare.py` reports the last check for every unit of a sorting
as `clock_lock`, and the coverage of detected events in clean and flagged
time as `results/bad_periods_coverage_<E>_<tag>.csv`.
