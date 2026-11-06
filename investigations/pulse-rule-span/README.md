---
title: The pulse rule's span after the offset
updated: 2026-11-06
---

# The pulse rule's span after the offset

The carry's pulse rule (`docs/sorting_pipeline.md`, step 5.4) admits a
sustained responder's fits under lower floors inside each laser pulse. A
sustained responder keeps firing through the pulse; a transient responder
answers only its onset. Until 2026-11-06 the rule's span ran from onset to
100 ms after the offset, the same span that keeps pulse times out of the
relative floor's reference and the template updates. The 100 ms entered
the carry bake-off as that exclusion
(`investigations/carry-bakeoff/_bake.py`, `post_ms`), and the pulse rule
reused it. The bake-off varied the rule's
floors but not its span, and its hybrid's sustained responders fired
half-size spikes only inside the pulse. The question: what did the 100 ms
after the offset admit in TTM-003-48H (207 ms pulses at 9.8 mW), which
was sorted with it?

`post_offset.py` reads the saved fits of TTM-003-48H's carry
(`sortings_carry/main/carry/`) and, for each window's sustained
responders, splits the accepted spikes by time from the pulse into those
that pass the normal rules (snr 4.5, bad periods, cosine 0.5, relative
floor 0.5 of the unit's median; the refractory rule aside) and those only
the pulse rule accepts (relaxed-only). `readmit_check.py` replays
admission on every window with `ttm.carry.readmit_window` without writing
anything. `template_effect.py` counts where the fits the replay dropped
could have changed a template update. The first two ran before
TTM-003-48H was replayed (`ttm-sort readmit`); since then the carry's own
decisions are in each window's `reason_carry` and its parameters under
`carried_with` in `carry_params.json`, and the scripts read them there.

## Findings

- Both kinds of spike stop together after the offset. For TT4's unit 7,
  over all 1,432 pulses, in 2 ms bins after the offset (spikes passing
  the normal rules, then relaxed-only): 234 and 159 from 0 to 2 ms, 31
  and 40 from 2 to 4 ms, 2 and 6 from 4 to 6 ms, and after 6 ms one
  relaxed-only fit (8 to 10 ms) and one spike passing the normal rules
  (14 to 16 ms) up to 30 ms (`offset_profile.csv`). Just before the
  offset it has about 105 and 56 per second.
- The units then stay nearly silent. Unit 7's spikes passing the normal
  rules come at 0.07, 0.12 and 0.24 per second from 10 to 20, 20 to 50
  and 50 to 100 ms after the offset (0.07, 0.12 and 0.34 with its
  relaxed-only fits), against 15.3 per second more than 2 s from any
  pulse.
- From 10 to 100 ms after the offset, over the 48 h, the rule admitted
  105 relaxed-only fits. Every sustained responder with a spike there:

  | Unit | Tetrode | Windows as a sustained responder | Passing the normal rules | Relaxed-only |
  |---|---|---|---|---|
  | 29 | TT11 | 24 | 24 | 23 |
  | 7 | TT4 | 24 | 23 | 7 |
  | 26 | TT10 | 13 | 129 | 33 |
  | 30 | TT11 | 9 | 11 | 4 |
  | 16 | TT6 | 4 | 1 | 3 |
  | 2005 | TT5 | 4 | 1 | 0 |
  | 1 | TT2 | 3 | 197 | 27 |
  | 44 | TT16 | 2 | 40 | 0 |
  | 43 | TT15 | 1 | 38 | 7 |
  | 0 | TT1 | 1 | 70 | 1 |
  | 28 | TT10 | 1 | 17 | 0 |

- Whether these fits are the units' spikes cannot be told from these
  data. Their cosines (0.50 to 0.58 for the first five units) are like
  those of the relaxed-only fits inside the pulse (medians 0.52 to 0.54).
  Their share of each unit's spikes from 10 to 100 ms (20 to 75 % for
  the first five) is near its share from 15 ms after onset to the offset
  (25 to 65 %). For units 7 and 26, which have most of the spikes there,
  it is also near the share of the unit's fits more than 2 s from any
  pulse that pass only the lower floors (23 and 20 % after the offset,
  20 and 19 % away from pulses). What the rule does after the offset is
  admit spikes there under a lower bar than in the rest of the time
  outside pulses, so rates after the offset and away from pulses would
  be measured differently.
- Inside the pulse, relaxed-only fits come at 5.2 to 15.5 times the rate
  of the same kind of fit (passing the pulse rule's floors, failing the
  normal ones) more than 2 s from any pulse for the units sustained in 4
  or more windows, and at 1.1 to 3.4 times for the five sustained in 1 to
  3 windows. They start and stop with the spikes that pass the normal
  rules. For the TT6 cell they are 65 % of its accepted spikes from 15 ms
  after onset to the offset.
- Replayed with the parameters the carry ran with, `readmit_window`
  reproduces all 845.9 M saved decisions except 61 fits, all in window
  28 and all inside a bad period of their tetrode: the pulse rule had
  accepted them because its candidates excluded bad periods only for fits
  of snr 4.5 or more (`readmit_check.csv`). With the span ending 5 ms
  after the offset, 112 more fits are dropped and none added.
- The dropped fits reach the template updates in one place
  (`template_effect.csv`). Of the 41 unit-windows with a dropped fit,
  none crosses the 100 accepted spikes an update needs. In one (unit 7,
  window 48) a spike just after the excluded span becomes a lead spike,
  and the pool the update draws 500 spikes from grows from 29,302 to
  29,303. A carry run with the 5 ms span would draw a different 500
  there, so its unit 7 template from window 49 on would differ slightly
  from the replayed run's.

## Outcome

The span splits in two (`ttm.carry.CarryParams`): the exclusion keeps
100 ms after the offset (`pulse_exclude_post_ms`), the pulse rule ends
5 ms after it (`pulse_relax_post_ms`), and `pulse_relax` (`ttm-sort
--no-pulse-relax`) keeps the exclusion without the rule, for pulses under
which no cell's spikes shrink. The 5 ms rests on the units' silence after
the offset, not on the fits there being false. The pulse rule's candidates
exclude bad periods for every fit. TTM-003-48H was replayed with
`ttm-sort readmit` on 2026-11-06, and TTM-003-48H-2 is to be sorted with
`--no-pulse-relax`.
