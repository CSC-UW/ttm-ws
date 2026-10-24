---
title: Small-spike capture on a full-size sorting's residual (ttm.staged_capture)
updated: 2026-10-23
status: needs_verification
---

What `ttm.staged_capture` does, which settings were measured and where, and
what it does not do. It is stage 3 of the staged pipeline of
`docs/plans/2026-10-22_ttm003_48h_sorting_plan.md` (track F). The driver is
`subjects/ttm-003/scripts/mp_staged.py`; the settings were measured with
`subjects/ttm-003/scripts/staged_settings.py` and are reported in
`docs/reports/2026-10-23_ttm003_48h_staged_pipeline.md`. It reuses the
scoring, peeling and refractory pieces of `ttm.capture`
(`docs/developer_notes/shape_capture.md`).

## What it does, per window and tetrode

1. `subtract_sorted` fits every spike of the full-size sorting (the June
   matcher plus the June residual pass) with its unit's window template and
   subtracts it. The fit is the least-squares scale `a = <snippet, t> / |t|^2`,
   floored at zero, at the shift within one sample that maximises the inner
   product, computed on the residual as it stands, deepest spike first. On
   the 48 h windows the median fitted scale per window and tetrode is 0.94
   (0.87 to 0.96 for 80 % of them), and the subtraction of one tetrode
   window of 1800 s takes seconds; with detection and the offer, a median
   of 75 s.
2. `residual_events` detects the residual's troughs below the floor times
   the spike-excised noise of the original traces
   (`ttm.capture.noise_mad_excised`, so the floor is the same number of noise
   units as in `ttm.capture`), and drops those within 0.3 ms (the detector's
   dead time) of a sorted spike, which are the remains of an imperfect fit.
3. `offer_events` scores every event against every template of the tetrode
   (full-window cosine, -1 to +2 ms, amplitude free) and offers it to the
   best one when the noise-corrected cosine reaches the gate. Offered events
   are subtracted from a working copy and the rest scored again, three passes
   in all. Every template competes, admitted or not, so an event that looks
   most like a unit that is not admitted stays unassigned.
4. `light_response` counts, per pulse, a unit's spikes from 10 ms after
   onset to offset ("in pulse") and in a window of the same length ending at
   onset ("pre-pulse"; inside a train this lies between pulses and is not a
   baseline), and gives the rate ratio, a one-sided exact test that the
   in-pulse rate is higher given the total count, and the per-pulse sign
   test. `admit` applies the rule: ratio at least `min_ratio`, count test
   below `max_p`, in-pulse rate at least `min_in_hz`.
5. `governing_train` and `accept_offered`: an offered event is kept only if
   its unit was admitted at the most recent delivered train before it
   (before the first train, the first train decides), is farther than
   1.5 ms from every sorted spike of that unit, and has the best raw cosine
   in its chain of the unit's other new events closer than 1.5 ms.

The light test runs on the unit's sorted spikes plus the events it was
offered in that train's window as if it were admitted, so a template whose
cell fires small spikes through the pulse passes even when the full-size
sorting kept only the first spike after onset.

## Settings chosen (2026-10-23)

Defaults of `mp_staged.py`, measured on the two October 16 sessions and on
48 h blocks 1 to 4 and 25 to 27:

- floor 3.5 excised MAD. At 4.5 the TT6 cell keeps 26 to 27 % fewer
  in-pulse spikes over the 4.2 mW 207 ms pulses (99.6 against 136.1 per s
  in OPTO1, 79.5 against 106.7 in OPTO2) while its refractory conflict
  ratio changes by under 0.01;
- gate 0.8 on the noise-corrected cosine, the lowest value at which the TT6
  cell's stage-3 spikes are at least as clean as with the candidate's raw
  cosine 0.6 by the refractory conflict ratio and by the rate added outside
  the trains, in both sessions;
- light test: in-pulse rate at least twice the pre-pulse rate and the count
  test below 0.001, no rate floor and no template floor.

## Limits

- No new templates: a cell without a unit in the full-size sorting stays
  unassigned, as in `ttm.capture`.
- The cross-correlogram of a unit's stage-3 events against its own sorted
  spikes is depleted within about 1.5 ms by the fit residue in the scored
  window, against other units' spikes as much as against the unit's own,
  so it cannot be read as the refractory signature of one cell; the
  refractory conflict ratio and the amplitude-split test of the whole unit
  are the measures used.
- A unit whose own full-size spikes are silenced by the light can still pass
  the light test if the events it is offered in the pulse come from a driven
  neighbour whose shape is closer to its template than to the neighbour's.
  On the 48 h this is the reading the data allow for TT6 units 14 and 15
  (findings and report); the test cannot tell the two cases apart.
- Stage 3 runs on whole windows, so most of what it adds lies between the
  pulses: on the 48 h, 0.14 M of 34.5 M stage-3 spikes are inside a pulse,
  and the refractory conflict ratio of the added spikes is 0.08 to 0.68
  across the 28 admitted units (median 0.27). `mp_staged.py --pulses-only`
  keeps only the in-pulse ones, leaving every unit its full-size self
  elsewhere.
- The pre-pulse rate rises with the animal's state and with noise; at a
  train in a noisy stretch (block 27, 21:11, after lights off) fewer units
  pass and those that fail lose stage 3 until the next train.
