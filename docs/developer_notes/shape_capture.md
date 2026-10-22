---
title: Shape-based spike capture (ttm.capture)
updated: 2026-10-22
status: needs_verification
---

What `ttm.capture` does, the choices that were measured on the TTM-003
optogenetic sessions of 2026-10-16, and what it does not do. The subject
driver is `subjects/ttm-003/scripts/opto_capture.py`; the numbers behind
each choice are in `subjects/ttm-003/optogenetic_sessions_findings.md`,
section 8, and `subjects/ttm-003/results/capture_*.json`.

## The problem it answers

MountainSort5 decides which events exist with a 5.5 MAD threshold on the
whitened signal and which cluster they join with amplitude-weighted PCA
features. On TTM-003 that drops the TT6 cell's in-pulse spikes (half
amplitude, same shape) and never isolates the TT10 cell (about 4.3 MAD).
The module keeps the sorter for one job only, proposing templates, and
redoes detection and assignment:

1. Noise per wire with spikes excised (`noise_mad_excised`): samples within
   1 ms of any crossing below 4 MAD on any wire are left out of the median,
   twice from the plain estimate. On a 32 min tetrode the median over every
   tenth sample differs from the full one by under 0.15 %, so the estimate
   is subsampled.
2. Detection (`detect_troughs`): every sample at which the four-wire minimum
   (each wire in its own noise units) is below 3.5 and is the minimum within
   0.3 ms on either side. One event per trough; the 0.3 ms is the detector's
   dead time and is passed to the arbiter below.
3. Scoring (`score_events`): the cosine of the event's four-wire snippet over
   the full template window (-1.0 to +2.0 ms, the window the June work kept
   after testing tighter ones) to every same-tetrode template, maximised over
   a one-sample shift either way. Amplitude plays no part. The best template's
   amplitude scale `a = conv / |t|^2` and a noise-corrected cosine (the
   snippet norm with the expected noise energy removed) are kept per event.
4. Assignment with peeling (`capture_tetrode`): the best template takes the
   event if its raw cosine reaches the gate (0.6). Accepted templates, scaled
   by `a`, are subtracted from a working copy of the traces and the remaining
   events are scored again, three passes in all. Without the peel a spike
   1 to 2 ms from a larger one sits inside a contaminated 3 ms snippet, fails
   the gate, and every neighbour cross-correlogram is hollow at the centre
   (every pair reads as a duplicate). With it, about 13 % of the accepted
   events on TT6 come from the second and third passes.
5. Refractory rule: among the events one template accepted, those closer
   than 1.5 ms form a chain and only the best cosine survives. The number of
   events a template lost this way, relative to what two independent trains
   of its size would produce, is the refractory conflict ratio; it replaces
   refractory contamination, which the rule makes zero by construction.
6. Everything not accepted is the tetrode's pool: unexplained events for a
   multi-unit bucket or a new template, and noise.

`match_wobble` runs SpikeInterface's wobble matcher on the same four
channels and templates (objective threshold 0.35 of the median template
energy, the June admit level) so the detector-plus-assigner alternative can
be scored on the same traces with the same gate.

## Purity

`unit_purity` and `plot_unit_purity` give, per unit: refractory
contamination (Llobet), the refractory conflict ratio, the amplitude-split
cross-correlogram (large versus small spikes of the same unit outside the
stimulation; one cell gives a zero-lag ratio near zero, two cells near one,
and a 2 to 3 ms peak means the small spikes follow the large ones in
bursts), cross-correlograms with every same-tetrode unit with the June
arbiter's verdict (duplicate below 0.3, distinct above 0.7, lags inside the
detector's dead time excluded from the central window), waveform overlays
on the template, amplitude against cosine and amplitude over time. No scalar
"purity bar" is encoded; the user's criterion is still an open question
(`investigations/spike-sorting/OPEN_CONCERNS.md`).

## Measured on TTM-003 (2026-10-22)

See the findings document for the tables. The short version: the
noise-corrected cosine and any raw gate above 0.6 discard the TT6 cell's own
large spikes on this busy tetrode (a 3 ms window holds other events even
after the peel), so the raw gate at 0.6 is the operating point; the cell's
in-pulse train is captured at 144 per second with an amplitude-split ratio
of 0.04; three small-template units on the same tetrode take 4 to 6 per
second of in-pulse events each, which the arbiter marks as ambiguous
collectors.

## Limits

- No new templates: a cell with no seed template stays in the pool. Seeds
  can come from a lower-threshold re-sort of one tetrode (`--seed`) or from
  another session's templates (`--templates-from`).
- Two troughs closer than 0.3 ms are one event; the peel recovers the
  second spike only when it is farther than that.
- The gate is noise-limited below about 5 MAD: a small spike of another
  shape and a small spike of the template's shape are not separable by
  cosine there, so small-template units collect hash. The refractory
  conflict ratio, the amplitude-split test and the neighbour verdicts are
  what show it.
- Runtime on tononi-2, one tetrode of a 32 min session (546 k events):
  load 6 s, noise 15 s, detection 8 s, three scoring passes 25 s, purity
  3 s; wobble on the same tetrode is reported in the findings table.
