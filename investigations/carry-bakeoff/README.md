---
title: Carry-forward bake-off on hybrid ground truth, TTM-003 48 h
updated: 2026-10-28
status: needs_verification
---

# Carry-forward bake-off

Two questions on the first 12 windows (6 h) of TTM-003's 48 h recording:
can circus-omp run on each tetrode on its own match or beat circus-omp on all
channels at once (June), and how should a carry-forward rebuild its
templates from one window to the next? The per-tetrode carries of
`investigations/per-tetrode-matching/` ran away on tetrodes with few or small
units. That README names two things not tried: a floor in noise units
instead of a fraction of the template, and re-estimation from full-size
spikes only.

## Ground truth

`make_hybrid.py` adds 21 known units to a copy of the first 12 windows of
`mp_baseline/binary` (`carry_bakeoff/v2/hybrid_binary` on the session's
sortings, 155 GB).
- The units are seed-bank templates moved to other tetrodes with their wires
  permuted, at 5 to 14 MAD; two are mixes with a destination unit, at cosine
  about 0.87.
- Each has rate swings, bursts at 0.85, 0.7 and 0.6 of full size, a slow
  change of size and wire ratios, size jitter and sub-sample timing.
- Two stop firing after window 5.
- Three answer the real laser pulses with a full-size spike, then a train at
  half size (the TT6 cell's pattern). The other units fall silent through each
  pulse.

`gt.npz` and `gt_plan.json` record every spike, its class and size factor,
and the injected template per window. A first version (v1) had driven units
of 5 and 6 MAD, whose half-size spikes no rule could find; it was replaced.

## Code

- `_bake.py`: the carry loop.
  - Matching: circus-omp, joint or per tetrode, with a per-template amplitude
    floor.
  - Per spike: fit amplitude; matched-filter signal-to-noise (projection on
    the unit-norm template over that projection's noise, measured per
    template and window); cosine of the raw snippet and of the peeled
    snippet (residual plus the spike's own fit); residual energy; trough.
  - Admission rules and re-estimation recipes, with optional spatial
    whitening per tetrode.
- `run_arm.py`: the arms. Every arm starts from the same bank: mp_baseline's
  seed bank plus one template per ground-truth unit, from its window-0
  spikes.
- `admission.py`, `admission_units.py`: stage A. Rules applied to permissive
  per-tetrode fits on frozen templates, against the joint runs.
- `score.py`, `compare_carries.py`: stage B. Carries over 12 windows, scored
  per window on the ground truth and on proxies for the real units.

Paths are literal tononi-2 paths.
