---
title: Wobble in place of circus-omp for the full-size stages, TTM-001 span test
updated: 2026-10-23
---

# Wobble for the full-size stages on TTM-001

Track G of `docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`, items 1 to 3:
the June TTM-001 pipeline with wobble in place of circus-omp, over 12 of the
95 windows of session `2026-05-27_09-07-52`, judged head to head against
June's assembly on the same windows. The report is
`docs/reports/2026-10-23_ttm001_wobble_fullsize_span.md`.

Paths in the scripts are literal `tononi-2` paths. Bulk products are under
`/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52/sortings_seed42_pcafix/track_eval/wobble_span/`
and figures under
`/nvme/neuropixels/tetrode_data/TTM-001/TTM-NOD/figures/wobble_span/`;
nothing under the June directory `mp_long_s2000_d170000/` is written.

## Span and settings

- Windows 55 to 66 of the June carry-forward (1800 s each, 28.1 to 34.1 h
  of the recording). June re-seeded after matching windows 12, 24, ..., 84
  (`reseed_sorts_c12/reseed_w*`, births in `reseed_cadence_sweep.npz`);
  window 60's re-seed added ten units, the most of any boundary, so the
  span has six windows on each side of it.
- The bank entering window 55 is rebuilt from June's own spike trains
  (each unit's template from its last window with at least 100 spikes,
  computed as June's loop did), and the re-seed units are June's rule
  applied to June's MountainSort5 re-sort of window 60. The rule gives the
  same ten units with the same at-birth cosines to four decimals.
- Arms: `circus` (amplitudes >= 0.8, June's matcher, re-run as a check),
  `w35` (wobble, objective threshold 0.35 times the window bank's median
  template energy, cosine to the assigned template >= 0.60; June's
  documented operating point), `w25` (the same with threshold 0.25, the
  lowest admit June mapped, which asks less of small templates) and
  `w35a` (`w35` with templates re-estimated only from spikes at amplitude
  scale 0.8 or more; a diagnostic added when `w35`'s bank shrank).
- Downstream unchanged (`97_production_pipeline.py` imported with its
  output directory redirected): CCG-guarded merge, residual single-unit
  pass, per-tetrode multi-unit bucket, three-axis scoreboard, all on the
  span's windows. The residual pass has no matcher in it.
- Amplitude classes: a stage-1 template is large when its trough on its
  deepest wire at the span start is at least 10 excised MAD (37 units),
  small otherwise (53).

## Outcome

Wobble as June configured it loses on this span: assignment purity 0.855
against June's 0.884 (tight window, single units; a re-run of June gives
0.886), its bank's median trough falls from 103 to 72 µV in six windows,
and the raw-snippet cosine gate deletes coincident spikes, which hollows
every same-tetrode correlogram (median central-to-flank ratio over all
pairs: 0.61 for circus-omp, 0.93 for wobble before the gate, 0.35 after
it) so the CCG-guarded merge joins independent cells. A gate scored on the residual
with anchored re-estimation is untested. Details and numbers:
`docs/reports/2026-10-23_ttm001_wobble_fullsize_span.md`.

## Files

| file | purpose |
|---|---|
| `_span.py` | paths, span, arms; the June loop body over a list of windows; bank and re-seed reconstruction; amplitude classes |
| `01_span_bank.py` | bank entering window 55 and the window-60 re-seed units, checked against June's births |
| `02_span_carry.py` | stage 1 over the span for one arm; per-window raw matcher output, timing |
| `03_span_downstream.py` | June's merge, residual, multi-unit and scoring stages over the span for one arm |
| `04_stage1_compare.py` | stage 1 per unit against June: shared, new, dropped, reassigned; galleries; bank drift |
| `05_purity_pages.py` | purity tables and pages (`ttm.purity`) of each arm's single units; pairing across arms |
| `06_summary.py` | scoreboards, coverage per tetrode, per-class aggregates, timing |
| `07_center_hollowing.py` | same-tetrode correlogram centres and short intervals at stage 1, per arm |
| `08_gate_or_matcher.py` | the same before and after the cosine gate on the raw matcher output; merges checked against ungated correlograms |
| `09_new_spike_shape.py` | waveform, footprint and size of what wobble adds and moves, every unit, window 61 |
| `results/` | the small tables these write |
