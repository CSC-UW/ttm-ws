---
title: TTM-003 48 h, the June baseline against the shape-capture candidate
updated: 2026-10-23
status: needs_verification
---

# TTM-003 48 h: June baseline against the shape-capture candidate

Track C and D of `docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`, run on
2026-10-22 and 23 on the 48 h recording `2026-10-20_09-02-12` (manifest
`TTM-003-48H`, 47.98 h, 24 delivered laser blocks: trains 1 to 4 with 58
pulses in the first 6 h and trains 25 to 44 with 60 pulses every 2 h from
8.15 h, all 207 ms at 0.995 Hz and 9.8 mW). Measurements are stated as
such; sentences marked *conjecture* are the analyst's reading. Everything
here is a decision dashboard for the user, not a verdict. Small results are
in `subjects/ttm-003/results/` (`tagged_blocks_*`, `tagged_match_*`,
`tagged_xsession_*`, `purity_*`), the sortings and analyzers under
`<session>/sortings_seed42_pcafix/{mp_baseline,mp_candidate,mp_candidate_xseed,mp_candidate_xfixed}/`,
the figures under `<data_root>/TTM-003/TTM-003-48H/figures/`
(`tagged_blocks_*.png`, `tagged_xsession.png`, `compare/purity_*.pdf`).

## 1. What was compared

Both sortings use the same preprocessed binary (bandpass 300 to 6000 Hz,
global common median reference, 1.33 TB float32) and the same 96 windows of
1800 s.

The baseline (`mp_baseline.py`) is the June TTM-001 production pipeline
unchanged: a MountainSort5 seed sort of the first window (55 units, 45 after
the snr and dedup bars), circus-omp carry-forward at amplitudes of 0.8 of
the template or more with per-window re-estimation and re-seeding every 12
windows (88 units, 43 of them born by re-seeding), then the CCG-guarded
merge (88 to 84), residual capture of unclaimed events of 10 MAD or more
(127 k recovered, 4.6 % of them), and the per-tetrode multi-unit bucket at
7 MAD or more (6.35 M events into 16 units): `assembled_prod`, 84 single
plus 16 multi-unit units, 85.5 M spikes in the carry-forward sorting.

The candidate (`mp_candidate.py`) keeps the baseline's windows, seed units
and identities and replaces detection and assignment only: per window and
tetrode, templates from the baseline units with at least 100 spikes there,
then `ttm.capture` (every trough at 3.5 times the spike-excised noise,
full-window cosine to the best same-tetrode template with the amplitude
free, gate 0.6, three peeling passes, the refractory rule), the leftover
events of a tetrode pooled into one multi-unit unit, then the same
CCG-guarded merge (81 non-empty single units to 70) and the same analyzer:
`assembled_candidate`, 70 single plus 16 multi-unit units. Templates with
a trough shallower than 40 uV are not captured, and a window with no
template above that floor is not detected at all. 447 M events were
detected, 156 M accepted by the 88 units (the baseline's 85.5 M) and 291 M
pooled.

Two further candidate variants (`--xseed`, section 2.5) add the October 16
templates of the four tagged or driven cells to their tetrodes, carried or
fixed.

## 2. The tagged cells block by block

`tagged_cells_blocks.py` matches each sorting's units to the OPTO2 capture
templates and measures every unit at every delivered block (baseline -50 to
0 ms, onset 1 to 10 ms, in pulse 10 ms to offset, post 0 to 50 ms after
offset, rates per second over the block's pulses).
`tagged_cells_xsession.py` is the control that does not depend on either
sorting: the traces around each block captured with the OPTO2 templates
themselves.

### 2.1 The TT6 cell (OPTO2 unit 22)

The cell is in the bank as the baseline's seed unit 16 (template trough
-112 uV, cosine 0.91 to the October cell template and 0.98 to the
silenced CH49 neighbour, OPTO2 unit 23, which shares its wire profile), and
it is present in the recording for the first 6 h only. The headline table
gives, at every delivered block, what each pipeline makes of it: the
baseline's unit 16, the candidate's unit 16 (same seed), the capture around
the block with the cell's own October template (`tagged_xsession`, no 48 h
sorting involved), and the cross-seeded variant's unit 16 and its extra
October template. Rates per second over the block's pulses.

| block | hour | baseline unit 16 base / onset / in pulse | candidate unit 16 base / onset / in pulse | October template alone onset / in pulse | cross-seeded: unit 16 in pulse, extra in pulse |
|---|---|---|---|---|---|
| 1 | 0.45 | 0.7 / 15 / 3.7 | 2.1 / 69 / 124.9 | 59 / 110.9 | 124.9, 0.0 |
| 2 | 1.98 | 3.1 / 13 / 0.1 | 5.9 / 90 / 94.7 | 84 / 94.8 | 94.7, 0.0 |
| 3 | 4.33 | 0.7 / 4 / 0.0 | 3.4 / 69 / 74.8 | 98 / 81.2 | 74.8, 0.0 |
| 4 | 5.96 | 2.1 / 4 / 0.0 | 9.3 / 56 / 70.5 | 50 / 73.1 | 70.5, 0.0 |
| 25 | 8.15 | 3.0 / 4 / 0.0 | 11.3 / 28 / 19.8 | 35 / 26.0 | 8.4, 35.5 |
| 26 | 10.15 | 5.0 / 2 / 0.0 | 16.3 / 26 / 10.6 | 17 / 13.7 | 7.7, 20.3 |
| 27 | 12.15 | 4.7 / 0 / 0.0 | 17.7 / 11 / 4.0 | 9 / 5.2 | 3.4, 5.6 |
| 28 | 14.15 | 5.0 / 2 / 0.0 | 16.7 / 17 / 4.9 | 4 / 6.4 | 4.7, 5.7 |
| 29 | 16.15 | 2.0 / 2 / 0.0 | 13.7 / 22 / 7.8 | 9 / 9.1 | 7.5, 6.1 |
| 30 | 18.15 | 2.3 / 0 / 0.0 | 17.3 / 9 / 2.5 | 6 / 2.8 | 2.4, 3.5 |
| 31 | 20.15 | 2.0 / 0 / 0.0 | 8.7 / 2 / 4.3 | 0 / 1.8 | 4.1, 3.3 |
| 32 | 22.15 | 4.0 / 2 / 0.0 | 10.0 / 11 / 4.7 | 6 / 2.6 | 4.2, 4.7 |
| 33 | 24.15 | 12.0 / 7 / 0.0 | 25.0 / 11 / 2.3 | 4 / 1.3 | 1.9, 1.5 |
| 34 | 26.15 | 4.3 / 2 / 0.0 | 20.0 / 9 / 3.0 | 6 / 2.3 | 2.3, 2.4 |
| 35 | 28.15 | 6.0 / 6 / 0.0 | 17.3 / 20 / 6.1 | 17 / 17.0 | 5.3, 4.0 |
| 36 | 30.15 | 3.7 / 4 / 0.0 | 8.0 / 24 / 6.5 | 13 / 12.1 | 5.7, 7.6 |
| 37 | 32.15 | 3.7 / 4 / 0.0 | 13.0 / 11 / 1.4 | 0 / 4.5 | 1.3, 4.8 |
| 38 | 34.15 | 5.3 / 4 / 0.0 | 18.0 / 9 / 1.9 | 2 / 2.1 | 1.5, 2.7 |
| 39 | 36.15 | 6.7 / 2 / 0.0 | 13.7 / 11 / 0.8 | 4 / 1.8 | 0.7, 1.9 |
| 40 | 38.15 | 4.7 / 0 / 0.0 | 10.3 / 6 / 1.0 | 0 / 0.9 | 1.0, 1.7 |
| 41 | 40.15 | 7.0 / 4 / 0.0 | 11.7 / 13 / 0.8 | 4 / 0.4 | 0.5, 3.2 |
| 42 | 42.15 | 4.7 / 2 / 0.0 | 11.0 / 7 / 0.8 | 7 / 0.3 | 0.6, 2.6 |
| 43 | 44.15 | 8.3 / 0 / 0.1 | 10.3 / 6 / 1.1 | 6 / 0.7 | 0.9, 1.4 |
| 44 | 46.15 | 14.0 / 4 / 0.0 | 25.0 / 6 / 1.1 | 7 / 0.3 | 0.8, 2.2 |

Over the first four blocks (medians): the baseline's unit 16 fires 1.4 at
base, 9 at onset and nothing in the pulse (reliability 0.08); the
candidate's unit 16 4.7, 69 and 85 (reliability 0.53); the October
template alone 3.1, 72 and 88 (0.53). So the baseline keeps a few onset
spikes of the cell and drops its train, exactly as on October 16, and the
candidate gives the same seed unit the whole train, matching the control
block for block (125 against 111, 95 against 95, 75 against 81, 71 against
73 per s). The medians over all 24 blocks (baseline 4.5, 3, 0; candidate
12.3, 11, 4) hide this because the cell is gone for 20 of them; the role
match in `tagged_blocks_*.json` also understates it, since unit 16's
template is nearer OPTO2 unit 23 than unit 22.

From block 25 (17:11 on October 20) the onset response and the train shrink
together in every column, and from block 27 (21:11) nothing on TT6 answers
the light beyond 10 per s, apart from a weak return at blocks 35 and 36
(28 to 30 h). The October template's spikes keep a median cosine of 0.70 to
0.73 and a trough of 7 to 9 MAD throughout, so what it collects later is
other activity of the same shape on the same wires, and the candidate's
unit 16 base rate rises from 2 to 9 per s to 10 to 25 as the capture fills
it with that activity. *Conjecture:* the cell moved out of reach of the
tetrode, or stopped responding, between 15:00 and 21:00 on the first day;
the first is the simpler reading since its baseline firing disappears too.
Whether the candidate's unit 16 in the first 6 h is the cell alone or the
cell plus the silenced CH49 neighbour (which the main sort's unit 23 was on
October 16) is for its purity page (section 5): the neighbour is silenced
by the light, so the in-pulse train is the cell's in either case.

### 2.2 The TT10 driven cell (OPTO2 re-sort cluster 4, template 1004)

Present at every block. Baseline seed 26 (the carry-forward's unit that
matches the driven cell's template best, cosine 0.87) and the candidate's
unit 23 (seeds 26 and 65 merged); the third column is the capture with the
October template alone:

| block | hour | baseline onset / in pulse | candidate onset / in pulse | October template onset / in pulse | baseline base | candidate base |
|---|---|---|---|---|---|---|
| 1 | 0.45 | 54 / 3.0 | 21 / 5.1 | 34 / 21.7 | 7.6 | 10.0 |
| 2 | 1.98 | 44 / 5.3 | 61 / 32.6 | 48 / 31.7 | 6.6 | 12.1 |
| 3 | 4.33 | 69 / 9.1 | 67 / 45.7 | 38 / 31.0 | 3.1 | 10.0 |
| 4 | 5.96 | 75 / 11.4 | 84 / 55.1 | 21 / 29.7 | 6.2 | 20.3 |
| 25 | 8.15 | 104 / 13.6 | 122 / 61.5 | 54 / 30.2 | 6.7 | 21.3 |
| 26 | 10.15 | 115 / 12.9 | 130 / 65.6 | 69 / 33.2 | 8.0 | 18.7 |
| 27 | 12.15 | 143 / 14.0 | 143 / 47.6 | 50 / 22.9 | 19.3 | 31.0 |
| 28 | 14.15 | 117 / 17.2 | 115 / 35.5 | 52 / 22.3 | 12.3 | 16.3 |
| 29 | 16.15 | 113 / 12.3 | 133 / 48.7 | 65 / 17.1 | 16.7 | 34.3 |
| 30 | 18.15 | 87 / 10.5 | 100 / 30.6 | 44 / 12.0 | 14.7 | 23.0 |
| 31 | 20.15 | 109 / 11.5 | 119 / 63.7 | 31 / 31.0 | 7.3 | 25.3 |
| 32 | 22.15 | 117 / 15.1 | 128 / 66.6 | 35 / 19.0 | 5.3 | 16.7 |
| 33 | 24.15 | 98 / 12.6 | 98 / 32.8 | 37 / 20.2 | 4.7 | 12.0 |
| 34 | 26.15 | 107 / 13.3 | 119 / 38.7 | 65 / 25.3 | 8.0 | 17.3 |
| 35 | 28.15 | 83 / 13.0 | 113 / 49.4 | 57 / 22.6 | 8.3 | 24.0 |
| 36 | 30.15 | 48 / 8.6 | 143 / 59.0 | 50 / 7.9 | 4.0 | 23.3 |
| 37 | 32.15 | 52 / 7.0 | 154 / 70.3 | 52 / 9.3 | 1.0 | 20.3 |
| 38 | 34.15 | 46 / 6.4 | 100 / 83.8 | 50 / 18.5 | 4.7 | 29.0 |
| 39 | 36.15 | 72 / 8.0 | 144 / 78.8 | 13 / 10.1 | 3.3 | 30.7 |
| 40 | 38.15 | 67 / 4.6 | 178 / 105.7 | 93 / 41.6 | 0.7 | 18.0 |
| 41 | 40.15 | 22 / 4.3 | 107 / 90.3 | 59 / 22.1 | 4.0 | 34.3 |
| 42 | 42.15 | 48 / 4.4 | 126 / 117.1 | 65 / 26.9 | 3.0 | 24.0 |
| 43 | 44.15 | 46 / 6.9 | 120 / 101.5 | 31 / 12.6 | 2.7 | 20.3 |
| 44 | 46.15 | 28 / 5.8 | 120 / 90.4 | 48 / 20.6 | 5.0 | 49.7 |

Medians over the 24 blocks: baseline 5.8 at base, 73 at onset, 9.8 in the
pulse, 0.3 after, reliability 0.53; candidate 20.8, 119, 60, 2.0, 0.70;
October template alone 3.3, 50, 22, 0.7, 0.37. The baseline keeps the onset
spikes and a tenth of the train, as the plan expected; the candidate keeps
the train, and its unit also carries three to four times the baseline
firing (the TT10 hash the capture assigns to this template, as on October
16, where this tetrode was already a multi-unit one). The candidate's
in-pulse rate climbs over the second day (30 to 60 per s on day one, 80 to
117 from 34 h) while the baseline's falls; *conjecture:* the unit's share
of the tetrode's driven activity grows as the carried template drifts,
which the purity pages for unit 23 should show as a widening amplitude
spread.

### 2.3 The TT4 cells (OPTO2 units 3 and 5)

Two TT4 units respond at every block the sorting has them. Baseline seed 49
(matches OPTO2 unit 3, cosine 0.92, mutual best) is absent from the
sorting for the first four blocks: it was born by re-seeding at window 12
(6 h). Candidate unit 43 is the same seed.

| block | hour | baseline onset / in pulse | candidate onset / in pulse | October template onset / in pulse | baseline base | candidate base |
|---|---|---|---|---|---|---|
| 1 | 0.45 | absent | absent | 96 / 75.0 | | |
| 2 | 1.98 | absent | absent | 115 / 72.4 | | |
| 3 | 4.33 | absent | absent | 119 / 88.5 | | |
| 4 | 5.96 | absent | absent | 92 / 67.8 | | |
| 25 | 8.15 | 83 / 11.7 | 119 / 18.3 | 180 / 118.2 | 6.0 | 6.3 |
| 26 | 10.15 | 272 / 32.0 | absent | 187 / 110.3 | 13.7 | |
| 27 | 12.15 | 230 / 39.1 | 131 / 32.7 | 163 / 104.6 | 32.0 | 16.3 |
| 28 | 14.15 | 244 / 24.6 | 122 / 36.3 | 144 / 112.5 | 12.3 | 9.7 |
| 29 | 16.15 | 235 / 12.1 | 144 / 42.7 | 128 / 98.1 | 2.3 | 6.7 |
| 30 | 18.15 | 220 / 15.6 | 122 / 46.2 | 124 / 82.5 | 4.0 | 6.0 |
| 31 | 20.15 | 239 / 15.5 | 146 / 64.3 | 111 / 90.4 | 2.7 | 8.0 |
| 32 | 22.15 | 281 / 10.5 | 176 / 76.6 | 139 / 87.8 | 3.7 | 10.7 |
| 33 | 24.15 | 200 / 9.4 | 131 / 59.5 | 78 / 52.8 | 6.7 | 12.7 |
| 34 | 26.15 | 196 / 6.5 | 137 / 69.7 | 104 / 92.0 | 3.7 | 13.0 |
| 35 | 28.15 | 185 / 3.5 | 131 / 66.7 | 119 / 90.2 | 3.3 | 11.7 |
| 36 | 30.15 | 237 / 7.7 | 176 / 89.4 | 96 / 49.8 | 2.0 | 11.7 |
| 37 | 32.15 | 289 / 7.0 | 165 / 90.0 | 120 / 62.6 | 1.3 | 11.7 |
| 38 | 34.15 | 224 / 8.3 | 141 / 74.6 | 76 / 70.8 | 5.0 | 14.0 |
| 39 | 36.15 | 239 / 4.5 | 174 / 79.1 | 102 / 93.8 | 3.3 | 13.3 |
| 40 | 38.15 | 200 / 4.0 | 120 / 68.1 | 80 / 95.9 | 0.7 | 8.0 |
| 41 | 40.15 | 185 / 3.4 | 113 / 46.9 | 81 / 65.2 | 3.0 | 9.7 |
| 42 | 42.15 | 130 / 2.4 | 135 / 54.7 | 87 / 88.8 | 2.7 | 9.7 |
| 43 | 44.15 | 146 / 2.2 | 150 / 51.7 | 100 / 101.7 | 1.0 | 7.3 |
| 44 | 46.15 | 170 / 1.4 | 159 / 49.6 | 115 / 79.2 | 2.3 | 14.3 |

Medians: baseline 2.8 at base, 200 at onset, 6.8 in the pulse, reliability
0.95; candidate 9.7, 131, 51, 0.82; October template 6.8, 113, 89, 0.69.
The baseline's onset volley is the largest of the three (the sorter keeps
the full-size first spike of every pulse) and its in-pulse rate decays from
39 per s at 12 h to 1.4 at 46 h while the October template's stays at 50 to
120; the candidate's in-pulse rate is 50 to 90 from 20 h on. The candidate
loses the unit at block 26 (10.15 h; the window had fewer than 100
baseline spikes for the template or the template was under the floor).

The second TT4 responder, baseline seed 7 (OPTO2 unit 5, cosine 0.88), is
present at every block: baseline 69 at onset and 16 in the pulse at the
median block, candidate 101 and 100, October template 72 and 42; the
candidate's in-pulse rate exceeds the October template's by two to three
times in the first day. *Conjecture:* on TT4 the capture's two templates
divide the driven activity between them differently from the October
templates, and the sum is what is stable; the purity pages decide whether
unit 7 and unit 43 are two cells.

### 2.4 The TT11 cell (OPTO2 unit 39)

Baseline seed 28 (template trough only -39 uV, cosine 0.55 to the October
template; the match by template is weak but the response identifies it)
and candidate unit 25 (seed 28):

| block | hour | baseline onset / in pulse | candidate onset / in pulse | October template onset / in pulse | baseline base | candidate base |
|---|---|---|---|---|---|---|
| 1 | 0.45 | 157 / 26.7 | 88 / 61.6 | 71 / 48.8 | 2.4 | 6.2 |
| 2 | 1.98 | 59 / 10.6 | 109 / 56.9 | 36 / 46.4 | 0.7 | 4.8 |
| 3 | 4.33 | 71 / 2.4 | 126 / 45.3 | 54 / 45.7 | 0.3 | 4.1 |
| 4 | 5.96 | 86 / 2.2 | 144 / 47.0 | 75 / 56.7 | 1.4 | 7.2 |
| 25 | 8.15 | 41 / 1.4 | 94 / 35.6 | 37 / 56.0 | 1.3 | 7.3 |
| 26 | 10.15 | 39 / 1.8 | 106 / 35.7 | 63 / 51.1 | 3.0 | 10.7 |
| 27 | 12.15 | 56 / 1.2 | 124 / 37.0 | 67 / 40.5 | 4.7 | 23.0 |
| 28 | 14.15 | 30 / 1.7 | 89 / 37.8 | 81 / 53.4 | 1.3 | 12.0 |
| 29 | 16.15 | 26 / 2.4 | 89 / 39.3 | 63 / 43.6 | 0.7 | 12.3 |
| 30 | 18.15 | 31 / 2.4 | 104 / 37.0 | 93 / 41.2 | 1.0 | 12.7 |
| 31 | 20.15 | 37 / 1.4 | 100 / 32.3 | 50 / 55.3 | 0.7 | 4.7 |
| 32 | 22.15 | 46 / 1.2 | 83 / 22.8 | 61 / 49.3 | 1.0 | 10.0 |
| 33 | 24.15 | 54 / 3.0 | 91 / 18.5 | 61 / 47.0 | 1.7 | 6.3 |
| 34 | 26.15 | 41 / 3.9 | absent | 61 / 46.4 | 2.0 | |
| 35 | 28.15 | 52 / 5.1 | absent | 54 / 53.6 | 2.3 | |
| 36 | 30.15 | 41 / 0.5 | 69 / 14.0 | 28 / 61.7 | 0.7 | 6.0 |
| 37 | 32.15 | 24 / 0.8 | 20 / 14.0 | 48 / 58.3 | 2.0 | 4.0 |
| 38 | 34.15 | 20 / 1.1 | 26 / 16.1 | 48 / 45.5 | 0.0 | 1.0 |
| 39 | 36.15 | 13 / 0.6 | 26 / 9.6 | 39 / 50.6 | 2.7 | 4.7 |
| 40 | 38.15 | 28 / 0.6 | 19 / 8.5 | 33 / 44.8 | 3.3 | 3.3 |
| 41 | 40.15 | 13 / 2.2 | 11 / 8.3 | 56 / 42.2 | 1.3 | 1.3 |
| 42 | 42.15 | 26 / 1.5 | 20 / 6.2 | 41 / 48.4 | 1.0 | 2.7 |
| 43 | 44.15 | 31 / 1.7 | 24 / 5.9 | 39 / 41.7 | 2.7 | 2.3 |
| 44 | 46.15 | 20 / 1.4 | 17 / 6.3 | 89 / 45.7 | 3.3 | 2.7 |

Medians: baseline 1.4 at base, 38 at onset, 1.7 in the pulse, reliability
0.30; candidate 4.7, 86, 21, 0.57; October template 3.8, 55, 48, 0.37. The
October template holds 40 to 60 per s in the pulse at every block; the
candidate's unit has the train for the first day and loses most of it
after 30 h (the carried template was skipped in 30 windows as shallower
than 40 uV, section 3), the baseline's never has it.

### 2.5 The cross-seeded candidates

`mp_candidate.py --xseed` offers the four October templates on their
tetrodes in every window on top of the bank. Two variants were run.

Carried and deduplicated (`mp_candidate_xseed`, the rule as specified): an
extra template is dropped in a window where the bank holds one within
cosine 0.95, and otherwise re-estimated from the spikes it claimed in the
previous window when there were 100 or more. The TT6 extra was
deduplicated in 15 windows, the first 6 h among them, so unit 16 keeps the
cell's train there (the last column of the table in 2.1) and the extra
captured nothing; once the cell faded the extra was carried on its own
claimed spikes and drifted into a hash collector (trough 51 uV, cosine
0.80 to the original, an onset volley of 33 to 69 per s at every later
block and no train). The TT10 extra (re-estimated in all 96 windows) holds
29 at onset and 17 in the pulse at the median block against 7 at base,
alongside the bank's unit 23 (119 and 60), so the two divide the driven
activity; the TT4 extra was merged with baseline seed 72 as a CCG duplicate
(cosine 0.975) and the merged unit holds 121 at onset and 76 in the pulse
at every block; the TT11 extra was merged with baseline seeds 29 and 30 on
cosine alone (0.973, verdict distinct) and holds 81 and 54. The merge rule
that accepts a cosine of 0.95 over a "distinct" verdict is the June
fallback and is worth revisiting here.

Fixed and offered everywhere (`mp_candidate_xfixed`): the October templates
as given, in every window. (Numbers filled in when the run finishes.)

What the two variants show is that cross-session templates work as seeds
when a cell has none in the bank (TT10, TT4 and TT11 extras carry their
cells' trains at every block, on tetrodes where the bank's own units also
had them), that the dedup rule correctly yields to a bank unit that holds
the cell, and that re-estimating a template from its own claimed spikes
after the cell has gone is how a template turns into a hash collector; a
carried template should be re-estimated only while its claimed spikes
keep the cell's signature, or kept fixed.

## 3. Population, presence and identity

Per unit across the 24 blocks, single units only, the ratio of the median
in-pulse rate to the median baseline rate:

| | baseline carry (88) | baseline production (84) | candidate (70) |
|---|---|---|---|
| suppressed (ratio under 0.5) | 22 | 22 | 12 |
| activated (ratio over 2) | 3 | 3 | 9 |
| unchanged | 18 | 19 | 17 |
| no baseline spikes at the median block | 45 | 40 | 32 |
| present with 10 or more spikes at every block | 32 | 32 | 27 |
| median presence fraction over blocks | 0.46 | 0.58 | 0.65 |
| suppressed units: silence onset (ms, median) | 7.5 | 7.5 | 7.5 |
| in-pulse fraction of baseline (median) | 0.022 | 0.026 | 0.087 |
| silence after offset (ms, median) | 128 | 130 | 121 |
| rebound (median) | 2.5 | 2.5 | 1.75 |
| per-unit IQR of silence onset over blocks (median) | 5.0 | 5.0 | 2.5 |

The suppression profile is stable across the 48 h for the units that show
it: the per-unit interquartile range of the silence latency over 24 blocks
is one to two bins of 2.5 ms, and the in-pulse fraction's IQR is 0.075
(baseline) and 0.18 (candidate). As on October 16 (findings section 8e),
the capture lifts small units' in-pulse rates with hash and halves the
count classed suppressed; the baseline classes 22 of 43 responsive units
suppressed, the candidate 12 of 38.

Presence over the 96 windows (20 or more spikes in a window), counted from
each unit's birth window:

| | baseline carry | candidate |
|---|---|---|
| median presence fraction from birth | 0.97 | 0.92 |
| units present in 90 % or more of their windows | 56 of 88 | 45 of 88 |
| units with more / fewer / equal windows present | | 0 / 36 / 52 |
| spikes per unit, candidate over baseline (median, quartiles) | | 1.73 (1.04 to 3.35) |

The candidate never adds a window and loses windows on the quiet tetrodes
and on TT7 and TT9, where the amplitude floor bites: the baseline tracks
units with template troughs of 28 to 38 uV (TT9 unit 23 at 8.7 Hz, skipped
in 94 of 96 windows; TT7 units 19 and 18 at 14 and 11 Hz, skipped in 82
and 43), which the 40 uV floor excludes, so TT7's median presence falls
from 1.0 to 0.55 and TT9's from 0.99 to 0.40. On the busy tetrodes (TT4 to
TT6, TT8, TT10 to TT12) presence is the same in both. The floor was set
from the October sessions' noise-trough problem; on this recording it needs
to be lower or expressed in noise units. Seven baseline units with 12 to
660 spikes in 48 h were never captured (no window reached 100 spikes).

The merge: the baseline merged 4 pairs (88 to 84), the candidate 12 (81 to
70), of which 6 are CCG duplicates on TT7, TT10 and TT11 (four TT7 units
17 to 20 collapse into one) and 6 are segregated pairs merged on cosine
alone (0.96 to 0.99) on TT1, TT2, TT13 and TT15. *Conjecture:* the TT7
collapse is the capture giving four small templates the same hash, which
the arbiter then reads as one cell; the purity pages for that unit will
show whether it is one.

## 4. The June scoreboard

Script 97's three axes, the same peaks (5.5 MAD, locally exclusive) for
both. Axis A is the fraction of detected events within 0.5 ms of a sorted
spike of the same tetrode, by event amplitude.

| | baseline carry (88) | baseline single (84) | baseline with multi-unit (100) | candidate single (70) | candidate with multi-unit (86) |
|---|---|---|---|---|---|
| A: 5.5 to 7 MAD claimed | 35.3 % | 35.3 % | 35.3 % | 61.2 % | 93.1 % |
| 7 to 9 | 67.0 | 67.0 | 84.5 | 77.9 | 96.0 |
| 9 to 10 | 81.3 | 81.3 | 90.4 | 86.2 | 97.4 |
| 10 to 12 | 84.0 | 85.1 | 92.1 | 87.0 | 98.0 |
| 12 to 16 | 85.2 | 85.8 | 91.9 | 86.5 | 98.2 |
| 16 to 24 | 70.3 | 71.3 | 77.5 | 73.3 | 97.9 |
| over 24 | 30.6 | 32.4 | 52.4 | 43.2 | 99.7 |
| overall | 52.3 | 52.5 | 59.0 | 69.9 | 94.8 |
| B: units flagged cross-contaminated | 4 of 82 | 4 of 79 | | 0 of 70 | |
| B: spike-weighted purity, full window | 0.953 | 0.953 | | 0.910 | |
| B: spike-weighted purity, tight window | 0.825 | 0.826 | | 0.804 | |
| B: median best-match fraction | 0.966 | 0.967 | | 0.926 | |
| C: same-tetrode pairs, CCG duplicate / distinct | 248: 13 / 119 | 227: 14 / 115 | | 164: 24 / 31 | |

The candidate's single units claim more of the small events (61 against 35
% at 5.5 to 7 MAD) and of the largest (43 against 31 % above 24 MAD), and
its multi-unit units take nearly everything else; the baseline leaves
48 % of events above 24 MAD unclaimed even with its bucket. On axis B the
baseline's units are more self-consistent by the June cosine measure (0.95
against 0.91 spike-weighted), the candidate having admitted spikes whose
best template is a neighbour's in 7 to 15 % of cases on the busy tetrodes;
no candidate unit is flagged because the CCG and held-out votes do not
join the cosine vote. On axis C the candidate has fewer distinct pairs (31
of 164 against 119 of 248): most of its same-tetrode pairs are ambiguous,
the hash they share filling the centre of the cross-correlogram.

## 5. Purity per unit

`sorting_purity.py` scores every spike of both sortings (trough depth in
excised-noise units and cosine to the unit's window template) and gives
per unit the `ttm.capture` purity numbers from the whole 48 h train;
single units only (`results/purity_TTM-003-48H_{candidate,baseline}.csv`):

| | candidate (70) | baseline production (84) |
|---|---|---|
| refractory contamination (Llobet, 1 ms): median | 0.000 | 1.000 |
| units under 0.05 / over 0.5 | 69 / 0 | 6 / 57 |
| fraction of ISIs under 1 ms: median (quartiles) | 0 | 0.95 % (0.34 to 2.3 %) |
| amplitude-split zero-lag ratio: median | 0.105 | 0.312 |
| units under 0.3 (one cell) / over 0.7 (two) | 45 / 19 | 39 / 28 |
| 2 to 3 ms burst ratio: median | 2.1 | 2.0 |
| neighbour verdicts: duplicate / ambiguous / distinct / segregated | 28 / 183 / 145 / 42 | 77 / 110 / 242 / 109 |
| units with no duplicate and no ambiguous neighbour | 11 | 8 |
| units with at least one duplicate neighbour | 16 | 51 |

The candidate's refractory period is enforced by its assignment rule, so
its contamination says nothing; the baseline's one per cent of sub-ms
intervals is enough to saturate the Llobet formula on 48 h trains, which
is the June concern 1 in numbers (the candidate's refractory conflict
ratio, the replacement measure, is not computed over 96 windows yet).
The amplitude-split test, which neither rule touches, favours the
candidate: 45 of 70 units look like one cell by it against 39 of 84, and
its median ratio is a third of the baseline's. The neighbour verdicts cut
the other way: the baseline has more distinct pairs and the candidate more
ambiguous ones, the hash the capture gives to small templates filling the
centre of their cross-correlograms (findings section 8e), and 16 candidate
units still have a duplicate neighbour after the merge.

The tagged units of the candidate by these numbers: TT6 unit 16
amplitude-split 0.040, neighbours 3 distinct and 2 ambiguous, no
duplicate; TT10 unit 23 split 0.105, contamination 0.025 (the merge of
seeds 26 and 65 let a few conflicts through), all 9 neighbours ambiguous;
TT4 unit 7 split 0.093 (5 ambiguous, 3 distinct) and unit 43 split 0.086
with a 2 to 3 ms burst ratio of 4.0 (5 ambiguous, 1 distinct); TT11 unit
25 split 0.052, burst ratio 4.0, 2 duplicate neighbours (units 29 and
58). The baseline's counterparts: TT6 unit 16 contamination 0.58 and
split 0.31; TT4 unit 49 contamination 1.0 and split 0.79; TT4 unit 7
contamination 1.0 and split 0.38.

Pairing the candidate's units to the baseline's by template cosine on the
same tetrode (`purity_match_*.csv`): 88 pairs, median spike agreement
0.65, 32 pairs under 0.5, 18 baseline units without a candidate partner
(the candidate merged 12 pairs and the floor dropped the low-amplitude
units of section 3). The pages, `figures/compare/purity_{candidate,baseline}.pdf`,
hold 20 random units of each plus every unit of a disagreeing pair:
waveform overlays from the unit's busiest window on its template,
amplitude against cosine, amplitude over the 48 h, the ISI histogram, the
amplitude-split cross-correlogram and the neighbour cross-correlograms
with verdicts. They are the evidence the user asked for on concern 1 and
are not summarised further here.

## 6. Wall-clock, memory and disk

| stage | baseline | candidate |
|---|---|---|
| materialise (96 jobs) | 0.77 h, 1.33 TB | shared |
| seed sort and carry-forward (16 jobs) | 6.8 h, 9.4 GB | |
| capture (80 processes, one window of one tetrode each) | | 0.62 h, 58 GB, 29 GB on disk |
| peaks (16 jobs) | 0.07 h | shared |
| merge, residual, multi-unit bucket, scoring | 2.9 h, 37 GB | merge 0.09 h, 85 GB; scoring 1.6 h |
| analyzer (16 jobs) | 1.9 h | 2.7 h (156 M spikes) |
| total after the binary | 11.7 h | 5.1 h |

The capture is tetrode-parallel and would take 8 h in one process; the
carry-forward is sequential by construction. The candidate's analyzer and
scoring are slower in proportion to its spike count.

## 7. Reading the dashboard

Measurements, then the analyst's reading, marked.

- The baseline keeps the onset spikes and drops the in-pulse trains of
  every driven cell it has (TT10 9.8 per s, TT4 6.8, TT11 1.7 in the
  pulse at the median block against 22, 89 and 48 with the October
  templates), exactly the June amplitude gate's expected failure.
- The candidate keeps the trains for the cells its bank has: the TT6
  cell's for the 6 h it is present (85 per s in the pulse against the
  baseline's 0 and the control's 88), and the TT10, TT4 and TT11 cells'
  at every block (60, 51 and 21 per s), at the cost of a baseline rate two
  to four times the baseline's on those units and of small units lifted by
  hash. A cell the bank lacks it cannot capture; the cross-seeded variants
  show cross-session templates fill that gap when it occurs.
- The TT6 cell is present for the first 6 h of the recording and not
  after; this is a property of the recording, not of either pipeline, and
  it means the 48 h tracking question for that cell is a 6 h question
  unless it returns in a form the templates do not match.
- *Conjecture:* the remaining gaps are in the seed bank's timing, not in
  the matcher: a cell the seed sort isolates late (TT4 seed 49, born at
  6 h) is missed by both pipelines until then. A pass that proposes new
  templates from the pool (the capture's unexplained events) at the pulse
  blocks, or cross-session templates where a cell was tagged before, is the
  fix; section 2.5 tests the second.
- *Conjecture:* the amplitude floor of 40 uV should be replaced by a floor
  in noise units (about 4 MAD of the excised noise) so that the baseline's
  low-amplitude units on TT7 and TT9 keep their windows.

Open questions for the user: whether units 23 (TT10), 7 and 43 (TT4) and
25 (TT11) of the candidate pass the purity pages as single cells; whether
the 4.5 MAD detection of findings section 8e should be the 48 h setting
(fewer hash events on small units, 8 % of the TT6 train lost); the cosine
purity bar (still not encoded); and whether TT6's loss after 15:00 on
October 20 matches anything in the lab log (handling, a cable, the fibre).
