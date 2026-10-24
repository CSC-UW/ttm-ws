---
title: TTM-003 48 h, the staged pipeline
updated: 2026-10-23
status: needs_verification
---

# TTM-003 48 h: the staged pipeline (track F)

Track F of `docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`, run on
2026-10-23: the June pipeline unchanged for full-size spikes, and shape
capture only on what it leaves, only for units that answer the light. The
settings were measured on the two October 16 sessions and on 48 h blocks 1
to 4 and 25 to 27, then the pipeline ran on the whole 48 h recording
(`TTM-003-48H`, 47.98 h, 24 delivered trains of 207 ms pulses at 0.995 Hz
and 9.8 mW; trains 5 to 24 reached no tissue and count as no light). It is
judged as track D was, against the June baseline and the shape-capture
candidate of `docs/reports/2026-10-23_ttm003_48h_baseline_vs_candidate.md`.
Measurements are stated as such; sentences marked *conjecture* are the
analyst's reading and have not been confirmed by the user.

## 1. Terms

- June stages, or stages 1 and 2: the June TTM-001 production pipeline as
  `mp_baseline.py` ran it on this recording (circus-omp carry-forward at
  amplitudes of 0.8 of the template or more, with per-window template
  re-estimation and re-seeding every 12 windows of 1800 s; the CCG-guarded
  merge; the residual pass that offers unclaimed events of 10 MAD or more to
  their best template at cosine 0.8). Their output, `assembled_prod_sua`
  (84 single units), is read, not re-run. "June spikes" of a unit are its
  spikes there.
- Stage 3: `ttm.staged_capture` (`docs/developer_notes/staged_capture.md`).
  Per window and tetrode the June spikes are fitted with their unit's window
  template and subtracted; what is left is the residual. Its troughs below
  the floor are scored against every template of the tetrode.
- Floor: the detection threshold of stage 3, in units of the spike-excised
  noise of each wire (the median absolute deviation of the samples farther
  than 1 ms from any crossing, `ttm.capture.noise_mad_excised`; "excised
  MAD"). On the 48 h windows that noise is 10 to 16 uV, so 3.5 excised MAD is
  35 to 56 uV.
- Raw cosine: the four-wire cosine between an event's -1 to +2 ms snippet
  and a template, amplitude free. Noise-corrected cosine: the same with the
  expected noise energy removed from the snippet's norm, so a small spike of
  exactly the template's shape scores near one instead of being pulled down
  by the noise in its window (`ttm.capture.noise_corrected_cosine`).
- Offered: an event whose best template (by raw cosine) reaches the gate on
  the noise-corrected cosine. Every template competes, so an event that
  looks most like a unit that is not admitted is not given to one that is.
- Admitted: a unit stage 3 runs for. A unit is admitted at a train when its
  June spikes plus the events it was offered around that train (its "trial"
  stage-3 spikes, kept as if it were admitted) pass the light test; the
  decision holds until the next delivered train, and before the first train
  the first train decides.
- Light test: per pulse, spikes from 10 ms after onset to offset ("in
  pulse", 197 ms here) against spikes in a window of the same length ending
  at onset ("pre-pulse"). Inside a train the pre-pulse window lies between
  pulses, 600 ms after the previous offset, and is not a baseline. The test
  passes when the in-pulse rate is at least twice the pre-pulse rate and a
  one-sided exact test on the two counts gives p below 0.001.
- Left: an event whose best template is the unit's but whose
  noise-corrected cosine is under the gate.
- Sorting-free control: the traces around each block captured with the
  October 16 (OPTO2) templates themselves at 3.5 MAD and raw cosine 0.6
  (`tagged_cells_xsession.py`), no 48 h sorting involved; on the October 16
  sessions the same capture of the whole session.
- Refractory conflict ratio: the events a unit was offered that its
  refractory rule removed (within 1.5 ms of one of its June spikes or of a
  better-matching new event), over the number two independent trains of
  those sizes would give. Near zero: the new spikes respect the unit's
  refractory period; near one: they come from another cell as often as
  chance allows.
- Amplitude-split ratio: coincidences within 1.5 ms between a unit's largest
  and smallest spikes over the independent expectation (one cell gives near
  zero, two cells near one); neighbour verdicts are the June arbiter's
  cross-correlogram verdicts against every same-tetrode unit (duplicate,
  ambiguous, distinct, segregated).
- Block: a delivered train, numbered by its TTL train (1 to 4, 25 to 44).

## 2. The pipeline as built

`subjects/ttm-003/scripts/mp_staged.py`, writing to
`<session>/sortings_seed42_pcafix/mp_staged/`:

1. and 2. The June stages, read from `mp_baseline/assembled_prod_sua`.
3. Offer: for each 1800 s window that holds a delivered train, every
   tetrode; the June spikes subtracted (the median fitted scale per window
   and tetrode is 0.94, 0.87 to 0.96 for 80 % of them, so the fit removes
   them at about their own size), the residual's troughs at the floor, the
   offer with three peeling passes (`offer/TT<n>_w<NNN>.npz`, every event with
   its best template, both cosines and whether it was offered). Then the
   light test at every delivered train (`admission.csv`), and the offer for
   the remaining windows on tetrodes with a unit admitted there.
   Admit: offered events are kept for units admitted at the governing train,
   outside 1.5 ms of the unit's June spikes and winning their refractory
   chain (`assembled_stage3`).
4. The June CCG-guarded merge of the 84 single units (propose cosine 0.90,
   accept on a duplicate verdict or cosine 0.95), then the June multi-unit
   bucket on what remains (5.5 MAD peaks of 7 MAD or more that no unit
   claims, with best template cosine 0.55 to 0.8, one unit per tetrode), the
   analyzer and the June scoreboard (`assembled_staged`,
   `analyzer_staged.zarr`, `production_scoreboard.json`). The bucket comes
   after the merge, as in June; with stage 3 empty the pipeline reduces to
   the June assembly up to a second merge pass, which merged nothing here.

The settings are arguments with the chosen defaults: `--floor 3.5`,
`--gate-rc 0.8`, `--min-ratio 2`, `--max-p 0.001`. The light test and the
offer are in the library; the train table (delivered trains only) and the
paths are the subject's.

## 3. Settings, measured on the October 16 sessions

`staged_settings.py` ran stage 3 on the June stages of each session (OPTO1
whole, 32 min, 44 units with a template; OPTO2 whole, 20.6 min, 42 units; the
June pipeline was run on OPTO2 for this, 51 single units) at floors 3.5 and
4.5 and at noise-corrected gates 0.5 to 0.95, with the raw cosine 0.6 of the
candidate as a reference arm, and kept for every unit its offered, trial and
left events (`results/staged_settings_<E>.csv`, summary
`results/staged_settings_summary.csv`). The June unit matched to the TT6
cell is OPTO1 unit 18 and OPTO2 unit 21 (template cosine 0.99 to the OPTO2
cell). Rates below are over the pattern-A trains (207 ms pulses; in OPTO1 at
0.5, 0.55 and 4.2 mW, in OPTO2 at 4.2 mW; the response saturates at the
lowest power, findings section 4) unless marked as the 64 ms 9.8 mW
pattern-B pulses.

### 3.1 The gate on the noise-corrected cosine

The TT6 cell's June spikes fire 0.2 to 0.9 per s in the pulse. What stage 3
adds, at floor 3.5:

| session | gate | in pulse, admitted | in pulse, left | 9.8 mW 64 ms, admitted | pre-pulse, admitted | added outside trains (per s) | refractory conflict ratio | offered vs June spikes, CCG ratio |
|---|---|---|---|---|---|---|---|---|
| OPTO1 | raw 0.6 | 101 | 54 | 139 | 7.6 | 14.8 | 0.126 | 0.47 |
| OPTO1 | 0.7 | 147 | 8 | 205 | 11.7 | 21.2 | 0.124 | 0.52 |
| OPTO1 | 0.75 | 143 | 12 | 195 | 9.5 | 17.5 | 0.106 | 0.50 |
| OPTO1 | 0.8 | 136 | 19 | 183 | 7.4 | 14.0 | 0.090 | 0.48 |
| OPTO1 | 0.85 | 127 | 28 | 164 | 5.4 | 10.7 | 0.071 | 0.41 |
| OPTO1 | 0.9 | 114 | 41 | 140 | 3.7 | 7.7 | 0.054 | 0.34 |
| OPTO2 | raw 0.6 | 74 | 55 | 111 | 12.0 | 15.1 | 0.188 | 0.64 |
| OPTO2 | 0.7 | 119 | 10 | 169 | 18.2 | 22.1 | 0.174 | 0.63 |
| OPTO2 | 0.75 | 114 | 15 | 161 | 14.2 | 17.8 | 0.144 | 0.61 |
| OPTO2 | 0.8 | 107 | 22 | 150 | 10.8 | 13.8 | 0.122 | 0.59 |
| OPTO2 | 0.85 | 97 | 32 | 133 | 7.9 | 10.3 | 0.089 | 0.42 |
| OPTO2 | 0.9 | 84 | 46 | 111 | 4.7 | 7.3 | 0.074 | 0.30 |

The sorting-free control gives the cell 103 and 87 per s in the pulse
(pattern A) and 144 and 131 (pattern B); the June unit's own spikes outside
the trains fire at 14.2 and 12.3 per s. The last column is the
cross-correlogram of the events the unit was offered against its own June
spikes outside the trains, 0.3 to 1.5 ms over 5 to 25 ms; the June unit's
autocorrelogram on the same measure is 0.37 and 0.55.

What the pages show (`figures/staged_settings/TTM-003-OPTO{1,2}_full_TT6_cell_f{3.5,4.5}.png`,
one column per gate: PSTH of admitted, left and June spikes around the
pulses, overlays of admitted events in and outside the pulses and of left
events, trough histograms; `summary_TTM-003-OPTO{1,2}.png` plots the table
against the gate for all four tagged units):

- The events left under the raw cosine 0.6 have the cell's in-pulse time
  course in the PSTH (54 and 55 per s in the pulse, decaying like the
  admitted ones) and the cell's shape in the overlays: the raw gate drops a
  third (OPTO1) to two fifths (OPTO2) of the train because a half-size spike
  in a 3 ms window of noise has a low raw cosine, which is what the
  noise-corrected cosine corrects.
- At every noise-corrected gate the admitted and left events in the pulse
  sum to the same 155 (OPTO1) and 129 (OPTO2) per s; the gate only divides
  them. In the pulse the tetrode's other units are silenced, so this pool
  is the cell's spikes plus at most the pre-pulse rate of foreign events.
- At 0.85 and above the left events form an onset peak in the first 10 to
  20 ms of the pulse: the cell's spikes that overlap the population volley
  at onset fail a strict gate.
- Outside the trains the gate trades the cell's small spikes against foreign
  ones: from 0.7 to 0.9 the rate added there halves or more and the
  refractory conflict ratio falls from 0.12 to 0.05 (OPTO1).

The choice: 0.8, the lowest gate at which the TT6 cell's stage-3 spikes are
at least as clean as with the candidate's raw 0.6 by the two measures that
do not depend on the fit (refractory conflict ratio 0.090 against 0.126 and
0.122 against 0.188; rate added outside the trains 14.0 against 14.8 and
13.8 against 15.1), in both sessions, while keeping 35 % (OPTO1) and 45 %
(OPTO2) more of the in-pulse train than the raw gate. The neighbours, for
the user to overrule: 0.75 keeps 5 to 7 % more of the train and adds 25 to
29 % more outside the trains; 0.85 loses 7 to 9 % of the train (mostly at
onset) and adds 24 to 25 % less outside them. The cross-correlogram column
is not used to choose: the fit residue inside the 3 ms scoring window
depletes offered events within 1.5 ms of any June spike, other units' as
much as the unit's own (for the TT6 cell at 0.8 in OPTO1, at lags of 0.75
to 1.5 ms, 0.51 of the flank against its own June spikes and 0.47 against
the other TT6 units'), so it understates foreign events at every gate
alike.

The other tagged units' June clusters are not single units by this measure
(autocorrelogram ratios 1.3 to 1.7 for the TT10 and TT4 units, 0.8 to 0.9
for TT11), so they cannot judge purity; at 0.8 they gain in the pulse over
the raw gate (TT10 53 against 41 and 24 against 23; TT4 43 against 23 and
61 against 33; TT11 85 against 76 and 104 against 100 per s).

### 3.2 The detection floor, 3.5 against 4.5 excised MAD

At gate 0.8, the TT6 cell's in-pulse train at 4.5 against 3.5: 100 against
136 per s (OPTO1) and 80 against 107 (OPTO2) over the pattern-A pulses, 153
against 183 and 124 against 150 over pattern B, so 4.5 loses 25 to 27 % of
the 207 ms train and 16 to 17 % of the 64 ms one. Its refractory conflict ratio is 0.086 against 0.090 and
0.115 against 0.122, and the rate added outside the trains 10.6 against
14.0 and 10.6 against 13.8. The TT4 and TT11 units lose more: 9 against 43
and 13 against 61 per s in the pulse on TT4, 40 against 85 and 59 against 104
on TT11. Findings section 9 showed, window by window, that the class only
3.5 admits shares the cell's shape, refractory period and pulse response on
October 16 and in 48 h blocks 1 and 4. The choice: 3.5, because 4.5 removes a
quarter or more of every driven unit's in-pulse train while the TT6
cell's refractory conflict ratio moves by under 0.01. In the staged pipeline the cost 4.5 was meant
to avoid (findings 8e: hash lifting the in-pulse rates of small units that
the light silences) does not arise, because stage 3 does not run for units
that fail the light test.

### 3.3 No template floor

None is used. A floor in noise units would drop admitted units that answer
the light: by the trough of the window template over the noisiest wire's
excised noise, the admitted units under 5 are, on October 16, the TT4 driven
unit (OPTO1 unit 3 at 4.2, OPTO2 unit 5 at 4.4), and in 48 h blocks 1 to 4
and 25 to 27, TT11 units 29 (3.3 to 4.9) and 32 (4.7 to 4.8), TT12 unit 38
(4.2 to 4.7), TT9 unit 23 at block 3 (3.0) and TT4 unit 49 at block 26 (2.4).
A floor at the detection floor itself (3.5) would drop the last two and TT11
unit 29 in some windows; a floor at 5 would also drop the TT4 cell of
October 16. The units the candidate's 40 uV floor protected against (quiet
tetrodes, TT14 on October 16) never pass the light test.

### 3.4 How a template earns and keeps its place

The test and its threshold, from the per-train measurements of every unit
(`results/staged_settings_rule.csv`, `figures/staged_settings/rule.png`):
in-pulse rate at least twice the pre-pulse rate, one-sided count test below
0.001, on the June spikes plus the trial stage-3 spikes of that train. On
October 16 the untagged units sit at 0.4 to 1.0 times (quartiles; 90 % of
them under 1.7), the tagged cells at 1.4 to 9.7 (the TT6 cell 4.4 to 9.7,
the TT11 cell 2.5 to 7.4, the TT4 cell 1.9 to 3.8, the TT10 driven unit 1.4
to 3.7). The admitted sets:

- OPTO1 (7 trains, 0.5 to 9.8 mW): the TT6 cell and the TT11 cell at all 7,
  the TT4 cell at 6 (1.94 at the 0.5 mW train), the TT10 driven unit at 5
  (1.83 and 1.98 at trains 1 and 6), TT5 units 9 and 10 at all 7 (the weakly
  activated TT5 units of findings 8e, 21 to 41 per s in the pulse), and TT10
  unit 28 and TT11 units 33 to 35 at one to three trains each; 48 of 308
  unit-trains.
- OPTO2 (5 trains): the TT6, TT4 and TT11 cells and TT5 units 11 and 13 at
  all 5, TT5 unit 14 at 4, the TT10 driven unit at 2 and TT10 unit 30 at the
  other 3, TT11 unit 36 at one; 35 of 210.
- No large unit that the light silences is admitted in either session.

On the 48 h windows, the TT6 cell's unit (June unit 16) at blocks 1, 2, 3,
4, 25, 26, 27: in pulse over pre-pulse 25, 21, 22, 11, 4.4, 1.8 and 0.5
(its June and trial spikes together fire 196, 158, 129, 106, 42, 25 and 8
per s in the pulse against 6 to 17 before it), so it is admitted early,
admitted weakly at block 25 and dropped
from block 26, which matches the traces (trains 1 to 4 have the cell, train
25 weakly, later trains not). Its June spikes alone fire 0 to 3.8 per s in
the pulse at every block, which is why the test needs the trial spikes. A
threshold of 1.5 would keep it at block 26; one of 3 would drop the TT10
driven unit at most October 16 trains.

Two TT6 units pass at trains where the cell was present although their own
June spikes are silenced in the pulse (0 per s at every block): unit 15 at
blocks 2, 3, 4, 25 and 26 (2.2 to 3.6 times, 25 to 46 per s of trial spikes
in the pulse) and unit 14 at blocks 3 and 4 (2.6 and 2.1). *Conjecture:*
their trial spikes are the TT6 cell's attenuated spikes, or after block 25
what remains of it, landing on templates that drift closer to them (unit
15's trough grows from -78 to -132 uV between blocks 1 and 25); the test
cannot tell this from a light-driven cell of their own, and the purity pages
of section 4.3 are where the user can judge.

### 3.5 Confirmed on 48 h blocks 1 to 4

`mp_staged.py --max-windows 12` (the first 6 h) admitted 12, 15, 19 and 16
units at blocks 1 to 4, the same sets as the settings run, in 5.3 min for
offer and admit, and gave the TT6 cell 195, 159, 129 and 107 per s in the
pulse against the control's 111, 95, 81 and 73. The full run reproduces
these windows exactly (same offers, same admissions), and its numbers are
the ones reported below.

## 4. The 48 h run

`mp_staged.py` with the defaults, all stages, on the baseline's binary and
June stages (`<session>/sortings_seed42_pcafix/mp_staged/`). Of the 84
June single units, 75 have a template in the window of at least one
delivered train; 28 were admitted at one train or more, 282 of the 1140
unit-trains tested (`results/staged_admission_TTM-003-48H.csv`,
`figures/staged_admission.png`). Per train 5 to 19 units are admitted; fewest
at the trains in noisy or active stretches (block 27 at 21:11, after lights
off: 7; blocks 28, 34 and 44: 5 or 6), where the pre-pulse rates of most
units are raised. Stage 3 added 34.5 M spikes to those 28 units (the June
stages have 85.3 M in all); 0.14 M of them (0.4 %) fall inside a pulse and
the rest between pulses and between trains, in the windows the units were
admitted for. The merge found no pair to merge (84 units stay 84), and the
multi-unit bucket took 5.50 M events (the baseline's 6.35 M, less what
stage 3 now claims).

Because 99.6 % of the added spikes lie outside the pulses, a second variant
was run from the same offers with stage 3 kept inside the delivered pulses
only (`--pulses-only --offer-from ""`, `mp_staged_pulses/`, offer reused,
admit, merge and bucket again): 0.14 M stage-3 spikes, every unit identical
to its June self outside the pulses, the same admissions. Both are reported.

### 4.1 The tagged cells block by block

`tagged_cells_blocks.py` on the staged single units
(`results/tagged_blocks_TTM-003-48H_staged.csv`; the pulses-only variant
`..._staged_pulses.csv`), next to the baseline, the candidate and the
sorting-free control of the earlier report. Rates per second; pre-pulse
here is the table's -50 to 0 ms window.

The TT6 cell (June unit 16, the same seed unit as in the baseline and the
candidate):

| block | hour | staged pre-pulse / onset / in pulse | June in pulse | candidate pre-pulse / in pulse | control onset / in pulse | admitted |
|---|---|---|---|---|---|---|
| 1 | 0.45 | 4.1 / 98 / 195.4 | 3.7 | 2.1 / 124.9 | 59 / 110.9 | yes |
| 2 | 1.98 | 6.9 / 88 / 158.6 | 0.1 | 5.9 / 94.7 | 84 / 94.8 | yes |
| 3 | 4.33 | 6.2 / 58 / 128.9 | 0.0 | 3.4 / 74.8 | 98 / 81.2 | yes |
| 4 | 5.96 | 12.1 / 46 / 106.6 | 0.0 | 9.3 / 70.5 | 50 / 73.1 | yes |
| 25 | 8.15 | 9.3 / 35 / 42.3 | 0.0 | 11.3 / 19.8 | 35 / 26.0 | yes |
| 26 | 10.15 | 5.0 / 2 / 0.0 | 0.0 | 16.3 / 10.6 | 17 / 13.7 | no |
| 27 | 12.15 | 4.7 / 0 / 0.0 | 0.0 | 17.7 / 4.0 | 9 / 5.2 | no |
| 28 to 44 (median) | | 4.7 / 2 / 0.0 | 0.0 | 13.7 / 2.3 | 6 / 2.1 | no |

Over blocks 1 to 4 the staged unit fires a median 144 per s in the pulse
against the control's 88, the candidate's 85 and the baseline's 0. Its
in-pulse spikes have a median trough of -59, -54 and -49 uV at blocks 1, 4
and 25 against -96, -138 and -142 uV for its June spikes there (0.61, 0.39
and 0.34 of their size). From block 26 the unit is the June unit unchanged:
its pre-pulse rate stays at 2 to 14 per s where the candidate's unit rises to
10 to 25 as its template collects look-alike activity (the collector of the
earlier report, section 2.1), and nothing is captured in the pulse. The weak
return the control shows at blocks 35 and 36 (17 and 12 per s in the pulse)
is not captured: the unit's light test there gives 0.5 and 1.4. The
pulses-only variant gives the same in-pulse rates (195.3, 158.3, 128.6,
106.2, 42.1) with the June pre-pulse rates (0.7, 3.1, 0.7, 2.1, 3.0).

The other tagged and driven units (medians over the 24 blocks; staged unit =
June unit):

| unit | staged pre-pulse / onset / in pulse | pulses-only pre-pulse / in pulse | June pre-pulse / onset / in pulse | candidate pre-pulse / onset / in pulse | control onset / in pulse | trains admitted |
|---|---|---|---|---|---|---|
| TT10 driven, unit 26 | 15.0 / 91 / 50.0 | 5.8 / 49.9 | 5.8 / 73 / 9.8 | 20.8 / 119 / 60.3 | 50 / 22.2 | 18 of 24 |
| TT4, unit 49 (born at 6 h) | 14.8 / 194 / 68.8 | 2.8 / 68.7 | 2.8 / 180 / 6.8 | 9.7 / 131 / 50.6 | 113 / 88.7 | 18 of 20 present |
| TT4, unit 7 | 11.5 / 107 / 142.8 | 3.5 / 142.8 | 3.5 / 69 / 16.4 | 5.8 / 101 / 99.8 | 72 / 41.6 | 24 |
| TT11, unit 28 | 6.7 / 46 / 27.9 | 1.4 / 27.9 | 1.4 / 38 / 1.7 | 4.7 / 86 / 20.6 | 55 / 47.7 | 21 |

The staged pipeline gives every driven unit its in-pulse train at the trains
where it passes the light test, five to sixteen times the June rate at the
median block (and from nothing for the TT6 cell). Where the
test fails the unit falls back to June for two hours: the TT10 unit at block
1 (where June unit 27 carries the driven activity and is admitted instead),
at blocks 27 to 30 (12 to 20 h, ratios 1.3 to 1.7, in-pulse 10 to 17 per s)
and at block 44; the TT11 unit at blocks 34, 36 and 44. The candidate,
without a test, keeps these trains at every block. In the whole-window
version the pre-pulse rate of the admitted units is two to five times June's
(the small spikes added between pulses); in the pulses-only version it is
June's.

### 4.2 The units the capture should not have touched

The TT7 and TT9 small units whose windows the candidate's 40 uV floor
dropped (TT9 unit 23 at -37 uV, TT7 units 17 to 19 at -42 to -49 uV) are June
units in the staged sorting, with June's presence: `results/staged_presence_TTM-003-48H.csv`
gives the median fraction of windows with 20 spikes or more, from each
unit's first window: 1.00 on TT7 and 1.00 on TT9 in the staged sorting and
in the baseline, against 0.60 and 0.41 in the candidate. On every tetrode
the staged value equals the baseline's (the quiet ones too: TT1 0.37, TT14
0.23, TT16 0.79, against the candidate's 0.14, 0.03 and 0.44), since stage 3
only adds spikes to units that exist. Stage 3 touched them only where they
passed the light test: TT7 unit 18 at 4 trains (0.38 M stage-3 spikes, 19 %
of its June count), TT7 unit 17 at block 36 (0.10 M), TT9 unit 23 at block 3
(0.10 M); in the pulses-only variant those additions are 899, 294 and 89
spikes in the pulses.

After the TT6 cell faded, no TT6 template keeps collecting: unit 16 leaves
stage 3 at block 26 and unit 14 at block 4; unit 15 stays admitted at blocks
25, 26, 29 and 36 (7 trains in all, 1.63 M stage-3 spikes, 3144 of them in
pulses), the case section 3.4 flags.

### 4.3 Purity

Per unit (`results/purity_TTM-003-48H_staged.csv`,
`results/staged_units_TTM-003-48H.csv`) and as pages
(`figures/compare/purity_staged.pdf`: the 28 units stage 3 touched, the TT7
and TT9 units and 10 others; `purity_staged_baseline.pdf`: the same units in
the baseline).

The refractory conflict ratio of the stage-3 spikes, per admitted unit over
the windows it was admitted in: median 0.27 (quartiles 0.18 to 0.45, range
0.08 to 0.68). The cleanest are TT6 unit 15 (0.08), TT7 unit 17 and TT11
unit 33 (0.10), TT11 unit 28 (0.16) and TT6 unit 14 (0.16); the TT6 cell's
unit 16 is at 0.28 (0.09 and 0.12 on October 16); the highest are TT10
units 27 and 56 (0.68, 0.62), TT5 unit 77 (0.55), TT4 unit 7 (0.54), TT10
unit 24 (0.53) and TT11 unit 32 (0.50). The rule removes the offered events
that conflict (1.0 M of the 35.5 M offered to admitted units) but not the
foreign ones that happen not to. *Conjecture:* if the foreign events fire
independently of the unit and its own never conflict, the ratio is their
share, so about a quarter of what stage 3 adds at the median unit and a
third to two thirds on the busiest TT4, TT5, TT10 and TT11 units.

The measures of the whole unit do not separate the pipelines well, for
reasons of construction: refractory contamination (Llobet) of a staged unit
falls as stage 3 adds spikes that cannot violate (TT11 unit 28 from 0.52 to
0.03), and the amplitude-split ratio of a staged unit compares mostly June
spikes (large) with stage-3 spikes (small), which the refractory rule keeps
apart by construction (TT4 unit 49 from 0.79 to 0.09), as it does for every
candidate unit. Neither is evidence of purity here. The neighbour verdicts
are not touched by the rule: over the 84 single units, duplicate /
ambiguous / distinct / segregated are 82 / 143 / 204 / 109 for the staged
sorting, 77 / 110 / 242 / 109 for the baseline and 28 / 183 / 145 / 42 for
the candidate (70 units); stage 3 turns 38 distinct pairs into ambiguous or
duplicate ones, the hash it adds to admitted units filling the centre of
their cross-correlograms (the candidate has 97 fewer distinct pairs than
the baseline, on 70 units). In the pulses-only variant every unit is its
June self outside the pulses, so its purity numbers are the baseline's.

### 4.4 The June scoreboard

Script 97's three axes on the same 5.5 MAD peaks
(`results/mp_staged_TTM-003-48H_scoreboard.json`; baseline and candidate
from the earlier report). Axis A is the fraction of detected events within
0.5 ms of a sorted spike of the same tetrode, by event amplitude; axis B the
cosine self-consistency of each unit's spikes and the units flagged by two
of three purity signals; axis C the CCG verdicts of same-tetrode pairs.

| | baseline single (84) | baseline with multi-unit (100) | staged single (84) | staged with multi-unit (100) | candidate single (70) | candidate with multi-unit (86) |
|---|---|---|---|---|---|---|
| A: 5.5 to 7 MAD claimed | 35.3 % | 35.3 % | 42.0 % | 42.0 % | 61.2 % | 93.1 % |
| 7 to 9 | 67.0 | 84.5 | 70.5 | 85.4 | 77.9 | 96.0 |
| 10 to 12 | 85.1 | 92.1 | 86.1 | 92.3 | 87.0 | 98.0 |
| over 24 | 32.4 | 52.4 | 32.4 | 52.4 | 43.2 | 99.7 |
| overall | 52.5 | 59.0 | 57.2 | 62.9 | 69.9 | 94.8 |
| B: units flagged cross-contaminated | 4 of 79 | | 5 of 79 | | 0 of 70 | |
| B: spike-weighted purity, full window | 0.953 | | 0.914 | | 0.910 | |
| B: spike-weighted purity, tight window | 0.826 | | 0.797 | | 0.804 | |
| B: median best-match fraction | 0.967 | | 0.945 | | 0.926 | |
| C: same-tetrode pairs, CCG duplicate / distinct | 227: 14 / 115 | | 227: 15 / 97 | | 164: 24 / 31 | |

Stage 3 claims more of the small events (42 against 35 % at 5.5 to 7 MAD)
and nothing of the largest (32.4 % above 24 MAD, as June). On axis B the
staged units are as self-consistent by the June cosine measure as the
candidate's (0.914 against 0.910 spike-weighted) and less so than June's
(0.953): the stage-3 spikes of an admitted unit have its window template
as best match less often, which is the same cost the conflict ratios show.
One more unit is flagged than in June (5 against 4). On axis C, 18 of
June's 115 distinct pairs are no longer distinct (97); the candidate has
31 distinct pairs among its 164.

### 4.5 Compute

Wall-clock on tononi-2, shared with another agent's jobs for part of the
run; the staged pipeline used at most 48 offer workers and 16 SpikeInterface
jobs (`results/mp_staged_TTM-003-48H_timing.json`).

| stage | baseline | candidate | staged |
|---|---|---|---|
| materialise (96 jobs) | 0.77 h, 1.33 TB | shared | shared |
| seed sort and carry-forward (stage 1) | 6.8 h | baseline's | baseline's, not re-run |
| peaks | 0.07 h | shared | shared |
| June merge and residual pass (stage 2) | 0.32 h | | baseline's, not re-run |
| capture, or offer and admit | | 0.62 h | 0.45 + 0.05 h |
| merge | (in stage 2) | 0.09 h | 0.05 h |
| multi-unit bucket | 0.39 h | (in the capture) | 0.26 h |
| scoreboard | 2.15 h | 1.6 h | 1.16 h |
| analyzer | 1.9 h | 2.7 h | 1.89 h |
| total after the binary | 11.7 h | 5.1 h | 3.86 h on top of stages 1 and 2 and the peaks (7.2 h), 11.1 h from the binary |

Peak resident memory 27 GB (scoreboard), disk 20 GB under `mp_staged/`
(6 GB of offer files, 3.5 GB of purity scores). The offer is
tetrode-parallel: 711 window-tetrode tasks with a median of 75 s each,
19.7 h of single-process work. The pulses-only variant, reusing the offers,
took 0.42 h for admit, merge and bucket. The confirmation on the first 12
windows took 5.3 min for offer and admit.

## 5. Reading the comparison

Measurements first, then the analyst's reading, marked.

- Settings: floor 3.5 excised MAD, gate 0.8 on the noise-corrected cosine,
  light test at twice the pre-pulse rate with p below 0.001, no template
  floor. At 0.8 the TT6 cell keeps 136 and 107 per s of its in-pulse train
  on October 16 (raw 0.6: 101 and 74), with a lower refractory conflict
  ratio than the raw gate and no more spikes added outside the trains; 4.5
  would cost a quarter of every driven unit's train.
- Against the June baseline: the in-pulse trains. The TT6 cell 144 per s in
  the pulse over blocks 1 to 4 (June 0); the TT10 unit, the two TT4 units
  (49 and 7) and the TT11 unit 50, 69, 143 and 28 per s at the median block
  (June 10, 7, 16 and 2). The
  56 units that never pass the light test are June units, unchanged.
- Against the candidate: the units the capture should never have touched.
  TT7 and TT9 keep June's presence (1.00 and 1.00 against 0.60 and 0.41),
  no unit is lost to a template floor, and after the TT6 cell fades its
  unit is the June unit again (pre-pulse 2 to 14 per s, nothing in the
  pulse) instead of a collector (10 to 25 per s). The candidate keeps more
  of the trains where the light test fails (TT10 at blocks 1, 27 to 30 and
  44; TT11 at 34, 36 and 44) and catches some of the TT6 return at blocks 35
  and 36 (6 per s against the control's 12 to 17), which the staged pipeline
  does not.
- The cost, in the whole-window version: stage 3 adds 34.5 M spikes to 28
  units, 99.6 % of them between the pulses, with a refractory conflict
  ratio of 0.27 at the median unit, and 38 same-tetrode pairs that were
  distinct in the baseline become ambiguous or duplicate. On the June
  scoreboard the staged single units claim 57.2 % of the detected events
  (June 52.5, candidate 69.9) at a spike-weighted cosine purity of 0.914
  (June 0.953, candidate 0.910). The pulses-only version adds 0.14 M spikes,
  all inside the pulses, and leaves every unit its June self elsewhere; its
  tagged-cell trains are the same.
- Compute: 3.86 h on top of the June stages (offer and admit 0.49 h, merge
  0.05 h, bucket 0.26 h, analyzer 1.89 h, scoreboard 1.16 h), 11.1 h from
  the binary counting the June stages it reads, against 11.7 h for the June
  pipeline and 5.1 h for the candidate after the binary.
- *Conjecture:* the spikes stage 3 adds between the pulses are partly the
  admitted cells' own small spikes (the burst-attenuated spikes the
  amplitude-split test found on October 16) and partly other cells' spikes
  of similar shape at 3.5 to 7 MAD, in a proportion that differs by tetrode
  (conflict ratio 0.08 on TT6 unit 15, 0.5 to 0.7 on the busiest TT4 and
  TT10 units). Inside the pulses the tetrode's other cells are silenced, so
  the same gate admits few foreign events there. If the deliverable must keep
  every unit offered as single at least as clean as June's, the pulses-only
  version does that and still carries the trains; the whole-window version
  is the one to use if the cells' small spikes between pulses matter more
  than that cost. That choice is the user's.
- *Conjecture:* the remaining gaps are in the June bank, not in stage 3:
  the TT4 cell of October 16 has no unit before 6 h (re-seeded at window
  12), and the TT10 driven activity moves between June units 26 and 27, so a
  train at which the wrong one passes the test leaves the other's share out.

## 6. Open questions for the user

- Whole-window or pulses-only stage 3 for the deliverable (section 5).
- Whether TT6 units 14 and 15, admitted while their own full-size spikes
  are silenced in the pulse, should be kept out of stage 3 (their pages are
  in `figures/compare/purity_staged.pdf`, pages 16 and 17).
- Whether a unit that fails the light test at one train in a noisy or
  active stretch (block 27, after lights off) should keep its place until it
  fails twice; the plan's rule, used here, drops it for two hours.
- The light test's threshold of 2: at 1.5 the TT6 cell would stay admitted
  at block 26, at 3 the TT10 driven unit would drop at most October 16
  trains (section 3.4).
