---
title: circus-omp per tetrode against circus-omp on all channels at once, TTM-003 48 h
updated: 2026-10-27
status: needs_verification
---

# Per-tetrode matching on TTM-003

circus-omp, the June pipeline's matcher, fits every template to all 64
channels of a 1 s chunk together. In each pass it keeps the best template at
each sample over the whole probe and one candidate per 6 ms of the probe, and
it stops the chunk after five passes that add no valid spike anywhere
(`docs/developer_notes/bad_periods.md`, "Why it matters to tetrodes without
bad periods"). The question: does matching each tetrode on its own channels
(`ttm.matching.find_spikes_per_group`) remove the loss on clean tetrodes
during artifacts, and what does it change in clean data?

Scripts run on 2026-10-27 on the `mp_baseline` binary of session
`2026-10-20_09-02-12`. Templates are those the carry-forward would use: every
unit of `mp_baseline/assembled_reseed_c12` with 100 spikes or more in the
window before the stretch, built from those spikes (44 units for window 0).
Matcher settings are June's (amplitudes 0.8 and up, `max_failures` 5, 1 s
chunks, 16 jobs). Paths are literal `tononi-2` paths. Each script's printed
output is the `_output.txt` beside it, named by stretch.

| Script | Question | Answer |
|---|---|---|
| `compare_stretches.py 2400 2460 --masked` | During the 2287 to 2700 s even-channel event, does per-tetrode matching keep TT5 to TT12's spikes? | Yes. Unmasked, joint finds 4,943 TT5 to TT12 spikes and per tetrode 42,730, as many as per tetrode masked; they fall evenly through each chunk. Joint masked finds 33,372. Claimed 10 MAD events in clean time: 60.7 % joint unmasked, 96.1 % joint masked, 96.7 % per tetrode either way. |
| `compare_stretches.py 2100 2220` | What changes on a clean stretch? | Per tetrode finds 74,124 spikes against 59,222 and claims 94.6 % of the 10 MAD events against 91.8 %. Most of the increase is on tetrodes with small or few units: TT9 5,434 against 1,145, TT1 1,083 against 213, TT16 2,470 against 918. |
| `extra_spikes.py 2100 2220` | Are the added spikes good? | Mostly not. Their median cosine to their unit's template is 0.40 to 0.64, against 0.58 to 0.81 for spikes both arms find, and their intervals under 1.5 ms multiply on the small units (TT9's unit 23: 30 joint, 911 per tetrode; TT16's unit 44: 132 and 933). On white noise at the clean noise level both arms find nothing. |
| `settings_sweep.py 2100 2220` | Does a higher amplitude floor or fewer failures make per-tetrode matching as precise as joint? | No. A floor of 0.9 already lowers TT5 to TT12's claimed share from 94.7 to 90.8 % while the few-unit tetrodes keep 12 to 21 % of intervals under 1.5 ms (joint: 3.4 %); `max_failures` barely matters. |
| `snr_and_coincidence.py 2100 2220` | Are the added spikes small fits, or events shared with other tetrodes? | Both. On TT9 they are smaller (median matched-filter signal-to-noise 9.8 against 14.0); on the few-unit tetrodes they are large (median 21.8, about that of joint's spikes on TT5 to TT12) and 64 % of them have a spike on another tetrode within 0.2 ms, against 34 % of joint's. No signal-to-noise floor brings the few-unit tetrodes' short intervals near joint's. |

Conclusion so far, conjecture until the carry-forward runs below are in: the
joint matcher's competition across tetrodes, which costs clean tetrodes their
spikes during multi-tetrode artifacts, also does two things the June
pipeline's precision relies on in clean data. A signal seen on several
tetrodes at once is fitted on one of them, and a weak candidate must beat the
best fit anywhere on the probe within 3 ms. Per-tetrode matching as it stands
gives up both.

Pending: 14-window carry-forwards with `mp_baseline.py --per-tetrode-matching`,
masked (`mp_baseline_masked_pertet_w14`, against `mp_baseline_masked_w14`) and
unmasked (`mp_baseline_pertet_w14`, against the first 14 windows of
`mp_baseline`).
