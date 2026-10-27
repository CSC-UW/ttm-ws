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
| `carry_compare.py joint=mp_baseline pertet=mp_baseline_pertet_w14 joint_masked=mp_baseline_masked_w14 pertet_masked=mp_baseline_masked_pertet_w14` | Over the first 14 windows of the carry-forward (`mp_baseline.py --max-windows 14`, with and without `--per-tetrode-matching` and `--bad-periods`), what does per-tetrode matching give in clean time? | Many more spikes, mostly bad. TT5 to TT12's single units have 35.8 M clean spikes per tetrode against 7.4 M joint, 59 % of them with an interval under 1.5 ms against 2.4 %, median cosine 0.29 against 0.72 (masked: 35.4 M against 7.6 M, 58.7 % against 2.5 %). The few-unit tetrodes have 29.5 M against 1.8 M unmasked and 55.1 M against 1.0 M masked. Tetrodes whose units are large change little (TT6: 0.87 % short intervals against 0.76 %). Masked, the claimed share on the few-unit tetrodes rises from 48.9 to 95.3 %, because their units fire at hundreds of hertz. At those rates chance alone puts a spike of another tetrode within 0.2 ms of about 70 % (unmasked) to 80 % (masked) of instants, so the share of spikes coincident with another tetrode (61 to 86 % per tetrode, 12 to 23 % joint) says little here. |
| `carry_windows.py <same runs> --tetrodes 1 6 9 13 14 16` | Window by window, how does that happen? | The small units' templates run away. Per tetrode, TT9 has 16,100 clean spikes in window 0 (joint 5,652) and 2.6 M masked, 3.0 M unmasked in window 6 (joint 15,643 and 18,766); their median trough falls from 4.6 to 3.0 MAD and their median cosine from 0.52 to about 0.21. TT13 and TT14 do the same in both per-tetrode runs, TT1 and TT16 in the masked one. Unmasked, TT1's and TT16's units take 275,000 to 407,000 spikes in window 1, most of them in flagged time, and have at most 18 clean spikes a window from window 2 to 12. TT6 looks the same in all four runs. |

The carries took 2.55 h per tetrode masked and 1.88 h unmasked, against
1.30 h for the joint masked carry, but the two per-tetrode carries ran at the
same time as each other and as the masked 48 h run's production stage, so
these are not a speed comparison. `bad_periods_compare.py --max-windows 14`
against the baseline wrote `subjects/ttm-003/results/bad_periods_*_w14_pertet.*`
and `*_w14_masked_pertet.*`; its per-unit medians (`rv_frac` 0.024 per
tetrode against 0.010) hide the runaway, which sits in a few units.

Conclusion, conjecture: per-tetrode matching as it stands cannot go into the
carry-forward. On a tetrode with few or small units it admits fits to noise
that the joint matcher keeps out, the carry rebuilds each template from its
window's spikes, and the smaller template admits more noise in the next
window. The amplitude floor is relative to the template, so it follows the
template down. In clean data the joint matcher's competition across
tetrodes holds this back: a signal seen on several tetrodes is fitted on one
of them, and a weak candidate must beat the best fit anywhere on the probe
within 3 ms. Masking removes the cross-tetrode loss without giving that up.
Not tried: per-tetrode matching with an amplitude floor in MAD rather than
relative to the template, or with templates re-estimated only from
full-size spikes.
