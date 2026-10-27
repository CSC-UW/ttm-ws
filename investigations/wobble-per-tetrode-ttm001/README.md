---
title: wobble and circus-omp, each on all channels at once and per tetrode, on clean stretches of TTM-001
updated: 2026-10-27
status: needs_verification
---

# Joint against per-tetrode matching for both matchers on TTM-001

The wobble span test (`investigations/wobble-fullsize-ttm001/`) ran both
matchers as spikeinterface does, fitting every template to all 64 channels
of a 1 s chunk together. On TTM-003, circus-omp run on each tetrode's
channels on its own admits about 25 % more spikes in clean data, mostly
poor fits (`investigations/per-tetrode-matching/`). The question: does
wobble change as much between joint and per tetrode as circus-omp does, and
if both change, does the order of the two matchers on the measures below
change between the joint pair and the per-tetrode pair?

Scripts run on 2026-10-27 on June's binary of session `2026-05-27_09-07-52`
(`mp_long_s2000_d170000/binary`, which starts 2000 s into the recording;
times here are in the binary's frame). The bank is the one entering window
55, rebuilt from June's spike trains by `wobble-fullsize-ttm001/01_span_bank.py`
(80 units, 2 to 9 per tetrode, 45 with a template trough under 10 MAD). The
two stretches are 300 s each inside window 55 (99,000 to 100,800 s), with no
flagged time on any tetrode in `bad_periods_crop_s2000_d170000/bad_periods.csv`.
Settings are June's: circus-omp with amplitudes 0.8 and up; wobble with an
objective threshold of 0.35 times the bank's median template energy
(1,409,713 in raw units squared, the same absolute threshold for every
tetrode in the per-tetrode arm), then a cosine gate of 0.60 on each spike's
raw snippet against its template. 1 s chunks, 16 jobs. Per tetrode is
`ttm.matching.find_spikes_per_group`. Each script's printed output is the
`_output.txt` beside it, named by stretch; spike arrays are in
`/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52/sortings_seed42_pcafix/track_eval/wobble_pertet/`.
Detected events are those of June's `spike_coverage.npz` (negative peaks at
5.5 MAD per tetrode), claimed by a spike of the same tetrode within 0.5 ms.

| Script | Question | Answer |
|---|---|---|
| `check_per_group_wobble.py 99600` | Does `find_spikes_per_group` run wobble as wobble runs separately on each tetrode's channels with that tetrode's templates? | Yes, with nothing special. Per group equals the separate runs spike for spike on the synthetic case of `packages/ttm/tests/test_matching.py` (540 spikes) and on 20 s of TTM-001 (18,502 spikes), and 16 jobs equal one. Every group gets the full bank's threshold, because the same method kwargs go to each group. |
| `compare_arms.py 99600 99900` | On a clean stretch, how far does each matcher move from joint to per tetrode? | circus-omp moves and wobble does not. circus-omp finds 172,662 spikes joint and 194,570 per tetrode; the 32,693 found only per tetrode have median cosine 0.670 (spikes both arms find: 0.844), median trough 5.9 MAD (11.4) and 1.64 % intervals under 1.5 ms (0.41 %). Most are on tetrodes with small templates: TT2 goes from 3,029 to 6,045 (four units at 3.5 to 4.4 MAD), TT12 and TT14 gain 6,106 and 6,911. Wobble finds 280,469 joint and 280,484 per tetrode before the gate (62 spikes per tetrode only, 47 joint only) and 260,301 and 260,315 after it. |
| same | Does the order of the two matchers change between the joint pair and the per-tetrode pair? | On no measure here. Per tetrode, circus-omp's precision moves toward ungated wobble's and its claimed share rises by about 2 points. Median cosine 0.837 joint and 0.821 per tetrode for circus-omp against 0.801 for gated wobble; intervals under 1.5 ms 0.36 and 0.61 % against 1.82 % before the gate and 0.16 % after it; claimed events at 10 MAD and up 84.6 and 86.7 % against 95.0 %, from 5.5 to 10 MAD 21.3 and 25.0 % against 42.6 %. The circus-omp spikes that gated wobble lacks go from 13,026 to 32,104, their median cosine from 0.679 to 0.626 and trough from 9.7 to 6.4 MAD. The share of spikes with a spike on another tetrode within 0.2 ms is at most 1.4 points above chance at the arm's rate in every arm on both stretches (circus-omp's per-tetrode-only spikes: 25.1 % against 22.9 %). On white noise at the clean noise level no arm finds a spike. |
| `compare_arms.py 100200 100500` | The same on a second stretch | The same. circus-omp 152,068 joint and 168,310 per tetrode (25,270 per tetrode only, median cosine 0.678, trough 6.2 MAD); wobble 249,241 in both arms before the gate (30 spikes differ each way), 232,288 and 232,294 after it. |

Wall time for 300 s on the shared machine (a production run was going):
circus-omp 21 s joint and 34 s per tetrode, wobble 9 s joint and 7 s per
tetrode, plus 4 s for the cosine gate.

In spikeinterface 0.105.2 wobble subtracts a found spike only from the
objective of templates that share a channel with it
(`WobbleSparsity.unit_overlap`), enforces the refractory period per unit and
iterates until no peak is above threshold, with no failure counter. With
group-sparse templates, nothing found on one tetrode changes another
tetrode's objective. Reading, conjecture: a peak that a larger peak on
another tetrode within 89 samples suppresses in one iteration is still
above threshold in the next, once the larger one is subtracted, so joint and
per tetrode differ only in the order in which a tetrode's own peaks are
taken. circus-omp's joint pass keeps one candidate per vicinity across the
whole probe and stops a chunk after five passes that add no valid spike
anywhere, which is where joint and per tetrode part.

Conclusion, conjecture: the span test compared wobble that is in effect
per tetrode already against joint circus-omp. Running both per tetrode
changes only the circus-omp arm, and on these stretches it adds small,
poorly fitting spikes rather than closing the gap to wobble's coverage. A
per-tetrode re-run of the span's wobble arms would repeat the joint result.
The one new arm would be a per-tetrode circus-omp carry over the span; the
TTM-003 carry suggests its small units would shrink window by window, TT2
here being the likeliest place. Not run.
