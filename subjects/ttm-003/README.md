---
title: Subject TTM-003
updated: 2026-10-22
---

# TTM-003

16-tetrode Open Ephys subject with an optic fibre for optogenetic
stimulation. Raw acquisitions: `/Volumes/tetrode/TTM3/<session>/`.

| Experiment | Session | What it is |
|---|---|---|
| `TTM-003-OPTO1` | `2026-10-16_14-06-44` | 32 min, seven laser trains |
| `TTM-003-OPTO2` | `2026-10-16_14-44-22` | 20.6 min, five laser trains |
| `TTM-003-48H` (name provisional) | `2026-10-20_09-02-12` | 47.98 h, 44 laser trains: 20 blocks every 2 h from 8.15 h, four earlier blocks, manual tests near 8 h |

## The 48 h recording

`/Volumes/tetrode/TTM3/2026-10-20_09-02-12/Record Node 102/experiment1/recording1/`,
one experiment and one recording, 671 GB, finished the morning of
2026-10-22. `continuous.dat` is 663,292,502,656 bytes, which is 47.98 h at
64 channels, 30 kHz, int16. Software Time 1792504932633 ms is 2026-10-20
14:02:12 UTC, 09:02:12 America/Chicago. Streams: `acquisition_board`
(64 channels, 0.195 uV per bit), `memory_usage` (1 channel, 100 Hz),
`IMU_port_A` (18 channels, 100 Hz). Events: `acquisition_board/TTL` and
`MessageCenter`.

Its `structure.oebin` has all 64 headstage channel names blank, the
Acquisition Board plugin bug in
`docs/developer_notes/openephys_blank_channel_names.md`; the two October
16 sessions were repaired on 2026-10-21 and this one has not been repaired
in place yet. On 2026-10-22 it was converted (blosc-zstd) from a symlink
mirror with a repaired copy of the oebin (`docs/SETUP.md`, tononi-2); the
session directory lives on `/ssd-raid0` behind a symlink under
`data_root`. The manifest `TTM-003-48H` carries no light schedule and
no deprivation; its bad channel is CH22 (TT15, the 400 Hz comb at 80 to
260 times its flanks), from `channel_qc.py` on this recording, which also
found TT14 and CH2, bad on October 16, to be ordinary quiet channels
here, and the spike-band noise of every channel two and a half times
higher at 7.4 h, 25.8 h and from 37 h to the end than in the quiet
windows (`figures/channel_qc.png`). The planned protocol (`stimulation_log.toml`)
was pattern A (200 ms on, 1 pulse/s) at laser level 4.5 (9.8 mW), 60
pulses per train, one train every 2 h; there is no train-by-train hand
log, so the block times and measured parameters come from the TTL,
`results/stimulation_blocks_2026-10-20_09-02-12.csv` (`stim_trains.py`,
`as_planned` column). All 1,584 real pulses are pattern A (207 ms at
0.995 Hz). Blocks 1 to 4 (58 pulses each) fall at 0.45, 1.98, 4.33 and
5.96 h; TTL trains 5 to 23 between 7.69 and 8.09 h are single pulses and
runs of 3 to 5 with one 483 ms pulse and one 164 s high interval, during
which the laser was not connected to the animal (`delivered` is False for
them); trains 24 to 44 (60 pulses each) start at 8.12 h and then every
2.000 h from 8.15 h to 46.15 h, the 21 trains that match the plan. Train
24 at 17:09:28 did not reach the animal either: in `opto_blocks.py` its
evoked spike-band response is 8 uV against 50 to 80 for every delivered
block, its slow potential a fifth of block 25's, and TT6, TT8 and TT12
show no suppression; the experimenter's fifth train is train 25 at 17:11.
So 24 trains were delivered: four of 58 pulses in the first 6 h and
twenty of 60 pulses every 2 h from 8.15 h to 46.15 h. The plan for it is
`docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`.

Per-block checks before any sort (`opto_blocks.py --include-undelivered`,
`results/opto_blocks_TTM-003-48H.{csv,json}`, `figures/blocks/`,
`figures/opto_blocks_summary.png`, 2026-10-22): the light-edge amplitude
of the 300 to 6000 Hz pulse-triggered average is 3 to 20 uV in the median
channel and 7 to 27 uV in the largest, against a noise floor of 1.4 to
11 uV that follows the animal's state, the same proportions as on October
16, so no spike-band artifact was found. The evoked volley at 1 to 6 ms
is present at every delivered block (34 to 80 uV in the largest channel,
peaking on TT4) and the slow evoked potential as well. By 5 MAD crossings
the in-pulse rate on TT6, TT8 and TT12 falls to 0 to 0.3 of the 50 ms
baseline in the blocks with a low noise floor and only to 0.5 to 1.0 in
the noisy ones (blocks 1, 2, 27 to 30, 41, 44), where the crossings are
dominated by noise; the sorted units will give the cleaner per-block
measure. TT4 and TT11 show in-pulse rates above baseline in many blocks
(up to 7 and 2 times), consistent with the driven candidates of the
findings' section 8e. The quiet tetrodes (TT1 to TT3, TT13 to TT16) have
no spikes, and their crossings are noise.

The stimulation protocol, the experimenter's train-by-train log and what
the laser TTL actually recorded are in `optogenetic_stimulation.md`
(structured copy: `stimulation_log.toml`; joined with the TTL by
`scripts/stim_trains.py` into `results/`). The manifests carry no light
schedule: both sessions are daytime recordings of under an hour and the
colony schedule for this subject has not been recorded here.

Both sessions are converted (`ttm-convert`, blosc-zstd) and have
`experiment_params.json`. `badChannels`/`badTetrodes` in the manifests come
from the raw-trace check in `optogenetic_sessions_findings.md`, which also
holds the first analysis of the sessions (artifact check, MountainSort5
sort, tagging candidates, population suppression).

`scripts/` holds the analysis scripts for these sessions; they take paths
from the manifests and write bulk outputs (figures, sortings) to the data
volume under `<data_root>/<session>/` and `<data_root>/TTM-003/<experiment>/`,
and small summaries to `results/`. Run them from this workspace's
environment, in this order:

| Script | Writes |
|---|---|
| `stim_trains.py` | `results/stimulation_trains.csv`, `results/stimulation_pulses_<session>.csv`; `results/stimulation_blocks_<session>.csv` for a session without a hand log |
| `opto_blocks.py --subject TTM-003 --experiment <E> [--blocks 24,25]` | per-block artifact and MUA suppression from windows around each train: `results/opto_blocks_<E>.{csv,json}`, `figures/blocks/opto_block<NN>.png`, `figures/opto_blocks_summary.png` |
| `channel_qc.py ... [--n-windows 12]` | per-channel std, rail fraction, MAD, crossing rate and 400 Hz comb in 60 s windows: `results/channel_qc_<E>{,_summary}.csv`, `figures/channel_qc.png` |
| `mp_baseline.py ... [--stages ...] [--max-windows N] [--tag T]` | the June TTM-001 production matching-pursuit pipeline unchanged: `<session>/sortings_seed42_pcafix/mp_baseline/{binary,assembled_reseed_c12,spike_coverage.npz,assembled_prod,analyzer_prod.zarr,production_scoreboard.json,baseline_timing.json}` |
| `opto_explore.py --subject TTM-003 --experiment <E>` | `results/opto_explore_<E>.json`, `figures/opto_*.png`, `<session>/mua_crossings_5mad.npz` |
| `opto_onset_zoom.py ...` | `results/opto_onset_zoom_<E>.json`, `figures/opto_onset_*.png` |
| `opto_sort.py ... --variant asis [--blank-ms PRE POST]` | `<session>/sortings_opto/<variant>/{preprocessed,by_group,sorting,analyzer}`, `results/units_<E>_<variant>.csv` |
| `opto_tag.py ... --variant asis` | `results/optotag_<E>_<variant>.csv`, `..._candidates.json`, `figures/optotag_*.png` |
| `opto_population.py ...` | `results/opto_population_<E>_<variant>.json`, `figures/opto_population_<variant>.png` |
| `opto_pulse_gallery.py ... [--pulses first,middle,last\|all] [--resort THR] [--thresh-mad 4.5]` | `figures/galleries/{channels,tetrodes,units,<TT>_traces}_<variant>[_resort<THR>][_<thresh>mad].pdf`, `<session>/sortings_opto/<variant>/crossings_*.npz` |
| `opto_tt6_inpulse.py ... --unit U [--tetrode TT6]` | `results/tt_inpulse_<E>_<TT>.json`, `figures/opto_<TT>_inpulse_<variant>.png` |
| `opto_tt_resort.py ... --unit U --thresholds 4.0 3.5 [--ms5 KEY=VALUE ...] [--tag T]` | `<session>/sortings_opto/<variant>/resort_<TT>_thr<thr>[_T]/`, `results/resort_<E>_<TT>[_T].csv`, `figures/opto_<TT>_resort_thr*.png` |
| `opto_tt_amplitude_split.py ... --threshold 3.5 --cluster C` | `results/ampsplit_<E>_<TT>_thr<thr>_c<C>.json`, `figures/opto_<TT>_ampsplit_*.png` |
| `opto_capture.py ... [--tetrodes TT6 TT10\|all] [--seed TT10=resort_TT10_thr4] [--templates-from DIR] [--wobble] [--noise-check] [--analyzer]` | `<session>/sortings_opto/<out>/{sorting,templates_<TT>.npz,capture_<TT>.npz,wobble_<TT>.npz,analyzer}`, `results/capture_<E>_<out>.json`, `figures/capture/<out>_<TT>_purity.pdf`, `..._noise_vs_rate.png` |
| `opto_capture_compare.py --experiments E1 E2 --out capture` | `results/capture_xsession_<E1>_<E2>_<out>.csv` |

`results/unit_metrics_<E>_asis.csv` (template and quality metrics per
unit) was written from the analyzers by hand on 2026-10-21.

`opto_capture.py` is the driver for `ttm.capture` (shape-based capture of
low-amplitude spikes from the main sort's templates;
`docs/developer_notes/shape_capture.md`). Its output variant (default
`capture`) has the same layout as `asis`, so `opto_tag.py` and
`opto_population.py` run on it with `--variant capture` once `--analyzer`
has built the analyzer.
