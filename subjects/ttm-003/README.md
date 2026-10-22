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
| (no manifest yet) | `2026-10-20_09-02-12` | 48 h, optotagging pulses every 2 h |

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
16 sessions were repaired on 2026-10-21 and this one was not. Repair it
as that note describes before `ttm-convert`. It has not been converted,
has no manifest and no `experiment_params.json`; the pulse pattern,
power and block times are not logged here and must come from the TTL and
the experimenter. The plan for it is
`docs/plans/2026-10-22_ttm003_48h_sorting_plan.md`.

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
| `stim_trains.py` | `results/stimulation_trains.csv`, `results/stimulation_pulses_<session>.csv` |
| `opto_explore.py --subject TTM-003 --experiment <E>` | `results/opto_explore_<E>.json`, `figures/opto_*.png`, `<session>/mua_crossings_5mad.npz` |
| `opto_onset_zoom.py ...` | `results/opto_onset_zoom_<E>.json`, `figures/opto_onset_*.png` |
| `opto_sort.py ... --variant asis [--blank-ms PRE POST]` | `<session>/sortings_opto/<variant>/{preprocessed,by_group,sorting,analyzer}`, `results/units_<E>_<variant>.csv` |
| `opto_tag.py ... --variant asis` | `results/optotag_<E>_<variant>.csv`, `..._candidates.json`, `figures/optotag_*.png` |
| `opto_population.py ...` | `results/opto_population_<E>_<variant>.json`, `figures/opto_population_<variant>.png` |
| `opto_pulse_gallery.py ... [--pulses first,middle,last\|all] [--resort THR] [--thresh-mad 4.5]` | `figures/galleries/{channels,tetrodes,units,<TT>_traces}_<variant>[_resort<THR>][_<thresh>mad].pdf`, `<session>/sortings_opto/<variant>/crossings_*.npz` |
| `opto_tt6_inpulse.py ... --unit U [--tetrode TT6]` | `results/tt_inpulse_<E>_<TT>.json`, `figures/opto_<TT>_inpulse_<variant>.png` |
| `opto_tt_resort.py ... --unit U --thresholds 4.0 3.5 [--ms5 KEY=VALUE ...] [--tag T]` | `<session>/sortings_opto/<variant>/resort_<TT>_thr<thr>[_T]/`, `results/resort_<E>_<TT>[_T].csv`, `figures/opto_<TT>_resort_thr*.png` |
| `opto_tt_amplitude_split.py ... --threshold 3.5 --cluster C` | `results/ampsplit_<E>_<TT>_thr<thr>_c<C>.json`, `figures/opto_<TT>_ampsplit_*.png` |

`results/unit_metrics_<E>_asis.csv` (template and quality metrics per
unit) was written from the analyzers by hand on 2026-10-21.
