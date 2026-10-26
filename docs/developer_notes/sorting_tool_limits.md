---
title: Limits of the sorting tools on tetrode data
updated: 2026-10-22
---

Facts about the sorters and SpikeInterface that shaped the pipeline, each
with the record that established it. They hold on any machine.

- MountainSort5 stores detected spike times as int32 in
  `mountainsort5/core/detect_spikes.py`, so any single recording passed to
  a scheme (a scheme 3 block, or a scheme 2 training or full recording)
  longer than 2^31 samples (19.9 h at 30 kHz) crashes. This, not memory,
  caps block and training duration. `investigations/spike-sorting/`,
  `SORTING_COMPARISON_FINDINGS.md` and script `23_`.
- MountainSort5 is bit-reproducible only with both random sources
  controlled: the whitening seed (`whitening_seed`, SpikeInterface
  0.105.0, PR 4600) and the PCA solver (deterministic in mountainsort5
  0.5.9, upstream PR 65). `SORTING_COMPARISON_FINDINGS.md`.
- MountainSort5's default detection threshold is 5.5 and is applied to
  the whitened signal; the SpikeInterface wrapper keeps it. Burst-attenuated
  and small spikes at 4 to 5 MAD are below it. `ms5_parameter_comparison.md`;
  `subjects/ttm-003/optogenetic_sessions_findings.md` sections 5a to 5c.
- `BinaryFolderRecording` loads its full float64 time vector into memory
  on every reconstruction (no `mmap_mode`), so `run_sorter_by_property` on
  a frame slice of a long recording multiplies that vector by the number
  of workers. Crop to a real recording and `reset_times()` before sorting;
  the drivers do. `investigations/si-frame-slice-memory/`.
- With a time vector set in memory, `to_dict` carries `times_kwargs` (the
  whole vector) at every nesting level; the Zarr writer dumped it into
  `.zattrs` (4.9 GB). `ttm.spikeinterface.convert_recording` strips it
  (`packages/ttm/tests/test_convert_provenance.py`); upstream #4885 strips
  the top level only.
- Kilosort4's final clustering allocates a dense (spikes by clusters)
  tensor, about 42 GiB for one tetrode's 48 h, which does not fit a 32 GB
  GPU. The `kilosort` fork pin chunks it (`clustering_chunk_size`, default
  2**18 rows; the chunk size does not change the sort).
  `investigations/spike-sorting/KS4_FEASIBILITY_FINDINGS.md`.
- Tetrode geometry in the ProbeGroup is fictional (tetrodes 300 um apart,
  wires 10 um apart), so sorters and drift models that localise spikes
  (DARTsort, DREDge, pitch-based templates) do not apply. Drift is handled
  by re-estimating templates in time windows. `TRACKING_FINDINGS.md`,
  `MATCHING_PURSUIT_FINDINGS.md`.
- The amplitude-gated matcher of the matching-pursuit production sort
  (circus-omp, `amplitudes=[0.8, inf]`) and its detection floors (events
  at 5.5 MAD; residual capture at 10, multi-unit bucket at 7) drop
  spikes below 0.8 of the template amplitude and below 5.5 MAD. The
  shape-gated wobble matcher that keeps them was built and tested but
  never adopted; the full 48 h comparison was not run.
  `MATCHING_PURSUIT_FINDINGS.md` ("Wobble as primary", "Phase B ladder"),
  `OPEN_CONCERNS.md`.
- UnitRefine's stock noise/neural classifier is near chance on tetrodes;
  the retrained advisory ranking in `ttm.unitrefine_advisory` orders units
  usefully but its probabilities are not calibrated. Sort by it, do not
  threshold it. `investigations/unitrefine-advisory/`.
- Open Ephys streams are selected by name (`oe_stream` in
  `ttm.spikeinterface.get_recording`), not by index: neo orders streams by
  name, and a rig with the IMU stream puts `IMU_port_A` first.
  `openephys_blank_channel_names.md`, "Stream selection on the IMU rig".
- The SpikeInterface `group` property is 0-based; the lab's tetrode names
  come from the 1-based `tetrode` property (TT1 to TT16). Label figures
  and tables `TT{group + 1}`; several TTM-003 figures were off by one
  before this was caught.
- "Hash" in the sorting documents means multi-unit activity: real spikes
  of several unresolved neurons. Where a document writes "hash/noise",
  the two were not separated.
