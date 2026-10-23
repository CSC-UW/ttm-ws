---
title: Single-sample timestamp glitches in an Open Ephys 48 h recording
updated: 2026-10-22
status: needs_verification
---

The TTM-003 48 h acquisition (`TTM3/2026-10-20_09-02-12`, Record Node 102,
Acquisition Board plugin 2.0.0) has three samples in the
`acquisition_board` stream's `timestamps.npy` whose time is 85,899.379 s
too large. The sample before and the sample after sit on the clock, and
`sample_numbers.npy` runs on without a break, so the data are continuous
and only the time vector is wrong at those three samples.

| sample index (from the first sample) | sample number | time of the sample before (s) | glitched value (s) |
|---|---|---|---|
| 756,781,655 | 757,632,230 | 25,254.4077 | 25,340.3070 |
| 2,367,394,391 | 2,368,244,966 | 78,941.4989 | 79,027.3982 |
| 3,978,007,127 | 3,978,857,702 | 132,628.5901 | 132,714.4894 |

The three are exactly 3 x 2^29 samples (1,610,612,736; 14.9 h) apart, and
their times are 2^29 / 10^4 s apart; the excess is close to 2^32 x 20 us.
A counter artefact in the plugin is the likely cause (conjecture; the
code was not traced). The same scan over the TTM-001 48 h recording (two
experiments, 28.0 and 20.3 h, Record Node 121), the TTM-002 48.23 h
recording of 2026-09-24 on Record Node 102 and TTM-003 OPTO1 found no
glitch, so it is not a property of every long recording on that rig.

The scan: for each 200 M-sample block of `timestamps.npy`, every
difference between consecutive samples that departs from one period by
more than half a period, with the matching difference of
`sample_numbers.npy`. It reads the whole file (41 GB over CIFS) in a few
minutes.

## Consequences

- spikeinterface's `resample` (0.105) reads the recording's time vector
  and refuses to run on it: "Detected 3 timestamp gap(s) ... Gap sizes
  (ms): 85899.379". `ttm-make-lfp --gap-tolerance-ms 1e8` runs the
  resampling through, treating the jumps as continuous, which is right
  for these glitches and wrong for a real gap. Check the LFP store's own
  time vector afterwards.
- `ttm.spikeinterface.get_recording` now repairs such samples on load
  (`repair_timestamp_glitches`: a forward jump of over 1.5 periods that
  the next step brings back to within one period of where it belongs is
  replaced by the mean of its neighbours), prints the count and records
  the indices in the slice table (`timestamp_glitches_repaired`). A real
  gap is not touched. A store converted before 2026-10-22 carries the
  glitches in its `times_seg0` array; the fix is to write the repaired
  values into that array in place, which needs no re-conversion. The
  TTM-003 48 h store was patched that way on 2026-10-23 (samples
  756781655, 2367394391 and 3978007127 set to the mean of their
  neighbours; the chunk rewrite is atomic in zarr 2's directory store,
  so the sorting jobs reading the store at the time were not disturbed).
- Anything that maps sample indices to times through the store's time
  vector (`sample_index_to_time`, the viz bundle's `spikes.parquet`)
  gives those three samples a time 23.9 h ahead. The sorters themselves
  do not use the time vector.
