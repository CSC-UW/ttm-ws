---
title: Subject TTM-002
updated: 2026-10-21
status: needs_verification
---

# TTM-002

Skeleton. `experiments/TTM-002-SD.toml` carries only what is verifiable from
the archive: the session `2026-09-24_09-03-22` at
`/Volumes/tetrode/TTM2/2026-09-24_09-03-22/Record Node 102` (moved there on
2026-10-19 from the old share, where it sat one directory deeper, next to
TTM2's October sessions). The session carries three streams (IMU, the
64-channel `acquisition_board`, memory usage); `get_recording` picks the
headstage stream by name. The light schedule, deprivation window and bad
channels are commented out until confirmed; `ttm-experiment-params` will
write an empty schedule until then.

The archived `structure.oebin` originally named all 64 acquisition-board
channels with a blank string (an Acquisition Board plugin bug, documented
in `docs/developer_notes/openephys_blank_channel_names.md`); it was repaired
in place on 2026-10-21 with the original kept as `structure.oebin.orig`, and
spikeinterface now reads the session. The products currently under
`<data_root>/2026-09-24_09-03-22/` and `<data_root>/TTM-002/TTM-002-SD/` were
made outside this repository from a renamed mirror and are not inputs to
anything here; they will be superseded when the session is re-processed.
