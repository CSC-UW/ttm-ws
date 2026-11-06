---
title: Subject TTM-002
updated: 2026-11-06
---

# TTM-002

`experiments/TTM-002-48H.toml` (named `TTM-002-SD` until 2026-11-03) describes the session `2026-09-24_09-03-22` at
`/Volumes/tetrode/TTM2/2026-09-24_09-03-22/Record Node 102` (moved there on
2026-10-19 from the old share, where it sat one directory deeper, next to
TTM2's October sessions). The session carries three streams (IMU, the
64-channel `acquisition_board`, memory usage); `get_recording` picks the
headstage stream by name. The light schedule (09:00 on, 21:00 off), the
novel-objects window (day two, 09:00 to 14:00) and the bad channel CH22 were
confirmed by Graham Findlay on 2026-10-21.

The archived `structure.oebin` originally named all 64 acquisition-board
channels with a blank string (an Acquisition Board plugin bug, documented
in `docs/developer_notes/openephys_blank_channel_names.md`); it was repaired
in place on 2026-10-21 with the original kept as `structure.oebin.orig`, and
spikeinterface now reads the session. The products under
`<data_root>/TTM-002/TTM-002-48H/` (`delta.idelta.zarr` and
`experiment_params.json` on 2026-11-06) were made outside this repository
from a renamed mirror and are not inputs to anything here; they will be
superseded when the session is re-processed. The session has not been
converted to a store here.
