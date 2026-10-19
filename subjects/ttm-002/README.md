---
title: Subject TTM-002
updated: 2026-10-19
status: needs_verification
---

# TTM-002

Skeleton. `experiments/TTM-002-SD.toml` carries only what is verifiable from
the archive: the session `2026-09-24_09-03-22`, whose Open Ephys directory is
nested one level deeper than TTM-001's
(`/Volumes/neuropixel_archive/tetrode_data/2026-09-24_09-03-22/2026-09-24_09-03-22/Record Node 102`,
the `neuropixel_archive` remote; TTM2's October sessions are on
`/Volumes/tetrode/TTM2/`), with
the 64-channel acquisition-board stream at index 0 (streams 1 and 2 are
memory usage and an IMU). The light schedule, deprivation window and bad
channels are commented out until confirmed; `ttm-experiment-params` will
write an empty schedule until then.

Known blocker for re-processing: the archived `structure.oebin` names all 64
acquisition-board channels with a single blank string, and neo refuses to
parse the session (`signal_channels do not have unique ids`). `get_recording`
needs a channel-name repair (for example a patched copy of the oebin) before
`ttm-convert` can run. The products currently under
`<data_root>/2026-09-24_09-03-22/` and `<data_root>/TTM-002/TTM-002-SD/` were
made outside this repository from a renamed mirror and are not inputs to
anything here; they will be superseded when the session is re-processed.
