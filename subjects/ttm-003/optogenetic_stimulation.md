---
title: TTM-003 optogenetic stimulation sessions (2026-10-16)
updated: 2026-10-21
---

# TTM-003 optogenetic stimulation, 2026-10-16

Two short Open Ephys sessions on `/Volumes/tetrode/TTM3/`, recorded on
Record Node 102 with the IMU rig (streams `IMU_port_A`, `acquisition_board`,
`memory_usage`; GUI 1.1.0, Acquisition Board plugin 2.0.0). Both
`structure.oebin` files carried the blank-channel-name bug and were repaired
on 2026-10-21 (`docs/developer_notes/openephys_blank_channel_names.md`).
The laser TTL is line 1 of the acquisition-board TTL events; line 2 holds
one sub-millisecond toggle per session that coincides with a line-1 edge
(crosstalk, ignored).

| Experiment | Session | Duration | Trains |
|---|---|---|---|
| `TTM-003-OPTO1` | `2026-10-16_14-06-44` | 32.0 min | 5 x pattern A, 2 x pattern B |
| `TTM-003-OPTO2` | `2026-10-16_14-44-22` | 20.6 min | 3 x pattern B, 2 x pattern A |

## Hand log (Graham Findlay, recorded 2026-10-21)

Verbatim from the experimenter. BC = behaviour change, NR = neural
response. The structured copy that scripts read is `stimulation_log.toml`.

```
# Recording 2026-10-16_14-06-44

## Pulse pattern A: 200ms laser on, 1/s, for 120s.

1. Laser level 4.0 (0.5 mW intensity), delivered at ~1m30s. Subject awake. No obvious behavior change (BC) or neural response (NR).
2. Laser level 4.1 (0.55 mW), ~5m. Awake. No BC, maybe NR.
3. Laser level 4.2 (4.2 mW), ~9m. Awake. No BC, clear NR.
4. Laser level 4.2 (4.2 mW), ~18m. Drowsy. No BC, clear NR.
5. Laser level 4.2 (4.2 mw), ~21m. Drowsy. No BC, clear NR.

## Pulse pattern B: 50ms laser on, 3/s, for 120s.

6. Laser level 4.5 (9.8 mW), ~26m. Awake. No BC, clear NR.
7. Laser level 4.5 (9.8 mW), ~29m. Awake. No BC, clear NR.

# Recording 2026-10-16_14-44-22

## Pulse pattern B (see above)

1. Laser level 4.5 (9.8 mW), ~1m. Awake. No BC, clear NR.
2. Laser level 4.5 (9.8 mW), ~4m. Drowsy. No BC, clear NR.

## Pulse pattern A (see above)

3. Laser level 4.2 (4.2 mW), ~8m30s. Drowsy. No BC, clear NR.
4. Laser level 4.2 (4.2 mw), ~11m30s. Drowsy. No BC, clear NR.
5. Laser level 4.5 (9.8 mW), ~15m. Drowsy. No BC, clear NR.
```

## What the TTL record says

`scripts/stim_trains.py` groups the laser TTL pulses into trains and joins
them, in order, with the hand log; the result is
`results/stimulation_trains.csv` (one row per train) and
`results/stimulation_pulses_<session>.csv` (one row per pulse, seconds from
the first acquisition-board sample). Every train in the hand log has a TTL
train at the logged time (all twelve onsets within 5 s of the estimates).
Measured pulse parameters differ from the nominal ones:

| Pattern | Nominal | Measured (TTL) |
|---|---|---|
| A | 200 ms, 1 Hz, 120 s | 207 ms pulses, 1004 ms period (0.996 Hz), 118 pulses over 118 s |
| B | 50 ms, 3 Hz, 120 s | 64 ms pulses, 349.5 ms period (2.86 Hz), 355 pulses over 124 s |

Each train is preceded (by 0.4 s for pattern A, within the first period for
pattern B) by two to five sub-millisecond toggles on the laser line; the
scripts drop them (`glitch` column) and they are not counted above.

The fifth train of `2026-10-16_14-44-22` is logged under pattern A at
laser level 4.5, but its TTL pulses are pattern B (355 x 64 ms at 2.86 Hz).
The hand log keeps the experimenter's entry; the CSV carries both
`logged_pattern` and `ttl_pattern`.
