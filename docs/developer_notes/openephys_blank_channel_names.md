---
title: Blank channel names from the Open Ephys Acquisition Board plugin
updated: 2026-10-21
---

Every recording made on the Record Node 102 rig between 2026-06-16 and
2026-10-16 (TTM2, TTM3, and the June and September test sessions) has
`"channel_name": " "` for all headstage channels of the `acquisition_board`
stream in `structure.oebin`. Nothing else is wrong with those recordings:
`continuous.dat`, `sample_numbers.npy`, `timestamps.npy`, `sync_messages.txt`,
the events, `bit_volts` (0.195), units and sample rates are all intact and
identical in form to the TTM1 recordings, which have `CH1` to `CH64`.

neo's `OpenEphysBinaryRawIO` uses `channel_name` as the channel id, so 64
identical ids make it refuse the stream with
`signal_channels do not have unique ids`, and spikeinterface's
`read_openephys` fails with it.

## Cause

A bug in the Acquisition Board plugin (`open-ephys-plugins/acquisition-board`,
libraryVersion 2.0.0 on both rigs; code references are to `main` at
`efa4a464c9cf9d61adec73003dc104f2dc191cc1`, 2026-07-21, the v2.0.0 release).

- The naming-scheme enum has no member with value 0
  (`Source/devices/Headstage.h`, lines 28 to 32):

  ```cpp
  enum ChannelNamingScheme
  {
      GLOBAL_INDEX = 1,
      STREAM_INDEX = 2
  };
  ```

- The editor reads the setting with a default of 0 when the attribute is
  absent, and casts it to the enum unchecked (`Source/DeviceEditor.cpp`,
  lines 632 and 712; also 558 for the no-board path):

  ```cpp
  board->setNamingScheme ((ChannelNamingScheme) xml->getIntAttribute ("Channel_Naming_Scheme", 0));
  ```

- Name generation is a `switch` over the two enum members with no `default`
  (`Source/devices/oni/HeadstageONI.cpp`, lines 135 to 155). With scheme 0 it
  clears `channelNames` and adds nothing, so the GUI records every headstage
  channel with an empty name.

`settings.xml` matches this on every recording we hold: all 11 blank
recordings carry `Channel_Naming_Scheme="0"` in the Acquisition Board
editor block; all 11 named recordings carry `"1"`. It follows the rig, not the
GUI version (1.0.2 and 1.1.0 both affected) or the IMU stream (the June test
sessions have no IMU stream and are blank too).

| Rig | Recordings | Channel_Naming_Scheme | Headstage names |
|---|---|---|---|
| Record Node 101 / 121, May 2026 | 11 (TTM1, Unknown 05-26) | 1 | CH1 to CH64 |
| Record Node 102, June to October 2026 | 11 (TTM2, TTM3, Unknown 06-16, 06-24, 09-18, 09-22) | 0 | blank |

How that rig's saved configuration came to hold 0 is not verifiable from the
data. The likely route, unconfirmed, is a signal chain created or saved under
a build that did not yet write the attribute, after which the default of 0 was
carried forward in the auto-saved configuration. No human metadata step is
involved: the GUI writes these names itself from the headstage position.

There is no upstream issue or pull request about this as of 2026-10-21
(searched the acquisition-board and plugin-GUI trackers).

## Prevention

On the rig, open the Acquisition Board editor's channel list, set the
channel naming dropdown to the global index scheme, and save the signal
chain; the next `settings.xml` should carry `Channel_Naming_Scheme="1"` and
the next `structure.oebin` real names. Check the first recording afterwards.
Until that is done every new recording on the rig will have the problem.

Upstream, the fix is a `default:` case in `generateChannelNames` (or a
default of `GLOBAL_INDEX` in the three `getIntAttribute` calls, and a range
check in `setNamingScheme`). Draft issue text, should we open one:

> **Blank channel names when `Channel_Naming_Scheme` is 0**
>
> `DeviceEditor::loadVisualizerEditorParameters` reads
> `Channel_Naming_Scheme` with `getIntAttribute(..., 0)` and casts it to
> `ChannelNamingScheme`, whose members are `GLOBAL_INDEX = 1` and
> `STREAM_INDEX = 2`. `HeadstageONI::generateChannelNames` switches over
> those two values with no default, so a scheme of 0 produces an empty name
> list and every headstage channel is recorded with `channel_name` `" "`.
> We have 11 recordings (GUI 1.0.2 and 1.1.0, plugin 2.0.0) whose
> `settings.xml` carries `Channel_Naming_Scheme="0"` and whose oebin names
> are all blank; recordings from a rig with `"1"` are fine. Blank duplicate
> names make neo (and therefore spikeinterface) refuse the stream. Suggest
> defaulting to `GLOBAL_INDEX` on load and adding a `default:` case, and
> possibly rejecting out-of-range values in `setNamingScheme`.

## Repair of the existing recordings (done 2026-10-21)

The GLOBAL_INDEX scheme names channel *i* (zero-based, across headstages)
`CH{i+1}`, so the repair is deterministic. Each affected `structure.oebin`
under `/Volumes/tetrode` was edited in place with a byte-level replacement of
the N occurrences of `"channel_name": " "` by `CH1` to `CHN` in file order,
leaving formatting, CRLF line endings and every other value untouched; the
original is kept beside it as `structure.oebin.orig`. Checks before writing:
exactly one `acquisition_board` stream, blank names only in that stream,
N equals `num_channels`, and the edited JSON equals the original except for
the names. Eleven files were changed (64 channels each, 80 for
Unknown/2026-06-16): TTM2 2026-09-24, 2026-10-09 (two recordings),
2026-10-10, 2026-10-11; TTM3 2026-10-16 (both sessions); Unknown 2026-06-16,
2026-06-24, 2026-09-18, 2026-09-22. The patched copy someone else made under
`TTM2/Analysis/2026-10-09_09-33-42/oe_mirror/` was left alone. TTM3
2026-10-20 (the 48 h recording) followed on 2026-10-23; its file is owned by
root, so the backup is a plain byte copy (`copyfile`, not `copy2`, which
fails setting the timestamps). The script is kept in this session's
scratchpad only; the repair is not expected to be needed again once the
rig is fixed.

After the repair spikeinterface reads TTM2 2026-09-24 as 64 channels
`CH1`..`CH64`, 30 kHz, one segment, 48.23 h, gain 0.195.

## Stream selection on the IMU rig

neo orders streams by name, so a session with the IMU stream has
`0 = IMU_port_A`, `1 = acquisition_board`, `2 = memory_usage`, whereas the
TTM1 sessions have `0 = acquisition_board`. A manifest written with the
positional `oe_stream_index = 0` therefore pointed TTM-002 at the IMU. Since
2026-10-21 `get_recording` selects the stream by its Open Ephys name
(`Experiment.oe_stream`, default `acquisition_board`, matched as a suffix of
neo's `<node>#<processor>.<stream_name>`; exactly one match required), and
the manifests no longer carry a stream index.
