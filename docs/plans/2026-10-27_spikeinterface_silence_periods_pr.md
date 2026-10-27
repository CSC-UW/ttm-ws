---
title: A spikeinterface PR to let silence_periods do what ttm.bad_periods' mask does
updated: 2026-10-27
---

# spikeinterface PR: per-channel noise fill in silence_periods

Not started. The user wants it opened eventually, not now.

`ttm.bad_periods.mask_bad_periods` (a `MaskedPeriodsRecording`) replaces a
tetrode's four channels inside its bad periods with white noise at their
clean noise level. spikeinterface 0.105.2 has most of this in
`spikeinterface.preprocessing.silence_periods(mode="noise")`: Gaussian noise
per channel at `get_noise_levels`, generated in seeded 1 s blocks so it does
not depend on chunking. It lacks three things the mask needs:

1. Channel scope. A period silences every channel; the periods array has
   `segment_index`, `start_sample_index` and `end_sample_index` only. The
   mask needs each period to name its channels (a tetrode's four).
2. Edges. Noise mode switches abruptly; only `apodization` tapers, and it
   tapers to zeros. The mask fades the signal out and the noise in over
   1 ms (raised cosine).
3. Noise level. Without `noise_levels`, the levels come from random chunks
   that can fall inside the periods being silenced. The mask estimates each
   channel's MAD on chunks that avoid that channel's periods
   (`clean_noise_levels`). Passing `noise_levels` covers this, so it may
   only need a note or a helper.

Until then, per-tetrode scope is possible with built-ins alone:
`split_by("group")`, `silence_periods` per group with that group's periods
and precomputed clean noise levels, then `aggregate_channels`. That loses
only the edge fade. Before the PR, check that this composition gives the
same masked samples as `mask_bad_periods` apart from the noise draws and
the 1 ms edges, and whether the abrupt edge changes detection.

The PR would add an optional channel field (or a list of channel ids per
period) and a taper for noise mode, with tests that a chunked read equals
a whole read. Check the current spikeinterface main for either first.
