---
title: Acquisition noise in the Open Ephys tetrode recordings
updated: 2026-10-19
---

# Acquisition noise

The ~400 Hz comb, the ~9.87 Hz transient and the period-2 sample-order
pattern seen in the first tetrode recordings, characterised in May 2026.

| File | What it is |
|---|---|
| `2026-05-30_tetrode-acquisition-noise-diagnosis.md` | The report: signatures, hypotheses tested (board LEDs, clock-locked divider), the sample-rate sweep, remaining suspects. Status: needs verification. |
| `5-19.ipynb` | First look at the 2026-05-19 test session (TTM-001/TTM-Test): loading, channel map, tetrode grouping on the SpikeInterface extractor. |
| `spectral_and_artifact.ipynb` | Per-channel PSDs, harmonic grouping and the transient's periodicity on the 2026-05-19 session. |
| `subj01-artifact-investigation.ipynb` | The same analysis on the first 60 s of the 2026-05-27 session (TTM-NOD), all 64 channels. |

Paths in the notebooks are as run on `tononi-2`; the raw sessions are under
the `archive` remote of the per-machine config.
