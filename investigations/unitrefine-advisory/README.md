---
title: UnitRefine noise/neural classifier on tetrodes
updated: 2026-10-19
---

# UnitRefine advisory

Scripts 46-50, run in June 2026 against the chunk-tracked TTM-001 analyzer
(`sortings_seed42_pcafix/tracked_48h/analyzer_clustered.zarr`). The outcome
is the `tetrode_analyses.unitrefine_advisory` module and its two persisted
unit properties; its docstring states what the ranking is and is not.

| Script | Step |
|---|---|
| `46_apply_unitrefine_noise_classifier.py` | The stock 37-feature model on the tetrode analyzer |
| `47_unitrefine_reduced_feature_retrain.py` | Retrain on UnitRefine's own data without geometry features |
| `48_unitrefine_drop_duration_confounded.py` | Drop the recording-length-confounded features |
| `49_unitrefine_drop_firing_range.py` | Drop `firing_range` |
| `50_persist_unitrefine_advisory.py` | Write `unitrefine_neural_prob` and `unitrefine_label` onto the analyzer |

Paths are as run on `tononi-2`.
