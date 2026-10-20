---
title: Derived-product baselines
updated: 2026-10-20
---

`scripts/snapshot_derived.py` summarises one experiment's derived products
(`experiment_params.json`, the production sorting's unit set and spike
counts, and the SWA and EMG Zarr stores: shape, attributes, per-tetrode
statistics, and a hash of the values) into a JSON file, and `compare`
diffs two of them. Run it before and after a change to the library, or
after re-running an entry point, to see what moved:

```bash
uv run python scripts/snapshot_derived.py snapshot \
    --subject TTM-001 --experiment TTM-NOD -o docs/baselines/<name>.json
uv run python scripts/snapshot_derived.py compare a.json b.json [--rtol 1e-5]
```

Exact comparison is the default. `--rtol` lets float statistics differ
within a relative tolerance and downgrades a differing value hash to a
note, for runs that are numerically but not bit-for-bit equivalent.

## TTM-001 / TTM-NOD, 2026-10-20

`TTM-001_TTM-NOD_pre-rerun.json` is the state after the restructure
(Phases 1 to 4) with the products as they stood: `experiment_params.json`
rewritten by `ttm-experiment-params` on 2026-10-19 from the manifest,
and the sorting, `delta.idelta.zarr`, `stft_bandpowers.zarr` and
`synthetic_emg_methods.zarr` as produced by the pre-restructure code in
May and June 2026 and only moved since. The sorting is the registry's
`tetrode_analyses.mp_production_sorting` (87 units, 71 SUA + 16 MUA,
91,707,740 spikes). No snapshot was taken before Phase 2, so the
pre-restructure state is represented by these untouched products.

`TTM-001_TTM-NOD_post-rerun.json` is the state after `ttm-swa all` was
re-run on TTM-001 from this workspace (10 minutes; numpy 2.4.6,
scipy 1.18.1, dask 2026.8.0, xarray 2026.7.0, zarr 2.18.6).

Result of `compare pre post`:

- `experiment_params.json`, the sorting, `stft_bandpowers.zarr` and
  `synthetic_emg_methods.zarr`: identical, including value hashes.
- `delta.idelta.zarr`: same shape, time axis, attributes, min and max;
  per-tetrode means and standard deviations agree to about 1e-12
  relative. Elementwise, 26,301 of 347,668,880 values (0.008 %) differ,
  with a maximum relative difference of 7.2e-6 (median 1.1e-7) and a
  maximum absolute difference of 0.016 µV² on values up to 1.4e7. This is
  float32 rounding of the dask bandpass/Hilbert overlap computation, not a
  change in the method; `compare --rtol 1e-5` passes.

The pre-rerun stores are kept at
`/nvme/neuropixels/tmp/ttm_baseline_2026-10-20/` on `tononi-2` for an
elementwise comparison.

Not re-run: `ttm-convert`, `ttm-make-lfp`, `ttm-emg --full`, the sorting
pipeline. Those products are the May/June outputs, moved unchanged.
