"""A small SortingAnalyzer for looking at the few-unit tetrodes' units by eye.

Window 48 (24 to 24.5 h) of the plain preprocessed binary (``mp_baseline/binary``, the store both
sortings were matched on; the per-tetrode carry matched on a whitened view of it), channels of TT1,
TT6 and TT16. Units, as string ids:
- ``ptN``: unit N's accepted spikes in the per-tetrode whitened carry (``R48_W_s45_c50_lnb_b2``);
- ``ptN_removed``: unit N's fits that the 1.5 ms refractory rule removed (the lower-snr fit of each
  pair), which the export leaves out;
- ``prodN``: production's unit N (``mp_baseline_masked/assembled_reseed_c12``) in the same window.
TT1 unit 0 and TT16 unit 44 are the units in question (production's TT16 unit is 43); TT6 unit 16
is a well-isolated reference. The analyzer is dense over the twelve channels.

Written to ``<BAKE_ROOT>/review/TT1_TT6_TT16_w48.zarr``.

Run: .venv/bin/python investigations/carry-bakeoff/review_analyzer.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

W = 48
TETRODES = (0, 5, 15)
PT_UNITS = {0: 0, 16: 5, 44: 15}
PROD_UNITS = {0: 0, 16: 5, 43: 15}


def main():
    a, b = W * bk.WIN, (W + 1) * bk.WIN
    rec = si.load(bk.BINARY)
    groups = np.asarray(rec.get_property("group")).astype(int)
    rec = rec.select_channels(list(np.asarray(rec.channel_ids)[np.isin(groups, TETRODES)])).frame_slice(a, b)
    sp = np.load(bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2" / f"spikes_w{W:02d}.npz")
    prod = si.load(bk.SORTINGS / "mp_baseline_masked" / "assembled_reseed_c12")
    trains, unit_groups = {}, {}
    for u, g in PT_UNITS.items():
        for reason, suffix in ((0, ""), (3, "_removed")):
            m = (sp["unit"] == u) & (sp["reason"] == reason)
            trains[f"pt{u}{suffix}"] = np.sort(sp["frame"][m]) - a
            unit_groups[f"pt{u}{suffix}"] = g
    for u, g in PROD_UNITS.items():
        trains[f"prod{u}"] = np.sort(prod.get_unit_spike_train(u, start_frame=a, end_frame=b)) - a
        unit_groups[f"prod{u}"] = g
    n = rec.get_num_frames()
    trains = {k: v[(v >= bk.NB) & (v < n - bk.NA)] for k, v in trains.items()}
    for k, v in trains.items():
        print(f"{k}: {v.size:,} spikes")
    sorting = si.NumpySorting.from_unit_dict(trains, sampling_frequency=rec.sampling_frequency)
    sorting.set_property("group", np.array([unit_groups[u] for u in sorting.unit_ids]))
    out = bk.BAKE_ROOT / "review" / f"TT1_TT6_TT16_w{W}.zarr"
    out.parent.mkdir(exist_ok=True)
    an = si.create_sorting_analyzer(sorting, rec, format="zarr", folder=out, sparse=False, overwrite=True)
    an.compute("random_spikes", max_spikes_per_unit=1000, seed=0)
    an.compute("waveforms", ms_before=1.0, ms_after=2.0, n_jobs=16)
    an.compute(["templates", "noise_levels"])
    an.compute("spike_amplitudes", n_jobs=16)
    an.compute("correlograms", window_ms=20.0, bin_ms=0.25)
    an.compute("isi_histograms", window_ms=20.0, bin_ms=0.25)
    an.compute("principal_components", n_components=3, mode="by_channel_local", n_jobs=16)
    an.compute("template_similarity")
    print(f"analyzer: {out}")


if __name__ == "__main__":
    main()
