"""A small SortingAnalyzer for the TT11 driven cell's in-pulse spikes that the two matchers treat differently.

Window 84 (42 to 42.5 h, 60 pulses) of the plain preprocessed binary (``mp_baseline/binary``),
channels of TT11. Spikes come from the per-tetrode carry (``R48_W_s45_c50_lnb_b2``, "pt") and the
joint carry with the same rules and rebuild (``R48_W_joint_masked_s45_c50_lnb_b2``, "joint"); two
spikes are the same within 10 samples. In-pulse is onset to offset; outside is more than 100 ms
from every pulse. Units, as string ids:
- ``u28_in_both``: in-pulse unit 28 spikes in both carries;
- ``u28pt_u29joint``: in-pulse unit 28 spikes per tetrode that joint gives to unit 29;
- ``u28pt_jointnone``: in-pulse unit 28 spikes per tetrode that joint leaves unfit (neither 28 nor 29);
- ``u29_in_both``: in-pulse unit 29 spikes in both carries;
- ``u28_out_both``, ``u29_out_both``: the two units' spikes outside the pulses in both carries.

Written to ``<BAKE_ROOT>/review/TT11_driven_w84.zarr``.

Run: .venv/bin/python investigations/carry-bakeoff/review_driven.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402

W = 84
G = 10
PT = bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2"
JR = bk.BAKE_ROOT / "real48" / "R48_W_joint_masked_s45_c50_lnb_b2"


def accepted(path: pathlib.Path, u: int) -> np.ndarray:
    z = np.load(path / f"spikes_w{W:02d}.npz")
    return np.sort(z["frame"][(z["reason"] == 0) & (z["unit"] == u)])


def main():
    a, b = W * bk.WIN, (W + 1) * bk.WIN
    rec = si.load(bk.BINARY)
    groups = np.asarray(rec.get_property("group")).astype(int)
    rec = rec.select_channels(list(np.asarray(rec.channel_ids)[groups == G])).frame_slice(a, b)
    on, off = bk.laser_pulses()
    p28, p29, j28, j29 = accepted(PT, 28), accepted(PT, 29), accepted(JR, 28), accepted(JR, 29)
    inp = lambda f: bk.in_intervals(f, on, off)  # noqa: E731
    out = lambda f: ~bk.in_intervals(f, on - int(0.1 * bk.FS), off + int(0.1 * bk.FS))  # noqa: E731
    trains = {
        "u28_in_both": p28[inp(p28) & near(p28, j28, 10)],
        "u28pt_u29joint": p28[inp(p28) & near(p28, j29, 10) & ~near(p28, j28, 10)],
        "u28pt_jointnone": p28[inp(p28) & ~near(p28, j28, 10) & ~near(p28, j29, 10)],
        "u29_in_both": p29[inp(p29) & near(p29, j29, 10)],
        "u28_out_both": p28[out(p28) & near(p28, j28, 10)],
        "u29_out_both": p29[out(p29) & near(p29, j29, 10)],
    }
    n = rec.get_num_frames()
    trains = {k: v[(v >= a + bk.NB) & (v < b - bk.NA)] - a for k, v in trains.items()}
    for k, v in trains.items():
        print(f"{k}: {v.size:,} spikes")
    assert all(v.max() < n for v in trains.values() if v.size)
    sorting = si.NumpySorting.from_unit_dict(trains, sampling_frequency=rec.sampling_frequency)
    sorting.set_property("group", np.full(len(trains), G))
    folder = bk.BAKE_ROOT / "review" / f"TT11_driven_w{W}.zarr"
    folder.parent.mkdir(exist_ok=True)
    an = si.create_sorting_analyzer(sorting, rec, format="zarr", folder=folder, sparse=False, overwrite=True)
    an.compute("random_spikes", max_spikes_per_unit=1000, seed=0)
    an.compute("waveforms", ms_before=1.0, ms_after=2.0, n_jobs=16)
    an.compute(["templates", "noise_levels"])
    an.compute("spike_amplitudes", n_jobs=16)
    an.compute("correlograms", window_ms=20.0, bin_ms=0.25)
    an.compute("isi_histograms", window_ms=20.0, bin_ms=0.25)
    an.compute("principal_components", n_components=3, mode="by_channel_local", n_jobs=16)
    an.compute("template_similarity")
    print(f"analyzer: {folder}")


if __name__ == "__main__":
    main()
