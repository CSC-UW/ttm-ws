"""Synthetic EMG from a converted 30 kHz tetrode store.

Both estimators of ``ecephys.xrsig.synthetic_emg`` (backed by ``emg_from_lfp``)
are computed on the same 1500 Hz decimation of a few single channels, one per
widely separated tetrode: ``per_window`` (exact per-window mean pairwise Pearson
correlation) and ``global`` (amplitude-weighted global normalisation). Each
Open Ephys experiment is processed on its own and the results concatenated, so
the wall-clock gap between experiments survives on the session timeline.
"""

from __future__ import annotations

import pathlib
import time

import dask.array as da
import numpy as np
import pandas as pd
import spikeinterface as si
import xarray as xr
import zarr
from ecephys.xrsig import core as xrsig

DEFAULT_TETRODES = ("TT1", "TT8", "TT16")
DECIM_Q = 20  # 30 kHz -> 1500 Hz, the most decimation synthetic_emg allows
EMG_TARGET_SF = 20  # synthetic_emg default output rate


def _segment_uv(zg, gains, col_idx, start_sample, end_sample, fs):
    traces = da.from_zarr(zg["traces_seg0"])
    seg = traces[start_sample:end_sample][:, col_idx]
    seg_uv = seg.astype("float32") * gains[col_idx].astype("float32")
    return xr.DataArray(seg_uv, dims=("time", "channel"), attrs={"fs": fs, "units": "microvolts"})


def segment_emgs(seg_uv: xr.DataArray, t_start: float, t_end: float) -> xr.Dataset:
    """Decimate once, then compute both EMG methods on a session-time axis."""
    t0 = time.perf_counter()
    dec = xrsig.decimate_timeseries(seg_uv, q=DECIM_Q).load()
    dec = dec.assign_coords(time=np.linspace(t_start, t_end, dec.sizes["time"]))
    print(f"  decimated -> {dec.attrs['fs']:.0f} Hz, {dec.sizes['time']} samples "
          f"({dec.sizes['time'] / dec.attrs['fs'] / 3600:.2f} h) in {time.perf_counter() - t0:.1f}s")
    t0 = time.perf_counter()
    ds = xrsig.synthetic_emg(dec)
    print(f"  emg (both) -> {ds.sizes['time']} samples in {time.perf_counter() - t0:.1f}s")
    return ds


def synthetic_emg_from_store(
    store: str | pathlib.Path,
    slice_table: pd.DataFrame,
    *,
    tetrodes: tuple[str, ...] = DEFAULT_TETRODES,
    window: tuple[float, float] | None = None,
) -> xr.Dataset:
    """Both synthetic-EMG estimators for a whole session, or for one window.

    ``window=(start_s, duration_s)`` restricts the computation to that span of
    the first segment (for a quick validation run); ``None`` processes every
    segment of ``slice_table`` and concatenates them.
    """
    store = pathlib.Path(store)
    rec = si.read_zarr(str(store))
    fs = float(rec.sampling_frequency)
    tetrode = np.asarray(rec.get_property("tetrode"))
    gains = np.asarray(rec.get_property("gain_to_uV"))
    col_idx = np.array([int(np.flatnonzero(tetrode == tt)[0]) for tt in tetrodes])
    print(f"channels (tetrodes {tetrodes}): cols {col_idx.tolist()}")
    zg = zarr.open(str(store), mode="r")

    if window is None:
        segs = []
        for row in slice_table.itertuples():
            print(f"segment {row.oe_experiment_name}/{row.oe_recording_name}: samples "
                  f"{int(row.start_sample)}:{int(row.end_sample)}  t {row.t_start:.0f}..{row.t_end:.0f}s")
            seg = _segment_uv(zg, gains, col_idx, int(row.start_sample), int(row.end_sample), fs)
            segs.append(segment_emgs(seg, row.t_start, row.t_end))
        return xr.concat(segs, dim="time")

    start_s, dur_s = window
    row = slice_table.iloc[0]
    s0 = int(row.start_sample) + int(start_s * fs)
    s1 = s0 + int(dur_s * fs)
    print(f"WINDOW run: {row.oe_experiment_name}, {dur_s / 60:.0f} min from {start_s / 3600:.1f} h")
    seg = _segment_uv(zg, gains, col_idx, s0, s1, fs)
    return segment_emgs(seg, row.t_start + start_s, row.t_start + start_s + dur_s)
