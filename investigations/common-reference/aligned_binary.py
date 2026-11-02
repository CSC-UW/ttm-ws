"""Per-window wire timing for TTM-003's 48 h session and a binary of its first 12 windows with the aligned reference.

``ttm.reference.estimate_wire_timing_window`` on each 1800 s window of the Zarr store, bandpassed 300
to 6000 Hz without a reference: fits on the samples near the median's large excursions in six blocks
of 20 s spread over the window, the per-wire median over the blocks with at least 3000 such samples,
and windows with none filled from the nearest window (16 windows at a time). Fits on all samples
gave gains up to 1.9 in quiet blocks, where the median follows smaller common signals. A single block of 60 s at the middle put window
1's block inside the 2287 to 2696 s event on seven tetrodes and gave offsets up to 119 us off. Writes
``<session>/sortings_seed42_pcafix/mp_aligned/wire_timing.npz`` (starts, shifts, gains) and
``results/wire_timing_ttm003.csv`` (window, channel, shift in us, gain), and prints per tetrode how
far the window-to-window shifts move.

``--binary`` then writes ``mp_aligned/binary_w00_11``: the first 12 windows (6 h) of the store,
bandpassed and referenced with ``preprocess_for_sorting(cmr="global_aligned")``, float32, the same
frames as ``mp_baseline/binary``.

Run: .venv/bin/python investigations/common-reference/aligned_binary.py [--binary] [--jobs 64]
"""
from __future__ import annotations

import argparse
import pathlib
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import spikeinterface as si
import spikeinterface.preprocessing as spre

from ttm.reference import estimate_wire_timing_window, fill_missing_windows
from ttm.sorting import preprocess_for_sorting

HERE = pathlib.Path(__file__).parent
SESSION = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12")
STORE = SESSION / f"{SESSION.name}.blosc-zstd.zarr"
OUT = SESSION / "sortings_seed42_pcafix" / "mp_aligned"
FS = 30000
WIN = 1800 * FS
BLOCK = 20 * FS
N_WIN_BINARY = 12


def bandpassed():
    return spre.bandpass_filter(si.read_zarr(str(STORE)), freq_min=300.0, freq_max=6000.0, dtype="float32")


def one_window(start: int):
    rec = bandpassed()
    return estimate_wire_timing_window(rec, start, min(start + WIN, rec.get_num_frames()), BLOCK, n_blocks=6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", action="store_true")
    ap.add_argument("--jobs", type=int, default=64)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "wire_timing.npz"
    rec = bandpassed()
    names = [str(c) for c in rec.channel_ids]
    groups = np.asarray(rec.get_property("group")).astype(int)
    if path.exists():
        z = np.load(path)
        timing = {k: z[k] for k in ("starts", "shifts", "gains", "filled")}
        print(f"reuse {path}")
    else:
        starts = np.arange(0, rec.get_num_frames(), WIN, dtype=np.int64)
        t0 = time.time()
        with ProcessPoolExecutor(16) as ex:
            res = list(ex.map(one_window, starts.tolist()))
        timing = fill_missing_windows({"starts": starts, "shifts": np.array([r[0] for r in res]),
                                       "gains": np.array([r[1] for r in res])})
        print(f"windows filled from a neighbour: {np.flatnonzero(timing['filled']).tolist()}")
        np.savez(path, **timing, channel_ids=np.array(names))
        print(f"{starts.size} windows in {time.time() - t0:.0f} s -> {path}")
        rows = [(w, names[c], timing["shifts"][w, c] / FS * 1e6, timing["gains"][w, c])
                for w in range(starts.size) for c in range(len(names))]
        pd.DataFrame(rows, columns=["window", "channel", "shift_us", "gain"]).round(3).to_csv(
            HERE / "results" / "wire_timing_ttm003.csv", index=False)
    sh = timing["shifts"] / FS * 1e6
    step = np.abs(np.diff(sh, axis=0))
    print(f"shift range {sh.min():.1f} to {sh.max():.1f} us; window-to-window change: median {np.median(step):.2f} us, "
          f"99th percentile {np.percentile(step, 99):.2f} us; gains {timing['gains'].min():.3f} to {timing['gains'].max():.3f}")
    print("per tetrode: largest change of a wire from window 0 over the 48 h, and within windows 0 to 11 (us)")
    for g in range(16):
        d = np.abs(sh[:, groups == g] - sh[0, groups == g])
        print(f"  TT{g + 1:<3d} {d.max():5.1f}  {d[:N_WIN_BINARY].max():5.1f}")
    if a.binary:
        out = OUT / "binary_w00_11"
        pp = preprocess_for_sorting(si.read_zarr(str(STORE)), cmr="global_aligned", wire_timing=timing)
        crop = pp.frame_slice(0, N_WIN_BINARY * WIN)
        t0 = time.time()
        crop.save(format="binary", folder=str(out), dtype="float32", n_jobs=a.jobs, chunk_duration="1s",
                  progress_bar=False, overwrite=True)
        print(f"{out} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
