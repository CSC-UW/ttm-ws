"""What the global and the cross-tetrode median reference cost, per 60 s of TTM-003's 48 h store.

One process, the whole lazy chain (Zarr read and decompression, bandpass 300 to 6000 Hz, reference)
on 60 s of all 64 channels, three repeats at different times, for: bandpass only; global median
(``common_reference(reference="global", operator="median")``); cross-tetrode median
(``ttm.sorting.preprocess_for_sorting(cmr="cross_tetrode")``, spikeinterface's ``out_of_group``
reference). The 48 h materialize with the global reference took 0.765 h with 96 jobs
(``mp_baseline/baseline_timing.json``).

Run: .venv/bin/python investigations/common-reference/cmr_cost.py
"""
from __future__ import annotations

import time

import numpy as np
import spikeinterface as si
import spikeinterface.preprocessing as spre

from ttm.sorting import preprocess_for_sorting

STORE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/2026-10-20_09-02-12.blosc-zstd.zarr"
FS = 30000


def main():
    raw = si.read_zarr(STORE)
    chains = {
        "bandpass only": spre.bandpass_filter(raw, freq_min=300.0, freq_max=6000.0, dtype="float32"),
        "global median": preprocess_for_sorting(raw, cmr="global"),
        "cross-tetrode median": preprocess_for_sorting(raw, cmr="cross_tetrode"),
    }
    starts = [int(h * 3600 * FS) for h in (2, 20, 40)]
    res = {k: [] for k in chains}
    for s in starts:
        for k, rec in chains.items():
            t0 = time.perf_counter()
            rec.get_traces(start_frame=s, end_frame=s + 60 * FS)
            res[k].append(time.perf_counter() - t0)
    for k, v in res.items():
        print(f"{k}: {np.median(v):.2f} s per 60 s (repeats {', '.join(f'{x:.2f}' for x in v)})")
    ref = {k: np.median(v) - np.median(res["bandpass only"]) for k, v in res.items()}
    print(f"reference step alone: global {ref['global median']:.2f} s, cross-tetrode {ref['cross-tetrode median']:.2f} s "
          f"per 60 s; whole chain ratio {np.median(res['cross-tetrode median']) / np.median(res['global median']):.2f}")


if __name__ == "__main__":
    main()
