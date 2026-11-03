"""Detected events of the aligned-reference 48 h binary, for the coverage measure of ``bad_periods_compare.py``.

``mp_baseline.py``'s peak stage (locally exclusive, negative, 5.5 MAD, radius 40 um) on
``mp_aligned/binary`` instead of ``mp_baseline/binary``. Writes ``mp_aligned/spike_coverage.npz``.

Run: BAKE_REFERENCE=aligned48 .venv/bin/python investigations/carry-bakeoff/aligned_peaks.py [--jobs 24]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np
import spikeinterface as si
from spikeinterface.core import get_noise_levels
from spikeinterface.sortingcomponents.peak_detection import detect_peaks

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

DETECT_THRESH = 5.5


def main():
    assert bk.REFERENCE == "aligned48", "run with BAKE_REFERENCE=aligned48"
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=24)
    a = ap.parse_args()
    out = bk.SORTINGS / "mp_aligned" / "spike_coverage.npz"
    rec = si.load(bk.BINARY)
    groups = np.asarray(rec.get_property("group"))
    noise = get_noise_levels(rec, return_in_uV=False)
    t0 = time.time()
    peaks = detect_peaks(rec, method="locally_exclusive",
                         method_kwargs={"peak_sign": "neg", "detect_threshold": DETECT_THRESH, "radius_um": 40.0,
                                        "noise_levels": noise},
                         job_kwargs={"n_jobs": a.jobs, "chunk_duration": "1s", "progress_bar": False})
    s = peaks["sample_index"].astype(np.int64)
    ch = peaks["channel_index"].astype(np.int64)
    amp = np.abs(peaks["amplitude"]) / noise[ch]
    np.savez(out, peak_sample=s, peak_group=groups[ch].astype(np.int64), amp_mad=amp)
    bk.log(f"{len(s):,} events ({len(s) / (rec.get_num_frames() / bk.FS):.0f}/s) in {(time.time() - t0) / 60:.1f} min "
           f"-> {out}")


if __name__ == "__main__":
    main()
