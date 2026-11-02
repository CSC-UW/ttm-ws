"""MountainSort5 re-seed sorts of the aligned-reference 48 h binary, for the A48 arms of ``run_arm.py``.

June's re-seed sorts window 12, 24, ..., 84 with ``ttm.tracking.sort_chunk`` (per tetrode, scheme 2),
on the binary (``mp_baseline/reseed_sorts``) and on the binary with each tetrode's bad periods masked
(``mp_baseline_masked/reseed_sorts``). This makes both from ``mp_aligned/binary``:
``mp_aligned/reseed_sorts/reseed_w<w>`` and ``mp_aligned_masked/reseed_sorts/reseed_w<w>``, one process
per sort. ``_bake.reseed`` reads a sort when its ``by_group`` exists, so these must finish before an
arm reaches the window. A sort whose ``by_group`` exists is skipped.

Run: BAKE_REFERENCE=aligned48 .venv/bin/python investigations/carry-bakeoff/reseed_aligned.py [--workers 14]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

WINDOWS = range(12, 96, 12)


def one(job):
    masked, w = job
    from ttm.bad_periods import mask_bad_periods, read_table
    from ttm.tracking import sort_chunk

    out = bk.SORTINGS / ("mp_aligned_masked" if masked else "mp_aligned") / "reseed_sorts" / f"reseed_w{w}"
    if (out / "by_group").exists():
        return f"{out} exists"
    rec = si.load(bk.BINARY)
    if masked:
        rec = mask_bad_periods(rec, read_table(bk.BAD_PERIODS))
    win = rec.frame_slice(w * bk.WIN, min((w + 1) * bk.WIN, rec.get_num_frames()))
    win.reset_times()
    t0 = time.time()
    s = sort_chunk(win, out)
    return f"{out}: {s.get_num_units()} units in {(time.time() - t0) / 60:.1f} min"


def main():
    assert bk.REFERENCE == "aligned48", "run with BAKE_REFERENCE=aligned48"
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=14)
    a = ap.parse_args()
    jobs = [(m, w) for w in WINDOWS for m in (False, True)]
    with ProcessPoolExecutor(a.workers) as ex:
        for msg in ex.map(one, jobs):
            bk.log(msg)


if __name__ == "__main__":
    main()
