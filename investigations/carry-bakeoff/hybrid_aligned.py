"""The hybrid on the aligned-reference binary, with the same injected ground truth.

``make_hybrid.py`` added the ground-truth units to the first 12 windows of ``mp_baseline/binary``
(global median). This adds the same injected signal, the old hybrid minus those frames of the old
binary, to the first 12 windows referenced with the aligned median
(``mp_aligned/binary_w00_11``), 60 s at a time with 32 processes, and writes it to
``v2_aligned/hybrid_binary`` with the aligned binary's metadata; ``gt.npz`` and ``gt_plan.json``
are copied unchanged. Every ground-truth spike, template and size is the same in both hybrids, so
the two differ only in the background.

Run: .venv/bin/python investigations/carry-bakeoff/hybrid_aligned.py
"""
from __future__ import annotations

import json
import pathlib
import shutil
from concurrent.futures import ProcessPoolExecutor

import numpy as np

SORTINGS = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix")
OLD_BINARY = SORTINGS / "mp_baseline" / "binary"
OLD_HYBRID = SORTINGS / "carry_bakeoff" / "v2" / "hybrid_binary"
NEW_BINARY = SORTINGS / "mp_aligned" / "binary_w00_11"
NEW_BAKE = SORTINGS / "carry_bakeoff" / "v2_aligned"
NEW_HYBRID = NEW_BAKE / "hybrid_binary"
RAW = "traces_cached_seg0.raw"
N_CH = 64
STEP = 60 * 30000


def meta(folder: pathlib.Path) -> dict:
    return json.loads((folder / "binary.json").read_text())


def add_chunk(args):
    start, n_frames = args
    shape = (n_frames, N_CH)
    old_b = np.memmap(OLD_BINARY / RAW, np.float32, "r", shape=(meta_frames(OLD_BINARY), N_CH))
    old_h = np.memmap(OLD_HYBRID / RAW, np.float32, "r", shape=shape)
    new_b = np.memmap(NEW_BINARY / RAW, np.float32, "r", shape=shape)
    out = np.memmap(NEW_HYBRID / RAW, np.float32, "r+", shape=shape)
    e = min(start + STEP, n_frames)
    out[start:e] = new_b[start:e] + (old_h[start:e] - old_b[start:e])
    out.flush()
    return start


def meta_frames(folder: pathlib.Path) -> int:
    return (folder / RAW).stat().st_size // (4 * N_CH)


def main():
    for f in (OLD_BINARY, OLD_HYBRID, NEW_BINARY):
        k = meta(f)["kwargs"]
        assert k["num_channels"] == N_CH and np.dtype(k["dtype"]) == np.float32, f
    ids = [meta(f)["kwargs"]["channel_ids"] for f in (OLD_BINARY, OLD_HYBRID, NEW_BINARY)]
    assert ids[0] == ids[1] == ids[2], "channel order differs"
    n = meta_frames(NEW_BINARY)
    assert n == meta_frames(OLD_HYBRID), (n, meta_frames(OLD_HYBRID))
    NEW_BAKE.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(NEW_HYBRID, ignore_errors=True)
    shutil.copytree(NEW_BINARY, NEW_HYBRID, ignore=shutil.ignore_patterns(RAW, "provenance.json"))
    np.memmap(NEW_HYBRID / RAW, np.float32, "w+", shape=(n, N_CH)).flush()
    with ProcessPoolExecutor(32) as ex:
        list(ex.map(add_chunk, [(s, n) for s in range(0, n, STEP)]))
    for f in ("gt.npz", "gt_plan.json"):
        shutil.copy2(OLD_HYBRID.parent / f, NEW_BAKE / f)
    print(f"{NEW_HYBRID}: {n} frames")


if __name__ == "__main__":
    main()
