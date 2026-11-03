"""Pooling blocks across neighbouring windows, against the window's own blocks, on held-out blocks.

``timing_validation.py`` found that an estimate from the window before, or from 4 h before, does as
well on held-out blocks as the window's own, and that windows where one block decides the estimate
can leave a wire 1.2 to 2.7 times worse than the global median (TTM-002 windows 32 and 84). Here,
from its per-block fits (``mp_aligned/timing_validation_blocks.npz`` beside the session's sortings), each window's estimate is the
per-wire median of the qualifying blocks (3000 event samples or more) among:

- own: blocks 0, 2, 4 of the window (as in ``timing_validation.py``);
- own, at least 2: the same, NaN (not scored) with fewer than 2 qualifying blocks;
- pooled +-2: blocks 0, 2, 4 of the window and all blocks of the two windows (TTM-002: sampled
  windows) on each side;
- pooled +-2, at least 3: the same, falling back to the global median (shift 0, gain 1) with fewer
  than 3 qualifying blocks.

Each is applied to blocks 1, 3, 5 of the window and scored as in ``timing_validation.py``: per wire,
RMS within the event mask, aligned over global. Writes ``results/timing_pooling.csv``.

Run: .venv/bin/python investigations/common-reference/timing_pooling.py
"""
from __future__ import annotations

import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from timing_validation import MIN_EV, N_BLOCKS, held_out  # noqa: E402

HERE = pathlib.Path(__file__).parent
BLOCKS = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix/mp_aligned") / "timing_validation_blocks.npz"
RECS = {"TTM-003 48H": range(96), "TTM-001 NOD": range(96), "TTM-002 SD": range(0, 96, 4)}


def pool(blocks, min_n):
    t = [tb for tb, gb, nb in blocks if nb >= MIN_EV]
    g = [gb for tb, gb, nb in blocks if nb >= MIN_EV]
    if len(t) < min_n:
        return None
    return np.median(t, axis=0), np.median(g, axis=0)


def main():
    Z = np.load(BLOCKS)
    jobs = []
    for name, wins in RECS.items():
        key = name.split()[0]
        wins = list(wins)
        F = {w: tuple(Z[f"{key}_w{w:02d}_{k}"] for k in ("tau", "gain", "nev")) for w in wins}

        def blk(w, idx):
            t, g, e = F[w]
            return [(t[i], g[i], e[i]) for i in idx]

        for j, w in enumerate(wins):
            own = blk(w, (0, 2, 4))
            nb = [blk(v, range(N_BLOCKS)) for v in wins[max(0, j - 2):j] + wins[j + 1:j + 3]]
            pooled = own + [b for x in nb for b in x]
            n_ch = F[w][0].shape[1]
            identity = (np.zeros(n_ch), np.ones(n_ch))
            jobs.append((name, w, {"own": pool(own, 1), "own, at least 2": pool(own, 2),
                                   "pooled +-2": pool(pooled, 1),
                                   "pooled +-2, at least 3": pool(pooled, 3) or identity}))
    with ProcessPoolExecutor(16) as ex:
        res = list(ex.map(held_out, jobs))
    rows = []
    for name, w, ratios in res:
        for k, r in ratios.items():
            rows.append({"recording": name, "window": w, "estimate": k, "median": float(np.median(r)),
                         "p90": float(np.percentile(r, 90)), "worst": float(r.max()),
                         "wires_worse_5pct": int((r > 1.05).sum())})
    df = pd.DataFrame(rows).round(4)
    df.to_csv(HERE / "results" / "timing_pooling.csv", index=False)
    pd.set_option("display.width", 250)
    for name, d in df.groupby("recording", sort=False):
        print(f"== {name}")
        s = d.groupby("estimate", sort=False).agg(
            windows=("window", "size"), median=("median", "median"), p90=("p90", "median"),
            worst_median=("worst", "median"), worst_max=("worst", "max"),
            windows_with_a_wire_5pct_worse=("wires_worse_5pct", lambda x: int((x > 0).sum())))
        print(s.round(3).to_string())


if __name__ == "__main__":
    main()
