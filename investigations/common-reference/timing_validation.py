"""Whether per-window wire timing estimated as ``ttm.reference`` does carries to held-out data, in three animals.

Per 1800 s window, the six blocks of 20 s that ``estimate_wire_timing_window`` reads (bandpassed 300
to 6000 Hz, unreferenced), each fitted on its own with ``estimate_wire_timing``. Per window: how many
blocks have 3000 event samples or more; the spread of a block's offset between blocks (median over
wires of the MAD over qualifying blocks); and a held-out test: offsets and gains from the qualifying
blocks among 0, 2, 4 (the per-wire median) applied to blocks 1, 3, 5, against the global median.
The measure is each wire's RMS within the event mask (``event_mask``) of the held-out block, as the
ratio aligned / global, and its median and 90th percentile over wires. The same with the estimate of
the window before (what ``fill_missing_windows`` would borrow), and of the window 8 before (4 h).

Recordings and windows as ``common_events.py``: TTM-003's 48 h (96 windows), TTM-001's NOD (96),
TTM-002's SD (every fourth). Writes ``results/timing_validation.csv`` (per window) and
``mp_aligned/timing_validation_blocks.npz`` beside the session's sortings (per block offsets, gains, event counts).

Run: .venv/bin/python investigations/common-reference/timing_validation.py
"""
from __future__ import annotations

import pathlib
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from spikeinterface.preprocessing.phase_shift import apply_frequency_shift

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from common_events import BLOCK, FS, N_BLOCKS, WIN, open_rec  # noqa: E402
from ttm.reference import event_mask, estimate_wire_timing  # noqa: E402

HERE = pathlib.Path(__file__).parent
BLOCKS = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix/mp_aligned") / "timing_validation_blocks.npz"
MIN_EV = 3000
EDGE = 3000


def blocks_of(name, w):
    rec = open_rec(name)
    n = rec.get_num_frames()
    s, e = w * WIN, min((w + 1) * WIN, n)
    gain = float(rec.get_channel_gains()[0])
    out = []
    for i in range(N_BLOCKS):
        c = s + (e - s) * (i + 1) // (N_BLOCKS + 1)
        out.append(rec.get_traces(start_frame=c - BLOCK // 2, end_frame=c + BLOCK // 2).astype(np.float64) * gain)
    return out


def fit_window(args):
    name, w = args
    bl = blocks_of(name, w)
    taus, gains, nev = [], [], []
    for x in bl:
        t, g, n = estimate_wire_timing(x, min_event_samples=MIN_EV)
        taus.append(t)
        gains.append(g)
        nev.append(n)
    return name, w, np.array(taus), np.array(gains), np.array(nev)


def held_out(args):
    name, w, estimates = args
    bl = blocks_of(name, w)
    out = {}
    for i in (1, 3, 5):
        x = bl[i]
        m = np.median(x, axis=1)
        ev = event_mask(m)
        ev[:EDGE] = ev[-EDGE:] = False
        if ev.sum() < MIN_EV:
            continue
        glob = np.sqrt(((x - m[:, None])[ev] ** 2).mean(0))
        for key, est in estimates.items():
            if est is None:
                continue
            tau, g = est
            y = apply_frequency_shift(x, -tau, axis=0)
            r = np.sqrt(((y - g[None, :] * np.median(y, axis=1, keepdims=True))[ev] ** 2).mean(0))
            out.setdefault(key, []).append(r / glob)
    return name, w, {k: np.mean(v, axis=0) for k, v in out.items()}


def median_of(taus, gains, nev, idx):
    ok = [i for i in idx if nev[i] >= MIN_EV]
    if not ok:
        return None
    return np.median(taus[ok], axis=0), np.median(gains[ok], axis=0)


def main():
    jobs = [("TTM-003 48H", w) for w in range(96)] + [("TTM-001 NOD", w) for w in range(96)] \
        + [("TTM-002 SD", w) for w in range(0, 96, 4)]
    with ProcessPoolExecutor(16) as ex:
        fits = list(ex.map(fit_window, jobs))
    F = {(n, w): (t, g, e) for n, w, t, g, e in fits}
    np.savez(BLOCKS,
             **{f"{n.split()[0]}_w{w:02d}_{k}": v for (n, w), tv in F.items()
                for k, v in zip(("tau", "gain", "nev"), tv)})
    step = {"TTM-003 48H": 1, "TTM-001 NOD": 1, "TTM-002 SD": 4}
    ho_jobs = []
    for (n, w), (t, g, e) in F.items():
        est = {"same window": median_of(t, g, e, (0, 2, 4))}
        for lag, key in ((1, "window before"), (8, "8 windows before")):
            p = F.get((n, w - lag * step[n]))
            est[key] = median_of(*p, range(N_BLOCKS)) if p else None
        ho_jobs.append((n, w, est))
    with ProcessPoolExecutor(16) as ex:
        ho = list(ex.map(held_out, ho_jobs))
    rows = []
    for (n, w, ratios) in ho:
        t, g, e = F[(n, w)]
        ok = e >= MIN_EV
        spread = float(np.median(np.median(np.abs(t[ok] - np.median(t[ok], axis=0)), axis=0)) / FS * 1e6) \
            if ok.sum() >= 3 else np.nan
        row = {"recording": n, "window": w, "blocks_ok": int(ok.sum()), "block_spread_us": spread}
        for k, r in ratios.items():
            row[f"{k}: median"] = float(np.median(r))
            row[f"{k}: p90"] = float(np.percentile(r, 90))
            row[f"{k}: worst"] = float(r.max())
        rows.append(row)
    df = pd.DataFrame(rows).sort_values(["recording", "window"]).round(4)
    df.to_csv(HERE / "results" / "timing_validation.csv", index=False)
    pd.set_option("display.width", 250)
    for n, d in df.groupby("recording", sort=False):
        print(f"== {n}: {len(d)} windows; blocks with enough events per window: "
              f"{d.blocks_ok.value_counts().sort_index().to_dict()}")
        print(f"  block-to-block spread of a wire's offset, median over windows {d.block_spread_us.median():.2f} us "
              f"(90th percentile {d.block_spread_us.quantile(0.9):.2f})")
        cols = [c for c in d.columns if ": " in c]
        print("  held-out RMS at events, aligned over global (per window, then median [10th, 90th percentile] "
              "over windows):")
        for c in cols:
            v = d[c].dropna()
            print(f"    {c:28s} {v.median():.3f} [{v.quantile(0.1):.3f}, {v.quantile(0.9):.3f}]  n={len(v)}")


if __name__ == "__main__":
    main()
