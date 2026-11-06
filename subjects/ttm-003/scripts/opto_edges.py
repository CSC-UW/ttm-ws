"""Spike-band average over every delivered pulse of a session at light onset and offset, for the artifact check.

``opto_blocks.py`` and ``opto_off_periods.py`` average per train, where a train's noise floor is a few uV.
This reads 100 ms around every onset and offset from the store, bandpasses it 300 to 6000 Hz, and averages
over all pulses, so the floor drops by the square root of the number of trains. A photoelectric artifact
is locked to the light edge (0 to 0.7 ms); light-evoked spikes start about 1 ms later and jitter. The same
average at sham times (each onset minus 500 ms, mid-way between pulses at 1 Hz) gives the floor's spread.

Writes ``results/opto_edges_<E>.json`` (per channel: the largest absolute value from -0.2 to 0.7 ms at
onset and offset, the same at sham times, and the noise of the average before -1 ms) and
``<experiment>/figures/opto_edges.png``.

Usage: python subjects/ttm-003/scripts/opto_edges.py --subject TTM-003 --experiment TTM-003-48H-2 [--workers 16]
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
import numpy as np
import pandas as pd
import scipy.signal as sps
import spikeinterface as si
from joblib import Parallel, delayed

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ttm.cli._common import add_experiment_args, experiment_from_args
from ttm.sorting import bad_channel_ids

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
HALF = int(0.05 * FS)
LO, HI = int(0.005 * FS), int(0.012 * FS)


def chunk_sum(store: str, drop: list[str], times: np.ndarray) -> tuple[np.ndarray, int]:
    rec = si.read_zarr(store)
    if drop:
        rec = rec.remove_channels(drop)
    sos = sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos")
    acc = None
    n = 0
    for t in times:
        f = round(t * FS)
        if f - HALF < 0 or f + HALF > rec.get_num_frames():
            continue
        x = rec.get_traces(start_frame=f - HALF, end_frame=f + HALF, return_in_uV=True).astype(np.float64)
        y = sps.sosfiltfilt(sos, x, axis=0)[HALF - LO:HALF + HI]
        acc = y if acc is None else acc + y
        n += 1
    return acc, n


def pooled(store, drop, times, workers):
    parts = Parallel(n_jobs=workers)(delayed(chunk_sum)(store, drop, c) for c in np.array_split(times, workers * 4))
    tot = sum(p[0] for p in parts if p[0] is not None)
    n = sum(p[1] for p in parts)
    return tot / n, n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--workers", type=int, default=16)
    a = ap.parse_args()
    exp = experiment_from_args(a)
    rec = si.read_zarr(str(exp.store_path()))
    drop = [str(c) for c in bad_channel_ids(rec, exp.badChannels, exp.badTetrodes)]
    ch = [str(c) for c in rec.channel_ids if str(c) not in drop]
    blocks = pd.read_csv(RESULTS / f"stimulation_blocks_{exp.session}.csv")
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    pulses = pulses[~pulses.glitch & pulses.train.isin(blocks.loc[blocks.delivered, "ttl_train"])]
    on, off = pulses.onset_s.to_numpy(), pulses.offset_s.to_numpy()
    store = str(exp.store_path())
    avg = {k: pooled(store, drop, t, a.workers) for k, t in (("onset", on), ("offset", off), ("sham", on - 0.5))}
    tt = (np.arange(LO + HI) - LO) / FS * 1000
    edge = (tt >= -0.2) & (tt < 0.7)
    pre = tt < -1
    out = {"experiment": exp.experiment, "pulses": int(on.size), "channels": ch}
    for k, (m, n) in avg.items():
        out[f"{k}_n"] = int(n)
        out[f"{k}_edge_uV"] = np.round(np.abs(m[edge]).max(0), 2).tolist()
        out[f"{k}_noise_uV"] = np.round(m[pre].std(0), 2).tolist()
    for k in ("onset", "offset", "sham"):
        e, z = np.array(out[f"{k}_edge_uV"]), np.array(out[f"{k}_noise_uV"])
        out[f"{k}_edge_over_noise_max"] = float((e / z).max())
        out[f"{k}_edge_uV_max"] = float(e.max())
        out[f"{k}_edge_uV_median"] = float(np.median(e))
    m_on = avg["onset"][0]
    ev = (tt >= 1) & (tt < 6)
    out["evoked_1_6ms_uV_max"] = float(np.abs(m_on[ev]).max())
    out["evoked_1_6ms_channel"] = ch[int(np.abs(m_on[ev]).max(0).argmax())]
    (RESULTS / f"opto_edges_{exp.experiment}.json").write_text(json.dumps(out, indent=1) + "\n")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    for ax, k in zip(axes, ("onset", "offset", "sham")):
        ax.plot(tt, avg[k][0], lw=0.5, alpha=0.7)
        ax.axvspan(-0.2, 0.7, color="y", alpha=0.2)
        ax.set_title(f"{k} (n={avg[k][1]}): edge max {out[f'{k}_edge_uV_max']:.1f} uV, "
                     f"max edge/noise {out[f'{k}_edge_over_noise_max']:.1f}", fontsize=9)
        ax.set_xlabel("ms")
    axes[0].set_ylabel("uV (300-6000 Hz, every channel)")
    fig.suptitle(f"{exp.experiment}: spike-band average over every delivered pulse")
    fig.tight_layout(); fig.savefig(exp.figures_dir / "opto_edges.png", dpi=100); plt.close(fig)
    print(json.dumps({k: v for k, v in out.items() if not isinstance(v, list)}, indent=1))


if __name__ == "__main__":
    main()
