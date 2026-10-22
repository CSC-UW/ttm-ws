"""Per-block artifact and multi-unit suppression checks for a session with many laser blocks.

``opto_explore.py`` loads a whole short session into memory; this reads a
window around each train (block) from the store instead, so it runs on the
48 h recording. For every block: the raw and 300 to 6000 Hz pulse-triggered
averages at light onset and offset (the artifact check of
``opto_onset_zoom.py``), the LFP evoked potential, and multi-unit threshold
crossings per tetrode (5 MAD on each wire, MAD from the minute before the
block) as PSTHs and as rates before, during and after the block. Writes
``figures/blocks/opto_block<NN>.png`` per block, a summary figure across
blocks, ``results/opto_blocks_<E>.json`` and ``results/opto_blocks_<E>.csv``
(one row per block and tetrode).

Trains come from ``results/stimulation_blocks_<session>.csv`` when the
session has no hand log, otherwise from ``results/stimulation_trains.csv``.
On OPTO1 the artifact numbers agree with ``opto_onset_zoom.py`` to within
0.5 uV; the crossing rates differ in level from ``opto_explore.py`` (which
takes its MAD from the whole session) by up to a factor of 1.7 either way,
while the in-pulse to baseline ratios agree, so compare rates only within
this script's output.

Usage: python subjects/ttm-003/scripts/opto_blocks.py --subject TTM-003 --experiment TTM-003-48H
           [--blocks 24,25,26] [--margin-s 60]  (default: delivered, non-anomalous trains of 10 or more pulses) [--max-pulses 200] [--workers 8]
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as sps
import spikeinterface as si
from joblib import Parallel, delayed

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
PAD = 3000
THRESH_MAD = 5.0
DEAD_MS = 1.0


def load_blocks(session: str) -> pd.DataFrame:
    path = RESULTS / f"stimulation_blocks_{session}.csv"
    if path.exists():
        df = pd.read_csv(path)
    else:
        df = pd.read_csv(RESULTS / "stimulation_trains.csv")
        df = df[df.session == session]
        df["onset_h"] = df.onset_s / 3600
        df["anomalous"] = False
    return df.reset_index(drop=True)


def load_pulses(session: str) -> pd.DataFrame:
    df = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    return df[~df.glitch].reset_index(drop=True)


def filt(x: np.ndarray, sos: np.ndarray) -> np.ndarray:
    return sps.sosfiltfilt(sos, x, axis=0).astype(np.float32)


def crossings(xf: np.ndarray, mad: np.ndarray, groups: np.ndarray) -> dict[int, np.ndarray]:
    """Negative crossings of -THRESH_MAD * mad per tetrode (sample indices), merged within DEAD_MS."""
    dead = int(DEAD_MS * FS / 1000)
    out = {}
    for g in np.unique(groups):
        idx = []
        for c in np.flatnonzero(groups == g):
            below = xf[:, c] < -THRESH_MAD * mad[c]
            idx.append(np.flatnonzero(below[1:] & ~below[:-1]) + 1)
        idx = np.sort(np.concatenate(idx))
        if idx.size:
            idx = idx[np.r_[True, np.diff(idx) > dead]]
        out[int(g)] = idx
    return out


def snippets(x: np.ndarray, idx: np.ndarray, pre: int, post: int) -> np.ndarray:
    keep = (idx - pre >= 0) & (idx + post <= x.shape[0])
    return np.stack([x[i - pre:i + post] for i in idx[keep]])


def psth(spk: np.ndarray, events: np.ndarray, pre: int, post: int, bin_s: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(-pre, post + 1, int(bin_s * FS))
    counts = np.zeros(len(edges) - 1)
    for e in events:
        rel = spk[(spk >= e - pre) & (spk < e + post)] - e
        counts += np.histogram(rel, edges)[0]
    return edges[:-1] / FS + bin_s / 2, counts / (len(events) * bin_s)


def one_block(store: str, blk: pd.Series, pulses: pd.DataFrame, groups: np.ndarray, ch_names: list[str],
              margin_s: float, max_pulses: int, figdir: pathlib.Path, experiment: str) -> dict:
    rec = si.read_zarr(store)
    pp = pulses[pulses.train == blk.ttl_train]
    on_s, off_s = pp.onset_s.values[:max_pulses], pp.offset_s.values[:max_pulses]
    pulse_s = blk.pulse_ms_median / 1000
    a = max(int((blk.onset_s - margin_s) * FS) - PAD, 0)
    b = min(int((blk.offset_s + margin_s) * FS) + PAD, rec.get_num_frames())
    x = rec.get_traces(start_frame=a, end_frame=b, return_in_uV=True).astype(np.float32)
    xf = filt(x, sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos"))[PAD:-PAD]
    xl = filt(x, sps.butter(4, 300, btype="lowpass", fs=FS, output="sos"))[PAD:-PAD]
    x = x[PAD:-PAD]
    a += PAD
    on = np.round(on_s * FS).astype(np.int64) - a
    off = np.round(off_s * FS).astype(np.int64) - a
    first = max(on[0] - int(margin_s * FS), 0)
    pre_block = xf[first:on[0]]
    mad = np.median(np.abs(pre_block), axis=0) / 0.6745
    mua = crossings(xf, mad, groups)
    tetrodes = sorted(mua)

    pre, post = int(0.05 * FS), int((pulse_s + 0.1) * FS)
    raw_on = snippets(x, on, pre, post).mean(0)
    tt = (np.arange(raw_on.shape[0]) - pre) / FS
    hf_pre, hf_post = int(0.005 * FS), int(0.012 * FS)
    hf_on = snippets(xf, on, hf_pre, hf_post).mean(0)
    hf_off = snippets(xf, off, hf_pre, hf_post).mean(0)
    hf_t = (np.arange(hf_on.shape[0]) - hf_pre) / FS * 1000
    lfp_on = snippets(xl, on, pre, post).mean(0)
    base = raw_on[tt < -0.005].mean(0)
    sustained = raw_on[(tt >= 0.02) & (tt < pulse_s)].mean(0) - base
    noise = hf_on[hf_t < -1].std(0)
    edge_on = np.abs(hf_on[(hf_t >= -0.2) & (hf_t < 0.7)]).max(0)
    edge_off = np.abs(hf_off[(hf_t >= -0.2) & (hf_t < 0.7)]).max(0)
    evoked = np.abs(hf_on[(hf_t >= 1) & (hf_t < 6)]).max(0)

    t0, t1 = on[0], off[-1]
    rows = []
    for g in tetrodes:
        spk = mua[g]
        tb, r = psth(spk, on, pre, post, 0.002)
        rows.append({
            "block": int(blk.ttl_train), "onset_h": round(float(blk.onset_h), 4), "tetrode": f"TT{g + 1}",
            "baseline_hz": float(r[tb < 0].mean()),
            "early_0_10ms_hz": float(r[(tb >= 0) & (tb < 0.010)].mean()),
            "pulse_10ms_end_hz": float(r[(tb >= 0.010) & (tb < pulse_s)].mean()),
            "post_0_50ms_hz": float(r[(tb >= pulse_s) & (tb < pulse_s + 0.05)].mean()),
            "post_50_100ms_hz": float(r[(tb >= pulse_s + 0.05) & (tb < pulse_s + 0.1)].mean()),
            "pre60_hz": float(((spk >= t0 - 60 * FS) & (spk < t0)).sum() / min(60, t0 / FS)),
            "during_hz": float(((spk >= t0) & (spk < t1)).sum() / ((t1 - t0) / FS)),
            "post60_hz": float(((spk >= t1) & (spk < t1 + 60 * FS)).sum() / 60),
            "mad_uV_median": float(np.median(mad[groups == g])),
        })
    for row in rows:
        row["pulse_over_baseline"] = row["pulse_10ms_end_hz"] / row["baseline_hz"] if row["baseline_hz"] > 0 else float("nan")

    label = f"block{int(blk.ttl_train):02d}"
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    ax = axes[0, 0]
    v = np.percentile(np.abs(raw_on), 99)
    im = ax.imshow(raw_on.T, aspect="auto", extent=[tt[0] * 1000, tt[-1] * 1000, len(ch_names) - 0.5, -0.5],
                   cmap="RdBu_r", vmin=-v, vmax=v)
    ax.axvline(0, color="k", lw=0.5); ax.axvline(pulse_s * 1000, color="k", lw=0.5)
    ax.set_title(f"raw pulse-triggered average (uV), n={len(on)}"); ax.set_xlabel("ms from onset"); ax.set_ylabel("channel (tetrode order)")
    plt.colorbar(im, ax=ax)
    ax = axes[0, 1]
    for g in tetrodes:
        ax.plot(tt * 1000, lfp_on[:, groups == g].mean(1), lw=0.8, label=f"TT{g + 1}")
    ax.axvspan(0, pulse_s * 1000, color="c", alpha=0.15)
    ax.set_title("LFP (<300 Hz) evoked potential, mean per tetrode"); ax.set_xlabel("ms"); ax.set_ylabel("uV"); ax.legend(ncol=4, fontsize=7)
    ax = axes[1, 0]
    ax.plot(hf_t, hf_on, lw=0.5, alpha=0.6)
    ax.set_title("300-6000 Hz trial average at onset (all channels)"); ax.set_xlabel("ms from onset"); ax.set_ylabel("uV")
    ax2 = ax.inset_axes([0.6, 0.6, 0.38, 0.38]); ax2.plot(hf_t, hf_off, lw=0.5, alpha=0.6); ax2.set_title("at offset", fontsize=8)
    ax = axes[1, 1]
    for g in tetrodes:
        tb, r = psth(mua[g], on, pre, post, 0.002)
        ax.plot(tb * 1000, r, lw=0.8, label=f"TT{g + 1}")
    ax.axvspan(0, pulse_s * 1000, color="c", alpha=0.15)
    ax.set_title(f"MUA PSTH per tetrode (2 ms bins, {THRESH_MAD:g} MAD crossings)"); ax.set_xlabel("ms from onset"); ax.set_ylabel("crossings/s"); ax.legend(ncol=4, fontsize=7)
    fig.suptitle(f"{experiment} {label} at {blk.onset_h:.2f} h ({blk.get('onset_local', '')}): {len(pp)} pulses of "
                 f"{blk.pulse_ms_median:.0f} ms at {blk.rate_hz:.2f} Hz, pattern {blk.ttl_pattern}")
    fig.tight_layout(); fig.savefig(figdir / f"opto_{label}.png", dpi=100); plt.close(fig)

    return {
        "block": int(blk.ttl_train), "onset_h": float(blk.onset_h), "onset_local": str(blk.get("onset_local", "")),
        "n_pulses": int(len(pp)), "n_averaged": int(len(on)), "pulse_ms": float(blk.pulse_ms_median), "pattern": str(blk.ttl_pattern),
        "hf_noise_uV_median": float(np.median(noise)),
        "onset_edge_uV_median": float(np.median(edge_on)), "onset_edge_uV_max": float(edge_on.max()),
        "offset_edge_uV_median": float(np.median(edge_off)), "offset_edge_uV_max": float(edge_off.max()),
        "evoked_1_6ms_uV_max": float(evoked.max()),
        "evoked_1_6ms_uV_by_tetrode": {f"TT{g + 1}": float(evoked[groups == g].max()) for g in tetrodes},
        "sustained_shift_uV_median": float(np.median(sustained)), "sustained_shift_uV_max": float(np.abs(sustained).max()),
        "sustained_shift_uV_by_channel": dict(zip(ch_names, np.round(sustained, 1).tolist())),
        "tetrodes": rows,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--blocks", default=None, help="comma-separated ttl_train numbers; default all non-anomalous")
    p.add_argument("--margin-s", type=float, default=60.0)
    p.add_argument("--max-pulses", type=int, default=200)
    p.add_argument("--workers", type=int, default=8)
    a = p.parse_args()
    exp = experiment_from_args(a)
    figdir = exp.figures_dir / "blocks"
    figdir.mkdir(parents=True, exist_ok=True)
    blocks = load_blocks(exp.session)
    pulses = load_pulses(exp.session)
    if a.blocks:
        want = [int(s) for s in a.blocks.split(",")]
        blocks = blocks[blocks.ttl_train.isin(want)]
    else:
        blocks = blocks[(~blocks.anomalous) & (blocks.n_pulses >= 10)]
        if "delivered" in blocks:
            blocks = blocks[blocks.delivered]
    blocks = blocks.reset_index(drop=True)
    rec = si.read_zarr(str(exp.store_path()))
    groups = np.asarray(rec.get_property("group")).astype(int)
    ch_names = [str(c) for c in rec.get_channel_ids()]
    print(f"{len(blocks)} blocks, {len(pulses)} pulses, {a.workers} workers")
    out = Parallel(n_jobs=a.workers)(delayed(one_block)(str(exp.store_path()), blk, pulses, groups, ch_names,
                                                        a.margin_s, a.max_pulses, figdir, exp.experiment)
                                     for _, blk in blocks.iterrows())
    out = sorted(out, key=lambda d: d["block"])
    rows = pd.DataFrame([r for d in out for r in d["tetrodes"]])
    rows.round(4).to_csv(RESULTS / f"opto_blocks_{exp.experiment}.csv", index=False)
    (RESULTS / f"opto_blocks_{exp.experiment}.json").write_text(json.dumps(
        {"experiment": exp.experiment, "session": exp.session, "thresh_mad": THRESH_MAD, "blocks": out}, indent=1) + "\n")

    # Summary across blocks: artifact edge, suppression per tetrode, baseline MUA per tetrode.
    tets = sorted(rows.tetrode.unique(), key=lambda s: int(s[2:]))
    hours = [d["onset_h"] for d in out]
    sup = rows.pivot(index="tetrode", columns="block", values="pulse_over_baseline").loc[tets]
    base = rows.pivot(index="tetrode", columns="block", values="baseline_hz").loc[tets]
    fig, axes = plt.subplots(3, 1, figsize=(16, 13), gridspec_kw={"height_ratios": [1, 2, 2]})
    ax = axes[0]
    ax.plot(hours, [d["onset_edge_uV_max"] for d in out], "o-", label="onset edge, max channel")
    ax.plot(hours, [d["offset_edge_uV_max"] for d in out], "s-", label="offset edge, max channel")
    ax.plot(hours, [d["hf_noise_uV_median"] for d in out], "k--", label="noise floor of the average, median channel")
    ax.plot(hours, [d["evoked_1_6ms_uV_max"] for d in out], "^-", label="evoked 1-6 ms, max channel")
    ax.set_ylabel("uV"); ax.set_xlabel("session hour"); ax.legend(fontsize=8, ncol=4); ax.set_title("300-6000 Hz pulse-triggered average: light-edge amplitude per block (artifact check)")
    for ax, mat, title, vmax in ((axes[1], sup, "MUA rate 10 ms to pulse end, fraction of the 50 ms baseline", 1.5),
                                 (axes[2], base, "baseline MUA rate (crossings/s, -50 to 0 ms)", None)):
        im = ax.imshow(mat.values, aspect="auto", cmap="viridis", vmin=0, vmax=vmax)
        ax.set_yticks(range(len(tets))); ax.set_yticklabels(tets, fontsize=7)
        ax.set_xticks(range(len(hours))); ax.set_xticklabels([f"{h:.1f}" for h in hours], fontsize=6, rotation=90)
        ax.set_xlabel("block onset (session hour)"); ax.set_title(title); plt.colorbar(im, ax=ax)
    fig.suptitle(f"{exp.experiment}: {len(out)} laser blocks")
    fig.tight_layout(); fig.savefig(exp.figures_dir / "opto_blocks_summary.png", dpi=110); plt.close(fig)
    print(rows.groupby("tetrode")[["baseline_hz", "pulse_over_baseline"]].median().loc[tets].round(3).to_string())
    print(f"wrote {RESULTS / f'opto_blocks_{exp.experiment}.csv'} and figures in {figdir}")


if __name__ == "__main__":
    main()
