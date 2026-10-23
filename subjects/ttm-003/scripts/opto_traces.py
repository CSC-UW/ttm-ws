"""Minimally processed traces of one tetrode around the laser pulses of chosen blocks.

No sorting: the stored traces of the tetrode's four wires, bandpassed 300 to
6000 Hz, in uV, either as they are (``--reference none``) or with the median
of every channel of the other tetrodes (bad channels excluded) subtracted
(``--reference other-tetrodes``, figure suffix ``_ref``). Per block one
figure under ``<experiment>/figures/opto_traces/``:

* left, one row per example pulse (``--pulses``, 1-based positions in the
  block; negative counts from the end): the four wires stacked, -50 ms to
  pulse end + 100 ms, the pulse shaded, a fixed uV scale across blocks;
* right, top: every pulse of the block as one row of an image of the wire
  given by ``--wire`` (default: the wire with the deepest pulse-triggered
  trough in the first block plotted), same time axis, so a train of spikes
  that recurs at every pulse shows as vertical structure;
* right, bottom: the pulse-triggered mean of each wire, and the per-pulse
  count of troughs below ``--thresh-uv`` (or ``--thresh-mad`` times the MAD of the second before the block) on the
  image wire, inside the pulse against the 50 ms before it, over the quietest third of pulses and over all.

Usage: python subjects/ttm-003/scripts/opto_traces.py --subject TTM-003 --experiment TTM-003-48H --tetrode TT6 --blocks 1,4,25,26
"""

from __future__ import annotations

import argparse
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
PRE_MS, POST_MS = 50, 100


def load_blocks(session: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"stimulation_blocks_{session}.csv")


def load_pulses(session: str) -> pd.DataFrame:
    df = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    return df[~df.glitch].reset_index(drop=True)


def tetrode_index(name: str) -> int:
    return int(name.upper().replace("TT", "")) - 1


def block_traces(store: str, blk: pd.Series, pulses: pd.DataFrame, group: int, margin_s: float,
                 reference: str = "none", bad: tuple[str, ...] = ()):
    """Bandpassed (n_samples, 4) traces of the block with ``margin_s`` before and after, and the
    pulse onsets and offsets as frames into them. With ``reference="other-tetrodes"`` the median
    over the bandpassed channels of every other tetrode (``bad`` excluded) is subtracted."""
    rec = si.read_zarr(store)
    groups = np.asarray(rec.get_property("group")).astype(int)
    names = np.asarray([str(c) for c in rec.get_channel_ids()])
    ids = list(names[groups == group])
    a = max(int((blk.onset_s - margin_s) * FS) - PAD, 0)
    b = min(int((blk.offset_s + margin_s) * FS) + PAD, rec.get_num_frames())
    sos = sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos")
    if reference == "none":
        x = rec.get_traces(start_frame=a, end_frame=b, channel_ids=ids, return_in_uV=True).astype(np.float32)
        xf = sps.sosfiltfilt(sos, x, axis=0).astype(np.float32)[PAD:-PAD]
    elif reference == "other-tetrodes":
        x = rec.get_traces(start_frame=a, end_frame=b, return_in_uV=True).astype(np.float32)
        xf_all = sps.sosfiltfilt(sos, x, axis=0).astype(np.float32)[PAD:-PAD]
        others = (groups != group) & ~np.isin(names, bad)
        ref = np.median(xf_all[:, others], axis=1, keepdims=True)
        xf = (xf_all[:, groups == group] - ref).astype(np.float32)
    else:
        raise ValueError(reference)
    a += PAD
    pp = pulses[pulses.train == blk.ttl_train]
    on = np.round(pp.onset_s.values * FS).astype(np.int64) - a
    off = np.round(pp.offset_s.values * FS).astype(np.int64) - a
    return xf, ids, on, off


def troughs(x: np.ndarray, thresh: float, dead: int = 30) -> np.ndarray:
    """Indices of local minima below ``-thresh`` with ``dead`` samples between them."""
    below = np.flatnonzero((x[1:-1] < -thresh) & (x[1:-1] <= x[:-2]) & (x[1:-1] < x[2:])) + 1
    keep = []
    last = -dead
    for i in below:
        if i - last >= dead:
            keep.append(i)
            last = i
    return np.asarray(keep, dtype=np.int64)


def one_block(store: str, blk: pd.Series, pulses: pd.DataFrame, group: int, tetrode: str,
              examples: list[int], wire: int | None, thresh_mad: float, thresh_uv: float | None,
              scale_uv: float, image_uv: float, figdir: pathlib.Path, experiment: str,
              reference: str = "none", bad: tuple[str, ...] = ()) -> dict:
    xf, ids, on, off = block_traces(store, blk, pulses, group, margin_s=1.0, reference=reference, bad=bad)
    suffix = "_ref" if reference == "other-tetrodes" else ""
    ref_label = ("bandpass 300 to 6000 Hz, minus the median of the other tetrodes" if reference == "other-tetrodes"
                 else "bandpass 300 to 6000 Hz only")
    pre = int(PRE_MS * FS / 1000)
    post = int(POST_MS * FS / 1000)
    pulse = int(np.median(off - on))
    span = pre + pulse + post
    t_ms = (np.arange(span) - pre) / FS * 1000
    seg = np.stack([xf[o - pre:o - pre + span] for o in on if o - pre >= 0 and o - pre + span <= len(xf)])
    mean = seg.mean(0)
    mad = np.asarray([round(float(m), 1) for m in np.median(np.abs(xf[:on[0]]), axis=0) / 0.6745])
    if wire is None:
        wire = int(np.argmin(mean[pre:pre + 6 * 30].min(0)))
    thr = thresh_uv if thresh_uv else thresh_mad * mad[wire]
    thr_label = f"{thr:.0f} uV" if thresh_uv else f"{thresh_mad:g} MAD = {thr:.0f} uV"
    n_in, n_pre = [], []
    for s in seg:
        tr = troughs(s[:, wire], thr)
        n_pre.append(int(((tr >= 0) & (tr < pre)).sum()))
        n_in.append(int(((tr >= pre + 30) & (tr < pre + pulse)).sum()))
    n_in, n_pre = np.asarray(n_in), np.asarray(n_pre)
    rate_in = n_in / ((pulse - 30) / FS)
    rate_pre = n_pre / (pre / FS)
    pre_rms = np.sqrt((seg[:, :pre, :] ** 2).mean(axis=(1, 2)))
    quiet = pre_rms <= np.percentile(pre_rms, 100 / 3)

    pos = [(e - 1) if e > 0 else (len(on) + e) for e in examples]
    pos = [p for p in pos if 0 <= p < len(seg)]
    fig = plt.figure(figsize=(18, 3.2 * max(len(pos), 3)))
    gs = fig.add_gridspec(max(len(pos), 3), 2, width_ratios=[1.4, 1])
    for r, pidx in enumerate(pos):
        ax = fig.add_subplot(gs[r, 0])
        for w in range(4):
            ax.plot(t_ms, seg[pidx][:, w] - w * scale_uv, lw=0.5, color="k" if w == wire else "0.4")
            ax.text(t_ms[0], -w * scale_uv + scale_uv * 0.25, ids[w], fontsize=8, va="bottom")
        ax.axvspan(0, pulse / FS * 1000, color="c", alpha=0.15, lw=0)
        ax.axhline(-wire * scale_uv - thr, color="r", lw=0.5, ls="--")
        ax.set_ylim(-3.5 * scale_uv, scale_uv)
        ax.set_yticks([-w * scale_uv for w in range(4)], ["0"] * 4)
        ax.set_xlim(t_ms[0], t_ms[-1])
        ax.set_title(f"pulse {pidx + 1} of {len(seg)}{' (quiet)' if quiet[pidx] else ''}: {n_pre[pidx]} troughs below "
                     f"{thr_label} on {ids[wire]} in the 50 ms before, {n_in[pidx]} in the pulse after 1 ms", fontsize=9)
        if r == len(pos) - 1:
            ax.set_xlabel("ms from pulse onset")
    ax = fig.add_subplot(gs[: max(1, len(pos) - 1), 1])
    im = ax.imshow(seg[:, :, wire], aspect="auto", cmap="RdBu_r", vmin=-image_uv, vmax=image_uv,
                   extent=[t_ms[0], t_ms[-1], len(seg) + 0.5, 0.5], interpolation="nearest")
    ax.axvline(0, color="k", lw=0.5)
    ax.axvline(pulse / FS * 1000, color="k", lw=0.5)
    for q in np.flatnonzero(quiet):
        ax.plot(t_ms[0] + 2, q + 1, "k>", ms=3)
    ax.set_ylabel("pulse in block (arrow: quietest third by the 50 ms before)")
    ax.set_title(f"{ids[wire]}: every pulse of the block (colour = uV, clipped at {image_uv:.0f})", fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    ax = fig.add_subplot(gs[max(1, len(pos) - 1):, 1])
    for w in range(4):
        ax.plot(t_ms, mean[:, w], lw=0.8, label=ids[w], color="k" if w == wire else None)
    ax.axvspan(0, pulse / FS * 1000, color="c", alpha=0.15, lw=0)
    ax.set_ylabel("pulse-triggered mean, uV")
    ax.set_xlabel("ms from pulse onset")
    ax.legend(fontsize=7, ncol=4, loc="lower right")
    ax.set_title(f"troughs below {thr_label} on {ids[wire]}, quietest third of pulses: median {np.median(rate_in[quiet]):.0f} /s "
                 f"in the pulse, {np.median(rate_pre[quiet]):.0f} /s in the 50 ms before; all pulses {np.median(rate_in):.0f} and "
                 f"{np.median(rate_pre):.0f}", fontsize=9)
    fig.suptitle(f"{experiment} {tetrode} block {blk.ttl_train} at {blk.onset_h:.2f} h ({blk.onset_local}): "
                 f"{len(seg)} pulses, {ref_label}\nMAD of each wire in the second before the block: {mad.tolist()} uV",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(figdir / f"{experiment}_{tetrode}_b{blk.ttl_train}{suffix}.png", dpi=110)
    plt.close(fig)
    return {"block": int(blk.ttl_train), "onset_h": float(blk.onset_h), "onset_local": blk.onset_local, "reference": reference,
            "n_pulses": int(len(seg)), "wire": ids[wire], "wire_index": int(wire),
            "mad_uV": mad.tolist(), "thresh_uV": round(float(thr), 1),
            "troughs_in_pulse_per_s_median": round(float(np.median(rate_in)), 1),
            "troughs_pre_per_s_median": round(float(np.median(rate_pre)), 1),
            "quiet_troughs_in_pulse_per_s_median": round(float(np.median(rate_in[quiet])), 1),
            "quiet_troughs_pre_per_s_median": round(float(np.median(rate_pre[quiet])), 1),
            "quiet_pulses_with_in_pulse_troughs_frac": round(float((n_in[quiet] > 0).mean()), 3),
            "pulses_with_in_pulse_troughs_frac": round(float((n_in > 0).mean()), 3),
            "mean_trough_1_6ms_uV": [round(float(v), 1) for v in mean[pre:pre + 180].min(0)]}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--tetrode", default="TT6")
    p.add_argument("--blocks", required=True, help="comma-separated ttl_train numbers")
    p.add_argument("--pulses", default="1,10,30,-1", help="1-based pulse positions to draw; negative from the end")
    p.add_argument("--wire", type=int, default=None, help="wire index 0 to 3 for the image; default deepest in the first block")
    p.add_argument("--thresh-mad", type=float, default=4.0, help="trough threshold in MADs of the second before the block")
    p.add_argument("--thresh-uv", type=float, default=None, help="fixed trough threshold in uV instead, the same for every block")
    p.add_argument("--scale-uv", type=float, default=100.0, help="uV between stacked wires in the trace panels")
    p.add_argument("--image-uv", type=float, default=60.0, help="colour clip of the all-pulses image")
    p.add_argument("--reference", choices=["none", "other-tetrodes"], default="none")
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    exp = experiment_from_args(a)
    figdir = exp.figures_dir / "opto_traces"
    figdir.mkdir(parents=True, exist_ok=True)
    blocks = load_blocks(exp.session)
    want = [int(s) for s in a.blocks.split(",")]
    blocks = blocks[blocks.ttl_train.isin(want)].reset_index(drop=True)
    pulses = load_pulses(exp.session)
    examples = [int(s) for s in a.pulses.split(",")]
    group = tetrode_index(a.tetrode)
    store = str(exp.store_path())
    bad = tuple(exp.badChannels)
    wire = a.wire
    if wire is None:
        first = one_block(store, blocks.iloc[0], pulses, group, a.tetrode, examples, None, a.thresh_mad,
                          a.thresh_uv, a.scale_uv, a.image_uv, figdir, exp.experiment, a.reference, bad)
        wire = first["wire_index"]
        rest = blocks.iloc[1:]
        rows = [first]
    else:
        rest = blocks
        rows = []
    rows += Parallel(n_jobs=a.workers)(delayed(one_block)(store, blk, pulses, group, a.tetrode, examples, wire,
                                                          a.thresh_mad, a.thresh_uv, a.scale_uv, a.image_uv, figdir,
                                                          exp.experiment, a.reference, bad)
                                       for _, blk in rest.iterrows())
    df = pd.DataFrame(rows).sort_values("block")
    out = RESULTS / f"opto_traces_{exp.experiment}_{a.tetrode}{'_ref' if a.reference == 'other-tetrodes' else ''}.csv"
    df.to_csv(out, index=False)
    print(df[["block", "onset_h", "wire", "thresh_uV", "quiet_troughs_pre_per_s_median", "quiet_troughs_in_pulse_per_s_median",
              "quiet_pulses_with_in_pulse_troughs_frac", "troughs_in_pulse_per_s_median"]].to_string(index=False))
    print("figures:", figdir, "table:", out)


if __name__ == "__main__":
    main()
