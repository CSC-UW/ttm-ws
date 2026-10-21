"""Pulse-triggered look at one TTM-003 stimulation session.

For every laser train: the raw (unfiltered) pulse-triggered average per
channel, the same after 300-6000 Hz filtering (where a photoelectric
artifact survives averaging and jittered spikes do not), multi-unit
threshold-crossing PSTHs per tetrode, and the per-tetrode LFP evoked
potential. Writes figures to ``<experiment>/figures/opto_*.png`` on the
data volume and a metrics JSON to ``subjects/ttm-003/results/``.

Usage: python subjects/ttm-003/scripts/opto_explore.py --subject TTM-003 --experiment TTM-003-OPTO1
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
from joblib import Parallel, delayed
import spikeinterface as si

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000


def load_pulses(session: str) -> pd.DataFrame:
    df = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    return df[~df.glitch].reset_index(drop=True)


def load_trains(session: str) -> pd.DataFrame:
    df = pd.read_csv(RESULTS / "stimulation_trains.csv")
    return df[df.session == session].reset_index(drop=True)


def load_traces(rec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Full session as float32 microvolts (samples x channels), channel groups, times."""
    x = rec.get_traces(return_in_uV=True).astype(np.float32)
    groups = np.asarray(rec.get_property("group")).astype(int)
    t = rec.get_times()
    return x, groups, t - t[0]


def _sosfilt_channels(x: np.ndarray, sos: np.ndarray) -> np.ndarray:
    cols = Parallel(n_jobs=16, prefer="threads")(delayed(sps.sosfiltfilt)(sos, x[:, c]) for c in range(x.shape[1]))
    return np.stack(cols, axis=1).astype(np.float32)


def bandpass(x: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return _sosfilt_channels(x, sps.butter(3, [lo, hi], btype="bandpass", fs=FS, output="sos"))


def lowpass(x: np.ndarray, hi: float) -> np.ndarray:
    return _sosfilt_channels(x, sps.butter(4, hi, btype="lowpass", fs=FS, output="sos"))


def snippets(x: np.ndarray, t: np.ndarray, events_s: np.ndarray, pre_s: float, post_s: float) -> np.ndarray:
    """(events x samples x channels) around event times."""
    idx = np.searchsorted(t, events_s)
    pre, post = int(pre_s * FS), int(post_s * FS)
    keep = (idx - pre >= 0) & (idx + post <= x.shape[0])
    return np.stack([x[i - pre:i + post] for i in idx[keep]])


def detect_mua(xf: np.ndarray, groups: np.ndarray, thresh_mad: float = 5.0, dead_ms: float = 1.0) -> dict[int, np.ndarray]:
    """Negative threshold crossings per tetrode (sample indices), merged within dead_ms."""
    out = {}
    dead = int(dead_ms * FS / 1000)
    for g in np.unique(groups):
        chans = np.where(groups == g)[0]
        idx = []
        for c in chans:
            v = xf[:, c]
            mad = np.median(np.abs(v)) / 0.6745
            below = v < -thresh_mad * mad
            on = np.where(below[1:] & ~below[:-1])[0] + 1
            idx.append(on)
        idx = np.sort(np.concatenate(idx))
        if idx.size:
            keep = np.r_[True, np.diff(idx) > dead]
            idx = idx[keep]
        out[int(g)] = idx
    return out


def psth(spike_s: np.ndarray, events_s: np.ndarray, pre: float, post: float, bin_s: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(-pre, post + bin_s / 2, bin_s)
    counts = np.zeros(len(edges) - 1)
    for e in events_s:
        rel = spike_s[(spike_s >= e - pre) & (spike_s < e + post)] - e
        counts += np.histogram(rel, edges)[0]
    rate = counts / (len(events_s) * bin_s)
    return edges[:-1] + bin_s / 2, rate


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    a = p.parse_args()
    exp = experiment_from_args(a)
    figdir = exp.figures_dir
    figdir.mkdir(parents=True, exist_ok=True)
    pulses = load_pulses(exp.session)
    trains = load_trains(exp.session)

    rec = si.read_zarr(str(exp.store_path()))
    x, groups, t = load_traces(rec)
    ch_names = list(rec.get_channel_ids())
    tetrodes = sorted(set(groups.tolist()))
    print(f"loaded {x.shape} ({x.nbytes / 1e9:.1f} GB), {len(pulses)} pulses, {len(trains)} trains")

    xf = bandpass(x, 300, 6000)
    xl = lowpass(x, 300)
    mua = detect_mua(xf, groups)
    mua_s = {g: t[i] for g, i in mua.items()}
    np.savez(exp.session_dir / "mua_crossings_5mad.npz", **{f"TT{g + 1}": v for g, v in mua_s.items()})
    print("MUA rates (Hz):", {g: round(len(v) / t[-1], 1) for g, v in mua_s.items()})

    metrics = {"experiment": exp.experiment, "session": exp.session, "trains": []}
    for tr in trains.itertuples():
        pp = pulses[pulses.train == tr.ttl_train]
        on = pp.onset_s.values
        off = pp.offset_s.values
        pulse_s = tr.pulse_ms_median / 1000
        pre, post = 0.05, pulse_s + 0.1
        label = f"train{tr.n}_{tr.ttl_pattern}_{tr.power_mW:g}mW"

        raw_on = snippets(x, t, on, pre, post).mean(0)
        raw_off = snippets(x, t, off, 0.02, 0.05).mean(0)
        hf_on = snippets(xf, t, on, 0.01, 0.02).mean(0)
        hf_off = snippets(xf, t, off, 0.01, 0.02).mean(0)
        lfp_on = snippets(xl, t, on, pre, post).mean(0)
        tt = np.arange(raw_on.shape[0]) / FS - pre

        base = raw_on[tt < -0.005]
        # Onset step within the first 2 ms, and the fast-transient amplitude in the
        # high-passed trial average (artifact signature) per channel.
        step = raw_on[(tt >= 0) & (tt < 0.002)].mean(0) - base.mean(0)
        hf_t = np.arange(hf_on.shape[0]) / FS - 0.01
        hf_noise = hf_on[hf_t < -0.002].std(0)
        hf_peak_on = np.abs(hf_on[(hf_t >= -0.0005) & (hf_t < 0.002)]).max(0)
        hf_peak_off = np.abs(hf_off[(hf_t >= -0.0005) & (hf_t < 0.002)]).max(0)
        sustained = raw_on[(tt >= 0.02) & (tt < pulse_s)].mean(0) - base.mean(0)

        # MUA per tetrode: baseline (-50..0 ms) vs early (0..10 ms) vs during pulse.
        mua_rows = {}
        for g in tetrodes:
            sp = mua_s[g]
            tb, r = psth(sp, on, pre, post, 0.002)
            mua_rows[f"TT{g + 1}"] = {
                "baseline_hz": float(r[tb < 0].mean()),
                "early_0_10ms_hz": float(r[(tb >= 0) & (tb < 0.010)].mean()),
                "pulse_10ms_end_hz": float(r[(tb >= 0.010) & (tb < pulse_s)].mean()),
                "post_0_50ms_hz": float(r[(tb >= pulse_s) & (tb < pulse_s + 0.05)].mean()),
            }
        # Whole-train MUA time course: rate in the 60 s before, during, and 60 s after.
        t0, t1 = on[0], off[-1]
        course = {}
        for g in tetrodes:
            sp = mua_s[g]
            course[f"TT{g + 1}"] = {
                "pre60_hz": float(((sp > t0 - 60) & (sp < t0)).sum() / 60),
                "during_hz": float(((sp > t0) & (sp < t1)).sum() / (t1 - t0)),
                "post60_hz": float(((sp > t1) & (sp < t1 + 60)).sum() / 60),
            }

        metrics["trains"].append({
            "n": int(tr.n), "pattern": tr.ttl_pattern, "power_mW": float(tr.power_mW),
            "state": tr.state, "n_pulses": int(len(on)),
            "onset_step_uV_by_channel": dict(zip(ch_names, np.round(step, 2).tolist())),
            "sustained_shift_uV_by_channel": dict(zip(ch_names, np.round(sustained, 2).tolist())),
            "hf_avg_onset_peak_uV_by_channel": dict(zip(ch_names, np.round(hf_peak_on, 2).tolist())),
            "hf_avg_offset_peak_uV_by_channel": dict(zip(ch_names, np.round(hf_peak_off, 2).tolist())),
            "hf_avg_noise_uV_median": float(np.median(hf_noise)),
            "mua_by_tetrode": mua_rows, "mua_course_by_tetrode": course,
        })

        # --- figures -------------------------------------------------------
        fig, axes = plt.subplots(2, 2, figsize=(16, 10))
        ax = axes[0, 0]
        im = ax.imshow(raw_on.T, aspect="auto", extent=[tt[0] * 1000, tt[-1] * 1000, 63.5, -0.5],
                       cmap="RdBu_r", vmin=-np.percentile(np.abs(raw_on), 99), vmax=np.percentile(np.abs(raw_on), 99))
        ax.axvline(0, color="k", lw=0.5); ax.axvline(pulse_s * 1000, color="k", lw=0.5)
        ax.set_title(f"raw pulse-triggered average (uV), n={len(on)}"); ax.set_xlabel("ms from onset"); ax.set_ylabel("channel (tetrode order)")
        for g in tetrodes[1:]:
            ax.axhline(np.where(groups == g)[0][0] - 0.5, color="k", lw=0.3)
        plt.colorbar(im, ax=ax)

        ax = axes[0, 1]
        for g in tetrodes:
            ax.plot(tt * 1000, lfp_on[:, groups == g].mean(1) + 0 , lw=0.8, label=f"TT{g + 1}")
        ax.axvspan(0, pulse_s * 1000, color="c", alpha=0.15)
        ax.set_title("LFP (<300 Hz) evoked potential, mean per tetrode"); ax.set_xlabel("ms"); ax.set_ylabel("uV")
        ax.legend(ncol=4, fontsize=7)

        ax = axes[1, 0]
        ax.plot(hf_t * 1000, hf_on, lw=0.5, alpha=0.6)
        ax.set_title("300-6000 Hz trial average at onset (all channels)"); ax.set_xlabel("ms from onset"); ax.set_ylabel("uV")
        ax2 = ax.inset_axes([0.6, 0.6, 0.38, 0.38])
        ax2.plot(hf_t * 1000, hf_off, lw=0.5, alpha=0.6); ax2.set_title("at offset", fontsize=8)

        ax = axes[1, 1]
        for g in tetrodes:
            tb, r = psth(mua_s[g], on, pre, post, 0.002)
            ax.plot(tb * 1000, r, lw=0.8, label=f"TT{g + 1}")
        ax.axvspan(0, pulse_s * 1000, color="c", alpha=0.15)
        ax.set_title("MUA PSTH per tetrode (2 ms bins)"); ax.set_xlabel("ms from onset"); ax.set_ylabel("crossings/s")
        ax.legend(ncol=4, fontsize=7)
        fig.suptitle(f"{exp.experiment} {label} ({tr.state}, {tr.n_pulses} pulses of {tr.pulse_ms_median:.0f} ms at {tr.rate_hz:.2f} Hz)")
        fig.tight_layout()
        fig.savefig(figdir / f"opto_{label}.png", dpi=110)
        plt.close(fig)

    # Session overview: MUA rate per tetrode over time with trains shaded.
    fig, ax = plt.subplots(figsize=(18, 8))
    edges = np.arange(0, t[-1] + 5, 5)
    for k, g in enumerate(tetrodes):
        r = np.histogram(mua_s[g], edges)[0] / 5
        ax.plot(edges[:-1] / 60, r / max(r.max(), 1e-9) + k, lw=0.6)
        ax.text(-0.5, k + 0.5, f"TT{g + 1}", fontsize=7, ha="right")
    for tr in trains.itertuples():
        ax.axvspan(tr.onset_s / 60, tr.offset_s / 60, color="c", alpha=0.2)
        ax.text(tr.onset_s / 60, len(tetrodes) + 0.2, f"{tr.n}:{tr.ttl_pattern} {tr.power_mW:g}mW", fontsize=7, rotation=90)
    ax.set_xlabel("min"); ax.set_ylabel("MUA rate per tetrode (normalised, 5 s bins)")
    ax.set_title(f"{exp.experiment}: MUA over the session; shaded = laser trains")
    fig.tight_layout(); fig.savefig(figdir / "opto_session_mua.png", dpi=110); plt.close(fig)

    out = RESULTS / f"opto_explore_{exp.experiment}.json"
    out.write_text(json.dumps(metrics, indent=1) + "\n")
    print(f"wrote {out} and figures in {figdir}")


if __name__ == "__main__":
    main()
