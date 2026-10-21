"""Spike-band trial averages and single trials around laser onset.

A photoelectric artifact is locked to the light edge (0 ms, and again at
offset) and scales with power; light-evoked spikes start after ~1 ms with
trial-to-trial jitter. This plots, for the lowest- and highest-power
trains, the 300-6000 Hz pulse-triggered average per tetrode from -5 to
+12 ms around onset and offset, and 20 single trials on the channel with
the largest onset response. Reads pulse windows straight from the store.

Usage: python subjects/ttm-003/scripts/opto_onset_zoom.py --subject TTM-003 --experiment TTM-003-OPTO1
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

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
PAD = 3000
PRE_MS, POST_MS = 5, 12


def snippets(rec, t, events_s, pre, post, sos, max_n=200):
    idx = np.searchsorted(t, events_s)[:max_n]
    out = np.stack([rec.get_traces(start_frame=i - pre - PAD, end_frame=i + post + PAD, return_in_uV=True) for i in idx])
    return sps.sosfiltfilt(sos, out, axis=1)[:, PAD:-PAD]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    a = p.parse_args()
    exp = experiment_from_args(a)
    figdir = exp.figures_dir
    figdir.mkdir(parents=True, exist_ok=True)
    rec = si.read_zarr(str(exp.store_path()))
    t = rec.get_times(); t = t - t[0]
    groups = np.asarray(rec.get_property("group")); ids = np.asarray(rec.channel_ids)
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv"); pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv"); trains = trains[trains.session == exp.session]
    sos = sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos")
    pre, post = PRE_MS * 30, POST_MS * 30
    tt = (np.arange(pre + post) - pre) / 30
    pick = [trains.loc[trains.power_mW.idxmin()], trains.loc[trains.power_mW.idxmax()]]
    tetrodes = sorted(set(groups.tolist()))
    summary = {}
    for tr in pick:
        pp = pulses[pulses.train == tr.ttl_train]
        avg_on = snippets(rec, t, pp.onset_s.values, pre, post, sos).mean(0)
        avg_off = snippets(rec, t, pp.offset_s.values, pre, post, sos).mean(0)
        noise = avg_on[tt < -1].std(0)
        edge_on = np.abs(avg_on[(tt >= -0.2) & (tt < 0.7)]).max(0)
        late_on = np.abs(avg_on[(tt >= 1) & (tt < 6)]).max(0)
        edge_off = np.abs(avg_off[(tt >= -0.2) & (tt < 0.7)]).max(0)
        peak_lat = tt[(tt >= 1) & (tt < 6)][np.abs(avg_on[(tt >= 1) & (tt < 6)]).argmax(0)]
        summary[f"train{tr.n}_{tr.power_mW:g}mW"] = {
            "noise_uV_median": float(np.median(noise)),
            "onset_edge_0ms_max_uV": float(edge_on.max()), "onset_edge_0ms_median_uV": float(np.median(edge_on)),
            "offset_edge_0ms_max_uV": float(edge_off.max()), "offset_edge_0ms_median_uV": float(np.median(edge_off)),
            "evoked_1_6ms_max_uV": float(late_on.max()), "evoked_1_6ms_by_tetrode_uV": {f"TT{g + 1}": float(late_on[groups == g].max()) for g in tetrodes},
            "evoked_peak_latency_ms_by_tetrode": {f"TT{g + 1}": float(peak_lat[groups == g][late_on[groups == g].argmax()]) for g in tetrodes},
        }
        fig, axes = plt.subplots(4, 8, figsize=(26, 12))
        for k, g in enumerate(tetrodes):
            ax = axes[(k // 8) * 2, k % 8]
            for i in np.where(groups == g)[0]:
                ax.plot(tt, avg_on[:, i], lw=0.8, label=str(ids[i]))
            ax.axvline(0, color="k", lw=0.5); ax.set_title(f"TT{g + 1} onset", fontsize=9); ax.legend(fontsize=5)
            ax = axes[(k // 8) * 2 + 1, k % 8]
            for i in np.where(groups == g)[0]:
                ax.plot(tt, avg_off[:, i], lw=0.8)
            ax.axvline(0, color="k", lw=0.5); ax.set_title(f"TT{g + 1} offset", fontsize=9)
        fig.suptitle(f"{exp.experiment} train {tr.n} ({tr.ttl_pattern}, {tr.power_mW:g} mW, n={min(len(pp), 200)}): 300-6000 Hz pulse-triggered average (uV) at onset and offset")
        fig.tight_layout(); fig.savefig(figdir / f"opto_onset_zoom_train{tr.n}_{tr.power_mW:g}mW.png", dpi=90); plt.close(fig)
        # single trials on the most responsive channel
        ci = int(late_on.argmax())
        fig, ax = plt.subplots(figsize=(12, 5))
        for i in np.searchsorted(t, pp.onset_s.values)[50:70]:
            s = rec.get_traces(start_frame=i - 60 - PAD, end_frame=i + 300 + PAD, channel_ids=[ids[ci]], return_in_uV=True)[:, 0]
            ax.plot((np.arange(360) - 60) / 30, sps.sosfiltfilt(sos, s)[PAD:-PAD], lw=0.6, alpha=0.7)
        ax.axvline(0, color="k"); ax.set_xlabel("ms from onset"); ax.set_ylabel("uV")
        ax.set_title(f"{exp.experiment} train {tr.n} ({tr.power_mW:g} mW): {ids[ci]} (TT{groups[ci] + 1}) 20 single trials, 300-6000 Hz")
        fig.tight_layout(); fig.savefig(figdir / f"opto_onset_single_trials_train{tr.n}_{tr.power_mW:g}mW.png", dpi=90); plt.close(fig)
    out = RESULTS / f"opto_onset_zoom_{exp.experiment}.json"
    out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if not isinstance(vv, dict)} for k, v in summary.items()}, indent=1))


if __name__ == "__main__":
    main()
