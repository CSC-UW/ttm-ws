"""Are the small and the large spikes of a merged cluster one cell?

Takes a cluster of a re-sort (``opto_tt_resort.py``), measures every
spike's trough on the dominant wire, splits the cluster into large and
small spikes by amplitude, and tests whether the two populations respect
a common refractory period: their cross-correlogram outside the laser
trains (0.5 ms bins) is compared with the rate expected from
independence. One cell cannot fire a small spike within ~2 ms of a large
one; two cells can. Also reports the ISI floor of the small spikes alone
and their amplitude distribution inside vs outside the pulses. Writes
``results/ampsplit_<experiment>_<TT>_thr<thr>_c<cluster>.json`` and a figure.

Usage: python subjects/ttm-003/scripts/opto_tt_amplitude_split.py --subject TTM-003 --experiment TTM-003-OPTO1 --threshold 3.5 --cluster 10 --large -85 --small -70
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
import spikeinterface as si

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000


def ccg_counts(a, b, window_s, bin_s):
    nb = int(round(window_s / bin_s)); edges = (np.arange(-nb, nb + 2) - 0.5) * bin_s
    lags = []
    i0 = np.searchsorted(b, a - window_s - bin_s); i1 = np.searchsorted(b, a + window_s + bin_s)
    for ta, s, e in zip(a, i0, i1):
        if e > s: lags.append(b[s:e] - ta)
    lags = np.concatenate(lags) if lags else np.array([])
    return edges[:-1] + bin_s / 2, np.histogram(lags, edges)[0]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis")
    p.add_argument("--tetrode", default="TT6")
    p.add_argument("--threshold", type=float, required=True)
    p.add_argument("--cluster", type=int, required=True)
    p.add_argument("--large", type=float, default=-85.0, help="uV: spikes at or below this trough are 'large'")
    p.add_argument("--small", type=float, default=-70.0, help="uV: spikes at or above this trough are 'small'")
    a = p.parse_args()
    exp = experiment_from_args(a)
    out = exp.sorting_dir(a.variant, subdir="sortings_opto")
    rec = si.load(str(out / "preprocessed"))
    t = rec.get_times(); t = t - t[0]
    ids = np.asarray(rec.channel_ids); groups = np.asarray(rec.get_property("group")).astype(int)
    g = int(a.tetrode.removeprefix("TT")) - 1; chans = np.flatnonzero(groups == g)
    gain = float(rec.get_property("gain_to_uV")[chans[0]])
    sorting = si.load(str(out / f"resort_{a.tetrode}_thr{a.threshold:g}" / "sorting"))
    sp = sorting.get_unit_spike_train(a.cluster)
    sp = sp[(sp > 30) & (sp < rec.get_num_frames() - 30)]
    # trough on the dominant wire for every spike (read in blocks)
    x = rec.get_traces(channel_ids=ids[chans])
    dom_tmpl = np.stack([x[f - 3:f + 4].min(0) for f in sp[:2000]]).mean(0)
    dom = int(np.argmin(dom_tmpl))
    amp = np.array([x[f - 3:f + 4, dom].min() for f in sp]) * gain
    del x

    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv"); pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv"); trains = trains[trains.session == exp.session]
    ts = t[sp]
    outside = np.ones(len(sp), bool)
    for tr in trains.itertuples():
        outside &= ~((ts > tr.onset_s - 1) & (ts < tr.offset_s + 1))
    inpulse = np.zeros(len(sp), bool)
    for a_, b_ in zip(pulses.onset_s.values, pulses.offset_s.values):
        inpulse |= (ts > a_ + 0.01) & (ts < b_)
    large = amp <= a.large; small = amp >= a.small
    L = ts[outside & large]; Sm = ts[outside & small]
    T_out = (outside.sum() / len(sp)) * t[-1]
    lags, c = ccg_counts(L, Sm, 0.02, 0.0005)
    expected = len(L) * len(Sm) * 0.0005 / T_out
    near = np.abs(lags) <= 0.0015
    obs_near = int(c[near].sum()); exp_near = float(expected * near.sum())
    isi_small = np.diff(Sm) * 1000; isi_large = np.diff(L) * 1000
    isi_small_in = np.diff(ts[inpulse & small]) * 1000; isi_small_in = isi_small_in[isi_small_in < 50]
    core = np.abs(lags) <= 0.00075
    res = {"experiment": exp.experiment, "tetrode": a.tetrode, "threshold": a.threshold, "cluster": a.cluster,
           "ccg_large_small_within_0.75ms": {"observed": int(c[core].sum()), "expected_if_independent": round(float(expected * core.sum()), 1),
                                             "ratio": round(float(c[core].sum() / (expected * core.sum())), 3)},
           "ccg_large_small_peak_2_3ms_ratio": round(float(c[(lags > 0.00175) & (lags < 0.00325)].mean() / expected), 2),
           "dominant_wire": str(ids[chans[dom]]), "n_spikes": int(len(sp)),
           "amp_uV_quantiles_outside": [round(float(v), 1) for v in np.percentile(amp[outside], [5, 25, 50, 75, 95])],
           "amp_uV_quantiles_in_pulse": [round(float(v), 1) for v in np.percentile(amp[inpulse], [5, 25, 50, 75, 95])] if inpulse.any() else None,
           "n_large_outside": int(len(L)), "n_small_outside": int(len(Sm)), "seconds_outside": round(float(T_out), 1),
           "ccg_large_small_within_1.5ms": {"observed": obs_near, "expected_if_independent": round(exp_near, 1),
                                            "ratio": round(obs_near / exp_near, 3) if exp_near else None},
           "ccg_large_small_5_20ms_ratio": round(float(c[np.abs(lags) > 0.005].sum() / (expected * (np.abs(lags) > 0.005).sum())), 3),
           "isi_small_outside_ms_pct_1_5_50": [round(float(v), 2) for v in np.percentile(isi_small, [1, 5, 50])],
           "isi_large_outside_ms_pct_1_5_50": [round(float(v), 2) for v in np.percentile(isi_large, [1, 5, 50])],
           "isi_small_in_pulse_ms_pct_1_5_50": [round(float(v), 2) for v in np.percentile(isi_small_in, [1, 5, 50])] if isi_small_in.size > 10 else None,
           "frac_small_isi_lt_2ms_outside": round(float((isi_small < 2).mean()), 4),
           "frac_large_isi_lt_2ms_outside": round(float((isi_large < 2).mean()), 4)}
    (RESULTS / f"ampsplit_{exp.experiment}_{a.tetrode}_thr{a.threshold:g}_c{a.cluster}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    ax = axes[0]
    ax.hist(amp[outside], bins=np.arange(-250, -20, 2.5), histtype="step", color="k", density=True, label=f"outside trains (n={outside.sum()})")
    if inpulse.any(): ax.hist(amp[inpulse], bins=np.arange(-250, -20, 2.5), histtype="step", color="tab:green", density=True, label=f"in pulse (n={inpulse.sum()})")
    ax.axvline(a.large, color="tab:red", lw=0.8); ax.axvline(a.small, color="tab:blue", lw=0.8); ax.legend(fontsize=8)
    ax.set_title(f"cluster {a.cluster} trough on {ids[chans[dom]]} (uV)")
    ax = axes[1]
    ax.bar(lags * 1000, c, width=0.5, color="k"); ax.axhline(expected, color="gray", lw=0.8)
    ax.set_title(f"CCG large -> small, outside trains (0.5 ms bins); within +/-1.5 ms: {obs_near} observed vs {exp_near:.0f} expected"); ax.set_xlabel("ms (small relative to large)")
    ax = axes[2]
    ax.hist(isi_small[isi_small < 50], bins=np.arange(0, 50, 0.5), histtype="step", color="tab:blue", density=True, label="small, outside")
    ax.hist(isi_large[isi_large < 50], bins=np.arange(0, 50, 0.5), histtype="step", color="tab:red", density=True, label="large, outside")
    if isi_small_in.size > 10: ax.hist(isi_small_in, bins=np.arange(0, 50, 0.5), histtype="step", color="tab:green", density=True, label="small, in pulse")
    ax.legend(fontsize=8); ax.set_title("ISI (ms)")
    fig.suptitle(f"{exp.experiment} {a.tetrode} re-sort thr {a.threshold:g}, cluster {a.cluster}: large (<= {a.large:g} uV) vs small (>= {a.small:g} uV) spikes")
    fig.tight_layout(); fig.savefig(exp.figures_dir / f"opto_{a.tetrode}_ampsplit_thr{a.threshold:g}_c{a.cluster}.png", dpi=100); plt.close(fig)


if __name__ == "__main__":
    main()
