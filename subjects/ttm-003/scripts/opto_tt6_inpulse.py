"""Do the TT6 cell's spikes continue through the pulse at reduced amplitude?

For one tetrode and one reference cluster, detects every negative crossing
(default 3.5 MAD, any wire) on the preprocessed recording around the pulses
of the highest-power trains, and for each event records the time in the
pulse, the trough amplitude on each wire, and the cosine similarity of its
four-wire amplitude vector to the reference template's (the tetrode
"fingerprint" that is independent of overall amplitude). Events are
labelled sorted (claimed by any cluster of the tetrode) or unassigned.
Writes ``results/tt_inpulse_<experiment>_<tetrode>.json`` and
``figures/opto_<tetrode>_inpulse_<variant>.png``.

Usage: python subjects/ttm-003/scripts/opto_tt6_inpulse.py --subject TTM-003 --experiment TTM-003-OPTO1 --unit 21
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


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis")
    p.add_argument("--tetrode", default="TT6")
    p.add_argument("--unit", type=int, required=True, help="reference cluster id on that tetrode")
    p.add_argument("--thresh-mad", type=float, default=3.5)
    a = p.parse_args()
    exp = experiment_from_args(a)
    out = exp.sorting_dir(a.variant, subdir="sortings_opto")
    rec = si.load(str(out / "preprocessed")); sorting = si.load(str(out / "sorting"))
    analyzer = si.load_sorting_analyzer(str(out / "analyzer"))
    t = rec.get_times(); t = t - t[0]
    ids = np.asarray(rec.channel_ids); groups = np.asarray(rec.get_property("group")).astype(int)
    g = int(a.tetrode.removeprefix("TT")) - 1
    chans = np.flatnonzero(groups == g)
    gain = float(rec.get_property("gain_to_uV")[chans[0]])
    ug = pd.Series(sorting.get_property("group"), index=sorting.unit_ids).astype(int)
    tet_units = [u for u in sorting.unit_ids if ug[u] == g]
    sorted_frames = np.sort(np.concatenate([sorting.get_unit_spike_train(u) for u in tet_units]))
    ref_frames = sorting.get_unit_spike_train(a.unit)
    tmpl = analyzer.get_extension("templates").get_unit_template(a.unit, operator="average")[:, chans]
    ref_vec = tmpl.min(axis=0); ref_vec = ref_vec / np.linalg.norm(ref_vec)

    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv"); pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv"); trains = trains[trains.session == exp.session]
    pmax = trains.power_mW.max()
    sel = pulses[pulses.train.isin(trains[trains.power_mW == pmax].ttl_train)]
    dur = float(np.median(sel.offset_s - sel.onset_s))

    ref = rec.get_traces(start_frame=0, end_frame=FS * 120, channel_ids=ids[chans])
    mad = np.median(np.abs(ref), axis=0) / 0.6745
    rows = []; wf_in, wf_base, isis = [], [], []
    for on in sel.onset_s.values:
        ev_in = []
        i0 = np.searchsorted(t, on - 0.3); i1 = np.searchsorted(t, on + dur + 0.1)
        x = rec.get_traces(start_frame=i0, end_frame=i1, channel_ids=ids[chans])
        below = (x < -a.thresh_mad * mad).any(1)
        starts = np.flatnonzero(below[1:] & ~below[:-1]) + 1
        starts = starts[np.r_[True, np.diff(starts) > 15]]
        for s in starts:
            if s < 10 or s > len(x) - 20: continue
            seg = x[s - 10:s + 20]
            k = np.unravel_index(np.argmin(seg), seg.shape)[0]
            troughs = seg[max(k - 3, 0):k + 4].min(axis=0) * gain
            frame = i0 + s - 10 + k
            d = np.abs(sorted_frames - frame).min() if sorted_frames.size else 1e9
            dref = np.abs(ref_frames - frame).min() if ref_frames.size else 1e9
            vec = troughs / max(np.linalg.norm(troughs), 1e-9)
            tm = (t[frame] - on) * 1000
            rows.append({"t_ms": tm, "amp_uV": float(troughs.min()),
                         "cos_ref": float(vec @ ref_vec), "sorted": bool(d <= 15), "is_ref": bool(dref <= 15),
                         **{str(ids[c]): float(v) for c, v in zip(chans, troughs)}})
            if vec @ ref_vec > 0.95 and d > 15 and 30 <= s - 10 + k and s - 10 + k + 60 <= len(x):
                w = x[s - 10 + k - 30:s - 10 + k + 60] * gain
                if 10 < tm < dur * 1000:
                    wf_in.append(w); ev_in.append(tm)
                elif tm < -0.5:
                    wf_base.append(w)
        isis += list(np.diff(sorted(ev_in)))
    df = pd.DataFrame(rows)
    df["phase"] = pd.cut(df.t_ms, [-300, -0.5, 10, dur * 1000, dur * 1000 + 100], labels=["baseline", "onset_0_10ms", "in_pulse", "post_100ms"])

    summ = {}
    for ph, d in df.groupby("phase", observed=True):
        n_ev = len(d); secs = len(sel) * {"baseline": 0.3, "onset_0_10ms": 0.0105, "in_pulse": dur - 0.01, "post_100ms": 0.1}[ph]
        look = d[d.cos_ref > 0.95]
        summ[ph] = {"events_per_s": round(n_ev / secs, 1), "sorted_per_s": round(d["sorted"].sum() / secs, 1),
                    "ref_unit_per_s": round(d.is_ref.sum() / secs, 1),
                    "lookalike_cos>0.95_per_s": round(len(look) / secs, 1),
                    "lookalike_unsorted_per_s": round((~look["sorted"]).sum() / secs, 1),
                    "lookalike_amp_uV_quantiles": [round(float(v), 1) for v in np.percentile(look.amp_uV, [10, 50, 90])] if len(look) else None,
                    "all_amp_uV_median": round(float(d.amp_uV.median()), 1)}
    wf_in = np.array(wf_in); wf_base = np.array(wf_base); isis = np.array(isis)
    def shape(w):
        return w / np.abs(w).max()
    c_ref = int(np.argmin(ref_vec))
    m_in, m_base = wf_in.mean(0), wf_base.mean(0)
    def ttp(w):
        k = int(np.argmin(w)); return (k + int(np.argmax(w[k:])) - k) / FS * 1000
    lookalike = {
        "n_in_pulse_unsorted": int(len(wf_in)), "n_baseline_unsorted": int(len(wf_base)),
        "in_pulse_isi_ms_quantiles": [round(float(v), 1) for v in np.percentile(isis, [5, 25, 50, 75, 95])] if isis.size else None,
        "in_pulse_isi_cv": round(float(isis.std() / isis.mean()), 2) if isis.size else None,
        "mean_trough_uV_by_wire": {"in_pulse": {str(ids[c]): round(float(v), 1) for c, v in zip(chans, m_in.min(0))},
                                   "baseline": {str(ids[c]): round(float(v), 1) for c, v in zip(chans, m_base.min(0))},
                                   "ref_template": {str(ids[c]): round(float(v), 1) for c, v in zip(chans, tmpl.min(0))}},
        "dominant_wire_shape_corr": {"in_pulse_vs_template": round(float(np.corrcoef(shape(m_in[:, c_ref]), shape(tmpl[:, c_ref]))[0, 1]), 3),
                                     "baseline_vs_template": round(float(np.corrcoef(shape(m_base[:, c_ref]), shape(tmpl[:, c_ref]))[0, 1]), 3)},
        "trough_to_peak_ms": {"in_pulse": round(ttp(m_in[:, c_ref]), 2), "baseline": round(ttp(m_base[:, c_ref]), 2), "template": round(ttp(tmpl[:, c_ref]), 2)},
    }
    res = {"experiment": exp.experiment, "unsorted_lookalikes": lookalike, "tetrode": a.tetrode, "ref_unit": a.unit, "power_mW": float(pmax), "n_pulses": int(len(sel)),
           "pulse_s": round(dur, 3), "thresh_mad": a.thresh_mad, "ref_fingerprint": {str(ids[c]): round(float(v), 3) for c, v in zip(chans, ref_vec)},
           "by_phase": summ}
    (RESULTS / f"tt_inpulse_{exp.experiment}_{a.tetrode}.json").write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))

    fig, axes = plt.subplots(3, 2, figsize=(16, 15))
    ax = axes[0, 0]
    for lab, d, c in [("unassigned", df[~df["sorted"]], "0.6"), ("other sorted", df[df["sorted"] & ~df.is_ref], "tab:blue"), (f"unit {a.unit}", df[df.is_ref], "tab:red")]:
        ax.scatter(d.t_ms, d.amp_uV, s=4, c=c, alpha=0.5, label=lab)
    ax.axvspan(0, dur * 1000, color="c", alpha=0.15); ax.set_xlabel("ms from onset"); ax.set_ylabel("trough amplitude (uV, min over wires)"); ax.legend(fontsize=8)
    ax.set_title(f"{a.tetrode} crossings (>{a.thresh_mad:g} MAD) around {len(sel)} pulses at {pmax:g} mW")
    ax = axes[0, 1]
    ax.scatter(df.t_ms, df.cos_ref, s=4, c=np.where(df.is_ref, "tab:red", np.where(df["sorted"], "tab:blue", "0.6")), alpha=0.5)
    ax.axvspan(0, dur * 1000, color="c", alpha=0.15); ax.axhline(0.95, color="k", lw=0.5, ls="--"); ax.set_xlabel("ms from onset"); ax.set_ylabel(f"cosine to unit {a.unit} four-wire fingerprint")
    ax = axes[1, 0]
    for ph, c in [("baseline", "k"), ("onset_0_10ms", "tab:orange"), ("in_pulse", "tab:green"), ("post_100ms", "tab:purple")]:
        d = df[df.phase == ph]
        ax.hist(d.cos_ref, bins=np.linspace(-1, 1, 81), histtype="step", density=True, color=c, label=f"{ph} (n={len(d)})")
    ax.legend(fontsize=8); ax.set_xlabel("cosine to fingerprint"); ax.set_title("fingerprint similarity by phase")
    ax = axes[1, 1]
    look = df[df.cos_ref > 0.95]
    for ph, c in [("baseline", "k"), ("in_pulse", "tab:green")]:
        d = look[look.phase == ph]
        ax.hist(d.amp_uV, bins=np.arange(-200, -20, 5), histtype="step", density=True, color=c, label=f"{ph} look-alikes (n={len(d)})")
    ax.hist(df[df.is_ref & (df.phase == "baseline")].amp_uV, bins=np.arange(-200, -20, 5), histtype="step", density=True, color="tab:red", label=f"unit {a.unit} baseline spikes")
    ax.legend(fontsize=8); ax.set_xlabel("trough amplitude (uV)"); ax.set_title("amplitude of fingerprint look-alikes (cos > 0.95)")
    tt = (np.arange(90) - 30) / FS * 1000
    ax = axes[2, 0]
    for c in range(len(chans)):
        ax.plot(tt, m_in[:, c], label=f"{ids[chans[c]]} in-pulse"); ax.plot(tt, m_base[:, c], ls="--", color=ax.lines[-1].get_color(), label=f"{ids[chans[c]]} baseline")
    ax.legend(fontsize=6, ncol=2); ax.set_title(f"mean waveform of unsorted look-alikes: in-pulse (n={len(wf_in)}, solid) vs baseline (n={len(wf_base)}, dashed)"); ax.set_xlabel("ms")
    ax = axes[2, 1]
    for c in range(len(chans)):
        ax.plot(tt, tmpl[:, c], label=str(ids[chans[c]]))
    ax.legend(fontsize=7); ax.set_title(f"unit {a.unit} template; in-pulse look-alike ISI median {np.median(isis):.1f} ms, CV {isis.std() / isis.mean():.2f}" if isis.size else f"unit {a.unit} template")
    fig.suptitle(f"{exp.experiment}: does unit {a.unit} ({a.tetrode}) keep firing during the pulse at lower amplitude?")
    fig.tight_layout(); fig.savefig(exp.figures_dir / f"opto_{a.tetrode}_inpulse_{a.variant}.png", dpi=100); plt.close(fig)


if __name__ == "__main__":
    main()
