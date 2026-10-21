"""Re-sort one tetrode at a lower MountainSort5 threshold and compare.

Runs MountainSort5 scheme 2 on the four wires of ``--tetrode`` taken from
the session's preprocessed (bandpass + CMR) recording, at each
``--thresholds`` value, under
``<session>/sortings_opto/<variant>/resort_<TT>_thr<thr>/``. For every new
cluster it reports spikes, rate, amplitude, ISI violations, the four-wire
fingerprint cosine to a reference cluster of the main sort (``--unit``),
the fraction of the reference cluster's spikes it reproduces (within
0.5 ms), and its firing around the highest-power pulses (baseline, onset
1 to 10 ms, rest of pulse, 100 ms after). Writes
``results/resort_<experiment>_<TT>.csv`` and a figure per threshold.

Usage: python subjects/ttm-003/scripts/opto_tt_resort.py --subject TTM-003 --experiment TTM-003-OPTO1 --unit 21 --thresholds 4.0 3.5
"""

from __future__ import annotations

import argparse
import pathlib
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spikeinterface as si
import spikeinterface.sorters as ss

from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
MATCH = 15


def frac_matched(ref: np.ndarray, other: np.ndarray) -> float:
    if ref.size == 0 or other.size == 0:
        return 0.0
    j = np.searchsorted(other, ref)
    lo = other[np.clip(j - 1, 0, len(other) - 1)]; hi = other[np.clip(j, 0, len(other) - 1)]
    return float((np.minimum(np.abs(ref - lo), np.abs(ref - hi)) <= MATCH).mean())


def rate(frames, t, events, lo, hi):
    ts = t[frames]
    return sum(((ts >= e + lo) & (ts < e + hi)).sum() for e in events) / (len(events) * (hi - lo))


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis")
    p.add_argument("--tetrode", default="TT6")
    p.add_argument("--unit", type=int, required=True)
    p.add_argument("--thresholds", nargs="+", type=float, default=[4.0, 3.5])
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()
    exp = experiment_from_args(a)
    out = exp.sorting_dir(a.variant, subdir="sortings_opto")
    rec = si.load(str(out / "preprocessed")); main_sorting = si.load(str(out / "sorting"))
    main_an = si.load_sorting_analyzer(str(out / "analyzer"))
    t = rec.get_times(); t = t - t[0]
    ids = np.asarray(rec.channel_ids); groups = np.asarray(rec.get_property("group")).astype(int)
    g = int(a.tetrode.removeprefix("TT")) - 1; chans = np.flatnonzero(groups == g)
    gain = float(rec.get_property("gain_to_uV")[chans[0]])
    rec_tt = rec.select_channels(ids[chans])
    ref_frames = main_sorting.get_unit_spike_train(a.unit)
    ref_tmpl = main_an.get_extension("templates").get_unit_template(a.unit, operator="average")[:, chans]
    ref_vec = ref_tmpl.min(0); ref_vec = ref_vec / np.linalg.norm(ref_vec)

    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv"); pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv"); trains = trains[trains.session == exp.session]
    pmax = trains.power_mW.max(); sel = pulses[pulses.train.isin(trains[trains.power_mW == pmax].ttl_train)]
    on = sel.onset_s.values; dur = float(np.median(sel.offset_s - on))
    outside = np.ones(len(t), bool)
    for tr in trains.itertuples():
        outside[np.searchsorted(t, tr.onset_s - 1):np.searchsorted(t, tr.offset_s + 1)] = False

    rows = []
    for thr in a.thresholds:
        folder = out / f"resort_{a.tetrode}_thr{thr:g}"
        if folder.exists() and a.overwrite:
            shutil.rmtree(folder)
        if (folder / "sorting").exists():
            sorting = si.load(str(folder / "sorting"))
        else:
            sorting = ss.run_sorter("mountainsort5", rec_tt, folder=str(folder / "ms5"), verbose=False,
                                    scheme="2", filter=False, whiten=True, whitening_seed=42,
                                    detect_threshold=thr, detect_sign=-1, remove_existing_folder=True)
            sorting = sorting.save(folder=str(folder / "sorting"), overwrite=True)
        an = si.create_sorting_analyzer(sorting, rec_tt, format="memory", sparse=False)
        an.compute({"random_spikes": {"max_spikes_per_unit": 500, "seed": 42}, "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
                    "templates": {"operators": ["average"]}, "noise_levels": {},
                    "quality_metrics": {"metric_names": ["snr", "isi_violation", "rp_violation"]}}, n_jobs=8, progress_bar=False)
        qm = an.get_extension("quality_metrics").get_data(); tm = an.get_extension("templates")
        scale = 1.0 if an.return_in_uV else gain
        print(f"threshold {thr:g}: {len(sorting.unit_ids)} clusters on {a.tetrode}")
        fig, axes = plt.subplots(len(sorting.unit_ids), 3, figsize=(16, 2.6 * len(sorting.unit_ids)), squeeze=False)
        for k, u in enumerate(sorting.unit_ids):
            sp = sorting.get_unit_spike_train(u)
            tmpl = tm.get_unit_template(u, operator="average")
            vec = tmpl.min(0); vec = vec / np.linalg.norm(vec)
            dom = int(np.argmin(tmpl.min(0)))
            r = {"threshold": thr, "cluster": int(u), "n_spikes": int(len(sp)), "rate_hz": round(len(sp) / t[-1], 2),
                 "rate_outside_trains_hz": round(outside[sp].sum() / (outside.sum() / FS), 2),
                 "trough_uV": round(float(tmpl.min() * scale), 1), "dominant_wire": str(ids[chans[dom]]),
                 "wire_ratio_dom_to_2nd": round(float(np.sort(tmpl.min(0))[0] / np.sort(tmpl.min(0))[1]), 2),
                 "cos_to_ref": round(float(vec @ ref_vec), 3),
                 "frac_of_ref_spikes_matched": round(frac_matched(ref_frames, sp), 3),
                 "frac_of_own_spikes_in_ref": round(frac_matched(sp, ref_frames), 3),
                 "snr": round(float(qm.loc[u, "snr"]), 1), "isi_viol_ratio": round(float(qm.loc[u, "isi_violations_ratio"]), 3),
                 "rp_contamination": round(float(qm.loc[u, "rp_contamination"]), 3),
                 "base_hz": round(rate(sp, t, on, -0.05, 0), 1), "onset_1_10ms_hz": round(rate(sp, t, on, 0.001, 0.010), 1),
                 "in_pulse_hz": round(rate(sp, t, on, 0.010, dur), 1), "post_100ms_hz": round(rate(sp, t, on + dur, 0, 0.1), 1)}
            rows.append(r)
            ax = axes[k, 0]
            for c in range(len(chans)):
                ax.plot((np.arange(tmpl.shape[0]) - 30) / 30, tmpl[:, c] * scale, label=str(ids[chans[c]]))
            ax.set_title(f"thr {thr:g} cluster {u}: {len(sp)} spikes, {r['rate_hz']} Hz, cos {r['cos_to_ref']}, matches {r['frac_of_ref_spikes_matched']:.0%} of unit {a.unit}", fontsize=8)
            ax.legend(fontsize=6)
            ax = axes[k, 1]
            edges = np.arange(-0.05, dur + 0.1, 0.002)
            rel = np.concatenate([t[sp[np.searchsorted(sp, np.searchsorted(t, e - 0.05)):np.searchsorted(sp, np.searchsorted(t, e + dur + 0.1))]] - e for e in on])
            ax.hist(rel * 1000, bins=edges * 1000, weights=np.full(len(rel), 1 / (len(on) * 0.002)), color="k")
            ax.axvspan(0, dur * 1000, color="c", alpha=0.2); ax.set_title(f"PSTH around {len(on)} pulses at {pmax:g} mW (Hz)", fontsize=8)
            ax = axes[k, 2]
            amps = an.get_extension("waveforms").get_waveforms_one_unit(u)[:, :, dom].min(1) * scale
            sel_idx = an.get_extension("random_spikes").get_data()
            spk = an.sorting.to_spike_vector()[sel_idx]
            fr = spk["sample_index"][spk["unit_index"] == list(sorting.unit_ids).index(u)]
            ax.scatter(t[fr] / 60, amps, s=2, color="k")
            for tr in trains.itertuples():
                ax.axvspan(tr.onset_s / 60, tr.offset_s / 60, color="c", alpha=0.2)
            ax.set_title(f"trough amplitude on {ids[chans[dom]]} over the session (uV)", fontsize=8); ax.set_xlabel("min")
        fig.suptitle(f"{exp.experiment} {a.tetrode} re-sorted at threshold {thr:g} (reference: main-sort unit {a.unit})")
        fig.tight_layout(); fig.savefig(exp.figures_dir / f"opto_{a.tetrode}_resort_thr{thr:g}.png", dpi=100); plt.close(fig)
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / f"resort_{exp.experiment}_{a.tetrode}.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
