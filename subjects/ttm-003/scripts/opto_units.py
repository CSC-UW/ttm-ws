"""Sorted units' responses to the laser pulses over a stretch of a 48 h recording, and the pulse rule's test.

Sorts a stretch of the run's preprocessed binary (``ttm-sort binary``) per tetrode with MountainSort5 as the
seed sort does (``ttm.tracking.sort_chunk``: scheme 2, threshold 5.5, whitening seed 42), training on
``--train-s`` of it, and keeps units with template snr 5 or more and 100 spikes or more (the seed's
criteria). For every kept unit, around every pulse of the stretch and around sham onsets (the same period,
ending 2 s before each train's first pulse):

- rates in the 500 ms before onset (baseline), 1 to 10 ms after onset (the onset volley), 10 ms to pulse
  end, and 0 to 50, 50 to 100, 100 to 200 and 200 to 400 ms after offset, and a 5 ms PSTH;
- the silence after the light: the first 5 ms bin after offset where the unit's pulse PSTH, smoothed over
  three bins, regains half its baseline (ms from onset);
- the pulse rule's driven test (``docs/sorting_pipeline.md``, step 5.4) on its spikes, per 1800 s window of
  the recording: from 15 ms after each onset to the offset, 20 spikes or more, at least twice the count
  expected from its rate outside the pulses (onset to offset plus 100 ms left out), one-sided Poisson
  p under 1e-4;
- the median trough of its spikes on its main wire inside pulses (the first 500) against outside (500 spread
  over the stretch), also by position in the pulse as a share of the median outside, with the share under
  half; its spikes per pulse where it fires in the pulse, and the median interval between them;
- per train, its in-pulse rate over its baseline.

Writes ``results/opto_units_<E>_<name>.csv`` (one row per unit), ``results/opto_units_<E>_<name>_trains.csv``
(one row per unit and train), ``results/opto_units_<E>_<name>_windows.csv`` (the driven test per unit and window)
and ``results/opto_units_<E>_<name>.json``, and a PSTH figure per tetrode in
``<experiment>/figures/opto_units_<name>/``. The sort is kept in ``<experiment>/opto_units/<name>/``.

Usage: python subjects/ttm-003/scripts/opto_units.py --subject TTM-003 --experiment TTM-003-48H-2
           --start-h 23.97 --end-h 28.97 --name sd [--train-s 1800]
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
import numpy as np
import pandas as pd
import scipy.stats as st
import spikeinterface as si

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ttm import carry as cr
from ttm.cli._common import add_experiment_args, experiment_from_args
from ttm.sorting import bad_channel_ids

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
PRE_S, POST_S = 0.5, 0.95
WINDOW_S = 1800.0


def rate(rel: np.ndarray, n: int, a: float, b: float) -> float:
    return float(((rel >= a) & (rel < b)).sum() / (n * (b - a))) if n else float("nan")


def driven_test(t: np.ndarray, on: np.ndarray, off: np.ndarray, a: float, b: float) -> dict:
    """The pulse rule's driven test on spike times ``t`` (s) for the pulses in [a, b)."""
    m = (on >= a) & (on < b)
    on, off = on[m], off[m]
    if on.size == 0:
        return {"pulses": 0}
    tt = t[(t >= a) & (t < b)]
    inside = sum(int(((tt >= o + 0.015) & (tt < f)).sum()) for o, f in zip(on, off))
    excl = sum(int(((tt >= o) & (tt < f + 0.1)).sum()) for o, f in zip(on, off))
    dur_in = float(np.sum(off - (on + 0.015)))
    dur_out = (b - a) - float(np.sum(off + 0.1 - on))
    out_rate = (tt.size - excl) / dur_out
    expected = out_rate * dur_in
    p = float(st.poisson.sf(inside - 1, expected)) if expected > 0 else (0.0 if inside else 1.0)
    return {"start_s": float(a), "pulses": int(on.size), "count": inside, "expected": float(expected), "p": p,
            "driven": bool(inside >= 20 and inside >= 2 * expected and p < 1e-4)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--start-h", type=float, required=True)
    ap.add_argument("--end-h", type=float, required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--train-s", type=float, default=1800.0)
    ap.add_argument("--run", default="main")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    run_dir = exp.sorting_dir(a.run, subdir="sortings_carry")
    store = si.read_zarr(str(exp.store_path()))
    drop = bad_channel_ids(store, exp.badChannels, exp.badTetrodes)
    rec = cr.open_sorting_recording(run_dir / "binary", drop)
    fs = rec.get_sampling_frequency()
    fa, fb = int(a.start_h * 3600 * fs), min(int(a.end_h * 3600 * fs), rec.get_num_frames())
    win = rec.frame_slice(fa, fb)
    win.reset_times()
    out_dir = exp.experiment_dir / "opto_units" / a.name
    out_dir.mkdir(parents=True, exist_ok=True)
    if (out_dir / "by_group").exists():
        sorting = cr.sort_window(win, out_dir, cr.CarryParams())
    else:
        from ttm.tracking import sort_chunk, to_int_numpy_sorting
        sorting = to_int_numpy_sorting(sort_chunk(win, out_dir, whitening_seed=42, detect_threshold=5.5,
                                                  scheme2_training_duration_sec=a.train_s, sort_n_jobs=16))
    dense, snr = cr.unit_templates(sorting, win, n_jobs=32, with_snr=True, seed=0)
    nsp = np.array([len(sorting.get_unit_spike_train(u)) for u in sorting.unit_ids])
    keep = (snr >= 5) & (nsp >= 100)
    groups = np.asarray(sorting.get_property("group")).astype(int)
    print(f"{sorting.get_num_units()} units, {int(keep.sum())} with snr >= 5 and >= 100 spikes")

    blocks = pd.read_csv(RESULTS / f"stimulation_blocks_{exp.session}.csv")
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    pulses = pulses[~pulses.glitch & pulses.train.isin(blocks.loc[blocks.delivered, "ttl_train"])]
    t_a, t_b = fa / fs, fb / fs
    pp = pulses[(pulses.onset_s >= t_a + PRE_S) & (pulses.onset_s < t_b - POST_S)].sort_values("onset_s")
    on, off, train = pp.onset_s.to_numpy(), pp.offset_s.to_numpy(), pp.train.to_numpy()
    sham = []
    for tr in np.unique(train):
        o = on[train == tr]
        per = float(np.median(np.diff(o))) if o.size > 1 else 1.0
        s = o[0] - 2.0 - per * np.arange(min(40, o.size))[::-1]
        sham.append(s[s >= t_a + PRE_S])
    sham = np.concatenate(sham) if sham else np.zeros(0)
    pulse_s = float(np.median(off - on))
    bins = np.arange(-PRE_S, POST_S + 1e-9, 0.005)
    centers = bins[:-1] + 0.0025

    rows, train_rows, window_rows, psths = [], [], [], {}
    for i, u in enumerate(sorting.unit_ids):
        if not keep[i]:
            continue
        g = int(groups[i])
        frames = sorting.get_unit_spike_train(u)
        t = frames / fs + t_a
        rel_list = [t[(t >= e - PRE_S) & (t < e + POST_S)] - e for e in on]
        rel = np.concatenate(rel_list) if rel_list else np.zeros(0)
        rel_s = np.concatenate([t[(t >= e - PRE_S) & (t < e + POST_S)] - e for e in sham]) if sham.size else np.zeros(0)
        h = np.histogram(rel, bins)[0] / (on.size * 0.005)
        hs = np.histogram(rel_s, bins)[0] / (max(sham.size, 1) * 0.005)
        psths[int(u)] = (g, h, hs)
        base = rate(rel, on.size, -PRE_S, -0.005)
        sm = np.convolve(h, np.ones(3) / 3, mode="same")
        after = np.flatnonzero((centers >= pulse_s) & (sm >= 0.5 * base))
        silence_end = float(centers[after[0]] * 1000) if after.size and base > 0 else float("nan")
        tests = [driven_test(t, on, off, max(w0, t_a), min(w0 + WINDOW_S, t_b))
                 for w0 in np.arange(np.floor(t_a / WINDOW_S) * WINDOW_S, t_b, WINDOW_S)]
        tests = [d for d in tests if d.get("pulses", 0)]
        window_rows += [{"unit": int(u), "tetrode": f"TT{g + 1}", **d} for d in tests]
        chans = cr.cp.tetrode_channel_indices(win, g)
        main = int(np.argmin(dense[i][:, chans].min(0)))
        tr_in, tr_out = [], []
        in_pulse = np.zeros(t.size, bool)
        j = np.searchsorted(on, t, "right") - 1
        ok = j >= 0
        in_pulse[ok] = t[ok] < off[j[ok]]
        sel_in = np.flatnonzero(in_pulse)[:500]
        sel_out = np.flatnonzero(~in_pulse)
        if sel_out.size:
            sel_out = sel_out[np.linspace(0, sel_out.size - 1, min(500, sel_out.size)).astype(int)]
        ch_id = win.channel_ids[chans[main]]
        lat_in = []
        for sel, acc in ((sel_in, tr_in), (sel_out, tr_out)):
            for k_s in sel:
                f = frames[k_s]
                if f - 5 >= 0 and f + 6 <= win.get_num_frames():
                    x = win.get_traces(start_frame=int(f) - 5, end_frame=int(f) + 6, channel_ids=[ch_id], return_in_uV=True)
                    acc.append(float(x.min()))
                    if acc is tr_in:
                        lat_in.append(float(t[k_s] - on[j[k_s]]))
        size_in = np.array(tr_in) / np.median(tr_out) if tr_in and tr_out else np.zeros(0)
        lat_in = np.array(lat_in)

        def size_between(lo, hi, size_in=size_in, lat_in=lat_in):
            m = (lat_in >= lo) & (lat_in < hi)
            return float(np.median(size_in[m])) if m.any() else float("nan")

        per_pulse = np.bincount(j[in_pulse], minlength=on.size)
        per_pulse = per_pulse[per_pulse > 0]
        same_pulse = np.diff(j[in_pulse]) == 0
        isi_in = np.diff(t[in_pulse])[same_pulse]
        rows.append({
            "unit": int(u), "tetrode": f"TT{g + 1}", "spikes": int(t.size), "snr": round(float(snr[i]), 2),
            "baseline_hz": base, "onset_1_10ms_hz": rate(rel, on.size, 0.001, 0.010),
            "pulse_10ms_end_hz": rate(rel, on.size, 0.010, pulse_s),
            "post_0_50ms_hz": rate(rel, on.size, pulse_s, pulse_s + 0.05),
            "post_50_100ms_hz": rate(rel, on.size, pulse_s + 0.05, pulse_s + 0.1),
            "post_100_200ms_hz": rate(rel, on.size, pulse_s + 0.1, pulse_s + 0.2),
            "post_200_400ms_hz": rate(rel, on.size, pulse_s + 0.2, pulse_s + 0.4),
            "sham_10ms_end_hz": rate(rel_s, sham.size, 0.010, pulse_s),
            "silence_end_ms": silence_end,
            "driven_windows": int(sum(d["driven"] for d in tests)), "tested_windows": len(tests),
            "driven_count_max": int(max((d["count"] for d in tests), default=0)),
            "trough_in_pulse_uV_median": float(np.median(tr_in)) if tr_in else float("nan"),
            "trough_outside_uV_median": float(np.median(tr_out)) if tr_out else float("nan"),
            "spikes_in_pulse": len(tr_in),
            "size_in_pulse_0_15ms": size_between(0, 0.015), "size_in_pulse_15_30ms": size_between(0.015, 0.03),
            "size_in_pulse_30ms_off": size_between(0.03, 1.0),
            "size_in_pulse_under_half": float(np.mean(size_in < 0.5)) if size_in.size else float("nan"),
            "spikes_in_pulse_all": int(in_pulse.sum()),
            "spikes_per_responding_pulse_mean": float(per_pulse.mean()) if per_pulse.size else 0.0,
            "spikes_per_responding_pulse_max": int(per_pulse.max()) if per_pulse.size else 0,
            "isi_in_pulse_ms_median": float(np.median(isi_in) * 1000) if isi_in.size else float("nan"),
        })
        for tr in np.unique(train):
            o = on[train == tr]
            rl = np.concatenate([t[(t >= e - PRE_S) & (t < e + POST_S)] - e for e in o])
            b0 = rate(rl, o.size, -PRE_S, -0.005)
            train_rows.append({"unit": int(u), "tetrode": f"TT{g + 1}", "ttl_train": int(tr),
                               "onset_h": float(o[0] / 3600), "baseline_hz": b0,
                               "pulse_10ms_end_hz": rate(rl, o.size, 0.010, pulse_s),
                               "post_0_100ms_hz": rate(rl, o.size, pulse_s, pulse_s + 0.1)})
    df = pd.DataFrame(rows)
    df["pulse_over_baseline"] = df.pulse_10ms_end_hz / df.baseline_hz
    df["onset_over_baseline"] = df.onset_1_10ms_hz / df.baseline_hz
    df["class"] = np.where(df.driven_windows > 0, "driven",
                           np.where(df.pulse_over_baseline < 0.5, "suppressed",
                                    np.where(df.pulse_over_baseline > 2, "activated", "unaffected")))
    tag = f"{exp.experiment}_{a.name}"
    df.round(4).to_csv(RESULTS / f"opto_units_{tag}.csv", index=False)
    tdf = pd.DataFrame(train_rows)
    tdf["pulse_over_baseline"] = tdf.pulse_10ms_end_hz / tdf.baseline_hz
    tdf.round(4).to_csv(RESULTS / f"opto_units_{tag}_trains.csv", index=False)
    pd.DataFrame(window_rows).round(6).to_csv(RESULTS / f"opto_units_{tag}_windows.csv", index=False)
    trend = []
    for u, q in tdf.groupby("unit"):
        q = q[np.isfinite(q.pulse_over_baseline)]
        if len(q) >= 5:
            rho, p = st.spearmanr(q.onset_h, q.pulse_over_baseline)
            trend.append({"unit": int(u), "tetrode": q.tetrode.iloc[0], "spearman_rho": float(rho), "p": float(p),
                          "first_third": float(q.pulse_over_baseline.iloc[:len(q) // 3].mean()),
                          "last_third": float(q.pulse_over_baseline.iloc[-(len(q) // 3):].mean())})
    summary = {"experiment": exp.experiment, "name": a.name, "start_h": a.start_h, "end_h": a.end_h,
               "pulses": int(on.size), "sham_onsets": int(sham.size), "pulse_ms": pulse_s * 1000,
               "units_sorted": int(sorting.get_num_units()), "units_kept": int(keep.sum()),
               "classes": df["class"].value_counts().to_dict(),
               "driven_units": df[df["class"] == "driven"][["unit", "tetrode", "driven_windows", "tested_windows",
                                                              "onset_over_baseline", "pulse_over_baseline"]].to_dict(orient="records"),
               "suppressed_silence_end_ms_median": float(df.loc[df["class"] == "suppressed", "silence_end_ms"].median()),
               "trend_per_unit": trend}
    (RESULTS / f"opto_units_{tag}.json").write_text(json.dumps(summary, indent=1) + "\n")

    figdir = exp.figures_dir / f"opto_units_{a.name}"
    figdir.mkdir(parents=True, exist_ok=True)
    for g in sorted({v[0] for v in psths.values()}):
        us = [u for u, v in psths.items() if v[0] == g]
        fig, axes = plt.subplots(len(us), 1, figsize=(10, 1.8 * len(us) + 0.6), sharex=True, squeeze=False)
        for ax, u in zip(axes[:, 0], us):
            _, h, hs = psths[u]
            r = df[df.unit == u].iloc[0]
            ax.plot(centers * 1000, h, "k", lw=0.8)
            ax.plot(centers * 1000, hs, color="0.6", lw=0.8)
            ax.axvspan(0, pulse_s * 1000, color="c", alpha=0.2)
            ax.set_title(f"unit {u} ({r['class']}): baseline {r.baseline_hz:.1f} Hz, in-pulse/baseline "
                         f"{r.pulse_over_baseline:.2f}, silence ends {r.silence_end_ms:.0f} ms", fontsize=8)
        axes[-1, 0].set_xlabel("ms from onset (black: pulses, grey: sham onsets)")
        fig.suptitle(f"{exp.experiment} {a.name}: TT{g + 1}, {on.size} pulses", fontsize=9)
        fig.tight_layout(); fig.savefig(figdir / f"TT{g + 1}.png", dpi=90); plt.close(fig)
    print(df[["unit", "tetrode", "spikes", "snr", "baseline_hz", "onset_over_baseline", "pulse_over_baseline",
              "silence_end_ms", "driven_windows", "class"]].round(2).to_string(index=False))
    print(json.dumps({k: v for k, v in summary.items() if k != "trend_per_unit"}, indent=1))


if __name__ == "__main__":
    main()
