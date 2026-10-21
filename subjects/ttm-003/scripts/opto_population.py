"""Population effects of the laser trains on the sorted units of one session.

Pools the units that the tagging step classed as suppressed and, per train,
measures the time course of their summed rate around the pulse: latency to
silence, duration of the silent period after light offset (time until the
pooled rate recovers to half its baseline), and the rebound after it. Also
the pooled onset volley (1 to 10 ms) across all units. Writes
``results/opto_population_<experiment>.json`` and a per-train figure.

Usage: python subjects/ttm-003/scripts/opto_population.py --subject TTM-003 --experiment TTM-003-OPTO1
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
BIN = 0.005


def pooled_psth(trains_spikes, events, lo, hi, bin_s=BIN):
    edges = np.arange(lo, hi + bin_s / 2, bin_s)
    counts = np.zeros(len(edges) - 1)
    for sp in trains_spikes:
        i0 = np.searchsorted(sp, events + lo); i1 = np.searchsorted(sp, events + hi)
        rel = np.concatenate([sp[a:b] - e for a, b, e in zip(i0, i1, events)]) if len(events) else np.array([])
        counts += np.histogram(rel, edges)[0]
    return edges[:-1] + bin_s / 2, counts / (len(events) * bin_s)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis")
    a = p.parse_args()
    exp = experiment_from_args(a)
    out_dir = exp.sorting_dir(a.variant, subdir="sortings_opto")
    sorting = si.load(str(out_dir / "sorting"))
    rec = si.load(str(out_dir / "preprocessed"))
    t = rec.get_times(); t0 = t[0]
    spikes = {u: t[sorting.get_unit_spike_train(u)] - t0 for u in sorting.unit_ids}
    tag = pd.read_csv(RESULTS / f"optotag_{exp.experiment}_{a.variant}.csv")
    tag = tag[tag.set == "B"].set_index("unit")
    suppressed = [u for u in sorting.unit_ids if tag.loc[u, "class"] == "suppressed"]
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv"); pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv"); trains = trains[trains.session == exp.session]

    rows = []
    fig, axes = plt.subplots(len(trains), 1, figsize=(14, 2.6 * len(trains)), sharex=False)
    for ax, tr in zip(np.atleast_1d(axes), trains.itertuples()):
        pp = pulses[pulses.train == tr.ttl_train]
        on, off = pp.onset_s.values, pp.offset_s.values
        dur = float(np.median(off - on))
        tb, r = pooled_psth([spikes[u] for u in suppressed], on, -0.1, dur + 0.3)
        base = r[tb < 0].mean()
        half = base / 2
        # latency to silence: first bin after onset with rate < half baseline
        after = tb >= 0
        below = (r < half) & after
        t_sil = float(tb[below][0]) if below.any() else np.nan
        # recovery: first bin after offset where the rate is back above half baseline and stays for 3 bins
        post = np.where(tb > dur)[0]
        t_rec = np.nan
        for i in post:
            if i + 3 <= len(r) and np.all(r[i:i + 3] >= half):
                t_rec = float(tb[i]); break
        off_period = t_rec if np.isnan(t_rec) else t_rec  # from onset
        reb = r[(tb > dur + 0.05) & (tb < dur + 0.25)].max() / base if base > 0 else np.nan
        min_rate = r[(tb > 0.02) & (tb < dur)].mean() / base if base > 0 else np.nan
        _, r_all = pooled_psth([spikes[u] for u in sorting.unit_ids], on, -0.05, 0.05, 0.001)
        tb1 = np.arange(-0.05, 0.05, 0.001) + 0.0005
        b_all = r_all[tb1 < 0].mean()
        volley = r_all[(tb1 >= 0.001) & (tb1 < 0.010)].max() / b_all if b_all > 0 else np.nan
        rows.append({"train": int(tr.n), "pattern": tr.ttl_pattern, "power_mW": float(tr.power_mW), "state": tr.state,
                     "pulse_s": round(dur, 3), "n_suppressed_units": len(suppressed), "pooled_base_hz": float(base),
                     "silence_onset_ms": None if np.isnan(t_sil) else round(t_sil * 1000, 1),
                     "in_pulse_rate_frac_20ms_to_end": None if np.isnan(min_rate) else round(float(min_rate), 3),
                     "recovery_from_onset_ms": None if np.isnan(t_rec) else round(t_rec * 1000, 1),
                     "silence_after_offset_ms": None if np.isnan(t_rec) else round((t_rec - dur) * 1000, 1),
                     "rebound_peak_frac": None if np.isnan(reb) else round(float(reb), 2),
                     "onset_volley_peak_frac_all_units": None if np.isnan(volley) else round(float(volley), 2)})
        ax.plot(tb * 1000, r, color="k", lw=1)
        ax.axhline(base, color="gray", lw=0.5); ax.axhline(half, color="gray", lw=0.5, ls="--")
        ax.axvspan(0, dur * 1000, color="c", alpha=0.2)
        if not np.isnan(t_rec): ax.axvline(t_rec * 1000, color="r", lw=0.8)
        ax.set_title(f"train {tr.n} ({tr.ttl_pattern}, {tr.power_mW:g} mW, {tr.state}): pooled rate of {len(suppressed)} suppressed units, 5 ms bins", fontsize=9)
        ax.set_ylabel("Hz")
    np.atleast_1d(axes)[-1].set_xlabel("ms from onset")
    fig.tight_layout(); fig.savefig(exp.figures_dir / f"opto_population_{a.variant}.png", dpi=100); plt.close(fig)
    df = pd.DataFrame(rows)
    (RESULTS / f"opto_population_{exp.experiment}_{a.variant}.json").write_text(json.dumps(rows, indent=1) + "\n")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
