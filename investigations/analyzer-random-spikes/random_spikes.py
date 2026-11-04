"""PC-based quality metrics, SNR and templates of the final TTM-003-48H sorting against the
number of random spikes per unit the SortingAnalyzer draws (500, 2,000, 5,000; 20,000 as the
reference), three seeds each.

    .venv/bin/python investigations/analyzer-random-spikes/random_spikes.py <scratch dir>

Each (n, seed) run is kept in <scratch dir> (metrics CSV, average templates npz) and reused.
"""

import pathlib
import sys
import time

import numpy as np
import pandas as pd
import spikeinterface as si

from ttm.sorting import create_group_sparse_analyzer

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = pathlib.Path(__file__).parent / "results"
N_LIST = (500, 2000, 5000, 20000)
REF = 20000
SEEDS = (0, 1, 2)
METRICS = ["snr", "mahalanobis", "d_prime", "nearest_neighbor", "silhouette"]
COLS = ["snr", "isolation_distance", "l_ratio", "d_prime", "nn_hit_rate", "nn_miss_rate", "silhouette"]


def one_run(sorting, rec, n, seed, scratch):
    csv, npz = scratch / f"n{n}_s{seed}.csv", scratch / f"n{n}_s{seed}.npz"
    if csv.exists() and npz.exists():
        return pd.read_csv(csv), np.load(npz)["T"]
    t0 = time.time()
    az = create_group_sparse_analyzer(sorting, rec, format="memory", return_in_uV=True, n_jobs=48, seed=0)
    t1 = time.time()
    az.compute({
        "random_spikes": {"method": "uniform", "max_spikes_per_unit": n, "seed": seed},
        "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
        "templates": {"operators": ["average", "std"]},
        "noise_levels": {"random_slices_kwargs": {"seed": 0}},
        "principal_components": {"n_components": 5, "mode": "by_channel_local"},
    })
    t2 = time.time()
    az.compute("quality_metrics", metric_names=METRICS)
    t3 = time.time()
    qm = az.get_extension("quality_metrics").get_data()
    T = np.asarray(az.get_extension("templates").get_data(operator="average"), dtype=np.float32)
    drawn = np.bincount(sorting.to_spike_vector()["unit_index"][az.get_extension("random_spikes").get_data()],
                        minlength=len(sorting.unit_ids))
    d = pd.DataFrame({"unit": sorting.unit_ids, "n_drawn": drawn})
    for c in COLS:
        d[c] = qm.loc[sorting.unit_ids, c].to_numpy()
    d["n"], d["seed"] = n, seed
    d["s_main_channel"], d["s_extensions"], d["s_metrics"] = t1 - t0, t2 - t1, t3 - t2
    d.to_csv(csv, index=False)
    np.savez(npz, T=T)
    print(f"n={n} seed={seed}: main channels {t1 - t0:.0f} s, extensions {t2 - t1:.0f} s, "
          f"metrics {t3 - t2:.0f} s", flush=True)
    return d, T


def main():
    scratch = pathlib.Path(sys.argv[1])
    scratch.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(exist_ok=True)
    si.set_global_job_kwargs(n_jobs=48, chunk_duration="1s", progress_bar=False)
    rec = si.load(f"{RUN}/binary")
    sorting = si.load(f"{RUN}/sorting")
    counts = sorting.count_num_spikes_per_unit(outputs="array")
    tetrode = np.asarray(sorting.get_property("group")) + 1

    runs, temps = [], {}
    for n in N_LIST:
        for seed in SEEDS:
            d, T = one_run(sorting, rec, n, seed, scratch)
            runs.append(d)
            temps[(n, seed)] = T
    d = pd.concat(runs, ignore_index=True)
    d["tetrode"] = d["unit"].map(dict(zip(sorting.unit_ids, tetrode)))
    d["n_spikes"] = d["unit"].map(dict(zip(sorting.unit_ids, counts)))

    ref_T = np.mean([temps[(REF, s)] for s in SEEDS], axis=0)
    ref_flat = ref_T.reshape(len(sorting.unit_ids), -1)
    ref_peak = np.abs(ref_T).max(axis=(1, 2))
    cos, rel = {}, {}
    for (n, seed), T in temps.items():
        flat = T.reshape(len(sorting.unit_ids), -1)
        c = (flat * ref_flat).sum(1) / np.linalg.norm(flat, axis=1) / np.linalg.norm(ref_flat, axis=1)
        r = np.abs(T - ref_T).max(axis=(1, 2)) / ref_peak
        for u, cu, ru in zip(sorting.unit_ids, c, r):
            cos[(n, seed, u)], rel[(n, seed, u)] = cu, ru
    key = list(zip(d["n"], d["seed"], d["unit"]))
    d["template_cos_to_ref"] = [cos[k] for k in key]
    d["template_maxdiff_rel_to_ref"] = [rel[k] for k in key]
    d["peak_uv_ref"] = d["unit"].map(dict(zip(sorting.unit_ids, ref_peak)))
    d.to_csv(OUT / "random_spikes_units.csv", index=False)

    # Only units with more spikes than the largest cap are compared across caps: a unit drawn in
    # full at some cap keeps its own count while its neighbours' counts grow, which changes the
    # PC metrics by their definition, not by sampling. Those units are listed separately.
    full = d.groupby("unit")["n_spikes"].first()
    capped = full.index[full > REF]
    rows = []
    ref = d[d.n == REF].groupby("unit")[COLS].mean()
    for n in N_LIST:
        sub = d[d.n == n]
        means = sub.groupby("unit")[COLS].mean()
        for c in COLS:
            wide = sub.pivot(index="unit", columns="seed", values=c)
            units = wide.index[wide.index.isin(capped) & wide.notna().all(axis=1)]
            w = wide.loc[units]
            between_units = float(ref.loc[units, c].std())
            seed_sd = w.std(axis=1)
            dev = w.sub(ref.loc[units, c], axis=0).abs()
            shift = (means.loc[units, c] - ref.loc[units, c]) / ref.loc[units, c].abs()
            rank = [float(pd.concat([w[s], ref.loc[units, c]], axis=1).corr(method="spearman").iloc[0, 1])
                    for s in SEEDS]
            rows.append({
                "n": n, "metric": c, "units": len(units),
                "seed_sd_median": float(seed_sd.median()),
                "seed_sd_max": float(seed_sd.max()),
                "seed_sd_over_unit_sd": float(seed_sd.median() / between_units) if between_units else np.nan,
                "shift_rel_median": float(shift.median()),
                "shift_up_fraction": float((shift > 0).mean()),
                "absdev_to_ref_median": float(dev.stack().median()),
                "absdev_to_ref_p90": float(dev.stack().quantile(0.9)),
                "spearman_to_ref_min": min(rank),
            })
        units = sub["unit"].isin(capped)
        rows.append({"n": n, "metric": "template_cos_to_ref", "units": int(units.sum() / len(SEEDS)),
                     "absdev_to_ref_median": float(1 - sub.loc[units, "template_cos_to_ref"].median()),
                     "absdev_to_ref_p90": float(1 - sub.loc[units, "template_cos_to_ref"].quantile(0.1))})
        rows.append({"n": n, "metric": "template_maxdiff_rel_to_ref", "units": int(units.sum() / len(SEEDS)),
                     "absdev_to_ref_median": float(sub.loc[units, "template_maxdiff_rel_to_ref"].median()),
                     "absdev_to_ref_p90": float(sub.loc[units, "template_maxdiff_rel_to_ref"].quantile(0.9))})
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "random_spikes_summary.csv", index=False)
    small = (d[~d["unit"].isin(capped)].groupby(["unit", "tetrode", "n_spikes", "n"])[["n_drawn", *COLS]]
             .mean().reset_index())
    small.to_csv(OUT / "random_spikes_small_units.csv", index=False)
    print(small.round(3).to_string(index=False))
    t = d.groupby(["n", "seed"])[["s_main_channel", "s_extensions", "s_metrics"]].first().reset_index()
    t.to_csv(OUT / "random_spikes_timing.csv", index=False)
    pd.set_option("display.width", 250)
    print(s.to_string(index=False))
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
