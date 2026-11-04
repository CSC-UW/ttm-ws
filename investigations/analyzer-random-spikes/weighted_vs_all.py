"""L-ratio and isolation distance of the final TTM-003-48H sorting from capped draws, plain (as
SpikeInterface computes them) and with each neighbour's drawn spikes weighted by its total count over
its drawn count, against the values on every spike.

    .venv/bin/python investigations/analyzer-random-spikes/weighted_vs_all.py <scratch dir on /nvme>

One PCA model (5 components per channel, fitted on up to 5,000 spikes per unit, seed 0) projects
every spike once (``run_for_all_spikes``, about 11 GB, kept in <scratch dir>); every capped draw is
then a subset of those projections, so the comparison isolates the estimator from the PCA fit. The
analyzer is in ADC counts because ``run_for_all_spikes`` reads raw traces whatever ``return_in_uV``
says; ``check_uv`` measures that on TT9.
"""

import pathlib
import sys
import time

import numpy as np
import pandas as pd
import scipy.stats
import spikeinterface as si

from ttm.sorting import create_group_sparse_analyzer

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = pathlib.Path(__file__).parent / "results"
CAPS = (500, 2000, 5000, 20000)
SEEDS = (0, 1, 2, 3, 4)
EXT = {"waveforms": {"ms_before": 1.0, "ms_after": 2.0},
       "principal_components": {"n_components": 5, "mode": "by_channel_local"}}


def d2(X, mu, VI, step=1_000_000):
    out = np.empty(len(X))
    for i in range(0, len(X), step):
        d = X[i:i + step] - mu
        out[i:i + step] = np.einsum("ij,jk,ik->i", d, VI, d)
    return out


def isolation_l_ratio(own, other, w_other=None, n_own_total=None):
    """Isolation distance and L-ratio of ``own`` against ``other`` (rows are spikes). With
    ``w_other``, each other spike stands for ``w_other`` spikes and the unit for ``n_own_total``."""
    mu = own.mean(0)
    VI = np.linalg.inv(np.cov(own.T))
    D2 = d2(other, mu, VI)
    sf = scipy.stats.chi2.sf(D2, own.shape[1])
    if w_other is None:
        n_own, n_other = len(own), len(other)
        lr = sf.sum() / n_own
        n = min(n_own, n_other)
        iso = np.partition(D2, n - 1)[n - 1] if n >= 2 else np.nan
        return iso, lr
    lr = (w_other * sf).sum() / n_own_total
    n = min(n_own_total, w_other.sum())
    o = np.argsort(D2)
    k = np.searchsorted(np.cumsum(w_other[o]), n - 0.5)
    return D2[o][min(k, len(o) - 1)], lr


def check_uv(sorting, rec, scratch):
    """Max difference between the analyzer's own projections of its random spikes and
    ``run_for_all_spikes`` at the same spikes, in µV and in counts, on TT9."""
    tt9 = sorting.select_units(sorting.unit_ids[np.asarray(sorting.get_property("group")) == 8])
    rows = []
    for uv in (True, False):
        az = create_group_sparse_analyzer(tt9, rec, format="memory", return_in_uV=uv, n_jobs=48, seed=0)
        az.compute({"random_spikes": {"method": "uniform", "max_spikes_per_unit": 500, "seed": 0}, **EXT})
        pc = az.get_extension("principal_components")
        f = scratch / f"tt9_uv{int(uv)}.npy"
        pc.run_for_all_spikes(file_path=f)
        allp = np.load(f)
        idx = az.get_extension("random_spikes").get_data()
        own = pc.data["pca_projection"]
        rows.append({"return_in_uV": uv, "max_abs_diff": float(np.abs(allp[idx] - own).max()),
                     "own_abs_max": float(np.abs(own).max()), "all_abs_max": float(np.abs(allp[idx]).max())})
    return pd.DataFrame(rows)


def compute(scratch):
    si.set_global_job_kwargs(n_jobs=48, chunk_duration="1s", progress_bar=False)
    rec = si.load(f"{RUN}/binary")
    sorting = si.load(f"{RUN}/sorting")

    uv = check_uv(sorting, rec, scratch)
    uv.to_csv(OUT / "run_for_all_spikes_uv.csv", index=False)
    print(uv.to_string(index=False), flush=True)

    allp = scratch / "all_pcs.npy"
    az = create_group_sparse_analyzer(sorting, rec, format="memory", return_in_uV=False, n_jobs=48, seed=0)
    az.compute({"random_spikes": {"method": "uniform", "max_spikes_per_unit": 5000, "seed": 0}, **EXT})
    pc = az.get_extension("principal_components")
    if not allp.exists():
        t0 = time.time()
        pc.run_for_all_spikes(file_path=allp)
        print(f"all-spike projections in {(time.time() - t0) / 60:.1f} min", flush=True)
    P = np.load(allp)
    idx = az.get_extension("random_spikes").get_data()
    print(f"random spikes against all-spike projections: max |diff| "
          f"{np.abs(P[idx] - pc.data['pca_projection']).max():.2e}", flush=True)

    labels = sorting.to_spike_vector()["unit_index"]
    group = np.asarray(sorting.get_property("group"))
    nch = {g: int(np.sum(az.sparsity.mask[np.flatnonzero(group == g)[0]])) for g in np.unique(group)}
    rows = []
    t0 = time.time()
    for g in np.unique(group):
        units = np.flatnonzero(group == g)
        if len(units) < 2:
            continue
        sel = np.flatnonzero(np.isin(labels, units))
        X = P[sel][:, :, :nch[g]].reshape(len(sel), -1).astype(np.float64)
        lab = labels[sel]
        by_unit = {u: np.flatnonzero(lab == u) for u in units}
        N = {u: len(by_unit[u]) for u in units}
        for u in units:
            others = np.concatenate([by_unit[v] for v in units if v != u])
            iso, lr = isolation_l_ratio(X[by_unit[u]], X[others])
            row = {"unit": sorting.unit_ids[u], "tetrode": int(g) + 1, "n_spikes": N[u],
                   "n_other": len(others), "iso_all": iso, "lr_all": lr}
            for cap in CAPS:
                for seed in SEEDS:
                    rng = np.random.default_rng([seed, cap])
                    draw = {v: rng.choice(by_unit[v], size=min(cap, N[v]), replace=False) for v in units}
                    o = np.concatenate([draw[v] for v in units if v != u])
                    w = np.concatenate([np.full(len(draw[v]), N[v] / len(draw[v])) for v in units if v != u])
                    row[f"iso_cap_{cap}_s{seed}"], row[f"lr_cap_{cap}_s{seed}"] = isolation_l_ratio(X[draw[u]], X[o])
                    row[f"iso_w_{cap}_s{seed}"], row[f"lr_w_{cap}_s{seed}"] = isolation_l_ratio(
                        X[draw[u]], X[o], w_other=w, n_own_total=N[u])
            rows.append(row)
        print(f"TT{int(g) + 1}: {len(units)} units, {len(sel):,} spikes, {time.time() - t0:.0f} s", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "weighted_vs_all_units.csv", index=False)
    return d


def main():
    scratch = pathlib.Path(sys.argv[1])
    scratch.mkdir(parents=True, exist_ok=True)
    units_csv = OUT / "weighted_vs_all_units.csv"
    d = pd.read_csv(units_csv) if units_csv.exists() else compute(scratch)

    # Isolation distance is undefined when a unit has at least as many spikes as all the other units on
    # its tetrode; SpikeInterface then returns the farthest other spike, a maximum that no draw estimates.
    defined = d["n_spikes"] < d["n_other"]
    print("isolation distance left out (own spikes >= other units' spikes):\n"
          + d.loc[~defined, ["unit", "tetrode", "n_spikes", "n_other", "iso_all"]].to_string(index=False))
    summ = []
    for m in ("iso", "lr"):
        keep = defined if m == "iso" else np.ones(len(d), bool)
        exact = d.loc[keep, f"{m}_all"]
        for cap in CAPS:
            for kind in ("cap", "w"):
                est = d.loc[keep, [f"{m}_{kind}_{cap}_s{s}" for s in SEEDS]]
                diff = est.sub(exact, axis=0)
                rel = diff.div(exact.abs(), axis=0)
                med = rel.median(axis=1)
                rank = min(float(pd.concat([est.iloc[:, i], exact], axis=1).corr(method="spearman").iloc[0, 1])
                           for i in range(len(SEEDS)))
                summ.append({"metric": {"iso": "isolation_distance", "lr": "l_ratio"}[m], "cap": cap,
                             "estimate": {"cap": "capped", "w": "weighted"}[kind], "units": int(exact.notna().sum()),
                             "median_rel_error": float(med.median()),
                             "median_abs_rel_error": float(med.abs().median()),
                             "max_abs_rel_error": float(med.abs().max()),
                             "median_abs_error": float(diff.median(axis=1).abs().median()),
                             "seed_rel_sd_median": float(rel.std(axis=1).median()),
                             "spearman_to_all_min": rank})
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "weighted_vs_all_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(s.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
