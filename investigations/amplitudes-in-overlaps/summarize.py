"""Summary tables of amp_le0.py, amp_le0_sub.py, overlap_fits.py, overlap_timing.py and amp_cutoff_check.py.

Arguments: the folder with amp_le0_units.csv, amp_le0_snippets.csv and amp_le0_sub.csv; the folder with
overlap_events_*_wNN.csv, overlap_lag_hist.csv and overlap_by_v_size.csv; the results folder, which must
already hold amp_cutoff_check.csv."""
import glob
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, poisson

A, B, R = sys.argv[1:4]


def tset(tt):
    tt = np.asarray(tt)
    return np.where(tt == 15, "TT15", np.where(np.isin(tt, [1, 2, 16]), "TT1/2/16",
                    np.where((tt >= 4) & (tt <= 12), "TT4-12", "TT3/13")))


# amplitudes of 0 or less per unit group
U = pd.read_csv(f"{A}/amp_le0_units.csv")
U["set"] = tset(U.tetrode)
g = U.groupby(["exp", "set"]).agg(units=("unit", "size"), units_with_le0=("n_le0", lambda v: int((v > 0).sum())),
                                  spikes=("spikes", "sum"), n_le0=("n_le0", "sum"),
                                  unit_frac_median=("frac_le0", "median"), unit_frac_max=("frac_le0", "max"))
g["pooled_frac"] = g.n_le0 / g.spikes
g.reset_index().to_csv(f"{R}/le0_by_group.csv", index=False)
U[U.peak_shift.abs() >= 10][["exp", "unit", "tetrode", "main_ch", "template_extremum_ch", "peak_shift", "spikes",
                              "frac_le0"]].to_csv(f"{R}/le0_far_peak_shift_units.csv", index=False)

# context of the sampled spikes
N = pd.read_csv(f"{A}/amp_le0_snippets.csv").merge(U[["exp", "unit", "tetrode", "median_amp", "peak_shift"]],
                                                    on=["exp", "unit"])
N["set"] = tset(N.tetrode)
N.loc[N.peak_shift.abs() >= 10, "set"] = "far peak shift"
N["shift_ge_3"] = N.lag_best.abs() >= 3
N["shift_at_limit"] = N.lag_best.abs() == 10
N["cos_raw_lt_0.307"] = N.cos_raw < 0.307
flags = ["in_bad", "in_pulse", "nb_sorted", "nb_mua", "nb_self", "shift_ge_3", "shift_at_limit", "cos_raw_lt_0.307",
         "tet_argmax_center_is_main"]
meds = ["snr_pct", "cos_pl_pct", "cos_raw_pct", "a_pct", "trough_pct", "cos0", "cos_best"]
c = N.groupby(["set", "kind"])[flags].mean().join(N.groupby(["set", "kind"])[meds].median())
c["n"] = N.groupby(["set", "kind"]).size()
c["lag_best_abs_median"] = N.assign(x=N.lag_best.abs()).groupby(["set", "kind"]).x.median()
c.reset_index().to_csv(f"{R}/le0_context.csv", index=False)

# after subtracting the neighbouring sorted spikes
D = pd.read_csv(f"{A}/amp_le0_sub.csv").merge(U[["exp", "unit", "tetrode", "median_amp", "peak_shift"]],
                                               on=["exp", "unit"])
D["set"] = tset(D.tetrode)
D.loc[D.peak_shift.abs() >= 10, "set"] = "far peak shift"
D["rel"] = D.corrected / D.median_amp
D["has_other"] = D.n_nbr_other_unit > 0


def sub(x):
    return pd.Series({"n": len(x), "other_unit_within_snippet": x.has_other.mean(), "corrected_gt_0": (x.corrected > 0).mean(),
                      "corrected_gt_half_median": (x.rel > 0.5).mean(), "corrected_over_median": x.rel.median(),
                      "own_scale": x.self_scale.median()})


s1 = D.groupby(["set", "kind"]).apply(sub).reset_index()
D2 = D[(D.set == "TT4-12") & (D.kind == "le0")]
s2 = D2.groupby(D2.has_other.rename("other_unit_present")).apply(sub).reset_index()
s2.insert(0, "set", "TT4-12 le0 by other unit present")
pd.concat([s1, s2]).to_csv(f"{R}/le0_subtraction.csv", index=False)

# trace test
E = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{B}/overlap_events_*_w*.csv"))])
E["set"] = tset(E.group + 1)
thr = float(E.none_dU.quantile(0.99))
E["size"] = pd.cut(E.aU_rel, [0, 0.7, 1.0, np.inf], labels=["small (<0.7)", "mid (0.7-1)", "large (>=1)"])
E["lag"] = pd.cut(E.d, [-46, -22, -10, 9, 21, 46], labels=["U 0.7-1.5 ms before V", "U 0.33-0.7 ms before V",
                                                           "within 0.3 ms", "U 0.33-0.7 ms after V", "U 0.7-1.5 ms after V"])
E["v_size"] = pd.cut(E.aV_rel, [0, 0.8, 0.95, 1.05, 1.2, np.inf])


def tr(x):
    return pd.Series({"n": len(x), "obs_dU": x.obs_dU.median(), "twin_dU": x.twin_dU.median(), "none_dU": x.none_dU.median(),
                      "obs_over_twin": (x.obs_dU / x.twin_dU).median(), "obs_over_twin_p10": (x.obs_dU / x.twin_dU).quantile(0.1),
                      "obs_below_none_p99": (x.obs_dU < thr).mean(),
                      "twin_below_none_p99": (x.twin_dU < thr).mean(), "obs_u_ls": x.obs_u_ls.median(),
                      "twin_u_ls": x.twin_u_ls.median(), "aU_rel": x.aU_rel.median()})


parts = [E.groupby(["exp", "set"]).apply(tr).reset_index().assign(split="set")]
W = E[E.set == "TT4-12"]
for col in ("size", "lag", "v_size"):
    parts.append(W.groupby(["exp", col], observed=True).apply(tr).reset_index().rename(columns={col: "stratum"}).assign(split=col))
sm = W[W.aU_rel < 0.7]
parts.append(sm.groupby(["exp", "lag"], observed=True).apply(tr).reset_index().rename(columns={"lag": "stratum"}).assign(split="small by lag"))
q05 = W.cos_rawU.quantile(0.05)
parts.append(W[W.cos_rawU < q05].groupby("exp").apply(tr).reset_index().assign(split=f"raw cosine < {q05:.3f} (lowest 5 %)"))
pd.concat(parts).to_csv(f"{R}/overlap_trace_summary.csv", index=False)
tail = pd.DataFrame([{"energy": x, "none_ge": (W.none_dU >= x).mean(), "obs_lt": (W.obs_dU < x).mean(),
                      "twin_lt": (W.twin_dU < x).mean()} for x in (20, thr, 100, 150)])
tail.to_csv(f"{R}/overlap_trace_tails.csv", index=False)
pp = W.groupby(["exp", "group", "U", "V"]).agg(n=("obs_dU", "size"), obs_over_twin=("obs_dU", lambda v: np.nan),
                                                below=("obs_dU", lambda v: float((v < thr).mean()))).reset_index()
pp["obs_over_twin"] = W.assign(r=W.obs_dU / W.twin_dU).groupby(["exp", "group", "U", "V"]).r.median().to_numpy()
pp.to_csv(f"{R}/overlap_trace_per_pair.csv", index=False)

# timing
H = pd.read_csv(f"{B}/overlap_lag_hist.csv")
H = H[(H.group + 1).between(4, 12)].reset_index(drop=True)
L = np.arange(-150, 151)
cols = [f"l{v}" for v in L]
X = H[cols].to_numpy(float)
far, near = np.abs(L) >= 60, np.abs(L) <= 45
base = X[:, far].mean(axis=1)
ok = base * 91 >= 20
bins = [(-45, -22), (-21, -10), (-9, -4), (-3, 3), (4, 9), (10, 21), (22, 45)]
rows = []
for (e, cl), idx in H[ok].groupby(["exp", "u_class"]).groups.items():
    r = {"exp": e, "u_class": cl, "pairs": len(idx)}
    for a0, a1 in bins:
        m = (L >= a0) & (L <= a1)
        r[f"{a0}..{a1}"] = float(np.median(X[idx][:, m].mean(axis=1) / base[idx]))
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{R}/overlap_lag_profile.csv", index=False)


def peak(x):
    xs = np.convolve(x, np.ones(21) / 21, mode="same")
    p = poisson.sf(x[near] - 1, np.maximum(xs[near], 1e-9))
    return min(1.0, p.min() * 91)


rng = np.random.default_rng(0)
rows = []
for (e, g_, u, v), z in H.groupby(["exp", "group", "U", "V"]):
    z = z.set_index("u_class")
    if not {"small", "large"} <= set(z.index):
        continue
    xs_, xl = z.loc["small", cols].to_numpy(float), z.loc["large", cols].to_numpy(float)
    if xs_[near].sum() < 500 or xs_.sum() >= xl.sum():
        continue
    xl_t = rng.binomial(xl.astype(np.int64), xs_.sum() / xl.sum()).astype(float)
    b = np.add.reduceat(np.vstack([xs_[near], xl_t[near]]), np.arange(0, 91, 5), axis=1)
    b = b[:, b.sum(axis=0) > 0]
    rows.append({"exp": e, "group": g_, "U": u, "V": v, "small_near": xs_[near].sum(), "peak_p_small": peak(xs_),
                 "peak_p_large_thinned": peak(xl_t), "shape_p": chi2_contingency(b)[1]})
pk = pd.DataFrame(rows)
pk.to_csv(f"{R}/overlap_lag_peaks.csv", index=False)
C = pd.read_csv(f"{B}/overlap_by_v_size.csv")
C = C[(C.group + 1).between(4, 12)]
t = C.groupby(["exp", "v_bin"])[["near", "far"]].sum()
smc = C[C.u_class == "small"].groupby(["exp", "v_bin"])[["near", "far", "n_v"]].sum()
t["near_over_far_all"] = (t.near / 91) / (t.far / 182)
t["near_over_far_small"] = (smc.near / 91) / (smc.far / 182)
t["excess_per_1000_V_fit_and_unit_all"] = (t.near - t.far / 2) / smc.n_v * 1000
t["excess_per_1000_V_fit_and_unit_small"] = (smc.near - smc.far / 2) / smc.n_v * 1000
t.reset_index().to_csv(f"{R}/overlap_by_v_size_summary.csv", index=False)
print("written; null 99th percentile", round(thr, 1))
print(pk.groupby("exp").apply(lambda z: pd.Series({"pairs": len(z), "small_peak": (z.peak_p_small < 0.01).mean(),
                                                   "large_thinned_peak": (z.peak_p_large_thinned < 0.01).mean(),
                                                   "shape_differs": (z.shape_p < 1e-3).mean()})).round(3).to_string())

# amplitude_cutoff (amp_cutoff_check.py's output, in the results folder)
Q = pd.read_csv(f"{R}/amp_cutoff_check.csv")
Q = Q[Q.stored.notna()]
Q["width_all"] = (Q.max_rel - Q.all_min_rel) / 100
Q["width_pos"] = (Q.max_rel - Q.pos_min_rel) / 100
Q["bins_above_p999_all"] = (Q.max_rel - Q.p999_rel) / (Q.max_rel - Q.all_min_rel)
Q["bins_above_p999_pos"] = (Q.max_rel - Q.p999_rel) / (Q.max_rel - Q.pos_min_rel)
Q["max_over_p999"] = Q.max_rel / Q.p999_rel
agg = Q.groupby("exp").agg(units=("unit", "size"), stored_median=("stored", "median"), stored_max=("stored", "max"),
                           all_vs_stored_max_diff=("all_cutoff", lambda v: np.nan),
                           lowest_bin_all=("all_lowest_bin", "median"), lowest_bin_pos=("pos_lowest_bin", "median"),
                           min_rel_all=("all_min_rel", "median"), max_rel_median=("max_rel", "median"),
                           max_rel_max=("max_rel", "max"), p999_rel=("p999_rel", "median"),
                           max_over_p999=("max_over_p999", "median"), width_all=("width_all", "median"),
                           width_pos=("width_pos", "median"), bins_above_p999_all=("bins_above_p999_all", "median"),
                           bins_above_p999_pos=("bins_above_p999_pos", "median"), pos_median=("pos_cutoff", "median"),
                           pos_max=("pos_cutoff", "max"))
agg["all_vs_stored_max_diff"] = Q.assign(x=(Q.all_cutoff - Q.stored).abs()).groupby("exp").x.max()
agg.reset_index().to_csv(f"{R}/amp_cutoff_summary.csv", index=False)
Qw = Q[(Q.tetrode >= 4) & (Q.tetrode <= 12)]
bands = ["<=0", "0-0.1", "0.1-0.25", "0.25-0.5", ">0.5"]
brow = []
for e, z in Qw.groupby("exp"):
    r = {"exp": e}
    for b in bands:
        r[f"n_{b}"] = int(z[f"n_{b}"].sum())
        r[f"other_unit_within_1ms_{b}"] = float((z[f"nbr_{b}"].fillna(0) * z[f"n_{b}"]).sum() / z[f"n_{b}"].sum())
    r["other_unit_within_1ms_all"] = float((z.nbr_all * z.spikes).sum() / z.spikes.sum())
    brow.append(r)
pd.DataFrame(brow).to_csv(f"{R}/amp_cutoff_bands.csv", index=False)
