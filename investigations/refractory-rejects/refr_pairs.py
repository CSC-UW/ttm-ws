"""Accepted pairs of one unit at lags 1.5-8 ms (every 8th window of TTM-003-48H): the ACG to 8 ms, and the
first and second spike's amplitude, cosine and snr as percentiles within the unit's accepted spikes with no
same-unit accepted spike within 10 ms."""
import sys
import numpy as np
import pandas as pd
RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = sys.argv[1]
fs = 30000.0
bins = [(1.5, 2.0), (2.0, 2.6), (2.6, 3.4), (3.4, 4.2), (4.2, 5.0), (5.0, 6.0), (6.0, 8.0)]
acc_rows, feat = [], []
for w in range(0, 96, 8):
    z = np.load(f"{RUN}/carry/spikes_w{w:02d}.npz")
    m = z["reason"] == 0
    f, u, a, cos, snr = (z[k][m] for k in ("frame", "unit", "a", "cos_pl", "snr"))
    for k in np.unique(u):
        i = np.flatnonzero(u == k)
        if i.size < 2000:
            continue
        i = i[np.argsort(f[i], kind="stable")]
        F = f[i]
        dprev = np.r_[np.inf, np.diff(F)]
        dnext = np.r_[np.diff(F), np.inf]
        iso = (dprev > 300) & (dnext > 300)
        refs = {n: np.sort(x[i][iso]) for n, x in (("a", a), ("cos", cos), ("snr", snr))}
        base = (np.searchsorted(F, F + 3000) - np.searchsorted(F, F + 1500)).sum() / 50.0
        r = {"w": w, "unit": int(k)}
        for lo, hi in bins:
            sel = np.flatnonzero((dnext >= lo * 30) & (dnext < hi * 30) & (dprev > 300))
            sel2 = sel[np.r_[np.diff(F), np.inf][sel + 1 if sel.size and sel[-1] + 1 < F.size else sel] > 300] if False else sel
            cnt = ((np.searchsorted(F, F + hi * 30, "left") - np.searchsorted(F, F + lo * 30, "left")).sum())
            r[f"{lo}-{hi}"] = cnt / (hi - lo) / base
            if sel.size >= 20:
                for n in ("a", "cos", "snr"):
                    x = {"a": a, "cos": cos, "snr": snr}[n][i]
                    feat.append({"w": w, "unit": int(k), "bin": f"{lo}-{hi}", "measure": n, "n": sel.size,
                                 "first": np.median(np.searchsorted(refs[n], x[sel]) / refs[n].size),
                                 "second": np.median(np.searchsorted(refs[n], x[sel + 1]) / refs[n].size)})
        acc_rows.append(r)
A = pd.DataFrame(acc_rows)
Fe = pd.DataFrame(feat)
A.to_csv(f"{OUT}/refr_pairs_acg.csv", index=False)
Fe.to_csv(f"{OUT}/refr_pairs_features.csv", index=False)
print("accepted ACG, median over unit-windows [25-75 %]:")
for lo, hi in bins:
    c = A[f"{lo}-{hi}"]
    print(f"  {lo}-{hi} ms  {c.median():.2f} [{c.quantile(.25):.2f}-{c.quantile(.75):.2f}]")
print("\npairs (first spike with no same-unit spike 10 ms before): percentile within isolated spikes, weighted by pairs")
for b in [f"{lo}-{hi}" for lo, hi in bins]:
    s = []
    for n in ("a", "cos", "snr"):
        d = Fe[(Fe.bin == b) & (Fe.measure == n)]
        s.append(f"{n} {np.average(d['first'], weights=d.n):.2f}/{np.average(d['second'], weights=d.n):.2f}")
    print(f"  {b:8s} first/second: " + "   ".join(s) + f"   (pairs {Fe[(Fe.bin == b) & (Fe.measure == 'a')].n.sum():,})")
