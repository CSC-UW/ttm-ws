"""Tables and numbers for the README from the outputs of refr_rejects.py, refr_rejects2.py, refr_acg.py and
refr_pairs.py (argv[1] holds the bulk outputs, argv[2] the results folder)."""
import sys

import numpy as np
import pandas as pd

S, R = sys.argv[1], sys.argv[2]
FS, WIN_S, RF = 30000.0, 1800.0, 45

lag = np.abs(np.load(f"{S}/refr_rejects_lag.npy"))
h = np.bincount(lag, minlength=RF + 1)[:RF + 1]
pd.DataFrame({"samples": np.arange(RF + 1), "ms": np.arange(RF + 1) / FS * 1e3, "n": h,
              "share": h / lag.size}).to_csv(f"{R}/lag_hist.csv", index=False)
print(f"rejected fits {lag.size:,}; |lag| median {np.median(lag) / 30:.2f} ms; under 0.3 ms {np.mean(lag < 9):.3f}; "
      f"under 0.1 ms {np.mean(lag < 3):.4f}")

pr = pd.read_csv(f"{S}/refr_rejects_pairs.csv.gz")
pr["d"] = pr.lag.abs()
by = (pr.groupby(pd.cut(pr.d, [0, 3, 6, 9, 12, 15, 21, 30, 46], right=False), observed=True)
      .agg(n=("d", "size"), a_surv=("a_surv", "median"), a_rej=("a_rej", "median"), cos_rej=("cos_rej", "median")))
by.index = [f"{iv.left / 30:.1f}-{iv.right / 30:.2f} ms" for iv in by.index]
by.to_csv(f"{R}/pairs_by_lag.csv", index_label="lag")
print("\nchains of two, by lag:\n" + by.round(2).to_string())

q = pr[pr.d.between(9, RF)].copy()
q["first_larger"] = np.where(q.lag > 0, q.a_surv > q.a_rej, q.a_rej > q.a_surv)
rows = []
for u, g in q.groupby("unit"):
    if len(g) < 300:
        continue
    c = np.bincount(g.d, minlength=RF + 2)
    mode = int(np.argmax(c))
    rows.append({"unit": u, "n": len(g), "mode_ms": mode / 30, "share_mode_pm1": c[mode - 1:mode + 2].sum() / len(g),
                 "first_larger": g.first_larger.mean()})
lm = pd.DataFrame(rows)
lm.to_csv(f"{R}/lag_mode_per_unit.csv", index=False)
print(f"\nunits with 300+ pairs at 0.3-1.5 ms: {len(lm)}; share within 1 sample of the mode median "
      f"{lm.share_mode_pm1.median():.3f} (10-90 %: {lm.share_mode_pm1.quantile(.1):.3f}-"
      f"{lm.share_mode_pm1.quantile(.9):.3f}; even spread {3 / (RF - 9 + 1):.3f}); highest:")
print(lm.sort_values("share_mode_pm1", ascending=False).head(3).round(3).to_string(index=False))
print(f"first fit larger: median {lm.first_larger.median():.3f}, range {lm.first_larger.min():.3f}-"
      f"{lm.first_larger.max():.3f}")

un = pd.read_csv(f"{R}/refr_rejects_units.csv")
n_u = un.n_iso + un.n_surv
contam = un.n_rej * WIN_S / (2 * RF / FS * n_u ** 2)
print(f"\nindependent contaminant needed, as a multiple of the unit's accepted count: median {contam.median():.2f}; "
      f"share of rejected fits in unit-windows needing more than 1: "
      f"{un.n_rej[contam > 1].sum() / un.n_rej.sum():.2f}")
print(f"another unit's accepted fit within 0.5 ms: rejected {np.average(un.near_other_rej, weights=un.n_rej):.3f}, "
      f"isolated accepted {np.average(un.near_other_iso, weights=un.n_iso):.3f} (weighted by count)")

co = pd.read_csv(f"{R}/refr_rejects_controls.csv")
print("\npercentiles within the unit's isolated accepted spikes, median over unit-windows:")
for m in ("a", "snr", "cos"):
    print(f"  {m:4s} rejected {co[m + '_rej_pct'].median():.2f}  survivor {co[m + '_surv_pct'].median():.2f}  "
          f"overlap control {co[m + '_iso_overlap_pct'].median():.2f}")

acg = pd.read_csv(f"{R}/refr_acg.csv")
bins = [f"{6 * i / 30:.1f}-{6 * (i + 1) / 30:.1f}" for i in range(25)]
t = acg.groupby("kind")[[str(i) for i in range(25)]].median()
t.columns = bins
print(f"\nACG over its 50-100 ms level, median over {acg.groupby(['w', 'unit']).ngroups} unit-windows "
      f"({acg.unit.nunique()} units), 0.2 ms bins:")
print(t.iloc[:, :10].round(2).to_string())
pa = pd.read_csv(f"{R}/refr_pairs_acg.csv")
print("\naccepted ACG 1.5-8 ms, median and quartiles over unit-windows:")
print(pa.drop(columns=["w", "unit"]).describe().loc[["50%", "25%", "75%"]].round(2).to_string())
pf = pd.read_csv(f"{R}/refr_pairs_features.csv")
print("\naccepted pairs, percentiles of the first and second spike, median over unit-windows:")
print(pf.groupby(["bin", "measure"])[["first", "second"]].median().unstack().round(2).to_string())
