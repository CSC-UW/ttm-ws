"""Fits the carry rejected under the refractory rule (reason 3) on TTM-003-48H, against criteria the rule
does not use: lag to the chain's survivor, fitted amplitude, shape cosine, and coincidence with another
unit's accepted spike on the same tetrode. Every 8th window."""
import sys

import numpy as np
import pandas as pd

from ttm import carry as cr

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = sys.argv[1]
WINDOWS = range(0, 96, 8)
rec = cr.open_sorting_recording(f"{RUN}/binary", ["CH22"])
fs = rec.get_sampling_frequency()
RF = int(1.5e-3 * fs)
CO = int(0.5e-3 * fs)


def others_near(fq, uq, gq, fa, ua, ga):
    """Per query fit, whether an accepted fit of another unit on its tetrode lies within CO samples."""
    out = np.zeros(fq.size, bool)
    for g in np.unique(gq):
        qa = np.flatnonzero(gq == g)
        aa = np.flatnonzero(ga == g)
        F = fa[aa]
        n_all = np.searchsorted(F, fq[qa] + CO, "right") - np.searchsorted(F, fq[qa] - CO)
        n_own = np.zeros(qa.size, np.int64)
        for u in np.unique(uq[qa]):
            m = uq[qa] == u
            Fu = F[ua[aa] == u]
            n_own[m] = np.searchsorted(Fu, fq[qa][m] + CO, "right") - np.searchsorted(Fu, fq[qa][m] - CO)
        out[qa] = n_all > n_own
    return out


rows, lag_all, pairs = [], [], []
for w in WINDOWS:
    z = np.load(f"{RUN}/carry/spikes_w{w:02d}.npz")
    bk = np.load(f"{RUN}/carry/bankin_w{w:02d}.npz")
    gmap = dict(zip(bk["ids"].tolist(), bk["group"].tolist()))
    keep = np.isin(z["reason"], (0, 3))
    f, u, a, cos, snr, reason = (z[k][keep] for k in ("frame", "unit", "a", "cos_pl", "snr", "reason"))
    g = np.array([gmap[x] for x in u.tolist()])
    o = np.lexsort((f, u))
    f, u, a, cos, snr, reason, g = f[o], u[o], a[o], cos[o], snr[o], reason[o], g[o]
    close = np.r_[False, (np.diff(f) < RF) & (u[1:] == u[:-1])]
    chain = np.cumsum(~close) - 1
    size = np.bincount(chain)[chain]
    n0 = np.bincount(chain, weights=reason == 0)[chain]
    acc = reason == 0
    surv_f = np.full(chain.max() + 1, -1, np.int64)
    surv_a = np.full(chain.max() + 1, np.nan)
    surv_f[chain[acc]] = f[acc]
    surv_a[chain[acc]] = a[acc]
    rej = (reason == 3) & (n0 == 1)
    iso = acc & (size == 1)
    sv = acc & (size > 1)
    lag = f[rej] - surv_f[chain[rej]]
    lag_all.append(lag)
    oa = np.argsort(f[acc], kind="stable")
    fa, ua, ga = f[acc][oa], u[acc][oa], g[acc][oa]
    near_rej = others_near(f[rej], u[rej], g[rej], fa, ua, ga)
    near_iso = others_near(f[iso], u[iso], g[iso], fa, ua, ga)
    two = rej & (size == 2)
    pairs.append(pd.DataFrame({"w": w, "unit": u[two], "lag": f[two] - surv_f[chain[two]], "a_rej": a[two],
                               "a_surv": surv_a[chain[two]], "cos_rej": cos[two]}))
    for k in np.unique(u[rej]):
        mi, mr, ms = iso & (u == k), rej & (u == k), sv & (u == k)
        if mi.sum() < 100:
            continue
        ai = np.sort(a[mi])
        ci = np.sort(cos[mi])
        rows.append({"w": w, "unit": int(k), "group": gmap[int(k)], "n_iso": int(mi.sum()), "n_surv": int(ms.sum()),
                     "n_rej": int(mr.sum()), "a_iso": np.median(a[mi]), "a_surv": np.median(a[ms]),
                     "a_rej": np.median(a[mr]), "a_rej_pct": np.median(np.searchsorted(ai, a[mr]) / ai.size),
                     "cos_iso": np.median(cos[mi]), "cos_rej": np.median(cos[mr]),
                     "cos_rej_pct": np.median(np.searchsorted(ci, cos[mr]) / ci.size),
                     "snr_iso": np.median(snr[mi]), "snr_rej": np.median(snr[mr]),
                     "near_other_rej": near_rej[u[rej] == k].mean(), "near_other_iso": near_iso[u[iso] == k].mean()})
    print(f"w{w:02d}: {rej.sum():,} rejected, {iso.sum():,} isolated accepted, {sv.sum():,} survivors; "
          f"chains >2: {(np.bincount(chain) > 2).sum():,}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/refr_rejects_units.csv", index=False)
pr = pd.concat(pairs)
pr.to_csv(f"{OUT}/refr_rejects_pairs.csv.gz", index=False)
lag = np.concatenate(lag_all)
np.save(f"{OUT}/refr_rejects_lag.npy", lag)
ms = np.abs(lag) / fs * 1e3
print(f"\nrejected fits: {lag.size:,}; |lag| to survivor (ms) quantiles 10/25/50/75/90:",
      np.round(np.percentile(ms, [10, 25, 50, 75, 90]), 3))
edges = np.arange(0, RF + 3, 3)
h = np.histogram(np.abs(lag), edges)[0]
print("|lag| histogram, 0.1 ms bins (share of rejected):")
for e0, n in zip(edges[:-1], h):
    print(f"  {e0 / fs * 1e3:4.1f}-{(e0 + 3) / fs * 1e3:4.1f} ms  {n / lag.size:6.3f}")
print("share after the survivor:", np.mean(lag > 0).round(3))
w = df.n_rej
for c in ("a_iso", "a_surv", "a_rej", "a_rej_pct", "cos_iso", "cos_rej", "cos_rej_pct", "snr_iso", "snr_rej",
          "near_other_rej", "near_other_iso"):
    print(f"{c:15s} median over unit-windows {df[c].median():.3f}; weighted by rejected {np.average(df[c], weights=w):.3f}")
pr["ms"] = pr.lag.abs() / fs * 1e3
pr["sum_rel"] = (pr.a_rej + pr.a_surv)
for lo, hi in ((0, 0.3), (0.3, 0.7), (0.7, 1.1), (1.1, 1.5)):
    m = (pr.ms >= lo) & (pr.ms < hi)
    print(f"pairs {lo}-{hi} ms: n {m.sum():,}; a_surv median {pr.a_surv[m].median():.2f}, a_rej {pr.a_rej[m].median():.2f}, "
          f"sum {pr.sum_rel[m].median():.2f}, cos_rej {pr.cos_rej[m].median():.2f}")
