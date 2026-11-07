"""Controls for refr_rejects.py: percentiles (within the unit's isolated accepted spikes) of the rejected
fits and of their survivors, for amplitude, shape cosine and snr; and the cosine of isolated accepted spikes
with another unit's accepted spike within 1.5 ms, the overlap control."""
import sys
import numpy as np
import pandas as pd
from ttm import carry as cr
RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = sys.argv[1]
fs = cr.open_sorting_recording(f"{RUN}/binary", ["CH22"]).get_sampling_frequency()
RF = int(1.5e-3 * fs)

def other_within(fq, uq, gq, fa, ua, ga, d):
    out = np.zeros(fq.size, bool)
    for g in np.unique(gq):
        qa, aa = np.flatnonzero(gq == g), np.flatnonzero(ga == g)
        F = fa[aa]
        n_all = np.searchsorted(F, fq[qa] + d, "right") - np.searchsorted(F, fq[qa] - d)
        n_own = np.zeros(qa.size, np.int64)
        for u in np.unique(uq[qa]):
            m = uq[qa] == u
            Fu = F[ua[aa] == u]
            n_own[m] = np.searchsorted(Fu, fq[qa][m] + d, "right") - np.searchsorted(Fu, fq[qa][m] - d)
        out[qa] = n_all > n_own
    return out

rows = []
for w in range(0, 96, 8):
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
    rej = (reason == 3) & (n0 == 1) & (size == 2)
    sv = acc & (size == 2)
    iso = acc & (size == 1)
    oa = np.argsort(f[acc], kind="stable")
    ov = np.zeros(f.size, bool)
    ov[iso] = other_within(f[iso], u[iso], g[iso], f[acc][oa], u[acc][oa], g[acc][oa], RF)
    for k in np.unique(u[rej]):
        mi = iso & (u == k) & ~ov
        if mi.sum() < 100 or (rej & (u == k)).sum() < 20:
            continue
        r = {"w": w, "unit": int(k), "n_rej": int((rej & (u == k)).sum())}
        for name, x in (("a", a), ("cos", cos), ("snr", snr)):
            ref = np.sort(x[mi])
            for lab, m in (("rej", rej), ("surv", sv), ("iso_overlap", iso & ov)):
                mm = m & (u == k)
                r[f"{name}_{lab}_pct"] = np.median(np.searchsorted(ref, x[mm]) / ref.size) if mm.any() else np.nan
        rows.append(r)
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/refr_rejects_controls.csv", index=False)
print(f"{len(df)} unit-windows, {df.n_rej.sum():,} rejected fits in pairs")
print("median percentile within the unit's isolated accepted spikes with no other unit within 1.5 ms")
print("(weighted by rejected count); two iid draws ordered by the measure would put the lower at 0.33, the higher at 0.67:")
for c in [c for c in df if c.endswith("_pct")]:
    m = df[c].notna()
    print(f"  {c:22s} {np.average(df[c][m], weights=df.n_rej[m]):.3f}")
