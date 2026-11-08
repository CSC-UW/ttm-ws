"""Timing and size of accepted fits of unit U around accepted fits of another unit V of the same tetrode,
over every window of both sessions.

Per ordered pair (U, V) and U-size class (a over U's median a in the window): the histogram of U's lags to
V's fits within 5 ms in 1-sample bins. Per V-size class (V's a over V's median): U fits per V fit at lags
of 1.5 ms or less (near, 91 samples) and at 2-5 ms either side (far, 182 samples), by U-size class."""
import sys

import numpy as np
import pandas as pd

TD = "/nvme/neuropixels/tetrode_data"
S = {"TTM-003-48H": "2026-10-20_09-02-12", "TTM-003-48H-2": "2026-11-03_09-01-38"}
OUT = sys.argv[1]
M = 150
U_EDGES = [0, 0.7, 1.0, np.inf]
U_NAMES = ["small", "mid", "large"]
V_EDGES = [0, 0.8, 0.95, 1.05, 1.2, np.inf]
hist_rows, cond_rows = [], []
for exp, ses in S.items():
    run = f"{TD}/{ses}/sortings_carry/main"
    H = {}
    Cn = {}
    for w in range(96):
        z = np.load(f"{run}/carry/spikes_w{w:02d}.npz")
        bk = np.load(f"{run}/carry/bankin_w{w:02d}.npz")
        gmap = dict(zip(bk["ids"].tolist(), bk["group"].tolist()))
        m = z["reason"] == 0
        f, u, a = z["frame"][m], z["unit"][m], z["a"][m]
        o = np.argsort(f, kind="stable")
        f, u, a = f[o], u[o], a[o]
        med = pd.Series(a).groupby(u).median()
        rel = a / med.reindex(u).to_numpy()
        g = np.array([gmap[int(k)] for k in u])
        for gg in np.unique(g):
            mg = g == gg
            fg, ug, rg = f[mg], u[mg], rel[mg]
            units = np.unique(ug)
            for U in units:
                iu = np.flatnonzero(ug == U)
                cu = np.digitize(rg[iu], U_EDGES) - 1
                for V in units:
                    if V == U:
                        continue
                    iv = np.flatnonzero(ug == V)
                    fv = fg[iv]
                    cv = np.digitize(rg[iv], V_EDGES) - 1
                    key = (int(gg), int(U), int(V))
                    h = H.setdefault(key, np.zeros((3, 2 * M + 1), np.int64))
                    c = Cn.setdefault(key, np.zeros((len(V_EDGES) - 1, 3, 2), np.int64))
                    lo = np.searchsorted(fv, fg[iu] - M)
                    hi = np.searchsorted(fv, fg[iu] + M, "right")
                    cnt = hi - lo
                    tot = int(cnt.sum())
                    if tot:
                        iu_rep = np.repeat(np.arange(iu.size), cnt)
                        vidx = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt) + np.repeat(lo, cnt)
                        lags = fg[iu][iu_rep] - fv[vidx]
                        cc = cu[iu_rep]
                        h += np.bincount(cc * (2 * M + 1) + lags + M, minlength=3 * (2 * M + 1)).reshape(3, -1)
                        cvs = cv[vidx]
                        nbin = len(V_EDGES) - 1
                        near = np.abs(lags) <= 45
                        far = np.abs(lags) >= 60
                        c[:, :, 0] += np.bincount(cvs[near] * 3 + cc[near], minlength=nbin * 3).reshape(nbin, 3)
                        c[:, :, 1] += np.bincount(cvs[far] * 3 + cc[far], minlength=nbin * 3).reshape(nbin, 3)
                    nv = np.bincount(cv, minlength=len(V_EDGES) - 1)
                    Cn.setdefault(key + ("nv",), np.zeros(len(V_EDGES) - 1, np.int64))
                    Cn[key + ("nv",)] += nv
        if w % 12 == 0:
            print(exp, w, flush=True)
    for key, h in H.items():
        for ci, cname in enumerate(U_NAMES):
            hist_rows.append({"exp": exp, "group": key[0], "U": key[1], "V": key[2], "u_class": cname,
                              **{f"l{L}": int(h[ci, L + M]) for L in range(-M, M + 1)}})
        c = Cn[key]
        nv = Cn[key + ("nv",)]
        for vi in range(len(V_EDGES) - 1):
            for ci, cname in enumerate(U_NAMES):
                cond_rows.append({"exp": exp, "group": key[0], "U": key[1], "V": key[2], "v_bin": vi,
                                  "u_class": cname, "n_v": int(nv[vi]), "near": int(c[vi, ci, 0]),
                                  "far": int(c[vi, ci, 1])})
pd.DataFrame(hist_rows).to_csv(f"{OUT}/overlap_lag_hist.csv", index=False)
pd.DataFrame(cond_rows).to_csv(f"{OUT}/overlap_by_v_size.csv", index=False)
print("all done", flush=True)
