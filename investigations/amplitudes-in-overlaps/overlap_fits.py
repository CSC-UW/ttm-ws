"""Are the carry's accepted fits that sit within 1.5 ms of another unit's accepted fit real spikes?

Trace test, on the whitened traces with every fit subtracted except the event's own (all reasons), per
window in WINDOWS, N_BLOCKS blocks of BLOCK_S:
  event   an accepted fit of unit U (tested) and one accepted fit of another unit V of the tetrode at lag
          d = f_U - f_V, |d| <= 45 samples, and no other fit with snr >= 3 reaching the snippet
  twin    a lone accepted fit of V (no other fit with snr >= 3 within 9 ms) plus a lone accepted fit of U
          scaled by a_U(event) / a_U(lone), placed at lag d; both snippets are traces, noise included
  none    the same lone V snippet alone, U tested at lag d
Models, least squares on the 270-sample snippet (V's frame at sample 80), energies in whitened units:
  E_V     V with its first and second time derivatives at V's frame +-2 samples
  E_VU    the same plus U's template at U's frame +-2 samples
  dU      E_V - E_VU: what U explains beyond a rescaled, shifted, widened or narrowed V
  dU0     the same with V's plain template at V's frame and U's at U's frame, as the matcher fitted
  u_ls    U's coefficient in E_VU over a_U
Outputs ``overlap_events_<session>_wNN.csv``."""
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

TD = "/nvme/neuropixels/tetrode_data"
S = {"TTM-003-48H": ("2026-10-20_09-02-12", ["CH22"]), "TTM-003-48H-2": ("2026-11-03_09-01-38", ["A1_CH22"])}
OUT = sys.argv[1]
WINDOWS = [8, 20, 32, 44, 56, 68, 80, 92]
N_BLOCKS, BLOCK_S = 30, 10.0
NB, NA, TL = 30, 60, 90
PRE, LEN = 80, 270
DMAX = 45
CAP = 60
ISO_HALF = 270


def basis(T):
    d1 = np.gradient(T, axis=0)
    return np.stack([T, d1, np.gradient(d1, axis=0)])


def cross(A, B):
    """cross[L] = sum over overlap of A placed at 0 and B placed at L, L in [-TL+1, TL-1]."""
    out = np.zeros(2 * TL - 1)
    for i, L in enumerate(range(-TL + 1, TL)):
        if L >= 0:
            out[i] = (A[L:] * B[:TL - L]).sum()
        else:
            out[i] = (A[:TL + L] * B[-L:]).sum()
    return out


def fit(P, Bv, Tu, d):
    """Per snippet in P (n, LEN, C): E0, E_V, E_VU, dU0, u coefficient at the best E_VU."""
    n = P.shape[0]
    E0 = (P ** 2).sum(axis=(1, 2))
    Gvv = np.einsum("itc,jtc->ij", Bv, Bv)
    Guu = (Tu ** 2).sum()
    cr_ = np.stack([cross(Bv[k], Tu) for k in range(3)])  # (3, 2TL-1): Bv[k] at 0, Tu at L
    EV = np.full(n, np.inf)
    EVU = np.full(n, np.inf)
    ucoef = np.zeros(n)
    dU0 = np.zeros(n)
    for j in range(n):
        p = P[j]
        bv = {s: np.einsum("ktc,tc->k", Bv, p[PRE - NB + s: PRE - NB + s + TL]) for s in range(-2, 3)}
        bu = {t: (Tu * p[PRE - NB + d[j] + t: PRE - NB + d[j] + t + TL]).sum() for t in range(-2, 3)}
        for s in range(-2, 3):
            ev = E0[j] - bv[s] @ np.linalg.solve(Gvv, bv[s])
            EV[j] = min(EV[j], ev)
            for t in range(-2, 3):
                L = d[j] + t - s
                g = cr_[:, L + TL - 1] if abs(L) < TL else np.zeros(3)
                G = np.zeros((4, 4))
                G[:3, :3] = Gvv
                G[:3, 3] = G[3, :3] = g
                G[3, 3] = Guu
                b = np.r_[bv[s], bu[t]]
                c = np.linalg.solve(G, b)
                e = E0[j] - b @ c
                if e < EVU[j]:
                    EVU[j], ucoef[j] = e, c[3]
        # plain templates at the fitted frames
        L = d[j]
        g0 = cr_[0, L + TL - 1] if abs(L) < TL else 0.0
        G2 = np.array([[Gvv[0, 0], g0], [g0, Guu]])
        b2 = np.array([bv[0][0], bu[0]])
        e_v0 = E0[j] - bv[0][0] ** 2 / Gvv[0, 0]
        dU0[j] = e_v0 - (E0[j] - b2 @ np.linalg.solve(G2, b2))
    return E0, EV, EVU, ucoef, dU0


def run(job):
    exp, w = job
    from ttm import capture as cp
    from ttm import carry as cr

    ses, drop = S[exp]
    RUN = f"{TD}/{ses}/sortings_carry/main"
    rng = np.random.default_rng(w)
    rec = cr.open_sorting_recording(f"{RUN}/binary", drop)
    W = np.load(f"{RUN}/whitening.npz")["W"]
    wrec = cr.whitened_recording(rec, W)
    fs = rec.get_sampling_frequency()
    groups = np.asarray(rec.get_property("group")).astype(int)
    z = np.load(f"{RUN}/carry/spikes_w{w:02d}.npz")
    bk = np.load(f"{RUN}/carry/bankin_w{w:02d}.npz")
    ids, bgrp, Tall = bk["ids"], bk["group"], bk["T"]
    idx_of = {int(k): i for i, k in enumerate(ids.tolist())}
    f, u, a, snr, reason = (z[k] for k in ("frame", "unit", "a", "snr", "reason"))
    cos_raw, cos_pl = z["cos_raw"], z["cos_pl"]
    o = np.argsort(f, kind="stable")
    f, u, a, snr, reason, cos_raw, cos_pl = f[o], u[o], a[o], snr[o], reason[o], cos_raw[o], cos_pl[o]
    acc = reason == 0
    med_a = pd.Series(a[acc]).groupby(u[acc]).median().to_dict()
    bad = cr.read_bad_periods(f"{TD}/{ses}/bad_periods/bad_periods.csv", fs)
    win0 = w * int(1800 * fs)
    starts = np.linspace(win0, win0 + int(1800 * fs) - int((BLOCK_S + 1) * fs), N_BLOCKS).astype(np.int64)
    rows = []
    for b0 in starts:
        b1 = b0 + int(BLOCK_S * fs)
        lo, hi = b0 - 600, b1 + 600
        for g in np.unique(bgrp).tolist():
            chn = np.flatnonzero(groups == g)
            C = chn.size
            if C == 0 or any(s0 < hi and e0 > lo for s0, e0 in bad.get(g, [])):
                continue
            units_g = ids[bgrp == g]
            if units_g.size < 2:
                continue
            m_g = np.isin(u, units_g) & (f >= lo + NB) & (f < hi - NA)
            fg, ug, ag, sg, rg = f[m_g], u[m_g], a[m_g], snr[m_g], reason[m_g]
            crg, cpg = cos_raw[m_g], cos_pl[m_g]
            if fg.size == 0:
                continue
            x = cp.load_tetrode_traces(wrec, g, lo, hi).astype(np.float64)
            T = Tall[[idx_of[int(k)] for k in units_g]][:, :, :C].astype(np.float64)
            lab = np.searchsorted(units_g, ug)
            model = np.zeros_like(x)
            cr.add_spikes(model, (fg - lo).astype(np.int64), lab.astype(np.int64), ag.astype(np.float32),
                          T.astype(np.float32), NB)
            X = x - model
            strong = np.flatnonzero(sg >= 3)
            fs_strong = fg[strong]

            def n_strong(s0, s1):
                return np.searchsorted(fs_strong, s1) - np.searchsorted(fs_strong, s0)

            def snippet_with(idx_list, start):
                p = X[start - lo: start - lo + LEN].copy()
                for i in idx_list:
                    one = np.zeros((LEN + TL * 2, C))
                    cr.add_spikes(one, np.array([fg[i] - start + TL], np.int64), np.array([lab[i]], np.int64),
                                  np.array([ag[i]], np.float32), T.astype(np.float32), NB)
                    p += one[TL:TL + LEN]
                return p

            acc_i = np.flatnonzero((rg == 0) & (fg >= b0) & (fg < b1))
            # lone accepted fits per unit
            lone = {}
            for i in acc_i:
                if n_strong(fg[i] - ISO_HALF, fg[i] + ISO_HALF) == (1 if sg[i] >= 3 else 0):
                    lone.setdefault(int(lab[i]), []).append(i)
            acc_all = np.flatnonzero(rg == 0)
            fa = fg[acc_all]
            events = {}
            for i in acc_i:
                j0, j1 = np.searchsorted(fa, fg[i] - DMAX), np.searchsorted(fa, fg[i] + DMAX, "right")
                nb = [acc_all[j] for j in range(j0, j1) if acc_all[j] != i and lab[acc_all[j]] != lab[i]]
                if len(nb) != 1:
                    continue
                v = nb[0]
                start = fg[v] - PRE
                if start - lo < 0 or start + LEN - lo > X.shape[0]:
                    continue
                expect = int(sg[i] >= 3) + int(sg[v] >= 3)
                if n_strong(start - NA, start + LEN + NB) != expect:
                    continue
                events.setdefault((int(lab[i]), int(lab[v])), []).append((i, v))
            for (ki, kv), evs in events.items():
                if ki not in lone or kv not in lone:
                    continue
                if len(evs) > CAP:
                    evs = [evs[t] for t in rng.choice(len(evs), CAP, replace=False)]
                d = np.array([fg[i] - fg[v] for i, v in evs])
                P = np.stack([snippet_with([i, v], fg[v] - PRE) for i, v in evs])
                Pt, Pn, sc = [], [], []
                for (i, v), dd in zip(evs, d):
                    lv = lone[kv][rng.integers(len(lone[kv]))]
                    lu = lone[ki][rng.integers(len(lone[ki]))]
                    pv = snippet_with([lv], fg[lv] - PRE)
                    s_ = ag[i] / ag[lu]
                    pu = snippet_with([lu], fg[lu] - PRE - dd)
                    Pt.append(pv + s_ * pu)
                    Pn.append(pv)
                    sc.append(s_)
                Pt, Pn = np.stack(Pt), np.stack(Pn)
                Bv = basis(T[kv])
                Tu = T[ki]
                res = {name: fit(Q, Bv, Tu, d) for name, Q in (("obs", P), ("twin", Pt), ("none", Pn))}
                for j, (i, v) in enumerate(evs):
                    r = {"exp": exp, "w": w, "group": g, "U": int(units_g[ki]), "V": int(units_g[kv]),
                         "frame": int(fg[i]), "d": int(d[j]), "aU": float(ag[i]), "aV": float(ag[v]),
                         "aU_rel": float(ag[i] / med_a[int(units_g[ki])]), "aV_rel": float(ag[v] / med_a[int(units_g[kv])]),
                         "snrU": float(sg[i]), "snrV": float(sg[v]), "cos_rawU": float(crg[i]), "cos_plU": float(cpg[i]),
                         "twin_scale": float(sc[j]), "TuNorm2": float((Tu ** 2).sum())}
                    for name, (E0, EV, EVU, uc, dU0) in res.items():
                        r[f"{name}_E0"] = E0[j]
                        r[f"{name}_dU"] = EV[j] - EVU[j]
                        r[f"{name}_dU0"] = dU0[j]
                        r[f"{name}_u_ls"] = uc[j] / ag[i]
                    rows.append(r)
    d = pd.DataFrame(rows)
    d.to_csv(f"{OUT}/overlap_events_{S[exp][0]}_w{w:02d}.csv", index=False)
    return exp, w, len(d)


if __name__ == "__main__":
    jobs = [(e, w) for e in S for w in WINDOWS]
    with ProcessPoolExecutor(16) as ex:
        for r in ex.map(run, jobs):
            print(*r, flush=True)
    print("all done", flush=True)
