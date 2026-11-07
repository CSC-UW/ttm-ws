"""Trace check of the carry's close same-unit pairs on TTM-003-48H.

Events, on the whitened traces with every other unit's fits subtracted (all reasons):
  rej   a refractory-rejected fit and its survivor (chain of exactly two), lag 0.3-1.5 ms
  far   two accepted fits of the unit 2.0-2.6 ms apart (two separate spikes, control)
  iso   an accepted fit with no accepted or rejected fit of any unit on the tetrode within 3 ms
No other fit of the unit with snr >= 3 may lie in an event's window. Each event gets a synthetic twin:
the matcher's own explanation (its fits' templates at their frames and amplitudes) added to a window of
the same traces, same unit, with no fit of the unit with snr >= 3.

Models fitted by least squares on each window (270 samples, anchor = first fit at sample 80):
  E1    one copy of the unit's template, best lag
  E1d   one copy plus its first and second time derivatives (shifted, broadened or narrowed), best lag
  E1xd  as E1d with the best template of any unit on the tetrode
  E2    two copies at the fitted frames, each refined by up to 2 samples
  E2f   two copies, the first within 2 samples of the anchor, the second 0.3-1.5 ms from it, best pair
Energies are in whitened-noise units (sum of squares)."""
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
BAD = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
OUT = sys.argv[1]
WINDOWS = [8, 20, 32, 44, 56, 68, 80, 92]
N_BLOCKS, BLOCK_S = 30, 10.0
NB, NA, TL = 30, 60, 90
PRE, LEN = 80, 270
TAUS = np.arange(LEN - TL + 1) - (PRE - NB)  # fit frames whose template lies inside the window
RF = 45


def chains(f, u, reason):
    """Per fit (sorted by unit then frame): size of its refractory chain among the fits with reason 0 or 3
    and the chain's survivor frame (0 and -1 for other fits)."""
    m = np.flatnonzero(np.isin(reason, (0, 3)))
    fm, um, rm = f[m], u[m], reason[m]
    close = np.r_[False, (np.diff(fm) < RF) & (um[1:] == um[:-1])]
    ch = np.cumsum(~close) - 1
    sf = np.full(ch.max() + 1, -1, np.int64)
    sf[ch[rm == 0]] = fm[rm == 0]
    size = np.zeros(f.size, np.int64)
    surv = np.full(f.size, -1, np.int64)
    size[m] = np.bincount(ch)[ch]
    surv[m] = sf[ch]
    return size, surv


def basis(T):
    d1 = np.gradient(T, axis=0)
    return np.stack([T, d1, np.gradient(d1, axis=0)])


def fit_models(P, Tu, Bx, ku, d, iso):
    """Energies per event. P (n, LEN, C) windows; Tu (TL, C) the unit's template; Bx (K, 3, TL, C) every
    unit's basis on the tetrode; d (n,) the second fit's offset from the anchor (ignored for iso)."""
    n = P.shape[0]
    E0 = (P ** 2).sum(axis=(1, 2))
    view = np.lib.stride_tricks.sliding_window_view(P, TL, axis=1)  # (n, S, C, TL)
    c = np.einsum("nsct,kjtc->nkjs", view, Bx, optimize=True)  # (n, K, 3, S)
    G = np.einsum("kitc,kjtc->kij", Bx, Bx)
    Gi = np.linalg.inv(G)
    gain3 = np.einsum("nkis,kij,nkjs->nks", c, Gi, c)  # (n, K, S)
    out = {"E0": E0}
    cu = c[:, ku, 0]
    n2 = (Tu ** 2).sum()
    out["E1"] = E0 - (cu ** 2).max(axis=1) / n2
    out["tau1"] = TAUS[(cu ** 2).argmax(axis=1)]
    out["E1d"] = E0 - gain3[:, ku].max(axis=1)
    out["E1xd"] = E0 - gain3.reshape(n, -1).max(axis=1)
    out["k1xd"] = gain3.max(axis=2).argmax(axis=1) == ku
    lags = np.arange(-TL + 1, TL)
    ac = np.array([(Tu[max(0, -L):TL - max(0, L)] * Tu[max(0, L):TL - max(0, -L)]).sum() for L in lags])

    def two(t1, t2):
        i1, i2 = t1 + 50, t2 + 50
        ok = (i1 >= 0) & (i1 < TAUS.size) & (i2 >= 0) & (i2 < TAUS.size)
        i1c, i2c = np.clip(i1, 0, TAUS.size - 1), np.clip(i2, 0, TAUS.size - 1)
        b1, b2 = cu[np.arange(n), i1c], cu[np.arange(n), i2c]
        r = ac[np.clip(t2 - t1, -TL + 1, TL - 1) + TL - 1]
        det = n2 * n2 - r * r
        g = (n2 * b1 ** 2 - 2 * r * b1 * b2 + n2 * b2 ** 2) / det
        a1 = (n2 * b1 - r * b2) / det
        a2 = (n2 * b2 - r * b1) / det
        return np.where(ok, E0 - g, np.inf), a1, a2

    best = np.full(n, np.inf)
    a1b, a2b = np.zeros(n), np.zeros(n)
    s1b, s2b = np.zeros(n, np.int64), np.zeros(n, np.int64)
    if not iso:
        for s1 in range(-2, 3):
            for s2 in range(-2, 3):
                e, a1, a2 = two(np.full(n, s1), d + s2)
                m = e < best
                best[m], a1b[m], a2b[m], s1b[m], s2b[m] = e[m], a1[m], a2[m], s1, d[m] + s2
    out["E2"], out["a1_ls"], out["a2_ls"] = best, a1b, a2b
    bf = np.full(n, np.inf)
    for s1 in range(-2, 3):
        for s in list(range(-45, -8)) + list(range(9, 46)):
            e, _, _ = two(np.full(n, s1), np.full(n, s1 + s))
            bf = np.minimum(bf, e)
    out["E2f"] = bf
    # residual energy over the fitted spikes' span against the rest of the window, per sample
    inside = np.zeros((n, LEN), bool)
    for j in range(n):
        inside[j, PRE - NB: PRE + NA + (0 if iso else d[j])] = True
    sd = gain3[:, ku].argmax(axis=1)
    coef = np.stack([Gi[ku] @ c[j, ku, :, sd[j]] for j in range(n)])
    R1d = P.copy()
    Rm = P.copy()
    mod = np.zeros(n)
    for j in range(n):
        R1d[j, sd[j]: sd[j] + TL] -= np.tensordot(coef[j], Bx[ku], axes=1)
        if iso:
            s = (cu[j] ** 2).argmax()
            amp = cu[j, s] / n2
            Rm[j, s: s + TL] -= amp * Tu
            mod[j] = amp ** 2 * n2
        else:
            i1, i2 = s1b[j] + PRE - NB, s2b[j] + PRE - NB
            Rm[j, i1: i1 + TL] -= a1b[j] * Tu
            Rm[j, i2: i2 + TL] -= a2b[j] * Tu
            mod[j] = (a1b[j] ** 2 + a2b[j] ** 2) * n2 + 2 * a1b[j] * a2b[j] * ac[s2b[j] - s1b[j] + TL - 1]
    C = P.shape[2]
    for name, R in (("m", Rm), ("m1d", R1d)):
        e = (R ** 2).sum(axis=2)
        e_in = (e * inside).sum(axis=1)
        e_out = (e * ~inside).sum(axis=1) / (~inside).sum(axis=1)
        out[f"misfit_{name}"] = e_in - e_out * inside.sum(axis=1)
        out["noise_out"] = e_out / C
    out["model_energy"] = mod
    return out


def run_window(w):
    from ttm import capture as cp
    from ttm import carry as cr

    rng = np.random.default_rng(w)
    rec = cr.open_sorting_recording(f"{RUN}/binary", ["CH22"])
    W = np.load(f"{RUN}/whitening.npz")["W"]
    wrec = cr.whitened_recording(rec, W)
    fs = rec.get_sampling_frequency()
    groups = np.asarray(rec.get_property("group")).astype(int)
    gain = rec.get_channel_gains()
    z = np.load(f"{RUN}/carry/spikes_w{w:02d}.npz")
    bk = np.load(f"{RUN}/carry/bankin_w{w:02d}.npz")
    ids, bgrp, Tall = bk["ids"], bk["group"], bk["T"]
    idx_of = {int(k): i for i, k in enumerate(ids.tolist())}
    f, u, a, snr, reason = (z[k] for k in ("frame", "unit", "a", "snr", "reason"))
    o = np.lexsort((f, u))
    f, u, a, snr, reason = f[o], u[o], a[o], snr[o], reason[o]
    size, sf = chains(f, u, reason)
    bad = cr.read_bad_periods(BAD, fs)
    win0 = w * int(1800 * fs)
    starts = np.linspace(win0, win0 + int(1800 * fs) - int((BLOCK_S + 1) * fs), N_BLOCKS).astype(np.int64)
    rows, avgs = [], []
    for b0 in starts:
        b1 = b0 + int(BLOCK_S * fs)
        lo, hi = b0 - 400, b1 + 400
        for g in np.unique(bgrp).tolist():
            chn = np.flatnonzero(groups == g)
            C = chn.size
            if C == 0 or any(s < hi and e > lo for s, e in bad.get(g, [])):
                continue
            x = cp.load_tetrode_traces(wrec, g, lo, hi).astype(np.float64)
            units_g = ids[bgrp == g]
            m_g = np.isin(u, units_g) & (f >= lo + NB) & (f < hi - NA)
            fg, ug, ag, sg, rg = f[m_g], u[m_g], a[m_g], snr[m_g], reason[m_g]
            szg, sfg = size[m_g], sf[m_g]
            T = Tall[[idx_of[int(k)] for k in units_g]][:, :, :C].astype(np.float64)
            Bx = np.stack([basis(t) for t in T])
            model = np.zeros_like(x)
            lab = np.searchsorted(units_g, ug)
            cr.add_spikes(model, (fg - lo).astype(np.int64), lab.astype(np.int64), ag.astype(np.float32),
                          T.astype(np.float32), NB)
            Winv = np.linalg.inv(W[np.ix_(chn, chn)].astype(np.float64))
            acc_any = np.sort(fg[np.isin(rg, (0, 3))])
            for ki, k in enumerate(units_g.tolist()):
                mu = ug == k
                if (mu & (rg == 0)).sum() < 20:
                    continue
                own = np.zeros_like(x)
                cr.add_spikes(own, (fg[mu] - lo).astype(np.int64), np.full(mu.sum(), ki, np.int64),
                              ag[mu].astype(np.float32), T.astype(np.float32), NB)
                X = x - model + own
                fu, ru, au, su, szu, sfu = fg[mu], rg[mu], ag[mu], sg[mu], szg[mu], sfg[mu]
                strong = np.sort(fu[su >= 3])
                inb = (fu >= b0) & (fu < b1)
                ev = []
                rj = np.flatnonzero(inb & (ru == 3) & (szu == 2) & (sfu >= 0))
                for i in rj:
                    s = sfu[i]
                    first, second = min(fu[i], s), max(fu[i], s)
                    d = second - first
                    if d < 9:
                        continue
                    a1 = au[i] if fu[i] == first else au[np.searchsorted(fu, s)]
                    a2 = au[i] if fu[i] == second else au[np.searchsorted(fu, s)]
                    ev.append(("rej", first, d, a1, a2, int(fu[i] == first)))
                acc_i = np.flatnonzero(ru == 0)
                fa = fu[acc_i]
                dn = np.diff(fa)
                for j in np.flatnonzero((dn >= 60) & (dn <= 78)):
                    if fa[j] < b0 or fa[j] >= b1:
                        continue
                    ev.append(("far", fa[j], int(dn[j]), au[acc_i[j]], au[acc_i[j + 1]], -1))
                cand = np.flatnonzero(inb & (ru == 0))
                rng.shuffle(cand)
                n_iso = 0
                for i in cand:
                    if n_iso >= 30:
                        break
                    q = fu[i]
                    if (np.searchsorted(acc_any, q + 90, "right") - np.searchsorted(acc_any, q - 90)) > 1:
                        continue
                    ev.append(("iso", q, 0, au[i], 0.0, -1))
                    n_iso += 1
                keep = []
                for e in ev:
                    s0, s1 = e[1] - PRE - NA, e[1] - PRE + LEN + NB
                    n_strong = np.searchsorted(strong, s1) - np.searchsorted(strong, s0)
                    expect = sum(su[np.searchsorted(fu, q)] >= 3 for q in ([e[1], e[1] + e[2]] if e[0] != "iso" else [e[1]]))
                    if n_strong == expect and e[1] - PRE - lo >= 0 and e[1] - PRE + LEN - lo <= X.shape[0]:
                        keep.append(e)
                if not keep:
                    continue
                P = np.stack([X[e[1] - PRE - lo: e[1] - PRE + LEN - lo] for e in keep])
                d = np.array([e[2] for e in keep])
                Psyn = np.empty_like(P)
                tries = 0
                for j, e in enumerate(keep):
                    while True:
                        tries += 1
                        q = int(rng.integers(b0, b1 - LEN))
                        if np.searchsorted(strong, q + LEN + NB) - np.searchsorted(strong, q - NA) == 0:
                            break
                        if tries > 100000:
                            raise RuntimeError("no noise window")
                    nz = X[q - lo: q - lo + LEN].copy()
                    nz[PRE - NB: PRE + NA] += e[3] * T[ki]
                    if e[0] != "iso":
                        nz[PRE - NB + e[2]: PRE + NA + e[2]] += e[4] * T[ki]
                    Psyn[j] = nz
                for kind in ("rej", "far", "iso"):
                    sel = np.array([e[0] == kind for e in keep])
                    if not sel.any():
                        continue
                    real = fit_models(P[sel], T[ki], Bx, ki, d[sel], kind == "iso")
                    syn = fit_models(Psyn[sel], T[ki], Bx, ki, d[sel], kind == "iso")
                    kk = [e for e in keep if e[0] == kind]
                    for j, e in enumerate(kk):
                        r = {"w": w, "group": g, "unit": k, "kind": kind, "frame": e[1], "d": e[2], "a1": e[3],
                             "a2": e[4], "rej_first": e[5]}
                        r.update({f"{m}": v[j] for m, v in real.items()})
                        r.update({f"syn_{m}": v[j] for m, v in syn.items()})
                        rows.append(r)
                    mainc = int(np.argmax(np.ptp(T[ki] @ Winv, axis=0)))
                    plain = P[sel] @ Winv * gain[chn][None, None, :]
                    plain_syn = Psyn[sel] @ Winv * gain[chn][None, None, :]
                    for j, e in enumerate(kk):
                        avgs.append({"w": w, "unit": k, "kind": kind, "d": e[2], "rej_first": e[5],
                                     "trace": plain[j, :, mainc].astype(np.float32),
                                     "trace_syn": plain_syn[j, :, mainc].astype(np.float32)})
    np.save(f"{OUT}/trace_pairs_traces_w{w:02d}.npy", np.array(avgs, dtype=object), allow_pickle=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    with ProcessPoolExecutor(len(WINDOWS)) as ex:
        dfs = list(ex.map(run_window, WINDOWS))
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(f"{OUT}/trace_pairs.csv.gz", index=False)
    print(df.groupby("kind").size().to_string())
