"""Ceiling of sub-sample template fitting on the carry's accepted and rejected fits.

For one arm's windows on the aligned hybrid, per tetrode: rebuild the residual of all
fits from the arm's spikes file and the bank it matched with (bankin_wNN), take each
fit's peeled snippet (residual plus its own fit) and evaluate the template at
fractional delays d in [-0.5, 0.5] (FFT shift, step 1/16). The peeled snippet does
not depend on the fit's own shift, so snr and cos_pl at d are its projection on the
shifted unit-norm template over sigma and over the snippet's norm. Residual energy
compares the stored fit (a * T) with the least-squares fit at the integer and at the
best delay. Ground-truth fits are hits when a spike of the unit lies within 10
samples; their injected quarter-sample delay checks the estimate.

Template blur: per unit, the mean of up to 500 accepted peeled snippets as stored
(integer aligned) against the mean after shifting each by its own best delay and,
for ground-truth units, by its true injected delay.

Reads only; writes per-fit and per-unit tables to the scratchpad.

Run: BAKE_REFERENCE=aligned_seed .venv/bin/python investigations/template-superres/ceiling.py --windows 3 8
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "carry-bakeoff"))
import _bake as bk  # noqa: E402
from ttm import capture as cp  # noqa: E402

ARM = bk.BAKE / "arms" / "W_s45_c50_lnb_b2"
GT = bk.BAKE_ROOT / "v2_aligned" / "gt.npz"
PLAN = bk.BAKE_ROOT / "v2_aligned" / "gt_plan.json"
OUT = pathlib.Path("/tmp/claude-295375638/-Users-gfindlay-ad-wisc-edu-projects-ttm-ws/"
                   "e980bd61-79fb-4043-a37f-61959b628a69/scratchpad/superres")
D = np.linspace(-0.5, 0.5, 17)
NPAD = 256


def shift_bank(T: np.ndarray, d: np.ndarray) -> np.ndarray:
    """(n_units, n_d, N_SAMP, 4): each template delayed by each d (band-limited, zero padded)."""
    pad = np.zeros((T.shape[0], NPAD, T.shape[2]))
    pad[:, 64:64 + T.shape[1]] = T
    f = np.fft.rfftfreq(NPAD)
    ph = np.exp(-2j * np.pi * f[None, :] * d[:, None])
    F = np.fft.rfft(pad, axis=1)
    out = np.fft.irfft(F[:, None] * ph[None, :, :, None], n=NPAD, axis=2)
    return out[:, :, 64:64 + T.shape[1]]


def shift_snips(S: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Each snippet S[i] (N_SAMP, 4) delayed by d[i]."""
    pad = np.zeros((S.shape[0], NPAD, S.shape[2]))
    pad[:, 64:64 + S.shape[1]] = S
    f = np.fft.rfftfreq(NPAD)
    out = np.fft.irfft(np.fft.rfft(pad, axis=1) * np.exp(-2j * np.pi * f[None, :] * d[:, None])[:, :, None],
                       n=NPAD, axis=1)
    return out[:, 64:64 + S.shape[1]]


def near(ref: np.ndarray, sorted_other: np.ndarray, tol: int = 10) -> tuple[np.ndarray, np.ndarray]:
    """Whether each ref has an element of sorted_other within tol, and the index of the nearest."""
    if sorted_other.size == 0:
        return np.zeros(ref.size, bool), np.zeros(ref.size, np.int64)
    j = np.clip(np.searchsorted(sorted_other, ref), 1, sorted_other.size - 1)
    jl = j - 1
    pick = np.where(np.abs(sorted_other[jl] - ref) <= np.abs(sorted_other[j] - ref), jl, j)
    return np.abs(sorted_other[pick] - ref) <= tol, pick


def tetrode(w: int, g: int, W: np.ndarray, max_per_unit: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    t0 = time.time()
    a, b = w * bk.WIN, (w + 1) * bk.WIN
    rec = bk.si.load(bk.HYBRID)
    groups = np.asarray(rec.get_property("group")).astype(int)
    ch = np.flatnonzero(groups == g)
    Wg = W[np.ix_(ch, ch)].astype(np.float64)
    xp = cp.load_tetrode_traces(rec, g, a, b)
    noise_plain, _ = cp.noise_mad_excised(xp)
    x = (xp.astype(np.float64) @ Wg).astype(np.float32)
    del xp
    noise, _ = cp.noise_mad_excised(x)
    t_load = time.time() - t0

    bank = np.load(ARM / f"bankin_w{w:02d}.npz")
    gi = np.flatnonzero(bank["group"] == g)
    T = bank["T"][gi].astype(np.float32)
    ids = bank["ids"][gi]
    sigma = bank["sigma"][gi]
    Wi = np.linalg.inv(Wg)
    trough_mad = np.array([-(t.astype(np.float64) @ Wi / noise_plain).min() for t in T])

    sp = np.load(ARM / f"spikes_w{w:02d}.npz")
    m = np.isin(sp["unit"], ids)
    fr = sp["frame"][m] - a
    lab = np.searchsorted(ids, sp["unit"][m])
    assert np.array_equal(ids[lab], sp["unit"][m])
    amp = sp["a"][m]
    reason = sp["reason"][m]
    snr_st, cos_st = sp["snr"][m], sp["cos_pl"][m]

    t1 = time.time()
    model = np.zeros_like(x)
    bk.add_spikes(model, fr, lab, amp.astype(np.float32), T, bk.NB)
    resid = x - model
    del model
    t_resid = time.time() - t1

    gt = np.load(GT)
    mw = (gt["frame"] >= a) & (gt["frame"] < b)
    gfr, gun, gup = gt["frame"][mw] - a, gt["unit"][mw], gt["up"][mw]

    rng = np.random.default_rng(seed + 100 * w + g)
    ok = (fr >= bk.NB + 1) & (fr + bk.NA + 1 <= x.shape[0])
    TD = shift_bank(T.astype(np.float64), D)  # (k, nd, ns, 4)
    tdn = np.sqrt((TD ** 2).sum(axis=(2, 3)))
    rows, urows = [], []
    t_score = 0.0
    for k in range(len(ids)):
        uid = int(ids[k])
        idx = np.flatnonzero(ok & (lab == k) & (reason <= 2))
        if idx.size == 0:
            continue
        if idx.size > max_per_unit:
            idx = np.sort(rng.choice(idx, max_per_unit, replace=False))
        t2 = time.time()
        r = cp.snippets(resid, fr[idx]).astype(np.float64)
        pl = r + amp[idx, None, None] * T[k][None].astype(np.float64)
        proj = np.einsum("ntc,dtc->nd", pl, TD[k]) / tdn[k][None]
        j0 = int(np.argmin(np.abs(D)))
        jb = proj.argmax(axis=1)
        # parabolic refinement of the best delay on the grid
        jc = np.clip(jb, 1, D.size - 2)
        y0, y1, y2 = proj[np.arange(idx.size), jc - 1], proj[np.arange(idx.size), jc], proj[np.arange(idx.size), jc + 1]
        den = y0 - 2 * y1 + y2
        dref = D[jc] + np.where(den < 0, 0.5 * (y0 - y2) / np.where(den < 0, den, -1), 0) * (D[1] - D[0])
        dref = np.clip(dref, -0.5, 0.5)
        pn = np.sqrt((pl ** 2).sum(axis=(1, 2)))
        p0, pb = proj[:, j0], proj[np.arange(idx.size), jb]
        Tb = TD[k][jb]
        a_b = pb / tdn[k][jb]
        a_0 = p0 / tdn[k][j0]
        inv2 = 1.0 / noise[None, None, :] ** 2
        e_st = ((r ** 2) * inv2).mean(axis=(1, 2))
        T0 = T[k][None].astype(np.float64)
        e_0 = (((pl - a_0[:, None, None] * T0) ** 2) * inv2).mean(axis=(1, 2))
        e_b = (((pl - a_b[:, None, None] * Tb) ** 2) * inv2).mean(axis=(1, 2))
        pk_st = (np.abs(r) / noise[None, None, :]).max(axis=(1, 2))
        pk_b = (np.abs(pl - a_b[:, None, None] * Tb) / noise[None, None, :]).max(axis=(1, 2))
        t_score += time.time() - t2
        df = pd.DataFrame({
            "window": w, "group": g, "unit": uid, "frame": fr[idx] + a, "reason": reason[idx], "a": amp[idx],
            "snr_stored": snr_st[idx], "cos_stored": cos_st[idx],
            "snr0": p0 / sigma[k], "snr_best": pb / sigma[k], "cos0": p0 / pn, "cos_best": pb / pn,
            "d_best": dref, "e_stored": e_st, "e_ls0": e_0, "e_best": e_b, "peak_stored": pk_st, "peak_best": pk_b,
            "trough_mad": trough_mad[k],
        })
        if uid >= 1000:
            i = uid - 1000
            gsel = np.flatnonzero(gun == i)
            gf = gfr[gsel]
            hit, pick = near(fr[idx], gf)
            df["hit"] = hit
            df["gt_off"] = (np.where(hit, gf[pick] + gup[gsel][pick] / 4.0 - fr[idx], np.nan) if gf.size
                            else np.nan)
        rows.append(df)

        # template blur: accepted spikes, as the rebuild would see them (peeled, integer aligned)
        acc = np.flatnonzero(ok & (lab == k) & (reason == 0))
        if acc.size >= 100:
            if acc.size > 500:
                acc = np.sort(rng.choice(acc, 500, replace=False))
            r = cp.snippets(resid, fr[acc]).astype(np.float64)
            pl = r + amp[acc, None, None] * T[k][None].astype(np.float64)
            proj = np.einsum("ntc,dtc->nd", pl, TD[k])
            db = D[proj.argmax(axis=1)]
            M0 = pl.mean(axis=0)
            M1 = shift_snips(pl, -db).mean(axis=0)
            u = {"window": w, "group": g, "unit": uid, "n": acc.size, "trough_mad": trough_mad[k],
                 "E0": float((M0 ** 2).sum()), "E1": float((M1 ** 2).sum()),
                 "min0": float(M0.min()), "min1": float(M1.min()), "sd_d": float(db.std())}
            # same procedure on noise: residual snippets at the spikes' times shifted by 45 samples
            rn = cp.snippets(resid, np.clip(fr[acc] + 45, bk.NB, x.shape[0] - bk.NA - 1)).astype(np.float64)
            pn_ = rn + amp[acc, None, None] * T[k][None].astype(np.float64)
            dn = D[np.einsum("ntc,dtc->nd", pn_, TD[k]).argmax(axis=1)]
            u["E1_noise_ctrl"] = float((shift_snips(pn_, -dn).mean(axis=0) ** 2).sum())
            u["E0_noise_ctrl"] = float((pn_.mean(axis=0) ** 2).sum())
            if uid >= 1000 and np.any(gun == uid - 1000):
                i = uid - 1000
                gsel = np.flatnonzero(gun == i)
                hit, pick = near(fr[acc], gfr[gsel])
                true = gt["T_win"][i, w] * gt["size_win"][i, w] @ Wg
                off = gfr[gsel][pick] + gup[gsel][pick] / 4.0 - fr[acc]
                h = np.flatnonzero(hit)
                c = np.median(off[h]) if h.size else 0.0
                # oracle: shift each hit by its true offset (relative to the unit's median offset)
                M2 = shift_snips(pl[h], -(off[h] - c)).mean(axis=0) if h.size else np.full_like(M0, np.nan)
                # truth delayed by the median offset so its sample NB lines up with the template's
                tr = shift_snips(true[None], np.array([-c + 0.0]))[0] if np.isfinite(c) else true

                def cos(p, q):
                    return float((p * q).sum() / np.sqrt((p * p).sum() * (q * q).sum()))
                u.update({"E2_oracle": float((M2 ** 2).sum()), "E_true": float((true ** 2).sum()),
                          "cos0_true": cos(M0, tr), "cos1_true": cos(M1, tr), "cos2_true": cos(M2, tr),
                          "hit_share": float(hit.mean())})
            urows.append(u)
    timing = {"window": w, "group": g, "load_s": t_load, "resid_s": t_resid, "score_s": t_score,
              "n_fits_tetrode": int(fr.size), "total_s": time.time() - t0}
    return (pd.concat(rows) if rows else pd.DataFrame()), pd.DataFrame(urows), timing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", type=int, nargs="+", default=[3, 8])
    ap.add_argument("--groups", type=int, nargs="*", default=None)
    ap.add_argument("--max-per-unit", type=int, default=20000)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    from joblib import Parallel, delayed

    OUT.mkdir(parents=True, exist_ok=True)
    W = np.load(bk.BAKE / "whiten_W.npz")["W"]
    groups = args.groups if args.groups else list(range(16))
    jobs = [(w, g) for w in args.windows for g in groups]
    res = Parallel(n_jobs=args.workers)(delayed(tetrode)(w, g, W, args.max_per_unit, 0) for w, g in jobs)
    fits = pd.concat([r[0] for r in res if len(r[0])])
    units = pd.concat([r[1] for r in res if len(r[1])])
    tag = "_".join(f"w{w:02d}" for w in args.windows)
    fits.to_parquet(OUT / f"fits_{tag}.parquet")
    units.to_csv(OUT / f"units_{tag}.csv", index=False)
    (OUT / f"timing_{tag}.json").write_text(json.dumps([r[2] for r in res], indent=1))
    print(fits.shape, units.shape)


if __name__ == "__main__":
    main()
