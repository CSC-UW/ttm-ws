"""Score arms of the carry-forward bake-off, window by window.

Ground-truth units (ids 1000 + i): a ground-truth spike is found when the
arm accepted a spike of the unit seeded from it within 10 samples (0.33 ms),
taken by another unit when another unit of the tetrode has one there
instead. Recall is split by spike class (isolated, burst lead, burst tail,
laser first spikes, laser half-size). Precision is the share of the unit's
accepted spikes that have a ground-truth spike of that unit within 10
samples. Template tracking is the bank's template after the window against
the injected template at the window's middle: its scale (projection on the
truth over the truth's energy) and cosine.

Real units (the seed bank): accepted spikes, share with an interval under
1.5 ms to the previous spike, share with a peeled trough under 3.5 excised
MAD, median cosine of the peeled snippet, median residual energy, spikes
within 10 samples of a ground-truth spike of the tetrode (taken from a
ground-truth unit), and the template's energy against window 0's.

Writes ``<arm>/score_gt.csv`` and ``<arm>/score_real.csv``.

Run: .venv/bin/python investigations/carry-bakeoff/score.py ARM [ARM ...]
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402

TOL = 10
SHORT = int(1.5e-3 * bk.FS)
CLASSES = ["iso", "lead", "tail", "laser_first", "laser_half"]
EVENT = (2287.0, 2696.0)  # s; the multi-tetrode event of window 1 (docs/reports/2026-10-26_ttm003_bad_periods.md)


def near(ref: np.ndarray, other: np.ndarray, tol: int = TOL) -> np.ndarray:
    """Whether each of ``ref`` has an element of sorted ``other`` within ``tol``."""
    if other.size == 0:
        return np.zeros(ref.size, bool)
    j = np.searchsorted(other, ref)
    hi = other[np.clip(j, 0, other.size - 1)]
    lo = other[np.clip(j - 1, 0, other.size - 1)]
    d = np.minimum(np.abs(hi - ref), np.abs(lo - ref))
    return d <= tol


def xtet_excess(f: np.ndarray, others: np.ndarray, tol: int = 6) -> float:
    """Share of ``f`` with a spike of ``others`` (sorted) within ``tol`` samples, minus the mean share at
    offsets of 2, 3, -2 and -3 ms (chance at the same rates)."""
    if f.size == 0 or others.size == 0:
        return np.nan
    obs = near(f, others, tol).mean()
    chance = np.mean([near(f + int(d * 1e-3 * bk.FS), others, tol).mean() for d in (2, 3, -2, -3)])
    return float(obs - chance)


def inpulse_rate(f: np.ndarray, on: np.ndarray, off: np.ndarray, a: int, b: int) -> float:
    """Spikes per second from 10 ms after each pulse onset to its offset, over the pulses in [a, b)."""
    m = (on >= a) & (on < b)
    if not m.any():
        return np.nan
    return float(bk.in_intervals(f, on[m], off[m]).sum() / ((off[m] - on[m]).sum() / bk.FS))


def score_arm(arm_dir) -> tuple[pd.DataFrame, pd.DataFrame]:
    gt = np.load(bk.BAKE / "gt.npz")
    b0 = bk.Bank.load(bk.BAKE / "bank_start.npz")
    e0 = {u: float((t.astype(np.float64) ** 2).sum()) for u, t in zip(b0.ids, b0.T)}
    grp = dict(zip(b0.ids.tolist(), b0.group.tolist()))
    rows_gt, rows_real = [], []
    import json

    cfg = json.loads((arm_dir / "arm.json").read_text())
    Wb = None
    if cfg.get("whiten"):
        W = bk.whitening_matrix(cfg.get("hybrid", True))
        groups = np.repeat(np.arange(16), 4)
        Wb = {g: W[np.ix_(groups == g, groups == g)] for g in range(16)}
    on, off = bk.laser_pulses()
    on = on + int(0.010 * bk.FS)
    w = 0
    while (arm_dir / f"spikes_w{w:02d}.npz").exists():
        sp = np.load(arm_dir / f"spikes_w{w:02d}.npz")
        bank = bk.Bank.load(arm_dir / f"bank_w{w:02d}.npz")
        acc = sp["reason"] == 0
        fr, un = sp["frame"][acc], sp["unit"][acc]
        a, b = w * bk.WIN, (w + 1) * bk.WIN
        mw = (gt["frame"] >= a) & (gt["frame"] < b)
        gfr, gun, gcl = gt["frame"][mw], gt["unit"][mw], gt["cls"][mw]
        tet_frames = {}
        for g in np.unique(b0.group):
            m = np.isin(un, b0.ids[b0.group == g])
            tet_frames[g] = (fr[m], un[m])
        for i in range(len(gt["group"])):
            u, g = 1000 + i, int(gt["group"][i])
            mine = np.sort(fr[un == u])
            gi = gfr[gun == i]
            found = near(gi, mine)
            tf, tu = tet_frames[g]
            other = np.sort(tf[tu != u])
            taken = ~found & near(gi, other)
            prec_hit = near(mine, np.sort(gi))
            row = {"window": w, "gt": i, "tetrode": f"TT{g + 1}", "n_true": gi.size, "n_unit": mine.size,
                   "recall": found.mean() if gi.size else np.nan,
                   "taken_by_other": taken.mean() if gi.size else np.nan,
                   "precision": prec_hit.mean() if mine.size else np.nan}
            for c, name in enumerate(CLASSES):
                mc = gcl[gun == i] == c
                row[f"recall_{name}"] = found[mc].mean() if mc.sum() >= 20 else np.nan
            ev = (gi >= EVENT[0] * bk.FS) & (gi < EVENT[1] * bk.FS)
            row["recall_event"] = found[ev].mean() if ev.sum() >= 20 else np.nan
            k = np.flatnonzero(bank.ids == u)[0]
            est = bank.T[k].astype(np.float64)
            true = gt["T_win"][i, w] * gt["size_win"][i, w]
            if Wb is not None:
                true = true @ Wb[g]
            row["tmpl_scale"] = float((est * true).sum() / (true * true).sum())
            row["tmpl_cos"] = float((est * true).sum() / np.sqrt((est * est).sum() * (true * true).sum()))
            rows_gt.append(row)
        all_fr = {g: np.sort(fr[np.isin(un, b0.ids[b0.group != g])]) for g in np.unique(b0.group)}
        for u in b0.ids[b0.ids < 1000]:
            m = un == u
            f = np.sort(fr[m])
            sel = acc & (sp["unit"] == u)
            g = grp[int(u)]
            gi = np.sort(gfr[np.isin(gun, np.flatnonzero(gt["group"] == g))])
            k = np.flatnonzero(bank.ids == u)[0]
            rows_real.append({
                "window": w, "unit": int(u), "tetrode": f"TT{g + 1}", "n": f.size,
                "short_isi": float((np.diff(f) < SHORT).mean()) if f.size > 1 else np.nan,
                "small": float((sp["trough"][sel] < 3.5).mean()) if f.size else np.nan,
                "cos_pl": float(np.nanmedian(sp["cos_pl"][sel])) if f.size else np.nan,
                "resid": float(np.nanmedian(sp["resid"][sel])) if f.size else np.nan,
                "snr": float(np.nanmedian(sp["snr"][sel])) if f.size else np.nan,
                "took_gt": int(near(f, gi).sum()) if f.size else 0,
                "inpulse_hz": inpulse_rate(f, on, off, a, b),
                "xtet_excess": xtet_excess(f, all_fr[g]),
                "energy_rel": float((bank.T[k].astype(np.float64) ** 2).sum() / e0[int(u)])})
        w += 1
    return pd.DataFrame(rows_gt), pd.DataFrame(rows_real)


def main():
    for name in sys.argv[1:]:
        d = bk.BAKE / "arms" / name
        g, r = score_arm(d)
        g.to_csv(d / "score_gt.csv", index=False)
        r.to_csv(d / "score_real.csv", index=False)
        last = g.window.max()
        print(f"== {name}: {last + 1} windows")
        print(g.groupby("window")[["recall", "precision", "taken_by_other", "tmpl_scale"]].median().round(3)
              .to_string())
        print(r.groupby("window")[["n", "short_isi", "small", "cos_pl", "energy_rel", "took_gt"]]
              .agg({"n": "sum", "short_isi": "median", "small": "median", "cos_pl": "median",
                    "energy_rel": "median", "took_gt": "sum"}).round(3).to_string())


if __name__ == "__main__":
    main()
