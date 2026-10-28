"""Stage A: what each admission rule keeps, on frozen templates, windows 0 and 1 of the hybrid binary.

``A_pt_open`` matched every tetrode on its own with a permissive floor (fit
amplitude 0.3 of the template and matched-filter signal-to-noise 2.5 or more).
Rules are applied to its fits afterwards; the joint runs (``A_joint``,
``A_joint_masked``) and per-tetrode circus-omp at June's floor (``A_pt``) are
the reference points. For each rule and window:

- ground-truth units: recall (all spikes and by class) and precision of the
  unit seeded from each, pooled over units in two size classes (trough 6 MAD
  or less, more than 6 MAD), and the share of ground-truth spikes a real unit
  took;
- real units: accepted spikes, share with an interval under 1.5 ms, share with
  a peeled trough under 3.5 MAD.

Window 1 holds the 2287 to 2696 s event on seven tetrodes; ``--event`` adds the
ground-truth recall inside it, per tetrode, for each run.

Writes ``admission_<w>.csv`` next to the arms and prints the table.

Run: .venv/bin/python investigations/carry-bakeoff/admission.py [--windows 0 1]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import CLASSES, SHORT, near  # noqa: E402

REF = ["A_joint", "A_joint_masked", "A_pt"]


def rules():
    """Matched-filter snr floor x relative floor on the fit amplitude x shape gate."""
    gates = {"": lambda z: True, " cos_pl>=0.5": lambda z: z["cos_pl"] >= 0.5,
             " cos_pl>=0.6": lambda z: z["cos_pl"] >= 0.6, " resid<=1.5xmed": lambda z: z["resid_rel"] <= 1.5}
    out = {"a>=0.8": lambda z: z["a"] >= 0.8}
    for k in (3.5, 4.0, 4.5, 5.0):
        for r in (0.0, 0.4, 0.5, 0.6):
            for gname, gate in gates.items():
                out[f"snr>={k} a>={r}{gname}"] = (lambda k, r, gate: lambda z: (z["snr"] >= k) & (z["a"] >= r)
                                                   & gate(z))(k, r, gate)
    return out


def refractory(z, keep):
    """Drop, within each unit, the lower-snr spike of any pair closer than 1.5 ms among kept spikes."""
    from ttm import capture as cp

    out = keep.copy()
    for u in np.unique(z["unit"][keep]):
        idx = np.flatnonzero(keep & (z["unit"] == u))
        k = cp.resolve_refractory(z["frame"][idx], np.nan_to_num(z["snr"][idx], nan=-1), SHORT)
        out[idx[~k]] = False
    return out


def evaluate(z, keep, gt, w, sizes):
    a, b = w * bk.WIN, (w + 1) * bk.WIN
    mw = (gt["frame"] >= a) & (gt["frame"] < b)
    gfr, gun, gcl = gt["frame"][mw], gt["unit"][mw], gt["cls"][mw]
    fr, un = z["frame"][keep], z["unit"][keep]
    res = {}
    for cls_name, members in sizes.items():
        tp = n_true = n_unit = 0
        by_cls = np.zeros((len(CLASSES), 2))
        for i in members:
            mine = np.sort(fr[un == 1000 + i])
            gi = gfr[gun == i]
            f = near(gi, mine)
            n_true += gi.size
            n_unit += mine.size
            tp += near(mine, np.sort(gi)).sum()
            for c in range(len(CLASSES)):
                m = gcl[gun == i] == c
                by_cls[c] += (f[m].sum(), m.sum())
            res.setdefault(f"recall_{cls_name}", 0)
            res[f"recall_{cls_name}"] += f.sum()
        res[f"recall_{cls_name}"] /= max(n_true, 1)
        res[f"precision_{cls_name}"] = tp / max(n_unit, 1)
        for c, name in enumerate(CLASSES):
            res[f"{name}_{cls_name}"] = by_cls[c, 0] / by_cls[c, 1] if by_cls[c, 1] >= 20 else np.nan
    ev = (gfr >= 2287 * bk.FS) & (gfr < 2696 * bk.FS) & np.isin(gun, np.flatnonzero((gt["group"] >= 4) & (gt["group"] <= 11)))
    if ev.sum() >= 100:
        hit = np.zeros(ev.sum(), bool)
        for i in np.unique(gun[ev]):
            m = gun[ev] == i
            hit[m] = near(gfr[ev][m], np.sort(fr[un == 1000 + i]))
        res["recall_event_TT5_12"] = hit.mean()
    real = un < 1000
    gall = np.sort(gfr)
    res["real_spikes"] = int(real.sum())
    res["real_took_gt"] = int(near(fr[real], gall).sum())
    sh, sm = [], []
    for u in np.unique(un[real]):
        f = np.sort(fr[un == u])
        if f.size >= 100:
            sh.append((np.diff(f) < SHORT).mean())
            sm.append((z["trough"][keep][un == u] < 3.5).mean())
    res["real_short_isi_median"] = float(np.median(sh)) if sh else np.nan
    res["real_small_median"] = float(np.median(sm)) if sm else np.nan
    return res


def load(name, w):
    z = dict(np.load(bk.BAKE / "arms" / name / f"spikes_w{w:02d}.npz"))
    # residual energy relative to the median over the tetrode's spikes with snr >= 4 (a correct fit's level)
    b0 = bk.Bank.load(bk.BAKE / "bank_start.npz")
    grp = dict(zip(b0.ids.tolist(), b0.group.tolist()))
    g = np.array([grp[int(u)] for u in z["unit"]])
    rel = np.full(z["resid"].size, np.nan)
    for gg in np.unique(g):
        m = g == gg
        ref = np.nanmedian(z["resid"][m & (z["snr"] >= 4)]) if (m & (z["snr"] >= 4)).any() else np.nanmedian(z["resid"][m])
        rel[m] = z["resid"][m] / ref
    z["resid_rel"] = rel
    z["group"] = g
    return z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--open-arm", default="A_pt_open")
    a = ap.parse_args()
    gt = np.load(bk.BAKE / "gt.npz")
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    sizes = {"small": [p["gt"] for p in plan if p["trough_mad"] <= 6],
             "large": [p["gt"] for p in plan if p["trough_mad"] > 6]}
    for w in a.windows:
        rows = []
        for name in REF:
            z = load(name, w)
            keep = z["reason"] == 0
            rows.append({"run": name, "rule": "as run", **evaluate(z, keep, gt, w, sizes)})
        z = load(a.open_arm, w)
        base = z["reason"] != 1
        for rname, rule in rules().items():
            keep = base & rule(z)
            if rname == "a>=0.8":
                rows.append({"run": a.open_arm, "rule": rname, **evaluate(z, keep, gt, w, sizes)})
            rows.append({"run": a.open_arm, "rule": rname + " +refr", **evaluate(z, refractory(z, keep), gt, w, sizes)})
        df = pd.DataFrame(rows)
        df.to_csv(bk.BAKE / f"admission_{a.open_arm}_w{w:02d}.csv", index=False)
        cols = ["run", "rule", "recall_small", "precision_small", "recall_large", "precision_large",
                "tail_small", "tail_large", "laser_half_large", "real_spikes", "real_took_gt",
                "real_short_isi_median", "real_small_median"] + (["recall_event_TT5_12"] if w == 1 else [])
        print(f"\n== window {w}")
        print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
