"""Stricter admission applied after the fact to the per-tetrode carry, against the joint carry.

On the aligned-reference hybrid with the re-sorted seed bank (``BAKE_REFERENCE=aligned_seed``),
windows 0 to 11: the accepted spikes of the ground-truth units in ``W_s45_c50_lnb_b2`` (per tetrode)
are filtered again by floors on snr, the relative fit amplitude ``a`` and ``cos_pl``, and by an
overlap rule (drop a spike within 10 samples of a spike of another unit of the tetrode with larger
snr, when its own snr is under a cap), and by a floor relative to the unit's median snr over its
accepted spikes in the window (snr >= max(4.5, f x median); in a carry it would be the previous
window's median). Recall and precision are pooled per window over small units
(trough 6 MAD or less) and large ones and averaged over the windows, as in ``compare_reference.py``;
the units that stop after window 5 are left out from window 6. ``W_joint_masked_s45_c50_lnb_b2``
(joint) is scored the same way, unfiltered. The filter does not feed back into the carry: templates,
peeling and the next window's bank are those of the run, so this measures only what each rule would
remove, not how the carry would continue with it.

Writes ``results/rule_sweep.csv``.

Run: BAKE_REFERENCE=aligned_seed .venv/bin/python investigations/carry-bakeoff/rule_sweep.py
"""
from __future__ import annotations

import itertools
import json
import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

HERE = pathlib.Path(__file__).parent
PT, JT = "W_s45_c50_lnb_b2", "W_joint_masked_s45_c50_lnb_b2"
TOL = 10


def near(ref, other, tol=TOL):
    if other.size == 0 or ref.size == 0:
        return np.zeros(ref.size, bool)
    i = np.searchsorted(other, ref - tol, side="left")
    j = np.searchsorted(other, ref + tol, side="right")
    return j > i


def overlap_larger(fr, snr, unit, group_of):
    """Per spike: another unit of the same tetrode has a spike with larger snr within TOL samples."""
    out = np.zeros(fr.size, bool)
    g = np.array([group_of.get(int(u), -1) for u in unit])
    for gg in np.unique(g):
        m = np.flatnonzero(g == gg)
        o = m[np.argsort(fr[m], kind="stable")]
        f, s, u = fr[o], snr[o], unit[o]
        for d in range(1, 40):
            if d >= f.size:
                break
            close = (f[d:] - f[:-d]) <= TOL
            if not close.any():
                break
            diff = u[d:] != u[:-d]
            k = np.flatnonzero(close & diff)
            out[o[k]] |= s[k + d] > s[k]
            out[o[k + d]] |= s[k] > s[k + d]
    return out


def load(arm, gt, leaves):
    """Per window: dataframe of accepted ground-truth-unit spikes with features and hit flags, and n_true."""
    out = []
    for w in range(12):
        sp = np.load(bk.BAKE / "arms" / arm / f"spikes_w{w:02d}.npz")
        bank = bk.Bank.load(bk.BAKE / "arms" / arm / f"bank_w{w:02d}.npz")
        group_of = dict(zip(bank.ids.tolist(), bank.group.tolist()))
        acc = sp["reason"] == 0
        fr, un, snr = sp["frame"][acc], sp["unit"][acc], sp["snr"][acc]
        ov = overlap_larger(fr, snr, un, group_of)
        a, b = w * bk.WIN, (w + 1) * bk.WIN
        mw = (gt["frame"] >= a) & (gt["frame"] < b)
        gfr, gun = gt["frame"][mw], gt["unit"][mw]
        rows, ntrue = [], {}
        for i in range(len(gt["group"])):
            if i in leaves and w >= 6:
                continue
            own = np.sort(gfr[gun == i])
            ntrue[i] = own
            m = un == 1000 + i
            d = pd.DataFrame({"frame": fr[m], "snr": snr[m], "a": sp["a"][acc][m], "cos_pl": sp["cos_pl"][acc][m],
                              "ov": ov[m]})
            d["gt_unit"] = i
            d["hit"] = near(d.frame.to_numpy(), own)
            d["med_snr"] = float(np.median(d.snr)) if len(d) else np.nan
            rows.append(d)
        out.append((pd.concat(rows, ignore_index=True), ntrue))
    return out


def score(windows, keep_fn, small, large):
    res = {"small": [], "large": []}
    for d, ntrue in windows:
        k = keep_fn(d)
        for name, members in (("small", small), ("large", large)):
            found = n_true = hits = n_unit = 0
            for i in members:
                if i not in ntrue:
                    continue
                own = ntrue[i]
                sel = (d.gt_unit == i).to_numpy() & k
                f = np.sort(d.frame.to_numpy()[sel])
                found += int(near(own, f).sum())
                n_true += own.size
                hits += int(d.hit.to_numpy()[sel].sum())
                n_unit += int(sel.sum())
            res[name].append((found / max(n_true, 1), hits / max(n_unit, 1)))
    out = {}
    for name, v in res.items():
        r, p = np.mean([x for x, _ in v]), np.mean([y for _, y in v])
        out.update({f"{name}_recall": r, f"{name}_precision": p, f"{name}_f1": 2 * r * p / (r + p)})
    return out


def main():
    assert bk.REFERENCE == "aligned_seed", "run with BAKE_REFERENCE=aligned_seed"
    gt = np.load(bk.BAKE_ROOT / "v2_aligned" / "gt.npz")
    plan = json.loads((bk.BAKE_ROOT / "v2" / "gt_plan.json").read_text())
    leaves = {p["gt"] for p in plan if p["kind"] == "leaves"}
    small = [p["gt"] for p in plan if p["trough_mad"] <= 6]
    large = [p["gt"] for p in plan if p["trough_mad"] > 6]
    pt, jt = load(PT, gt, leaves), load(JT, gt, leaves)
    rows = [{"arm": "joint", "rule": "as run", **score(jt, lambda d: np.ones(len(d), bool), small, large)},
            {"arm": "per tetrode", "rule": "as run", **score(pt, lambda d: np.ones(len(d), bool), small, large)}]
    for s, a, c, ovcap in itertools.product((4.5, 5.0, 5.5, 6.0), (0.5, 0.6, 0.7), (0.5, 0.55, 0.6),
                                            (None, 6.0, 8.0, np.inf)):
        if (s, a, c, ovcap) == (4.5, 0.5, 0.5, None):
            continue

        def keep(d, s=s, a=a, c=c, ovcap=ovcap):
            k = (d.snr >= s) & (d.a >= a) & (d.cos_pl >= c)
            if ovcap is not None:
                k &= ~(d.ov & (d.snr < ovcap))
            return k.to_numpy()

        rows.append({"arm": "per tetrode", "rule": f"snr>={s} a>={a} cos_pl>={c}"
                     + ("" if ovcap is None else f" overlap<{ovcap}"), **score(pt, keep, small, large)})
    for f, ovcap in itertools.product((0.3, 0.35, 0.4, 0.45, 0.5, 0.55), (None, 8.0)):
        def keep(d, f=f, ovcap=ovcap):
            k = d.snr >= np.maximum(4.5, f * d.med_snr)
            if ovcap is not None:
                k &= ~(d.ov & (d.snr < ovcap))
            return k.to_numpy()

        rows.append({"arm": "per tetrode", "rule": f"snr>=max(4.5, {f} x unit median)"
                     + ("" if ovcap is None else f" overlap<{ovcap}"), **score(pt, keep, small, large)})
    df = pd.DataFrame(rows).round(4)
    df.to_csv(HERE / "results" / "rule_sweep.csv", index=False)
    j = df.iloc[0]
    pd.set_option("display.width", 250)
    print(df.iloc[:2].to_string(index=False))
    dom = df[(df.arm == "per tetrode") & (df.small_recall >= j.small_recall) & (df.small_precision >= j.small_precision)
             & (df.large_recall >= j.large_recall) & (df.large_precision >= j.large_precision)]
    print(f"\nper-tetrode rules at or above joint on all four: {len(dom)}")
    print(dom.sort_values("small_f1", ascending=False).head(15).to_string(index=False))
    print("\nbest per-tetrode rules by small F1 + large F1")
    print("\nrelative floors")
    print(df[df.rule.str.contains("median")].to_string(index=False))
    print(df.assign(t=df.small_f1 + df.large_f1).sort_values("t", ascending=False).head(15).drop(columns="t")
          .to_string(index=False))


if __name__ == "__main__":
    main()
