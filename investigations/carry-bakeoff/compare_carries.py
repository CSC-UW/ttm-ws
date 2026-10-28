"""Compare carry-forward arms window by window (run ``score.py`` on each first).

Per arm and window, pooled over ground-truth units in two size classes
(trough 6 MAD or less, more): recall, precision, recall of burst-tail and
laser half-size spikes, the share of ground-truth spikes another unit took,
and the median template scale against the truth (1 tracks it; below 1 the
template shrank more than the injected unit). For the two units that stop
firing after window 5, the spikes their templates collect afterwards. For
real units: accepted spikes, median share with an interval under 1.5 ms,
median share under 3.5 MAD, median template energy against window 0's, and
ground-truth spikes they took.

Writes ``results/carries_summary.csv`` and, under ``<_bake.BAKE>/figures/``,
``carries_gt.png`` (per ground-truth unit: template scale and recall over
windows, one line per arm) and ``carries_real.png`` (per real unit:
template energy over windows).

Run: .venv/bin/python investigations/carry-bakeoff/compare_carries.py ARM [ARM ...]
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

HERE = pathlib.Path(__file__).parent


def pooled(g: pd.DataFrame, members) -> pd.Series:
    d = g[g["gt"].isin(members)]
    found = (d.recall * d.n_true).sum()
    hits = (d.precision.fillna(0) * d.n_unit).sum()
    return pd.Series({"recall": found / max(d.n_true.sum(), 1), "precision": hits / max(d.n_unit.sum(), 1),
                      "taken": (d.taken_by_other * d.n_true).sum() / max(d.n_true.sum(), 1),
                      "tail": d.recall_tail.median(), "laser_half": d.recall_laser_half.median(),
                      "event": d.recall_event.median(),
                      "scale": d.tmpl_scale.median(), "scale_min": d.tmpl_scale.min()})


def main():
    arms = sys.argv[1:]
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    small = [p["gt"] for p in plan if p["trough_mad"] <= 6]
    large = [p["gt"] for p in plan if p["trough_mad"] > 6]
    leaves = [p["gt"] for p in plan if p["kind"] == "leaves"]
    rows = []
    gts, reals = {}, {}
    for a in arms:
        d = bk.BAKE / "arms" / a
        g = pd.read_csv(d / "score_gt.csv")
        r = pd.read_csv(d / "score_real.csv")
        gts[a], reals[a] = g, r
        for w, gw in g.groupby("window"):
            row = {"arm": a, "window": w}
            for cname, mem in (("small", small), ("large", large)):
                row.update({f"{k}_{cname}": v for k, v in pooled(gw[~gw["gt"].isin(leaves) | (w < 6)], mem).items()})
            row["left_collected"] = int(gw[gw["gt"].isin(leaves)].n_unit.sum()) if w >= 6 else 0
            rw = r[r.window == w]
            row.update({"real_n": int(rw.n.sum()), "real_short_isi": rw.short_isi.median(),
                        "real_small": rw.small.median(), "real_energy": rw.energy_rel.median(),
                        "real_energy_min": rw.energy_rel.min(), "real_took_gt": int(rw.took_gt.sum())})
            rows.append(row)
    df = pd.DataFrame(rows)
    (HERE / "results").mkdir(exist_ok=True)
    df.to_csv(HERE / "results" / "carries_summary.csv", index=False)
    pd.set_option("display.width", 250)
    for a in arms:
        print(f"\n== {a}")
        print(df[df.arm == a].drop(columns="arm").to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    figdir = bk.BAKE / "figures"
    figdir.mkdir(exist_ok=True)
    n = len(plan)
    nr = int(np.ceil(n / 10))
    fig, axes = plt.subplots(2 * nr, 10, figsize=(30, 5.5 * nr), sharex=True)
    for i, p in enumerate(plan):
        ax_s, ax_r = axes[2 * (i // 10), i % 10], axes[2 * (i // 10) + 1, i % 10]
        for a in arms:
            g = gts[a][gts[a]["gt"] == i]
            ax_s.plot(g.window, g.tmpl_scale, marker=".", label=a)
            ax_r.plot(g.window, g.recall, marker=".", label=a)
            ax_r.plot(g.window, g.precision, ls=":", color=ax_r.lines[-1].get_color())
        ax_s.axhline(1, color="k", lw=0.5)
        ax_s.set_title(f"GT{i} {p['tetrode']} {p['trough_mad']} MAD {p['course']} {p['kind']}", fontsize=8)
        ax_s.set_ylim(0, 1.6)
        ax_r.set_ylim(0, 1.02)
    for j in range(n, 10 * nr):
        axes[2 * (j // 10), j % 10].axis("off")
        axes[2 * (j // 10) + 1, j % 10].axis("off")
    axes[0, 0].set_ylabel("template scale / truth")
    axes[1, 0].set_ylabel("recall (solid), precision (dotted)")
    axes[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figdir / "carries_gt.png", dpi=80)
    units = sorted(reals[arms[0]].unit.unique())
    nc = 9
    fig, axes = plt.subplots(int(np.ceil(len(units) / nc)), nc, figsize=(27, 2.2 * np.ceil(len(units) / nc)),
                             sharex=True)
    for k, u in enumerate(units):
        ax = axes.flat[k]
        for a in arms:
            r = reals[a][reals[a].unit == u]
            ax.plot(r.window, r.energy_rel, marker=".", label=a)
        ax.axhline(1, color="k", lw=0.5)
        ax.set_title(f"unit {u} {reals[arms[0]][reals[arms[0]].unit == u].tetrode.iloc[0]}", fontsize=8)
    for k in range(len(units), axes.size):
        axes.flat[k].axis("off")
    axes.flat[0].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(figdir / "carries_real.png", dpi=80)
    print(f"figures in {figdir}")


if __name__ == "__main__":
    main()
