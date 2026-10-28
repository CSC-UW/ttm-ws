"""Figures of the carry-forward bake-off (run ``score.py`` and ``compare_carries.py`` on the arms first).

- ``ladder.png``: ground-truth recall and precision per window, small and large units, for June's
  carry (joint, masked), joint matching on whitened traces, and the per-tetrode whitened carry with
  frozen templates and with blend 0.2 re-estimation.
- ``template_scale.png``: each ground-truth unit's template scale against the injected template per
  window, June's carry against the per-tetrode whitened carries.
- ``laser.png``: recall of the driven units' half-size in-pulse spikes and of all their spikes in the
  four train windows, per arm.

Written to ``<_bake.BAKE>/figures/``.

Run: .venv/bin/python investigations/carry-bakeoff/figures.py
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

HERE = pathlib.Path(__file__).parent
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
INK, MUTED, GRID = "#1a1a19", "#6b6a63", "#e4e3dd"
ARMS = [("joint_june_masked", "June: joint, masked"), ("W_joint_masked_june", "joint, whitened, masked"),
        ("W_s45_c50_frozen", "per tetrode, whitened, frozen"), ("W_s45_c50_lnb_b2", "per tetrode, whitened, blend 0.2")]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def ladder(df, out):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True, sharey=True)
    for r, cls in enumerate(("small", "large")):
        for c, metric in enumerate(("recall", "precision")):
            ax = axes[r, c]
            style(ax)
            for k, (arm, label) in enumerate(ARMS):
                d = df[df.arm == arm]
                ax.plot(d.window, d[f"{metric}_{cls}"], color=SERIES[k], lw=2, marker="o", ms=4, label=label)
            ax.set_ylim(0, 1)
            ax.set_title(f"{metric}, {'trough 6 MAD or less' if cls == 'small' else 'trough over 6 MAD'}",
                         fontsize=9, color=INK, loc="left")
            if r == 1:
                ax.set_xlabel("1800 s window", fontsize=8, color=MUTED)
            ax.axvspan(0.5, 1.5, color=GRID, alpha=0.6, lw=0)
    axes[0, 0].text(1, 0.03, "event", ha="center", fontsize=7, color=MUTED)
    axes[0, 0].legend(fontsize=7, frameon=False, loc="lower left", bbox_to_anchor=(0, 0.05))
    fig.suptitle("Ground-truth units over 12 carry-forward windows (pooled over units)", fontsize=10, color=INK,
                 x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out / "ladder.png", dpi=110)


def template_scale(out):
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    arms = [("joint_june_masked", "June: joint, masked"), ("W_s45_c50_lnb", "per tetrode, whitened, full replace"),
            ("W_s45_c50_lnb_b2", "per tetrode, whitened, blend 0.2")]
    fig, axes = plt.subplots(1, len(arms), figsize=(12, 3.8), sharey=True)
    for ax, (arm, label) in zip(axes, arms):
        style(ax)
        g = pd.read_csv(bk.BAKE / "arms" / arm / "score_gt.csv")
        for p in plan:
            d = g[g["gt"] == p["gt"]]
            small = p["trough_mad"] <= 6
            ax.plot(d.window, d.tmpl_scale, color=SERIES[0] if small else SERIES[1], lw=1, alpha=0.7)
        ax.axhline(1, color=INK, lw=0.8)
        ax.set_yscale("log")
        ax.set_ylim(0.4, 4)
        ax.set_yticks([0.5, 1, 2, 4])
        ax.set_yticklabels(["0.5", "1", "2", "4"])
        ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_title(label, fontsize=9, color=INK, loc="left")
        ax.set_xlabel("1800 s window", fontsize=8, color=MUTED)
    axes[0].set_ylabel("template scale / injected", fontsize=8, color=MUTED)
    axes[0].plot([], [], color=SERIES[0], label="trough 6 MAD or less")
    axes[0].plot([], [], color=SERIES[1], label="trough over 6 MAD")
    axes[0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out / "template_scale.png", dpi=110)


def laser(out):
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    driven = [p for p in plan if p["kind"] == "driven"]
    arms = [("joint_june_masked", "June: joint, masked"), ("W_joint_masked_june", "joint, whitened"),
            ("D_s4_frozen", "per tetrode, plain, pulse rule"), ("W_s45_c50_lnb_b2", "per tetrode, whitened, pulse rule")]
    fig, ax = plt.subplots(figsize=(9, 3.6))
    style(ax)
    width = 0.8 / len(arms)
    for k, (arm, label) in enumerate(arms):
        g = pd.read_csv(bk.BAKE / "arms" / arm / "score_gt.csv")
        vals = [g[(g["gt"] == p["gt"]) & g.window.isin([0, 3, 8, 11])].recall_laser_half.mean() for p in driven]
        x = np.arange(len(driven)) + (k - (len(arms) - 1) / 2) * width
        ax.bar(x, vals, width * 0.92, color=SERIES[k], label=label)
    ax.set_xticks(np.arange(len(driven)))
    ax.set_xticklabels([f"GT{p['gt']} {p['tetrode']}, {p['trough_mad']:g} MAD" for p in driven], fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("half-size in-pulse spikes found", fontsize=8, color=MUTED)
    ax.legend(fontsize=7, frameon=False, loc="upper right")
    ax.set_title("Driven ground-truth units, the four train windows", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(out / "laser.png", dpi=110)


def main():
    out = bk.BAKE / "figures"
    out.mkdir(exist_ok=True)
    df = pd.read_csv(HERE / "results" / "carries_summary_final.csv")
    ladder(df, out)
    template_scale(out)
    laser(out)
    print(f"figures in {out}")


if __name__ == "__main__":
    main()
