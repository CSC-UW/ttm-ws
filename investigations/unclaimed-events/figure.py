"""Waveforms of the events the fold and multi-unit steps would take, next to the templates, from
replay.py's per-window files.

    .venv/bin/python investigations/unclaimed-events/figure.py <replay out dir> <figure path>

Top row: the four units with the most folded events; each unit's template averaged over the windows
(black) and the mean of its folded events (red), the four wires side by side, µV. Lower rows: TT5 to
TT12; the templates of the tetrode's units (grey, averaged over the windows), 30 random multi-unit events
(light red) and the mean of the multi-unit events sampled (up to 50 per window, red).
"""

import pathlib
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAIN = 0.195
N_SAMP = 90


def flat(w):
    """(..., 90, 4) -> (..., 4 * 92), the wires side by side with a gap."""
    pad = np.full(w.shape[:-2] + (2, 4), np.nan)
    return np.concatenate([w, pad], axis=-2).swapaxes(-1, -2).reshape(*w.shape[:-2], -1)


def main():
    d, fig_path = pathlib.Path(sys.argv[1]), sys.argv[2]
    tsum, tcnt, fold, mua, mua_g, group_of = {}, {}, {}, [], [], {}
    for f in sorted(d.glob("w*.npz")):
        z = np.load(f)
        for u, t in zip(z["tmpl_ids"].tolist(), z["tmpl_T"]):
            tsum[u] = tsum.get(u, 0) + t
            tcnt[u] = tcnt.get(u, 0) + 1
        ok = z["ev_best_uid"] >= 0
        group_of.update(dict(zip(z["ev_best_uid"][ok].tolist(), z["ev_group"][ok].tolist())))
        folded = z["ev_folded"][z["prop_idx"]]
        for u, sn in zip(z["ev_best_uid"][z["prop_idx"]][folded].tolist(), z["prop_snip"][folded]):
            fold.setdefault(u, []).append(sn)
        mua.append(z["mua_snip"])
        mua_g.append(z["ev_group"][z["mua_idx"]])
    T = {u: tsum[u] / tcnt[u] for u in tsum}
    mua, mua_g = np.concatenate(mua), np.concatenate(mua_g)
    top = sorted(fold, key=lambda u: -len(fold[u]))[:4]
    fig, axes = plt.subplots(3, 4, figsize=(16, 9))
    x = np.arange(4 * (N_SAMP + 2))
    for ax, u in zip(axes[0], top):
        ax.plot(x, flat(T[u]) * GAIN, "k", lw=1)
        ax.plot(x, flat(np.mean(fold[u], axis=0)) * GAIN, "r", lw=1)
        ax.set_title(f"unit {u} (TT{group_of[u] + 1}): {len(fold[u])} folded", fontsize=9)
    rng = np.random.default_rng(0)
    for ax, g in zip(axes[1:].ravel(), range(4, 12)):
        for u, t in T.items():
            if group_of.get(u) == g:
                ax.plot(x, flat(t) * GAIN, color="0.5", lw=1)
        sn = mua[mua_g == g]
        for i in rng.permutation(len(sn))[:30]:
            ax.plot(x, flat(sn[i]) * GAIN, color="r", alpha=0.1, lw=0.5)
        if len(sn):
            ax.plot(x, flat(sn.mean(axis=0)) * GAIN, "r", lw=1.2)
        ax.set_title(f"TT{g + 1}: multi-unit, {len(sn)} sampled", fontsize=9)
        ax.set_ylim(-400, 200)
    for ax in axes.ravel():
        ax.set_xticks([])
        ax.axhline(0, color="0.85", lw=0.5, zorder=0)
    for ax in axes[:, 0]:
        ax.set_ylabel("µV")
    fig.tight_layout()
    fig.savefig(fig_path, dpi=110)


if __name__ == "__main__":
    main()
