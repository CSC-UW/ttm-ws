"""Figures for the multi-unit gate: the cosine distributions per tetrode, and a gallery of events per
tetrode and cosine band.

    .venv/bin/python investigations/multi-unit-gate/figure.py <gate out dir> <figure dir> [--sub]

With --sub, the scores, best units, scales and gallery traces are those after subtracting the neighbouring
sorted spikes (the residual traces), from the results files with the suffix _sub, and every figure name
gets the suffix _sub.

cosine_distributions.png, from results/: per tetrode of TT5 to TT12, the best cosine of the reference
events (blue), the unclaimed events (black) and the background stretches (grey), each as a share per 0.02
bin, with the cutoffs losing 1 % (dotted), 5 % (dashed) and 10 % (dash-dot) and 97's 0.55 (red).

lags.png, from results/lag_hist.csv: per tetrode of TT5 to TT12, the sorted spikes of the tetrode at each
lag from the event (spike minus event, -3 to 3 ms), over their rate 20 to 50 ms away, for reference events and for
unclaimed events grouped by the tetrode's cutoffs losing 1, 5 and 10 %.

gallery_TTn.png, from gate.py's traces: columns are bands of best cosine; rows are random unclaimed
events of 7 to 9 times the noise, unclaimed events of 9 or more, and reference events of 7 to 9, the
tetrode's four wires side by side in µV (black, or blue for reference), with the best template scaled to
the event by least squares (red) and the times of the tetrode's other sorted spikes in the window (grey
ticks; a reference event's own claiming spike is left out). Each title gives the amplitude in units of noise, the best cosine, and
the number of other tetrodes with an event within 0.5 ms.
"""

import pathlib
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gate import AFTER, BEFORE, CLAIM, Setup

RESULTS = pathlib.Path(__file__).parent / "results"
GAIN = 0.195
N_SAMP = 90
BANDS = ((-1.0, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 1.01))
SUB = "--sub" in sys.argv
SFX = "_sub" if SUB else ""
ROWS = (("unclaimed, 7-9", 0, False, 4), ("unclaimed, 9+", 0, True, 2), ("reference, 7-9", 1, False, 3))


def flat(w):
    """(..., 90, 4) -> (..., 4 * 92), the wires side by side with a gap."""
    pad = np.full(w.shape[:-2] + (2, 4), np.nan)
    return np.concatenate([w, pad], axis=-2).swapaxes(-1, -2).reshape(*w.shape[:-2], -1)


def distributions(fig_dir):
    h = pd.read_csv(RESULTS / f"cosine_hist{SFX}.csv")
    t = pd.read_csv(RESULTS / f"tetrodes{SFX}.csv").set_index("tetrode")
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharex=True)
    for ax, tt in zip(axes.ravel(), range(5, 13)):
        for kind, color in (("reference", "tab:blue"), ("unclaimed", "k"), ("background", "0.6")):
            s = h[(h["tetrode"] == tt) & (h["kind"] == kind)]
            ax.step(s["bin_lo"] + 0.01, s["count"] / s["count"].sum(), where="mid", color=color, label=kind)
        ax.axvline(t.loc[tt, "cut_1"], color="k", ls=":", lw=1, label="loses 1 %")
        ax.axvline(t.loc[tt, "cut_5"], color="k", ls="--", lw=1, label="loses 5 %")
        ax.axvline(t.loc[tt, "cut_10"], color="k", ls="-.", lw=1, label="loses 10 %")
        ax.axvline(0.55, color="r", lw=1, label="0.55")
        ax.set_title(f"TT{tt}: cutoffs {t.loc[tt, 'cut_1']:.2f} (1 %), {t.loc[tt, 'cut_5']:.2f} (5 %), "
                     f"{t.loc[tt, 'cut_10']:.2f} (10 %)", fontsize=9)
        ax.set_xlim(-0.3, 1)
    axes[0, 0].legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("best template cosine")
    for ax in axes[:, 0]:
        ax.set_ylabel("share per 0.02")
    fig.tight_layout()
    fig.savefig(fig_dir / f"cosine_distributions{SFX}.png", dpi=110)
    plt.close(fig)


def lags(fig_dir):
    h = pd.read_csv(RESULTS / f"lag_hist{SFX}.csv")
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharex=True)
    colors = {"reference": "tab:blue", "unclaimed below 1 %": "0.75", "unclaimed 1 to 5 %": "0.5",
              "unclaimed 5 to 10 %": "0.25", "unclaimed at or above 10 %": "k",
              "unclaimed, residual below 7": "tab:orange"}
    for ax, tt in zip(axes.ravel(), range(5, 13)):
        for kind, color in colors.items():
            s = h[(h["tetrode"] == tt) & (h["kind"] == kind)]
            s = s[s["lag"].abs() > CLAIM].dropna()
            if s.empty:
                continue
            for side in (s[s["lag"] < 0], s[s["lag"] > 0]):
                ax.plot(side["lag"] / 30, side["ratio"], color=color, lw=1,
                        label=kind if side["lag"].iloc[0] > 0 else None)
        ax.axhline(1, color="0.85", lw=0.5, zorder=0)
        ax.set_title(f"TT{tt}", fontsize=9)
    axes[0, 0].legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("sorted spike minus event (ms)")
    for ax in axes[:, 0]:
        ax.set_ylabel("spikes per lag / rate at 20-50 ms")
    fig.tight_layout()
    fig.savefig(fig_dir / f"lags{SFX}.png", dpi=110)
    plt.close(fig)


def galleries(d, fig_dir):
    t = pd.read_csv(RESULTS / f"tetrodes{SFX}.csv").set_index("tetrode")
    st = Setup()
    keys = ("group", "claims", "amp", "best_cos", "others", "alpha", "best_uid", "frame")
    src = {"best_cos": "sub_cos", "alpha": "sub_alpha", "best_uid": "sub_uid"} if SUB else {}
    sn = {k: [] for k in keys + ("trace", "tmpl")}
    for f in sorted(d.glob("w*.npz")):
        z = np.load(f)
        i = z["snip_idx"]
        for k in keys:
            sn[k].append(z[f"ev_{src.get(k, k)}"][i])
        sn["trace"].append(z["snip_res" if SUB else "snip"].astype(np.float32))
        tid = {int(u): j for j, u in enumerate(z["tmpl_ids"])}
        uid, alpha = z[f"ev_{src.get('best_uid', 'best_uid')}"][i], z[f"ev_{src.get('alpha', 'alpha')}"][i]
        sn["tmpl"].append(z["tmpl_T"][[tid[int(u)] for u in uid]] * alpha[:, None, None])
    sn = {k: np.concatenate(v) for k, v in sn.items()}
    rng = np.random.default_rng(0)
    x = np.arange(4 * (N_SAMP + 2))
    for g in range(4, 12):
        tt = g + 1
        n_rows = sum(r[3] for r in ROWS)
        fig, axes = plt.subplots(n_rows, len(BANDS), figsize=(2.6 * len(BANDS), 1.45 * n_rows))
        for j, (lo, hi) in enumerate(BANDS):
            marks = [f"{q} % cut {t.loc[tt, f'cut_{q}']:.2f}" for q in (1, 5, 10)
                     if lo <= t.loc[tt, f"cut_{q}"] < hi]
            if lo <= 0.55 < hi:
                marks.append("0.55")
            head = f"cosine {max(lo, 0):.1f} to {min(hi, 1):.1f}" if lo > -1 else "cosine below 0.2"
            axes[0, j].set_title(head + ("\n" + ", ".join(marks) if marks else "\n"), fontsize=8)
            row = 0
            for label, claims, big, n in ROWS:
                m = np.flatnonzero((sn["group"] == g) & (sn["claims"] == claims) & (sn["best_cos"] >= lo)
                                   & (sn["best_cos"] < hi) & ((sn["amp"] >= 9) == big))
                pick = rng.permutation(m)[:n]
                for i in range(n):
                    ax = axes[row + i, j]
                    ax.set_xticks([])
                    ax.tick_params(labelsize=6)
                    if j == 0:
                        ax.set_ylabel(label, fontsize=7)
                    if i >= pick.size:
                        ax.set_yticks([])
                        continue
                    e = pick[i]
                    ax.plot(x, flat(sn["trace"][e]) * GAIN, color="tab:blue" if claims else "k", lw=0.7)
                    ax.plot(x, flat(sn["tmpl"][e]) * GAIN, color="r", lw=0.7)
                    p = st.pool_f[g]
                    lag = p[np.searchsorted(p, sn["frame"][e] - BEFORE):np.searchsorted(p, sn["frame"][e] + AFTER)] \
                        - sn["frame"][e]
                    if claims:
                        lag = np.delete(lag, np.argmin(np.abs(lag)))
                    for lg in lag:
                        for wire in range(4):
                            ax.axvline(wire * (N_SAMP + 2) + BEFORE + lg, color="0.6", lw=0.6, ymax=0.12)
                    ax.text(0.01, 0.02, f"{sn['amp'][e]:.1f}, cos {sn['best_cos'][e]:.2f}, +{sn['others'][e]}",
                            transform=ax.transAxes, fontsize=6)
                row += n
        what = "residual traces after subtracting the neighbouring sorted spikes, " if SUB else ""
        fig.suptitle(f"TT{tt}: unclaimed (black) and reference (blue) events by best cosine; {what}best template "
                     "scaled to the event (red); amplitude in noise units, cosine, other tetrodes", fontsize=9)
        fig.tight_layout()
        fig.savefig(fig_dir / f"gallery{SFX}_TT{tt}.png", dpi=100)
        plt.close(fig)


def main():
    d, fig_dir = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    fig_dir.mkdir(parents=True, exist_ok=True)
    distributions(fig_dir)
    lags(fig_dir)
    galleries(d, fig_dir)


if __name__ == "__main__":
    main()
