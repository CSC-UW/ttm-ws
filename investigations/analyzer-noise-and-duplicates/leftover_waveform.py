"""What the two mean spikes leave of the mean coincident waveform, for the pairs of consistent_waveform.py.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/leftover_waveform.py [figure.png]

Per pair, W is the mean of up to 5,000 coincident events at the most common lag (aligned on a's spike),
A and B the mean spikes of 2,000 spikes of a and of b with no partner within 45 samples, and R the
residual of W fitted by least squares with A and B (shifted by the lag). Reported for R: its main channel
and peak time against A's and B's, its peak over A's, how its energy spreads over the four channels
(cosine of the per-channel RMS against A's and B's), and its max-shift cosine with A, B and the mean spike
(1,000 random spikes) of every other unit of the tetrode.

Whether R is a change in a's or b's own spike when the unit fires in a burst: D is the mean of the unit's
spikes without a partner whose previous spike of the unit is under 10 ms before, minus the mean of those
whose previous spike is over 50 ms before. W is refitted with A, B, D_a and D_b (B and D_b shifted by the
lag), and with their time derivatives too. Also the share of coincident and of partnerless spikes whose
previous spike of the unit is under 10 ms before.

Whether R is a difference in width: W is refitted with A, B and their first and second time derivatives.
As a calibration, the mean of the 10 % largest (by their least-squares amplitude on A) of 20,000 spikes of
a without a partner is fitted with A, and with A and its first and second derivatives.

Pairs are grouped as in the README: sharp where the best lag within 3 samples is at least 3.5 times the
lowest, broad under 2 times, in between otherwise, and the controls at about chance.
"""

import pathlib
import sys

import matplotlib
import numpy as np
import pandas as pd
import spikeinterface as si

from consistent_waveform import snips
from double_fits import AFTER, BEFORE, ISO, NEAR, RUN, partners, shift
from ttm.tracking import cosine_from_templates

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = pathlib.Path(__file__).parent / "results"
NONE = np.iinfo(np.int64).min
N_TEMPLATE, N_EVENTS, N_OTHER, N_BIG = 2000, 5000, 1000, 20000
BURST, FIRST = 300, 1500


def fit(W, cols):
    M = np.stack([c.ravel() for c in cols], 1)
    coef, *_ = np.linalg.lstsq(M, W.ravel(), rcond=None)
    r = (W.ravel() - M @ coef).reshape(W.shape)
    return coef, r, float((r ** 2).sum() / (W ** 2).sum())


def prev_isi(t, idx):
    return np.where(idx > 0, t[idx] - t[np.maximum(idx - 1, 0)], np.iinfo(np.int64).max)


def burst_diff(raw, t, iso, cols, gain, rng):
    p = prev_isi(t, iso)
    burst, first = iso[p < BURST], iso[p > FIRST]
    mb = snips(raw, t[rng.choice(burst, min(N_TEMPLATE, len(burst)), replace=False)], cols, gain).mean(0)
    mf = snips(raw, t[rng.choice(first, min(N_TEMPLATE, len(first)), replace=False)], cols, gain).mean(0)
    return mb - mf, len(burst)


def group_of(row):
    lags = row[[f"lag_{l:+d}" for l in range(-NEAR, NEAR + 1)]].astype(float)
    r = lags.max() / lags.min()
    if row["kind"] == "about chance":
        return "control"
    return "sharp" if r >= 3.5 else "broad" if r < 2 else "in between"


def main():
    fig_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None
    rng = np.random.default_rng(3)
    rec = si.load(f"{RUN}/binary")
    s = si.load(f"{RUN}/sorting")
    n_frames = rec.get_num_frames()
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", dtype="<i2", mode="r").reshape(n_frames, -1)
    gains = rec.get_channel_gains().astype(np.float32)
    rec_group = np.asarray(rec.get_property("group"))
    group = dict(zip(s.unit_ids, np.asarray(s.get_property("group"))))
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in s.unit_ids}
    tr = {u: t[(t >= BEFORE + NEAR) & (t < n_frames - AFTER - NEAR)] for u, t in tr.items()}
    pairs = pd.read_csv(OUT / "consistent_waveform.csv")
    pairs["group"] = pairs.apply(group_of, axis=1)

    rows, shapes = [], {}
    for _, p in pairs.iterrows():
        a, b, mode, g = int(p.unit_a), int(p.unit_b), int(p.lag_mode), group[int(p.unit_a)]
        cols = np.flatnonzero(rec_group == g)
        gain = gains[cols]
        lag = partners(tr[a], tr[b], -NEAR, NEAR)
        at = np.flatnonzero(lag == mode)
        co = np.flatnonzero(lag != NONE)
        iso_a = np.flatnonzero(partners(tr[a], tr[b], -ISO, ISO) == NONE)
        iso_b = np.flatnonzero(partners(tr[b], tr[a], -ISO, ISO) == NONE)
        A = snips(raw, tr[a][rng.choice(iso_a, N_TEMPLATE, replace=False)], cols, gain).mean(0)
        B = snips(raw, tr[b][rng.choice(iso_b, N_TEMPLATE, replace=False)], cols, gain).mean(0)
        Da, n_burst_a = burst_diff(raw, tr[a], iso_a, cols, gain, rng)
        Db, n_burst_b = burst_diff(raw, tr[b], iso_b, cols, gain, rng)
        W = snips(raw, tr[a][rng.choice(at, min(N_EVENTS, len(at)), replace=False)], cols, gain).mean(0)
        Bm, Dbm = shift(B, mode), shift(Db, mode)
        (al, be), R, left = fit(W, [A, Bm])
        _, _, left_burst = fit(W, [A, Bm, Da, Dbm])
        dA, dB = np.gradient(A, axis=0), np.gradient(Bm, axis=0)
        _, _, left_all = fit(W, [A, Bm, dA, dB, Da, Dbm])
        _, Rd, left_deriv = fit(W, [A, Bm, dA, dB])
        d2A, d2B = np.gradient(dA, axis=0), np.gradient(dB, axis=0)
        _, R2, left_second = fit(W, [A, Bm, dA, dB, d2A, d2B])
        X = snips(raw, tr[a][rng.choice(iso_a, min(N_BIG, len(iso_a)), replace=False)], cols, gain)
        amp = X.reshape(len(X), -1) @ A.ravel() / (A ** 2).sum()
        big = X[amp >= np.quantile(amp, 0.9)].mean(0)
        (big_alpha,), Rbig, left_big = fit(big, [A])
        _, _, left_big_second = fit(big, [A, dA, d2A])

        def foot(x):
            return np.sqrt((x ** 2).mean(0))

        def main_ch(x):
            return int(np.argmax(np.abs(x).max(0)))

        def peak_t(x):
            c = main_ch(x)
            return int(np.argmax(np.abs(x[:, c]))) - BEFORE

        others = {}
        for u in s.unit_ids:
            if group[u] == g and u not in (a, b):
                t = tr[u]
                others[u] = snips(raw, t[rng.choice(len(t), min(N_OTHER, len(t)), replace=False)], cols, gain).mean(0)
        cos_other = {u: float(cosine_from_templates(R, T, max_shift_samples=10)) for u, T in others.items()}
        best_other = max(cos_other, key=cos_other.get) if cos_other else None
        burst_a = prev_isi(tr[a], np.arange(len(tr[a])))
        burst_b = prev_isi(tr[b], np.arange(len(tr[b])))
        b_idx = np.searchsorted(tr[b], tr[a][at] + mode)
        row = {
            "group": p.group, "unit_a": a, "unit_b": b, "tetrode": int(g) + 1, "lag_mode": mode, "events": len(at),
            "alpha": float(al), "beta": float(be), "left": left, "left_derivatives": left_deriv,
            "left_burst": left_burst, "left_burst_derivatives": left_all, "left_second_derivatives": left_second,
            "big_a_alpha": float(big_alpha), "big_a_left": left_big, "big_a_left_second_derivatives": left_big_second,
            "cos_big_a": float(cosine_from_templates(R, Rbig, max_shift_samples=10)),
            "r_main_channel": main_ch(R), "a_main_channel": main_ch(A), "b_main_channel": main_ch(Bm),
            "r_peak_sample": peak_t(R), "a_peak_sample": peak_t(A), "b_peak_sample": peak_t(Bm),
            "r_peak_uv": float(R[peak_t(R) + BEFORE, main_ch(R)]),
            "r_max_sample": int(np.argmax(R[:, main_ch(R)])) - BEFORE, "r_min_sample": int(np.argmin(R[:, main_ch(R)])) - BEFORE,
            "w_trough_over_fit": float(W[:, main_ch(R)].min() / (al * A + be * Bm)[:, main_ch(R)].min()),
            "a_peak_uv": float(A[peak_t(A) + BEFORE, main_ch(A)]),
            "r_peak_over_a_peak": float(np.abs(R).max() / np.abs(A).max()),
            "footprint_cos_a": float(foot(R) @ foot(A) / np.linalg.norm(foot(R)) / np.linalg.norm(foot(A))),
            "footprint_cos_b": float(foot(R) @ foot(B) / np.linalg.norm(foot(R)) / np.linalg.norm(foot(B))),
            "cos_a": float(cosine_from_templates(R, A, max_shift_samples=10)),
            "cos_b": float(cosine_from_templates(R, Bm, max_shift_samples=10)),
            "cos_burst_a": float(cosine_from_templates(R, Da, max_shift_samples=10)),
            "cos_burst_b": float(cosine_from_templates(R, Dbm, max_shift_samples=10)),
            "best_other_unit": best_other, "cos_best_other": cos_other.get(best_other, np.nan),
            "burst_share_a_coincident": float((burst_a[co] < BURST).mean()),
            "burst_share_a_alone": float((burst_a[iso_a] < BURST).mean()),
            "burst_share_b_coincident": float((burst_b[np.clip(b_idx, 0, len(tr[b]) - 1)] < BURST).mean()),
            "burst_share_b_alone": float((burst_b[iso_b] < BURST).mean()),
            "burst_spikes_a": n_burst_a, "burst_spikes_b": n_burst_b,
        }
        rows.append(row)
        shapes[(a, b)] = (W, al * A + be * Bm, R, p.group, mode, Rbig, main_ch(R))
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(OUT / "leftover_waveform.csv", index=False)
    pd.set_option("display.width", 300)
    print(d.round(3).to_string(index=False))

    if fig_path is not None:
        order = ["sharp", "in between", "broad", "control"]
        keys = sorted(shapes, key=lambda k: (order.index(shapes[k][3]), -d.loc[(d.unit_a == k[0]) & (d.unit_b == k[1]), "left"].iloc[0]))
        t_ms = np.arange(-BEFORE, AFTER) / 30.0
        fig, axes = plt.subplots(len(keys), 4, figsize=(12, 1.6 * len(keys)), sharex=True)
        for i, k in enumerate(keys):
            W, fitted, R, grp, mode, _, _ = shapes[k]
            lim = np.abs(W).max() * 1.05
            for c in range(4):
                ax = axes[i, c]
                ax.plot(t_ms, W[:, c], color="k", lw=1)
                ax.plot(t_ms, fitted[:, c], color="0.6", lw=1, ls="--")
                ax.plot(t_ms, 5 * R[:, c], color="tab:red", lw=1)
                ax.set_ylim(-lim, lim)
                ax.axvline(0, color="0.85", lw=0.5, zorder=0)
                ax.tick_params(labelsize=6)
                if c == 0:
                    ax.set_ylabel(f"{k[0]}/{k[1]}\n{grp}, lag {mode:+d}", fontsize=7)
                if i == 0:
                    ax.set_title(f"channel {c + 1}", fontsize=8)
        axes[-1, 0].set_xlabel("ms from a's spike", fontsize=7)
        fig.suptitle("mean coincident waveform (black), alpha A + beta B (grey), residual x5 (red); µV", fontsize=9)
        fig.tight_layout()
        fig.savefig(fig_path, dpi=110)

        zoom = [k for k in keys if shapes[k][3] == "sharp"] + [k for k in keys if shapes[k][3] == "control"][:4]
        fig, axes = plt.subplots(3, 4, figsize=(13, 8))
        for ax, k in zip(axes.ravel(), zoom):
            W, fitted, R, grp, mode, Rbig, c = shapes[k]
            win = (t_ms > -0.5) & (t_ms < 1.0)
            ax.plot(t_ms[win], W[win, c], color="k", lw=1, label="W")
            ax.plot(t_ms[win], fitted[win, c], color="0.6", lw=1, ls="--", label="alpha A + beta B")
            ax.plot(t_ms[win], 5 * R[win, c], color="tab:red", lw=1, label="residual x5")
            ax.plot(t_ms[win], 5 * Rbig[win, c], color="tab:blue", lw=1, label="large a spikes minus A, x5")
            ax.axvline(0, color="0.85", lw=0.5, zorder=0)
            ax.set_title(f"{k[0]}/{k[1]} {grp}, lag {mode:+d}, channel {c + 1}", fontsize=8)
            ax.tick_params(labelsize=6)
        axes[0, 0].legend(fontsize=6)
        fig.tight_layout()
        fig.savefig(fig_path.with_name(fig_path.stem + "_zoom.png"), dpi=110)


if __name__ == "__main__":
    main()
