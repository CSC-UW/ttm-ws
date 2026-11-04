"""Whether the zero-lag coincident events of same-tetrode pairs on TT5 to TT12 of the final TTM-003-48H
sorting look like one waveform of their own, which the matcher rebuilds from two templates, or like two
spikes whose amplitudes vary independently. Pairs with an excess above 1 % (zero_lag.py) are compared
with pairs whose coincidences within 3 samples are 0.85 to 1.15 times the expected rate.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/consistent_waveform.py

Snippets (1 ms before to 2 ms after, µV) are read from the int16 binary. A and B are the mean spikes of
1,000 spikes of a and of b with no partner within 45 samples. Each coincident event (up to 3,000, a's
spike with one of b within 3 samples, lag L) is fitted as alpha * A + beta * B shifted by L. Reported:

- the coincidence rate at each lag from -3 to +3 samples over the expected rate per sample;
- the correlation of alpha and beta after removing each 1800 s window's mean, for coincident events and,
  as the correlation that noise alone gives, for 1,000 other spikes of a without a partner;
- the mean residual energy of coincident events against the sum of those of a's and b's spikes without a
  partner minus that of noise snippets (1,000 random times outside the tetrode's bad intervals with no
  spike of its units within 90 samples), which is what two independent spikes would leave;
- the mean coincident waveform at the most common lag fitted with A and B, and with A, B and their time
  derivatives (sub-sample misalignment), and the cosine between the mean residuals of two halves of the
  events.
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from double_fits import AFTER, BEFORE, ISO, LIVE, NEAR, RUN, partners, shift

OUT = pathlib.Path(__file__).parent / "results"
BAD = f"{RUN.split('/sortings_carry')[0]}/bad_periods/bad_periods.csv"
NONE = np.iinfo(np.int64).min
N_TEMPLATE, N_TEST, N_CO = 1000, 1000, 3000
WIN = 1800 * 30000
CLEAR = 90
OFF = np.arange(-BEFORE, AFTER)


def snips(raw, frames, cols, gain):
    return raw[np.asarray(frames, np.int64)[:, None] + OFF][:, :, cols].astype(np.float32) * gain


def fit(X, cols):
    """Least-squares amplitudes of the templates ``cols`` for each snippet of ``X``, and residual energy."""
    M = np.stack([c.ravel() for c in cols], 1)
    Xf = X.reshape(len(X), -1)
    coef = Xf @ np.linalg.pinv(M).T
    return coef, ((Xf - coef @ M.T) ** 2).sum(1)


def window_corr(t, x, y):
    d = pd.DataFrame({"w": t // WIN, "x": x, "y": y})
    g = d.groupby("w")
    keep = g["x"].transform("size") >= 3
    dx, dy = (d["x"] - g["x"].transform("mean"))[keep], (d["y"] - g["y"].transform("mean"))[keep]
    return float(np.corrcoef(dx, dy)[0, 1])


def left(W, cols):
    M = np.stack([c.ravel() for c in cols], 1)
    coef, *_ = np.linalg.lstsq(M, W.ravel(), rcond=None)
    r = W.ravel() - M @ coef
    return float((r ** 2).sum() / (W ** 2).sum()), r


def main():
    rng = np.random.default_rng(2)
    rec = si.load(f"{RUN}/binary")
    s = si.load(f"{RUN}/sorting")
    n_frames = rec.get_num_frames()
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", dtype="<i2", mode="r").reshape(n_frames, -1)
    gains = rec.get_channel_gains().astype(np.float32)
    rec_group = np.asarray(rec.get_property("group"))
    group = dict(zip(s.unit_ids, np.asarray(s.get_property("group"))))
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in s.unit_ids}
    tr = {u: t[(t >= BEFORE + 3) & (t < n_frames - AFTER - 3)] for u, t in tr.items()}
    probe = rng.integers(BEFORE + CLEAR, n_frames - AFTER - CLEAR, 50000)
    check = rec.get_traces(start_frame=int(probe[0]) - BEFORE, end_frame=int(probe[0]) + AFTER, return_in_uV=True)
    assert np.allclose(snips(raw, probe[:1], slice(None), gains)[0], check)
    bad = pd.read_csv(BAD)

    z = pd.read_csv(OUT / "zero_lag_pairs.csv")
    z = z[(z["pair"] == "same tetrode") & z["tetrode_a"].isin(LIVE)]
    z = z[(z["excess_of_smaller_unit"] > 0.01) | z["near_over_expected"].between(0.85, 1.15)]
    rows = []
    for a, b, tt, ratio, excess, expected in z[["unit_a", "unit_b", "tetrode_a", "near_over_expected",
                                                "excess_of_smaller_unit", "expected"]].itertuples(index=False):
        g = group[a]
        cols = np.flatnonzero(rec_group == g)
        gain = gains[cols]
        lag = partners(tr[a], tr[b], -NEAR, NEAR)
        co = np.flatnonzero(lag != NONE)
        iso_a = np.flatnonzero(partners(tr[a], tr[b], -ISO, ISO) == NONE)
        iso_b = np.flatnonzero(partners(tr[b], tr[a], -ISO, ISO) == NONE)
        pa = rng.choice(iso_a, N_TEMPLATE + N_TEST, replace=False)
        pb = rng.choice(iso_b, N_TEMPLATE + N_TEST, replace=False)
        A = snips(raw, tr[a][pa[:N_TEMPLATE]], cols, gain).mean(0)
        B = snips(raw, tr[b][pb[:N_TEMPLATE]], cols, gain).mean(0)

        ev = np.sort(rng.choice(co, min(N_CO, len(co)), replace=False))
        L, t_ev = lag[ev], tr[a][ev]
        X = snips(raw, t_ev, cols, gain)
        coef, res_co = np.empty((len(ev), 2), np.float32), np.empty(len(ev), np.float32)
        for l in np.unique(L):
            m = L == l
            coef[m], res_co[m] = fit(X[m], [A, shift(B, int(l))])
        ta = tr[a][pa[N_TEMPLATE:]]
        coef_a, res_a = fit(snips(raw, ta, cols, gain), [A, B])
        _, res_b = fit(snips(raw, tr[b][pb[N_TEMPLATE:]], cols, gain), [B, A])

        unit_trains = [tr[u] for u in s.unit_ids if group[u] == g]
        quiet = np.ones(len(probe), bool)
        for t in unit_trains:
            quiet &= np.searchsorted(t, probe - CLEAR, "left") == np.searchsorted(t, probe + CLEAR, "right")
        for lo, hi in zip(bad.loc[bad["group"] == g, "start_s"], bad.loc[bad["group"] == g, "end_s"]):
            quiet &= (probe < lo * 30000 - CLEAR) | (probe > hi * 30000 + CLEAR)
        _, res_0 = fit(snips(raw, probe[quiet][:N_TEST], cols, gain), [A, B])

        counts = np.array([(lag[co] == l).sum() for l in range(-NEAR, NEAR + 1)])
        mode = int(np.argmax(counts)) - NEAR
        at = L == mode
        W = X[at].mean(0)
        Bm = shift(B, mode)
        left2, _ = left(W, [A, Bm])
        left4, _ = left(W, [A, Bm, np.gradient(A, axis=0), np.gradient(Bm, axis=0)])
        halves = []
        for h in (0, 1):
            _, r = left(X[at][h::2].mean(0), [A, Bm])
            halves.append(r)
        row = {"kind": "excess" if excess > 0.01 else "about chance", "unit_a": a, "unit_b": b, "tetrode": tt,
               "near_over_expected": ratio, "coincident": len(co), "events": len(ev), "lag_mode": mode,
               "events_at_mode": int(at.sum())}
        row.update({f"lag_{l:+d}": float(c / (expected / (2 * NEAR + 1))) for l, c in zip(range(-NEAR, NEAR + 1), counts)})
        row.update({
            "alpha_median": float(np.median(coef[:, 0])), "beta_median": float(np.median(coef[:, 1])),
            "corr_coincident": window_corr(t_ev, coef[:, 0], coef[:, 1]),
            "corr_a_alone": window_corr(ta, coef_a[:, 0], coef_a[:, 1]),
            "res_coincident": float(res_co.mean()), "res_a": float(res_a.mean()), "res_b": float(res_b.mean()),
            "res_noise": float(res_0.mean()), "noise_snippets": int(min(quiet.sum(), N_TEST)),
            "res_over_two_spikes": float(res_co.mean() / (res_a.mean() + res_b.mean() - res_0.mean())),
            "left_mean_ab": left2, "left_mean_ab_derivatives": left4,
            "half_cosine": float(halves[0] @ halves[1] / np.linalg.norm(halves[0]) / np.linalg.norm(halves[1])),
        })
        rows.append(row)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "consistent_waveform.csv", index=False)
    pd.set_option("display.width", 300)
    print(d.drop(columns=["res_coincident", "res_a", "res_b", "res_noise"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
