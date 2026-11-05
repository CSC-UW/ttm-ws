"""Whether unit 13's coincident events with its four partners on TT6 (12, 14, 15, 16) are one waveform or
13's spike plus each partner's.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/unit13_partners.py

For each partner x, the events are 13's spikes with a spike of x at the most common lag within 3 samples
and no spike of any other TT6 unit within 45 samples; up to 5,000 are averaged, aligned on 13's spike (W_x).
Mean spikes T_u are averages of 2,000 spikes of u with no spike of another TT6 unit within 45 samples.

If the events are one source that the matcher builds from 13 and whichever partner fits, W_x is the same
waveform for every x, apart from the chance coincidences among the events. If they are 13's spike plus
x's, W_x differs between partners as S_x = T_13 + T_x (shifted by the lag) does. For each two partners:

- d_w: the energy of W_x - W_y over the mean energy of the two, at the best whole-sample shift within 3
  samples (samples 5 to 84 of the 90), and noise, what the events' spread alone gives to d_w;
- d_s: the same for S_x and S_y;
- d_f: the same for alpha T_13 + beta T_x, each W's least-squares fit;
- left_w, left_s: the fraction of W_x's (S_x's) energy that W_y (S_y) and its time derivative do not fit,
  averaged over both directions, which allows for a scale and a sub-sample shift.

Two hubs on TT6 whose partners coincide with them at about the expected rate (16 with 12 and 14, 12 with
15 and 16) are positive controls: their events are two spikes, so d_w should follow d_s there.
"""

import itertools
import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from consistent_waveform import snips
from double_fits import AFTER, BEFORE, ISO, NEAR, RUN, partners, shift

OUT = pathlib.Path(__file__).parent / "results"
NONE = np.iinfo(np.int64).min
N_TEMPLATE, N_EVENTS = 2000, 5000
INNER = slice(5, 85)
HUBS = {13: [12, 14, 15, 16], 16: [12, 14], 12: [15, 16]}


def near(t, trains, lo, hi):
    m = np.zeros(len(t), bool)
    for u in trains:
        m |= np.searchsorted(u, t + lo, "left") < np.searchsorted(u, t + hi, "right")
    return m


def dist(x, y):
    """Energy of x - y at the best whole-sample shift of y within 3 samples, over their mean energy."""
    norm = ((x[INNER] ** 2).sum() + (y[INNER] ** 2).sum()) / 2
    return min(float(((x[INNER] - shift(y, s)[INNER]) ** 2).sum()) for s in range(-NEAR, NEAR + 1)) / norm


def left(x, y):
    """Fraction of x's energy that y and its time derivative leave, by least squares."""
    M = np.stack([y[INNER].ravel(), np.gradient(y, axis=0)[INNER].ravel()], 1)
    coef, *_ = np.linalg.lstsq(M, x[INNER].ravel(), rcond=None)
    r = x[INNER].ravel() - M @ coef
    return float((r ** 2).sum() / (x[INNER] ** 2).sum())


def main():
    rng = np.random.default_rng(4)
    rec = si.load(f"{RUN}/binary")
    s = si.load(f"{RUN}/sorting")
    n_frames = rec.get_num_frames()
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", dtype="<i2", mode="r").reshape(n_frames, -1)
    gains = rec.get_channel_gains().astype(np.float32)
    group = dict(zip(s.unit_ids, np.asarray(s.get_property("group"))))
    g = group[13]
    cols = np.flatnonzero(np.asarray(rec.get_property("group")) == g)
    gain = gains[cols]
    tt6 = [u for u in s.unit_ids if group[u] == g]
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in tt6}
    tr = {u: t[(t >= BEFORE + NEAR) & (t < n_frames - AFTER - NEAR)] for u, t in tr.items()}

    T = {}
    for u in tt6:
        clean = np.flatnonzero(~near(tr[u], [tr[v] for v in tt6 if v != u], -ISO, ISO))
        T[u] = snips(raw, tr[u][rng.choice(clean, N_TEMPLATE, replace=False)], cols, gain).mean(0)

    pair_rows, comp_rows = [], []
    for hub, xs in HUBS.items():
        W, S, F, noise_var = {}, {}, {}, {}
        for x in xs:
            lag = partners(tr[hub], tr[x], -NEAR, NEAR)
            others = [tr[v] for v in tt6 if v not in (hub, x)]
            clean = ~near(tr[hub], others, -ISO, ISO)
            co = (lag != NONE) & clean
            mode = int(pd.Series(lag[co]).mode().iloc[0])
            ev = np.flatnonzero(co & (lag == mode))
            ev = rng.choice(ev, min(N_EVENTS, len(ev)), replace=False)
            X = snips(raw, tr[hub][ev], cols, gain)
            W[x] = X.mean(0)
            noise_var[x] = X.var(0, ddof=1) / len(X)
            Tx = shift(T[x], mode)
            S[x] = T[hub] + Tx
            M = np.stack([T[hub].ravel(), Tx.ravel()], 1)
            (al, be), *_ = np.linalg.lstsq(M, W[x].ravel(), rcond=None)
            F[x] = al * T[hub] + be * Tx
            r = W[x] - F[x]
            pair_rows.append({"hub": hub, "partner": x, "lag_mode": mode, "events": len(ev), "alpha": float(al),
                              "beta": float(be), "left_fit": float((r ** 2).sum() / (W[x] ** 2).sum())})
        for x, y in itertools.combinations(xs, 2):
            norm = ((W[x][INNER] ** 2).sum() + (W[y][INNER] ** 2).sum()) / 2
            comp_rows.append({
                "hub": hub, "x": x, "y": y,
                "d_w": dist(W[x], W[y]), "noise": float((noise_var[x] + noise_var[y])[INNER].sum() / norm),
                "d_s": dist(S[x], S[y]), "d_f": dist(F[x], F[y]),
                "left_w": (left(W[x], W[y]) + left(W[y], W[x])) / 2,
                "left_s": (left(S[x], S[y]) + left(S[y], S[x])) / 2,
                "partner_cosine_raw": float((T[x] * T[y]).sum() / np.linalg.norm(T[x]) / np.linalg.norm(T[y])),
            })
    p = pd.DataFrame(pair_rows)
    c = pd.DataFrame(comp_rows)
    c["d_w_over_d_s"] = (c["d_w"] - c["noise"]) / c["d_s"]
    p.to_csv(OUT / "unit13_partners_pairs.csv", index=False)
    c.to_csv(OUT / "unit13_partners.csv", index=False)
    pd.set_option("display.width", 200)
    print(p.round(3).to_string(index=False))
    print(c.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
