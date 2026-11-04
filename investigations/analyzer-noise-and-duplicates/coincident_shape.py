"""Shape of the zero-lag coincident events of double_fits.py's pairs: is their mean waveform a scaled
copy of one unit's mean spike (one large spike whose residual a second template was fitted to), or
does it need both units' mean spikes?

    .venv/bin/python investigations/analyzer-noise-and-duplicates/coincident_shape.py

Per pair, the coincident events at the most common lag are averaged (aligned on a's spike). That mean
waveform W is fitted by least squares with free amplitudes as alpha * A, beta * B (B shifted by the
lag) and alpha * A + beta * B, where A and B are the units' mean spikes over spikes with no partner
within 1.5 ms; the fraction of W's energy each fit leaves is reported, and the mean number of the
tetrode's other units with a spike within 1.5 ms of a's spike in those events.
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from double_fits import ISO, NEAR, RUN, N_EVENTS, partners, shift, snippets

OUT = pathlib.Path(__file__).parent / "results"
NONE = np.iinfo(np.int64).min


def fit(W, cols):
    M = np.stack([c.ravel() for c in cols], 1)
    coef, *_ = np.linalg.lstsq(M, W.ravel(), rcond=None)
    r = W.ravel() - M @ coef
    return coef, float((r ** 2).sum() / (W ** 2).sum())


def main():
    rng = np.random.default_rng(1)
    rec = si.load(f"{RUN}/binary")
    s = si.load(f"{RUN}/sorting")
    group = dict(zip(s.unit_ids, np.asarray(s.get_property("group"))))
    rec_group = np.asarray(rec.get_property("group"))
    n = rec.get_num_frames()
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in s.unit_ids}
    tr = {u: t[(t >= 100) & (t < n - 100)] for u, t in tr.items()}
    pairs = pd.read_csv(OUT / "double_fits.csv")
    pairs = pairs[pairs["events"] == "coincident"][["pairs", "unit_a", "unit_b", "tetrode", "cosine"]]
    rows = []
    for kind, a, b, tt, cos in pairs.itertuples(index=False):
        ch = np.asarray(rec.channel_ids)[rec_group == group[a]]
        lag = partners(tr[a], tr[b], -NEAR, NEAR)
        co = lag != NONE
        mode = int(pd.Series(lag[co]).mode().iloc[0])
        ev = np.flatnonzero(lag == mode)
        ev = rng.choice(ev, min(N_EVENTS, len(ev)), replace=False)
        iso_a = partners(tr[a], tr[b], -ISO, ISO) == NONE
        iso_b = partners(tr[b], tr[a], -ISO, ISO) == NONE
        A = snippets(rec, rng.choice(tr[a][iso_a], N_EVENTS, replace=False), ch).mean(0)
        B = shift(snippets(rec, rng.choice(tr[b][iso_b], N_EVENTS, replace=False), ch).mean(0), mode)
        W = snippets(rec, tr[a][ev], ch).mean(0)
        (al,), ra = fit(W, [A])
        (be,), rb = fit(W, [B])
        (al2, be2), rab = fit(W, [A, B])
        others = [u for u in s.unit_ids if group[u] == group[a] and u not in (a, b)]
        third = sum(int((partners(tr[a][ev], tr[u], -ISO, ISO) != NONE).sum()) for u in others) / len(ev)
        rows.append({"pairs": kind, "unit_a": a, "unit_b": b, "tetrode": tt, "cosine": cos,
                     "lag_mode": mode, "share_at_mode": float((lag[co] == mode).mean()), "events": len(ev),
                     "alpha_only_a": al, "left_only_a": ra, "beta_only_b": be, "left_only_b": rb,
                     "alpha_both": al2, "beta_both": be2, "left_both": rab, "third_unit_within_1p5ms": third})
        print(rows[-1], flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "coincident_shape.csv", index=False)
    pd.set_option("display.width", 250)
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
