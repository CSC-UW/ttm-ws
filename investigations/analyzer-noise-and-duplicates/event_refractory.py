"""Whether the spikes of a and b in a zero-lag coincident event are spikes of a and of b, from each
unit's other spikes: a real spike of b is followed and preceded by b's refractory gap and any bursts,
a spike of a third cell that the matcher rebuilds as a plus b is not. Same pairs as
consistent_waveform.py.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/event_refractory.py

For the spikes of a unit in coincident events (a and b within 3 samples), and for its spikes with no
spike of the other unit within 45 samples, the rate of the unit's other spikes at 1.5 to 3 ms and 3 to
5 ms either side, over the rate at 20 to 50 ms. The carry keeps no two spikes of a unit closer than
1.5 ms, so the first bin starts there.
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from double_fits import ISO, NEAR, RUN, partners

OUT = pathlib.Path(__file__).parent / "results"
NONE = np.iinfo(np.int64).min
BINS = {"1.5_3ms": (45, 90), "3_5ms": (90, 150)}
BASE = (600, 1500)


def rate(train, centers, lo, hi):
    n = (np.searchsorted(train, centers + hi, "left") - np.searchsorted(train, centers + lo, "left")
         + np.searchsorted(train, centers - lo, "right") - np.searchsorted(train, centers - hi, "right"))
    return n.sum() / (2 * (hi - lo) * len(centers))


def indices(train, centers):
    base = rate(train, centers, *BASE)
    return {k: float(rate(train, centers, lo, hi) / base) for k, (lo, hi) in BINS.items()}


def main():
    s = si.load(f"{RUN}/sorting")
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in s.unit_ids}
    pairs = pd.read_csv(OUT / "consistent_waveform.csv")[["kind", "unit_a", "unit_b", "tetrode", "near_over_expected"]]
    rows = []
    for kind, a, b, tt, ratio in pairs.itertuples(index=False):
        lag = partners(tr[a], tr[b], -NEAR, NEAR)
        co = lag != NONE
        ev = {a: tr[a][co], b: np.unique(tr[a][co] + lag[co])}
        row = {"kind": kind, "unit_a": a, "unit_b": b, "tetrode": tt, "near_over_expected": ratio, "events": int(co.sum())}
        for u, v, name in ((a, b, "a"), (b, a, "b")):
            alone = tr[u][partners(tr[u], tr[v], -ISO, ISO) == NONE]
            for k, x in indices(tr[u], ev[u]).items():
                row[f"{name}_{k}_coincident"] = x
            for k, x in indices(tr[u], alone).items():
                row[f"{name}_{k}_alone"] = x
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "event_refractory.csv", index=False)
    pd.set_option("display.width", 250)
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
