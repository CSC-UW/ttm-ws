"""Cross-correlograms out to 1.5 ms at one-sample resolution for every same-tetrode pair on TT5 to TT12 of the
final TTM-003-48H sorting, and for zero_lag.py's 70 pairs of units on two different tetrodes among TT5 to
TT12 as a control: whether the coincidences within 0.1 ms sit on broader co-firing, and whether a peak
within 0.1 ms could be coincidences moved there from neighbouring lags.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/lag_profiles.py

For each pair, b's spikes at each lag from -45 to +45 samples around a's spikes, over the mean count per
lag 150 to 750 samples away on either side (the 5 to 25 ms that zero_lag.py takes the expected rate from).
Means over ranges of the absolute lag go to results/lag_profiles.csv, the ratio at each lag to
results/lag_profiles_by_lag.csv. Pairs carry the group leftover_waveform.py gave them (sharp, broad, in
between, control), "other", or "different tetrodes". Also, in units of the expected count per lag: the
excess within 3 samples
above the mean over 9 to 15 samples (excess_0_3), the shortfall at 4 to 8 samples below that mean
(shortfall_4_8, negative where there is none), and the excess over the expected rate from -45 to +45
samples (excess_45).
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from double_fits import LIVE, RUN

OUT = pathlib.Path(__file__).parent / "results"
FLANK = (150, 750)
RANGES = {"0_3": (0, 3), "4_8": (4, 8), "9_15": (9, 15), "16_30": (16, 30), "31_45": (31, 45)}


def lags(ta, tb, lo, hi):
    """b - a for every pair of spikes with lo <= b - a <= hi, for sorted frame arrays."""
    i0 = np.searchsorted(tb, ta + lo, "left")
    n = np.searchsorted(tb, ta + hi, "right") - i0
    keep = n > 0
    n, i0, a = n[keep], i0[keep], ta[keep]
    off = np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n)
    return tb[np.repeat(i0, n) + off] - np.repeat(a, n)


def main():
    s = si.load(f"{RUN}/sorting")
    group = dict(zip(s.unit_ids, np.asarray(s.get_property("group"))))
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in s.unit_ids}
    z = pd.read_csv(OUT / "zero_lag_pairs.csv")
    live = z["tetrode_a"].isin(LIVE) & z["tetrode_b"].isin(LIVE)
    same = z[(z["pair"] == "same tetrode") & live][["unit_a", "unit_b"]]
    cross = z[(z["pair"] == "different tetrodes") & live][["unit_a", "unit_b"]].assign(group="different tetrodes")
    grp = pd.read_csv(OUT / "leftover_waveform.csv")[["group", "unit_a", "unit_b"]]
    same = same.merge(grp, on=["unit_a", "unit_b"], how="left").fillna({"group": "other"})
    pairs = pd.concat([same, cross])[["group", "unit_a", "unit_b"]]
    rows, by_lag = [], []
    lag = np.arange(-45, 46)
    for grp, a, b in pairs.itertuples(index=False):
        L = lags(tr[a], tr[b], -FLANK[1], FLANK[1])
        per_lag = (np.abs(L) >= FLANK[0]).sum() / (2 * (FLANK[1] - FLANK[0] + 1))
        prof = np.bincount(L[np.abs(L) <= 45] + 45, minlength=91) / per_lag
        row = {"group": grp, "unit_a": a, "unit_b": b, "tetrode": int(group[a]) + 1,
               "peak": float(prof.max()), "peak_lag": int(lag[prof.argmax()])}
        row.update({k: float(prof[(np.abs(lag) >= lo) & (np.abs(lag) <= hi)].mean()) for k, (lo, hi) in RANGES.items()})
        level = row["9_15"]
        row["excess_0_3"] = float((prof[np.abs(lag) <= 3] - level).sum())
        row["shortfall_4_8"] = float((level - prof[(np.abs(lag) >= 4) & (np.abs(lag) <= 8)]).sum())
        row["excess_45"] = float((prof - 1).sum())
        rows.append(row)
        by_lag += [{"unit_a": a, "unit_b": b, "lag": int(l), "over_expected": float(p)} for l, p in zip(lag, prof)]
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "lag_profiles.csv", index=False)
    pd.DataFrame(by_lag).to_csv(OUT / "lag_profiles_by_lag.csv", index=False)
    pd.set_option("display.width", 200)
    print(d.sort_values(["group", "peak"], ascending=[False, False]).round(2).to_string(index=False))
    print(d.groupby("group")[list(RANGES)].agg(["min", "max"]).round(2).to_string())


if __name__ == "__main__":
    main()
