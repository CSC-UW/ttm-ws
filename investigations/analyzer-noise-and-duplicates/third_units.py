"""How often a third unit of the tetrode fires near a spike of unit a, for a's spikes that coincide with
unit b within 3 samples and for a's spikes without b, over every same-tetrode pair on TT5 to TT12 of
the final TTM-003-48H sorting.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/third_units.py

Spikes of a fall in three sets by where b's nearest spike is: within 3 samples (coincident), none
within 45 samples (lone), and none within 45 but one 150 to 750 samples away (flank, the 5 to 25 ms
range zero_lag.py takes the expected rate from). For each set, the share of a's spikes with a spike of
another unit of the tetrode (neither a nor b) within 45 samples, and within 3 samples.
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from double_fits import ISO, LIVE, NEAR, RUN

OUT = pathlib.Path(__file__).parent / "results"
FLANK = (150, 750)


def near(t, trains, lo, hi):
    m = np.zeros(len(t), bool)
    for u in trains:
        m |= np.searchsorted(u, t + lo, "left") < np.searchsorted(u, t + hi, "right")
    return m


def main():
    s = si.load(f"{RUN}/sorting")
    group = dict(zip(s.unit_ids, np.asarray(s.get_property("group"))))
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in s.unit_ids}
    z = pd.read_csv(OUT / "zero_lag_pairs.csv")
    z = z[(z["pair"] == "same tetrode") & z["tetrode_a"].isin(LIVE)]
    rows = []
    for a, b, tt, ratio, excess in z[["unit_a", "unit_b", "tetrode_a", "near_over_expected",
                                      "excess_of_smaller_unit"]].itertuples(index=False):
        ta, tb = tr[a], [tr[b]]
        others = [tr[u] for u in s.unit_ids if group[u] == group[a] and u not in (a, b)]
        within = near(ta, tb, -ISO, ISO)
        sets = {"coincident": near(ta, tb, -NEAR, NEAR), "lone": ~within,
                "flank": ~within & (near(ta, tb, *FLANK) | near(ta, tb, -FLANK[1], -FLANK[0]))}
        row = {"unit_a": a, "unit_b": b, "tetrode": tt, "near_over_expected": ratio, "excess_of_smaller_unit": excess,
               "other_units": len(others)}
        for name, m in sets.items():
            t = ta[m]
            row[f"n_{name}"] = len(t)
            row[f"third_45_{name}"] = float(near(t, others, -ISO, ISO).mean()) if len(t) else np.nan
            row[f"third_3_{name}"] = float(near(t, others, -NEAR, NEAR).mean()) if len(t) else np.nan
        rows.append(row)
    d = pd.DataFrame(rows)
    d["kind"] = np.where(d["excess_of_smaller_unit"] > 0.01, "excess",
                         np.where(d["near_over_expected"].between(0.85, 1.15), "about chance", "other"))
    d.to_csv(OUT / "third_units.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["unit_a", "unit_b", "tetrode", "kind", "near_over_expected", "n_coincident",
            "third_45_coincident", "third_45_flank", "third_45_lone", "third_3_coincident", "third_3_flank", "third_3_lone"]
    print(d.sort_values(["kind", "near_over_expected"], ascending=[True, False])[cols].round(3).to_string(index=False))
    print(d.groupby("kind")[[c for c in cols if c.startswith("third")]].median().round(3).to_string())


if __name__ == "__main__":
    main()
