"""Coincidences within 3 samples (0.1 ms) between same-tetrode units of the final TTM-003-48H sorting,
against the rate 5 to 25 ms away, as a count of excess coincidences and a fraction of the smaller unit;
the 200 pairs from different tetrodes in ccg_pairs.csv, which the per-tetrode matcher fits separately,
as a control.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/zero_lag.py
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = pathlib.Path(__file__).parent / "results"
NEAR, F0, F1 = 3, 150, 750


def count(ta, tb, lo, hi):
    return int((np.searchsorted(tb, ta + hi, "right") - np.searchsorted(tb, ta + lo, "left")).sum())


def main():
    s = si.load(f"{RUN}/sorting")
    uids = list(s.unit_ids)
    group = np.asarray(s.get_property("group"))
    tr = [np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in uids]
    cos = pd.read_csv(OUT / "ccg_pairs.csv").set_index(["unit_a", "unit_b"])["cosine"]
    ctrl = pd.read_csv(OUT / "ccg_pairs.csv")
    ctrl = {(a, b) for a, b in zip(ctrl.loc[ctrl["pair"] == "different tetrodes", "unit_a"],
                                   ctrl.loc[ctrl["pair"] == "different tetrodes", "unit_b"])}
    rows = []
    for i in range(len(uids)):
        for j in range(i + 1, len(uids)):
            same = group[i] == group[j]
            if not same and (uids[i], uids[j]) not in ctrl:
                continue
            near = count(tr[i], tr[j], -NEAR, NEAR)
            flank = count(tr[i], tr[j], F0, F1) + count(tr[i], tr[j], -F1, -F0)
            expected = flank / (2 * (F1 - F0 + 1)) * (2 * NEAR + 1)
            small = min(len(tr[i]), len(tr[j]))
            rows.append({"pair": "same tetrode" if same else "different tetrodes",
                         "unit_a": uids[i], "unit_b": uids[j], "tetrode_a": int(group[i]) + 1,
                         "tetrode_b": int(group[j]) + 1,
                         "spikes_a": len(tr[i]), "spikes_b": len(tr[j]),
                         "cosine": cos.get((uids[i], uids[j]), np.nan),
                         "near": near, "expected": expected, "near_over_expected": near / expected if expected else np.nan,
                         "excess": near - expected, "excess_of_smaller_unit": (near - expected) / small})
    d = pd.DataFrame(rows).sort_values("excess_of_smaller_unit", ascending=False)
    d.to_csv(OUT / "zero_lag_pairs.csv", index=False)
    pd.set_option("display.width", 250)
    for kind, g in d.groupby("pair"):
        print(kind)
        print(g.head(12).round(4).to_string(index=False))
        print(f"{kind}: excess above 1 % of the smaller unit in {(g.excess_of_smaller_unit > 0.01).sum()} of {len(g)} "
              f"pairs; near/expected median {g.near_over_expected.median():.2f}")


if __name__ == "__main__":
    main()
