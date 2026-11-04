"""Same-tetrode unit pairs of the final TTM-003-48H sorting whose cross-correlogram has a refractory
dip, as `investigations/spike-sorting/97_production_pipeline.py`'s merge step judged them, with
their template cosine; random pairs from different tetrodes, which the per-tetrode matcher never
fits together, as a control.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/ccg_pairs.py

The CCG ratio and verdict follow `_assignment_eval.ccg_verdict_pair` (co-active 1800 s windows with
at least 5 spikes of each unit; coincidences within 1.5 ms per ms over those 5 to 25 ms away per ms;
below 0.3 "duplicate", above 0.7 "distinct", fewer than 30 flank coincidences or 2 co-active windows
"SEGREGATED"), with the lag counts vectorised. Template cosine is `ttm.tracking.cosine_from_templates`
with a 10-sample shift, on average templates of up to 5,000 random spikes per unit in µV.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "spike-sorting"))
from _assignment_eval import DIP, FILL, FLANK_MS, MIN_CO_WIN, MIN_FLANK, REFR_MS, verdict_of  # noqa: E402

from ttm.sorting import create_group_sparse_analyzer  # noqa: E402
from ttm.tracking import cosine_from_templates  # noqa: E402

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = pathlib.Path(__file__).parent / "results"
FS = 30000
WIN = 1800 * FS
N_CONTROL = 200


def count_within(tb, tu, lo, hi):
    """Number of (b, u) pairs with lo <= u - b <= hi, for sorted frame arrays."""
    return int((np.searchsorted(tu, tb + hi, "right") - np.searchsorted(tu, tb + lo, "left")).sum())


def ccg(ta, tb):
    nwin = int(max(ta[-1], tb[-1]) // WIN) + 1
    na, nb = np.bincount(ta // WIN, minlength=nwin), np.bincount(tb // WIN, minlength=nwin)
    co = np.flatnonzero((na >= 5) & (nb >= 5))
    ca, cb = ta[np.isin(ta // WIN, co)], tb[np.isin(tb // WIN, co)]
    r, f0, f1 = int(REFR_MS * FS / 1000), int(FLANK_MS[0] * FS / 1000), int(FLANK_MS[1] * FS / 1000)
    central = count_within(cb, ca, -r, r)
    flank = count_within(cb, ca, f0, f1) + count_within(cb, ca, -f1, -f0)
    ratio = (central / (2 * REFR_MS)) / (flank / (2 * (FLANK_MS[1] - FLANK_MS[0]))) if flank >= MIN_FLANK else np.nan
    return {"n_co_windows": int(co.size), "central": central, "flank": flank, "ratio": ratio,
            "verdict": verdict_of(ratio, co.size, dip=DIP, fill=FILL, min_co_win=MIN_CO_WIN)}


def main():
    OUT.mkdir(exist_ok=True)
    si.set_global_job_kwargs(n_jobs=48, chunk_duration="1s", progress_bar=False)
    rec = si.load(f"{RUN}/binary")
    sorting = si.load(f"{RUN}/sorting")
    az = create_group_sparse_analyzer(sorting, rec, format="memory", return_in_uV=True, n_jobs=48, seed=0)
    az.compute({"random_spikes": {"method": "uniform", "max_spikes_per_unit": 5000, "seed": 0},
                "templates": {"ms_before": 1.0, "ms_after": 2.0, "operators": ["average"]}})
    dense = np.asarray(az.get_extension("templates").get_data(operator="average"), dtype=np.float32)
    mask = az.sparsity.mask
    uids = list(sorting.unit_ids)
    group = np.asarray(sorting.get_property("group"))
    trains = [np.sort(sorting.get_unit_spike_train(u)).astype(np.int64) for u in uids]

    rows = []
    for i in range(len(uids)):
        for j in range(i + 1, len(uids)):
            same = group[i] == group[j]
            rows.append((i, j, same))
    same_pairs = [(i, j) for i, j, s in rows if s]
    rng = np.random.default_rng(0)
    cross = [(i, j) for i, j, s in rows if not s]
    control = [cross[k] for k in rng.choice(len(cross), size=N_CONTROL, replace=False)]

    out = []
    for kind, pairs in (("same tetrode", same_pairs), ("different tetrodes", control)):
        for i, j in pairs:
            row = {"pair": kind, "unit_a": uids[i], "unit_b": uids[j], "tetrode_a": int(group[i]) + 1,
                   "tetrode_b": int(group[j]) + 1, "spikes_a": len(trains[i]), "spikes_b": len(trains[j])}
            if kind == "same tetrode":
                ch = np.flatnonzero(mask[i])
                row["cosine"] = float(cosine_from_templates(dense[i][:, ch], dense[j][:, ch], max_shift_samples=10))
            row.update(ccg(trains[i], trains[j]))
            out.append(row)
    d = pd.DataFrame(out)
    d.to_csv(OUT / "ccg_pairs.csv", index=False)

    same = d[d["pair"] == "same tetrode"].copy()
    same["cosine_bin"] = pd.cut(same["cosine"], [-1, 0.5, 0.9, 0.95, 1.0], labels=["<0.5", "0.5-0.9", "0.9-0.95", ">=0.95"])
    pd.set_option("display.width", 250)
    print(pd.crosstab(same["cosine_bin"], same["verdict"], margins=True).to_string())
    ctrl = d[d["pair"] == "different tetrodes"]
    print("different tetrodes:", ctrl["verdict"].value_counts().to_dict(),
          "ratio median", round(float(ctrl["ratio"].median()), 3),
          "5th percentile", round(float(ctrl["ratio"].quantile(0.05)), 3))
    merge97 = same[(same["cosine"] >= 0.90) & ((same["verdict"] == "duplicate") | (same["cosine"] >= 0.95))]
    print(f"pairs the 97 merge rule would join: {len(merge97)}")
    print(same.sort_values("ratio")[["unit_a", "unit_b", "tetrode_a", "spikes_a", "spikes_b", "cosine",
                                     "ratio", "central", "flank", "n_co_windows", "verdict"]]
          .head(25).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
