"""Summary tables for the report: scoreboards, coverage per tetrode, per-class unit aggregates, timing.

Reads what scripts 02, 03 and 05 wrote and the June products, and writes
``results/scoreboard.csv`` (the June three-axis scoreboard of each arm over
the span), ``results/coverage_by_tetrode.csv`` (detected events >= 10 MAD
claimed by any unit on the tetrode, single units alone and with the
multi-unit bucket), ``results/units_by_class.csv`` (per-unit medians and
counts by amplitude class for each arm), ``results/timing.json`` and
``results/summary.json``.

    uv run python investigations/wobble-fullsize-ttm001/06_summary.py
"""
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import spikeinterface as si

from _span import JUNE, OUT, SORTING_DIR, SPAN, WLEN, log

sys.path.insert(0, str(SORTING_DIR))
from _wobble_eval import coverage_by_band  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
ARMS = ["june", "w35", "w25", "w35a"]


def scoreboards():
    rows = []
    for arm in ARMS + ["circus"]:
        f = OUT / f"down_{arm}" / "production_scoreboard.json"
        if not f.exists():
            continue
        sb = json.loads(f.read_text())
        for var, res in sb.items():
            A, B, C = res.get("axis_A", {}), res.get("axis_B", {}), res.get("axis_C", {})
            rows.append({"arm": arm, "variant": var, "units": res["n_units"],
                         "cov_5.5-7": A.get("5.5-7"), "cov_7-9": A.get("7-9"), "cov_9-10": A.get("9-10"),
                         "cov_ge10": A.get(">=10_pooled"), "cov_ge12": A.get(">=12_pooled"), "cov_all": A.get("overall"),
                         "purity_full_sw": B.get("spikeweighted_purity"),
                         "purity_tight_sw": B.get("spikeweighted_purity_tight"),
                         "flagged": B.get("n_flagged"), "cross_contaminated": B.get("n_cross_contaminated"),
                         "oversplit": B.get("n_oversplit"), "ccg_dup_pairs": C.get("n_ccg_duplicate_pairs"),
                         "ccg_distinct_pairs": C.get("n_ccg_distinct_pairs"), "same_tet_pairs": C.get("n_same_tet_pairs")})
    return pd.DataFrame(rows)


def coverage_by_tetrode():
    z = np.load(JUNE / "spike_coverage.npz", mmap_mode="r")
    a0, b0 = SPAN[0] * WLEN, (SPAN[-1] + 1) * WLEN
    lo, hi = np.searchsorted(z["peak_sample"], [a0, b0])
    ps, pg, am = (np.asarray(z[k][lo:hi]) for k in ("peak_sample", "peak_group", "amp_mad"))
    rows = []
    for arm in ARMS:
        d = OUT / f"down_{arm}"
        if not (d / "assembled_prod").exists():
            continue
        for var in ("assembled_base", "assembled_prod_sua", "assembled_prod"):
            s = si.load(d / var)
            _, claimed = coverage_by_band(s, ps, pg, am)
            for g in range(16):
                m = pg == g
                big = m & (am >= 10)
                mid = m & (am >= 7) & (am < 10)
                rows.append({"arm": arm, "variant": var.replace("assembled_", ""), "tetrode": f"TT{g + 1}",
                             "events_ge10": int(big.sum()),
                             "cov_ge10": round(100 * claimed[big].mean(), 2) if big.any() else None,
                             "events_7_10": int(mid.sum()),
                             "cov_7_10": round(100 * claimed[mid].mean(), 2) if mid.any() else None})
    return pd.DataFrame(rows)


def by_class():
    rows = []
    for arm in ARMS:
        f = RESULTS / f"units_{arm}.csv"
        if not f.exists():
            continue
        t = pd.read_csv(f)
        for c, d in [("all", t), ("large", t[t.amp_class == "large"]), ("small", t[t.amp_class == "small"])]:
            rows.append({"arm": arm, "class": c, "units": len(d), "spikes": int(d.n_spikes.sum()),
                         "rate_hz_median": round(d.rate_hz.median(), 2),
                         "rp_contamination_median": round(d.rp_contamination.median(), 4),
                         "units_rp_ge_0.1": int((d.rp_contamination >= 0.1).sum()),
                         "purity_full_median": round(d.purity_full.median(), 4),
                         "purity_tight_median": round(d.purity_tight.median(), 4),
                         "purity_tight_spikeweighted": round(float((d.purity_tight * d.n_spikes).sum() / d.n_spikes[d.purity_tight.notna()].sum()), 4),
                         "units_purity_tight_lt_0.8": int((d.purity_tight < 0.8).sum()),
                         "ampsplit_ratio_median": round(d.ampsplit_ratio.median(), 3),
                         "units_ampsplit_ge_0.5": int((d.ampsplit_ratio >= 0.5).sum()),
                         "cos_median_median": round(d.cos_median.median(), 3),
                         "cos_p10_median": round(d.cos_p10.median(), 3),
                         "resid_median_median": round(d.resid_median.median(), 3),
                         "amp_cv_median": round(d.amp_cv.median(), 3),
                         "verdict_duplicate": int(d.n_duplicate.sum()), "verdict_ambiguous": int(d.n_ambiguous.sum()),
                         "verdict_distinct": int(d.n_distinct.sum()), "verdict_segregated": int(d.n_segregated.sum())})
    return pd.DataFrame(rows)


def timing():
    out = {}
    for arm in ["circus", "w35", "w25"]:
        f = OUT / f"span_{arm}" / "timing.json"
        if f.exists():
            t = json.loads(f.read_text())
            out[arm] = {k: t[k] for k in ("wall_s", "match_s", "fit_s", "reestimate_s", "pipeline_s", "n_spikes",
                                          "n_jobs", "started")}
        f = OUT / f"span_{arm}_timing" / "timing.json"
        if f.exists():
            t = json.loads(f.read_text())
            out[arm + "_alone_3w"] = {k: t[k] for k in ("windows", "wall_s", "match_s", "fit_s", "reestimate_s",
                                                        "pipeline_s", "n_jobs", "started")}
    for arm in ARMS + ["circus"]:
        f = OUT / f"down_{arm}" / "timing.json"
        if f.exists():
            out["downstream_" + arm] = json.loads(f.read_text())
    june_dirs = sorted((JUNE / "reseed_sorts_c12").iterdir())
    stamps = {p.name: p.stat().st_mtime for p in june_dirs}
    names = sorted(stamps, key=lambda n: int(n.removeprefix("reseed_w")))
    out["june_c12_segments_min"] = {f"{a}->{b}": round((stamps[b] - stamps[a]) / 60, 1) for a, b in zip(names, names[1:])}
    return out


def main():
    RESULTS.mkdir(exist_ok=True)
    sb = scoreboards()
    sb.to_csv(RESULTS / "scoreboard.csv", index=False)
    cov = coverage_by_tetrode()
    cov.to_csv(RESULTS / "coverage_by_tetrode.csv", index=False)
    bc = by_class()
    bc.to_csv(RESULTS / "units_by_class.csv", index=False)
    tm = timing()
    (RESULTS / "timing.json").write_text(json.dumps(tm, indent=1))
    extra = {}
    for arm in ARMS:
        for name in ("production_residual.json", "production_mua.json"):
            f = OUT / f"down_{arm}" / name
            if f.exists():
                extra.setdefault(arm, {})[name.removesuffix(".json")] = json.loads(f.read_text())
        f = OUT / f"down_{arm}" / "merges.json"
        if f.exists():
            m = json.loads(f.read_text())
            extra.setdefault(arm, {})["merges"] = {"n": len(m), "by_verdict": pd.Series([x["verdict"] for x in m]).value_counts().to_dict()}
    (RESULTS / "summary.json").write_text(json.dumps(extra, indent=1, default=int))
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(sb.to_string(index=False))
        print(bc.to_string(index=False))
        print(cov.groupby(["arm", "variant"]).apply(
            lambda d: round(100 * (d.cov_ge10 * d.events_ge10).sum() / 100 / d.events_ge10.sum(), 2)).to_string())
    print(json.dumps(tm, indent=1))
    log("done")


if __name__ == "__main__":
    main()
