"""Does the hollow correlogram centre come from the cosine gate or from wobble itself?

For each wobble arm and each span window, the ungated matcher output saved by
``02_span_carry.py`` (``raw/w<i>.npz``) is scored twice for every same-tetrode
pair: all matched spikes, and only those passing the gate (cosine >= 0.60).
The bank is the same in both, so the difference is the gate alone. Reported
per arm: the median central (|lag| <= 1.5 ms) to flank (5 to 25 ms) ratio
over pairs with at least 30 flank coincidences, pooled over windows, and the
same for the circus re-run (no gate, so its gated column is the gate applied
after the fact). Writes ``results/gate_or_matcher.json`` and the per-pair
ratios ``results/gate_pairs.csv``. Each arm's merges (``merges.json``) are
then looked up in June's stage-1 pair ratios (``results/center_ratios.csv``)
and in the arm's own ungated ratios: ``results/merges_vs_ccg.csv``.

    uv run python investigations/wobble-fullsize-ttm001/08_gate_or_matcher.py
"""
import json
import pathlib

import numpy as np
import pandas as pd

from _span import FS, OUT, SPAN, log
from ttm import capture as cp

HERE = pathlib.Path(__file__).resolve().parent
MAX_LAG = int(25 * FS / 1000)


def pair_counts(s, ci, ids, groups):
    out = {}
    tr = {k: np.sort(s[ci == k]) for k in range(len(ids))}
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if groups[i] != groups[j] or tr[i].size == 0 or tr[j].size == 0:
                continue
            lags = np.abs(cp.ccg_lags(tr[i], tr[j], MAX_LAG)) / FS * 1000
            out[(int(ids[i]), int(ids[j]))] = np.array([(lags <= 1.5).sum(), ((lags >= 5) & (lags <= 25)).sum()])
    return out


def ratios(acc):
    r = [(c / 3.0) / (f / 40.0) for c, f in acc.values() if f >= 30]
    return np.array(r)


def pair_ratio(v):
    c, f = v
    return (c / 3.0) / (f / 40.0) if f >= 30 else np.nan


def main():
    res, rows = {}, []
    for arm in ["circus", "w35", "w25", "w35a"]:
        acc_all, acc_gate = {}, {}
        for wi in SPAN:
            z = np.load(OUT / f"span_{arm}" / "raw" / f"w{wi}.npz")
            s, ci, r = z["sample_index"], z["cluster_index"], z["r"]
            for acc, keep in ((acc_all, np.ones(s.size, bool)), (acc_gate, r >= 0.6)):
                for k, v in pair_counts(s[keep], ci[keep], z["ids"], z["groups"]).items():
                    acc[k] = acc.get(k, 0) + v
        for k in acc_all:
            rows.append({"arm": arm, "a": k[0], "b": k[1], "ratio_ungated": pair_ratio(acc_all[k]),
                         "ratio_gated": pair_ratio(acc_gate[k]) if k in acc_gate else np.nan})
        ra, rg = ratios(acc_all), ratios(acc_gate)
        res[arm] = {"pairs_ungated": int(ra.size), "median_ratio_ungated": round(float(np.median(ra)), 3),
                    "pairs_below_0.3_ungated": int((ra < 0.3).sum()), "pairs_gated": int(rg.size),
                    "median_ratio_gated": round(float(np.median(rg)), 3), "pairs_below_0.3_gated": int((rg < 0.3).sum())}
        log(f"{arm}: {res[arm]}")
    (HERE / "results" / "gate_or_matcher.json").write_text(json.dumps(res, indent=1))
    gp = pd.DataFrame(rows)
    gp.to_csv(HERE / "results" / "gate_pairs.csv", index=False)
    june = pd.read_csv(HERE / "results" / "center_ratios.csv")
    june = june[june.arm == "june"].set_index(["a", "b"]).ratio
    out = []
    for arm in ["june", "circus", "w35", "w25", "w35a"]:
        own = gp[gp.arm == arm].set_index(["a", "b"]) if arm != "june" else None
        for m in json.loads((OUT / f"down_{arm}" / "merges.json").read_text()):
            a, b = sorted((m["a"], m["b"]))
            out.append({"arm": arm, "a": a, "b": b, "verdict": m["verdict"], "cosine": m["cosine"],
                        "june_ratio": june.get((a, b), np.nan),
                        "own_ungated_ratio": own.ratio_ungated.get((a, b), np.nan) if own is not None else np.nan,
                        "own_gated_ratio": own.ratio_gated.get((a, b), np.nan) if own is not None else np.nan})
    mv = pd.DataFrame(out)
    mv.to_csv(HERE / "results" / "merges_vs_ccg.csv", index=False)
    for arm, d in mv.groupby("arm", sort=False):
        log(f"{arm}: {len(d)} merges; June ratio > 0.7 on {int((d.june_ratio > 0.7).sum())}, own ungated > 0.7 on "
            f"{int((d.own_ungated_ratio > 0.7).sum())}, both on {int(((d.june_ratio > 0.7) & (d.own_ungated_ratio > 0.7)).sum())}")


if __name__ == "__main__":
    main()
