"""Is the cosine gate hollowing the centre of every correlogram?

The gate scores each spike by the cosine of its raw 3 ms snippet to the
assigned template, so a spike with another spike inside that window scores
lower and is deleted. If that is what wobble's near-zero refractory
contamination and its extra CCG "duplicate" verdicts come from, then at stage
1 (same unit ids in every arm) the coincidences between same-tetrode units
within 1.5 ms, relative to their 5 to 25 ms flanks, should fall in the wobble
arms against June for pairs of units June calls distinct, and the
autocorrelograms should lose their 1 to 3 ms intervals too.

Per same-tetrode pair and arm: the central (|lag| <= 1.5 ms) to flank
(5 to 25 ms) ratio per ms, the June arbiter's quantity without its
co-activity restriction. Pooled over units: the interval histogram from 0 to
6 ms and the same-tetrode cross-unit lag histogram from 0 to 6 ms, each per
1000 spikes. Writes ``results/center_ratios.csv``,
``results/center_hollowing.json`` and ``center_hollowing.png``.

    uv run python investigations/wobble-fullsize-ttm001/07_center_hollowing.py
"""
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spikeinterface as si

from _span import FIG, FS, JUNE, OUT, SPAN, WLEN, log, restrict
from ttm import capture as cp

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
ARMS = ["june", "circus", "w35", "w25", "w35a"]
COLOURS = {"june": "k", "circus": "0.6", "w35": "tab:blue", "w25": "tab:green", "w35a": "tab:purple"}
EDGES_MS = np.arange(0, 6.25, 0.25)


def stage1(arm):
    a0, b0 = SPAN[0] * WLEN, (SPAN[-1] + 1) * WLEN
    if arm == "june":
        return restrict(si.load(JUNE / "assembled_reseed_c12"), a0, b0)
    return si.load(OUT / f"span_{arm}" / "sorting")


def main():
    RESULTS.mkdir(exist_ok=True)
    max_lag = int(25 * FS / 1000)
    rows, isi_h, xlag_h, nspk = [], {}, {}, {}
    for arm in ARMS:
        s = stage1(arm)
        groups = dict(zip([int(u) for u in s.unit_ids], np.asarray(s.get_property("group")).astype(int)))
        tr = {u: np.sort(s.get_unit_spike_train(u).astype(np.int64)) for u in groups}
        n = sum(t.size for t in tr.values())
        nspk[arm] = n
        isi = np.concatenate([np.diff(t) for t in tr.values() if t.size > 1]) / FS * 1000
        isi_h[arm] = np.histogram(isi, EDGES_MS)[0] / n * 1000
        xl = []
        for u in tr:
            for v in tr:
                if v <= u or groups[v] != groups[u] or not tr[u].size or not tr[v].size:
                    continue
                lags = np.abs(cp.ccg_lags(tr[u], tr[v], max_lag)) / FS * 1000
                central = int((lags <= 1.5).sum())
                flank = int(((lags >= 5) & (lags <= 25)).sum())
                ratio = (central / 3.0) / (flank / 40.0) if flank >= 30 else np.nan
                rows.append({"arm": arm, "a": u, "b": v, "tetrode": groups[u] + 1, "central": central, "flank": flank,
                             "ratio": ratio})
                xl.append(lags[lags < 6])
        xlag_h[arm] = np.histogram(np.concatenate(xl), EDGES_MS)[0] / n * 1000
        log(f"{arm}: {n:,} spikes")
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "center_ratios.csv", index=False)
    piv = df.pivot_table(index=["a", "b"], columns="arm", values="ratio")
    june_distinct = piv["june"] > 0.7
    out = {"pairs": int(len(piv)), "june_distinct_pairs": int(june_distinct.sum()), "per_arm": {}}
    for arm in ARMS:
        r = piv[arm]
        out["per_arm"][arm] = {
            "median_ratio_all": round(float(r.median()), 3),
            "median_ratio_on_june_distinct": round(float(r[june_distinct].median()), 3),
            "june_distinct_now_below_0.3": int((r[june_distinct] < 0.3).sum()),
            "pairs_below_0.3": int((r < 0.3).sum()), "pairs_above_0.7": int((r > 0.7).sum()),
            "isi_per_1000_0_1ms": round(float(isi_h[arm][:4].sum()), 3),
            "isi_per_1000_1_3ms": round(float(isi_h[arm][4:12].sum()), 3),
            "isi_per_1000_3_6ms": round(float(isi_h[arm][12:].sum()), 3),
            "xlag_per_1000_0_0.25ms": round(float(xlag_h[arm][0]), 3),
            "xlag_per_1000_0.25_2ms": round(float(xlag_h[arm][1:8].sum()), 3),
            "xlag_per_1000_2_3.75ms": round(float(xlag_h[arm][8:15].sum()), 3),
            "xlag_per_1000_0_1.5ms": round(float(xlag_h[arm][:6].sum()), 3),
            "xlag_per_1000_3_6ms": round(float(xlag_h[arm][12:].sum()), 3)}
    (RESULTS / "center_hollowing.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))

    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    c = EDGES_MS[:-1] + 0.125
    for arm in ARMS:
        axes[0].plot(c, isi_h[arm], color=COLOURS[arm], label=arm, drawstyle="steps-mid")
        axes[1].plot(c, xlag_h[arm], color=COLOURS[arm], label=arm, drawstyle="steps-mid")
    axes[0].set_title("intervals within a unit, pooled (stage 1, span 55-66)", fontsize=9)
    axes[1].set_title("lags between same-tetrode units, pooled", fontsize=9)
    for ax in axes[:2]:
        ax.set_xlabel("ms"); ax.set_ylabel("per 1000 spikes per 0.25 ms"); ax.legend(fontsize=8)
    for arm in ("circus", "w35", "w35a"):
        axes[2].scatter(piv["june"], piv[arm], s=8, color=COLOURS[arm], label=arm, alpha=0.7)
    lim = [0, max(2.0, float(np.nanpercentile(piv.values, 99)))]
    axes[2].plot(lim, lim, color="k", lw=0.6)
    axes[2].axvline(0.7, color="k", lw=0.5, ls=":"); axes[2].axhline(0.3, color="k", lw=0.5, ls=":")
    axes[2].set_xlim(lim); axes[2].set_ylim(lim)
    axes[2].set_xlabel("June: central/flank ratio of the pair"); axes[2].set_ylabel("arm: same pair")
    axes[2].set_title("same-tetrode pairs (0.3 duplicate, 0.7 distinct)", fontsize=9); axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "center_hollowing.png", dpi=100)
    plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
