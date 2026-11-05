"""How the sorted spikes of a tetrode fall around its unclaimed and reference events, from gate.py's
per-window files.

    .venv/bin/python investigations/multi-unit-gate/lags.py <gate out dir>

For each event, the sorted spikes of its tetrode at each lag from -3 to 3 ms (lag = spike minus event,
so positive lags are spikes after the event) and at 20 to 50 ms on either side. Within 0.5 ms a reference
event has only its claiming spike and an unclaimed event none, by definition, so those lags are left
empty and the claiming spike is not counted. Events are reference, or unclaimed and grouped by
where their best cosine falls among the tetrode's cutoffs losing 1, 5 and 10 % (results/tetrodes.csv).

Writes to results/:
lags.csv: per tetrode and kind, the events, the share with a sorted spike 0.5 to 1.5 ms away, and the
    spikes 0.5 to 1.5 ms before and after the event, each per ms over the spikes 20 to 50 ms away per ms.
lag_hist.csv: per tetrode, kind and lag in samples, the spikes at that lag per event, over the spikes per
    sample per event at 20 to 50 ms (empty within 0.5 ms).
"""

import pathlib
import sys

import numpy as np
import pandas as pd
from gate import CLAIM, Setup

OUT = pathlib.Path(__file__).parent / "results"
MAX_LAG, NEAR = 90, 45
BASE = (600, 1500)


def main():
    d = pathlib.Path(sys.argv[1])
    st = Setup()
    cuts = pd.read_csv(OUT / "tetrodes.csv").set_index("tetrode")[["cut_1", "cut_5", "cut_10"]]
    keys = ("frame", "group", "claims", "best_cos")
    ev = {k: [] for k in keys}
    for f in sorted(d.glob("w*.npz")):
        z = np.load(f)
        m = z["ev_best_uid"] >= 0
        for k in keys:
            ev[k].append(z[f"ev_{k}"][m])
    ev = {k: np.concatenate(v) for k, v in ev.items()}
    rows, hist_rows = [], []
    for g in sorted(np.unique(ev["group"]).tolist()):
        p = st.pool_f[g]
        on = ev["group"] == g
        un = on & (ev["claims"] == 0)
        c1, c5, c10 = cuts.loc[g + 1]
        c = ev["best_cos"]
        kinds = {"reference": on & (ev["claims"] == 1), "unclaimed below 1 %": un & (c < c1),
                 "unclaimed 1 to 5 %": un & (c >= c1) & (c < c5), "unclaimed 5 to 10 %": un & (c >= c5) & (c < c10),
                 "unclaimed at or above 10 %": un & (c >= c10)}
        for kind, m in kinds.items():
            fr = ev["frame"][m]
            if fr.size == 0:
                continue
            h = np.zeros(2 * MAX_LAG + 1, np.int64)
            for c0 in range(0, fr.size, 1_000_000):
                q = fr[c0:c0 + 1_000_000]
                i0 = np.searchsorted(p, q - MAX_LAG, "left")
                i1 = np.searchsorted(p, q + MAX_LAG, "right")
                cnt = i1 - i0
                e = np.repeat(np.arange(q.size), cnt)
                j = np.repeat(i0, cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
                lag = p[j] - q[e]
                h += np.bincount(lag + MAX_LAG, minlength=2 * MAX_LAG + 1)
            base = (np.searchsorted(p, fr + BASE[1], "right") - np.searchsorted(p, fr + BASE[0], "left")
                    + np.searchsorted(p, fr - BASE[0], "right") - np.searchsorted(p, fr - BASE[1], "left")).sum()
            per_sample = base / (2 * (BASE[1] - BASE[0] + 1)) / fr.size
            lags = np.arange(-MAX_LAG, MAX_LAG + 1)
            near_before = h[(lags >= -NEAR) & (lags <= -CLAIM - 1)].sum() / 30 / fr.size
            near_after = h[(lags >= CLAIM + 1) & (lags <= NEAR)].sum() / 30 / fr.size
            i0 = np.searchsorted(p, fr - NEAR, "left")
            i1 = np.searchsorted(p, fr + NEAR, "right")
            n_near = i1 - i0 - (1 if kind == "reference" else 0)
            rows.append({"tetrode": g + 1, "kind": kind, "events": int(fr.size),
                         "share_with_spike_0p5_1p5ms": float(np.mean(n_near > 0)),
                         "before_ratio": near_before / per_sample, "after_ratio": near_after / per_sample})
            hist_rows += [{"tetrode": g + 1, "kind": kind, "lag": int(lg),
                           "ratio": c / fr.size / per_sample if abs(lg) > CLAIM else np.nan}
                          for lg, c in zip(lags, h)]
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "lags.csv", index=False)
    pd.DataFrame(hist_rows).to_csv(OUT / "lag_hist.csv", index=False)
    pd.set_option("display.width", 200)
    print(t[t["tetrode"].between(5, 12)].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
