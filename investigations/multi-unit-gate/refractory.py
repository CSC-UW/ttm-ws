"""Whether the events in each band of best cosine respect a unit's refractory period, from gate.py's
per-window files.

    .venv/bin/python investigations/multi-unit-gate/refractory.py <gate out dir>

For a reference event the unit is the one whose spike claims it; for an unclaimed event, the unit whose
template it matches best. The ratio is the unit's spikes 0.5 to 1.5 ms from the event (either side), per
ms, over its spikes 20 to 50 ms from the event, per ms. Within 0.5 ms a reference event has only its own
claiming spike and an unclaimed event none, by definition, so the window starts at 0.5 ms for both. A unit's
own spikes give a ratio near 0; events that fire independently of the unit give about 1.

Writes results/refractory_by_cosine.csv: per tetrode, kind and band of best cosine, the events, the
spikes counted in each window, the ratio, and for reference events the share whose claiming spike sits
within 2 samples of its usual place: the median, over the unit's reference events, of the claiming spike's
frame minus the event's.
"""

import pathlib
import sys

import numpy as np
import pandas as pd
from gate import CLAIM, Setup

OUT = pathlib.Path(__file__).parent / "results"
NEAR, BASE = (16, 45), (600, 1500)
ALIGN = 2
BANDS = ((-1.0, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5), (0.5, 0.6), (0.6, 0.7), (0.7, 1.01))


def count(train, frames, lo, hi):
    """Spikes of ``train`` at lags lo..hi samples (inclusive) on either side of each frame."""
    after = np.searchsorted(train, frames + hi, "right") - np.searchsorted(train, frames + lo, "left")
    before = np.searchsorted(train, frames - lo, "right") - np.searchsorted(train, frames - hi, "left")
    return after + before


def main():
    d = pathlib.Path(sys.argv[1])
    st = Setup()
    keys = ("frame", "group", "claims", "own", "best_uid", "best_cos")
    ev = {k: [] for k in keys}
    for f in sorted(d.glob("w*.npz")):
        z = np.load(f)
        m = z["ev_best_uid"] >= 0
        for k in keys:
            ev[k].append(z[f"ev_{k}"][m])
    ev = {k: np.concatenate(v) for k, v in ev.items()}
    unit = np.where(ev["claims"] == 1, ev["own"], ev["best_uid"])
    dt = np.zeros(unit.size, np.int64)
    for g, p in st.pool_f.items():
        i = np.flatnonzero((ev["group"] == g) & (ev["claims"] == 1))
        dt[i] = p[np.searchsorted(p, ev["frame"][i] - CLAIM, "left")] - ev["frame"][i]
    aligned = np.zeros(unit.size, bool)
    for u in np.unique(ev["own"][ev["claims"] == 1]).tolist():
        i = np.flatnonzero((ev["claims"] == 1) & (ev["own"] == u))
        aligned[i] = np.abs(dt[i] - np.median(dt[i])) <= ALIGN
    near = np.zeros(unit.size, np.int64)
    base = np.zeros(unit.size, np.int64)
    for u in np.unique(unit).tolist():
        i = np.flatnonzero(unit == u)
        near[i] = count(st.trains[u], ev["frame"][i], *NEAR)
        base[i] = count(st.trains[u], ev["frame"][i], *BASE)
    width = {w: 2 * (w[1] - w[0] + 1) / 30.0 for w in (NEAR, BASE)}
    rows = []
    for g in sorted(np.unique(ev["group"]).tolist()):
        for claims, kind in ((1, "reference"), (0, "unclaimed")):
            for lo, hi in BANDS:
                m = (ev["group"] == g) & (ev["claims"] == claims) & (ev["best_cos"] >= lo) & (ev["best_cos"] < hi)
                n_near, n_base = int(near[m].sum()), int(base[m].sum())
                rows.append({"tetrode": g + 1, "kind": kind, "cos_lo": lo, "cos_hi": min(hi, 1.0),
                             "events": int(m.sum()), "near": n_near, "base": n_base,
                             "ratio": (n_near / width[NEAR]) / (n_base / width[BASE]) if n_base else np.nan,
                             "aligned_within_2": float(aligned[m].mean()) if claims and m.any() else np.nan})
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "refractory_by_cosine.csv", index=False)
    live = t[t["tetrode"].between(5, 12)]
    print(live.pivot_table(index=["tetrode", "kind"], columns="cos_lo", values="ratio").round(2).to_string())
    print(live.pivot_table(index=["tetrode", "kind"], columns="cos_lo", values="events").to_string())
    print(live.pivot_table(index="tetrode", columns="cos_lo", values="aligned_within_2").round(3).to_string())


if __name__ == "__main__":
    main()
