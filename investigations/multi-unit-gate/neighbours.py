"""Whether the sorted spikes 0.5 to 1.5 ms after the unclaimed events that the subtracted score admits are
spikes at their own time, from gate.py's per-window files (run with the subtraction).

    .venv/bin/python investigations/multi-unit-gate/neighbours.py <gate out dir>

In every 8th window, for each unclaimed event of TT5 to TT12 at or above its tetrode's cutoff losing 5 % after
the subtraction (results/tetrodes_sub.csv), with a residual amplitude of 7 or more and a sorted spike 16 to 45
samples after it, and for up to 3,000 reference events per tetrode and window with such a spike: the nearest
such spike's trough (the most negative sample on any wire within 3 samples of it, over that wire's noise), as a
share of its unit's template trough measured the same way, and the cosine of its own 3 ms trace to its unit's
template. If a sorted spike were the event's own spike placed about 1 ms late, its trough at its own time would
be a small share of the template's.

Writes results/neighbours.csv: per tetrode, kind (admitted, reference) and whether the event's best unit after
the subtraction is the neighbour's unit, the events, the median trough share and cosine, the median lag, and
the share of the kind's events whose best unit is the neighbour's unit, observed and expected if the two were
unrelated (the sum over units of the unit's share of best units times its share of neighbours).
"""

import pathlib
import sys

import numpy as np
import pandas as pd
from gate import AFTER, BEFORE, CLAIM, FLOOR, RUN, Setup

RESULTS = pathlib.Path(__file__).parent / "results"
NEAR, HALF, N_REF, STEP = 45, 3, 3000, 8


def main():
    d = pathlib.Path(sys.argv[1])
    st = Setup()
    cut = pd.read_csv(RESULTS / "tetrodes_sub.csv").set_index("tetrode")["cut_5"]
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", np.int16, "r").reshape(-1, len(st.groups))
    off = np.arange(-BEFORE, AFTER)
    rows = []
    for w in range(0, 96, STEP):
        z = np.load(d / f"w{w:02d}.npz")
        rng = np.random.default_rng(w)
        tid = {int(u): j for j, u in enumerate(z["tmpl_ids"])}
        for g in range(4, 12):
            cols = np.flatnonzero(st.groups == g)
            nz = st.noise[cols]
            p, pu = st.pool_f[g], st.pool_u[g]
            m = (z["ev_group"] == g) & (z["ev_sub_uid"] >= 0)
            f, cl, su = z["ev_frame"][m], z["ev_claims"][m], z["ev_sub_uid"][m]
            j = np.searchsorted(p, f + CLAIM + 1, "left")
            has = (j < p.size) & (p[np.minimum(j, p.size - 1)] <= f + NEAR)
            admitted = (cl == 0) & (z["ev_sub_cos"][m] >= cut[g + 1]) & (z["ev_res_amp"][m] >= FLOOR) & has
            for kind, sel in (("admitted", admitted), ("reference", (cl == 1) & has)):
                ii = np.flatnonzero(sel)
                if kind == "reference" and ii.size > N_REF:
                    ii = rng.choice(ii, N_REF, replace=False)
                nf, nu = p[j[ii]], pu[j[ii]]
                ok = np.array([int(u) in tid for u in nu], bool)
                ii, nf, nu = ii[ok], nf[ok], nu[ok]
                if ii.size == 0:
                    continue
                x = raw[nf[:, None] + off][:, :, cols].astype(np.float32)
                T = z["tmpl_T"][[tid[int(u)] for u in nu]][:, :, :cols.size]
                cos = (x * T).sum((1, 2)) / np.sqrt((x ** 2).sum((1, 2)) * (T ** 2).sum((1, 2)))
                share = (-x[:, BEFORE - HALF:BEFORE + HALF + 1] / nz).max((1, 2)) \
                    / (-T[:, BEFORE - HALF:BEFORE + HALF + 1] / nz).max((1, 2))
                rows += [(g + 1, kind, int(b), int(n), c, a, lag)
                         for b, n, c, a, lag in zip(su[ii], nu, cos, share, nf - f[ii])]
    t = pd.DataFrame(rows, columns=["tetrode", "kind", "best_unit", "neighbour_unit", "cos", "trough_share", "lag"])
    t["same_unit"] = t["best_unit"] == t["neighbour_unit"]
    out = t.groupby(["tetrode", "kind", "same_unit"]).agg(
        events=("cos", "size"), trough_share_median=("trough_share", "median"), cos_median=("cos", "median"),
        lag_median=("lag", "median")).reset_index()
    for (tt, kind), q in t.groupby(["tetrode", "kind"]):
        pb, pn = q["best_unit"].value_counts(normalize=True), q["neighbour_unit"].value_counts(normalize=True)
        m = (out["tetrode"] == tt) & (out["kind"] == kind)
        out.loc[m, "same_unit_share"] = q["same_unit"].mean()
        out.loc[m, "same_unit_expected"] = float((pb * pn.reindex(pb.index, fill_value=0)).sum())
    out.to_csv(RESULTS / "neighbours.csv", index=False)
    pd.set_option("display.width", 200)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
