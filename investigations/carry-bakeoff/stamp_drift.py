"""Where each unit's spike times sit on its waveform, early and late in the 48 h.

For production (``mp_baseline_masked/assembled_reseed_c12``) and the two 48 h carries
(``carry_R48_pt``, ``carry_R48_joint``): per unit and window, the mean snippet of up to 300 of its
spikes on the plain binary (spike time at sample 30), and the sample of the trough on its deepest
wire. A carry whose templates wander inside their 90-sample window stamps its spikes away from the
trough by the same amount.

Writes ``results/stamp_drift.csv``.

Run: .venv/bin/python investigations/carry-bakeoff/stamp_drift.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from ttm import capture as cp  # noqa: E402

HERE = pathlib.Path(__file__).parent
SORTS = {"production": "mp_baseline_masked", "per_tetrode": "carry_R48_pt", "joint_whitened": "carry_R48_joint"}
WINDOWS = (0, 24, 48, 72, 94)


def main():
    rec = si.load(bk.BINARY)
    sorts = {n: si.load(bk.SORTINGS / p / "assembled_reseed_c12") for n, p in SORTS.items()}
    rng = np.random.default_rng(0)
    rows = []
    for w in WINDOWS:
        a, b = w * bk.WIN, (w + 1) * bk.WIN
        for g in range(16):
            x = cp.load_tetrode_traces(rec, g, a, b)
            for n, s in sorts.items():
                grp = np.asarray(s.get_property("group")).astype(int)
                for u in np.asarray(s.unit_ids)[grp == g]:
                    f = s.get_unit_spike_train(u, start_frame=a + bk.NB, end_frame=b - bk.NA)
                    if f.size < 50:
                        continue
                    f = np.sort(rng.choice(f, min(300, f.size), replace=False)) - a
                    T = cp.snippets(x, f).mean(axis=0)
                    c = int(T.min(axis=0).argmin())
                    rows.append({"sorting": n, "unit": int(u), "tetrode": f"TT{g + 1}", "window": w,
                                 "n": int(f.size), "trough_sample": int(T[:, c].argmin()),
                                 "trough": float(T[:, c].min())})
        print(f"window {w} done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "stamp_drift.csv", index=False)
    off = df.assign(off=(df.trough_sample - 30).abs())
    print(off.pivot_table(index="sorting", columns="window", values="off",
                          aggfunc=lambda v: f"{np.median(v):.0f} / {np.percentile(v, 90):.0f} / {(v >= 5).sum()}"
                          ).to_string())
    print("cells: median / 90th percentile of |trough - 30| samples / units at 5 or more")


if __name__ == "__main__":
    main()
