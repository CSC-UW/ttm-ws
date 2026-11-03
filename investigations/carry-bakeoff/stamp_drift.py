"""Where each unit's spike times sit on its waveform, early and late in the 48 h.

For production (``mp_baseline_masked/assembled_reseed_c12``) and the two 48 h carries
(``carry_R48_pt``, ``carry_R48_joint``): per unit and window, the mean snippet of up to 300 of its
spikes on the plain binary (spike time at sample 30), the sample of the trough on its deepest wire,
and on its largest wire (largest absolute value) the samples of the minimum and the maximum. A carry
whose templates wander inside their 90-sample window stamps its spikes away from the trough by the
same amount. For a unit whose maximum is more than twice the depth of its minimum (positive-going)
the trough is a small lobe. The summary also gives, per unit and window, the lag that best aligns
the mean snippet (all four wires) with the unit's mean snippet in its first window in the table,
which follows the waveform whatever its sign. The mean snippets are written to
``<BAKE_ROOT>/review/stamp_drift_means.npz``.

Sortings that are not on disk are skipped. With ``BAKE_REFERENCE=aligned48`` the binary is the aligned
48 h binary and the sortings are the carries made on it (``carry_A48_pt``, ``carry_A48_joint_rules``).

Writes ``results/stamp_drift.csv`` (``results/stamp_drift_aligned48.csv``).

Run: [BAKE_REFERENCE=aligned48] .venv/bin/python investigations/carry-bakeoff/stamp_drift.py
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
SORTS = {"production": "mp_baseline_masked", "per_tetrode": "carry_R48_pt", "joint_whitened": "carry_R48_joint",
         "joint_rules": "carry_R48_joint_rules", "per_tetrode_rc": "carry_R48_pt_rc"}
SUFFIX = ""
if bk.REFERENCE == "aligned48":
    SORTS = {"aligned_per_tetrode": "carry_A48_pt", "aligned_joint_rules": "carry_A48_joint_rules"}
    SUFFIX = "_aligned48"
WINDOWS = (0, 24, 48, 72, 94)


def main():
    rec = si.load(bk.BINARY)
    sorts = {n: si.load(bk.SORTINGS / p / "assembled_reseed_c12") for n, p in SORTS.items()
             if (bk.SORTINGS / p / "assembled_reseed_c12").exists()}
    print("sortings: " + ", ".join(sorts))
    rng = np.random.default_rng(0)
    rows, means = [], {}
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
                    means[(n, int(u), w)] = T
                    c = int(T.min(axis=0).argmin())
                    m = int(np.abs(T).max(axis=0).argmax())
                    rows.append({"sorting": n, "unit": int(u), "tetrode": f"TT{g + 1}", "window": w,
                                 "n": int(f.size), "trough_sample": int(T[:, c].argmin()),
                                 "trough": float(T[:, c].min()), "wire": m,
                                 "min_sample": int(T[:, m].argmin()), "min": float(T[:, m].min()),
                                 "max_sample": int(T[:, m].argmax()), "max": float(T[:, m].max())})
        print(f"window {w} done", flush=True)
    df = pd.DataFrame(rows)
    df["align"] = [align_shift(means, r.sorting, r.unit, r.window) for r in df.itertuples()]
    df.to_csv(HERE / "results" / f"stamp_drift{SUFFIX}.csv", index=False)
    out = bk.BAKE_ROOT / "review" / f"stamp_drift_means{SUFFIX}.npz"
    np.savez(out, keys=np.array([f"{n}|{u}|{w}" for n, u, w in means]), T=np.stack(list(means.values())))
    summarize(df)


def align_shift(means: dict, n: str, u: int, w: int, max_lag: int = 45) -> int:
    """Samples by which the unit's mean snippet in window ``w`` is late against its first window's."""
    w0 = min(k[2] for k in means if k[0] == n and k[1] == u)
    A, B = means[(n, u, w0)].astype(np.float64), means[(n, u, w)].astype(np.float64)
    lags = np.arange(-max_lag, max_lag + 1)
    L = A.shape[0]
    score = [(A[max(0, -d):L - max(0, d)] * B[max(0, d):L - max(0, -d)]).sum() for d in lags]
    return int(lags[int(np.argmax(score))])


def summarize(df: pd.DataFrame) -> None:
    pos = df["max"] > 2 * df["min"].abs()
    df = df.assign(ref=np.where(pos, df.max_sample, df.min_sample), positive=pos)
    first = df.sort_values("window").groupby(["sorting", "unit"]).first()
    df = df.join(first[["ref", "wire"]].rename(columns={"ref": "ref0", "wire": "wire0"}), on=["sorting", "unit"])
    df["off"] = (df.ref - df.ref0).abs()
    df["off_trough"] = (df.trough_sample - 30).abs()
    agg = lambda v: f"{np.median(v):.0f} / {np.percentile(v, 90):.0f} / {(v >= 5).sum()} / {v.max():.0f}"
    df["off_align"] = df["align"].abs()
    for col, what in (("off_trough", "|trough - 30|"), ("off", "shift of the larger-sign extremum from the first window"),
                      ("off_align", "|lag that aligns the mean snippet with the first window's|")):
        print(f"{what}: median / 90th percentile / units at 5 or more / largest, in samples")
        print(df.pivot_table(index="sorting", columns="window", values=col, aggfunc=agg).to_string())
    last = df[df.window == df.window.max()]
    print("positive-going units in the last window: "
          + ", ".join(f"{r.sorting} {r.tetrode} u{r.unit}" for r in last[last.positive].itertuples()))
    print("wire changed since the first window, last window: "
          + ", ".join(f"{r.sorting} {r.tetrode} u{r.unit}" for r in last[last.wire != last.wire0].itertuples()))
    big = last[last.off_align >= 5].sort_values("off_align", ascending=False)
    print(big[["sorting", "tetrode", "unit", "trough_sample", "min_sample", "max_sample", "min", "max", "off", "align"]]
          .round(0).to_string(index=False))


if __name__ == "__main__":
    main()
