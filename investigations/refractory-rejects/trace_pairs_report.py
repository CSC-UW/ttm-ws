"""Summaries and a figure from trace_pairs.py output (argv[1]); tables and the figure go to argv[2]."""
import glob
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = sys.argv[1]
R = sys.argv[2]
df = pd.read_csv(f"{S}/trace_pairs.csv.gz")
df["ms"] = df.d / 30.0
df["bin"] = df.kind
r = df.kind == "rej"
df.loc[r, "bin"] = pd.cut(df.ms[r], [0.29, 0.5, 0.7, 1.0, 1.5], labels=["rej 0.3-0.5", "rej 0.5-0.7", "rej 0.7-1.0",
                                                                        "rej 1.0-1.5"]).astype(str)
df.loc[df.kind == "far", "bin"] = "far 2.0-2.6"
df.loc[df.kind == "iso", "bin"] = "lone"
for p in ("", "syn_"):
    df[p + "gain2"] = df[p + "E1d"] - df[p + "E2"]
    df[p + "gain2x"] = df[p + "E1xd"] - df[p + "E2"]
    df[p + "one_wins"] = df[p + "E1xd"] < df[p + "E2"]
    df[p + "rel"] = df[p + "misfit_m"] / df[p + "model_energy"]
    df[p + "gain_d"] = df[p + "E1"] - df[p + "E1d"]
order = ["rej 0.3-0.5", "rej 0.5-0.7", "rej 0.7-1.0", "rej 1.0-1.5", "far 2.0-2.6", "lone"]
rows = []
for b in order:
    d = df[df.bin == b]
    pair = b != "lone"
    rows.append({
        "events": len(d), "units": d.unit.nunique(),
        "gain 2 copies vs 1 flexible: real": d.gain2.median() if pair else np.nan,
        "twin": d.syn_gain2.median() if pair else np.nan,
        "real/twin per event": (d.gain2 / d.syn_gain2).median() if pair else np.nan,
        "1 flexible (any unit) beats 2: real": d.one_wins.mean() if pair else np.nan,
        "twin ": d.syn_one_wins.mean() if pair else np.nan,
        "misfit / model energy: real": d.rel.median(), "twin  ": d.syn_rel.median(),
        "real - twin": (d.rel - d.syn_rel).median(),
        "derivative gain: real": d.gain_d.median(), "twin   ": d.syn_gain_d.median(),
        "background per sample: real": d.noise_out.median(), "twin    ": d.syn_noise_out.median(),
    })
tab = pd.DataFrame(rows, index=order).T
pd.set_option("display.width", 200)
print(tab.round(3).to_string())
tab.to_csv(f"{R}/trace_pairs_summary.csv")
df.to_csv(f"{S}/trace_pairs_scored.csv.gz", index=False)

per_unit = (df[df.bin.isin(order[:2])].assign(ratio=lambda x: x.gain2 / x.syn_gain2)
            .groupby("unit").agg(n=("ratio", "size"), ratio=("ratio", "median"), one_wins=("one_wins", "mean"),
                                 syn_one_wins=("syn_one_wins", "mean")))
per_unit = per_unit[per_unit.n >= 30]
per_unit.to_csv(f"{R}/trace_pairs_per_unit.csv")
print(f"\nrejected pairs 0.3-0.7 ms, units with 30+ events ({len(per_unit)}): real/twin gain median per unit")
print(per_unit.ratio.describe(percentiles=[.1, .25, .5, .75, .9]).round(3).to_string())
print("units where one flexible spike beats two copies in more than twice the twins' share:",
      int((per_unit.one_wins > 2 * per_unit.syn_one_wins + 0.02).sum()))

tr = np.concatenate([np.load(f, allow_pickle=True) for f in sorted(glob.glob(f"{S}/trace_pairs_traces_w*.npy"))])
meta = pd.DataFrame([{k: v for k, v in t.items() if not k.startswith("trace")} for t in tr])
top = meta[(meta.kind == "rej") & meta.d.between(9, 21)].unit.value_counts().index[:6]
t_ms = (np.arange(270) - 80) / 30.0
fig, axes = plt.subplots(2, 3, figsize=(12, 6.5), sharex=True)
blue, orange, gray = "#2a78d6", "#eb6834", "#7d7c78"
for ax, k in zip(axes.flat, top):
    m = (meta.unit == k) & (meta.kind == "rej") & meta.d.between(9, 21)
    dmode = int(meta.d[m].mode().iloc[0])
    sel = np.flatnonzero(m & meta.d.between(dmode - 1, dmode + 1))
    real = np.mean([tr[i]["trace"] for i in sel], axis=0)
    syn = np.mean([tr[i]["trace_syn"] for i in sel], axis=0)
    iso_i = np.flatnonzero((meta.unit == k) & (meta.kind == "iso"))
    iso = np.mean([tr[i]["trace"] for i in iso_i], axis=0)
    ax.axhline(0, color="#d6d5cf", lw=0.8)
    ax.plot(t_ms, iso, color=gray, lw=1.5, ls="--", label=f"lone spike (n={iso_i.size})")
    ax.plot(t_ms, syn, color=orange, lw=2, label="twin: the two fits + background")
    ax.plot(t_ms, real, color=blue, lw=2, label="data")
    ax.set_title(f"unit {k}: pairs {dmode / 30:.2f} ms apart (±1 sample), n={sel.size}", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
for ax in axes[-1]:
    ax.set_xlabel("ms from the first fit")
for ax in axes[:, 0]:
    ax.set_ylabel("µV, main wire")
axes.flat[0].legend(frameon=False, fontsize=8, loc="lower right")
fig.suptitle("Refractory-rejected pairs on TTM-003-48H: averaged traces, other units subtracted", fontsize=11)
fig.tight_layout()
fig.savefig(f"{R}/trace_pairs.png", dpi=130)
print("figure:", f"{R}/trace_pairs.png")
