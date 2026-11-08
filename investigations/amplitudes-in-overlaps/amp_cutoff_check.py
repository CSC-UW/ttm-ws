"""SpikeInterface's amplitude_cutoff on the sorted units of both analyzers, as stored and recomputed on amplitudes
above 0, the histogram it works from, and how often spikes of each amplitude band have another unit's sorted
spike within 1 ms (30 samples) on the same tetrode.

Amplitudes are the analyzer's spike_amplitudes over all spikes (the stored metric uses the periods), with the
sign flipped when a unit's median is negative, as amplitude_cutoff does. Bands are amplitude over the unit's
median amplitude."""
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.metrics.quality.misc_metrics import amplitude_cutoff

OUT = sys.argv[1]
TD = "/nvme/neuropixels/tetrode_data"
S = {"TTM-003-48H": "2026-10-20_09-02-12", "TTM-003-48H-2": "2026-11-03_09-01-38"}
BANDS = [-np.inf, 0, 0.1, 0.25, 0.5, np.inf]
NAMES = ["<=0", "0-0.1", "0.1-0.25", "0.25-0.5", ">0.5"]
rows = []
for exp, ses in S.items():
    a = si.load_sorting_analyzer(f"{TD}/{ses}/sortings_carry/main/analyzer.zarr")
    s = a.sorting
    mua = np.asarray(s.get_property("is_mua"))
    grp = np.asarray(s.get_property("group")).astype(int)
    ids = s.unit_ids
    qm = a.get_extension("quality_metrics").get_data()
    amps = a.get_extension("spike_amplitudes").get_data(outputs="by_unit", concatenated=True)
    sv = s.to_spike_vector()
    for k in np.flatnonzero(~mua):
        u, g = ids[k], grp[k]
        fr = s.get_unit_spike_train(u)
        x = np.asarray(amps[u], float)
        if np.median(x) < 0:
            x = -x
        med = np.median(x)
        others = np.isin(sv["unit_index"], np.flatnonzero((grp == g) & ~mua & (np.arange(ids.size) != k)))
        of = np.sort(sv["sample_index"][others])
        nb = (np.searchsorted(of, fr + 30, "right") > np.searchsorted(of, fr - 30, "left")) if of.size else np.zeros(fr.size, bool)
        r = {"exp": exp, "unit": int(u), "tetrode": int(g) + 1, "spikes": int(x.size), "stored": float(qm.loc[u, "amplitude_cutoff"])}
        for name, y in (("all", x), ("pos", x[x > 0])):
            h = np.histogram(y, 100)[0]
            r[f"{name}_lowest_bin"] = int(h[0])
            r[f"{name}_min_rel"] = y.min() / med
            r[f"{name}_cutoff"] = amplitude_cutoff(y)
        r["max_rel"] = x.max() / med
        r["p999_rel"] = np.percentile(x, 99.9) / med
        band = np.digitize(x / med, BANDS[1:-1], right=True)
        for i, name in enumerate(NAMES):
            m = band == i
            r[f"n_{name}"] = int(m.sum())
            r[f"nbr_{name}"] = float(nb[m].mean()) if m.any() else np.nan
        r["nbr_all"] = float(nb.mean())
        rows.append(r)
    print(exp, "done", flush=True)
d = pd.DataFrame(rows)
d.to_csv(f"{OUT}/amp_cutoff_check.csv", index=False)
pd.set_option("display.width", 250)
d["set"] = np.where(d.tetrode == 15, "TT15", np.where(d.tetrode.isin([1, 2, 16]), "TT1/2/16", np.where(d.tetrode.between(4, 12), "TT4-12", "other")))
print(d.groupby("exp")[["stored", "pos_cutoff", "all_lowest_bin", "pos_lowest_bin", "max_rel", "p999_rel"]].median().round(5).to_string())
print(d.groupby("exp")[["stored", "pos_cutoff", "max_rel"]].max().round(4).to_string())
w = d[d.set == "TT4-12"]
for e, z in w.groupby("exp"):
    pooled = {n: (z[f"nbr_{n}"] * z[f"n_{n}"]).sum() / z[f"n_{n}"].sum() for n in NAMES}
    print(e, "TT4-12 share with another unit's spike within 1 ms, pooled:", {k: round(v, 3) for k, v in pooled.items()},
          "all:", round((z.nbr_all * z.spikes).sum() / z.spikes.sum(), 3),
          "counts:", {n: int(z[f"n_{n}"].sum()) for n in NAMES})
