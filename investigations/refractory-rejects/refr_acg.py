"""Autocorrelograms 0-5 ms of the fits that passed admission rules 1-4 (accepted plus refractory-rejected)
and of the accepted fits alone, per unit and window (every 8th window of TTM-003-48H), each normalised to
its own count per ms at lags 50-100 ms (1 = no dependence on the unit's previous spike)."""
import sys
import numpy as np
import pandas as pd
RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = sys.argv[1]
fs = 30000.0
edges = np.arange(0, 151, 6)  # 0.2 ms bins to 5 ms
B0, B1 = 1500, 3000

def acg(F):
    h = np.zeros(edges.size - 1)
    k = 1
    while k < F.size:
        d = F[k:] - F[:-k]
        if d.min() > edges[-1]:
            break
        h += np.histogram(d, edges)[0]
        k += 1
    base = (np.searchsorted(F, F + B1) - np.searchsorted(F, F + B0)).sum() / ((B1 - B0) / fs * 1e3)
    return h / (np.diff(edges) / fs * 1e3) / base

rows = []
for w in range(0, 96, 8):
    z = np.load(f"{RUN}/carry/spikes_w{w:02d}.npz")
    f, u, reason = z["frame"], z["unit"], z["reason"]
    for k in np.unique(u):
        m = (u == k) & np.isin(reason, (0, 3))
        acc = (u == k) & (reason == 0)
        if acc.sum() < 2000:
            continue
        rows.append({"w": w, "unit": int(k), "n": int(acc.sum()), "kind": "rules1to4", **dict(enumerate(acg(np.sort(f[m]))))})
        rows.append({"w": w, "unit": int(k), "n": int(acc.sum()), "kind": "accepted", **dict(enumerate(acg(np.sort(f[acc]))))})
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/refr_acg.csv", index=False)
cols = list(range(edges.size - 1))
print(f"{df.unit.nunique()} units, {len(df) // 2} unit-windows with 2,000+ accepted spikes")
print("lag (ms)   rules 1-4: median [25-75 %]        accepted: median [25-75 %]")
for c in cols:
    a = df[df.kind == "rules1to4"][c]
    b = df[df.kind == "accepted"][c]
    print(f"{edges[c] / 30:4.1f}-{edges[c + 1] / 30:3.1f}   {a.median():5.2f} [{a.quantile(.25):5.2f}-{a.quantile(.75):5.2f}]"
          f"            {b.median():5.2f} [{b.quantile(.25):5.2f}-{b.quantile(.75):5.2f}]")
