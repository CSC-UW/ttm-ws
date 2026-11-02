"""Whether the rate of large common events follows the animal's state, in TTM-003 and TTM-001.

Per 1800 s window, the events per second of ``common_events.py`` (``results/common_events.csv``)
against two state measures taken over the same six blocks of 20 s: the synthetic EMG (``ttm.emg``,
``per_window``: correlation of 300 to 600 Hz between single unreferenced wires of TT1, TT8 and
TT16) and slow-wave activity (``delta.idelta.zarr``, median over tetrodes). The EMG is not
independent of the events, since a signal common to all wires raises the correlation it measures;
slow-wave activity is in another band. Spearman correlation over windows, the event rate by
quartile of each measure, and inside against outside the sleep deprivation (TTM-003 windows 48 to
57, 09:00 to 14:00 on 2026-10-21 in the manifest). Writes ``results/event_state.csv``.

Run: .venv/bin/python investigations/common-reference/event_state.py
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr

HERE = pathlib.Path(__file__).parent
WIN = 1800.0
BLOCK = 20.0
N_BLOCKS = 6
EXP = {"TTM-003 48H": pathlib.Path("/nvme/neuropixels/tetrode_data/TTM-003/TTM-003-48H"),
       "TTM-001 NOD": pathlib.Path("/nvme/neuropixels/tetrode_data/TTM-001/TTM-NOD")}


def block_median(t, v, w):
    keep = np.zeros(t.size, bool)
    for i in range(N_BLOCKS):
        c = w * WIN + WIN * (i + 1) / (N_BLOCKS + 1)
        keep |= (t >= c - BLOCK / 2) & (t < c + BLOCK / 2)
    return float(np.nanmedian(v[keep])) if keep.any() else np.nan


def main():
    ev = pd.read_csv(HERE / "results" / "common_events.csv")
    rows = []
    for name, d in EXP.items():
        emg = xr.open_zarr(d / "synthetic_emg_methods.zarr")
        te, ve = emg.time.values, emg.per_window.values
        swa = xr.open_zarr(d / "delta.idelta.zarr").pwr
        ts = swa.time.values
        vs = swa.median("tetrode").values
        de = ev[ev.recording == name]
        for _, r in de.iterrows():
            w = int(r.window)
            rows.append({"recording": name, "window": w, "events_per_s": r.events_per_s,
                         "over_100uV_per_s": r.over_100uV_per_s, "emg": block_median(te, ve, w),
                         "swa": block_median(ts, vs, w)})
    df = pd.DataFrame(rows)
    df.round(5).to_csv(HERE / "results" / "event_state.csv", index=False)
    for name, d in df.groupby("recording", sort=False):
        print(f"== {name}")
        for m in ("emg", "swa"):
            r1 = spearmanr(d.events_per_s, d[m], nan_policy="omit")[0]
            r2 = spearmanr(d.over_100uV_per_s, d[m], nan_policy="omit")[0]
            q = pd.qcut(d[m], 4, labels=False)
            by_q = d.groupby(q).events_per_s.median().round(1).tolist()
            print(f"  {m}: Spearman with events/s {r1:+.2f}, with events over 100 uV/s {r2:+.2f}; "
                  f"events/s by quartile, low to high: {by_q}")
        sd = d.window.between(48, 57)
        print(f"  events/s in windows 48 to 57: {d[sd].events_per_s.median():.1f}, "
              f"elsewhere: {d[~sd].events_per_s.median():.1f}")
        off = d[~sd]
        print(f"  outside windows 48 to 57, Spearman with EMG {spearmanr(off.events_per_s, off.emg)[0]:+.2f}, "
              f"with SWA {spearmanr(off.events_per_s, off.swa)[0]:+.2f}")


if __name__ == "__main__":
    main()
