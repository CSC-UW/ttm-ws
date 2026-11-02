"""How many large common events the median of all wires carries, and when, in three recordings.

The median of the 64 wires, bandpassed 300 to 6000 Hz without a reference, in six blocks of 20 s
spread over each 1800 s window (the blocks ``ttm.reference.estimate_wire_timing_window`` reads).
An event is a local extreme of the median beyond 5 of its MAD in the block (either sign), at least
1 ms from a larger one; also counted at fixed sizes (beyond 50 and 100 uV). Per window: events per
second, the share of samples within 0.5 ms of one (``ttm.reference.event_mask``, what the timing
fit uses), the median trace's MAD and the 99th percentile of the events' size.

Recordings: TTM-003's 48 h session (96 windows), TTM-001's NOD session (Zarr store, 96 windows),
TTM-002's SD session (Open Ephys binary, every fourth window). Writes ``results/common_events.csv``
and prints per recording the distribution over windows.

Run: .venv/bin/python investigations/common-reference/common_events.py
"""
from __future__ import annotations

import pathlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre

from ttm.reference import event_mask

HERE = pathlib.Path(__file__).parent
FS = 30000
WIN = 1800 * FS
BLOCK = 20 * FS
N_BLOCKS = 6
DATA = pathlib.Path("/nvme/neuropixels/tetrode_data")


def open_rec(name):
    if name == "TTM-003 48H":
        raw = si.read_zarr(str(DATA / "2026-10-20_09-02-12" / "2026-10-20_09-02-12.blosc-zstd.zarr"))
    elif name == "TTM-001 NOD":
        raw = si.read_zarr(str(DATA / "2026-05-27_09-07-52" / "2026-05-27_09-07-52.blosc-zstd.zarr"))
    else:
        folder = "/Volumes/tetrode/TTM2/2026-09-24_09-03-22"
        names, _ = se.get_neo_streams("openephysbinary", folder)
        stream = [n for n in names if n.endswith("acquisition_board")][0]
        raw = se.read_openephys(folder, stream_name=stream)
        raw = raw.select_channels([c for c in raw.channel_ids if str(c).startswith("CH")][:64])
    return spre.bandpass_filter(raw, freq_min=300.0, freq_max=6000.0, dtype="float32")


def peaks(m, thr, dead=30):
    a = np.abs(m - np.median(m))
    cand = np.flatnonzero((a[1:-1] > thr) & (a[1:-1] >= a[:-2]) & (a[1:-1] >= a[2:])) + 1
    keep, taken = [], np.zeros(m.size, bool)
    for c in cand[np.argsort(-a[cand])]:
        if not taken[max(0, c - dead):c + dead + 1].any():
            keep.append(c)
            taken[c] = True
    return np.array(keep, np.int64), a


def one_window(args):
    name, w = args
    rec = open_rec(name)
    n = rec.get_num_frames()
    s, e = w * WIN, min((w + 1) * WIN, n)
    gain = float(rec.get_channel_gains()[0])
    counts = {"n5": 0, "n50": 0, "n100": 0, "event_frac": [], "mad": [], "sizes": []}
    for i in range(N_BLOCKS):
        c = s + (e - s) * (i + 1) // (N_BLOCKS + 1)
        x = rec.get_traces(start_frame=c - BLOCK // 2, end_frame=c + BLOCK // 2).astype(np.float64) * gain
        m = np.median(x, axis=1)
        mad = np.median(np.abs(m - np.median(m))) / 0.6745
        k, a = peaks(m, 5 * mad)
        counts["n5"] += k.size
        counts["n50"] += int((a[k] > 50).sum())
        counts["n100"] += int((a[k] > 100).sum())
        counts["event_frac"].append(event_mask(m).mean())
        counts["mad"].append(mad)
        counts["sizes"].extend(a[k].tolist())
    secs = N_BLOCKS * BLOCK / FS
    sz = np.array(counts["sizes"]) if counts["sizes"] else np.zeros(1)
    return {"recording": name, "window": w, "events_per_s": counts["n5"] / secs,
            "over_50uV_per_s": counts["n50"] / secs, "over_100uV_per_s": counts["n100"] / secs,
            "event_sample_frac": float(np.mean(counts["event_frac"])), "median_mad_uV": float(np.median(counts["mad"])),
            "size_p99_uV": float(np.percentile(sz, 99))}


def main():
    jobs = [("TTM-003 48H", w) for w in range(96)] + [("TTM-001 NOD", w) for w in range(96)] \
        + [("TTM-002 SD", w) for w in range(0, 96, 4)]
    with ProcessPoolExecutor(24) as ex:
        rows = list(ex.map(one_window, jobs))
    df = pd.DataFrame(rows).round(4)
    df.to_csv(HERE / "results" / "common_events.csv", index=False)
    pd.set_option("display.width", 200)
    for name, d in df.groupby("recording", sort=False):
        print(f"== {name}: {len(d)} windows")
        print(d[["events_per_s", "over_50uV_per_s", "over_100uV_per_s", "event_sample_frac", "median_mad_uV",
                 "size_p99_uV"]].describe(percentiles=[0.1, 0.5, 0.9]).round(2).to_string())
        print("  by quarter of the session, events per s (median over windows): "
              + ", ".join(f"{q}: {d[d.window // 24 == q].events_per_s.median():.1f}" for q in range(4)))
        print("  windows with fewer than 1 event per s: " + str(d[d.events_per_s < 1].window.tolist()))


if __name__ == "__main__":
    main()
