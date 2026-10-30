"""Each wire's time offset and gain against the median of all wires, across time and animals.

Per recording and time point, two blocks of 60 s bandpassed 300 to 6000 Hz without a reference. Per
wire, a least-squares fit over the block of the wire to the median of all 64 wires and its time
derivative (central difference), x = g * m + b * m'. A wire that records the common signal tau
samples after the median has x(t) = g * m(t - tau), close to g * m - g * tau * m', so tau = -b / g.
Only differences between wires mean anything: the median has no sampling time of its own.

Recordings: TTM-003's 48 h session (Zarr store) at 0, 12, 24, 36 and 47 h; TTM-001's NOD session
(Zarr store) at 1, 24 and 47 h; TTM-002's SD session (Open Ephys binary on /Volumes/tetrode) at 1,
24 and 47 h.

Writes ``results/sampling_offsets.csv`` (recording, hour, channel, tau in us, g) and prints, per
recording, the offset against channel number within CH1-CH32 and CH33-CH64, its change across time
points, and the correlation of the offset profiles between recordings.

Run: .venv/bin/python investigations/common-reference/sampling_offsets.py
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre

HERE = pathlib.Path(__file__).parent
FS = 30000.0
DATA = pathlib.Path("/nvme/neuropixels/tetrode_data")
RECORDINGS = {
    "TTM-003 48H": (lambda: si.read_zarr(str(DATA / "2026-10-20_09-02-12" / "2026-10-20_09-02-12.blosc-zstd.zarr")),
                    (0.1, 12, 24, 36, 47)),
    "TTM-001 NOD": (lambda: si.read_zarr(str(DATA / "2026-05-27_09-07-52" / "2026-05-27_09-07-52.blosc-zstd.zarr")),
                    (1, 24, 47)),
    "TTM-002 SD": (None, (1, 24, 47)),
}


def open_ttm2():
    folder = "/Volumes/tetrode/TTM2/2026-09-24_09-03-22"
    names, _ = se.get_neo_streams("openephysbinary", folder)
    name = [n for n in names if n.endswith("acquisition_board")][0]
    rec = se.read_openephys(folder, stream_name=name)
    return rec.select_channels([c for c in rec.channel_ids if str(c).startswith("CH")])


def fit(rec, start: int, n: int):
    gain = rec.get_channel_gains().astype(np.float64)
    x = rec.get_traces(start_frame=start, end_frame=start + n).astype(np.float64) * gain
    m = np.median(x, axis=1)
    A = np.stack([m, np.gradient(m)], axis=1)
    coef, *_ = np.linalg.lstsq(A, x, rcond=None)
    return coef[0], -coef[1] / coef[0]


def main():
    rows = []
    for name, (opener, hours) in RECORDINGS.items():
        rec = open_ttm2() if name.startswith("TTM-002") else opener()
        bp = spre.bandpass_filter(rec, freq_min=300.0, freq_max=6000.0, dtype="float32")
        chno = np.array([int(str(c)[2:]) for c in bp.channel_ids])
        dur_h = bp.get_num_frames() / FS / 3600
        for h in hours:
            if h + 0.05 > dur_h:
                continue
            gs, ts = [], []
            for k in range(2):
                g, t = fit(bp, int((h * 3600 + 120 * k) * FS), int(60 * FS))
                gs.append(g)
                ts.append(t)
            g, t = np.mean(gs, 0), np.mean(ts, 0)
            for c in range(len(chno)):
                rows.append({"recording": name, "hour": h, "channel": int(chno[c]), "tau_us": t[c] / FS * 1e6,
                             "tau_us_block_diff": (ts[1][c] - ts[0][c]) / FS * 1e6, "g": g[c]})
            print(f"{name} {h} h done", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "sampling_offsets.csv", index=False)
    for name, d in df.groupby("recording", sort=False):
        p = d.pivot_table(index="channel", columns="hour", values="tau_us")
        mean = p.mean(axis=1)
        print(f"\n== {name}: offset range {mean.min():.1f} to {mean.max():.1f} us; "
              f"largest change of a wire across time points {(p.max(axis=1) - p.min(axis=1)).max():.1f} us; "
              f"median |difference| between a time point's two blocks {d.tau_us_block_diff.abs().median():.2f} us")
        for lo, hi in ((1, 32), (33, 64)):
            m = mean[(mean.index >= lo) & (mean.index <= hi)]
            sl, ic = np.polyfit(m.index, m.values, 1)
            res = m.values - (sl * m.index + ic)
            print(f"  CH{lo}-CH{hi}: slope {sl:.2f} us per channel, correlation {np.corrcoef(m.index, m.values)[0, 1]:.2f}, "
                  f"off the line by more than 5 us: " + ", ".join(f"CH{c} ({r:+.0f})" for c, r in zip(m.index, res) if abs(r) > 5))
        gg = d.groupby("channel").g.mean()
        print(f"  g: {gg.min():.2f} to {gg.max():.2f}; below 0.9: " + ", ".join(f"CH{c} {v:.2f}" for c, v in gg.items() if v < 0.9))
    prof = df.groupby(["recording", "channel"]).tau_us.mean().unstack("recording")
    print("\ncorrelation of the offset profiles between recordings")
    print(prof.corr().round(2).to_string())


if __name__ == "__main__":
    main()
