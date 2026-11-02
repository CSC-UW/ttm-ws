"""What a gain fitted on all samples of a quiet block does to the large common events.

TTM-003's 48 h session, windows 0 and 5, bandpassed without a reference. In each window the block of
20 s with the fewest event samples (``ttm.reference.event_mask``) is fitted twice with
``ttm.reference.estimate_wire_timing``: on all its samples, and on the event samples of the window's
other blocks (the estimator as used). Each fit is applied to the window's other blocks (shift each
wire, subtract its gain times the median of the shifted wires) and compared with the global median:
per wire, the RMS of what is left within 0.5 ms of the median's large excursions, and between them.

Run: .venv/bin/python investigations/common-reference/gain_choice.py
"""
from __future__ import annotations

import numpy as np
import spikeinterface as si
import spikeinterface.preprocessing as spre
from spikeinterface.preprocessing.phase_shift import apply_frequency_shift

from ttm.reference import event_mask, estimate_wire_timing

STORE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/2026-10-20_09-02-12.blosc-zstd.zarr"
FS = 30000
WIN = 1800 * FS
BLOCK = 20 * FS
EDGE = 3000


def apply(x, shift, gain):
    y = apply_frequency_shift(x, shift, axis=0)
    return y - gain[None, :] * np.median(y, axis=1, keepdims=True)


def main():
    rec = spre.bandpass_filter(si.read_zarr(STORE), freq_min=300.0, freq_max=6000.0, dtype="float32")
    names = [str(c) for c in rec.channel_ids]
    uv = float(rec.get_channel_gains()[0])
    for w in (0, 5):
        blocks = []
        for i in range(6):
            c = w * WIN + WIN * (i + 1) // 7
            blocks.append(rec.get_traces(start_frame=c - BLOCK // 2, end_frame=c + BLOCK // 2).astype(np.float64) * uv)
        ev = [event_mask(np.median(b, axis=1)) for b in blocks]
        q = int(np.argmin([e.sum() for e in ev]))
        tau_all, g_all, _ = estimate_wire_timing(blocks[q], event_mad=1e9)
        rest = [b for i, b in enumerate(blocks) if i != q]
        tau_ev, g_ev, n_ev = estimate_wire_timing(np.concatenate(rest), min_event_samples=0)
        print(f"== window {w}: quiet block {q} with {ev[q].sum()} event samples; other blocks {n_ev}")
        res = {"global median": [], "aligned, all-sample fit on the quiet block": [],
               "aligned, event-sample fit": []}
        for i, b in enumerate(blocks):
            if i == q:
                continue
            e = ev[i].copy()
            e[:EDGE] = e[-EDGE:] = False
            ne = ~ev[i]
            ne[:EDGE] = ne[-EDGE:] = False
            outs = {"global median": b - np.median(b, axis=1, keepdims=True),
                    "aligned, all-sample fit on the quiet block": apply(b, -tau_all, g_all),
                    "aligned, event-sample fit": apply(b, -tau_ev, g_ev)}
            for k, o in outs.items():
                res[k].append((np.sqrt((o[e] ** 2).mean(0)), np.sqrt((o[ne] ** 2).mean(0))))
        worst = np.argsort(-np.abs(g_all - g_ev))[:5]
        print("  wires with the largest gain difference: " + ", ".join(
            f"{names[c]} {g_all[c]:.2f} against {g_ev[c]:.2f}" for c in worst))
        for k, v in res.items():
            at_ev = np.mean([a for a, _ in v], axis=0)
            off = np.mean([b_ for _, b_ in v], axis=0)
            print(f"  {k:44s} RMS at events median {np.median(at_ev):6.1f} uV, worst wires "
                  + " ".join(f"{at_ev[c]:6.1f}" for c in worst) + f"; between events median {np.median(off):5.1f}")


if __name__ == "__main__":
    main()
