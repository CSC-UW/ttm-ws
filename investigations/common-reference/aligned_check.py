"""What the aligned reference leaves at the few-unit tetrodes' units' spikes, against the global median.

Windows 3, 6 and 9 (the aligned binary covers windows 0 to 11): each unit's accepted spikes in the
per-tetrode carry over the 48 h (``R48_W_s45_c50_lnb_b2``, whose frames index both binaries), the
mean on its own tetrode's wires in ``mp_baseline/binary`` (global median) and in
``mp_aligned/binary_w00_11`` (aligned median): minimum and maximum (uV), the energy of the mean in
the aligned binary over that in the global one, and the cosine between the two means. The few-unit
tetrodes' units are listed first, then well-isolated units for comparison.

Run: .venv/bin/python investigations/common-reference/aligned_check.py
"""
from __future__ import annotations

import pathlib

import numpy as np
import spikeinterface as si

from ttm import capture as cp

SORTINGS = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix")
ARM = SORTINGS / "carry_bakeoff" / "real48" / "R48_W_s45_c50_lnb_b2"
FS = 30000
WIN = 1800 * FS
WINDOWS = (3, 6, 9)
FEW = (0, 2, 3, 39, 40, 43, 44)
ISOLATED = (4, 5, 16, 12, 26, 28)


def main():
    recs = {"global": si.load(SORTINGS / "mp_baseline" / "binary"),
            "aligned": si.load(SORTINGS / "mp_aligned" / "binary_w00_11")}
    gain = recs["global"].get_channel_gains()[0]
    acc = {}
    for w in WINDOWS:
        sp = np.load(ARM / f"spikes_w{w:02d}.npz")
        bank = np.load(ARM / f"bank_w{w:02d}.npz")
        grp = dict(zip(bank["ids"].tolist(), bank["group"].tolist()))
        a, b = w * WIN, (w + 1) * WIN
        for u in FEW + ISOLATED:
            if u not in grp:
                continue
            f = np.sort(sp["frame"][(sp["reason"] == 0) & (sp["unit"] == u)])
            f = f[(f > a + 30) & (f < b - 60)]
            if f.size > 5000:
                f = np.sort(np.random.default_rng(w).choice(f, 5000, replace=False))
            for name, rec in recs.items():
                x = cp.load_tetrode_traces(rec, grp[u], a, b)
                s, n = acc.get((u, name), (0.0, 0))
                acc[(u, name)] = (s + cp.snippets(x, f - a).sum(0), n + f.size)
            acc[(u, "group")] = grp[u]
    print(f"{'unit':10s} {'n':>6s} {'global min/max':>16s} {'aligned min/max':>16s} {'energy ratio':>12s} {'cosine':>7s}")
    for u in FEW + ISOLATED:
        if (u, "global") not in acc:
            continue
        (sg, n), (sa, _) = acc[(u, "global")], acc[(u, "aligned")]
        G, A = sg / n * gain, sa / n * gain
        cos = float((G * A).sum() / np.sqrt((G * G).sum() * (A * A).sum()))
        tag = "TT%d u%d" % (acc[(u, "group")] + 1, u)
        print(f"{tag:10s} {n:6d} {G.min():7.1f} {G.max():7.1f}  {A.min():7.1f} {A.max():7.1f}  "
              f"{(A * A).sum() / (G * G).sum():12.2f} {cos:7.3f}")


if __name__ == "__main__":
    main()
