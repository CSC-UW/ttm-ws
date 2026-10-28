"""Does a per-tetrode whitening matrix from the first 300 s of the 48 h binary still whiten later on?

For each tetrode and a 120 s stretch at the start of windows 0, 12, 24, 48, 72 and 95, the spike-excised
noise is whitened with the window-0 matrix; printed are the largest absolute correlation between wires
and the range of the wires' variances after whitening (1 and 0 when the matrix still fits).

Run: .venv/bin/python investigations/carry-bakeoff/whitening_drift.py
"""
from __future__ import annotations

import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402
from ttm import capture as cp  # noqa: E402

W = bk.whitening_matrix(hybrid=False)
rec = si.load(bk.BINARY)
groups = np.asarray(rec.get_property("group")).astype(int)
for g in range(16):
    ch = groups == g
    Wg = W[np.ix_(ch, ch)]
    out = []
    for w in (0, 12, 24, 48, 72, 95):
        a = w * bk.WIN
        x = cp.load_tetrode_traces(rec, g, a, a + int(120 * bk.FS)).astype(np.float64)
        noise, _ = cp.noise_mad_excised(x)
        keep = cp.excision_mask(x, noise, detect_mad=4.0, pad_frames=30)
        C = np.cov((x[keep] @ Wg).T)
        r = C / np.sqrt(np.outer(np.diag(C), np.diag(C)))
        out.append(f"w{w}: |r|max {np.abs(r[np.triu_indices(4, 1)]).max():.2f} var {np.diag(C).min():.2f}-{np.diag(C).max():.2f}")
    print(f"TT{g + 1:2d} " + " | ".join(out), flush=True)
