"""Why the median reference leaves part of the common event on each wire.

Five blocks of 60 s of window 48 from the Zarr store, bandpassed 300 to 6000 Hz without a reference
(as ``reference_check.py``). Per wire, a least-squares fit of the wire to the median of all 64 wires
and its time derivative (central difference), x = g * m + b * m'. For a wire that receives the
common signal ``tau`` samples after the median, x(t) = g * m(t - tau), close to g * m - g * tau * m',
so tau = -b / g. Per wire: g and tau in each block, by the wire's channel number.

Then the mean at the few-unit tetrodes' units' spikes on their own tetrode after subtracting, from
each wire: the median of all 64 wires (the sorted binary); the median of the other tetrodes' 60
wires (cross-tetrode); g * m; g * m + b * m'; and, with spikeinterface's ``phase_shift`` moving each
wire by -tau first (and by +tau, the wrong way, as a control), the median of the shifted wires, and
each shifted wire's least-squares gain on that median times the median.

Run: .venv/bin/python investigations/carry-bakeoff/reference_skew.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import spikeinterface as si
import spikeinterface.preprocessing as spre

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from ttm import capture as cp  # noqa: E402

ARM = bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2"
STORE = bk.SESSION / f"{bk.SESSION.name}.blosc-zstd.zarr"
W = 48
UNITS = (0, 3, 39, 40, 43, 44, 2, 5, 16, 4)


def main():
    binr = si.load(bk.BINARY)
    bp = spre.bandpass_filter(si.read_zarr(str(STORE)), freq_min=300.0, freq_max=6000.0, dtype="float32")
    bp = bp.select_channels(list(binr.channel_ids))
    names = [str(c) for c in binr.channel_ids]
    chno = np.array([int(n[2:]) for n in names])
    groups = np.asarray(binr.get_property("group")).astype(int)
    gain = bp.get_channel_gains()
    a, L = W * bk.WIN, int(60 * bk.FS)
    starts = [a + int(s * bk.FS) for s in (60, 400, 750, 1100, 1450)]
    X = {s: bp.get_traces(start_frame=s, end_frame=s + L).astype(np.float64) * gain for s in starts}
    M = {s: np.median(X[s], axis=1) for s in starts}
    D = {s: np.gradient(M[s]) for s in starts}

    fits = []
    for s in starts:
        A = np.stack([M[s], D[s]], axis=1)
        coef, *_ = np.linalg.lstsq(A, X[s], rcond=None)
        fits.append(coef)
    fits = np.array(fits)
    g_blk = fits[:, 0]
    tau_blk = -fits[:, 1] / fits[:, 0]
    tau = tau_blk.mean(0)
    print(f"tau (samples, wire after the median): range {tau.min():.3f} to {tau.max():.3f}, "
          f"spread across blocks (median sd) {np.median(tau_blk.std(0)):.3f}")
    print(f"g: range {g_blk.mean(0).min():.3f} to {g_blk.mean(0).max():.3f}")
    o = np.argsort(chno)
    print("by channel number: CH tau_us g")
    for i in o:
        print(f"  CH{chno[i]:<3d} TT{groups[i] + 1:<3d} {tau[i] / bk.FS * 1e6:7.1f} us  g {g_blk[:, i].mean():.3f}")
    r = np.corrcoef(chno, tau)[0, 1]
    print(f"correlation of tau with channel number {r:.2f}")
    for lo, hi in ((1, 32), (33, 64)):
        m = (chno >= lo) & (chno <= hi)
        r_ = np.corrcoef(chno[m], tau[m])[0, 1]
        sl = np.polyfit(chno[m], tau[m] / bk.FS * 1e6, 1)[0]
        print(f"  CH{lo}-CH{hi}: correlation {r_:.2f}, slope {sl:.2f} us per channel")

    g_all, b_all = fits[:, 0].mean(0), fits[:, 1].mean(0)
    ref = {
        "median of 64": {s: np.repeat(M[s][:, None], 64, axis=1) for s in starts},
        "median of other 60": {},
        "g * median": {s: M[s][:, None] * g_all[None] for s in starts},
        "g * median + b * slope": {s: M[s][:, None] * g_all[None] + D[s][:, None] * b_all[None] for s in starts},
    }
    for s in starts:
        R = np.empty_like(X[s])
        for g in range(16):
            R[:, groups == g] = np.median(X[s][:, groups != g], axis=1, keepdims=True)
        ref["median of other 60"][s] = R
    for sign, tag in ((-1, "shifted"), (1, "shifted the wrong way")):
        ps = spre.phase_shift(bp, inter_sample_shift=sign * tau)
        Y = {s: ps.get_traces(start_frame=s, end_frame=s + L).astype(np.float64) * gain for s in starts}
        my = {s: np.median(Y[s], axis=1) for s in starts}
        ref[f"median, {tag}"] = {s: X[s] - Y[s] + my[s][:, None] for s in starts}
        if sign == -1:
            gy = sum((Y[s] * my[s][:, None]).sum(0) for s in starts) / sum(float((my[s] ** 2).sum()) for s in starts)
            ref["shifted, g * median"] = {s: X[s] - Y[s] + my[s][:, None] * gy[None] for s in starts}
    sp = dict(np.load(ARM / f"spikes_w{W:02d}.npz"))
    b = bk.Bank.load(ARM / f"bank_w{W:02d}.npz")
    grp = dict(zip(b.ids.tolist(), b.group.tolist()))
    acc = sp["reason"] == 0
    print("own tetrode, min and max of the mean at the unit's spikes, after subtracting:")
    print("unit      " + " | ".join(f"{k:>22s}" for k in ref))
    for u in UNITS:
        f = sp["frame"][acc & (sp["unit"] == u)]
        g = grp[u]
        cells = []
        for k, R in ref.items():
            T, n = 0.0, 0
            for s in starts:
                ff = f[(f > s + bk.NB) & (f < s + L - bk.NA)] - s
                if ff.size:
                    T, n = T + cp.snippets(X[s][:, groups == g] - R[s][:, groups == g], ff).sum(0), n + ff.size
            T = T / n
            cells.append(f"{T.min():8.1f} {T.max():8.1f}")
        print(f"u{u:<3d} TT{g + 1:<3d} " + " | ".join(f"{c:>22s}" for c in cells))


if __name__ == "__main__":
    main()
