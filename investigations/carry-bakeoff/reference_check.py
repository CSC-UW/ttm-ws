"""Do the few-unit tetrodes' events come through the global median reference?

The sorted binary (``mp_baseline/binary``) is the session's Zarr store bandpassed 300 to 6000 Hz
with the median of all 64 channels subtracted at each sample. Here, in five blocks of 60 s of
window 48, the same bandpass without the reference, and the median itself. At each unit's accepted
spikes in the per-tetrode whitened carry (``R48_W_s45_c50_lnb_b2``): the mean of the median trace
(minimum and maximum, uV), and per tetrode the minimum and maximum of the mean on its wires with
and without the reference. A deflection that the median carries and a tetrode's own wires do not
shows up on that tetrode, inverted, after referencing.

Then the troughs of the median trace below 5 of its MAD (dead time 1 ms): their rate, and per unit
the share of its accepted spikes within 1 ms of one. Last, each wire's least-squares gain on the
median, and the mean at each unit's spikes on its own tetrode with that gain times the median
subtracted instead of the median itself.

Run: .venv/bin/python investigations/carry-bakeoff/reference_check.py
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
UNITS = (0, 3, 39, 40, 43, 44, 2, 5, 21, 16, 26, 4, 12)


def main():
    binr = si.load(bk.BINARY)
    raw = si.read_zarr(str(STORE))
    bp = spre.bandpass_filter(raw, freq_min=300.0, freq_max=6000.0, dtype="float32")
    bp = bp.select_channels(list(binr.channel_ids))
    groups = np.asarray(binr.get_property("group")).astype(int)
    gain_b, gain_r = binr.get_channel_gains(), bp.get_channel_gains()
    a, L = W * bk.WIN, int(60 * bk.FS)
    starts = [a + int(s * bk.FS) for s in (60, 400, 750, 1100, 1450)]
    X, M, B = {}, {}, {}
    for s in starts:
        x = bp.get_traces(start_frame=s, end_frame=s + L).astype(np.float32) * gain_r
        X[s] = x
        M[s] = np.median(x, axis=1, keepdims=True)
        B[s] = binr.get_traces(start_frame=s, end_frame=s + L).astype(np.float32) * gain_b
    d = np.concatenate([(X[s] - M[s] - B[s])[::97] for s in starts])
    print(f"bandpass minus median against the binary: median |difference| {np.median(np.abs(d)):.2f} uV, "
          f"99th percentile {np.percentile(np.abs(d), 99):.2f} uV")
    med_noise = np.median(np.abs(np.concatenate([M[s][::10, 0] for s in starts]))) / 0.6745
    print(f"median trace noise {med_noise:.1f} uV")
    sp = dict(np.load(ARM / f"spikes_w{W:02d}.npz"))
    b = bk.Bank.load(ARM / f"bank_w{W:02d}.npz")
    grp = dict(zip(b.ids.tolist(), b.group.tolist()))
    acc = sp["reason"] == 0

    def mean_at(D, f):
        T, n = 0.0, 0
        for s in starts:
            ff = f[(f > s + bk.NB) & (f < s + L - bk.NA)] - s
            if ff.size:
                T, n = T + cp.snippets(D[s], ff, nbefore=bk.NB, nafter=bk.NA).sum(0), n + ff.size
        return T / max(n, 1), n

    print("per unit: median trace at its spikes (min, max uV); own tetrode with / without the reference "
          "(min, max of the mean on its wires); other tetrodes with the largest |min| without the reference")
    for u in UNITS:
        f = sp["frame"][acc & (sp["unit"] == u)]
        g = grp[u]
        tm, n = mean_at(M, f)
        tb, _ = mean_at(B, f)
        tx, _ = mean_at(X, f)
        own_b, own_x = tb[:, groups == g], tx[:, groups == g]
        others = sorted(((tx[:, groups == h].min(), h) for h in range(16) if h != g))[:3]
        print(f"u{u:<3d} TT{g + 1:<3d} n={n:5d}  median {tm.min():6.1f} {tm.max():6.1f} | own with "
              f"{own_b.min():6.1f} {own_b.max():6.1f}, without {own_x.min():6.1f} {own_x.max():6.1f} | "
              + ", ".join(f"TT{h + 1} {v:.0f}" for v, h in others))
    print("per unit, without the reference: the median over the 64 wires of each wire's mean trough "
          "(uV) and how many wires reach half the deepest")
    for u in UNITS:
        f = sp["frame"][acc & (sp["unit"] == u)]
        tx, _ = mean_at(X, f)
        tr = tx.min(axis=0)
        print(f"u{u:<3d} TT{grp[u] + 1:<3d} median wire trough {np.median(tr):6.1f}, deepest {tr.min():6.1f}, "
              f"wires at half the deepest or more {(tr <= tr.min() / 2).sum()}")
    ev = []
    for s in starts:
        m = M[s][:, 0]
        thr = -5 * med_noise
        cand = np.flatnonzero((m[1:-1] < thr) & (m[1:-1] <= m[:-2]) & (m[1:-1] <= m[2:])) + 1
        keep = []
        cand = cand[np.argsort(m[cand])]
        taken = np.zeros(m.size, bool)
        for c in cand:
            if not taken[max(0, c - 30):c + 31].any():
                keep.append(c)
                taken[c] = True
        ev.append(np.sort(np.array(keep, np.int64)) + s)
    ev = np.concatenate(ev)
    print(f"median-trace troughs below {5 * med_noise:.0f} uV: {ev.size} in {len(starts) * 60} s "
          f"({ev.size / (len(starts) * 60):.1f} per s)")
    from score import near
    for u in UNITS:
        f = np.sort(sp["frame"][acc & (sp["unit"] == u)])
        f = f[np.any([(f > s0 + bk.NB) & (f < s0 + L - bk.NA) for s0 in starts], axis=0)]
        print(f"u{u:<3d} TT{grp[u] + 1:<3d} within 1 ms of a median-trace trough {near(f, ev, 30).mean():.2f} (n={f.size})")
    # each wire's own gain on the median: least squares over the blocks
    num = sum((X[s] * M[s]).sum(0) for s in starts)
    den = sum(float((M[s] ** 2).sum()) for s in starts)
    beta = num / den
    print(f"each wire's gain on the median: {np.percentile(beta, 5):.3f} to {np.percentile(beta, 95):.3f} "
          f"(5th to 95th percentile); per tetrode " + " ".join(f"{np.median(beta[groups == g]):.2f}" for g in range(16)))
    G = {s: X[s] - M[s] * beta[None, :] for s in starts}
    print("own tetrode, min and max of the mean: with the median subtracted / with each wire's gain on it")
    for u in UNITS:
        f = sp["frame"][acc & (sp["unit"] == u)]
        g = grp[u]
        tb, _ = mean_at(B, f)
        tg, _ = mean_at(G, f)
        print(f"u{u:<3d} TT{g + 1:<3d} {tb[:, groups == g].min():6.1f} {tb[:, groups == g].max():6.1f} / "
              f"{tg[:, groups == g].min():6.1f} {tg[:, groups == g].max():6.1f}")


if __name__ == "__main__":
    main()
