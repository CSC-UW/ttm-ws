"""Are the few-unit tetrodes' units local cells or events shared across tetrodes? And the mask's noise
under whitening.

Window 48 of the per-tetrode whitened carry (``R48_W_s45_c50_lnb_b2``):
- the mean of each unit's accepted spikes on all 64 channels of the plain binary (five blocks of
  60 s), as the trough on each tetrode's deepest wire, in uV;
- the share of its accepted spikes that have an accepted spike on another tetrode within 0.2 ms,
  against the same share with the spikes shifted by 2 and 3 ms either way (chance), and the
  tetrodes that share most.

Then, per tetrode, the variance after the first 300 s' whitening of noise that is independent
across wires at each wire's own level, which is what ``ttm.bad_periods.mask_bad_periods`` puts in
place of a bad period (real noise comes out at 1).

Run: .venv/bin/python investigations/carry-bakeoff/shared_events.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402
from ttm import capture as cp  # noqa: E402

ARM = bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2"
W = 48
UNITS = (0, 1, 3, 39, 40, 43, 44, 16, 26, 4, 12)


def footprints(sp, grp) -> None:
    rec = si.load(bk.BINARY)
    groups = np.asarray(rec.get_property("group")).astype(int)
    gain = rec.get_channel_gains()
    a, L = W * bk.WIN, int(60 * bk.FS)
    blocks = [a + int(s * bk.FS) for s in (60, 400, 750, 1100, 1450)]
    X = {s: rec.get_traces(start_frame=s, end_frame=s + L).astype(np.float32) * gain for s in blocks}
    print(f"mean trough per tetrode (deepest wire, uV) at the unit's accepted spike times, window {W}; own tetrode *")
    print("unit  own   " + " ".join(f"TT{g + 1:<3d}" for g in range(16)))
    acc = sp["reason"] == 0
    for u in UNITS:
        f = sp["frame"][acc & (sp["unit"] == u)]
        T, n = 0.0, 0
        for s in blocks:
            ff = f[(f > s + bk.NB) & (f < s + L - bk.NA)] - s
            if ff.size:
                T, n = T + cp.snippets(X[s], ff).sum(0), n + ff.size
        T = T / n
        print(f"u{u:<3d} TT{grp[u] + 1:<3d}"
              + " ".join(f"{T[:, groups == g].min():5.0f}{'*' if g == grp[u] else ' '}" for g in range(16))
              + f"   n={n}")
    noise = np.median(np.abs(np.concatenate([X[s][::50] for s in blocks])), axis=0) / 0.6745
    print("noise (uV, median of wires): " + " ".join(f"{np.median(noise[groups == g]):.1f}" for g in range(16)))


def coincidence(sp, grp) -> None:
    acc = sp["reason"] == 0
    g_of = np.array([grp[int(u)] for u in sp["unit"]])
    print("share of accepted spikes with an accepted spike on another tetrode within 0.2 ms (chance); top tetrodes")
    for u in UNITS:
        g = grp[u]
        f = np.sort(sp["frame"][acc & (sp["unit"] == u)])
        other = np.sort(sp["frame"][acc & (g_of != g)])
        chance = np.mean([near(f + int(d * 1e-3 * bk.FS), other, 6).mean() for d in (2, 3, -2, -3)])
        part = {h: near(f, np.sort(sp["frame"][acc & (g_of == h)]), 6).mean() for h in range(16) if h != g}
        top = sorted(part.items(), key=lambda kv: -kv[1])[:4]
        print(f"u{u:<3d} TT{g + 1:<3d} {near(f, other, 6).mean():.2f} (chance {chance:.2f}); "
              + ", ".join(f"TT{h + 1} {v:.2f}" for h, v in top))


def mask_noise() -> None:
    Wm = bk.whitening_matrix(False).astype(np.float64)
    print("variance after whitening of noise independent across wires at each wire's level (real noise: 1)")
    for g in range(16):
        s = slice(4 * g, 4 * g + 4)
        Wg = Wm[s, s]
        C = np.linalg.inv(Wg @ Wg)
        v = np.diag(Wg.T @ np.diag(np.diag(C)) @ Wg)
        r = (C / np.sqrt(np.outer(np.diag(C), np.diag(C))))[np.triu_indices(4, 1)]
        print(f"TT{g + 1:<2d} wire correlation {r.min():.2f} to {r.max():.2f}: {v.min():.2f} to {v.max():.2f}")


def main():
    sp = dict(np.load(ARM / f"spikes_w{W:02d}.npz"))
    b = bk.Bank.load(ARM / f"bank_w{W:02d}.npz")
    grp = dict(zip(b.ids.tolist(), b.group.tolist()))
    footprints(sp, grp)
    coincidence(sp, grp)
    mask_noise()


if __name__ == "__main__":
    main()
