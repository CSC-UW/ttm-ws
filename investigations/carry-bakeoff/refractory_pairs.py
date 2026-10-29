"""Refractory pairs of the per-tetrode whitened carry: the procedure, or the units?

Per unit, the share of its fits that the 1.5 ms refractory rule removed (removed over accepted plus
removed):
- on the hybrid (``W_s45_c50_lnb_b2``, window 0), real units next to the ground-truth units put on
  the same tetrodes, which have no intervals under 1.5 ms;
- on the real 48 h (``R48_W_s45_c50_lnb_b2``, windows 0 and 48), with the lag from each removed fit
  to the nearest accepted spike of its unit, in 0.1 ms bins (double fits of one event would sit
  under 0.3 ms).

``--shapes`` adds, for the units with most removals, the cosine of the mean waveform of their
accepted spikes in four snr bins to that of their spikes of snr 12 or more, and the trough of
each wire over the deepest (plain traces).

Run: .venv/bin/python investigations/carry-bakeoff/refractory_pairs.py [--shapes]
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

REAL = bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2"
UNITS = (0, 1, 3, 39, 40, 42, 43, 44, 16, 26)
SNR_BINS = (4.5, 6, 8, 12, np.inf)


def removed_share(sp, u) -> tuple[float, int]:
    m = sp["unit"] == u
    n_acc, n_rem = int((m & (sp["reason"] == 0)).sum()), int((m & (sp["reason"] == 3)).sum())
    return n_rem / max(n_acc + n_rem, 1), n_acc


def hybrid() -> None:
    b = bk.Bank.load(bk.BAKE / "bank_start.npz")
    sp = np.load(bk.BAKE / "arms" / "W_s45_c50_lnb_b2" / "spikes_w00.npz")
    print("hybrid, window 0: share removed by the refractory rule (accepted spikes)")
    for g in np.unique(b.group):
        out = []
        for u in b.ids[b.group == g]:
            share, n = removed_share(sp, u)
            if n >= 200:
                out.append(f"{f'GT{u - 1000}' if u >= 1000 else f'u{u}'} {share:.2f} ({n})")
        print(f"  TT{g + 1}: " + ", ".join(out))


def lags() -> None:
    edges = np.arange(0, 46, 3)
    print("real 48 h: share removed; lag of removed fits to the nearest accepted spike, 0.1 ms bins to 1.5 ms")
    for w in (0, 48):
        sp = np.load(REAL / f"spikes_w{w:02d}.npz")
        for u in UNITS:
            m = sp["unit"] == u
            acc = np.sort(sp["frame"][m & (sp["reason"] == 0)])
            rem = sp["frame"][m & (sp["reason"] == 3)]
            j = np.clip(np.searchsorted(acc, rem), 1, acc.size - 1)
            d = np.minimum(np.abs(acc[j] - rem), np.abs(rem - acc[j - 1]))
            h = np.histogram(d, edges)[0]
            print(f"  w{w} u{u:<3d} removed {removed_share(sp, u)[0]:.2f} of {acc.size + rem.size:>7,}:",
                  " ".join(f"{v:.2f}" for v in h / max(h.sum(), 1)))


def shapes() -> None:
    import spikeinterface as si

    from ttm import capture as cp

    rec = si.load(bk.BINARY)
    b0 = bk.Bank.load(bk.BAKE_ROOT / "bank0.npz")
    grp = dict(zip(b0.ids.tolist(), b0.group.tolist()))
    rng = np.random.default_rng(0)
    rows = []
    for w in (0, 48):
        sp = np.load(REAL / f"spikes_w{w:02d}.npz")
        a = w * bk.WIN
        xs = {}
        for u in UNITS:
            g = grp[u]
            x = xs.setdefault(g, cp.load_tetrode_traces(rec, g, a, a + bk.WIN))
            acc = (sp["unit"] == u) & (sp["reason"] == 0)
            f, snr = sp["frame"][acc] - a, sp["snr"][acc]
            ok = (f > bk.NB) & (f < bk.WIN - bk.NA)
            f, snr = f[ok], snr[ok]
            top = f[snr >= 12] if (snr >= 12).sum() >= 30 else f[np.argsort(-snr)[:300]]
            T_top = cp.snippets(x, np.sort(top[:1000])).mean(0)
            for lo, hi in zip(SNR_BINS[:-1], SNR_BINS[1:]):
                fb = f[(snr >= lo) & (snr < hi)]
                if fb.size < 30:
                    continue
                T = cp.snippets(x, np.sort(rng.choice(fb, min(2000, fb.size), replace=False))).mean(0)
                pk = -T.min(0)
                rows.append({"window": w, "unit": u, "tetrode": f"TT{g + 1}", "snr": f"{lo:g} to {hi:g}",
                             "n": int(fb.size),
                             "cos_to_top": round(float((T * T_top).sum() / np.sqrt((T * T).sum() * (T_top * T_top).sum())), 3),
                             "wire_troughs": "/".join(f"{v:.2f}" for v in pk / pk.max())})
    pd.set_option("display.width", 200)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    hybrid()
    lags()
    if "--shapes" in sys.argv:
        shapes()
