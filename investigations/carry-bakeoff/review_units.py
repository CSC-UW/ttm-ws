"""Follow-ups to the review of ``review/TT1_TT6_TT16_w48.zarr``, window 48 of the per-tetrode
whitened carry (``R48_W_s45_c50_lnb_b2``) on the plain binary.

- TT1 unit 0 split in two by its CH36 waveform (a two-component Gaussian mixture on the first three
  principal components of CH36): per part, the mean on TT1's wires, the width of the CH36 trough,
  the share of its spikes with an accepted spike on another tetrode within 0.2 ms (and chance, with
  the spikes shifted by 2 and 3 ms either way), and its mean trough on each tetrode's deepest wire
  in five blocks of 60 s.
- TT16 unit 44: per spike, the sample of the minimum and of the maximum on its largest wire, and
  how much each wanders across spikes.
- TT6: production's unit 16 against per tetrode, by which per-tetrode fit (any reason) each of its
  spikes lands on, with the mean on TT6's wires of the spikes per-tetrode unit 16 did not accept.
- The refractory rule: per unit, the lags of the removed fits to the nearest accepted fit of the same
  unit, and the share of removals at 1.25 to 1.5 ms (what a 1.25 ms rule would keep).

Run: .venv/bin/python investigations/carry-bakeoff/review_units.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import spikeinterface as si
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402
from ttm import capture as cp  # noqa: E402

ARM = bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2"
W = 48
TOL = 6


def chans(rec, g):
    idx = cp.tetrode_channel_indices(rec, g)
    return [str(c) for c in np.asarray(rec.channel_ids)[idx]], rec.get_channel_gains()[idx].astype(np.float32)


def footprint(rec, frames_abs, blocks, L):
    groups = np.asarray(rec.get_property("group")).astype(int)
    T, n = 0.0, 0
    for s, X in blocks.items():
        ff = frames_abs[(frames_abs > s + bk.NB) & (frames_abs < s + L - bk.NA)] - s
        if ff.size:
            T, n = T + cp.snippets(X, ff).sum(0), n + ff.size
    T = T / max(n, 1)
    return [float(T[:, groups == g].min()) for g in range(16)], n


def coincidence(f, other):
    obs = near(f, other, TOL).mean()
    chance = np.mean([near(f + int(d * 1e-3 * bk.FS), other, TOL).mean() for d in (2, 3, -2, -3)])
    return obs, chance


def width(w):
    """Full width at half the trough depth, in ms."""
    i = int(w.argmin())
    half = w[i] / 2
    lo = i
    while lo > 0 and w[lo] < half:
        lo -= 1
    hi = i
    while hi < w.size - 1 and w[hi] < half:
        hi += 1
    return (hi - lo) / bk.FS * 1e3


def main():
    rec = si.load(bk.BINARY)
    gain = rec.get_channel_gains()
    a, b = W * bk.WIN, (W + 1) * bk.WIN
    sp = dict(np.load(ARM / f"spikes_w{W:02d}.npz"))
    bank = bk.Bank.load(ARM / f"bank_w{W:02d}.npz")
    grp = dict(zip(bank.ids.tolist(), bank.group.tolist()))
    g_of = np.array([grp[int(u)] for u in sp["unit"]])
    acc = sp["reason"] == 0
    L = int(60 * bk.FS)
    starts = [a + int(s * bk.FS) for s in (60, 400, 750, 1100, 1450)]
    blocks = {s: rec.get_traces(start_frame=s, end_frame=s + L).astype(np.float32) * gain for s in starts}

    print("== TT1 unit 0 split by its CH36 waveform")
    g = 0
    names, gn = chans(rec, g)
    x = cp.load_tetrode_traces(rec, g, a, b) * gn
    f = np.sort(sp["frame"][acc & (sp["unit"] == 0)])
    f = f[(f - a >= bk.NB) & (f - a + bk.NA <= x.shape[0])]
    sn = cp.snippets(x, f - a)
    c36 = names.index("CH36")
    pcs = PCA(3, random_state=0).fit_transform(sn[:, :, c36])
    gm = GaussianMixture(2, random_state=0, n_init=5).fit(pcs)
    lab = gm.predict(pcs)
    order = np.argsort([sn[lab == k, :, c36].mean(0).min() for k in range(2)])
    other = np.sort(sp["frame"][acc & (g_of != g)])
    for k in order:
        m = lab == k
        T = sn[m].mean(0)
        obs, chance = coincidence(f[m], other)
        fp, n = footprint(rec, f[m], blocks, L)
        top = sorted(((v, h) for h, v in enumerate(fp) if h != g))[:4]
        print(f"part {k}: {m.sum():,} spikes ({m.mean():.0%}); trough per wire "
              + ", ".join(f"{c} {T[:, i].min():.0f}" for i, c in enumerate(names))
              + f" uV; CH36 width {width(T[:, c36]):.2f} ms")
        print(f"   with another tetrode within 0.2 ms {obs:.2f} (chance {chance:.2f}); own trough {fp[g]:.0f} uV, "
              f"largest elsewhere " + ", ".join(f"TT{h + 1} {v:.0f}" for v, h in top) + f" (n={n})")
        isi = np.diff(f[m]) / bk.FS * 1e3
        print(f"   ISIs under 3 ms {np.mean(isi < 3):.3f}, rate {m.sum() / (b - a) * bk.FS:.1f} Hz")

    print("\n== TT16 unit 44: where the minimum and maximum sit on its largest wire")
    g = 15
    names, gn = chans(rec, g)
    x = cp.load_tetrode_traces(rec, g, a, b) * gn
    f = np.sort(sp["frame"][acc & (sp["unit"] == 44)])
    f = f[(f - a >= bk.NB) & (f - a + bk.NA <= x.shape[0])]
    sn = cp.snippets(x, f - a).astype(np.float64)
    T = sn.mean(0)
    c = int(np.abs(T).max(0).argmax())
    print(f"template on {names[c]}: min {T[:, c].min():.1f} at sample {T[:, c].argmin()}, "
          f"max {T[:, c].max():.1f} at sample {T[:, c].argmax()} (spike time at {bk.NB})")
    seg = slice(bk.NB - 15, bk.NB + 30)
    imin = sn[:, seg, c].argmin(1) + seg.start
    imax = sn[:, seg, c].argmax(1) + seg.start
    for nm, i in (("minimum", imin), ("maximum", imax)):
        q = np.percentile(i, [10, 25, 50, 75, 90])
        print(f"per-spike {nm}: samples at 10/25/50/75/90 % " + " ".join(f"{v:.0f}" for v in q)
              + f"; interquartile {q[3] - q[1]:.0f} samples")

    print("\n== TT6: production unit 16 against per tetrode")
    g = 5
    names, gn = chans(rec, g)
    prod = si.load(bk.SORTINGS / "mp_baseline_masked" / "assembled_reseed_c12")
    fp_ = np.sort(prod.get_unit_spike_train(16, start_frame=a + bk.NB, end_frame=b - bk.NA))
    x = cp.load_tetrode_traces(rec, g, a, b) * gn
    mine = np.flatnonzero(g_of == g)
    fr_g = sp["frame"][mine]
    j = np.clip(np.searchsorted(fr_g, fp_), 1, fr_g.size - 1)
    pick = np.where(np.abs(fr_g[j - 1] - fp_) <= np.abs(fr_g[j] - fp_), j - 1, j)
    hit = np.abs(fr_g[pick] - fp_) <= TOL
    u_hit = np.where(hit, sp["unit"][mine][pick], -1)
    r_hit = np.where(hit, sp["reason"][mine][pick], -1)
    print(f"production unit 16: {fp_.size:,} spikes")
    for (u, r), n in sorted(((k, int(v)) for k, v in zip(*np.unique(np.c_[u_hit, r_hit], axis=0, return_counts=True))),
                           key=lambda kv: -kv[1])[:8]:
        sel = fp_[(u_hit == u) & (r_hit == r)] - a
        T = cp.snippets(x, sel).mean(0)
        what = "no per-tetrode fit" if u < 0 else f"pt unit {u}, reason {r}"
        print(f"  {what}: {n:,} ({n / fp_.size:.0%}); trough per wire "
              + ", ".join(f"{c} {T[:, i].min():.0f}" for i, c in enumerate(names)))

    print("\n== refractory removals: lag to the nearest accepted fit of the same unit")
    rf = sp["reason"] == 3
    rows = []
    for u in np.unique(sp["unit"][rf]):
        fa = np.sort(sp["frame"][acc & (sp["unit"] == u)])
        fr_ = sp["frame"][rf & (sp["unit"] == u)]
        if fa.size < 1000 or fr_.size < 50:
            continue
        j = np.clip(np.searchsorted(fa, fr_), 1, fa.size - 1)
        lag = np.minimum(np.abs(fa[j] - fr_), np.abs(fa[j - 1] - fr_)) / bk.FS * 1e3
        lag = lag[lag < 1.5]
        rows.append((u, grp[int(u)], fa.size, fr_.size / (fa.size + fr_.size), lag))
    for u, g, n_acc, share, lag in sorted(rows, key=lambda r: -r[3]):
        q = np.percentile(lag, [25, 50, 75])
        print(f"u{u:<4d} TT{g + 1:<3d} accepted {n_acc:6,d} removed share {share:.2f}; lag quartiles "
              + " ".join(f"{v:.2f}" for v in q) + f" ms; at 1.25-1.5 ms {np.mean(lag >= 1.25):.2f}; "
              f"under 0.5 ms {np.mean(lag < 0.5):.2f}")


if __name__ == "__main__":
    main()
