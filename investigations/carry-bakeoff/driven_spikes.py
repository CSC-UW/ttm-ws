"""Where the TT11 driven cell's in-pulse spikes go under joint matching with the per-tetrode rules.

Unit 28 on TT11 (the OPTO2 session's unit 39) in the per-tetrode carry (``R48_W_s45_c50_lnb_b2``)
and the joint carry with the same rules and rebuild (``R48_W_joint_masked_s45_c50_lnb_b2``), in
the 24 windows with pulses. Two spikes are the same within 10 samples.

1. Each in-pulse unit-28 spike of the per-tetrode carry is classed by the joint carry: accepted on
   unit 28; fitted by unit 28 and rejected (by reason); accepted on another TT11 unit; not fitted.
   For the accepted, the not fitted and those on another unit: the share with an accepted spike on
   another tetrode within 3 ms whose snr exceeds the spike's per-tetrode snr. Shares by quarter of
   the session.
2. Windows 60, 72 and 84: the mean plain-trace snippet of each class against the OPTO2 templates
   of units 39 and 41 and the two units' means outside the pulses (more than 100 ms away), with the
   troughs on the four wires.
3. On ``review_driven.py``'s analyzer (window 84, up to 1000 spikes per group): per spike, the
   cosine to unit 28's and unit 29's means outside the pulses and the scale on each (projection
   over the template's energy), and for the in-pulse groups the scale on unit 28 by position in
   the pulse.

Writes ``results/driven_spikes.txt``.

Run: .venv/bin/python investigations/carry-bakeoff/driven_spikes.py
"""
from __future__ import annotations

import io
import pathlib
import sys
from contextlib import redirect_stdout

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402
from ttm import capture as cp  # noqa: E402

HERE = pathlib.Path(__file__).parent
PT = bk.BAKE_ROOT / "real48" / "R48_W_s45_c50_lnb_b2"
JR = bk.BAKE_ROOT / "real48" / "R48_W_joint_masked_s45_c50_lnb_b2"
OPTO2 = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-16_14-44-22/sortings_opto/capture/templates_TT11.npz")
U, U2, G = 28, 29, 10
REASONS = {1: "snr or floor", 2: "shape", 3: "refractory", 5: "bad period"}


def cos(a, b):
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def fates(on, off):
    counts = {}
    ratio = {"kept": [], "not fitted": [], "on another TT11 unit": []}
    for w in range(96):
        a, b = w * bk.WIN, (w + 1) * bk.WIN
        if not ((on >= a) & (on < b)).any():
            continue
        p, j = np.load(PT / f"spikes_w{w:02d}.npz"), np.load(JR / f"spikes_w{w:02d}.npz")
        bank = bk.Bank.load(JR / f"bank_w{w:02d}.npz")
        lut = np.full(int(bank.ids.max()) + 1, -1)
        lut[bank.ids] = bank.group
        m = (p["reason"] == 0) & (p["unit"] == U)
        fp, sp = p["frame"][m], p["snr"][m]
        k = bk.in_intervals(fp, on, off)
        o = np.argsort(fp[k])
        fp, sp = fp[k][o], sp[k][o]
        jf, ju, jr, js = j["frame"], j["unit"], j["reason"], j["snr"]
        jg, acc = lut[ju], jr == 0
        kept = near(fp, np.sort(jf[acc & (ju == U)]), 10)
        left = ~kept
        cls = np.full(fp.size, "", object)
        cls[kept] = "kept"
        for r, name in REASONS.items():
            f_ = np.sort(jf[(jr == r) & (ju == U)])
            hit = left & (near(fp, f_, 10) if f_.size else False)
            cls[hit] = f"fitted, rejected: {name}"
            left &= ~hit
        hit = left & near(fp, np.sort(jf[acc & (jg == G) & (ju != U)]), 10)
        cls[hit] = "on another TT11 unit"
        left &= ~hit
        cls[left] = "not fitted"
        mx = acc & (jg != G)
        o = np.argsort(jf[mx])
        xf, xs = jf[mx][o], js[mx][o]
        lo, hi = np.searchsorted(xf, fp - 90), np.searchsorted(xf, fp + 91)
        big = np.array([xs[l_:h].max() if h > l_ else 0.0 for l_, h in zip(lo, hi)])
        for name in ratio:
            ratio[name].append(np.c_[big[cls == name] / sp[cls == name], np.full((cls == name).sum(), w)])
        for c in np.unique(cls):
            counts[(w // 24, c)] = counts.get((w // 24, c), 0) + int((cls == c).sum())
    total = sum(counts.values())
    print(f"1. per-tetrode in-pulse unit {U} spikes: {total}")
    names = sorted({c for _, c in counts})
    for c in names:
        print(f"   {c:32s} {sum(v for (q, cc), v in counts.items() if cc == c):7d}")
    for name, r in ratio.items():
        r = np.concatenate(r)
        print(f"   {name:22s} with a larger other-tetrode spike within 3 ms {np.mean(r[:, 0] > 1):.2f}")
    print("   share by quarter of the session (windows):")
    for q in range(4):
        n = sum(v for (qq, _), v in counts.items() if qq == q) or 1
        print(f"   {q * 24:2d}-{q * 24 + 23:2d} (n={n}): " + ", ".join(
            f"{c} {counts.get((q, c), 0) / n:.2f}" for c in names))


def means(on, off):
    z = np.load(OPTO2)
    Tref, ids, wires = z["templates"], [int(u) for u in z["unit_ids"]], [str(w) for w in z["wires"]]
    rec = si.load(bk.BINARY)
    order = [wires.index(str(c)) for c in rec.channel_ids[cp.tetrode_channel_indices(rec, G)]]
    acc = {}
    for w in (60, 72, 84):
        a, b = w * bk.WIN, (w + 1) * bk.WIN

        def tr(path, u):
            z_ = np.load(path / f"spikes_w{w:02d}.npz")
            return np.sort(z_["frame"][(z_["reason"] == 0) & (z_["unit"] == u)])
        p28, p29, j28, j29 = tr(PT, U), tr(PT, U2), tr(JR, U), tr(JR, U2)
        inp = bk.in_intervals(p28, on, off)
        away = lambda f: ~bk.in_intervals(f, on - int(0.1 * bk.FS), off + int(0.1 * bk.FS))  # noqa: E731
        sets = {
            "in pulse, 28 in both": p28[inp & near(p28, j28, 10)],
            "in pulse, 28 pt / 29 joint": p28[inp & near(p28, j29, 10) & ~near(p28, j28, 10)],
            "in pulse, 28 pt / none joint": p28[inp & ~near(p28, j28, 10) & ~near(p28, j29, 10)],
            "outside, 28 in both": p28[away(p28) & near(p28, j28, 10)],
            "outside, 29 in both": p29[away(p29) & near(p29, j29, 10)],
        }
        x = cp.load_tetrode_traces(rec, G, a, b)
        for k, f in sets.items():
            f = f[(f > a + bk.NB) & (f < b - bk.NA)][:3000]
            if f.size:
                s, n = acc.get(k, (0, 0))
                acc[k] = (s + cp.snippets(x, f - a).sum(0), n + f.size)
    M = {k: s / n for k, (s, n) in acc.items()}
    refs = {"OPTO2 39": Tref[ids.index(39)][:, order], "OPTO2 41": Tref[ids.index(41)][:, order],
            "28 outside": M["outside, 28 in both"], "29 outside": M["outside, 29 in both"]}
    gain = rec.get_channel_gains()[0]
    print("\n2. windows 60, 72, 84: cosine of each class's mean; troughs (uV) on the four wires")
    print(f"   {'':30s} {'n':>6s} " + " ".join(f"{k:>10s}" for k in refs))
    for k, m in M.items():
        print(f"   {k:30s} {acc[k][1]:6d} " + " ".join(f"{cos(m, r):10.3f}" for r in refs.values())
              + "   " + " ".join(f"{v:6.1f}" for v in m.min(0) * gain))


def per_spike(on):
    an = si.load_sorting_analyzer(bk.BAKE_ROOT / "review" / "TT11_driven_w84.zarr")
    wf, sp = an.get_extension("waveforms"), an.get_extension("random_spikes").get_random_spikes()
    units = list(an.sorting.unit_ids)
    W = {u: wf.get_waveforms_one_unit(u) for u in units}
    T28, T29 = W["u28_out_both"].mean(0).ravel(), W["u29_out_both"].mean(0).ravel()

    def c_(X, t):
        X = X.reshape(len(X), -1)
        return X @ t / (np.linalg.norm(X, axis=1) * np.linalg.norm(t))

    def s_(X, t):
        return X.reshape(len(X), -1) @ t / (t @ t)
    print("\n3. window 84, per spike against the two units' means outside the pulses")
    print(f"   {'group':18s} {'n':>5s} {'closer to 28':>12s} {'scale on 28':>11s} {'scale on 29':>11s}")
    for u in units:
        X = W[u]
        print(f"   {u:18s} {len(X):5d} {np.mean(c_(X, T28) > c_(X, T29)):12.2f} "
              f"{np.median(s_(X, T28)):11.2f} {np.median(s_(X, T29)):11.2f}")
    print("   in-pulse groups: median scale on 28 by position in the pulse")
    a0 = 84 * bk.WIN
    for i, u in enumerate(units):
        if "out" in u:
            continue
        f = sp["sample_index"][sp["unit_index"] == i] + a0
        k = np.searchsorted(on, f, side="right") - 1
        _, off = bk.laser_pulses()
        pos = (f - on[k]) / (off[k] - on[k])
        sc = s_(W[u], T28)
        cells = [f"{lo:.2f}-{hi:.2f} {np.median(sc[(pos >= lo) & (pos < hi)]):.2f} (n {((pos >= lo) & (pos < hi)).sum()})"
                 for lo, hi in ((0, .1), (.1, .25), (.25, .5), (.5, 1.0))]
        print(f"   {u:18s} " + ", ".join(cells))


def main():
    on, off = bk.laser_pulses()
    buf = io.StringIO()
    with redirect_stdout(buf):
        fates(on, off)
        means(on, off)
        per_spike(on)
    print(buf.getvalue())
    (HERE / "results" / "driven_spikes.txt").write_text(buf.getvalue())


if __name__ == "__main__":
    main()
