"""Hybrid ground truth: known units added to the first 12 windows (6 h) of the 48 h binary.

Each ground-truth unit is a seed-bank template from another tetrode with its
four wires permuted (or, for the "similar" units, a destination template
mixed with another until their cosine is about 0.85), scaled so its trough
on its deepest wire is a target number of excised-MAD noise units of the
destination tetrode in window 0. Over the 6 h each unit has:

- a firing rate that swings by half around its mean over 30 to 90 min, with a
  2.5 ms dead time, and, for two units, no spikes from window 6 on (a cell
  that leaves; whatever its template collects afterwards is false);
- bursts: a share of spikes followed by one to three spikes 3 to 7 ms apart
  at 0.85, 0.7 and 0.6 of full size;
- a slow change of size (flat, falling to 0.6, rising to 1.3, a dip to 0.7,
  or a step to 0.75 at 2.5 h) and of the four wires' ratios (each wire's gain
  moves linearly by up to 20 %), the latter constant within a window;
- spike-to-spike size jitter (7 %) and sub-sample timing (quarter samples);
- during every delivered laser pulse: silence from 5 ms after onset to 60 ms
  after offset, or, for the three driven units, a full-size spike 3.5 to 8 ms
  after onset, two more at 0.85 and 0.7, then a train at about 140 per second
  at half size until offset; the share of pulses that drive the train falls
  over the four trains (1, 1, 0.6, 0.3).

Writes the hybrid binary (float32, 64 channels, frames 0 to 12 windows of the
original) to ``_bake.HYBRID`` and ``gt.npz`` beside it: per spike the unit,
frame, class (0 isolated, 1 burst lead, 2 burst tail, 3 laser first spikes,
4 laser half-size) and total size factor; per unit and window the injected
template (without the slow size factor) and that factor at the window's middle.

Run: .venv/bin/python investigations/carry-bakeoff/make_hybrid.py [--jobs 64]
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import spikeinterface as si
from spikeinterface.core.generate import InjectTemplatesRecording

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402
from ttm import capture as cp  # noqa: E402

N_WIN = 12
SPAN = N_WIN * bk.WIN
UP = 4
# (tetrode index, trough MAD, mean rate Hz, burst share, size course, kind)
PLAN = [
    (0, 5.0, 6.0, 0.10, "down", "plain"),
    (0, 9.0, 2.0, 0.15, "flat", "plain"),
    (2, 5.5, 15.0, 0.05, "up", "plain"),
    (2, 12.0, 1.0, 0.20, "dip", "plain"),
    (8, 5.0, 10.0, 0.10, "flat", "plain"),
    (8, 7.0, 3.0, 0.15, "step", "plain"),
    (12, 9.0, 4.0, 0.10, "dip", "driven"),
    (12, 10.0, 8.0, 0.10, "flat", "leaves"),
    (13, 5.0, 2.0, 0.10, "flat", "leaves"),
    (13, 8.0, 12.0, 0.05, "down", "plain"),
    (15, 5.5, 5.0, 0.10, "flat", "plain"),
    (15, 14.0, 3.0, 0.20, "up", "plain"),
    (3, 6.0, 6.0, 0.10, "down", "similar"),
    (5, 10.0, 3.0, 0.15, "flat", "driven"),
    (5, 9.0, 8.0, 0.10, "step", "plain"),
    (6, 5.0, 10.0, 0.10, "dip", "plain"),
    (9, 5.0, 4.0, 0.10, "flat", "driven"),
    (10, 7.0, 5.0, 0.10, "up", "similar"),
    (11, 5.5, 8.0, 0.10, "down", "plain"),
    (4, 12.0, 4.0, 0.15, "flat", "plain"),
    (7, 6.0, 6.0, 0.10, "down", "plain"),
]
DRIVE_P = [1.0, 1.0, 0.6, 0.3]


def size_course(kind: str, t_h: np.ndarray) -> np.ndarray:
    T = SPAN / bk.FS / 3600
    if kind == "flat":
        return np.ones_like(t_h)
    if kind == "down":
        return 1 - 0.4 * t_h / T
    if kind == "up":
        return 1 + 0.3 * t_h / T
    if kind == "dip":
        return 1 - 0.3 * np.sin(np.pi * t_h / T)
    if kind == "step":
        return np.where(t_h < 2.5, 1.0, 0.75)
    raise ValueError(kind)


def fractional_shift(t: np.ndarray, d: float) -> np.ndarray:
    """``t`` (N_SAMP, 4) delayed by ``d`` samples (band-limited, zero padded)."""
    n = 256
    pad = np.zeros((n, t.shape[1]))
    pad[64:64 + t.shape[0]] = t
    f = np.fft.rfftfreq(n)
    out = np.fft.irfft(np.fft.rfft(pad, axis=0) * np.exp(-2j * np.pi * f * d)[:, None], n=n, axis=0)
    return out[64:64 + t.shape[0]]


def taper(t: np.ndarray, k: int = 6) -> np.ndarray:
    w = np.ones(t.shape[0])
    ramp = 0.5 - 0.5 * np.cos(np.pi * np.arange(k) / k)
    w[:k], w[-k:] = ramp, ramp[::-1]
    return t * w[:, None]


def trains(rng, rate, burst, leaves, driven, pulses, nfr):
    """Frames, classes and class size factors of one unit over [0, nfr)."""
    on, off = pulses
    span_s = nfr / bk.FS
    P, ph = rng.uniform(1800, 5400), rng.uniform(0, 2 * np.pi)
    n = rng.poisson(1.5 * rate * span_s)
    t = np.sort(rng.uniform(0, span_s, n))
    keep = rng.uniform(0, 1.5, n) < 1 + 0.5 * np.sin(2 * np.pi * t / P + ph)
    t = t[keep]
    t = t[np.r_[True, np.diff(t) > 0.0025]]
    if leaves:
        t = t[t < 6 * 1800]
    fr = (t * bk.FS).astype(np.int64)
    cls = np.zeros(fr.size, np.int8)
    mult = np.ones(fr.size)
    # bursts
    lead = rng.uniform(size=fr.size) < burst
    add_f, add_c, add_m = [], [], []
    for f in fr[lead]:
        k = rng.integers(1, 4)
        isi = np.cumsum(rng.uniform(3, 7, k)) * 1e-3 * bk.FS
        add_f.append(f + isi.astype(np.int64))
        add_c.append(np.full(k, 2, np.int8))
        add_m.append(np.array([0.85, 0.7, 0.6][:k]))
    cls[lead] = 1
    if add_f:
        fr = np.r_[fr, np.concatenate(add_f)]
        cls = np.r_[cls, np.concatenate(add_c)]
        mult = np.r_[mult, np.concatenate(add_m)]
    # silence in and after pulses
    quiet = bk.in_intervals(fr, on + int(0.005 * bk.FS), off + int(0.060 * bk.FS))
    fr, cls, mult = fr[~quiet], cls[~quiet], mult[~quiet]
    if driven:
        lf, lc, lm = [], [], []
        train_of = np.searchsorted(np.r_[on[np.r_[True, np.diff(on) > 10 * bk.FS]], np.inf], on, side="right") - 1
        for o, e, tr in zip(on, off, train_of):
            if leaves and o >= 6 * bk.WIN or rng.uniform() > DRIVE_P[min(tr, 3)]:
                continue
            s = o + int(rng.uniform(3.5, 8) * 1e-3 * bk.FS)
            j = 0
            while s < e:
                lf.append(s)
                lc.append(3 if j < 3 else 4)
                lm.append([1.0, 0.85, 0.7][j] if j < 3 else rng.normal(0.5, 0.06))
                s += int(max(rng.normal(7.0, 1.5), 4.0) * 1e-3 * bk.FS)
                j += 1
        fr = np.r_[fr, np.array(lf, np.int64)]
        cls = np.r_[cls, np.array(lc, np.int8)]
        mult = np.r_[mult, np.array(lm)]
    o = np.argsort(fr, kind="stable")
    fr, cls, mult = fr[o], cls[o], mult[o]
    ok = np.r_[True, np.diff(fr) > int(0.002 * bk.FS)] & (fr > bk.N_SAMP) & (fr < nfr - bk.N_SAMP)
    return fr[ok], cls[ok], mult[ok]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=64)
    ap.add_argument("--seed", type=int, default=20261028)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rec = si.load(bk.BINARY)
    groups = np.asarray(rec.get_property("group")).astype(int)
    bank = bk.Bank.load(bk.BAKE_ROOT / "bank0.npz")
    noise0 = {g: cp.noise_mad_excised(cp.load_tetrode_traces(rec, g, 0, int(300 * bk.FS)))[0] for g in range(16)}
    pulses = bk.laser_pulses()
    pulses = (pulses[0][pulses[0] < SPAN], pulses[1][pulses[0] < SPAN])
    n_gt = len(PLAN)
    base = np.zeros((n_gt, bk.N_SAMP, 4))
    info = []
    for i, (g, mad, rate, burst, course, kind) in enumerate(PLAN):
        dest = bank.T[bank.group == g].astype(np.float64)
        for _ in range(200):
            if kind == "similar":
                t0 = dest[rng.integers(len(dest))]
                other = bank.T[rng.choice(np.flatnonzero(bank.group != g))][:, rng.permutation(4)].astype(np.float64)
                for beta in np.linspace(0.05, 2, 60):
                    t = t0 / np.linalg.norm(t0) + beta * other / np.linalg.norm(other)
                    if cp.template_cosine(t, t0) < 0.86:
                        break
            else:
                src = rng.choice(np.flatnonzero(bank.group != g))
                t = bank.T[src][:, rng.permutation(4)].astype(np.float64)
            cmax = max([cp.template_cosine(t, d) for d in dest] + [cp.template_cosine(t, base[j]) for j in range(i)
                                                                   if PLAN[j][0] == g])
            if (kind == "similar" and 0.8 <= cmax <= 0.9) or (kind != "similar" and cmax < 0.8):
                break
        t = taper(t)
        c = int(np.argmin(t.min(axis=0)))
        t = t * (mad * noise0[g][c] / -t[:, c].min())
        base[i] = t
        info.append({"gt": i, "group": g, "tetrode": f"TT{g + 1}", "trough_mad": mad, "rate_hz": rate, "burst": burst,
                     "course": course, "kind": kind, "max_cos_to_tetrode": round(float(cmax), 3)})
    drift = rng.uniform(-0.2, 0.2, (n_gt, 4))
    wmid_h = (np.arange(N_WIN) + 0.5) * 0.5
    T_win = np.stack([np.stack([base[i] * (1 + drift[i] * w / (N_WIN - 1))[None, :] for w in range(N_WIN)])
                      for i in range(n_gt)])
    size_win = np.stack([size_course(p[4], wmid_h) for p in PLAN])
    # spikes
    F, U, C, M = [], [], [], []
    for i, (g, mad, rate, burst, course, kind) in enumerate(PLAN):
        fr, cls, mult = trains(rng, rate, burst, kind == "leaves", kind == "driven", pulses, SPAN)
        F.append(fr), U.append(np.full(fr.size, i)), C.append(cls), M.append(mult)
    F, U, C, M = map(np.concatenate, (F, U, C, M))
    o = np.argsort(F, kind="stable")
    F, U, C, M = F[o], U[o], C[o], M[o]
    size = M.copy()
    for i, p in enumerate(PLAN):
        m = U == i
        size[m] *= size_course(p[4], F[m] / bk.FS / 3600)
    size *= rng.normal(1.0, 0.07, size.size)
    w = np.minimum(F // bk.WIN, N_WIN - 1)
    sub = U * N_WIN + w
    up = rng.integers(0, UP, F.size)
    # dense injection templates: (n_gt * N_WIN, N_SAMP, 64, UP)
    inj = np.zeros((n_gt * N_WIN, bk.N_SAMP, len(groups), UP), np.float32)
    for i in range(n_gt):
        ch = np.flatnonzero(groups == PLAN[i][0])
        for wi in range(N_WIN):
            for k in range(UP):
                inj[i * N_WIN + wi][:, ch, k] = fractional_shift(T_win[i, wi], k / UP)
    sorting = si.NumpySorting.from_samples_and_labels([F], [sub], sampling_frequency=bk.FS,
                                                      unit_ids=np.arange(n_gt * N_WIN))
    sv = sorting.to_spike_vector()
    assert np.array_equal(sv["sample_index"], F) and np.array_equal(sv["unit_index"], sub)
    span = rec.frame_slice(0, SPAN)
    hyb = InjectTemplatesRecording(sorting, inj, nbefore=bk.NB, amplitude_factor=size.astype(np.float32),
                                   parent_recording=span, upsample_vector=up)
    for key in span.get_property_keys():
        hyb.set_property(key, span.get_property(key))
    bk.BAKE.mkdir(parents=True, exist_ok=True)
    np.savez(bk.BAKE / "gt.npz", frame=F, unit=U, cls=C, size=size, up=up, T_win=T_win, size_win=size_win,
             group=np.array([p[0] for p in PLAN]), drift=drift)
    (bk.BAKE / "gt_plan.json").write_text(json.dumps(info, indent=2))
    for r in info:
        print(r, f"spikes {int((U == r['gt']).sum()):,}")
    print(f"{F.size:,} ground-truth spikes; classes {np.bincount(C).tolist()}")
    bk.log(f"writing {bk.HYBRID}")
    hyb.save(folder=bk.HYBRID, format="binary", dtype="float32", n_jobs=a.jobs, chunk_duration="1s",
             progress_bar=False)
    bk.log("done")


if __name__ == "__main__":
    main()
