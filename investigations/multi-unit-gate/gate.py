"""Where to set the multi-unit cosine gate on the final TTM-003-48H sorting: the best template cosine of
the detected events no sorted spike claims, next to that of the events exactly one sorted spike claims and
of random stretches of background, window by window on the sort's own binary.

    .venv/bin/python investigations/multi-unit-gate/gate.py <out dir> [workers] [windows...]

Events are negative peaks at 5.5 times each channel's noise (SpikeInterface "locally_exclusive", radius
40 µm, so one event per trough per tetrode), kept at 7 times the noise or more, outside their tetrode's bad
intervals and the first 30 and last 59 samples of each 1800 s window. Each channel's noise is fixed for the
whole recording: the mean MAD (as a standard deviation) of every 0.5 s chunk outside its tetrode's bad
intervals, from investigations/analyzer-noise-and-duplicates/results/noise_levels_channels.csv. The
analyzer stage draws 2,000 of those chunks, which stays within 1.0 % of that mean on every channel in 95 %
of draws.

An event is unclaimed when no sorted spike of its tetrode is within 0.5 ms of it, and a reference event
when exactly one is. Both are scored the same way: the full-window cosine (1 ms before to 2 ms after the
trough, the tetrode's wires, no shift) to each template of the tetrode, a template being per window the
mean of up to 500 random spikes of a unit with 100 or more spikes there. Background stretches are random
frames at least 1 ms from every detected event and every sorted spike of the tetrode, outside its bad
intervals, scored the same way.

Every event and background stretch is also scored after subtracting the tetrode's sorted spikes around
it. Each unit's extended template is the mean of the same spikes as its template, from 3.5 ms before to
3.5 ms after each (units with fewer than 100 spikes in the window get one from all their spikes, for
subtraction only). From the 3 ms trace, every sorted spike of the tetrode whose extended template reaches
it is subtracted, scaled by least squares on the spike's own samples -10 to +20, except a spike within
0.5 ms of the event (a reference event's claiming spike). The residual is scored as above, and its
amplitude is the most negative sample on any of the tetrode's wires within 3 samples of the event, over
that wire's noise.

<out dir>/wNN.npz holds, per window: ev_* for every unclaimed and reference event (amplitude, claims, the
claiming unit, the number of other tetrodes with an event within 0.5 ms, whether it falls in a laser pulse
or the 100 ms after one, the best and second-best cosine, the best unit, the least-squares scale of the
best template, the cosine to the claiming unit's template; after the subtraction the best cosine, best
unit and scale, the residual amplitude and the spikes subtracted and not subtracted, a unit with no spike
in the window), bg_* for the background stretches (cosine before and after the subtraction), counts of
events of 7 or more excluded or claimed by two or more spikes, the templates, and the raw and residual
traces of up to 3 random events per tetrode, kind, 0.05 bin of the cosine after the subtraction and
amplitude class (7 to 9, 9 or more).
"""

import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from ttm.carry import in_intervals, read_bad_periods, read_pulses

SESSION = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12"
RUN = f"{SESSION}/sortings_carry/main"
BAD = f"{SESSION}/bad_periods/bad_periods.csv"
ROOT = pathlib.Path(__file__).resolve().parents[2]
PULSES = ROOT / "subjects/ttm-003/results/laser_pulses_2026-10-20_09-02-12.csv"
NOISE = ROOT / "investigations/analyzer-noise-and-duplicates/results/noise_levels_channels.csv"
FS = 30000.0
WIN = 54_000_000
BEFORE, AFTER = 30, 60
CLAIM, BG_GAP = 15, 30
DETECT, FLOOR, AMP_SPLIT = 5.5, 7.0, 9.0
MIN_SPIKES, N_TEMPLATE = 100, 500
N_BG = 1000
PULSE_POST = 3000
COS_BINS, N_SNIP = 20, 3
EXT, FIT, HALF = 105, (-10, 21), 3
CHUNK = 50_000
N_GROUPS = 16


def near(q, ref, tol):
    """Bool per frame of ``q``: a frame of the sorted ``ref`` within ``tol``."""
    if ref.size == 0:
        return np.zeros(q.size, bool)
    j = np.searchsorted(ref, q)
    d1 = np.where(j > 0, q - ref[np.clip(j - 1, 0, ref.size - 1)], tol + 1)
    d2 = np.where(j < ref.size, ref[np.clip(j, 0, ref.size - 1)] - q, tol + 1)
    return np.minimum(d1, d2) <= tol


class Setup:
    def __init__(self):
        self.rec = si.load(f"{RUN}/binary")
        self.groups = np.asarray(self.rec.get_property("group")).astype(int)
        self.nfr = self.rec.get_num_frames()
        t = pd.read_csv(NOISE).set_index("channel")
        ids = [str(c) for c in self.rec.channel_ids]
        assert (t.loc[ids, "tetrode"].to_numpy() == self.groups + 1).all()
        self.noise = t.loc[ids, "mean_chunk_mad_clean"].to_numpy(np.float64)
        s = si.load(f"{RUN}/sorting")
        self.uids = np.asarray(s.unit_ids).astype(int)
        ug = np.asarray(s.get_property("group")).astype(int)
        self.ugroup = dict(zip(self.uids.tolist(), ug.tolist()))
        self.trains = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in self.uids.tolist()}
        self.pool_f, self.pool_u = {}, {}
        for g in np.unique(ug).tolist():
            us = self.uids[ug == g].tolist()
            f = np.concatenate([self.trains[u] for u in us])
            u = np.concatenate([np.full(self.trains[u].size, u, np.int32) for u in us])
            o = np.argsort(f, kind="stable")
            self.pool_f[g], self.pool_u[g] = f[o], u[o]
        bad = read_bad_periods(BAD, FS)
        self.bad = {g: (np.array([a for a, _ in v], np.int64), np.array([b for _, b in v], np.int64))
                    for g, v in bad.items() if len(v)}
        on, off = read_pulses(PULSES, FS)
        self.pulse = (on, off + PULSE_POST)

    def in_bad(self, frames, g):
        if g not in self.bad:
            return np.zeros(frames.size, bool)
        return in_intervals(frames, *self.bad[g])


S = None


def init():
    global S
    S = Setup()


def run_window(w, out_dir):
    t0 = time.time()
    rng = np.random.default_rng(w)
    a, b = w * WIN, min((w + 1) * WIN, S.nfr)
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", np.int16, "r").reshape(-1, len(S.groups))
    ra, rb = max(0, a - 2 * EXT), min(S.nfr, b + 2 * EXT)
    x = np.ascontiguousarray(raw[ra:rb])
    rng_snip = np.random.default_rng(10_000 + w)
    peaks = detect_peaks(S.rec.frame_slice(a, b), method="locally_exclusive",
                         method_kwargs={"peak_sign": "neg", "detect_threshold": DETECT, "radius_um": 40.0,
                                        "noise_levels": S.noise},
                         job_kwargs={"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False})
    f = peaks["sample_index"].astype(np.int64) + a
    ch = peaks["channel_index"].astype(np.int64)
    g = S.groups[ch]
    amp = -peaks["amplitude"].astype(np.float64) / S.noise[ch]
    detected = {k: np.sort(f[g == k]) for k in range(N_GROUPS)}
    keep = (f - BEFORE >= a) & (f + AFTER <= b)
    for k in np.unique(g).tolist():
        m = g == k
        keep[m] &= ~S.in_bad(f[m], k)
    n_excluded = np.bincount(g[~keep & (amp >= FLOOR)], minlength=N_GROUPS)
    f, ch, g, amp = f[keep], ch[keep], g[keep], amp[keep]
    others = np.zeros(f.size, np.int8)
    for k in range(N_GROUPS):
        others += near(f, np.sort(f[g == k]), CLAIM).astype(np.int8)
    others -= 1
    claims = np.zeros(f.size, np.int64)
    own = np.full(f.size, -1, np.int32)
    for k, pf in S.pool_f.items():
        m = np.flatnonzero(g == k)
        i0 = np.searchsorted(pf, f[m] - CLAIM, "left")
        i1 = np.searchsorted(pf, f[m] + CLAIM, "right")
        claims[m] = i1 - i0
        one = i1 - i0 == 1
        own[m[one]] = S.pool_u[k][i0[one]]
    big = amp >= FLOOR
    n_events = np.bincount(g[big], minlength=N_GROUPS)
    n_overlap = np.bincount(g[big & (claims >= 2)], minlength=N_GROUPS)
    sel = np.flatnonzero(big & (claims <= 1))
    ev = {"frame": f[sel], "group": g[sel].astype(np.int8), "chan": ch[sel].astype(np.int8),
          "amp": amp[sel].astype(np.float32), "claims": claims[sel].astype(np.int8), "own": own[sel],
          "others": others[sel], "in_pulse": in_intervals(f[sel], *S.pulse)}
    n = sel.size
    ev["best_uid"] = np.full(n, -1, np.int32)
    ev["sub_uid"] = np.full(n, -1, np.int32)
    for k in ("best_cos", "second_cos", "alpha", "own_cos", "sub_cos", "sub_alpha", "res_amp"):
        ev[k] = np.full(n, np.nan, np.float32)
    ev["n_sub"] = np.zeros(n, np.int16)
    ev["n_unsub"] = np.zeros(n, np.int16)

    off = np.arange(-BEFORE, AFTER)
    ext = np.arange(-EXT, EXT + 1)
    fit = np.arange(*FIT)
    tmpl_ids, tmpl_T, bg_g, bg_cos, bg_sub = [], [], [], [], []
    snip_idx, snip, snip_res = [], [], []
    for k in range(N_GROUPS):
        cols = np.flatnonzero(S.groups == k)
        xt = np.ascontiguousarray(x[:, cols])
        nz = S.noise[cols].astype(np.float32)

        def take(frames, rel, xt=xt):
            return xt[np.clip((frames - ra)[:, None] + rel, 0, xt.shape[0] - 1)].astype(np.float32)

        T, ids, E, eids = [], [], [], []
        for u in (u for u in S.uids.tolist() if S.ugroup[u] == k):
            t = S.trains[u]
            t = t[(t - BEFORE >= a) & (t + AFTER <= b)]
            if t.size == 0:
                continue
            scorable = t.size >= MIN_SPIKES
            if scorable and t.size > N_TEMPLATE:
                t = np.sort(rng.choice(t, N_TEMPLATE, replace=False))
            e = take(t, ext).mean(axis=0)
            E.append(e)
            eids.append(u)
            if scorable:
                T.append(e[EXT - BEFORE:EXT + AFTER])
                ids.append(u)
        if not ids:
            continue
        T, ids = np.stack(T), np.asarray(ids, np.int32)
        E = np.stack(E)
        ecol = np.full(S.uids.max() + 1, -1)
        ecol[eids] = np.arange(len(eids))
        pf, pu = S.pool_f[k], S.pool_u[k]
        Tp = np.zeros((ids.size, BEFORE + AFTER, 4), np.float32)
        Tp[:, :, :cols.size] = T
        tmpl_ids.append(ids)
        tmpl_T.append(Tp)
        tn = np.sqrt((T ** 2).sum(axis=(1, 2)))
        col = np.full(S.uids.max() + 1, -1)
        col[ids] = np.arange(ids.size)

        def cosines(s, T=T, tn=tn):
            dots = np.einsum("ntc,utc->nu", s, T)
            sn = np.sqrt(np.einsum("ntc,ntc->n", s, s))
            return dots, dots / np.maximum(sn[:, None] * tn[None], 1e-12)

        def score(frames, take=take, cosines=cosines):
            return cosines(take(frames, off))

        def residual(frames, s, take=take, E=E, ecol=ecol, pf=pf, pu=pu):
            i0 = np.searchsorted(pf, frames - BEFORE - EXT, "left")
            cnt = np.searchsorted(pf, frames + AFTER - 1 + EXT, "right") - i0
            pe = np.repeat(np.arange(frames.size), cnt)
            pj = np.repeat(i0, cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
            sp, su = pf[pj], pu[pj]
            keep = np.abs(sp - frames[pe]) > CLAIM
            known = ecol[su] >= 0
            n_unsub = np.bincount(pe[keep & ~known], minlength=frames.size)
            keep &= known
            pe, sp, ti = pe[keep], sp[keep], ecol[su[keep]]
            model = np.zeros_like(s)
            if sp.size:
                tf = E[ti][:, EXT + FIT[0]:EXT + FIT[1]]
                alpha = np.einsum("ntc,ntc->n", take(sp, fit), tf) / np.maximum(np.einsum("ntc,ntc->n", tf, tf), 1e-12)
                j = (frames[pe] - BEFORE - sp)[:, None] + np.arange(BEFORE + AFTER) + EXT
                inside = (j >= 0) & (j <= 2 * EXT)
                contrib = E[ti[:, None], np.clip(j, 0, 2 * EXT)] * (inside * alpha[:, None])[:, :, None]
                starts = np.flatnonzero(np.r_[True, pe[1:] != pe[:-1]])
                model[pe[starts]] = np.add.reduceat(contrib, starts, axis=0)
            return s - model, np.bincount(pe, minlength=frames.size), n_unsub

        def score_sub(frames, take=take, residual=residual, cosines=cosines):
            res, n_sub, n_unsub = residual(frames, take(frames, off))
            dots, c = cosines(res)
            return res, dots, c, n_sub, n_unsub

        idx = np.flatnonzero(ev["group"] == k)
        for c0 in range(0, idx.size, CHUNK):
            ii = idx[c0:c0 + CHUNK]
            dots, c = score(ev["frame"][ii])
            order = np.argsort(-c, axis=1)
            r = np.arange(ii.size)
            ev["best_uid"][ii] = ids[order[:, 0]]
            ev["best_cos"][ii] = c[r, order[:, 0]]
            if ids.size > 1:
                ev["second_cos"][ii] = c[r, order[:, 1]]
            ev["alpha"][ii] = dots[r, order[:, 0]] / tn[order[:, 0]] ** 2
            oc = np.where(ev["own"][ii] >= 0, col[np.maximum(ev["own"][ii], 0)], -1)
            ev["own_cos"][ii] = np.where(oc >= 0, c[r, np.maximum(oc, 0)], np.nan)
            res, dots, c, n_sub, n_unsub = score_sub(ev["frame"][ii])
            kb = c.argmax(axis=1)
            ev["sub_uid"][ii] = ids[kb]
            ev["sub_cos"][ii] = c[r, kb]
            ev["sub_alpha"][ii] = dots[r, kb] / tn[kb] ** 2
            ev["res_amp"][ii] = (-res[:, BEFORE - HALF:BEFORE + HALF + 1] / nz).max(axis=(1, 2))
            ev["n_sub"][ii] = n_sub
            ev["n_unsub"][ii] = n_unsub

        cand = rng.integers(a + BEFORE, b - AFTER, 4 * N_BG)
        ok = ~near(cand, detected[k], BG_GAP) & ~near(cand, S.pool_f.get(k, np.zeros(0, np.int64)), BG_GAP)
        cand = np.sort(cand[ok & ~S.in_bad(cand, k)][:N_BG])
        bg_g.append(np.full(cand.size, k, np.int8))
        bg_cos.append(score(cand)[1].max(axis=1).astype(np.float32))
        bg_sub.append(score_sub(cand)[2].max(axis=1).astype(np.float32))

        cbin = np.clip((np.nan_to_num(ev["sub_cos"][idx]) * COS_BINS).astype(int), 0, COS_BINS - 1)
        key = (ev["claims"][idx].astype(int) * COS_BINS + cbin) * 2 + (ev["amp"][idx] >= AMP_SPLIT)
        perm = rng_snip.permutation(idx.size)
        perm = perm[np.argsort(key[perm], kind="stable")]
        ks = key[perm]
        rank = np.arange(ks.size) - np.searchsorted(ks, ks, "left")
        si_k = np.sort(idx[perm[rank < N_SNIP]])
        raw_k = take(ev["frame"][si_k], off)
        res_k = residual(ev["frame"][si_k], raw_k)[0]
        pad = np.zeros((si_k.size, BEFORE + AFTER, 4), np.float32)
        snip_idx.append(si_k)
        snip.append(np.concatenate([raw_k, pad[:, :, cols.size:]], axis=2).astype(np.int16))
        snip_res.append(np.concatenate([res_k, pad[:, :, cols.size:]], axis=2).round().astype(np.int16))
    snip_idx = np.concatenate(snip_idx) if snip_idx else np.zeros(0, np.int64)
    o = np.argsort(snip_idx)
    snip_idx = snip_idx[o]
    snip = np.concatenate(snip)[o] if snip else np.zeros((0, BEFORE + AFTER, 4), np.int16)
    snip_res = np.concatenate(snip_res)[o] if snip_res else np.zeros((0, BEFORE + AFTER, 4), np.int16)

    np.savez(pathlib.Path(out_dir) / f"w{w:02d}.npz", a=a, b=b, noise=S.noise, n_events=n_events,
             n_excluded=n_excluded, n_overlap=n_overlap,
             tmpl_ids=np.concatenate(tmpl_ids) if tmpl_ids else np.zeros(0, np.int32),
             tmpl_T=np.concatenate(tmpl_T) if tmpl_T else np.zeros((0, BEFORE + AFTER, 4), np.float32),
             bg_group=np.concatenate(bg_g) if bg_g else np.zeros(0, np.int8),
             bg_cos=np.concatenate(bg_cos) if bg_cos else np.zeros(0, np.float32),
             bg_sub_cos=np.concatenate(bg_sub) if bg_sub else np.zeros(0, np.float32),
             snip_idx=snip_idx, snip=snip, snip_res=snip_res, **{f"ev_{k}": v for k, v in ev.items()})
    un = ev["claims"] == 0
    return w, int(un.sum()), int((~un).sum()), time.time() - t0


def main():
    out_dir = pathlib.Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    n_win = int(np.ceil(si.load(f"{RUN}/binary").get_num_frames() / WIN))
    todo = [int(v) for v in sys.argv[3:]] or list(range(n_win))
    todo = [w for w in todo if not (out_dir / f"w{w:02d}.npz").exists()]
    with ProcessPoolExecutor(workers, initializer=init) as ex:
        for w, nu, nr, dt in ex.map(run_window, todo, [out_dir] * len(todo)):
            print(f"[{time.strftime('%H:%M:%S')}] window {w}: {nu} unclaimed, {nr} reference, {dt:.0f} s",
                  flush=True)


if __name__ == "__main__":
    main()
