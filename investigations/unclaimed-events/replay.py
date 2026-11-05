"""What the fold and multi-unit steps of investigations/spike-sorting/97_production_pipeline.py would do
to the final TTM-003-48H sorting, replayed window by window on the sort's own binary.

    .venv/bin/python investigations/unclaimed-events/replay.py <out dir> [workers] [windows...]

Events are negative peaks at 5.5 times each channel's MAD noise (SpikeInterface "locally_exclusive",
radius 40 µm, so one event per trough per tetrode); an event is claimed when any sorted unit of its
tetrode fires within 0.5 ms of it. The fold takes the unclaimed events of 10 MAD or more whose best
full-window cosine (1 ms before to 2 ms after the trough, the tetrode's wires, no shift) to a template of
the same tetrode is 0.8 or more, and adds each to that unit unless the unit fires within 1.5 ms of it.
Templates are per 1800 s window: the mean of up to 500 random spikes of each unit with 100 or more spikes
there. The multi-unit step then pools the events of 7 MAD or more still unclaimed (folded spikes now
claim too) whose best cosine is 0.55 to 0.8 into one unit per tetrode. In both, an event within 0.3 ms
of the previous one of the same train (the unit, or the tetrode's multi-unit unit) is dropped, as in 97.

Departures from 97:
- 97 took its events from one detection over the whole recording with one SpikeInterface noise estimate
  (20 random chunks from anywhere). Here each window is detected on its own, with each channel's noise
  the mean MAD of 40 random 0.5 s chunks of the window outside its tetrode's bad intervals; the 7 and 10
  MAD floors use the same noise.
- Events inside a bad interval of their tetrode are left out, and so are events in the first 30 or last
  60 samples of a window (the sort has no spikes there, and 97 scores no event there either).
- 97 rebuilt the templates from the sorting with the folded spikes added before the multi-unit step;
  here both steps use the same templates.

<out dir>/wNN.npz holds, per window: ev_* for every unclaimed event of 7 MAD or more (scores, the
matcher's highest-snr fit on the tetrode within 0.5 ms, and that of the best-matching unit, with why
each was rejected, the number of other tetrodes with a kept event within 0.5 ms, the decisions); host_*
for the spikes each template was built from, scored the same way; counts of every kept event by
tetrode, claimed, amplitude band and number of other tetrodes; the noise and the templates; the traces
of every event the fold proposes and of up to 50 random multi-unit events per tetrode.
"""

import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import spikeinterface as si
from spikeinterface.sortingcomponents.peak_detection import detect_peaks
from ttm.carry import in_intervals, read_bad_periods, read_pulses

SESSION = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12"
RUN = f"{SESSION}/sortings_carry/main"
BAD = f"{SESSION}/bad_periods/bad_periods.csv"
PULSES = pathlib.Path(__file__).parents[2] / "subjects/ttm-003/results/laser_pulses_2026-10-20_09-02-12.csv"
FS = 30000.0
WIN = 54_000_000
BEFORE, AFTER, SHIFT = 30, 60, 3
CLAIM, REFR, COINC = 15, 45, 9
DETECT, MUA_MAD, FOLD_MAD = 5.5, 7.0, 10.0
MUA_COS, FOLD_COS = 0.55, 0.8
MIN_SPIKES, N_TEMPLATE = 100, 500
N_CHUNKS, CHUNK = 40, 15000
PULSE_POST = 3000
N_MUA_SNIPPETS = 50
BANDS = (5.5, 7.0, 10.0, np.inf)
N_GROUPS = 16
MAD_TO_SD = 0.6744897501960817


def near(q, ref, tol):
    """Bool per frame of ``q``: a frame of the sorted ``ref`` within ``tol``."""
    if ref.size == 0:
        return np.zeros(q.size, bool)
    j = np.searchsorted(ref, q)
    d1 = np.where(j > 0, q - ref[np.clip(j - 1, 0, ref.size - 1)], tol + 1)
    d2 = np.where(j < ref.size, ref[np.clip(j, 0, ref.size - 1)] - q, tol + 1)
    return np.minimum(d1, d2) <= tol


def strongest(evf, evg, fits, fit_g, unit=None):
    """Per event, the index into ``fits`` of the highest-snr matcher fit on its tetrode within CLAIM
    (only fits of ``unit[i]`` if given), or -1."""
    out = np.full(evf.size, -1, np.int64)
    snr = np.nan_to_num(fits["snr"], nan=-np.inf)
    for g in np.unique(evg):
        m = np.flatnonzero(evg == g)
        fm = np.flatnonzero(fit_g == g)
        if fm.size == 0:
            continue
        o = fm[np.argsort(fits["frame"][fm], kind="stable")]
        fr = fits["frame"][o]
        i0 = np.searchsorted(fr, evf[m] - CLAIM, "left")
        cnt = np.searchsorted(fr, evf[m] + CLAIM, "right") - i0
        e = np.repeat(np.arange(m.size), cnt)
        j = o[np.repeat(i0, cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)]
        if unit is not None:
            ok = fits["unit"][j] == unit[m][e]
            e, j = e[ok], j[ok]
        order = np.lexsort((-snr[j], e))
        e, j = e[order], j[order]
        first = np.r_[True, np.diff(e) > 0] if e.size else np.zeros(0, bool)
        out[m[e[first]]] = j[first]
    return out


def event_snippets(x, lo, ev, idx):
    """(n, BEFORE + AFTER, 4) traces (raw units) around the events ``idx`` on their tetrode's wires."""
    out = np.zeros((idx.size, BEFORE + AFTER, 4), np.float32)
    off = np.arange(-BEFORE, AFTER)
    for g in np.unique(ev["group"][idx]):
        pos = np.flatnonzero(ev["group"][idx] == g)
        cols = np.flatnonzero(S.groups == g)
        out[pos, :, :cols.size] = x[(ev["frame"][idx[pos]] - lo)[:, None] + off][:, :, cols]
    return out


def drop_close(frames):
    """Mask keeping a sorted frame unless it is within COINC of the previous frame (97's collapse)."""
    return np.r_[True, np.diff(frames) > COINC] if frames.size else np.zeros(0, bool)


def cosines(ext, T):
    """(n, units, shifts) cosine of each event's snippet, moved by -SHIFT..SHIFT samples, to each template."""
    tn = np.sqrt((T ** 2).sum(axis=(1, 2)))
    out = np.empty((ext.shape[0], len(T), 2 * SHIFT + 1), np.float32)
    for i, k in enumerate(range(-SHIFT, SHIFT + 1)):
        s = ext[:, SHIFT + k:SHIFT + k + BEFORE + AFTER]
        sn = np.sqrt(np.einsum("ntc,ntc->n", s, s))
        out[:, :, i] = np.einsum("ntc,utc->nu", s, T) / np.maximum(sn[:, None] * tn[None], 1e-12)
    return out


class Setup:
    def __init__(self):
        self.rec = si.load(f"{RUN}/binary")
        self.groups = np.asarray(self.rec.get_property("group")).astype(int)
        self.nfr = self.rec.get_num_frames()
        s = si.load(f"{RUN}/sorting")
        self.uids = np.asarray(s.unit_ids).astype(int)
        self.ugroup = dict(zip(self.uids.tolist(), np.asarray(s.get_property("group")).astype(int).tolist()))
        self.group_of = np.full(self.uids.max() + 1, -1)
        self.group_of[self.uids] = [self.ugroup[u] for u in self.uids.tolist()]
        self.trains = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in self.uids.tolist()}
        self.pooled = {g: np.sort(np.concatenate([t for u, t in self.trains.items() if self.ugroup[u] == g]))
                       for g in set(self.ugroup.values())}
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


def window_noise(x, lo, a, b, rng):
    """Per channel, the mean MAD (as a standard deviation) of N_CHUNKS random 0.5 s chunks of [a, b)
    outside the channel's tetrode's bad intervals; NaN where there are none."""
    starts = rng.integers(a, b - CHUNK, 20 * N_CHUNKS)
    noise = np.full(len(S.groups), np.nan)
    for g in np.unique(S.groups):
        ok = starts
        if g in S.bad:
            st, en = S.bad[g]
            i = np.searchsorted(st, starts + CHUNK, side="left")
            j = np.searchsorted(en, starts, side="right")
            ok = starts[i == j]
        ok = ok[:N_CHUNKS]
        if ok.size == 0:
            continue
        ch = np.flatnonzero(S.groups == g)
        m = []
        for s0 in ok:
            c = x[s0 - lo:s0 - lo + CHUNK][:, ch].astype(np.float32)
            m.append(np.median(np.abs(c - np.median(c, axis=0)), axis=0) / MAD_TO_SD)
        noise[ch] = np.mean(m, axis=0)
    return noise


def run_window(w, out_dir):
    t0 = time.time()
    rng = np.random.default_rng(w)
    a, b = w * WIN, min((w + 1) * WIN, S.nfr)
    lo, hi = max(a - SHIFT, 0), min(b + SHIFT, S.nfr)
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", np.int16, "r").reshape(-1, len(S.groups))
    x = np.ascontiguousarray(raw[lo:hi])
    noise = window_noise(x, lo, a, b, rng)
    thr_noise = np.where(np.isfinite(noise), noise, 1e9)
    peaks = detect_peaks(S.rec.frame_slice(a, b), method="locally_exclusive",
                         method_kwargs={"peak_sign": "neg", "detect_threshold": DETECT, "radius_um": 40.0,
                                        "noise_levels": thr_noise},
                         job_kwargs={"n_jobs": 1, "chunk_duration": "1s", "progress_bar": False})
    f_all = peaks["sample_index"].astype(np.int64) + a
    ch_all = peaks["channel_index"].astype(np.int64)
    g_all = S.groups[ch_all]
    amp_all = -peaks["amplitude"].astype(np.float64) / thr_noise[ch_all]
    keep = (f_all - a >= BEFORE) & (f_all + AFTER <= b) & (f_all - BEFORE - SHIFT >= lo) & (f_all + AFTER + SHIFT <= hi)
    for g in np.unique(g_all):
        m = g_all == g
        keep[m] &= ~S.in_bad(f_all[m], g)
    n_excluded = np.bincount(g_all[~keep], minlength=N_GROUPS)
    f_all, ch_all, g_all, amp_all = f_all[keep], ch_all[keep], g_all[keep], amp_all[keep]
    by_g = {g: np.sort(f_all[g_all == g]) for g in range(N_GROUPS)}
    others = np.zeros(f_all.size, np.int8)
    for g in range(N_GROUPS):
        others += near(f_all, by_g[g], CLAIM).astype(np.int8)
    others -= 1
    claimed = np.zeros(f_all.size, bool)
    for g, ref in S.pooled.items():
        m = g_all == g
        claimed[m] = near(f_all[m], ref, CLAIM)
    band = np.digitize(amp_all, BANDS[1:-1])
    counts = np.zeros((N_GROUPS, 2, len(BANDS) - 1, N_GROUPS), np.int64)
    np.add.at(counts, (g_all, claimed.astype(int), band, others), 1)

    with np.load(f"{RUN}/carry/spikes_w{w:02d}.npz") as z:
        fits = {k: z[k] for k in ("frame", "unit", "reason", "snr", "a", "cos_pl")}
    fit_g = S.group_of[fits["unit"]]

    sel = np.flatnonzero(~claimed & (amp_all >= MUA_MAD))
    ev = {"frame": f_all[sel], "group": g_all[sel].astype(np.int8), "chan": ch_all[sel].astype(np.int8),
          "amp_mad": amp_all[sel].astype(np.float32), "others": others[sel],
          "in_pulse": in_intervals(f_all[sel], *S.pulse)}
    n = sel.size
    for k, dt in (("best_uid", np.int32), ("best_shift_uid", np.int32), ("best_shift", np.int8)):
        ev[k] = np.full(n, -1, dt)
    for k in ("best_cos", "second_cos", "best_shift_cos", "alpha"):
        ev[k] = np.full(n, np.nan, np.float32)
    host = {k: [] for k in ("uid", "frame", "amp_mad", "own_cos", "best_cos", "best_uid", "own_shift_cos",
                            "best_shift_uid", "alpha")}
    tmpl_ids, tmpl_T = [], []
    off = np.arange(-BEFORE - SHIFT, AFTER + SHIFT)
    for g in range(N_GROUPS):
        m = ev["group"] == g
        cols = np.flatnonzero(S.groups == g)
        nz = noise[cols].astype(np.float32)
        xt = np.ascontiguousarray(x[:, cols])
        units = [u for u in S.uids.tolist() if S.ugroup[u] == g]
        T, ids, samples = [], [], []
        for u in units:
            t = S.trains[u]
            t = t[(t >= a) & (t < b)]
            if t.size < MIN_SPIKES:
                continue
            t = np.sort(rng.choice(t, N_TEMPLATE, replace=False)) if t.size > N_TEMPLATE else t
            sn = xt[(t - lo)[:, None] + off].astype(np.float32)
            T.append(sn[:, SHIFT:SHIFT + BEFORE + AFTER].mean(axis=0))
            ids.append(u)
            samples.append((u, t, sn))
        if not ids:
            continue
        T = np.stack(T)
        ids = np.asarray(ids)
        Tp = np.zeros((len(ids), BEFORE + AFTER, 4), np.float32)
        Tp[:, :, :len(cols)] = T
        tmpl_ids.append(ids)
        tmpl_T.append(Tp)
        tsq = (T ** 2).sum(axis=(1, 2))
        for u, t, sn in samples:
            c = cosines(sn, T)
            iu = int(np.flatnonzero(ids == u)[0])
            c0 = c[:, :, SHIFT]
            cs = c.max(axis=2)
            s0 = sn[:, SHIFT:SHIFT + BEFORE + AFTER]
            host["uid"].append(np.full(t.size, u, np.int32))
            host["frame"].append(t)
            host["amp_mad"].append((-(s0[:, BEFORE - 5:BEFORE + 6] / nz).min(axis=(1, 2))).astype(np.float32))
            host["own_cos"].append(c0[:, iu])
            host["best_cos"].append(c0.max(axis=1))
            host["best_uid"].append(ids[c0.argmax(axis=1)].astype(np.int32))
            host["own_shift_cos"].append(cs[:, iu])
            host["best_shift_uid"].append(ids[cs.argmax(axis=1)].astype(np.int32))
            host["alpha"].append((np.einsum("ntc,tc->n", s0, T[iu]) / tsq[iu]).astype(np.float32))
        idx = np.flatnonzero(m)
        for c0_ in range(0, idx.size, 200_000):
            ii = idx[c0_:c0_ + 200_000]
            ext = xt[(ev["frame"][ii] - lo)[:, None] + off].astype(np.float32)
            c = cosines(ext, T)
            c0 = c[:, :, SHIFT]
            order = np.argsort(-c0, axis=1)
            r = np.arange(ii.size)
            ev["best_uid"][ii] = ids[order[:, 0]]
            ev["best_cos"][ii] = c0[r, order[:, 0]]
            if len(ids) > 1:
                ev["second_cos"][ii] = c0[r, order[:, 1]]
            flat = c.reshape(ii.size, -1).argmax(axis=1)
            ev["best_shift_uid"][ii] = ids[flat // (2 * SHIFT + 1)]
            ev["best_shift"][ii] = (flat % (2 * SHIFT + 1)) - SHIFT
            ev["best_shift_cos"][ii] = c.reshape(ii.size, -1).max(axis=1)
            s0 = ext[:, SHIFT:SHIFT + BEFORE + AFTER]
            ev["alpha"][ii] = np.einsum("ntc,ntc->n", s0, T[order[:, 0]]) / tsq[order[:, 0]]

    for prefix, unit in (("fit", None), ("hfit", ev["best_uid"])):
        j = strongest(ev["frame"], ev["group"], fits, fit_g, unit)
        hit, jj = j >= 0, np.maximum(j, 0)
        ev[f"{prefix}_unit"] = np.where(hit, fits["unit"][jj], -1).astype(np.int32)
        ev[f"{prefix}_reason"] = np.where(hit, fits["reason"][jj], -1).astype(np.int8)
        ev[f"{prefix}_dt"] = np.where(hit, fits["frame"][jj] - ev["frame"], 0).astype(np.int32)
        for k in ("snr", "a", "cos_pl"):
            ev[f"{prefix}_{k}"] = np.where(hit, fits[k][jj], np.nan).astype(np.float32)

    valid = ev["best_uid"] >= 0
    propose = valid & (ev["amp_mad"] >= FOLD_MAD) & (ev["best_cos"] >= FOLD_COS)
    refr_ok = np.zeros(n, bool)
    for u in np.unique(ev["best_uid"][propose]):
        i = np.flatnonzero(propose & (ev["best_uid"] == u))
        refr_ok[i] = ~near(ev["frame"][i], S.trains[int(u)], REFR)
    folded = propose & refr_ok
    folded_kept = np.zeros(n, bool)
    for u in np.unique(ev["best_uid"][folded]):
        i = np.flatnonzero(folded & (ev["best_uid"] == u))
        i = i[np.argsort(ev["frame"][i], kind="stable")]
        folded_kept[i[drop_close(ev["frame"][i])]] = True
    claimed_after = folded_kept.copy()
    for g in range(N_GROUPS):
        m = (ev["group"] == g) & ~folded_kept
        fk = np.sort(ev["frame"][folded_kept & (ev["group"] == g)])
        claimed_after[m] = near(ev["frame"][m], fk, CLAIM)
    mua = valid & ~claimed_after & (ev["best_cos"] >= MUA_COS) & (ev["best_cos"] < FOLD_COS)
    mua_kept = np.zeros(n, bool)
    for g in range(N_GROUPS):
        i = np.flatnonzero(mua & (ev["group"] == g))
        i = i[np.argsort(ev["frame"][i], kind="stable")]
        mua_kept[i[drop_close(ev["frame"][i])]] = True
    ev.update(propose=propose, refr_ok=refr_ok, folded=folded_kept, claimed_after=claimed_after, mua=mua_kept)
    host = {k: np.concatenate(v) if v else np.zeros(0) for k, v in host.items()}
    prop_idx = np.flatnonzero(propose)
    mua_idx = np.sort(np.concatenate([rng.permutation(np.flatnonzero(mua_kept & (ev["group"] == g)))[:N_MUA_SNIPPETS]
                                      for g in range(N_GROUPS)]))
    np.savez(pathlib.Path(out_dir) / f"w{w:02d}.npz", noise=noise, counts=counts, n_excluded=n_excluded,
             tmpl_ids=np.concatenate(tmpl_ids) if tmpl_ids else np.zeros(0, int),
             tmpl_T=np.concatenate(tmpl_T) if tmpl_T else np.zeros((0, BEFORE + AFTER, 4), np.float32),
             prop_idx=prop_idx, prop_snip=event_snippets(x, lo, ev, prop_idx),
             mua_idx=mua_idx, mua_snip=event_snippets(x, lo, ev, mua_idx),
             **{f"ev_{k}": v for k, v in ev.items()}, **{f"host_{k}": v for k, v in host.items()})
    return w, n, int(folded_kept.sum()), int(mua_kept.sum()), time.time() - t0


def main():
    out_dir = pathlib.Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    n_win = int(np.ceil(si.load(f"{RUN}/binary").get_num_frames() / WIN))
    todo = [int(v) for v in sys.argv[3:]] or list(range(n_win))
    todo = [w for w in todo if not (out_dir / f"w{w:02d}.npz").exists()]
    with ProcessPoolExecutor(workers, initializer=init) as ex:
        for w, n, nf, nm, dt in ex.map(run_window, todo, [out_dir] * len(todo)):
            print(f"[{time.strftime('%H:%M:%S')}] window {w}: {n} unclaimed >= {MUA_MAD} MAD, "
                  f"{nf} folded, {nm} multi-unit, {dt:.0f} s", flush=True)


if __name__ == "__main__":
    main()
