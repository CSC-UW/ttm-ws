"""What the events with a sorted spike of their tetrode 0.5 to 1.5 ms away are: whether each survives
subtracting the waveforms of the sorted spikes around it, and how it scores once they are gone.

    .venv/bin/python investigations/multi-unit-gate/overlaps.py <gate out dir> <out dir> [workers] [--no-figure]

Per window and unit, the extended template is the mean of up to 500 of the unit's spikes from 3.5 ms
before to 3.5 ms after each. The events are those of TT5 to TT12 that gate.py scored and that have a sorted
spike of their tetrode 16 to 45 samples (0.5 to 1.5 ms) before or after them: every such unclaimed event,
and up to 2,000 reference events per tetrode and window as a control. From each event's window (1 ms before
to 2 ms after it, the tetrode's wires) every sorted spike of the tetrode whose extended template reaches the
window is subtracted, as its unit's extended template scaled by least squares on the spike's own samples
-10 to +20; a reference event's claiming spike is not subtracted. Then:

- the amplitude is the most negative sample on any of the tetrode's wires within 3 samples of the event,
  over that wire's noise, before (raw) and after (residual) the subtraction;
- the predicted amplitude is the subtracted waveforms' value at the event's sample and channel, over the
  noise, so the amplitude the neighbouring spikes alone would give there;
- the cosine is the best full-window cosine, as in gate.py, to the cores (-30 to +59) of the extended
  templates of units with 100 or more spikes in the window, before and after the subtraction.

<out dir>/wNN.npz holds these per event, with the nearest sorted spike 16 to 45 samples after and before
the event (lag and unit), the number of spikes subtracted and not subtracted (unit with no spike in the
window), and the raw and residual traces of up to 2 random unclaimed events per tetrode and window that score
below the tetrode's cutoff losing 10 % and have a sorted spike 16 to 45 samples after them.

Writes results/overlaps.csv: per tetrode, kind (reference; unclaimed grouped by the tetrode's cutoffs losing
1, 5 and 10 %, results/tetrodes.csv) and side (a sorted spike after only, before only, both), the events
of that kind, the events with a sorted spike on that side, the median predicted amplitude, the shares whose
residual amplitude is below 5.5 and 7 or more, the median cosine before and after, the share whose residual
cosine reaches the cutoff losing 10 %, and the share whose best unit after the subtraction is the unit of
the nearest spike on that side, against the share expected if the two units were unrelated (the sum over
units of the product of their shares among best units and among nearest spikes).

Writes results/overlaps_rescored.csv: per tetrode, the unclaimed events below the cutoff losing 10 % and
those of them with a sorted spike 0.5 to 1.5 ms away, the same-unit shares above for those with a spike after
only, and for the cutoffs losing 5 and 10 % the unclaimed events kept, those below the cutoff whose residual
reaches it with an amplitude of 7 or more (added), and those at or above it whose residual does not
(dropped). Draws overlaps.png into <out dir>.
"""

import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gate import AFTER, BEFORE, CLAIM, DETECT, FLOOR, MIN_SPIKES, N_TEMPLATE, RUN, Setup

RESULTS = pathlib.Path(__file__).parent / "results"
EXT = 105
NEAR = 45
FIT = (-10, 21)
HALF = 3
N_REF = 2000
N_SNIP = 2
LIVE = range(4, 12)
GAIN = 0.195
S = None


def init():
    global S
    S = Setup()


def nearest(p, pu, f, lo, hi):
    """Lag and unit of the sorted spike nearest each frame at lags lo..hi (by |lag|), or 0 and -1."""
    lag = np.zeros(f.size, np.int64)
    unit = np.full(f.size, -1, np.int32)
    if lo > 0:
        j = np.searchsorted(p, f + lo, "left")
        ok = (j < p.size) & (p[np.minimum(j, p.size - 1)] <= f + hi)
    else:
        j = np.searchsorted(p, f + hi, "right") - 1
        ok = (j >= 0) & (p[np.maximum(j, 0)] >= f + lo)
    lag[ok] = p[j[ok]] - f[ok]
    unit[ok] = pu[j[ok]]
    return lag, unit


def best(v, Tc, tn, ids):
    """Best cosine of each window in ``v`` to the templates ``Tc`` (norms ``tn``), and its unit."""
    dots = np.einsum("ntc,utc->nu", v, Tc)
    c = dots / np.maximum(np.sqrt(np.einsum("ntc,ntc->n", v, v))[:, None] * tn[None], 1e-12)
    k = c.argmax(axis=1)
    return c[np.arange(v.shape[0]), k], ids[k]


def trough(v, nz):
    """Most negative sample within HALF of the event on any wire, over that wire's noise, sign flipped."""
    return (-v[:, BEFORE - HALF:BEFORE + HALF + 1] / nz).max(axis=(1, 2))


def run_window(w, gate_dir, out_dir, cuts):
    t0 = time.time()
    z = np.load(pathlib.Path(gate_dir) / f"w{w:02d}.npz")
    a, b = int(z["a"]), int(z["b"])
    rng = np.random.default_rng(1000 + w)
    raw = np.memmap(f"{RUN}/binary/traces_cached_seg0.raw", np.int16, "r").reshape(-1, len(S.groups))
    ev = {k[3:]: z[k] for k in z.files if k.startswith("ev_")}
    ext = np.arange(-EXT, EXT + 1)
    off = np.arange(-BEFORE, AFTER)
    fit = np.arange(*FIT)
    out, snips = [], []
    for g in LIVE:
        cols = np.flatnonzero(S.groups == g)
        nz = S.noise[cols].astype(np.float32)
        p, pu = S.pool_f[g], S.pool_u[g]
        units, T, scorable = [], [], []
        for u in (u for u in S.uids.tolist() if S.ugroup[u] == g):
            t = S.trains[u]
            t = t[(t >= a) & (t < b) & (t - EXT >= 0) & (t + EXT < S.nfr)]
            if t.size == 0:
                continue
            scorable.append(t.size >= MIN_SPIKES)
            t = np.sort(rng.choice(t, N_TEMPLATE, replace=False)) if t.size > N_TEMPLATE else t
            T.append(raw[t[:, None] + ext][:, :, cols].astype(np.float32).mean(axis=0))
            units.append(u)
        if not any(scorable):
            continue
        T, units, scorable = np.stack(T), np.asarray(units, np.int32), np.asarray(scorable)
        col = np.full(S.uids.max() + 1, -1)
        col[units] = np.arange(units.size)
        Tc = T[scorable][:, EXT - BEFORE:EXT + AFTER]
        tn = np.sqrt((Tc ** 2).sum(axis=(1, 2)))

        on = np.flatnonzero((ev["group"] == g) & (ev["best_uid"] >= 0))
        f = ev["frame"][on]
        la, ua = nearest(p, pu, f, CLAIM + 1, NEAR)
        lb, ub = nearest(p, pu, f, -NEAR, -CLAIM - 1)
        has = (ua >= 0) | (ub >= 0)
        un = has & (ev["claims"][on] == 0)
        ref = np.flatnonzero(has & (ev["claims"][on] == 1))
        ref = rng.choice(ref, N_REF, replace=False) if ref.size > N_REF else ref
        sel = np.sort(np.concatenate([np.flatnonzero(un), ref]))
        if sel.size == 0:
            continue
        i_ev, f = on[sel], f[sel]
        la, ua, lb, ub = la[sel], ua[sel], lb[sel], ub[sel]
        x = raw[f[:, None] + off][:, :, cols].astype(np.float32)

        i0 = np.searchsorted(p, f - BEFORE - EXT, "left")
        i1 = np.searchsorted(p, f + AFTER - 1 + EXT, "right")
        cnt = i1 - i0
        pe = np.repeat(np.arange(f.size), cnt)
        pj = np.repeat(i0, cnt) + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        s, su = p[pj], pu[pj]
        keep = np.abs(s - f[pe]) > CLAIM
        known = col[su] >= 0
        n_unsub = np.bincount(pe[keep & ~known], minlength=f.size)
        keep &= known
        pe, s, su = pe[keep], s[keep], su[keep]
        n_sub = np.bincount(pe, minlength=f.size)
        ti = col[su]
        if s.size:
            us, inv = np.unique(s, return_inverse=True)
            xf = raw[us[:, None] + fit][:, :, cols].astype(np.float32)[inv]
            tf = T[ti][:, EXT + FIT[0]:EXT + FIT[1]]
            alpha = np.einsum("ntc,ntc->n", xf, tf) / np.maximum(np.einsum("ntc,ntc->n", tf, tf), 1e-12)
            idx = (f[pe] - BEFORE - s)[:, None] + np.arange(BEFORE + AFTER) + EXT
            valid = (idx >= 0) & (idx <= 2 * EXT)
            contrib = T[ti[:, None], np.clip(idx, 0, 2 * EXT)] * (valid * alpha[:, None])[:, :, None]
            model = np.zeros_like(x)
            np.add.at(model, pe, contrib)
        else:
            model = np.zeros_like(x)
        res = x - model

        cos_raw, uid_raw = best(x, Tc, tn, units[scorable])
        cos_res, uid_res = best(res, Tc, tn, units[scorable])
        ch = np.searchsorted(cols, ev["chan"][i_ev])
        r = {"frame": f, "group": np.full(f.size, g, np.int8), "claims": ev["claims"][i_ev],
             "amp": ev["amp"][i_ev], "best_cos": ev["best_cos"][i_ev], "raw_amp": trough(x, nz), "res_amp": trough(res, nz),
             "pred_amp": -model[np.arange(f.size), BEFORE, ch] / nz[ch], "cos_raw": cos_raw, "uid_raw": uid_raw,
             "cos_res": cos_res, "uid_res": uid_res, "lag_after": la, "unit_after": ua, "lag_before": lb,
             "unit_before": ub, "n_sub": n_sub, "n_unsub": n_unsub}
        out.append(r)
        pick = np.flatnonzero((r["claims"] == 0) & (ua >= 0) & (r["best_cos"] < cuts[g + 1]))
        for i in rng.permutation(pick)[:N_SNIP]:
            xp = np.zeros((BEFORE + AFTER, 4), np.float32)
            rp = np.zeros_like(xp)
            xp[:, :cols.size], rp[:, :cols.size] = x[i], res[i]
            snips.append((g, f[i], r["amp"][i], r["best_cos"][i], r["res_amp"][i], r["cos_res"][i], la[i], ua[i],
                          xp, rp))
    keys = out[0].keys()
    arrs = {k: np.concatenate([r[k] for r in out]) for k in keys}
    sn = {"snip_group": np.array([v[0] for v in snips], np.int8), "snip_frame": np.array([v[1] for v in snips]),
          "snip_info": np.array([v[2:8] for v in snips], np.float64).reshape(-1, 6),
          "snip_raw": np.array([v[8] for v in snips], np.float32).reshape(-1, BEFORE + AFTER, 4),
          "snip_res": np.array([v[9] for v in snips], np.float32).reshape(-1, BEFORE + AFTER, 4)}
    np.savez(pathlib.Path(out_dir) / f"w{w:02d}.npz", **arrs, **sn)
    return w, arrs["frame"].size, time.time() - t0


def summarize(gate_dir, out_dir, tet):
    files = sorted(pathlib.Path(out_dir).glob("w*.npz"))
    assert len(files) == 96
    d = {}
    for f in files:
        z = np.load(f)
        for k in z.files:
            if not k.startswith("snip"):
                d.setdefault(k, []).append(z[k])
    d = {k: np.concatenate(v) for k, v in d.items()}
    tot = {}
    for f in sorted(pathlib.Path(gate_dir).glob("w*.npz")):
        z = np.load(f)
        m = (z["ev_best_uid"] >= 0) & (z["ev_group"] >= 4) & (z["ev_group"] < 12)
        tot.setdefault("group", []).append(z["ev_group"][m])
        tot.setdefault("claims", []).append(z["ev_claims"][m])
        tot.setdefault("best_cos", []).append(z["ev_best_cos"][m])
    tot = {k: np.concatenate(v) for k, v in tot.items()}

    def kinds(gr, cl, bc, tt):
        c1, c5, c10 = tet.loc[tt, ["cut_1", "cut_5", "cut_10"]]
        on = gr == tt - 1
        u = on & (cl == 0)
        return {"reference": on & (cl == 1), "unclaimed below 1 %": u & (bc < c1),
                "unclaimed 1 to 5 %": u & (bc >= c1) & (bc < c5), "unclaimed 5 to 10 %": u & (bc >= c5) & (bc < c10),
                "unclaimed at or above 10 %": u & (bc >= c10)}

    def chance(best_u, near_u):
        return float(sum(np.mean(best_u == u) * np.mean(near_u == u) for u in np.unique(np.r_[best_u, near_u])))

    rows, re_rows = [], []
    for tt in range(5, 13):
        on, on_all = (d["group"] == tt - 1) & (d["claims"] == 0), (tot["group"] == tt - 1) & (tot["claims"] == 0)
        c10 = tet.loc[tt, "cut_10"]
        low = on & (d["best_cos"] < c10) & (d["unit_after"] >= 0) & (d["unit_before"] < 0)
        r = {"tetrode": tt, "unclaimed_below_10": int((on_all & (tot["best_cos"] < c10)).sum()),
             "with_neighbour_below_10": int((on & (d["best_cos"] < c10)).sum()),
             "same_unit_below_10_after_only": float(np.mean(d["uid_res"][low] == d["unit_after"][low])),
             "same_unit_chance_below_10_after_only": chance(d["uid_res"][low], d["unit_after"][low])}
        for q in (5, 10):
            c = tet.loc[tt, f"cut_{q}"]
            ok = (d["cos_res"] >= c) & (d["res_amp"] >= FLOOR)
            r[f"kept_{q}"] = int((on_all & (tot["best_cos"] >= c)).sum())
            r[f"added_{q}"] = int((on & (d["best_cos"] < c) & ok).sum())
            r[f"dropped_{q}"] = int((on & (d["best_cos"] >= c) & ~ok).sum())
        re_rows.append(r)
        k_ev = kinds(d["group"], d["claims"], d["best_cos"], tt)
        k_all = kinds(tot["group"], tot["claims"], tot["best_cos"], tt)
        for kind, m in k_ev.items():
            sides = {"after only": (d["unit_after"] >= 0) & (d["unit_before"] < 0),
                     "before only": (d["unit_after"] < 0) & (d["unit_before"] >= 0),
                     "both": (d["unit_after"] >= 0) & (d["unit_before"] >= 0)}
            for side, sm in sides.items():
                mm = m & sm
                if not mm.any():
                    continue
                nu = d["unit_after"] if side != "before only" else d["unit_before"]
                rows.append({"tetrode": tt, "kind": kind, "side": side, "events_of_kind": int(k_all[kind].sum()),
                             "events": int(mm.sum()), "pred_amp_median": float(np.median(d["pred_amp"][mm])),
                             "res_below_5p5": float(np.mean(d["res_amp"][mm] < DETECT)),
                             "res_7_or_more": float(np.mean(d["res_amp"][mm] >= FLOOR)),
                             "cos_raw_median": float(np.median(d["cos_raw"][mm])),
                             "cos_res_median": float(np.median(d["cos_res"][mm])),
                             "cos_res_cut_10": float(np.mean(d["cos_res"][mm] >= c10)),
                             "res_7_and_cut_10": float(np.mean((d["cos_res"][mm] >= c10) & (d["res_amp"][mm] >= FLOOR))),
                             "same_unit": float(np.mean(d["uid_res"][mm] == nu[mm])),
                             "same_unit_chance": chance(d["uid_res"][mm], nu[mm]),
                             "unsubtracted": int(d["n_unsub"][mm].sum())})
    t = pd.DataFrame(rows)
    t.to_csv(RESULTS / "overlaps.csv", index=False)
    re_t = pd.DataFrame(re_rows)
    re_t.to_csv(RESULTS / "overlaps_rescored.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.round(3).to_string(index=False))
    print(re_t.round(3).to_string(index=False))
    print(re_t.drop(columns="tetrode").sum().round(1).to_string())


def figure(out_dir):
    sn = {k: [] for k in ("snip_group", "snip_info", "snip_raw", "snip_res")}
    for f in sorted(pathlib.Path(out_dir).glob("w*.npz")):
        z = np.load(f)
        for k, v in sn.items():
            v.append(z[k])
    sn = {k: np.concatenate(v) for k, v in sn.items()}
    rng = np.random.default_rng(0)
    n_col = 4
    x = np.arange(4 * (BEFORE + AFTER + 2))
    fig, axes = plt.subplots(8, n_col, figsize=(4.2 * n_col, 1.7 * 8))
    for r, g in enumerate(LIVE):
        pick = rng.permutation(np.flatnonzero(sn["snip_group"] == g))[:n_col]
        for c in range(n_col):
            ax = axes[r, c]
            ax.set_xticks([])
            ax.tick_params(labelsize=6)
            if c == 0:
                ax.set_ylabel(f"TT{g + 1}", fontsize=8)
            if c >= pick.size:
                continue
            i = pick[c]
            for v, color in ((sn["snip_raw"][i], "k"), (sn["snip_raw"][i] - sn["snip_res"][i], "r"),
                             (sn["snip_res"][i], "tab:green")):
                pad = np.full((2, 4), np.nan)
                ax.plot(x, np.concatenate([v, pad]).T.ravel() * GAIN, color=color, lw=0.7)
            amp, cos, ramp, rcos, lag, unit = sn["snip_info"][i]
            for wire in range(4):
                ax.axvline(wire * (BEFORE + AFTER + 2) + BEFORE, color="0.7", lw=0.5, zorder=0)
            ax.set_title(f"amp {amp:.1f} -> {ramp:.1f}, cos {cos:.2f} -> {rcos:.2f}, unit {int(unit)} at "
                         f"+{lag / 30:.2f} ms", fontsize=7)
    fig.suptitle("Unclaimed events below the cutoff losing 10 % with a sorted spike 0.5 to 1.5 ms after: "
                 "trace (black), subtracted sorted spikes (red), residual (green); µV, wires side by side",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(pathlib.Path(out_dir) / "overlaps.png", dpi=100)
    plt.close(fig)


def main():
    gate_dir, out_dir = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
    workers = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else 12
    out_dir.mkdir(parents=True, exist_ok=True)
    tet = pd.read_csv(RESULTS / "tetrodes.csv").set_index("tetrode")
    cuts = tet["cut_10"].to_dict()
    todo = [w for w in range(96) if not (out_dir / f"w{w:02d}.npz").exists()]
    with ProcessPoolExecutor(workers, initializer=init) as ex:
        for w, n, dt in ex.map(run_window, todo, [gate_dir] * len(todo), [out_dir] * len(todo),
                               [cuts] * len(todo)):
            print(f"[{time.strftime('%H:%M:%S')}] window {w}: {n} events, {dt:.0f} s", flush=True)
    summarize(gate_dir, out_dir, tet)
    if "--no-figure" not in sys.argv:
        figure(out_dir)


if __name__ == "__main__":
    main()
