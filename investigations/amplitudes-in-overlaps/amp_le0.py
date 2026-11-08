"""Where the sorted units' spike amplitudes of 0 or less (sign-flipped so a unit's median is positive) come from.

SpikeInterface's amplitude is one sample: the trace at the spike frame plus the unit's template peak shift, on
the unit's main channel. Per unit: counts, the template peak and the main channel's noise, where the spikes sit
(bad intervals, pulses, other spikes of the tetrode, windows). For up to 100 of them per unit and 100 of the
unit's other spikes: carry scores (percentile within the unit's carry fits of the same window) and the tetrode
trace around them against the unit's template."""
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.core.template_tools import get_template_peak_shift_on_main_channel
from ttm import carry as cr

OUT = sys.argv[1]
TD = "/nvme/neuropixels/tetrode_data"
FS = 30000.0
WIN = int(1800 * FS)
S = {"TTM-003-48H": "2026-10-20_09-02-12", "TTM-003-48H-2": "2026-11-03_09-01-38"}
NS = 100
HALF = 30
rng = np.random.default_rng(0)


def near(t, ref, d):
    if ref.size == 0:
        return np.zeros(t.size, bool)
    return np.searchsorted(ref, t + d, "right") > np.searchsorted(ref, t - d, "left")


def in_iv(t, st, en):
    if st.size == 0:
        return np.zeros(t.size, bool)
    o = np.argsort(st)
    st, en = st[o], en[o]
    i = np.searchsorted(st, t, "right") - 1
    return (i >= 0) & (t < en[np.clip(i, 0, None)])


units_rows, snip_rows = [], []
for exp, ses in S.items():
    run = f"{TD}/{ses}/sortings_carry/main"
    a = si.load_sorting_analyzer(f"{run}/analyzer.zarr")
    s = a.sorting
    rec = a.recording
    chans = list(rec.channel_ids)
    ids = s.unit_ids
    mua = np.asarray(s.get_property("is_mua"))
    grp = np.asarray(s.get_property("group")).astype(int)
    main_id = np.asarray(s.get_property("main_channel_id"))
    rgrp = np.asarray(rec.get_property("group")).astype(int)
    sv = s.to_spike_vector()
    amps = a.get_extension("spike_amplitudes").get_data()
    text = a.get_extension("templates")
    T = text.get_templates(operator="average")
    nbefore = text.nbefore
    shifts = get_template_peak_shift_on_main_channel(a, with_dict=True)
    noise = a.get_extension("noise_levels").get_data()
    bad = cr.read_bad_periods(f"{TD}/{ses}/bad_periods/bad_periods.csv", FS)
    p = pd.read_csv(f"subjects/ttm-003/results/laser_pulses_{ses}.csv")
    pon = np.round(p.onset_s.to_numpy() * FS).astype(np.int64)
    poff = np.round(p.offset_s.to_numpy() * FS).astype(np.int64) + int(0.010 * FS)
    order = np.argsort(sv["unit_index"], kind="stable")
    bounds = np.searchsorted(sv["unit_index"][order], np.arange(ids.size + 1))
    fbu = {k: sv["sample_index"][order[bounds[k]:bounds[k + 1]]] for k in range(ids.size)}
    abu = {k: amps[order[bounds[k]:bounds[k + 1]]] for k in range(ids.size)}
    sorted_by_group, mua_by_group = {}, {}
    for k in range(ids.size):
        (mua_by_group if mua[k] else sorted_by_group).setdefault(grp[k], []).append(k)
    requests = []
    for k in np.flatnonzero(~mua):
        u, g = int(ids[k]), int(grp[k])
        fr = fbu[k]
        x = abu[k].astype(np.float64)
        sign = -1.0 if np.median(x) < 0 else 1.0
        x = sign * x
        le0 = x <= 0
        mc = chans.index(main_id[k])
        tpk = sign * T[k, nbefore + shifts[ids[k]], mc]
        ext_ch = chans[int(np.argmax(-sign * -T[k].min(axis=0) if sign < 0 else T[k].max(axis=0)))]
        sd = noise[mc]
        bl = np.array(bad.get(g, []), np.int64).reshape(-1, 2)
        in_bad = in_iv(fr, bl[:, 0], bl[:, 1])
        in_pulse = in_iv(fr, pon, poff)
        others = np.sort(np.concatenate([fbu[j] for j in sorted_by_group[g] if j != k] or [np.array([], np.int64)]))
        muaf = np.sort(np.concatenate([fbu[j] for j in mua_by_group.get(g, [])] or [np.array([], np.int64)]))
        nb_sorted = near(fr, others, 30)
        nb_mua = near(fr, muaf, 30)
        d_prev, d_next = np.diff(np.r_[-10**12, fr]), np.diff(np.r_[fr, 10**13])
        nb_self = (d_prev <= 90) | (d_next <= 90)
        win = np.minimum(fr // WIN, 95).astype(int)
        wc = np.bincount(win[le0], minlength=96)
        wn = np.bincount(win, minlength=96)
        r = {"exp": exp, "unit": u, "tetrode": g + 1, "main_ch": main_id[k], "template_extremum_ch": ext_ch,
             "spikes": fr.size, "peak_shift": int(shifts[ids[k]]), "n_le0": int(le0.sum()), "frac_le0": le0.mean(), "median_amp": np.median(x),
             "tmpl_peak": tpk, "noise_main": sd, "median_over_noise": np.median(x) / sd,
             "top_window_share": wc.max() / max(le0.sum(), 1), "windows_with_le0": int((wc > 0).sum()),
             "max_window_frac": (wc / np.maximum(wn, 1)).max()}
        for name, m in (("bad", in_bad), ("pulse", in_pulse), ("sorted_1ms", nb_sorted), ("mua_1ms", nb_mua),
                        ("self_3ms", nb_self)):
            r[f"{name}_le0"] = m[le0].mean() if le0.any() else np.nan
            r[f"{name}_all"] = m.mean()
        units_rows.append(r)
        i_le = np.flatnonzero(le0)
        i_le = rng.choice(i_le, min(NS, i_le.size), replace=False) if i_le.size else i_le
        i_ct = rng.choice(np.flatnonzero(~le0), min(NS, int((~le0).sum())), replace=False)
        for kind, idx in (("le0", i_le), ("control", i_ct)):
            for i in idx:
                requests.append({"exp": exp, "k": k, "unit": u, "g": g, "sign": sign, "mc": mc, "kind": kind,
                                 "frame": int(fr[i]), "w": int(win[i]), "amp": x[i], "in_bad": bool(in_bad[i]),
                                 "in_pulse": bool(in_pulse[i]), "nb_sorted": bool(nb_sorted[i]),
                                 "nb_mua": bool(nb_mua[i]), "nb_self": bool(nb_self[i])})
    R = pd.DataFrame(requests)
    for w, rw in R.groupby("w"):
        z = np.load(f"{run}/carry/spikes_w{w:02d}.npz")
        cu, cf = z["unit"], z["frame"]
        o = np.lexsort((cf, cu))
        cu, cf = cu[o], cf[o]
        cols = {c: z[c][o] for c in ("a", "snr", "cos_raw", "cos_pl", "trough", "reason")}
        for _, q in rw.iterrows():
            row = {c: q[c] for c in ("exp", "unit", "kind", "frame", "w", "amp", "in_bad", "in_pulse", "nb_sorted",
                                     "nb_mua", "nb_self")}
            lo, hi = np.searchsorted(cu, q.unit, "left"), np.searchsorted(cu, q.unit, "right")
            j = lo + np.searchsorted(cf[lo:hi], q.frame)
            if j < hi and cf[j] == q.frame:
                acc = cols["reason"][lo:hi] == 0
                for c in ("a", "snr", "cos_raw", "cos_pl", "trough"):
                    v = cols[c][j]
                    row[c] = v
                    row[c + "_pct"] = (cols[c][lo:hi][acc] < v).mean()
                row["reason"] = int(cols["reason"][j])
            k = int(q.k)
            sh = shifts[ids[k]]
            tet = np.flatnonzero(rgrp == q.g)
            c0 = q.frame + sh
            tr = q.sign * rec.get_traces(start_frame=c0 - HALF, end_frame=c0 + HALF + 1,
                                         channel_ids=[chans[t] for t in tet], return_in_uV=True).astype(np.float64)
            Tp = np.pad(T[k][:, tet].astype(np.float64), ((2 * HALF, 2 * HALF), (0, 0)))
            tm = q.sign * Tp[2 * HALF + nbefore + sh - HALF:2 * HALF + nbefore + sh + HALF + 1]
            mcl = list(tet).index(q.mc)
            row["amp_check"] = tr[HALF, mcl]
            row["main_max_05ms"] = tr[HALF - 15:HALF + 16, mcl].max()
            row["main_argmax_05ms"] = int(tr[HALF - 15:HALF + 16, mcl].argmax()) - 15
            row["tet_max_center"] = tr[HALF].max()
            row["tet_argmax_center_is_main"] = bool(tr[HALF].argmax() == mcl)
            tt = tm[15:46]
            best = (-2.0, 0, 0.0)
            for lag in range(-10, 11):
                seg = tr[15 + lag:46 + lag]
                cv = (seg * tt).sum() / max(np.sqrt((seg ** 2).sum() * (tt ** 2).sum()), 1e-9)
                if cv > best[0]:
                    best = (cv, lag, (seg * tt).sum() / max((tt ** 2).sum(), 1e-9))
            seg0 = tr[15:46]
            row["cos0"] = (seg0 * tt).sum() / max(np.sqrt((seg0 ** 2).sum() * (tt ** 2).sum()), 1e-9)
            row["scale0"] = (seg0 * tt).sum() / max((tt ** 2).sum(), 1e-9)
            row["cos_best"], row["lag_best"], row["scale_best"] = best
            row["tmpl_peak"] = tm[HALF, mcl]
            snip_rows.append(row)
    print(exp, "done", flush=True)

U = pd.DataFrame(units_rows)
N = pd.DataFrame(snip_rows)
U.to_csv(f"{OUT}/amp_le0_units.csv", index=False)
N.to_csv(f"{OUT}/amp_le0_snippets.csv", index=False)
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 200)
print(U.round(4).to_string(index=False))
