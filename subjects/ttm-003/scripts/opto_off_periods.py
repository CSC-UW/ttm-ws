"""Pulse-locked silences, slow waves and multi-unit responses at every laser train, against sham trains.

For each delivered train of an experiment (``results/stimulation_blocks_<session>.csv``) this reads the store
from 72 s before the first pulse to 11 s after the last and computes:

- the light-edge artifact check of ``opto_blocks.py``: the 300 to 6000 Hz pulse-triggered average per channel,
  without a reference, its largest value from -0.2 to 0.7 ms around onset and offset, and its noise;
- multi-unit crossings: 300 to 6000 Hz, the median of every good channel subtracted, negative crossings of
  5 MAD per wire (MAD from the 60 s before the train), merged per tetrode within 1 ms;
- the LFP: 0.5 to 30 Hz, the mean of each tetrode's wires, at 500 Hz;
- events: every pulse onset, and a sham train of up to 40 onsets at the train's pulse period, the last one
  2 s before the first pulse (no light, the same state); for each event the crossings of every tetrode and
  the LFP from -0.5 to +0.95 s;
- the 60 s before the train: each tetrode's crossings per second and the delta (0.5 to 4 Hz) share of its
  LFP power (0.5 to 30 Hz).

Each train is cached in ``<experiment>/opto_off_periods/train<NNN>.npz``. The summary then:

- classes a tetrode as suppressed when its crossing rate from 10 ms after onset to pulse end, over all
  pulses, is under half its rate in the 500 ms before onset; the population is the suppressed tetrodes
  pooled (or, with none, every tetrode whose pre-train rate is 5/s or more);
- calls an event an OFF period when the population has a gap between crossings of 50 ms or more that
  spans onset + 20 ms, and records where that gap ends (the OFF period's end, from onset);
- measures each event's slow wave as the mean over the population's tetrodes of the LFP from 100 to 300 ms
  after onset minus its mean in the 200 ms before onset, and the same from 10 to 100 ms (the early wave);
- tabulates per train the OFF fraction and median OFF end of the pulses and of the sham onsets, both LFP
  waves of each, and the pre-train state, flags the trains inside the sleep deprivation window, and draws
  the figures.

Writes ``results/opto_off_periods_<E>.csv`` (one row per train), ``results/opto_off_periods_<E>.json``
(tetrode classes and the sleep-deprivation trend) and figures in ``<experiment>/figures/opto_off_periods/``.

Usage: python subjects/ttm-003/scripts/opto_off_periods.py --subject TTM-003 --experiment TTM-003-48H-2
           [--workers 12] [--summary-only]
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
import numpy as np
import pandas as pd
import scipy.signal as sps
import scipy.stats as st
import spikeinterface as si
from joblib import Parallel, delayed

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ttm.cli._common import add_experiment_args, experiment_from_args
from ttm.experiment import ExperimentParams, dt2t, read_session_t0_unix

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
LFP_FS = 500
PRE_S, POST_S = 0.5, 0.95
THRESH_MAD = 5.0
DEAD_MS = 1.0
N_SHAM = 40
OFF_MIN_S = 0.05
PROBE_S = 0.02
SW_WIN = (0.1, 0.3)
EARLY_WIN = (0.01, 0.1)
N_BOOT = 5000


def crossings(xf: np.ndarray, mad: np.ndarray, groups: np.ndarray) -> dict[int, np.ndarray]:
    dead = int(DEAD_MS * FS / 1000)
    out = {}
    for g in np.unique(groups):
        idx = []
        for c in np.flatnonzero(groups == g):
            below = xf[:, c] < -THRESH_MAD * mad[c]
            idx.append(np.flatnonzero(below[1:] & ~below[:-1]) + 1)
        idx = np.sort(np.concatenate(idx))
        if idx.size:
            idx = idx[np.r_[True, np.diff(idx) > dead]]
        out[int(g)] = idx
    return out


def one_train(store: str, drop: list[str], blk: pd.Series, on_s: np.ndarray, off_s: np.ndarray, out: pathlib.Path):
    if out.exists():
        return
    rec = si.read_zarr(store)
    if drop:
        rec = rec.remove_channels(drop)
    groups = np.asarray(rec.get_property("group")).astype(int)
    tets = np.unique(groups)
    period = float(np.median(np.diff(on_s))) if on_s.size > 1 else 1.0
    a = max(int((on_s[0] - 72) * FS), 0)
    b = min(int((on_s[-1] + 11) * FS), rec.get_num_frames())
    x = rec.get_traces(start_frame=a, end_frame=b, return_in_uV=True).astype(np.float32)
    t0 = a / FS
    hf = sps.sosfiltfilt(sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos"), x, axis=0).astype(np.float32)
    lf = sps.sosfiltfilt(sps.butter(3, [0.5, 30], btype="bandpass", fs=FS, output="sos"), x, axis=0)
    del x
    lf = lf[::FS // LFP_FS].astype(np.float32)
    lfp = np.stack([lf[:, groups == g].mean(1) for g in tets], axis=1)
    del lf
    on = np.round((on_s - t0) * FS).astype(np.int64)
    off = np.round((off_s - t0) * FS).astype(np.int64)

    hp, hq = int(0.005 * FS), int(0.012 * FS)
    ok = (on - hp >= 0) & (off + hq <= hf.shape[0])
    hf_on = np.mean([hf[i - hp:i + hq] for i in on[ok]], axis=0)
    hf_off = np.mean([hf[i - hp:i + hq] for i in off[ok]], axis=0)
    ht = (np.arange(hp + hq) - hp) / FS * 1000
    edge = (ht >= -0.2) & (ht < 0.7)

    ref = hf - np.median(hf, axis=1, keepdims=True)
    del hf
    pre_a = max(on[0] - 60 * FS, 0)
    mad = np.median(np.abs(ref[pre_a:on[0]]), axis=0) / 0.6745
    mua = crossings(ref, mad, groups)
    del ref

    last_sham = on_s[0] - 2.0
    n_sham = int(min(N_SHAM, max(0, np.floor((last_sham - t0 - PRE_S - 1) / period) + 1)))
    sham_s = last_sham - period * np.arange(n_sham)[::-1]
    ev_s = np.r_[on_s, sham_s]
    is_sham = np.r_[np.zeros(on_s.size, bool), np.ones(n_sham, bool)]

    rel_t, rel_ev, rel_tet = [], [], []
    for j, g in enumerate(tets):
        spk = mua[int(g)] / FS + t0
        lo = np.searchsorted(spk, ev_s - PRE_S)
        hi = np.searchsorted(spk, ev_s + POST_S)
        for e in range(ev_s.size):
            seg = spk[lo[e]:hi[e]] - ev_s[e]
            rel_t.append(seg.astype(np.float32))
            rel_ev.append(np.full(seg.size, e, np.int32))
            rel_tet.append(np.full(seg.size, j, np.int16))
    n_l = int((PRE_S + POST_S) * LFP_FS)
    li = np.round((ev_s - t0 - PRE_S) * LFP_FS).astype(np.int64)
    lfp_ev = np.stack([lfp[i:i + n_l] for i in li]).astype(np.float32)

    pre_l = lfp[max(0, int((on_s[0] - t0 - 60) * LFP_FS)):int((on_s[0] - t0) * LFP_FS)]
    f, pxx = sps.welch(pre_l, fs=LFP_FS, nperseg=4 * LFP_FS, axis=0)
    band = (f >= 0.5) & (f <= 30)
    delta_share = pxx[(f >= 0.5) & (f <= 4)].sum(0) / pxx[band].sum(0)
    pre_dur = (on[0] - pre_a) / FS
    pre_rate = np.array([((mua[int(g)] >= pre_a) & (mua[int(g)] < on[0])).sum() / pre_dur for g in tets])
    np.savez_compressed(out, ttl_train=int(blk.ttl_train), onset_h=float(blk.onset_h), ev_s=ev_s, is_sham=is_sham,
                        pulse_s=float(np.median(off_s - on_s)), tets=tets, rel_t=np.concatenate(rel_t),
                        rel_ev=np.concatenate(rel_ev), rel_tet=np.concatenate(rel_tet), lfp=lfp_ev,
                        edge_on=np.abs(hf_on[edge]).max(0), edge_off=np.abs(hf_off[edge]).max(0),
                        hf_noise=hf_on[ht < -1].std(0), delta_share=delta_share, pre_rate=pre_rate,
                        mad=mad, groups=groups)


def off_period(times: np.ndarray, probe: float = PROBE_S) -> tuple[bool, float, float]:
    """Whether a gap of OFF_MIN_S or more between pooled crossings spans ``probe``; its end and length."""
    before = times[times <= probe]
    after = times[times > probe]
    start = before.max() if before.size else -PRE_S
    end = after.min() if after.size else POST_S
    return (end - start) >= OFF_MIN_S, float(end), float(end - start)


def light_dark_state(p, df: pd.DataFrame, sd) -> dict:
    """Trains grouped by light phase and day (the manifest's light times), and a least-squares fit of each
    pulse effect on session hour, delta share before the train, dark and sleep deprivation."""
    tr = sorted([(dt2t(p, x) / 3600, "light") for x in p.lightsOn] + [(dt2t(p, x) / 3600, "dark") for x in p.lightsOff])
    lab, day = [], []
    for h in df.onset_h:
        prev = [t for t in tr if t[0] <= h]
        lab.append(prev[-1][1] if prev else ("dark" if tr[0][1] == "light" else "light"))
        day.append(sum(1 for t in tr if t[1] == "light" and t[0] <= h))
    d = df.assign(phase=lab, day=day)
    d["group"] = np.where(d.in_sleep_deprivation, "sleep deprivation", d.phase + ", day " + d.day.astype(str))
    cols = ["off_frac_pulse", "off_frac_sham", "early_wave_uV_pulse", "slow_wave_uV_pulse", "pre_delta_share"]
    g = d.groupby("group", sort=False).agg(trains=("ttl_train", "size"), first_h=("onset_h", "min"),
                                           last_h=("onset_h", "max"), **{c: (c, "mean") for c in cols})
    X = np.column_stack([np.ones(len(d)), d.onset_h, d.pre_delta_share, (d.phase == "dark").astype(float),
                         d.in_sleep_deprivation.astype(float)])
    fits = {}
    for y in ("early_wave_uV_pulse", "off_frac_pulse", "slow_wave_uV_pulse"):
        v = d[y].to_numpy(float)
        b = np.linalg.lstsq(X, v, rcond=None)[0]
        r = v - X @ b
        se = np.sqrt(np.diag(r @ r / (len(v) - X.shape[1]) * np.linalg.inv(X.T @ X)))
        fits[y] = {n: {"coef": float(c), "ci95": [float(c - 1.96 * e), float(c + 1.96 * e)]}
                   for n, c, e in zip(("intercept", "per_hour", "delta_share", "dark", "sleep_deprivation"), b, se)}
    return {"light_transitions_h": [[round(t, 4), k] for t, k in tr], "groups": g.reset_index().to_dict(orient="records"),
            "fits": fits}


def summarize(exp, cache: pathlib.Path, figdir: pathlib.Path) -> None:
    files = sorted(cache.glob("train*.npz"))
    Z = [dict(np.load(f)) for f in files]
    tets = Z[0]["tets"]
    names = [f"TT{int(g) + 1}" for g in tets]
    bins = np.arange(-PRE_S, POST_S + 1e-9, 0.005)
    centers = bins[:-1] + 0.0025
    H = {k: np.zeros((len(tets), bins.size - 1)) for k in ("real", "sham")}
    n_ev = {"real": 0, "sham": 0}
    in_pulse = np.zeros(len(tets))
    early = np.zeros(len(tets))
    base = np.zeros(len(tets))
    from15, post_a, post_b = np.zeros(len(tets)), np.zeros(len(tets)), np.zeros(len(tets))
    train_from15 = np.full((len(Z), len(tets)), np.nan)
    for i_z, z in enumerate(Z):
        ps = float(z["pulse_s"])
        for kind, m in (("real", ~z["is_sham"]), ("sham", z["is_sham"])):
            evs = np.flatnonzero(m)
            n_ev[kind] += evs.size
            sel = np.isin(z["rel_ev"], evs)
            for j in range(len(tets)):
                tj = z["rel_t"][sel & (z["rel_tet"] == j)]
                H[kind][j] += np.histogram(tj, bins)[0]
                if kind == "real":
                    in_pulse[j] += ((tj >= 0.01) & (tj < ps)).sum() / (ps - 0.01)
                    early[j] += ((tj >= 0.001) & (tj < 0.01)).sum() / 0.009
                    base[j] += ((tj >= -PRE_S) & (tj < 0)).sum() / PRE_S
                    n15 = ((tj >= 0.015) & (tj < ps)).sum() / (ps - 0.015)
                    from15[j] += n15
                    post_a[j] += ((tj >= ps) & (tj < ps + 0.05)).sum() / 0.05
                    post_b[j] += ((tj >= ps + 0.05) & (tj < ps + 0.1)).sum() / 0.05
                    b_j = ((tj >= -PRE_S) & (tj < 0)).sum() / PRE_S
                    if b_j > 0:
                        train_from15[i_z, j] = n15 / b_j
    ratio = in_pulse / np.maximum(base, 1e-9)
    early_ratio = early / np.maximum(base, 1e-9)
    rate_before = base / max(n_ev["real"], 1)
    ps_all = float(np.median([z["pulse_s"] for z in Z]))
    after = centers >= ps_all + 0.1
    rebound = (H["real"][:, after] / (max(n_ev["real"], 1) * 0.005)).max(1) / np.maximum(rate_before, 1e-9)
    pre_rate = np.median(np.stack([z["pre_rate"] for z in Z]), axis=0)
    suppressed = (ratio < 0.5) & (pre_rate >= 5)
    pop = suppressed if suppressed.any() else pre_rate >= 5
    tet_table = pd.DataFrame({"tetrode": names, "pre_rate_hz": pre_rate.round(2), "pulse_over_baseline": ratio.round(3),
                              "onset_1_10ms_over_baseline": early_ratio.round(3),
                              "suppressed": suppressed, "in_population": pop,
                              "rate_500ms_before_hz": rate_before.round(2),
                              "pulse_15ms_off_over_baseline": (from15 / np.maximum(base, 1e-9)).round(3),
                              "post_0_50ms_over_baseline": (post_a / np.maximum(base, 1e-9)).round(3),
                              "post_50_100ms_over_baseline": (post_b / np.maximum(base, 1e-9)).round(3),
                              "rebound_max_5ms_over_baseline": rebound.round(3),
                              "onset_1_10ms_extra_per_pulse": ((early_ratio - 1) * rate_before * 0.009).round(3),
                              "trains_15ms_off_over_1p2": (train_from15 > 1.2).sum(0),
                              "train_15ms_off_ratio_median": np.nanmedian(train_from15, axis=0).round(3)})
    edge_ratio = np.stack([np.maximum(z["edge_on"], z["edge_off"]) / z["hf_noise"] for z in Z])
    ch_ids = [str(c) for c in si.read_zarr(str(exp.store_path())).channel_ids
              if str(c) not in {str(b) for b in exp.badChannels}]

    p = ExperimentParams.from_experiment(exp, read_session_t0_unix(exp.slice_table_path))
    sd = None
    if exp.novel_objects_start and exp.novel_objects_end:
        sd = (dt2t(p, exp.novel_objects_start) / 3600, dt2t(p, exp.novel_objects_end) / 3600)

    rows, per_event = [], []
    lfp_t = np.arange(int((PRE_S + POST_S) * LFP_FS)) / LFP_FS - PRE_S
    sw_m = (lfp_t >= SW_WIN[0]) & (lfp_t < SW_WIN[1])
    bl_m = (lfp_t >= -0.2) & (lfp_t < 0)
    early_m = (lfp_t >= EARLY_WIN[0]) & (lfp_t < EARLY_WIN[1])
    L = {"real": [], "sham": []}
    for z in Z:
        popj = np.flatnonzero(pop)
        pooled_sel = np.isin(z["rel_tet"], popj)
        lf = z["lfp"][:, :, popj].mean(2)
        sw = lf[:, sw_m].mean(1) - lf[:, bl_m].mean(1)
        early = lf[:, early_m].mean(1) - lf[:, bl_m].mean(1)
        res = {}
        for kind, m in (("real", ~z["is_sham"]), ("sham", z["is_sham"])):
            evs = np.flatnonzero(m)
            offs = [off_period(np.sort(z["rel_t"][pooled_sel & (z["rel_ev"] == e)])) for e in evs]
            is_off = np.array([o[0] for o in offs], bool)
            ends = np.array([o[1] for o in offs])
            res[kind] = (is_off, ends, sw[evs], early[evs])
            L[kind].append(lf[evs])
            for e, o, s, w in zip(evs, offs, sw[evs], early[evs]):
                per_event.append({"ttl_train": int(z["ttl_train"]), "kind": kind, "event": int(e), "off": o[0],
                                  "off_end_s": o[1], "gap_s": o[2], "slow_wave_uV": float(s), "early_wave_uV": float(w)})
        h = float(z["onset_h"])
        r_off, r_end, r_sw, r_ew = res["real"]
        s_off, s_end, s_sw, s_ew = res["sham"]
        rows.append({
            "ttl_train": int(z["ttl_train"]), "onset_h": round(h, 4),
            "in_sleep_deprivation": bool(sd is not None and sd[0] <= h < sd[1]),
            "pulse_ms": round(float(z["pulse_s"]) * 1000, 1), "n_pulses": int((~z["is_sham"]).sum()),
            "n_sham": int(z["is_sham"].sum()),
            "off_frac_pulse": float(r_off.mean()), "off_frac_sham": float(s_off.mean()) if s_off.size else np.nan,
            "off_end_ms_median_pulse": float(np.median(r_end[r_off]) * 1000) if r_off.any() else np.nan,
            "off_end_ms_median_sham": float(np.median(s_end[s_off]) * 1000) if s_off.any() else np.nan,
            "slow_wave_uV_pulse": float(r_sw.mean()), "slow_wave_uV_sham": float(s_sw.mean()) if s_sw.size else np.nan,
            "early_wave_uV_pulse": float(r_ew.mean()),
            "early_wave_uV_sham": float(s_ew.mean()) if s_ew.size else np.nan,
            "pre_delta_share": float(z["delta_share"][pop].mean()),
            "pre_rate_hz": float(z["pre_rate"][pop].sum()),
            "edge_on_uV_max": float(z["edge_on"].max()), "edge_off_uV_max": float(z["edge_off"].max()),
            "hf_noise_uV_median": float(np.median(z["hf_noise"])),
        })
    df = pd.DataFrame(rows).sort_values("onset_h")
    df["off_frac_excess"] = df.off_frac_pulse - df.off_frac_sham
    df["slow_wave_excess_uV"] = df.slow_wave_uV_pulse - df.slow_wave_uV_sham
    df["early_wave_excess_uV"] = df.early_wave_uV_pulse - df.early_wave_uV_sham
    df.round(4).to_csv(RESULTS / f"opto_off_periods_{exp.experiment}.csv", index=False)
    ev = pd.DataFrame(per_event)
    ev.to_csv(cache / "events.csv", index=False)

    trend = {}
    sdd = df[df.in_sleep_deprivation]
    if len(sdd) >= 4:
        for col in ("off_frac_pulse", "off_frac_sham", "off_frac_excess", "off_end_ms_median_pulse",
                    "slow_wave_uV_pulse", "slow_wave_excess_uV", "early_wave_uV_pulse", "early_wave_excess_uV",
                    "pre_delta_share", "pre_rate_hz"):
            v = sdd[col].to_numpy(float)
            okv = np.isfinite(v)
            rho, pv = st.spearmanr(sdd.onset_h.to_numpy()[okv], v[okv]) if okv.sum() >= 4 else (np.nan, np.nan)
            k = max(1, len(sdd) // 3)
            trend[col] = {"spearman_rho": float(rho), "p": float(pv), "first_third": float(np.nanmean(v[:k])),
                          "last_third": float(np.nanmean(v[-k:]))}
            x, idx = sdd.onset_h.to_numpy(), np.flatnonzero(okv)
            rng = np.random.default_rng(0)
            bs = [np.polyfit(x[b], v[b], 1)[0] for b in (rng.choice(idx, idx.size) for _ in range(N_BOOT)) if np.ptp(x[b]) > 0]
            span = sd[1] - sd[0]
            trend[col]["change_over_window"] = float(np.polyfit(x[idx], v[idx], 1)[0] * span)
            trend[col]["change_over_window_ci95"] = [float(q * span) for q in np.percentile(bs, [2.5, 97.5])]
        k = max(1, len(sdd) // 3)
        for name, trs in (("first", sdd.ttl_train.iloc[:k]), ("last", sdd.ttl_train.iloc[-k:])):
            e = ev[ev.ttl_train.isin(trs) & (ev.kind == "real")]
            ci = st.binomtest(int(e.off.sum()), len(e)).proportion_ci()
            trend[f"off_pulse_{name}_trains"] = {"trains": [int(t) for t in trs], "off": int(e.off.sum()), "pulses": len(e),
                                                 "ci95": [float(ci.low), float(ci.high)]}
        # pulse vs sham inside each third, events pooled
        thirds = np.array_split(sdd.ttl_train.to_numpy(), 3)
        trend["thirds_pooled"] = []
        for i, trs in enumerate(thirds):
            e = ev[ev.ttl_train.isin(trs)]
            rp, rs = e[e.kind == "real"], e[e.kind == "sham"]
            tab = [[int(rp.off.sum()), int((~rp.off).sum())], [int(rs.off.sum()), int((~rs.off).sum())]]
            trend["thirds_pooled"].append({
                "trains": [int(t) for t in trs], "off_frac_pulse": float(rp.off.mean()), "off_frac_sham": float(rs.off.mean()),
                "fisher_p": float(st.fisher_exact(tab)[1]),
                "off_end_ms_median_pulse": float(rp.off_end_s[rp.off].median() * 1000) if rp.off.any() else None,
                "slow_wave_uV_pulse": float(rp.slow_wave_uV.mean()), "slow_wave_uV_sham": float(rs.slow_wave_uV.mean()),
                "early_wave_uV_pulse": float(rp.early_wave_uV.mean()),
                "early_wave_uV_sham": float(rs.early_wave_uV.mean())})
    allp = ev[ev.kind == "real"]
    alls = ev[ev.kind == "sham"]
    state = light_dark_state(p, df, sd)
    lm = np.concatenate(L["real"]).mean(0)
    lm = lm - lm[bl_m].mean()
    w1, w2 = np.flatnonzero((lfp_t >= 0) & (lfp_t < 0.15)), np.flatnonzero((lfp_t >= 0.1) & (lfp_t < 0.4))
    i_pk, i_tr = w1[lm[w1].argmax()], w2[lm[w2].argmin()]
    i_zero = i_pk + int(np.argmax(lm[i_pk:] <= 0))
    i_back = i_tr + int(np.argmax(lm[i_tr:] >= 0))
    lfp_shape = {"peak_uV": float(lm[i_pk]), "peak_ms": float(lfp_t[i_pk] * 1000), "zero_ms": float(lfp_t[i_zero] * 1000),
                 "trough_uV": float(lm[i_tr]), "trough_ms": float(lfp_t[i_tr] * 1000), "back_ms": float(lfp_t[i_back] * 1000)}
    by_off = allp.groupby("off")[["early_wave_uV", "slow_wave_uV"]].mean()
    summary = {
        "experiment": exp.experiment, "trains": len(df), "pulses": len(allp), "sham_onsets": len(alls),
        "sleep_deprivation_h": sd, "tetrodes": tet_table.to_dict(orient="records"),
        "population": [n for n, q in zip(names, pop) if q],
        "off_frac_pulse": float(allp.off.mean()), "off_frac_sham": float(alls.off.mean()),
        "off_end_ms_median_pulse": float(allp.off_end_s[allp.off].median() * 1000),
        "off_end_ms_quartiles_pulse": [float(q * 1000) for q in allp.off_end_s[allp.off].quantile([0.25, 0.75])],
        "slow_wave_uV_pulse": float(allp.slow_wave_uV.mean()), "slow_wave_uV_sham": float(alls.slow_wave_uV.mean()),
        "early_wave_uV_pulse": float(allp.early_wave_uV.mean()), "early_wave_uV_sham": float(alls.early_wave_uV.mean()),
        "edge_on_uV_max": float(df.edge_on_uV_max.max()), "edge_off_uV_max": float(df.edge_off_uV_max.max()),
        "hf_noise_uV_median": float(df.hf_noise_uV_median.median()),
        "edge_over_noise_median_train_max_channel": float(np.median(edge_ratio.max(1))),
        "edge_over_noise_by_channel_median": dict(zip(ch_ids,
                                                      np.round(np.median(edge_ratio, axis=0), 2).tolist())),
        "sleep_deprivation_trend": trend,
        "lfp_mean_pulse": lfp_shape,
        "waves_by_off": {("off" if k else "no_off"): {c: float(v) for c, v in r.items()} for k, r in by_off.iterrows()},
        "state": state,
    }
    (RESULTS / f"opto_off_periods_{exp.experiment}.json").write_text(json.dumps(summary, indent=1) + "\n")

    # Figures: PSTH per tetrode, evoked LFP, per-train metrics against time.
    figdir.mkdir(parents=True, exist_ok=True)
    ps_med = float(df.pulse_ms.median()) / 1000
    nt = len(tets)
    fig, axes = plt.subplots(int(np.ceil(nt / 4)), 4, figsize=(18, 3 * int(np.ceil(nt / 4))), sharex=True)
    for j, ax in enumerate(axes.ravel()):
        if j >= nt:
            ax.axis("off")
            continue
        for kind, col in (("real", "k"), ("sham", "0.6")):
            ax.plot(centers * 1000, H[kind][j] / (max(n_ev[kind], 1) * 0.005), color=col, lw=0.8, label=kind)
        ax.axvspan(0, ps_med * 1000, color="c", alpha=0.2)
        ax.set_title(f"{names[j]}: in-pulse/baseline {ratio[j]:.2f}{' (population)' if pop[j] else ''}", fontsize=9)
    axes.ravel()[0].legend(fontsize=7)
    fig.supxlabel("ms from onset"); fig.supylabel("crossings/s")
    fig.suptitle(f"{exp.experiment}: multi-unit crossings around {n_ev['real']} pulses (black) and {n_ev['sham']} sham onsets (grey), 5 ms bins")
    fig.tight_layout(); fig.savefig(figdir / "psth_tetrodes.png", dpi=100); plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    lr, ls = np.concatenate(L["real"]), np.concatenate(L["sham"])
    ax = axes[0]
    ax.plot(lfp_t * 1000, lr.mean(0), "k", label=f"pulses (n={len(lr)})")
    ax.plot(lfp_t * 1000, ls.mean(0), color="0.6", label=f"sham (n={len(ls)})")
    ax.axvspan(0, ps_med * 1000, color="c", alpha=0.2)
    ax.axvspan(SW_WIN[0] * 1000, SW_WIN[1] * 1000, color="y", alpha=0.15)
    ax.axvspan(EARLY_WIN[0] * 1000, EARLY_WIN[1] * 1000, ymin=0.95, color="m", alpha=0.4)
    ax.set_xlabel("ms from onset"); ax.set_ylabel("uV"); ax.legend(); ax.set_title("LFP 0.5-30 Hz, mean of the population tetrodes")
    ax = axes[1]
    for kind, col in (("real", "k"), ("sham", "0.6")):
        e = ev[(ev.kind == kind) & ev.off]
        ax.hist(e.off_end_s * 1000, bins=np.arange(0, POST_S * 1000, 10), histtype="step", color=col,
                density=False, label=f"{kind}: OFF in {ev[ev.kind == kind].off.mean():.0%} of events")
    ax.axvline(ps_med * 1000, color="c")
    ax.set_xlabel("end of the OFF period spanning onset + 20 ms (ms from onset)"); ax.set_ylabel("events"); ax.legend()
    fig.suptitle(exp.experiment)
    fig.tight_layout(); fig.savefig(figdir / "lfp_and_off_end.png", dpi=100); plt.close(fig)

    fig, axes = plt.subplots(5, 1, figsize=(16, 16), sharex=True)
    for ax, (cols, title) in zip(axes, [
            (("off_frac_pulse", "off_frac_sham"), "share of events with an OFF period spanning onset + 20 ms"),
            (("off_end_ms_median_pulse", "off_end_ms_median_sham"), "median end of that OFF period (ms from onset)"),
            (("early_wave_uV_pulse", "early_wave_uV_sham"),
             f"LFP {EARLY_WIN[0]*1000:.0f}-{EARLY_WIN[1]*1000:.0f} ms minus pre-onset (uV)"),
            (("slow_wave_uV_pulse", "slow_wave_uV_sham"), f"LFP {SW_WIN[0]*1000:.0f}-{SW_WIN[1]*1000:.0f} ms minus pre-onset (uV)"),
            (("pre_delta_share",), "delta share of LFP power in the 60 s before the train")]):
        for c, col in zip(cols, ("k", "0.6")):
            ax.plot(df.onset_h, df[c], "o-", ms=3, lw=0.6, label=c, color=col)
        if sd:
            ax.axvspan(*sd, color="orange", alpha=0.15)
        ax.set_title(title, fontsize=9); ax.legend(fontsize=7)
    axes[-1].set_xlabel("session hour (orange: sleep deprivation)")
    fig.suptitle(f"{exp.experiment}: per-train pulse effects against sham onsets")
    fig.tight_layout(); fig.savefig(figdir / "trains_over_time.png", dpi=100); plt.close(fig)
    print(tet_table.to_string(index=False))
    print(json.dumps({k: v for k, v in summary.items() if k not in ("tetrodes",)}, indent=1))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--trains", default=None, help="comma-separated ttl_train numbers (default: every delivered train)")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    cache = exp.experiment_dir / "opto_off_periods"
    cache.mkdir(parents=True, exist_ok=True)
    figdir = exp.figures_dir / "opto_off_periods"
    if not a.summary_only:
        blocks = pd.read_csv(RESULTS / f"stimulation_blocks_{exp.session}.csv")
        blocks = blocks[blocks.delivered & (blocks.n_pulses >= 10)].reset_index(drop=True)
        if a.trains:
            blocks = blocks[blocks.ttl_train.isin([int(t) for t in a.trains.split(",")])].reset_index(drop=True)
        pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
        pulses = pulses[~pulses.glitch]
        rec = si.read_zarr(str(exp.store_path()))
        from ttm.sorting import bad_channel_ids
        drop = [str(c) for c in bad_channel_ids(rec, exp.badChannels, exp.badTetrodes)]
        jobs = []
        for _, blk in blocks.iterrows():
            pp = pulses[pulses.train == blk.ttl_train].sort_values("onset_s")
            jobs.append(delayed(one_train)(str(exp.store_path()), drop, blk, pp.onset_s.to_numpy(),
                                           pp.offset_s.to_numpy(), cache / f"train{int(blk.ttl_train):03d}.npz"))
        print(f"{len(jobs)} trains, {a.workers} workers, bad channels left out: {drop}")
        Parallel(n_jobs=a.workers)(jobs)
    summarize(exp, cache, figdir)


if __name__ == "__main__":
    main()
