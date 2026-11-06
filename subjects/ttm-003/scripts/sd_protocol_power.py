"""Precision of the first-and-last-hour stimulation protocol, from one 48 h opto experiment.

The outcome is the OFF share: the fraction of pulses followed by an OFF period as ``opto_off_periods.py``
defines it, per pulse in its cached ``events.csv``. This script measures:

- the per-pulse variance of the OFF outcome and the variance of the train OFF share beyond the binomial
  part (between-train variance), over all trains and the sleep deprivation's;
- the OFF share of the first and last five pulses of the trains, and the change in log-odds per pulse
  position from a logistic fit, with its 95 % interval;
- the covariance of the OFF outcome between pulses of one train by lag (1-3, 4-7, 8-12 and 13-19 s),
  after removing the mean of each light phase, dark phase and the sleep deprivation, a linear trend in
  each, and the mean of each pulse position, with 95 % intervals from 300 resamplings of trains;
- the covariance of the train OFF shares between trains 20 min to 6 h apart, with the same means and
  trends removed, and the standard error of each mean product;
- for spontaneous activity of the suppressed tetrodes (the multi-unit rate in 1 s bins and the delta,
  0.5 to 4 Hz, and sigma, 10 to 15 Hz, shares of 0.5 to 30 Hz LFP power in 4 s windows stepped by 1 s),
  over the sleep deprivation and over session hours 2 to 7, trains and the 2 s before and 12 s after
  them excluded: the share of variance at periods of 5 to 1,024 s that lies above 120 s, between 10 and
  120 s and between 5 and 10 s, and the frequency of the largest power times frequency (Welch, 1,024 s
  Hann segments with half overlap inside gap-free stretches, linear detrend);
- for three models of the state covariance of the OFF outcome, the smallest last-minus-first-hour
  difference in OFF share that one sleep deprivation detects with 80 % power at two-sided p < 0.05
  (2.8 times the standard error, with the hour means' variance computed exactly from the model and
  ``sd_stim_schedule.py`` schedules), the sleep deprivations needed for 0.10 and 0.05, the model's
  covariance in the lag bins above, and how much subtracting shams widens the difference.

The models share a fast component (variance 0.007, exponential with a 3 s time constant). A adds 0.010
with a 16 s time constant (a corner at 0.01 Hz); B splits that into 0.0085 at 16 s and 0.0015 at 400 s;
C into 0.0075 at 16 s and 0.0025 at 1,200 s; D puts all 0.010 at 1,200 s. The per-pulse binomial variance
is 0.205. For slow components of variance 0.0025, 0.005 and 0.010 with time constants of 400, 1,200 and
3,600 s, the chi-square of the covariances between trains against the component's prediction, beside the
chi-square against zero.

Spontaneous activity is cached in ``<experiment>/sd_protocol_power/``. Writes
``results/sd_protocol_power_<E>.json``.

Usage: python subjects/ttm-003/scripts/sd_protocol_power.py --subject TTM-003 --experiment TTM-003-48H-2
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
import scipy.signal as sps
import statsmodels.formula.api as smf
import zarr
from joblib import Parallel, delayed
from sd_stim_schedule import draw
from ttm.cli._common import add_experiment_args, experiment_from_args

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"
FS = 30000
CHUNK = 900000
LAG_BINS = {"1_3s": (1, 3), "4_7s": (4, 7), "8_12s": (8, 12), "13_19s": (13, 19)}
TRAIN_LAGS_H = (1 / 3, 2 / 3, 1, 4 / 3, 2, 3, 4, 5, 6)
FAST = (0.007, 3.0)
MODELS = {"A": [(0.010, 16.0)], "B": [(0.0085, 16.0), (0.0015, 400.0)], "C": [(0.0075, 16.0), (0.0025, 1200.0)],
          "D": [(0.010, 1200.0)]}
SLOW_GRID = [(v, tau) for tau in (400.0, 1200.0, 3600.0) for v in (0.0025, 0.005, 0.010)]
SW2 = 0.205
DESIGNS = ((240, 2), (480, 2), (960, 5))
SEG = 1024


def kernel(model: str, n: int) -> np.ndarray:
    k = np.arange(n + 1.0)
    return sum(v * np.exp(-k / tau) for v, tau in [FAST, *MODELS[model]])


def phase_groups(tr: pd.DataFrame, transitions: list) -> list[str]:
    out = []
    for h, sd in zip(tr.onset_h, tr.in_sleep_deprivation):
        prev = [t for t in transitions if t[0] <= h]
        ph = prev[-1][1] if prev else ("dark" if transitions[0][1] == "light" else "light")
        out.append("sd" if sd else f"{ph}{sum(1 for t in transitions if t[1] == 'light' and t[0] <= h)}")
    return out


def lag_cov(d: pd.DataFrame, y: str) -> dict[str, float]:
    by = d.set_index(["ttl_train", "event"])[y]
    out = {}
    for name, (a, b) in LAG_BINS.items():
        prods = []
        for k in range(a, b + 1):
            x = by[by.index.get_level_values(1) < 20 - k]
            z = by.reindex(pd.MultiIndex.from_arrays([x.index.get_level_values(0), x.index.get_level_values(1) + k]))
            prods.append(float(np.nanmean(x.to_numpy() * z.to_numpy())))
        out[name] = float(np.mean(prods))
    return out


def within_train(pp: pd.DataFrame) -> dict:
    pp = pp.copy()
    pp["r"] = pp.off - smf.ols("off ~ C(grp) * onset_h + C(event)", pp).fit().fittedvalues
    rng = np.random.default_rng(0)
    res = {}
    for name, d in (("all_trains", pp), ("sleep_deprivation", pp[pp.in_sleep_deprivation])):
        est = lag_cov(d, "r")
        g = {t: x for t, x in d.groupby("ttl_train")}
        tt = np.array(list(g))
        bs = []
        for _ in range(300):
            pick = rng.choice(tt, tt.size)
            bs.append(lag_cov(pd.concat([g[t].assign(ttl_train=i) for i, t in enumerate(pick)]), "r"))
        res[name] = {k: {"cov": round(v, 5), "ci95": [round(float(np.percentile([b[k] for b in bs], q)), 5) for q in (2.5, 97.5)]}
                     for k, v in est.items()}
        tm = d.groupby("ttl_train").off.agg(["mean", "var", "size"])
        res[name]["per_pulse_var"] = round(float(d.off.var()), 4)
        res[name]["between_train_var"] = round(float(tm["mean"].var() - tm["var"].mean() / tm["size"].mean()), 4)
        res[name]["off_share"] = round(float(d.off.mean()), 3)
        fit = smf.logit("off ~ event", d).fit(disp=0)
        res[name]["position"] = {"off_share_pulses_1_5": round(float(d.off[d.event < 5].mean()), 3),
                                 "off_share_pulses_16_20": round(float(d.off[d.event >= 15].mean()), 3),
                                 "log_odds_per_pulse": round(float(fit.params["event"]), 4),
                                 "log_odds_per_pulse_ci95": [round(float(x), 4) for x in fit.conf_int().loc["event"]]}
    return res


def between_trains(pp: pd.DataFrame) -> dict:
    t = pp.groupby("ttl_train").agg(off=("off", "mean"), h=("onset_h", "first"), grp=("grp", "first")).sort_values("h")
    r = (t.off - smf.ols("off ~ C(grp) * h", t).fit().fittedvalues).to_numpy()
    h = t.h.to_numpy()
    dh = h[None] - h[:, None]
    rr = r[:, None] * r[None]
    out = []
    for lag in TRAIN_LAGS_H:
        p = rr[np.abs(dh - lag) < 1 / 6]
        out.append({"lag_h": round(lag, 3), "cov": round(float(p.mean()), 5), "se": round(float(p.std() / np.sqrt(p.size / 2)), 5),
                    "pairs": int(p.size)})
    lag_s = np.array([o["lag_h"] for o in out]) * 3600
    c = np.array([o["cov"] for o in out])
    se = np.array([o["se"] for o in out])
    fit = [{"variance": v, "time_constant_s": tau, "chi2": round(float((((c - v * np.exp(-lag_s / tau)) / se) ** 2).sum()), 2)}
           for v, tau in SLOW_GRID]
    return {"lags": out, "mean_over_lags": round(float(np.mean(c)), 5), "chi2_none": round(float(((c / se) ** 2).sum()), 2),
            "slow_component_chi2": fit}


def spont_chunk(store: str, a: int, cols: list[int], groups: np.ndarray, tets: list[int], gain: float):
    tr = zarr.open(store, "r")["traces_seg0"]
    lo, hi = max(a - FS, 0), a + CHUNK + FS
    x = tr.get_orthogonal_selection((slice(lo, hi), cols)).astype(np.float32) * gain
    hf = sps.sosfiltfilt(sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos"), x, axis=0)
    hf = (hf - np.median(hf, axis=1, keepdims=True))[a - lo:a - lo + CHUNK]
    mad = np.median(np.abs(hf), axis=0) / 0.6745
    counts = np.zeros(CHUNK // FS, np.int64)
    for g in tets:
        idx = []
        for c in np.flatnonzero(groups == g):
            b = hf[:, c] < -5 * mad[c]
            idx.append(np.flatnonzero(b[1:] & ~b[:-1]) + 1)
        idx = np.sort(np.concatenate(idx))
        if idx.size:
            idx = idx[np.r_[True, np.diff(idx) > 30]]
        counts += np.bincount(idx // FS, minlength=CHUNK // FS)[:CHUNK // FS]
    lf = sps.sosfiltfilt(sps.butter(3, [0.5, 30], btype="bandpass", fs=FS, output="sos"), x, axis=0)[a - lo:a - lo + CHUNK:60]
    return a, counts, lf.mean(1).astype(np.float32)


def spontaneous(exp, h0: float, h1: float, tets: list[int], cache: pathlib.Path) -> dict[str, np.ndarray]:
    if not cache.exists():
        z = zarr.open(str(exp.store_path()), "r")
        ids = [str(c) for c in z["channel_ids"][:]]
        groups = np.asarray(z["properties"]["group"][:]).astype(int)
        bad = {str(b) for b in exp.badChannels}
        cols = [c for c in range(len(ids)) if groups[c] in tets and ids[c] not in bad]
        gain = float(z["properties"]["gain_to_uV"][0])
        f0, f1 = int(h0 * 3600 * FS) // CHUNK * CHUNK, int(np.ceil(h1 * 3600 * FS / CHUNK)) * CHUNK
        res = sorted(Parallel(n_jobs=16)(delayed(spont_chunk)(str(exp.store_path()), a, cols, groups[cols], tets, gain)
                                         for a in range(f0, f1, CHUNK)), key=lambda r: r[0])
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez(cache, start_s=f0 / FS, counts=np.concatenate([r[1] for r in res]), lfp=np.concatenate([r[2] for r in res]))
    z = np.load(cache)
    cnt, lfp, s0 = z["counts"].astype(float), z["lfp"], float(z["start_s"])
    n = cnt.size
    idx = np.arange(0, min(lfp.size - 4 * 500, (n - 4) * 500), 500)
    f, pxx = sps.welch(np.lib.stride_tricks.sliding_window_view(lfp, 4 * 500)[idx], fs=500, nperseg=1000, axis=-1)
    tot = pxx[:, (f >= 0.5) & (f <= 30)].sum(1)
    out = {"t": s0 + np.arange(n), "multi-unit rate (sqrt)": np.sqrt(cnt)}
    for name, (a, b) in {"delta share": (0.5, 4), "sigma share": (10, 15)}.items():
        v = np.full(n, np.nan)
        v[2:2 + idx.size] = (pxx[:, (f >= a) & (f <= b)].sum(1) / tot)[:n - 2]
        out[name] = v
    return out


def spectrum_shares(series: dict, pulses: pd.DataFrame, h0: float, h1: float) -> dict:
    t = series["t"]
    ok = (t >= h0 * 3600) & (t < h1 * 3600)
    for _, g in pulses.groupby("train"):
        ok &= ~((t >= g.onset_s.min() - 2) & (t < g.offset_s.max() + 12))
    for k, v in series.items():
        if k != "t":
            ok &= np.isfinite(v)
    edges = np.flatnonzero(np.diff(np.r_[0, ok.astype(int), 0]))
    out = {}
    for name, x in series.items():
        if name == "t":
            continue
        P, nseg = 0.0, 0
        for a, b in zip(edges[::2], edges[1::2]):
            for s in range(a, b - SEG + 1, SEG // 2):
                f, pp = sps.periodogram(sps.detrend(x[s:s + SEG]), fs=1.0, window="hann")
                P, nseg = P + pp, nseg + 1
        band = (f >= 1 / SEG) & (f <= 0.2)
        tot = P[band].sum()

        def share(lo, hi, f=f, P=P, tot=tot):
            return round(float(P[(f >= lo) & (f < hi)].sum() / tot), 3)

        pk = (f >= 0.005) & (f < 0.1)
        out[name] = {"segments": nseg, "above_120s": share(1 / SEG, 1 / 120), "10_120s": share(1 / 120, 0.1),
                     "5_10s": share(0.1, 0.2 + 1e-9), "peak_hz_power_x_f": round(float(f[pk][np.argmax((P * f)[pk])]), 4)}
    return out


def hour_mean_var(t_s: np.ndarray, c: np.ndarray) -> float:
    w = np.bincount(np.round(t_s).astype(int), minlength=c.size + 1) / t_s.size
    nfft = 1 << int(np.ceil(np.log2(2 * w.size)))
    F = np.fft.rfft(w, nfft)
    R = np.fft.irfft(F * np.conj(F), nfft)[:c.size]
    return SW2 / t_s.size + c[0] * R[0] + 2 * float((c[1:] * R[1:]).sum())


def precision() -> dict:
    out = {}
    for m in MODELS:
        c = kernel(m, 3700)
        rows = []
        for R, k in DESIGNS:
            sched = [draw(R * k // (k - 1), 3_600_000, 5000, 2000, k, np.random.default_rng(s)) for s in range(5)]
            v = float(np.mean([hour_mean_var(d.offset_ms[d.type == "pulse"].to_numpy() / 1000, c) for d in sched]))
            mdd = 2.8 * np.sqrt(2 * v)
            S = R // (k - 1)
            widen = {str(p): round(float(np.sqrt(1 + 2 * p * (1 - p) / S / (2 * v)) - 1), 3) for p in (0.005, 0.09)}
            rows.append({"pulses_per_hour": R, "shams_per_hour": S, "detectable_difference": round(float(mdd), 3),
                         "sleep_deprivations_for_0.10": int(np.ceil((mdd / 0.10) ** 2)),
                         "sleep_deprivations_for_0.05": int(np.ceil((mdd / 0.05) ** 2)), "sham_widening_by_sham_rate": widen})
        out[m] = {"designs": rows, "model_cov_in_lag_bins": {n: round(float(c[a:b + 1].mean()), 4) for n, (a, b) in LAG_BINS.items()}}
    return out


def intervals() -> list[dict]:
    rows = []
    for R, k in DESIGNS:
        ev, pp = [], []
        for s in range(20):
            d = draw(R * k // (k - 1), 3_600_000, 5000, 2000, k, np.random.default_rng(s))
            ev.append(np.diff(d.offset_ms) / 1000)
            pp.append(np.diff(d.offset_ms[d.type == "pulse"]) / 1000)
        ev, pp = np.concatenate(ev), np.concatenate(pp)
        rows.append({"pulses_per_hour": R, "events_per_hour": R * k // (k - 1), "event_interval_s_median": round(float(np.median(ev)), 2),
                     "event_interval_s_p1_p99": [round(float(q), 2) for q in np.percentile(ev, [1, 99])],
                     "pulse_interval_s_median": round(float(np.median(pp)), 2)})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    a = ap.parse_args()
    exp = experiment_from_args(a)
    js = json.loads((RESULTS / f"opto_off_periods_{exp.experiment}.json").read_text())
    tr = pd.read_csv(RESULTS / f"opto_off_periods_{exp.experiment}.csv")
    tr["grp"] = phase_groups(tr, js["state"]["light_transitions_h"])
    ev = pd.read_csv(exp.experiment_dir / "opto_off_periods" / "events.csv")
    pp = ev[ev.kind == "real"].merge(tr[["ttl_train", "onset_h", "grp", "in_sleep_deprivation"]], on="ttl_train")
    pp["off"] = pp.off.astype(float)
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    tets = [int(t[2:]) - 1 for t in js["population"]]
    sd0, sd1 = js["sleep_deprivation_h"]
    spont = {}
    for name, (h0, h1) in {"sleep_deprivation": (sd0, sd1), "session_hours_2_7": (2.0, 7.0)}.items():
        s = spontaneous(exp, h0, h1, tets, exp.experiment_dir / "sd_protocol_power" / f"spont_{name}.npz")
        spont[name] = spectrum_shares(s, pulses, h0, h1)
    out = {"experiment": exp.experiment, "population": js["population"], "within_train": within_train(pp),
           "between_trains": between_trains(pp), "spontaneous_spectra": spont, "models": precision(), "intervals": intervals()}
    (RESULTS / f"sd_protocol_power_{exp.experiment}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
