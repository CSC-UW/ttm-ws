"""Optotagging metrics and SOM-candidate characterization for one sorted session.

Per unit and per pulse pattern: PSTH, first-spike latency and jitter,
reliability, a SALT-style first-spike-latency test, a signed-rank test of
per-pulse counts, rate modulation during and after the pulse, and the
correlation between in-pulse and out-of-train waveforms. Units are then
classed as activated, suppressed or unaffected per pattern. For the
ranked tagging candidates: cross-correlograms with every other unit
outside the trains, spike-triggered LFP averages and band phase locking.

Writes ``results/optotag_<experiment>_<variant>.csv`` (one row per unit
and pattern) and ``results/optotag_<experiment>_<variant>_candidates.json``;
figures go to ``<experiment>/figures/optotag_*.png``.

Usage: python subjects/ttm-003/scripts/opto_tag.py --subject TTM-003 --experiment TTM-003-OPTO1 --variant asis
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as sps
import scipy.stats as sst
import spikeinterface as si

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
SORTINGS_SUBDIR = "sortings_opto"
FS = 30000
EARLY = (0.001, 0.010)        # s after onset: direct-activation window (skip the first ms)
BASE = (-0.050, 0.0)
BANDS = {"delta": (0.5, 4), "theta": (5, 10), "sigma": (10, 16), "gamma": (30, 80)}


# --- spike-train helpers -----------------------------------------------------

def rel_times(spikes: np.ndarray, events: np.ndarray, lo: float, hi: float) -> list[np.ndarray]:
    i0 = np.searchsorted(spikes, events + lo)
    i1 = np.searchsorted(spikes, events + hi)
    return [spikes[a:b] - e for a, b, e in zip(i0, i1, events)]


def counts(spikes, events, lo, hi) -> np.ndarray:
    return np.array([len(r) for r in rel_times(spikes, events, lo, hi)])


def psth(spikes, events, lo, hi, bin_s):
    edges = np.arange(lo, hi + bin_s / 2, bin_s)
    rel = np.concatenate(rel_times(spikes, events, lo, hi)) if len(events) else np.array([])
    return edges[:-1] + bin_s / 2, np.histogram(rel, edges)[0] / (len(events) * bin_s)


def first_latency(spikes, events, lo, hi) -> np.ndarray:
    out = np.full(len(events), np.nan)
    for k, r in enumerate(rel_times(spikes, events, lo, hi)):
        if r.size:
            out[k] = r[0]
    return out


def js_divergence(p, q):
    m = 0.5 * (p + q)
    def kl(a, b):
        mask = a > 0
        return np.sum(a[mask] * np.log2(a[mask] / b[mask]))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def salt(spikes, events, win=0.010, n_base=20, nbins=10):
    """SALT (Kvitsiani et al. 2013): first-spike latency distribution in the
    test window vs equal-length baseline windows, compared by Jensen-Shannon
    divergence. Returns (p, test_median_JS, baseline_median_JS)."""
    edges = np.linspace(0, win, nbins + 1)
    def dist(lo):
        lat = first_latency(spikes, events, lo, lo + win) - lo
        h = np.histogram(lat[~np.isnan(lat)], edges)[0].astype(float)
        h = np.r_[h, np.isnan(lat).sum()]
        return h / max(h.sum(), 1)
    base = [dist(-(k + 1) * win - 0.010) for k in range(n_base)]
    test = dist(EARLY[0])
    bb = [js_divergence(base[i], base[j]) for i in range(n_base) for j in range(i + 1, n_base)]
    tb = [js_divergence(test, b) for b in base]
    bb, tb = np.sqrt(np.array(bb)), np.sqrt(np.array(tb))
    p = (np.sum(bb >= np.median(tb)) + 1) / (len(bb) + 1)
    return float(p), float(np.median(tb)), float(np.median(bb))


def ccg(a, b, window=0.05, bin_s=0.001):
    """Counts of b spikes at lag (b - a) in bins centred on integer multiples of
    bin_s; exact zero lags are dropped (self pairs of an autocorrelogram)."""
    nb = int(round(window / bin_s))
    edges = (np.arange(-nb, nb + 2) - 0.5) * bin_s
    lags = []
    i0 = np.searchsorted(b, a - window - bin_s)
    i1 = np.searchsorted(b, a + window + bin_s)
    for t, s, e in zip(a, i0, i1):
        if e > s:
            lags.append(b[s:e] - t)
    lags = np.concatenate(lags) if lags else np.array([])
    lags = lags[lags != 0]
    return edges[:-1] + bin_s / 2, np.histogram(lags, edges)[0]


def outside_trains(spikes, trains, pad=1.0):
    keep = np.ones(len(spikes), bool)
    for tr in trains.itertuples():
        keep &= ~((spikes > tr.onset_s - pad) & (spikes < tr.offset_s + pad))
    return spikes[keep]


# --- main --------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis")
    p.add_argument("--n-candidates", type=int, default=4)
    a = p.parse_args()
    exp = experiment_from_args(a)
    figdir = exp.figures_dir
    figdir.mkdir(parents=True, exist_ok=True)
    tag = f"{exp.experiment}_{a.variant}"

    out = exp.sorting_dir(a.variant, subdir=SORTINGS_SUBDIR)
    analyzer = si.load_sorting_analyzer(str(out / "analyzer"))
    sorting = analyzer.sorting
    rec = analyzer.recording
    t = rec.get_times()
    t0 = t[0]
    groups = pd.Series(sorting.get_property("group"), index=sorting.unit_ids)
    qm = analyzer.get_extension("quality_metrics").get_data()
    templates = analyzer.get_extension("templates")
    sparsity = analyzer.sparsity

    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv")
    trains = trains[trains.session == exp.session]
    train_power = dict(zip(trains.ttl_train, trains.power_mW))
    train_pattern = dict(zip(trains.ttl_train, trains.ttl_pattern))
    pulses = pulses.assign(power_mW=pulses.train.map(train_power), pattern=pulses.train.map(train_pattern))

    # session-relative spike times (seconds from the first sample), like the pulse CSVs
    spikes = {u: t[sorting.get_unit_spike_train(u)] - t0 for u in sorting.unit_ids}

    # pattern sets: strongest power per pattern (A at 4.2 mW, B at 9.8 mW), and the sub-threshold A trains
    sets = {}
    for pat in ("A", "B"):
        pp = pulses[pulses.pattern == pat]
        if len(pp):
            pmax = pp.power_mW.max()
            sets[pat] = pp[pp.power_mW == pmax]
    low = pulses[(pulses.pattern == "A") & (pulses.power_mW < 1)]
    if len(low):
        sets["A_low"] = low

    rows = []
    for u in sorting.unit_ids:
        sp = spikes[u]
        base_row = {"unit": u, "group": int(groups[u]), "tetrode": f"TT{int(groups[u]) + 1}", "n_spikes": int(len(sp)),
                    "rate_hz": float(len(sp) / (t[-1] - t0)), "snr": float(qm.loc[u, "snr"]),
                    "isi_viol_ratio": float(qm.loc[u, "isi_violations_ratio"]),
                    "rp_contamination": float(qm.loc[u, "rp_contamination"])}
        for name, pp in sets.items():
            on, off = pp.onset_s.values, pp.offset_s.values
            dur = float(np.median(off - on))
            n_base = counts(sp, on, *BASE)
            n_early = counts(sp, on, *EARLY)
            n_pulse = counts(sp, on, EARLY[1], dur)
            n_post = counts(sp, off, 0.0, 0.050)
            n_post2 = counts(sp, off, 0.050, 0.200)
            base_hz = n_base.sum() / (len(on) * (BASE[1] - BASE[0]))
            early_hz = n_early.sum() / (len(on) * (EARLY[1] - EARLY[0]))
            pulse_hz = n_pulse.sum() / (len(on) * (dur - EARLY[1]))
            post_hz = n_post.sum() / (len(on) * 0.050)
            post2_hz = n_post2.sum() / (len(on) * 0.150)
            lat = first_latency(sp, on, EARLY[0], 0.020)
            ok = ~np.isnan(lat)
            # signed-rank: early window count vs an equal-length baseline window
            n_base_eq = counts(sp, on, BASE[0], BASE[0] + (EARLY[1] - EARLY[0]))
            try:
                p_early = float(sst.wilcoxon(n_early, n_base_eq, alternative="greater").pvalue) if np.any(n_early != n_base_eq) else 1.0
            except ValueError:
                p_early = 1.0
            n_base_eq2 = counts(sp, on, BASE[0], BASE[0] + (dur - EARLY[1])) if dur - EARLY[1] <= -BASE[0] else counts(sp, on - 1.0 + dur, 0, dur - EARLY[1])
            try:
                p_supp = float(sst.wilcoxon(n_pulse, n_base_eq2, alternative="less").pvalue) if np.any(n_pulse != n_base_eq2) else 1.0
            except ValueError:
                p_supp = 1.0
            p_salt, js_t, js_b = salt(sp, on)
            # in-pulse vs out-of-train waveform correlation
            in_pulse = np.concatenate([sp[(sp > s + EARLY[0]) & (sp < s + dur)] for s in on]) if len(on) else np.array([])
            wf_r = np.nan
            if len(in_pulse) >= 20:
                st_in = np.searchsorted(t - t0, in_pulse)
                out_sp = outside_trains(sp, trains)
                rng = np.random.default_rng(0)
                st_out = np.searchsorted(t - t0, rng.choice(out_sp, min(500, len(out_sp)), replace=False)) if len(out_sp) >= 20 else None
                if st_out is not None:
                    chans = sparsity.unit_id_to_channel_indices[u]
                    def mean_wf(idx):
                        idx = idx[(idx > 60) & (idx < rec.get_num_frames() - 60)]
                        acc = np.zeros((90, len(chans)))
                        for i in idx[:500]:
                            acc += rec.get_traces(start_frame=i - 30, end_frame=i + 60, channel_ids=rec.channel_ids[chans])
                        return acc / len(idx[:500])
                    w_in, w_out = mean_wf(st_in), mean_wf(st_out)
                    wf_r = float(np.corrcoef(w_in.ravel(), w_out.ravel())[0, 1])
            rows.append({**base_row, "set": name, "n_pulses": int(len(on)), "pulse_s": round(dur, 3),
                         "base_hz": base_hz, "early_hz": early_hz, "pulse_hz": pulse_hz,
                         "post0_50_hz": post_hz, "post50_200_hz": post2_hz,
                         "early_ratio": (early_hz + 0.1) / (base_hz + 0.1), "pulse_ratio": (pulse_hz + 0.1) / (base_hz + 0.1),
                         "post_ratio": (post_hz + 0.1) / (base_hz + 0.1),
                         "reliability": float(ok.mean()),
                         "latency_ms": float(np.nanmedian(lat) * 1000) if ok.any() else np.nan,
                         "jitter_ms": float(np.nanstd(lat) * 1000) if ok.sum() > 1 else np.nan,
                         "p_early": p_early, "p_supp": p_supp, "p_salt": p_salt, "salt_js_test": js_t, "salt_js_base": js_b,
                         "wf_r_in_vs_out": wf_r})
    df = pd.DataFrame(rows)

    def classify(r):
        if r.p_salt < 0.01 and r.early_ratio > 3 and r.reliability > 0.3 and (np.isnan(r.wf_r_in_vs_out) or r.wf_r_in_vs_out > 0.85):
            return "activated"
        if r.p_supp < 0.05 and r.pulse_ratio < 0.5:
            return "suppressed"
        if r.p_early < 0.01 and r.early_ratio > 2:
            return "activated_weak"
        return "unaffected"
    df["class"] = df.apply(classify, axis=1)
    df.to_csv(RESULTS / f"optotag_{tag}.csv", index=False)
    print(df.groupby(["set", "class"]).size().unstack(fill_value=0))

    # --- candidates: activated in B (short pulses) ranked by SALT JS and reliability ---
    cand = df[(df.set == "B") & (df["class"].str.startswith("activated"))].copy()
    cand["score"] = cand.reliability * (cand.salt_js_test - cand.salt_js_base).clip(lower=0) * (cand.early_ratio.clip(upper=50))
    cand = cand.sort_values("score", ascending=False)
    print("\ncandidates (pattern B):")
    print(cand[["unit", "tetrode", "rate_hz", "snr", "base_hz", "early_hz", "early_ratio", "reliability", "latency_ms", "jitter_ms", "p_salt", "wf_r_in_vs_out", "class"]].head(12).to_string(index=False))

    # LFP: raw store low-passed, one lead channel per tetrode plus the candidate's home tetrode
    raw = si.read_zarr(str(exp.store_path()))
    lfp_fs = 1000
    q = FS // lfp_fs
    cand_json = {"experiment": exp.experiment, "variant": a.variant, "candidates": []}
    for u in cand.unit.head(a.n_candidates):
        u = type(sorting.unit_ids[0])(u)
        g = int(groups[u])
        sp_all = spikes[u]
        sp_out = outside_trains(sp_all, trains)
        info = {"unit": int(u) if np.issubdtype(type(u), np.integer) else str(u), "group": g, "tetrode": f"TT{g + 1}",
                "n_spikes": int(len(sp_all)), "n_spikes_outside_trains": int(len(sp_out))}

        fig, axes = plt.subplots(3, 4, figsize=(20, 13))
        # template
        ax = axes[0, 0]
        chans = sparsity.unit_id_to_channel_indices[u]
        tmpl = templates.get_unit_template(u, operator="average")
        for c in chans:
            ax.plot(np.arange(tmpl.shape[0]) / FS * 1000 - 1.0, tmpl[:, c], label=str(rec.channel_ids[c]))
        ax.set_title(f"unit {u} TT{g + 1} template (uV)"); ax.legend(fontsize=7); ax.set_xlabel("ms")
        # ACG + ISI
        ax = axes[0, 1]
        lags, c = ccg(sp_out, sp_out, 0.05, 0.001)
        ax.bar(lags * 1000, c, width=1.0, color="k"); ax.set_title("ACG outside trains (1 ms bins)"); ax.set_xlabel("ms")
        ax = axes[0, 2]
        isi = np.diff(sp_out) * 1000
        ax.hist(isi[isi < 100], bins=np.arange(0, 100, 1), color="k"); ax.set_title("ISI (ms)")
        isi_viol = float((isi < 2).mean()) if isi.size else np.nan
        info["isi_lt2ms_fraction_outside"] = isi_viol
        # PSTHs per set
        for k, name in enumerate([n for n in ("B", "A", "A_low") if n in sets]):
            pp = sets[name]
            on = pp.onset_s.values
            dur = float(np.median(pp.offset_s.values - on))
            ax = axes[1, k]
            rel = rel_times(sp_all, on, -0.05, dur + 0.1)
            for i, r in enumerate(rel):
                ax.plot(r * 1000, np.full(len(r), i), "|", color="k", ms=2)
            ax.axvspan(0, dur * 1000, color="c", alpha=0.2); ax.set_title(f"{name}: raster ({len(on)} pulses)"); ax.set_xlabel("ms")
            ax = axes[2, k]
            tb, r = psth(sp_all, on, -0.05, dur + 0.1, 0.002)
            ax.plot(tb * 1000, r, color="k"); ax.axvspan(0, dur * 1000, color="c", alpha=0.2); ax.set_title(f"{name}: PSTH (2 ms bins, Hz)")
            ax2 = ax.inset_axes([0.55, 0.55, 0.43, 0.43])
            tb2, r2 = psth(sp_all, on, -0.005, 0.025, 0.0005)
            ax2.plot(tb2 * 1000, r2, color="k"); ax2.set_title("onset, 0.5 ms bins", fontsize=7)
        # latency histogram (B)
        ax = axes[0, 3]
        lat = first_latency(sp_all, sets["B"].onset_s.values, EARLY[0], 0.030) * 1000 if "B" in sets else np.array([])
        ax.hist(lat[~np.isnan(lat)], bins=np.arange(0, 30, 0.5), color="k"); ax.set_title("first-spike latency, pattern B (ms)")
        # spike-count time course within the 200 ms pulse (10 ms bins) across the three A trains
        ax = axes[1, 3]
        if "A" in sets:
            for trn, pp in sets["A"].groupby("train"):
                tb, r = psth(sp_all, pp.onset_s.values, -0.1, 0.4, 0.010)
                ax.plot(tb * 1000, r, label=f"train {trn} ({train_power[trn]:g} mW)")
            ax.axvspan(0, 207, color="c", alpha=0.2); ax.legend(fontsize=7); ax.set_title("A: PSTH per train (10 ms bins)")
        # pulse-by-pulse reliability through the B trains
        ax = axes[2, 3]
        if "B" in sets:
            for trn, pp in sets["B"].groupby("train"):
                n = counts(sp_all, pp.onset_s.values, *EARLY)
                ax.plot(np.convolve(n > 0, np.ones(20) / 20, mode="valid"), label=f"train {trn}")
            ax.set_title("B: reliability (20-pulse running mean) through the train"); ax.set_xlabel("pulse #"); ax.legend(fontsize=7)
        fig.suptitle(f"{tag} candidate unit {u} on TT{g + 1}")
        fig.tight_layout(); fig.savefig(figdir / f"optotag_{a.variant}_unit{u}_summary.png", dpi=100); plt.close(fig)

        # CCGs with all other units outside the trains
        others = [v for v in sorting.unit_ids if v != u]
        ncol = 8
        nrow = int(np.ceil(len(others) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.6 * ncol, 2.2 * nrow), squeeze=False)
        ccg_rows = []
        for k, v in enumerate(others):
            sv = outside_trains(spikes[v], trains)
            lags, c = ccg(sp_out, sv, 0.05, 0.001)
            base = c[(np.abs(lags) > 0.020)].mean()
            sd = max(np.sqrt(base), 1.0)
            # inhibitory signature: trough at +1..+5 ms (reference spikes of u precede)
            post = c[(lags > 0.0005) & (lags < 0.0055)]
            pre = c[(lags < -0.0005) & (lags > -0.0055)]
            ccg_rows.append({"other": int(v) if np.issubdtype(type(v), np.integer) else str(v), "other_group": int(groups[v]), "other_tetrode": f"TT{int(groups[v]) + 1}",
                             "n_ref": int(len(sp_out)), "n_other": int(len(sv)), "baseline_per_bin": float(base),
                             "post_1_5ms_z": float((post.mean() - base) / sd * np.sqrt(len(post))),
                             "pre_1_5ms_z": float((pre.mean() - base) / sd * np.sqrt(len(pre))),
                             "zero_bin_z": float((c[len(c) // 2] - base) / sd)})
            ax = axes[k // ncol, k % ncol]
            ax.bar(lags * 1000, c, width=1.0, color="k" if groups[v] != g else "tab:red")
            ax.axhline(base, color="gray", lw=0.5)
            ax.set_title(f"u{v} TT{groups[v] + 1} z+={ccg_rows[-1]['post_1_5ms_z']:.1f}", fontsize=7); ax.tick_params(labelsize=6)
        for k in range(len(others), nrow * ncol):
            axes[k // ncol, k % ncol].axis("off")
        fig.suptitle(f"{tag}: CCG of unit {u} (reference, TT{g + 1}) vs every other unit, outside trains; red = same tetrode")
        fig.tight_layout(); fig.savefig(figdir / f"optotag_{a.variant}_unit{u}_ccgs.png", dpi=90); plt.close(fig)
        info["ccg"] = ccg_rows

        # spike-triggered LFP + phase locking on the home tetrode (lead channel), outside trains
        home_ids = [cid for cid, gg in zip(raw.channel_ids, raw.get_property("group")) if int(gg) == g]
        lead = home_ids[0]
        x = raw.get_traces(channel_ids=[lead], return_in_uV=True)[:, 0].astype(np.float32)
        x = sps.decimate(x, q, ftype="fir", zero_phase=True)
        tl = (t - t0)[::q][: len(x)]
        sp_idx = np.searchsorted(tl, sp_out)
        sp_idx = sp_idx[(sp_idx > lfp_fs) & (sp_idx < len(x) - lfp_fs)]
        win = np.arange(-500, 501)
        sta = np.mean([x[i + win] for i in sp_idx[:20000]], axis=0) if len(sp_idx) else np.zeros(len(win))
        locking = {}
        for band, (lo, hi) in BANDS.items():
            sos = sps.butter(3, [lo, hi], btype="bandpass", fs=lfp_fs, output="sos")
            ph = np.angle(sps.hilbert(sps.sosfiltfilt(sos, x)))
            phs = ph[sp_idx]
            if len(phs) > 50:
                R = np.abs(np.mean(np.exp(1j * phs)))
                z = len(phs) * R ** 2
                locking[band] = {"n": int(len(phs)), "mrl": float(R), "mean_phase_deg": float(np.degrees(np.angle(np.mean(np.exp(1j * phs))))),
                                 "rayleigh_p": float(np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * len(phs))))}
        info["spike_lfp_home_lead"] = str(lead)
        info["phase_locking_outside_trains"] = locking
        fig, axes = plt.subplots(1, 2, figsize=(14, 4))
        axes[0].plot(win, sta, color="k"); axes[0].axvline(0, color="r", lw=0.5)
        axes[0].set_title(f"spike-triggered LFP, {lead} (TT{g + 1}), n={min(len(sp_idx), 20000)}, outside trains"); axes[0].set_xlabel("ms"); axes[0].set_ylabel("uV")
        for band, (lo, hi) in BANDS.items():
            if band in locking:
                sos = sps.butter(3, [lo, hi], btype="bandpass", fs=lfp_fs, output="sos")
                ph = np.angle(sps.hilbert(sps.sosfiltfilt(sos, x)))[sp_idx]
                axes[1].hist(np.degrees(ph), bins=36, histtype="step", label=f"{band} MRL={locking[band]['mrl']:.3f} p={locking[band]['rayleigh_p']:.2g}")
        axes[1].legend(fontsize=7); axes[1].set_title("spike phase distributions"); axes[1].set_xlabel("deg")
        fig.tight_layout(); fig.savefig(figdir / f"optotag_{a.variant}_unit{u}_lfp.png", dpi=100); plt.close(fig)
        cand_json["candidates"].append(info)

    (RESULTS / f"optotag_{tag}_candidates.json").write_text(json.dumps(cand_json, indent=1) + "\n")

    # population figure: PSTH of every unit, z-scored, sorted by tetrode, for B and A
    for name in [n for n in ("B", "A") if n in sets]:
        pp = sets[name]; on = pp.onset_s.values; dur = float(np.median(pp.offset_s.values - on))
        mats = []
        for u in sorting.unit_ids:
            tb, r = psth(spikes[u], on, -0.1, dur + 0.2, 0.005)
            b = r[tb < 0]
            mats.append((r - b.mean()) / (b.std() + 1e-3))
        order = np.argsort(groups.values, kind="stable")
        fig, ax = plt.subplots(figsize=(12, 0.25 * len(order) + 2))
        im = ax.imshow(np.array(mats)[order], aspect="auto", extent=[tb[0] * 1000, tb[-1] * 1000, len(order) - 0.5, -0.5], cmap="RdBu_r", vmin=-5, vmax=5)
        ax.axvline(0, color="k", lw=0.5); ax.axvline(dur * 1000, color="k", lw=0.5)
        ax.set_yticks(range(len(order))); ax.set_yticklabels([f"u{sorting.unit_ids[i]} TT{groups.values[i] + 1}" for i in order], fontsize=6)
        ax.set_title(f"{tag} pattern {name}: z-scored PSTH per unit (5 ms bins)"); ax.set_xlabel("ms from onset")
        plt.colorbar(im, ax=ax); fig.tight_layout(); fig.savefig(figdir / f"optotag_{a.variant}_population_{name}.png", dpi=100); plt.close(fig)
    print(f"wrote results and figures in {figdir}")


if __name__ == "__main__":
    main()
