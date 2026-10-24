"""Stage 1 head to head: what wobble admits that June's circus-omp does not, unit by unit, and the reverse.

Both start from the same bank with the same unit ids, so stage-1 units pair
by id. For each unit, a spike of June's train (``assembled_reseed_c12`` cut
to the span) is shared when the wobble train of the same unit has a spike
within 10 frames (0.33 ms); otherwise it is reassigned when a wobble spike of
another unit on the tetrode is that close, and dropped when none is. The
wobble spikes are split the same way into shared, taken from a neighbour and
new. Each spike's trough depth (excised-MAD units) and cosine to the unit's
own window template come from ``ttm.purity.score_sorting`` on each sorting.

Arm ``circus`` (the re-run of June's matcher from the same bank) goes
through the same comparison and gives its noise floor.

Amplitude classes come from the span-start bank (``template_amp_w55.npy``,
written by ``01_span_bank.py``; ``_span.amp_classes``): a template is large when its trough
on its deepest wire is at least 10 excised MAD, small otherwise; ``q`` is its
energy over the bank's median energy, and ``bar35`` the amplitude scale
wobble's objective threshold at 0.35 asks of it, ``(1 + 0.35 / q) / 2``.

Writes ``results/stage1_units_<arm>.csv``, ``results/stage1_summary.json``,
``results/bank_drift.csv`` and figures under the TTM-NOD figures directory
(``wobble_span/stage1_*``).

    OMP_NUM_THREADS=1 uv run python investigations/wobble-fullsize-ttm001/04_stage1_compare.py
"""
import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spikeinterface as si

from _span import FIG, FS, JUNE, LARGE_MAD, OUT, SPAN, WLEN, amp_classes, log, recording, restrict
from ttm import capture as cp
from ttm import purity as pu

HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
TOL = 10
ARMS = ["w35", "w25", "w35a"]
COLOURS = {"w35": "tab:blue", "w25": "tab:green", "w35a": "tab:purple", "circus": "0.3"}
GALLERY_WINDOW = 61


def within(q, ref, tol=TOL):
    if q.size == 0 or ref.size == 0:
        return np.zeros(q.size, bool)
    j = np.searchsorted(ref, q)
    lo = ref[np.clip(j - 1, 0, ref.size - 1)]
    hi = ref[np.clip(j, 0, ref.size - 1)]
    return np.minimum(np.abs(q - lo), np.abs(q - hi)) <= tol


def lookup(frames_scored, values, q):
    """Values of the scored spikes at frames q (exact match); NaN for spikes the scoring skipped."""
    j = np.searchsorted(frames_scored, q)
    j = np.clip(j, 0, max(frames_scored.size - 1, 0))
    ok = frames_scored.size > 0
    hit = (frames_scored[j] == q) if ok else np.zeros(q.size, bool)
    out = np.full(q.size, np.nan)
    out[hit] = values[j[hit]]
    return out


def med(x):
    x = x[np.isfinite(x)]
    return round(float(np.median(x)), 3) if x.size else None


def trains_of(s):
    return {int(u): np.sort(s.get_unit_spike_train(u).astype(np.int64)) for u in s.unit_ids}


def compare(cls, june, arm_s, sc_j, sc_w):
    tj, tw = trains_of(june), trains_of(arm_s)
    gj = dict(zip([int(u) for u in june.unit_ids], np.asarray(june.get_property("group")).astype(int)))
    rows, groups = [], {}
    for u in sorted(set(tj) | set(tw)):
        C = tj.get(u, np.zeros(0, np.int64))
        W = tw.get(u, np.zeros(0, np.int64))
        g = gj.get(u, cls.get(u, {}).get("tetrode"))
        oth_w = np.sort(np.concatenate([tw[v] for v in tw if v != u and gj.get(v) == g] or [np.zeros(0, np.int64)]))
        oth_c = np.sort(np.concatenate([tj[v] for v in tj if v != u and gj.get(v) == g] or [np.zeros(0, np.int64)]))
        c_sh = within(C, W)
        c_re = ~c_sh & within(C, oth_w)
        c_dr = ~c_sh & ~c_re
        w_sh = within(W, C)
        w_re = ~w_sh & within(W, oth_c)
        w_new = ~w_sh & ~w_re
        aj = lookup(sc_j[0].get(u, np.zeros(0, np.int64)), sc_j[1].get(u, np.zeros(0)), C)
        rj = lookup(sc_j[0].get(u, np.zeros(0, np.int64)), sc_j[2].get(u, np.zeros(0)), C)
        aw = lookup(sc_w[0].get(u, np.zeros(0, np.int64)), sc_w[1].get(u, np.zeros(0)), W)
        rw = lookup(sc_w[0].get(u, np.zeros(0, np.int64)), sc_w[2].get(u, np.zeros(0)), W)
        groups[u] = {"C": C, "W": W, "c_sh": c_sh, "c_re": c_re, "c_dr": c_dr, "w_sh": w_sh, "w_re": w_re,
                     "w_new": w_new, "aj": aj, "aw": aw}
        row = {"unit": u, **cls.get(u, {"tetrode": g}), "n_june": int(C.size), "n_wobble": int(W.size),
               "ratio": round(W.size / C.size, 3) if C.size else None,
               "june_shared": int(c_sh.sum()), "june_reassigned": int(c_re.sum()), "june_dropped": int(c_dr.sum()),
               "wobble_shared": int(w_sh.sum()), "wobble_from_neighbour": int(w_re.sum()), "wobble_new": int(w_new.sum()),
               "amp_shared": med(aw[w_sh]), "amp_wobble_new": med(aw[w_new]), "amp_wobble_from_neighbour": med(aw[w_re]),
               "amp_june_dropped": med(aj[c_dr]), "amp_june_reassigned": med(aj[c_re]),
               "cos_shared_wobble": med(rw[w_sh]), "cos_wobble_new": med(rw[w_new]),
               "cos_shared_june": med(rj[c_sh]), "cos_june_dropped": med(rj[c_dr]),
               "new_ge10mad": int((aw[w_new] >= 10).sum()), "dropped_ge10mad": int((aj[c_dr] >= 10).sum())}
        rows.append(row)
    return pd.DataFrame(rows), groups


def summarize(df):
    out = {}
    for name, d in [("all", df), ("large", df[df.amp_class == "large"]), ("small", df[df.amp_class == "small"]),
                    ("bar35_over_0.8", df[df.bar35 > 0.8])]:
        nj, nw = int(d.n_june.sum()), int(d.n_wobble.sum())
        out[name] = {"units": int(len(d)), "n_june": nj, "n_wobble": nw, "ratio": round(nw / nj, 3) if nj else None,
                     "june_dropped_frac": round(d.june_dropped.sum() / nj, 4) if nj else None,
                     "june_reassigned_frac": round(d.june_reassigned.sum() / nj, 4) if nj else None,
                     "wobble_new_frac": round(d.wobble_new.sum() / nw, 4) if nw else None,
                     "wobble_from_neighbour_frac": round(d.wobble_from_neighbour.sum() / nw, 4) if nw else None,
                     "median_unit_ratio": med(d.ratio.to_numpy(dtype=float)),
                     "units_ratio_over_1.5": int((d.ratio > 1.5).sum()), "units_ratio_under_0.67": int((d.ratio < 0.67).sum())}
    return out


def bank_drift():
    rows = []
    for arm in ["circus", *[a for a in ARMS if a != "circus"]]:
        for wi in SPAN:
            f = OUT / f"span_{arm}" / "raw" / f"w{wi}.npz"
            if not f.exists():
                continue
            z = np.load(f)
            T = z["templates4"].astype(float)
            e = (T ** 2).sum(axis=(1, 2))
            for k, u in enumerate(z["ids"]):
                rows.append({"arm": arm, "window": wi, "unit": int(u), "trough_raw": float(-T[k].min()),
                             "energy": float(e[k]), "threshold": float(z["threshold"]),
                             "n_matched": int((z["cluster_index"] == k).sum()),
                             "n_kept": int(((z["cluster_index"] == k) & (z["r"] >= 0.6)).sum()) if arm != "circus"
                             else int((z["cluster_index"] == k).sum())})
    return pd.DataFrame(rows)


def gallery(rec, gains, u, grp, g, tmpl, arm, path):
    """Rows: wobble shared, wobble new, wobble taken from a neighbour, June dropped, June reassigned;
    up to 40 snippets each from window GALLERY_WINDOW (or the unit's busiest), template in black."""
    C, W = grp["C"], grp["W"]
    allf = np.concatenate([C, W])
    wi = GALLERY_WINDOW if ((allf // WLEN) == GALLERY_WINDOW).any() else int(np.bincount(allf // WLEN).argmax())
    a, b = wi * WLEN, (wi + 1) * WLEN
    x = cp.load_tetrode_traces(rec, g, a, b)
    sets = [("wobble and June", W[grp["w_sh"]], grp["aw"][grp["w_sh"]], "0.5"),
            ("wobble only, new", W[grp["w_new"]], grp["aw"][grp["w_new"]], "tab:blue"),
            ("wobble only, a neighbour's in June", W[grp["w_re"]], grp["aw"][grp["w_re"]], "tab:purple"),
            ("June only, dropped by wobble", C[grp["c_dr"]], grp["aj"][grp["c_dr"]], "tab:red"),
            ("June only, a neighbour's in wobble", C[grp["c_re"]], grp["aj"][grp["c_re"]], "tab:orange")]
    tt = (np.arange(cp.N_SAMP) - cp.NBEFORE) / FS * 1000
    fig, axes = plt.subplots(len(sets), 2, figsize=(13, 2.6 * len(sets)), gridspec_kw={"width_ratios": [3, 1]})
    rng = np.random.default_rng(0)
    for k, (lab, fr, am, col) in enumerate(sets):
        ax, axh = axes[k, 0], axes[k, 1]
        loc = fr[(fr >= a + cp.NBEFORE + 2) & (fr + cp.NAFTER + 2 < b)] - a
        pick = np.sort(rng.choice(loc, min(40, loc.size), replace=False)) if loc.size else loc
        if pick.size:
            sn = cp.snippets(x, pick) * gains
            for s in sn:
                for c in range(4):
                    ax.plot(tt + c * 3.2, s[:, c], color=col, lw=0.4, alpha=0.45)
            for c in range(4):
                ax.plot(tt + c * 3.2, sn[:, :, c].mean(axis=0), color=col, lw=1.6, ls="--")
        for c in range(4):
            ax.plot(tt + c * 3.2, tmpl[:, c] * gains, color="k", lw=1.2)
        ax.set_title(f"{lab}: {fr.size:,} in the span, {pick.size} shown from window {wi}"
                     f" (median {np.nanmedian(am) if np.isfinite(am).any() else float('nan'):.1f} MAD)", fontsize=8)
        ax.set_ylabel("uV")
        ok = np.isfinite(am)
        if ok.any():
            axh.hist(am[ok], bins=np.linspace(0, max(30, np.nanpercentile(am, 99)), 60), color=col)
        axh.axvline(10, color="k", lw=0.6)
        axh.set_xlabel("trough (MAD)", fontsize=7)
    axes[-1, 0].set_xlabel("ms, wires offset by 3.2 ms")
    fig.suptitle(f"TTM-001 span 55-66, unit {u} (TT{g + 1}), {arm} against June circus-omp; black: June template "
                 f"in window {wi}, dashed: mean of the snippets shown", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def main():
    global ARMS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", nargs="+", default=ARMS)
    ARMS = ap.parse_args().arms
    RESULTS.mkdir(exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rec = recording()
    binary = str(JUNE / "binary")
    a0, b0 = SPAN[0] * WLEN, (SPAN[-1] + 1) * WLEN
    cls = amp_classes()
    june = restrict(si.load(JUNE / "assembled_reseed_c12"), a0, b0)
    june = june.select_units([u for u in june.unit_ids if june.get_unit_spike_train(u).size])
    log(f"June stage 1 over the span: {june.get_num_units()} units")
    sc_j = pu.score_sorting(binary, june, 16, cache=OUT / "scores_stage1_june.npz", frame_range=(a0, b0))
    drift = bank_drift()
    drift.to_csv(RESULTS / "bank_drift.csv", index=False)
    summary = {}
    gains = float(rec.get_property("gain_to_uV")[0])
    for arm in ARMS:
        s = si.load(OUT / f"span_{arm}" / "sorting")
        sc_w = pu.score_sorting(binary, s, 16, cache=OUT / f"scores_stage1_{arm}.npz", frame_range=(a0, b0))
        df, groups = compare(cls, june, s, sc_j, sc_w)
        df.to_csv(RESULTS / f"stage1_units_{arm}.csv", index=False)
        summary[arm] = summarize(df)
        log(f"{arm}: {json.dumps(summary[arm]['all'])}")
        # galleries: per class, the units where the two matchers differ most
        df["diff"] = (df.wobble_new + df.june_dropped + df.wobble_from_neighbour + df.june_reassigned) / df.n_june.clip(lower=1)
        T, _ = pu.unit_templates(binary, trains_of(june), dict(zip([int(u) for u in june.unit_ids],
                                 np.asarray(june.get_property("group")).astype(int))), rec.get_num_frames(),
                                 workers=16, frame_range=(GALLERY_WINDOW * WLEN, (GALLERY_WINDOW + 1) * WLEN))
        for c in ("large", "small"):
            d = df[(df.amp_class == c) & (df.n_june >= 500)].sort_values("diff", ascending=False)
            for u in d.unit.head(4):
                gallery(rec, gains, int(u), groups[int(u)], int(cls[int(u)]["tetrode"]), T[int(u)], arm,
                        FIG / f"stage1_gallery_{arm}_{c}_u{int(u)}.png")
    old = json.loads((RESULTS / "stage1_summary.json").read_text()) if (RESULTS / "stage1_summary.json").exists() else {}
    old.update(summary)
    (RESULTS / "stage1_summary.json").write_text(json.dumps(old, indent=1))

    # overview: count ratio against template amplitude; amplitude of what each admits alone
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for arm in [a for a in ARMS if a != "circus"]:
        col = COLOURS[arm]
        df = pd.read_csv(RESULTS / f"stage1_units_{arm}.csv")
        axes[0].scatter(df.template_mad, df.ratio, s=14, color=col, label=arm)
    axes[0].axhline(1, color="k", lw=0.6); axes[0].axvline(LARGE_MAD, color="k", lw=0.6, ls=":")
    axes[0].set_xscale("log"); axes[0].set_yscale("log")
    axes[0].set_xlabel("template trough at window 55 (excised MAD)"); axes[0].set_ylabel("spikes wobble / spikes June")
    axes[0].legend(fontsize=8); axes[0].set_title("stage-1 spike count per unit, span 55-66", fontsize=9)
    d = drift.pivot_table(index="window", columns="arm", values="trough_raw", aggfunc="median") * gains
    for arm in d:
        axes[1].plot(d.index, d[arm], marker="o", label=arm, color=COLOURS.get(arm, "k"))
    axes[1].axvline(60.5, color="k", lw=0.6, ls=":")
    axes[1].set_xlabel("window (re-seed after 60)"); axes[1].set_ylabel("median template trough (uV)")
    axes[1].set_title("bank templates entering each window", fontsize=9); axes[1].legend(fontsize=8)
    piv = drift[drift.arm != "circus"].merge(drift[drift.arm == "circus"][["window", "unit", "trough_raw"]],
                                             on=["window", "unit"], suffixes=("", "_circus"))
    last = piv[piv.window == SPAN[-1]].copy()
    last["shrink"] = last.trough_raw / last.trough_raw_circus
    last["template_mad"] = last.unit.map({u: v["template_mad"] for u, v in cls.items()})
    for arm in [a for a in ARMS if a != "circus"]:
        col = COLOURS[arm]
        m = last[last.arm == arm]
        axes[2].scatter(m.template_mad, m.shrink, s=14, color=col, label=arm)
    axes[2].axhline(1, color="k", lw=0.6); axes[2].set_xscale("log")
    axes[2].set_xlabel("template trough at window 55 (excised MAD)")
    axes[2].set_ylabel("template trough, wobble bank / circus bank")
    axes[2].set_title(f"bank entering window {SPAN[-1]}", fontsize=9); axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "stage1_overview.png", dpi=100)
    plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
