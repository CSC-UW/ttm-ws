"""Purity numbers and pictures for every unit of a long sorting, and the disagreements with another.

Per 1800 s window and tetrode the preprocessed binary is read once, the
spike-excised noise estimated, and every spike of every unit of the
tetrode gets its trough depth in MAD units and its cosine to the unit's
window-local template. From the whole trains ``ttm.capture.unit_purity``
then gives per unit the refractory contamination, the amplitude-split
zero-lag ratio and 2 to 3 ms burst ratio, and the cross-correlogram
verdict against every same-tetrode unit (``results/purity_<E>_<name>.csv``).
``figures/compare/purity_<name>.pdf`` holds one page per unit for a random
sample of ``--n-pages`` units and, with ``--other``, every unit on which the
two sortings disagree: units are paired by template cosine on the same
tetrode and the agreement is the fraction of the smaller train's spikes
within 0.5 ms of a spike of the other (``results/purity_match_<E>_<name>_vs_<other-name>.csv``);
pairs under ``--agree`` and unpaired units get a page in both files.

Usage: python subjects/ttm-003/scripts/sorting_purity.py --subject TTM-003 --experiment TTM-003-48H
           --sorting <dir>/assembled_candidate --name candidate [--dead-ms 0.3]
           [--other <dir>/assembled_prod --other-name baseline] [--n-pages 20] [--workers 48]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from tagged_cells_blocks import unit_templates  # noqa: E402

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
WIN_S = 1800.0


def _window_scores(binary: str, group: int, a: int, b: int, trains: dict) -> dict:
    """Per unit in this window: (local frames, amp_mad, cosine to the window template)."""
    rec = si.load(binary)
    x = cp.load_tetrode_traces(rec, group, a, b)
    noise, _ = cp.noise_mad_excised(x)
    local = {u: fr[(fr >= a + cp.NBEFORE) & (fr + cp.NAFTER <= b)] - a for u, fr in trains.items()}
    local = {u: fr for u, fr in local.items() if fr.size}
    out = {"group": group, "start": a, "noise": noise, "units": {}}
    if not local:
        return out
    T, ids = cp.templates_from_trains(x, {u: fr for u, fr in local.items() if fr.size >= 20} or local)
    for u, fr in local.items():
        sn = cp.snippets(x, fr)
        amp = -(sn / noise[None, None, :]).min(axis=(1, 2))
        if u in ids:
            r = cp.cosine_to_templates(x, fr, T[ids.index(u)][None], s_max=1)[:, 0]
        else:
            r = np.full(fr.size, np.nan)
        out["units"][u] = (fr + a, amp, r)
    return out


def score_sorting(binary: str, sorting, workers: int) -> tuple[dict, dict, dict, dict]:
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    units = [int(u) for u in sorting.unit_ids]
    groups = np.asarray(sorting.get_property("group")).astype(int)
    ugroup = dict(zip(units, groups.tolist()))
    trains = {u: np.sort(sorting.get_unit_spike_train(u).astype(np.int64)) for u in units}
    wf = int(WIN_S * FS)
    tasks = [(g, s, min(s + wf, nfr)) for s in range(0, nfr, wf) for g in sorted(set(groups.tolist()))]
    res = Parallel(n_jobs=workers, batch_size=1)(
        delayed(_window_scores)(binary, g, s, e, {u: trains[u] for u in units if ugroup[u] == g}) for g, s, e in tasks)
    amps = {u: [] for u in units}
    rs = {u: [] for u in units}
    frames = {u: [] for u in units}
    for r in res:
        for u, (fr, am, rr) in r["units"].items():
            frames[u].append(fr); amps[u].append(am); rs[u].append(rr)
    cat = lambda d, u, dtype: np.concatenate(d[u]) if d[u] else np.zeros(0, dtype)
    frames = {u: cat(frames, u, np.int64) for u in units}
    amps = {u: cat(amps, u, float) for u in units}
    rs = {u: cat(rs, u, float) for u in units}
    return frames, amps, rs, ugroup


def frac_within(ref: np.ndarray, other: np.ndarray, tol: int = 15) -> float:
    if ref.size == 0 or other.size == 0:
        return 0.0
    j = np.searchsorted(other, ref)
    lo = other[np.clip(j - 1, 0, other.size - 1)]
    hi = other[np.clip(j, 0, other.size - 1)]
    return float((np.minimum(np.abs(ref - lo), np.abs(ref - hi)) <= tol).mean())


def purity_table(frames: dict, amps: dict, ugroup: dict, nfr: int, dead_ms: float) -> pd.DataFrame:
    rows = []
    for g in sorted(set(ugroup.values())):
        members = [u for u in frames if ugroup[u] == g]
        pur = cp.unit_purity({u: frames[u] for u in members}, {u: amps[u] for u in members}, nfr, dead_ms=dead_ms)
        for u in members:
            d = pur[u]
            verdicts = list(d["neighbours"].values())
            rows.append({"unit": u, "tetrode": f"TT{g + 1}", "n_spikes": d["n_spikes"], "rate_hz": round(d["n_spikes"] / (nfr / FS), 3),
                         "amp_mad_median": round(float(np.median(amps[u])), 2) if amps[u].size else None,
                         "rp_contamination": None if np.isnan(d["rp_contamination"]) else round(d["rp_contamination"], 4),
                         "rp_violations": d["rp_violations"],
                         "ampsplit_ratio": None if d["ampsplit_ratio"] is None or np.isnan(d["ampsplit_ratio"]) else round(d["ampsplit_ratio"], 3),
                         "ampsplit_burst_ratio": None if d.get("ampsplit_burst_ratio") is None or np.isnan(d["ampsplit_burst_ratio"]) else round(d["ampsplit_burst_ratio"], 3),
                         "n_duplicate": verdicts.count("duplicate"), "n_ambiguous": verdicts.count("ambiguous"),
                         "n_distinct": verdicts.count("distinct"), "n_segregated": verdicts.count("segregated"),
                         "neighbours": " ".join(f"{v}:{d['neighbours'][v][:3]}" for v in sorted(d["neighbours"]))})
    return pd.DataFrame(rows)


def unit_page(axes, binary: str, u: int, frames: dict, amps: dict, rs: dict, ugroup: dict, templates: dict, gain: float,
              nfr: int, dead_ms: float, title: str) -> None:
    """One purity page from whole-recording trains; the overlays come from the window with most spikes."""
    fr, am, rr = frames[u], amps[u], rs[u]
    g = ugroup[u]
    wf = int(WIN_S * FS)
    tt = (np.arange(cp.N_SAMP) - cp.NBEFORE) / FS * 1000
    ax = axes[0, 0]
    if fr.size:
        w = int(np.bincount(fr // wf).argmax())
        a_, b_ = w * wf, min((w + 1) * wf, nfr)
        rec = si.load(binary)
        x = cp.load_tetrode_traces(rec, g, a_, b_)
        loc = fr[(fr >= a_ + cp.NBEFORE) & (fr + cp.NAFTER <= b_)] - a_
        rng = np.random.default_rng(0)
        pick = np.sort(rng.choice(loc, min(120, loc.size), replace=False)) if loc.size else loc
        if pick.size:
            sn = cp.snippets(x, pick) * gain
            for k in range(pick.size):
                for c in range(4):
                    ax.plot(tt + c * 3.2, sn[k, :, c], color="0.6", lw=0.3, alpha=0.5)
        ax.set_title(f"{title} unit {u} (TT{g + 1}): {fr.size} spikes; {pick.size} overlaid from the {a_ / FS / 3600:.1f} h window", fontsize=8)
    for c in range(4):
        ax.plot(tt + c * 3.2, templates[u][:, c] * gain, color="k", lw=1.2)
    ax.set_xlabel("ms, wires offset by 3.2 ms"); ax.set_ylabel("uV")
    ax = axes[0, 1]
    ax.scatter(rr, am, s=2, c="0.4", alpha=0.4)
    ax.set_xlabel("cosine to the unit's window template"); ax.set_ylabel("trough amplitude (MAD)"); ax.set_title("amplitude by shape", fontsize=8)
    ax = axes[0, 2]
    ax.scatter(fr / FS / 3600, am, s=1, c="0.4", alpha=0.4)
    ax.set_xlabel("session hour"); ax.set_ylabel("trough amplitude (MAD)"); ax.set_title("amplitude over time", fontsize=8)
    ax = axes[1, 0]
    isi = np.diff(fr) / FS * 1000
    c, nv = cp.rp_contamination(fr, nfr)
    ax.hist(isi[isi < 50], bins=np.arange(0, 50, 0.25), color="k")
    ax.axvline(1.0, color="r", lw=0.6)
    ax.set_title(f"ISI; rp contamination {c:.3f} ({nv} < 1 ms)", fontsize=8); ax.set_xlabel("ms")
    ax = axes[1, 1]
    sp = cp.amplitude_split(fr, am, nfr, dead_ms=dead_ms)
    if "lag_ms" in sp:
        ax.bar(sp["lag_ms"], sp["counts"], width=0.5, color="k")
        ax.axhline(sp["expected_per_bin"], color="gray", lw=0.8)
        ax.set_title(f"large (>= {sp['large_min']:.1f} MAD, n={sp['n_large']}) vs small (<= {sp['small_max']:.1f} MAD, n={sp['n_small']}): "
                     f"zero-lag ratio {sp['ratio']:.2f}, 2-3 ms {sp['burst_2_3ms_ratio']:.2f}", fontsize=7)
        ax.set_xlabel("ms (small relative to large)")
    ax = axes[1, 2]
    others = [v for v in sorted(frames) if v != u and ugroup[v] == g]
    max_lag = int(25 * FS / 1000)
    for k, v in enumerate(others):
        lags = cp.ccg_lags(fr, frames[v], max_lag) / FS * 1000
        h, e = np.histogram(lags, np.arange(-25, 25.5, 0.5))
        verdict = cp.ccg_verdict(fr, frames[v], dead_ms=dead_ms)
        base = h.mean() if h.size else 1
        lab = f"unit {v}: {verdict['verdict']}" + (f" ({verdict['ratio']:.2f})" if np.isfinite(verdict["ratio"]) else "")
        ax.plot(e[:-1] + 0.25, h / max(base, 1e-9) + 1.5 * k, lw=0.7, label=lab)
    ax.set_title("cross-correlograms with same-tetrode units (stacked, normalised)", fontsize=8)
    ax.set_xlabel("ms")
    if others:
        ax.legend(fontsize=6)


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--sorting", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--binary", default=None)
    p.add_argument("--dead-ms", type=float, default=0.0, help="detector dead time of this sorting (0.3 for the capture)")
    p.add_argument("--other", default=None, help="a second sorting to pair units with")
    p.add_argument("--other-name", default="other")
    p.add_argument("--other-dead-ms", type=float, default=0.0)
    p.add_argument("--agree", type=float, default=0.5)
    p.add_argument("--n-pages", type=int, default=20)
    p.add_argument("--workers", type=int, default=48)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    exp = experiment_from_args(a)
    binary = a.binary or str(exp.sorting_dir("mp_baseline") / "binary")
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    gain = float(rec.get_property("gain_to_uV")[0])
    figdir = exp.figures_dir / "compare"
    figdir.mkdir(parents=True, exist_ok=True)

    def one(path, name, dead):
        s = si.load(path)
        print(f"{name}: {s.get_num_units()} units; scoring every spike ...", flush=True)
        frames, amps, rs, ugroup = score_sorting(binary, s, a.workers)
        tab = purity_table(frames, amps, ugroup, nfr, dead)
        if "is_mua" in s.get_property_keys():
            mua = dict(zip([int(u) for u in s.unit_ids], np.asarray(s.get_property("is_mua")).astype(bool)))
            tab["is_mua"] = tab.unit.map(mua)
        tab.to_csv(RESULTS / f"purity_{exp.experiment}_{name}.csv", index=False)
        T, _ = unit_templates(binary, frames, ugroup, nfr, workers=min(a.workers, 16))
        return frames, amps, rs, ugroup, tab, T

    f1, a1, r1, g1, t1, T1 = one(a.sorting, a.name, a.dead_ms)
    pages = {a.name: set()}
    rng = np.random.default_rng(a.seed)
    sua1 = [u for u in f1 if not (("is_mua" in t1) and bool(t1.set_index("unit").loc[u, "is_mua"]))]
    pages[a.name] |= set(rng.choice(sua1, min(a.n_pages, len(sua1)), replace=False).tolist())
    if a.other:
        f2, a2, r2, g2, t2, T2 = one(a.other, a.other_name, a.other_dead_ms)
        pages[a.other_name] = set()
        sua2 = [u for u in f2 if not (("is_mua" in t2) and bool(t2.set_index("unit").loc[u, "is_mua"]))]
        rows = []
        used = set()
        for u in sua1:
            cands = [v for v in sua2 if g2[v] == g1[u]]
            if not cands:
                rows.append({a.name: u, a.other_name: None, "cosine": None, "agreement": 0.0}); continue
            cos = np.array([cp.template_cosine(T1[u], T2[v]) for v in cands])
            v = cands[int(np.argmax(cos))]
            used.add(v)
            small, big = (f1[u], f2[v]) if f1[u].size <= f2[v].size else (f2[v], f1[u])
            rows.append({a.name: u, a.other_name: v, "cosine": round(float(cos.max()), 4), "agreement": round(frac_within(small, big), 4),
                         f"n_{a.name}": int(f1[u].size), f"n_{a.other_name}": int(f2[v].size),
                         f"in_{a.other_name}": round(frac_within(f1[u], f2[v]), 4), f"in_{a.name}": round(frac_within(f2[v], f1[u]), 4)})
        for v in sua2:
            if v not in used:
                rows.append({a.name: None, a.other_name: v, "cosine": None, "agreement": 0.0, f"n_{a.other_name}": int(f2[v].size)})
        m = pd.DataFrame(rows)
        m.to_csv(RESULTS / f"purity_match_{exp.experiment}_{a.name}_vs_{a.other_name}.csv", index=False)
        dis = m[m.agreement < a.agree]
        pages[a.name] |= {int(u) for u in dis[a.name].dropna()}
        pages[a.other_name] |= {int(v) for v in dis[a.other_name].dropna()}
        sample2 = rng.choice(sua2, min(a.n_pages, len(sua2)), replace=False).tolist()
        pages[a.other_name] |= set(sample2)
        print(f"pairs: {len(m)}; agreement under {a.agree}: {len(dis)} (median agreement {m.agreement.median():.3f})", flush=True)
    for name, (fr, am, rr, ug, T, dead) in {a.name: (f1, a1, r1, g1, T1, a.dead_ms),
                                            **({a.other_name: (f2, a2, r2, g2, T2, a.other_dead_ms)} if a.other else {})}.items():
        with PdfPages(figdir / f"purity_{name}.pdf") as pdf:
            for u in sorted(pages[name], key=lambda u: (ug[u], u)):
                if fr[u].size == 0:
                    continue
                fig, axes = plt.subplots(2, 3, figsize=(20, 10))
                unit_page(axes, binary, u, fr, am, rr, ug, T, gain, nfr, dead, f"{exp.experiment} {name}")
                fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        print(f"wrote {figdir / f'purity_{name}.pdf'} ({len(pages[name])} pages)", flush=True)
    summ = {}
    for name, tab in {a.name: t1, **({a.other_name: t2} if a.other else {})}.items():
        s = tab[~tab.is_mua] if "is_mua" in tab else tab
        summ[name] = {"n_units": int(len(s)), "rp_contamination_median": float(s.rp_contamination.median()),
                      "units_rp_under_0.1": int((s.rp_contamination < 0.1).sum()),
                      "ampsplit_ratio_median": float(s.ampsplit_ratio.median()),
                      "units_ampsplit_under_0.3": int((s.ampsplit_ratio < 0.3).sum()),
                      "neighbour_verdicts": {k: int(s[f"n_{k}"].sum()) for k in ("duplicate", "ambiguous", "distinct", "segregated")}}
    (RESULTS / f"purity_{exp.experiment}_{a.name}_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
