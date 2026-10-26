"""The staged 48 h sorting against the June baseline and the candidate: admissions, purity, presence.

Reads ``mp_staged/`` (``admission.csv``, ``stage3_counts.npz``,
``staged_merges.json``, ``assembled_staged``) and the baseline's
``assembled_prod`` and the candidate's ``assembled_candidate``. Writes

- ``results/staged_admission_<E>.csv``: the light test of every unit at every
  delivered train (a copy of ``admission.csv``), and
  ``figures/staged_admission.png``: which units stage 3 ran for, train by train;
- ``results/purity_<E>_staged.csv`` (``ttm.purity.purity_table`` of the staged
  single and multi-unit units) and ``results/staged_units_<E>.csv``: per June
  unit its stage-3 spike count, the trains it was admitted at, and the purity
  numbers of its staged unit next to the baseline's (``results/purity_<E>_baseline.csv``)
  and the refractory conflict ratio of its stage-3 spikes;
- ``figures/compare/purity_staged.pdf`` and ``purity_staged_baseline.pdf``: pages
  for every unit stage 3 touched and for the TT7 and TT9 units the candidate's
  template floor dropped, in the staged sorting and in the baseline;
- ``results/staged_presence_<E>.csv``: per tetrode the median fraction of
  windows (from each unit's first) in which a unit has 20 spikes or more, for
  the three sortings.

With ``--tag T`` the run ``mp_staged<T>`` is read and the outputs carry the
tag (``purity_<E>_staged<T>.csv``, ``staged_units_<E><T>.csv``,
``figures/compare/purity_staged<T>.pdf``, ...); for a run made from saved
offers with a post hoc filter, the per-unit conflict ratios come from
``staged_split.py``'s table, split into in-pulse and between-pulse spikes.

Usage: python subjects/ttm-003/scripts/staged_compare.py --subject TTM-003 --experiment TTM-003-48H [--workers 48]
           [--tag _split_f5.5_g0.95 --no-baseline-pages]
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
import spikeinterface as si

from ttm import purity as pu
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
WIN_S = 1800.0
FLOOR_TETRODES = ("TT7", "TT9")


def presence(sorting, nfr: int) -> pd.DataFrame:
    wf = int(WIN_S * FS)
    nwin = int(np.ceil(nfr / wf))
    groups = np.asarray(sorting.get_property("group")).astype(int)
    mua = (np.asarray(sorting.get_property("is_mua")).astype(bool) if "is_mua" in sorting.get_property_keys()
           else np.zeros(len(groups), bool))
    rows = []
    for u, g, m in zip(sorting.unit_ids, groups, mua):
        if m:
            continue
        c = np.bincount(sorting.get_unit_spike_train(u) // wf, minlength=nwin)[:nwin]
        present = np.flatnonzero(c >= 20)
        if present.size == 0:
            rows.append({"unit": int(u), "tetrode": f"TT{g + 1}", "presence_from_first": 0.0, "n_spikes": int(c.sum())})
            continue
        rows.append({"unit": int(u), "tetrode": f"TT{g + 1}", "presence_from_first": float((c[present[0]:] >= 20).mean()),
                     "n_spikes": int(c.sum())})
    return pd.DataFrame(rows)


def conflict_ratios(out: pathlib.Path, base_dir: pathlib.Path, session: str, adm: pd.DataFrame) -> pd.DataFrame:
    """Per admitted unit, over the windows it was admitted in: the events it was offered there, those its
    refractory rule removed (within 1.5 ms of one of its June spikes or of a better new event), and the ratio
    of that count to what independent trains of the offered and the kept sizes would give."""
    import glob

    from ttm import staged_capture as sc

    wf = int(WIN_S * FS)
    blocks = pd.read_csv(RESULTS / f"stimulation_blocks_{session}.csv")
    blocks = blocks[blocks.delivered & (blocks.n_pulses >= 10)].sort_values("onset_s")
    bounds = blocks.onset_s.values - 1.0
    order = list(blocks.ttl_train)
    admitted = {(int(r.train), int(r.unit)) for r in adm[adm.admitted].itertuples(index=False)}
    units = sorted({u for _, u in admitted})
    sua = si.load(str(base_dir / "assembled_prod_sua"))
    st = si.load(str(out / "assembled_stage3"))
    group = {int(u): int(g) for u, g in zip(sua.unit_ids, sua.get_property("group"))}
    s2 = {u: np.sort(sua.get_unit_spike_train(u)) for u in units}
    s3 = {u: np.setdiff1d(st.get_unit_spike_train(u), s2[u]) for u in units}
    acc = {u: [0, 0, 0.0] for u in units}
    for f in sorted(glob.glob(str(out / "offer" / "TT*_w*.npz"))):
        tt = pathlib.Path(f).stem.split("_")[0]
        wi = int(pathlib.Path(f).stem.split("_w")[1])
        z = np.load(f)
        if z["frames"].size == 0:
            continue
        k = sc.governing_train(z["frames"] / FS, bounds)
        for u in units:
            if f"TT{group[u] + 1}" != tt:
                continue
            m = z["offered"] & (z["best"] == u)
            if not m.any():
                continue
            el = np.array([(order[kk], u) in admitted for kk in k[m]])
            n_off = int(el.sum())
            if n_off == 0:
                continue
            a, b = wi * wf, (wi + 1) * wf
            n_kept = int(((s3[u] >= a) & (s3[u] < b)).sum())
            n_june = int(((s2[u] >= a) & (s2[u] < b)).sum())
            span = wf * el.mean()
            acc[u][0] += n_off
            acc[u][1] += n_off - n_kept
            acc[u][2] += n_off * (n_june * el.mean() + n_kept) * 2 * 45 / span
    return pd.DataFrame([{"june_unit": u, "stage3_offered_admitted": v[0], "stage3_refractory_lost": v[1],
                          "stage3_conflict_ratio": round(v[1] / v[2], 3) if v[2] > 0 else None} for u, v in acc.items()])


def split_conflicts(experiment: str, params: dict) -> pd.DataFrame:
    """Per admitted unit, the in-pulse and between-pulse stage-3 counts and conflict ratios of a run made from
    saved offers with a post hoc filter, from ``staged_split.py``'s table at that run's floor and gate."""
    f = params.get("between_floor") or 3.5
    g = params.get("between_gate") or 0.8
    t = pd.read_csv(RESULTS / f"staged_split_{experiment}.csv")
    t = t[(t.arm == "split") & (t.floor == f) & (t.gate == g)]
    if params.get("inpulse_floor") or params.get("inpulse_gate"):
        t = pd.read_csv(RESULTS / f"staged_split_{experiment}.csv")
        t = t[(t.arm == "both") & (t.floor == f) & (t.gate == g)]
    p = t.pivot_table(index="unit", columns="region", values=["n_kept", "conflict_ratio", "conflict_host", "conflict_internal"])
    p.columns = [f"stage3_{a}_{b}" for a, b in p.columns]
    if params.get("pulses_only"):
        for c in [c for c in p.columns if c.endswith("_between")]:
            p[c] = 0 if c.startswith("stage3_n_kept") else np.nan
    return p.reset_index().rename(columns={"unit": "june_unit"})


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--tag", default="")
    p.add_argument("--workers", type=int, default=48)
    p.add_argument("--n-random", type=int, default=10, help="untouched staged units given a page as a sample")
    p.add_argument("--no-baseline-pages", action="store_true",
                   help="skip the baseline pages (the same units as for the run without a tag)")
    p.add_argument("--no-pages", action="store_true",
                   help="stop after the per-unit table, reusing results/purity_<E>_staged.csv when it exists")
    a = p.parse_args()
    exp = experiment_from_args(a)
    base_dir = exp.sorting_dir("mp_baseline")
    out = exp.sorting_dir(f"mp_staged{a.tag}")
    binary = str(base_dir / "binary")
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    gain = float(rec.get_property("gain_to_uV")[0])
    figdir = exp.figures_dir / "compare"
    figdir.mkdir(parents=True, exist_ok=True)

    name = f"staged{a.tag}"
    params = json.loads((out / "stage3_summary.json").read_text())["params"]
    adm = pd.read_csv(out / "admission.csv")
    if not a.tag:
        adm.to_csv(RESULTS / f"staged_admission_{exp.experiment}.csv", index=False, float_format="%.4g")
    trains = sorted(adm.train.unique())
    ever = adm[adm.admitted].groupby(["tetrode", "unit"]).size()
    units_adm = sorted(ever.index.tolist(), key=lambda x: (int(x[0][2:]), x[1]))
    mat = np.full((len(units_adm), len(trains)), np.nan)
    ratio = np.full_like(mat, np.nan)
    for i, (tt, u) in enumerate(units_adm):
        for j, t in enumerate(trains):
            r = adm[(adm.unit == u) & (adm.train == t)]
            if len(r):
                mat[i, j] = float(r.admitted.iloc[0])
                ratio[i, j] = float(r.ratio.iloc[0])
    fig, ax = plt.subplots(figsize=(14, 0.32 * len(units_adm) + 2))
    ax.imshow(np.where(np.isnan(mat), 0.5, mat), aspect="auto", cmap="RdYlGn", vmin=-0.2, vmax=1.2)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(ratio[i, j]):
                ax.text(j, i, f"{ratio[i, j]:.1f}", ha="center", va="center", fontsize=5)
    hours = adm.groupby("train").onset_h.first()
    ax.set_xticks(range(len(trains))); ax.set_xticklabels([f"{t}\n{hours[t]:.1f} h" for t in trains], fontsize=6)
    ax.set_yticks(range(len(units_adm))); ax.set_yticklabels([f"{tt} u{u}" for tt, u in units_adm], fontsize=7)
    ax.set_title("stage 3: units admitted at each delivered train (green), with the in-pulse over pre-pulse rate; "
                 "pale yellow = no template in that window", fontsize=9)
    fig.tight_layout()
    if not a.tag:
        fig.savefig(exp.figures_dir / "staged_admission.png", dpi=110)
    plt.close(fig)

    staged = si.load(str(out / "assembled_staged"))
    merges = json.loads((out / "staged_merges.json").read_text())
    id_map = {int(k): int(v) for k, v in merges["id_map"].items()}
    counts = np.load(out / "stage3_counts.npz")
    s3 = dict(zip(counts["unit_ids"].tolist(), counts["counts"].sum(axis=1).tolist()))
    tab_path = RESULTS / f"purity_{exp.experiment}_{name}.csv"
    if a.no_pages and tab_path.exists():
        tab = pd.read_csv(tab_path)
    else:
        print(f"staged: {staged.get_num_units()} units; scoring every spike ...", flush=True)
        frames, amps, rs, resid, ugroup = pu.score_sorting(binary, staged, a.workers, cache=out / "purity_scores_assembled_staged.npz")
        tab = pu.purity_table(frames, amps, ugroup, nfr, 0.0, rs=rs, resid=resid)
        mua = dict(zip([int(u) for u in staged.unit_ids], np.asarray(staged.get_property("is_mua")).astype(bool)))
        tab["is_mua"] = tab.unit.map(mua)
        tab.to_csv(tab_path, index=False)

    prod = si.load(str(base_dir / "assembled_prod"))
    base_tab = pd.read_csv(RESULTS / f"purity_{exp.experiment}_baseline.csv").set_index("unit")
    st = tab.set_index("unit")
    june_group = {int(u): int(g) for u, g in zip(prod.unit_ids, prod.get_property("group"))}
    rows = []
    for ju, su in sorted(id_map.items()):
        adm_trains = sorted(adm[(adm.unit == ju) & adm.admitted].train.tolist())
        row = {"june_unit": ju, "tetrode": f"TT{june_group[ju] + 1}", "staged_unit": su,
               "merged_with": " ".join(str(k) for k, v in id_map.items() if v == su and k != ju),
               "n_june": int(prod.get_unit_spike_train(ju).size), "n_stage3": int(s3.get(ju, 0)),
               "admitted_trains": " ".join(str(t) for t in adm_trains), "n_admitted_trains": len(adm_trains)}
        for name, t, u in (("staged", st, su), ("baseline", base_tab, ju)):
            if u in t.index:
                for col in ("rp_contamination", "ampsplit_ratio", "ampsplit_burst_ratio", "n_duplicate", "n_ambiguous", "n_distinct"):
                    row[f"{name}_{col}"] = t.loc[u, col]
        rows.append(row)
    units = pd.DataFrame(rows)
    cr = (conflict_ratios(out, base_dir, exp.session, adm) if params.get("offer_from") is None
          else split_conflicts(exp.experiment, params))
    units = units.merge(cr, on="june_unit", how="left")
    units.to_csv(RESULTS / f"staged_units_{exp.experiment}{a.tag}.csv", index=False)
    if a.no_pages:
        return

    touched_june = sorted(units[units.n_stage3 > 0].june_unit.tolist())
    floor_june = sorted(units[units.tetrode.isin(FLOOR_TETRODES)].june_unit.tolist())
    rng = np.random.default_rng(0)
    rest = sorted(set(units.june_unit) - set(touched_june) - set(floor_june))
    sample = sorted(rng.choice(rest, min(a.n_random, len(rest)), replace=False).tolist())
    page_june = sorted(set(touched_june) | set(floor_june) | set(sample))
    page_staged = sorted({id_map[u] for u in page_june})
    T, _ = pu.unit_templates(binary, {u: frames[u] for u in page_staged}, {u: ugroup[u] for u in page_staged}, nfr,
                             workers=min(a.workers, 16))
    n = pu.write_pages(figdir / f"purity_{name}.pdf", page_staged, binary, frames, amps, rs, ugroup,
                       {**{u: T[u] for u in page_staged}}, gain, nfr, 0.0, f"{exp.experiment} {name}")
    print(f"wrote {figdir / f'purity_{name}.pdf'} ({n} pages)", flush=True)
    if not a.no_baseline_pages:
        # the same units in the baseline, with their same-tetrode neighbours so the neighbour CCGs are complete
        tts = {june_group[u] for u in page_june}
        sub_ids = [int(u) for u in prod.unit_ids if june_group.get(int(u)) in tts and int(u) in id_map]
        sub = prod.select_units(sub_ids)
        bf, ba, br, bres, bg = pu.score_sorting(binary, sub, a.workers, cache=out / "purity_scores_baseline_subset.npz")
        TB, _ = pu.unit_templates(binary, {u: bf[u] for u in page_june}, {u: bg[u] for u in page_june}, nfr,
                                  workers=min(a.workers, 16))
        n = pu.write_pages(figdir / f"purity_{name}_baseline.pdf", page_june, binary, bf, ba, br, bg, TB, gain, nfr, 0.0,
                           f"{exp.experiment} baseline")
        print(f"wrote {figdir / f'purity_{name}_baseline.pdf'} ({n} pages)", flush=True)

    pres = []
    cand = si.load(str(exp.sorting_dir("mp_candidate") / "assembled_candidate"))
    for nm, srt in (("baseline", prod), ("candidate", cand), ("staged", staged)):
        d = presence(srt, nfr)
        d["sorting"] = nm
        pres.append(d)
    pres = pd.concat(pres)
    tt = pres.groupby(["tetrode", "sorting"]).presence_from_first.median().unstack()
    tt["n_units_baseline"] = pres[pres.sorting == "baseline"].groupby("tetrode").size()
    tt["n_units_candidate"] = pres[pres.sorting == "candidate"].groupby("tetrode").size()
    tt["n_units_staged"] = pres[pres.sorting == "staged"].groupby("tetrode").size()
    tt = tt.reindex(sorted(tt.index, key=lambda x: int(x[2:])))
    tt.round(3).to_csv(RESULTS / f"staged_presence_{exp.experiment}{a.tag}.csv")
    summ = {"units_touched_by_stage3": len(touched_june), "stage3_spikes": int(sum(s3.values())),
            "june_spikes_of_touched_units": int(units[units.n_stage3 > 0].n_june.sum()),
            "merges": len(merges["merges"]), "n_single_staged": int((~tab.is_mua).sum()),
            "admitted_per_train": adm[adm.admitted].groupby("train").size().reindex(trains, fill_value=0).to_dict(),
            "presence_median": pres.groupby("sorting").presence_from_first.median().round(3).to_dict()}
    s = tab[~tab.is_mua]
    summ["staged_purity"] = {"ampsplit_ratio_median": float(s.ampsplit_ratio.median()),
                             "units_ampsplit_under_0.3": int((s.ampsplit_ratio < 0.3).sum()),
                             "units_ampsplit_over_0.7": int((s.ampsplit_ratio > 0.7).sum()),
                             "rp_contamination_median": float(s.rp_contamination.median()),
                             "neighbour_verdicts": {k: int(s[f"n_{k}"].sum()) for k in ("duplicate", "ambiguous", "distinct", "segregated")}}
    (RESULTS / f"staged_compare_{exp.experiment}{a.tag}.json").write_text(json.dumps(summ, indent=1, default=int) + "\n")
    print(json.dumps(summ, indent=1, default=int))
    print(tt.round(2).to_string())


if __name__ == "__main__":
    main()
