"""Purity numbers and pictures for every unit of a long sorting, and the disagreements with another.

The scoring, table, templates, pairing and pages are ``ttm.purity``: per
1800 s window and tetrode the preprocessed binary is read once, the
spike-excised noise estimated, and every spike of every unit of the
tetrode gets its trough depth in MAD units, its cosine to the unit's
window-local template and its residual energy in noise units. From the
whole trains ``ttm.capture.unit_purity`` then gives per unit the refractory
contamination, the amplitude-split zero-lag ratio and 2 to 3 ms burst
ratio, and the cross-correlogram verdict against every same-tetrode unit
(``results/purity_<E>_<name>.csv``, with the within-cluster spread columns).
``figures/compare/purity_<name>.pdf`` holds one page per unit for a random
sample of ``--n-pages`` units and, with ``--other``, every unit on which the
two sortings disagree: units are paired by template cosine on the same
tetrode and the agreement is the fraction of the smaller train's spikes
within 0.5 ms of a spike of the other (``results/purity_match_<E>_<name>_vs_<other-name>.csv``);
pairs under ``--agree`` and unpaired units get a page in both files. The
per-spike scores are cached next to each sorting
(``purity_scores_<sorting>.npz``) so a rerun only redraws; scatter pages
show at most 20,000 points, rasterised.

Usage: python subjects/ttm-003/scripts/sorting_purity.py --subject TTM-003 --experiment TTM-003-48H
           --sorting <dir>/assembled_candidate --name candidate [--dead-ms 0.3]
           [--other <dir>/assembled_prod --other-name baseline] [--n-pages 20] [--workers 48]
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import numpy as np
import spikeinterface as si

from ttm import purity as pu
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"


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
        cache = pathlib.Path(path).parent / f"purity_scores_{pathlib.Path(path).name}.npz"
        frames, amps, rs, resid, ugroup = pu.score_sorting(binary, s, a.workers, cache=cache)
        tab = pu.purity_table(frames, amps, ugroup, nfr, dead, rs=rs, resid=resid)
        if "is_mua" in s.get_property_keys():
            mua = dict(zip([int(u) for u in s.unit_ids], np.asarray(s.get_property("is_mua")).astype(bool)))
            tab["is_mua"] = tab.unit.map(mua)
        tab.to_csv(RESULTS / f"purity_{exp.experiment}_{name}.csv", index=False)
        T, _ = pu.unit_templates(binary, frames, ugroup, nfr, workers=min(a.workers, 16))
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
        m = pu.match_units((a.name, a.other_name), (f1, f2), (g1, g2), (T1, T2), (sua1, sua2))
        m.to_csv(RESULTS / f"purity_match_{exp.experiment}_{a.name}_vs_{a.other_name}.csv", index=False)
        dis = m[m.agreement < a.agree]
        pages[a.name] |= {int(u) for u in dis[a.name].dropna()}
        pages[a.other_name] |= {int(v) for v in dis[a.other_name].dropna()}
        sample2 = rng.choice(sua2, min(a.n_pages, len(sua2)), replace=False).tolist()
        pages[a.other_name] |= set(sample2)
        print(f"pairs: {len(m)}; agreement under {a.agree}: {len(dis)} (median agreement {m.agreement.median():.3f})", flush=True)
    for name, (fr, am, rr, ug, T, dead) in {a.name: (f1, a1, r1, g1, T1, a.dead_ms),
                                            **({a.other_name: (f2, a2, r2, g2, T2, a.other_dead_ms)} if a.other else {})}.items():
        pu.write_pages(figdir / f"purity_{name}.pdf", pages[name], binary, fr, am, rr, ug, T, gain, nfr, dead,
                       f"{exp.experiment} {name}")
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
