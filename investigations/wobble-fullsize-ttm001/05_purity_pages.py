"""Per-unit purity pages and tables for the single-unit sortings of each arm over the span.

For each arm's ``down_<arm>/assembled_prod_sua`` (merge plus residual pass)
``ttm.purity`` scores every spike over windows 55 to 66 and gives per unit
the refractory contamination, amplitude-split ratios, the CCG verdict against
every same-tetrode unit and the within-cluster spread (cosine median and
10th percentile, residual energy, trough-amplitude CV), with rates over the
span's 6 h. The June assignment purity (``_scoreboard.windowed_axis_b``:
fraction of a unit's spikes whose best same-tetrode template is its own,
full window and the tight trough window) is added per unit.

Units of the two arms pair by their stage-1 members: the merge keeps a record
of which stage-1 ids each unit holds (``merges.json``), and stage-1 ids are
shared, so a unit pairs with the unit of the other arm holding the same
members, or failing that the most of them. Each unit's amplitude class is
that of its member with the most spikes (``_span.amp_classes``).
``ttm.purity.match_units`` (template cosine and spike agreement, as the
TTM-003 pages pair) is written alongside.

Writes ``results/units_<arm>.csv``, ``results/pairs_june_<arm>.csv``,
``results/match_june_<arm>.csv`` and the PDFs ``purity_<arm>.pdf`` (one page
per single unit) under the TTM-NOD figures directory.

    OMP_NUM_THREADS=1 uv run python investigations/wobble-fullsize-ttm001/05_purity_pages.py
"""
import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import spikeinterface as si

from _span import FIG, JUNE, OUT, SPAN, WLEN, amp_classes, log, recording
from ttm import purity as pu


HERE = pathlib.Path(__file__).resolve().parent
RESULTS = HERE / "results"
ARMS = ["june", "w35", "w25", "w35a"]


def components(base_ids, merges):
    """Merged unit id -> sorted stage-1 member ids, reproducing ccg_guarded_merge's numbering
    (components ordered by their smallest member)."""
    parent = {u: u for u in base_ids}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    for m in merges:
        ra, rb = find(m["a"]), find(m["b"])
        if ra != rb:
            parent[rb] = ra
    comp = {}
    for u in sorted(base_ids):
        comp.setdefault(find(u), []).append(u)
    return {k: sorted(v) for k, v in enumerate(comp.values())}


def axis_b_per_unit(rec, sorting, cache):
    if cache.exists():
        d = json.loads(cache.read_text())
        return ({int(k): v for k, v in d["full"].items()}, {int(k): v for k, v in d["tight"].items()})
    import sys
    from _span import SORTING_DIR
    sys.path.insert(0, str(SORTING_DIR))
    from _scoreboard import windowed_axis_b
    full, tight = windowed_axis_b(rec, sorting, [w * 0.5 for w in SPAN], win_s=1800.0, n_jobs=16,
                                  min_spikes_template=100)
    cache.write_text(json.dumps({"full": full, "tight": tight}, default=float))
    return full, tight


def arm_tables(arm, rec, binary, cls, span):
    d = OUT / f"down_{arm}"
    base = si.load(d / "assembled_base")
    sua = si.load(d / "assembled_prod_sua")
    merges = json.loads((d / "merges.json").read_text())
    comp = components([int(u) for u in base.unit_ids], merges)
    assert len(comp) == sua.get_num_units(), (arm, len(comp), sua.get_num_units())
    nbase = {int(u): int(base.get_unit_spike_train(u).size) for u in base.unit_ids}
    frames, amps, rs, resid, ugroup = pu.score_sorting(binary, sua, 16, cache=d / "purity_scores_sua.npz",
                                                       frame_range=span)
    tab = pu.purity_table(frames, amps, ugroup, span[1] - span[0], 0.0, rs=rs, resid=resid)
    full, tight = axis_b_per_unit(rec, sua, d / "axis_b_units.json")
    tab["members"] = tab.unit.map(lambda u: " ".join(str(m) for m in comp[u]))
    lead = {u: max(comp[u], key=lambda m: nbase.get(m, 0)) for u in comp}
    tab["lead_member"] = tab.unit.map(lead)
    for k in ("template_mad", "q", "bar35", "amp_class"):
        tab[k] = tab.lead_member.map(lambda m, k=k: cls.get(m, {}).get(k))
    tab["purity_full"] = tab.unit.map(lambda u: full.get(u, {}).get("best_match_frac"))
    tab["purity_tight"] = tab.unit.map(lambda u: tight.get(u, {}).get("best_match_frac"))
    tab["top_neighbour_tight"] = tab.unit.map(lambda u: tight.get(u, {}).get("top_neighbor"))
    tab.to_csv(RESULTS / f"units_{arm}.csv", index=False)
    T, _ = pu.unit_templates(binary, frames, ugroup, rec.get_num_frames(), workers=16, frame_range=span)
    return {"tab": tab, "comp": comp, "frames": frames, "amps": amps, "rs": rs, "ugroup": ugroup, "T": T}


def pair_by_members(A, B, nameA, nameB):
    rows, used = [], set()
    setsB = {v: set(m) for v, m in B["comp"].items()}
    for u, mem in A["comp"].items():
        mu = set(mem)
        best, score = None, 0.0
        for v, mv in setsB.items():
            j = len(mu & mv) / len(mu | mv)
            if j > score:
                best, score = v, j
        if best is not None:
            used.add(best)
        fa = A["frames"][u]
        fb = B["frames"][best] if best is not None else np.zeros(0, np.int64)
        rows.append({nameA: u, nameB: best, "members_" + nameA: " ".join(map(str, mem)),
                     "members_" + nameB: " ".join(map(str, B["comp"][best])) if best is not None else "",
                     "member_jaccard": round(score, 3), "same_members": score == 1.0,
                     "n_" + nameA: int(fa.size), "n_" + nameB: int(fb.size),
                     f"in_{nameB}": round(pu.frac_within(fa, fb), 4), f"in_{nameA}": round(pu.frac_within(fb, fa), 4)})
    for v in B["comp"]:
        if v not in used:
            rows.append({nameA: None, nameB: v, "members_" + nameB: " ".join(map(str, B["comp"][v])),
                         "member_jaccard": 0.0, "same_members": False, "n_" + nameB: int(B["frames"][v].size)})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arms", nargs="+", default=ARMS)
    ap.add_argument("--no-pages", action="store_true")
    a = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    rec = recording()
    binary = str(JUNE / "binary")
    span = (SPAN[0] * WLEN, (SPAN[-1] + 1) * WLEN)
    gain = float(rec.get_property("gain_to_uV")[0])
    cls = amp_classes()
    res = {}
    for arm in a.arms:
        if not (OUT / f"down_{arm}" / "assembled_prod_sua").exists():
            log(f"{arm}: no downstream output yet, skipped")
            continue
        res[arm] = arm_tables(arm, rec, binary, cls, span)
        r = res[arm]
        if a.no_pages:
            continue
        n = pu.write_pages(FIG / f"purity_{arm}.pdf", list(r["frames"]), binary, r["frames"], r["amps"], r["rs"],
                           r["ugroup"], r["T"], gain, rec.get_num_frames(), 0.0,
                           f"TTM-001 span 55-66 {arm}", duration_frames=span[1] - span[0])
        log(f"{arm}: {len(r['tab'])} single units, {n} pages -> {FIG / f'purity_{arm}.pdf'}")
    for arm in [x for x in a.arms if x != "june" and x in res]:
        pair_by_members(res["june"], res[arm], "june", arm).to_csv(RESULTS / f"pairs_june_{arm}.csv", index=False)
        A, B = res["june"], res[arm]
        m = pu.match_units(("june", arm), (A["frames"], B["frames"]), (A["ugroup"], B["ugroup"]), (A["T"], B["T"]),
                           (list(A["frames"]), list(B["frames"])))
        m.to_csv(RESULTS / f"match_june_{arm}.csv", index=False)
        log(f"june vs {arm}: median template-pair agreement {m.agreement.median():.3f}")


if __name__ == "__main__":
    main()
