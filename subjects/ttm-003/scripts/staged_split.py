"""Stage-3 additions inside and between the pulses, and stricter settings between them, from saved offers.

Reads a staged run's offer files (``mp_staged<tag>/offer/``: every residual
event with its trough depth in excised MAD, ``amp_mad``, and the raw and
noise-corrected cosines to its best template, ``r`` and ``rc``, as offered at
floor 3.5 and gate 0.8) and its ``admission.csv``, and re-runs only the
acceptance (``ttm.staged_capture.accept_offered``) with stricter filters:
an event offered at 3.5 and 0.8 is kept under a setting when its depth and
noise-corrected cosine also clear that setting's floor and gate. Admissions
are the run's. Two arms per setting:

  split  in-pulse events (onset to offset of a delivered pulse) at floor 3.5
         and gate 0.8 as run, between-pulse events at the setting;
  both   the setting applied inside the pulses as well (the control arm).

For every arm, setting and admitted unit, separately for in-pulse and
between-pulse events over the time the unit was admitted: events offered,
kept, lost to the refractory rule, and the refractory conflict ratio (lost
over what independent trains of the offered and the unit's June plus kept
sizes would give in that time), and the kept count relative to floor 3.5
and gate 0.8. Per delivered block, each unit's kept events from 10 ms after
onset to offset and from 1 to 10 ms, for the tagged-cell table.

``--check`` measures how far the post hoc filter is from re-running the
offer at the stricter setting (the peeling passes subtract fewer events at
a stricter setting, so a later pass can score an event on a less clean
snippet): on a few window-tetrodes the offer is recomputed from the binary
and the two offered sets and kept counts compared.

Writes ``results/staged_split_<E><tag>.csv`` (arm, floor, gate, unit,
region), ``results/staged_split_<E><tag>_blocks.csv``,
``results/staged_split_<E><tag>_summary.json`` (median and quartiles over
the admitted units per setting, and the setting chosen by ``--target``),
``results/staged_split_<E><tag>_check.csv`` with ``--check``, and
``figures/staged_split.png`` (gate by floor panels, one dot per unit).

Usage: python subjects/ttm-003/scripts/staged_split.py --subject TTM-003 --experiment TTM-003-48H
           [--floors 3.5 4.0 4.5 5.0 5.5] [--gates 0.8 0.85 0.9 0.95] [--target 0.1] [--check] [--workers 48]
"""

from __future__ import annotations

import argparse
import itertools
import json
import pathlib
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm import staged_capture as sc
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from mp_staged import DEAD_FRAMES, MIN_SPIKES_TEMPLATE, TRAIN_MARGIN_S, delivered_trains  # noqa: E402

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
WIN_S = 1800.0
REFR = 45
BASE = (3.5, 0.8)


def intervals_mask(t: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Boolean per time: inside one of the sorted, non-overlapping intervals [starts, ends)."""
    k = np.searchsorted(starts, t, side="right") - 1
    m = np.zeros(t.size, bool)
    ok = k >= 0
    m[ok] = t[ok] < ends[k[ok]]
    return m


def overlap_s(a0: float, a1: float, starts: np.ndarray, ends: np.ndarray) -> float:
    """Total length of the intervals inside [a0, a1)."""
    return float(np.clip(np.minimum(ends, a1) - np.maximum(starts, a0), 0, None).sum())


def governed_intervals(adm: pd.DataFrame, trains: pd.DataFrame, nfr: int) -> dict:
    """{unit: (starts_s, ends_s)} of the time each unit is admitted: from a train it passed (1 s before
    onset) to the next delivered train, the first train also governing the time before it."""
    bounds = trains.onset_s.values - TRAIN_MARGIN_S
    order = list(trains.ttl_train)
    edges = np.r_[0.0, bounds[1:], nfr / FS]
    out = {}
    for u, g in adm[adm.admitted].groupby("unit"):
        ks = sorted(order.index(int(t)) for t in g.train)
        out[int(u)] = (np.array([edges[k] for k in ks]), np.array([edges[k + 1] for k in ks]))
    return out


def settings_list(floors: list, gates: list) -> list:
    """The split and both arms of every floor and gate, and always the split arm at floor 3.5 and gate 0.8
    (the run as made, which the survival fractions refer to)."""
    out = [("split", f, g) for f, g in itertools.product(floors, gates)]
    out += [("both", f, g) for f, g in itertools.product(floors, gates) if (f, g) != BASE]
    if ("split", *BASE) not in out:
        out.append(("split", *BASE))
    return out


def file_task(path: str, june: dict, settings: list, on: np.ndarray, off: np.ndarray, blocks: list,
              adm_mat: np.ndarray, all_u: np.ndarray, bounds: np.ndarray, gate_kind: str = "rc") -> dict:
    z = np.load(path)
    res = {"counts": [], "blocks": []}
    if z["frames"].size == 0:
        return res
    fr, amp, rc, r, best, offered = z["frames"], z["amp_mad"], z["rc"], z["r"], z["best"], z["offered"]
    ids = [int(u) for u in z["template_ids"]]
    ids_arr = np.asarray(ids, np.int64)
    srt = np.argsort(ids_arr)
    lab = np.where(best >= 0, srt[np.clip(np.searchsorted(ids_arr[srt], best), 0, ids_arr.size - 1)], -1)
    k = sc.governing_train(fr / FS, bounds)
    col = np.clip(np.searchsorted(all_u, np.clip(best, 0, None)), 0, all_u.size - 1)
    elig = (best >= 0) & adm_mat[k, col]
    if not (offered & elig).any():
        return res
    t = fr / FS
    inp = intervals_mask(t, on, off)
    units = [u for u in ids if (offered & elig & (best == u)).any()]
    clash = np.zeros(fr.size, bool)
    for u in units:
        mu = best == u
        clash[mu] = sc.near(fr[mu], june[u], REFR)
    blk = [(b, bon, boff) for b, bon, boff in blocks if bon[0] * FS < fr[-1] and boff[-1] * FS > fr[0]]
    for arm, f, g in settings:
        strict = (amp >= f) & ((rc if gate_kind == "rc" else r) >= g)
        if (f, g) == BASE:
            strict = (amp >= f) & (rc >= g)
        mask = offered & (strict | inp) if arm == "split" else offered & strict
        o = sc.Offer(frames=fr, amp_mad=amp, best=lab, r=r.astype(float), rc=rc.astype(float),
                     a=np.zeros(fr.size), offered=mask, pass_index=z["pass_index"])
        keep = sc.accept_offered(o, elig, {u: june[u] for u in ids}, ids, refr_frames=REFR)
        for u in units:
            mu = best == u
            for region, m in (("in_pulse", inp), ("between", ~inp)):
                res["counts"].append((arm, f, g, u, region, int((mask & elig & mu & m).sum()), int((keep & mu & m).sum()),
                                      int((mask & elig & mu & m & clash).sum())))
            for b, bon, boff in blk:
                ku = np.sort(t[keep & mu])
                n_in = int((np.searchsorted(ku, boff) - np.searchsorted(ku, bon + 0.010)).sum())
                n_on = int((np.searchsorted(ku, bon + 0.010) - np.searchsorted(ku, bon + 0.001)).sum())
                res["blocks"].append((arm, f, g, u, b, n_in, n_on))
    return res


def check_task(binary: str, group: int, wi: int, june: dict, offer_path: str, settings: list) -> list:
    """Recompute the offer of one window-tetrode at each setting and compare with the post hoc filter."""
    from threadpoolctl import threadpool_limits
    threadpool_limits(2)
    rec = si.load(binary)
    wf = int(WIN_S * FS)
    a, b = wi * wf, min((wi + 1) * wf, rec.get_num_frames())
    x = cp.load_tetrode_traces(rec, group, a, b)
    noise, _ = cp.noise_mad_excised(x)
    local = {u: fr[(fr >= a) & (fr < b)] - a for u, fr in june.items()}
    all_sorted = np.sort(np.concatenate(list(local.values())))
    bank = {u: fr for u, fr in local.items() if fr.size >= MIN_SPIKES_TEMPLATE}
    T, ids = cp.templates_from_trains(x, bank)
    lab = {u: i for i, u in enumerate(ids)}
    frames = np.concatenate([bank[u] for u in ids])
    labels = np.concatenate([np.full(bank[u].size, lab[u]) for u in ids])
    order = np.argsort(frames, kind="stable")
    residual, _ = sc.subtract_sorted(x, frames[order], labels[order], T)
    z = np.load(offer_path)
    rows = []
    for f, g in settings:
        ev, amp = sc.residual_events(residual, noise, all_sorted, threshold_mad=f, dead_frames=DEAD_FRAMES)
        off = sc.offer_events(residual, ev, amp, T, noise, gate_rc=g)
        exact = off.frames[off.offered] + a
        exact_best = np.asarray(ids)[off.best[off.offered]]
        ph = z["offered"] & (z["amp_mad"] >= f) & (z["rc"] >= g)
        post = z["frames"][ph]
        post_best = z["best"][ph]
        both = np.intersect1d(exact, post)
        same_unit = 0
        if both.size:
            e_map = dict(zip(exact.tolist(), exact_best.tolist()))
            p_map = dict(zip(post.tolist(), post_best.tolist()))
            same_unit = sum(e_map[v] == p_map[v] for v in both.tolist())
        # the refractory rule on both sets, every template treated as admitted
        lab_post = np.array([ids.index(int(v)) if v >= 0 else -1 for v in z["best"]])
        o_post = sc.Offer(frames=z["frames"] - a, amp_mad=z["amp_mad"], best=lab_post, r=z["r"].astype(float),
                          rc=z["rc"].astype(float), a=np.zeros(z["frames"].size), offered=ph, pass_index=z["pass_index"])
        k_post = sc.accept_offered(o_post, np.ones(ph.size, bool), bank, ids, refr_frames=REFR)
        k_exact = sc.accept_offered(off, np.ones(off.frames.size, bool), bank, ids, refr_frames=REFR)
        lost_post = 1 - k_post.sum() / max(ph.sum(), 1)
        lost_exact = 1 - k_exact.sum() / max(off.offered.sum(), 1)
        rows.append({"tetrode": f"TT{group + 1}", "window": wi, "floor": f, "gate": g,
                     "lost_frac_exact": round(float(lost_exact), 5), "lost_frac_post_hoc": round(float(lost_post), 5),
                     "n_exact": int(exact.size),
                     "n_post_hoc": int(post.size), "n_both": int(both.size), "n_both_same_unit": int(same_unit),
                     "exact_found_post_hoc": round(both.size / max(exact.size, 1), 4),
                     "post_hoc_found_exact": round(both.size / max(post.size, 1), 4)})
    return rows


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--tag", default="", help="the staged run whose offers and admissions are read (mp_staged<tag>)")
    p.add_argument("--floors", nargs="+", type=float, default=[3.5, 4.0, 4.5, 5.0, 5.5])
    p.add_argument("--gates", nargs="+", type=float, default=[0.8, 0.85, 0.9, 0.95])
    p.add_argument("--target", type=float, default=0.1, help="median between-pulse conflict ratio to reach")
    p.add_argument("--negligible", type=float, default=0.01,
                   help="a setting keeping under this fraction of the floor-3.5 gate-0.8 additions cannot be chosen")
    p.add_argument("--workers", type=int, default=48)
    p.add_argument("--check", action="store_true", help="also compare the post hoc filter with a recomputed offer")
    p.add_argument("--check-only", action="store_true", help="only the comparison with a recomputed offer")
    p.add_argument("--check-pairs", nargs="+", default=["TT6:0", "TT4:20", "TT10:20", "TT11:40", "TT5:20", "TT12:60"])
    p.add_argument("--out-suffix", default="", help="appended to the output names (an extended grid, for example)")
    p.add_argument("--gate-kind", choices=["rc", "raw"], default="rc",
                   help="the stricter gate on the noise-corrected cosine (default) or on the raw cosine")
    a = p.parse_args()
    exp = experiment_from_args(a)
    src = exp.sorting_dir(f"mp_staged{a.tag}")
    base_dir = exp.sorting_dir("mp_baseline")
    binary = str(base_dir / "binary")
    nfr = si.load(binary).get_num_frames()
    trains, pulses = delivered_trains(exp.session)
    pulses = pulses.sort_values("onset_s").reset_index(drop=True)
    on, off = pulses.onset_s.values, pulses.offset_s.values
    adm = pd.read_csv(src / "admission.csv")
    sua = si.load(str(base_dir / "assembled_prod_sua"))
    group = {int(u): int(g) for u, g in zip(sua.unit_ids, sua.get_property("group"))}
    june = {int(u): np.sort(sua.get_unit_spike_train(u).astype(np.int64)) for u in sua.unit_ids}
    order = list(trains.ttl_train)
    all_u = np.array(sorted(june), np.int64)
    adm_mat = np.zeros((len(order), all_u.size), bool)
    for r in adm[adm.admitted].itertuples(index=False):
        adm_mat[order.index(int(r.train)), np.searchsorted(all_u, int(r.unit))] = True
    bounds = trains.onset_s.values - TRAIN_MARGIN_S
    blocks = [(int(t), g.onset_s.values, g.offset_s.values) for t, g in pulses.groupby("train")]
    settings = settings_list(a.floors, a.gates)
    files = sorted((src / "offer").glob("TT*_w*.npz"))
    tag = a.tag + a.out_suffix
    if a.check_only:
        a.check = True
        files = []
    print(f"{len(files)} offer files, {len(settings)} arm-settings, {a.workers} workers", flush=True)
    t0 = time.perf_counter()
    res = Parallel(n_jobs=a.workers, batch_size=1)(
        delayed(file_task)(str(f), {u: june[u] for u in june if group[u] == int(f.stem.split("_")[0][2:]) - 1},
                           settings, on, off, blocks, adm_mat, all_u, bounds, a.gate_kind) for f in files)
    print(f"acceptance re-run in {(time.perf_counter() - t0) / 60:.1f} min", flush=True)
    if a.check_only:
        run_check(a, exp, binary, june, group, src, tag)
        return
    wf = int(WIN_S * FS)
    per_file = []
    for f, rr in zip(files, res):
        wi = int(f.stem.split("_w")[1])
        for arm, fl, g, u, region, n_off, n_kept, n_host in rr["counts"]:
            per_file.append((arm, fl, g, u, region, wi, n_off, n_kept, n_host))
    d = pd.DataFrame(per_file, columns=["arm", "floor", "gate", "unit", "region", "window", "n_offered", "n_kept",
                                        "n_lost_host"])
    gov = governed_intervals(adm, trains, nfr)
    # June spikes and time per unit, window and region, inside the admitted time
    jrows = []
    for u, (gs, ge) in gov.items():
        tu = june[u] / FS
        in_gov = intervals_mask(tu, gs, ge)
        in_p = intervals_mask(tu, on, off)
        for wi in sorted(d[d.unit == u].window.unique()):
            w0, w1 = wi * WIN_S, min((wi + 1) * WIN_S, nfr / FS)
            in_w = (tu >= w0) & (tu < w1)
            t_gov = sum(overlap_s(max(s, w0), min(e, w1), on, off) for s, e in zip(gs, ge) if e > w0 and s < w1)
            span = sum(max(0.0, min(e, w1) - max(s, w0)) for s, e in zip(gs, ge))
            sel = in_w & in_gov & ~in_p
            isi = np.diff(june[u][sel])
            viol = int(((isi > DEAD_FRAMES) & (isi < REFR)).sum())
            jrows.append((u, wi, "in_pulse", int((in_w & in_gov & in_p).sum()), t_gov, 0))
            jrows.append((u, wi, "between", int(sel.sum()), span - t_gov, viol))
    j = pd.DataFrame(jrows, columns=["unit", "window", "region", "n_june", "seconds", "june_viol"])
    d = d.merge(j, on=["unit", "window", "region"], how="left")
    per_frame = np.where(d.seconds > 0, 2 * REFR / (d.seconds * FS), 0.0)
    d["expected_host"] = d.n_offered * d.n_june * per_frame
    d["expected_internal"] = (d.n_offered - d.n_lost_host) * d.n_kept * per_frame
    d["expected"] = d.expected_host + d.expected_internal
    d["n_lost"] = d.n_offered - d.n_kept
    d["june_viol_expected"] = np.where(d.seconds > 0, d.n_june ** 2 * (REFR - DEAD_FRAMES) / (d.seconds * FS), 0.0)
    d["n_lost_internal"] = d.n_lost - d.n_lost_host
    agg = d.groupby(["arm", "floor", "gate", "unit", "region"]).agg(
        n_offered=("n_offered", "sum"), n_kept=("n_kept", "sum"), n_lost=("n_lost", "sum"), expected=("expected", "sum"),
        n_lost_host=("n_lost_host", "sum"), expected_host=("expected_host", "sum"),
        n_lost_internal=("n_lost_internal", "sum"), expected_internal=("expected_internal", "sum"),
        n_june=("n_june", "sum"), seconds=("seconds", "sum"), june_viol=("june_viol", "sum"),
        june_viol_expected=("june_viol_expected", "sum")).reset_index()
    agg["june_conflict_admitted"] = np.where(agg.june_viol_expected > 0, agg.june_viol / agg.june_viol_expected, np.nan)
    agg["conflict_ratio"] = np.where(agg.expected > 0, agg.n_lost / agg.expected, np.nan)
    agg["conflict_host"] = np.where(agg.expected_host > 0, agg.n_lost_host / agg.expected_host, np.nan)
    agg["conflict_internal"] = np.where(agg.expected_internal > 0, agg.n_lost_internal / agg.expected_internal, np.nan)
    base = agg[(agg.arm == "split") & (agg.floor == BASE[0]) & (agg.gate == BASE[1])].set_index(["unit", "region"]).n_kept
    agg["survival"] = [k / base.get((u, rg), np.nan) if base.get((u, rg), 0) else np.nan
                       for u, rg, k in zip(agg.unit, agg.region, agg.n_kept)]
    agg["tetrode"] = agg.unit.map(lambda u: f"TT{group[u] + 1}")
    jself = {}
    for u in agg.unit.unique():
        isi = np.diff(june[u])
        n = june[u].size
        jself[u] = float(((isi > DEAD_FRAMES) & (isi < REFR)).sum() / (n * n * 2 * (REFR - DEAD_FRAMES) / nfr / 2)) if n > 1 else np.nan
    agg["june_self_conflict"] = agg.unit.map(jself)
    agg.round(4).to_csv(RESULTS / f"staged_split_{exp.experiment}{tag}.csv", index=False)

    # blocks: kept counts per block with the June counts in the same windows, as rates
    brows = [x for rr in res for x in rr["blocks"]]
    b = pd.DataFrame(brows, columns=["arm", "floor", "gate", "unit", "block", "kept_in_pulse", "kept_onset"])
    b = b.groupby(["arm", "floor", "gate", "unit", "block"]).sum().reset_index()
    jb = []
    for blk, bon, boff in blocks:
        for u in sorted(b.unit.unique()):
            tu = june[u] / FS
            n_in = int((np.searchsorted(tu, boff) - np.searchsorted(tu, bon + 0.010)).sum())
            n_on = int((np.searchsorted(tu, bon + 0.010) - np.searchsorted(tu, bon + 0.001)).sum())
            jb.append((u, blk, n_in, n_on, len(bon), float((boff - bon - 0.010).sum())))
    jb = pd.DataFrame(jb, columns=["unit", "block", "june_in_pulse", "june_onset", "n_pulses", "in_pulse_s"])
    b = b.merge(jb, on=["unit", "block"], how="left")
    b["in_pulse_hz"] = (b.kept_in_pulse + b.june_in_pulse) / b.in_pulse_s
    b["onset_hz"] = (b.kept_onset + b.june_onset) / (b.n_pulses * 0.009)
    b["june_in_pulse_hz"] = b.june_in_pulse / b.in_pulse_s
    b.round(3).to_csv(RESULTS / f"staged_split_{exp.experiment}{tag}_blocks.csv", index=False)

    # summary over units per setting, between-pulse region of the split arm
    summ = {"settings": [], "n_units": int(agg.unit.nunique())}
    s_all = agg[(agg.arm == "split") & (agg.region == "between")]
    s = s_all[s_all.floor.isin(a.floors) & s_all.gate.isin(a.gates)]
    base_kept = s_all[(s_all.floor == BASE[0]) & (s_all.gate == BASE[1])].n_kept.sum()
    for (fl, g), q in s.groupby(["floor", "gate"]):
        cr = q.conflict_ratio.dropna()
        summ["settings"].append({
            "floor": fl, "gate": g, "conflict_median": round(float(cr.median()), 3),
            "conflict_q25": round(float(cr.quantile(0.25)), 3), "conflict_q75": round(float(cr.quantile(0.75)), 3),
            "units_over_0.2": int((cr > 0.2).sum()), "units_over_target": int((cr > a.target).sum()),
            "host_median": round(float(q.conflict_host.median()), 3),
            "internal_median": round(float(q.conflict_internal.median()), 3),
            "june_admitted_median": round(float(q.june_conflict_admitted.median()), 3),
            "between_kept": int(q.n_kept.sum()), "survival_median": round(float(q.survival.median()), 3),
            "survival_total": round(float(q.n_kept.sum() / base_kept), 4)})
    ins = agg[(agg.arm == "split") & (agg.region == "in_pulse") & (agg.floor == BASE[0]) & (agg.gate == BASE[1])]
    summ["in_pulse_as_run"] = {"conflict_median": round(float(ins.conflict_ratio.median()), 3),
                               "conflict_q25": round(float(ins.conflict_ratio.quantile(0.25)), 3),
                               "conflict_q75": round(float(ins.conflict_ratio.quantile(0.75)), 3),
                               "kept": int(ins.n_kept.sum())}
    ok = [x for x in summ["settings"] if x["conflict_median"] <= a.target and x["survival_total"] >= a.negligible]
    summ["negligible_below"] = a.negligible
    summ["target"] = a.target
    summ["gate_kind"] = a.gate_kind
    summ["chosen"] = max(ok, key=lambda x: x["between_kept"]) if ok else None
    (RESULTS / f"staged_split_{exp.experiment}{tag}_summary.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))

    # figure: gate (rows) by floor (columns), one dot per unit
    fig, axes = plt.subplots(len(a.gates), len(a.floors), figsize=(3.2 * len(a.floors), 2.8 * len(a.gates)),
                             sharex=True, sharey=True, squeeze=False)
    for i, g in enumerate(a.gates):
        for jx, fl in enumerate(a.floors):
            ax = axes[i, jx]
            q = s[(s.floor == fl) & (s.gate == g)]
            x = q.n_kept.clip(lower=1)
            y = q.conflict_ratio
            ax.scatter(x, q.conflict_internal.clip(upper=1.0), s=12, facecolors="none", edgecolors="0.55", lw=0.6)
            ax.scatter(x[y <= 1.0], y[y <= 1.0], s=12, c="k")
            ax.scatter(x[y > 1.0], np.full(int((y > 1.0).sum()), 0.98), s=18, c="k", marker="^")
            ax.set_ylim(0, 1.0)
            ax.axhline(a.target, color="tab:red", lw=0.8)
            ax.axhline(0.2, color="tab:orange", lw=0.6, ls="--")
            med = q.conflict_ratio.median()
            ax.axhline(med, color="tab:blue", lw=0.8, ls=":")
            ax.set_xscale("log")
            ax.set_title(f"floor {fl:g}, gate {g:g}: median {med:.2f}, {int(q.n_kept.sum()) / 1e6:.2f} M kept", fontsize=7)
            if i == len(a.gates) - 1:
                ax.set_xlabel("between-pulse spikes added")
            if jx == 0:
                ax.set_ylabel("conflict ratio")
    ins_med = summ["in_pulse_as_run"]["conflict_median"]
    fig.suptitle(f"{exp.experiment}: stage-3 additions between the pulses, one dot per admitted unit: conflict ratio "
                 f"(black; triangles = above 1)\nand its part among the added spikes alone (open grey). Red = target "
                 f"{a.target:g}, orange = 0.2, blue dotted = median; in-pulse additions as run: median {ins_med:.2f}",
                 fontsize=9)
    fig.tight_layout(); fig.savefig(exp.figures_dir / f"staged_split{tag}.png", dpi=110); plt.close(fig)

    if a.check:
        run_check(a, exp, binary, june, group, src, tag)


def run_check(a, exp, binary: str, june: dict, group: dict, src: pathlib.Path, tag: str) -> None:
    pairs = [(int(x.split(":")[0][2:]) - 1, int(x.split(":")[1])) for x in a.check_pairs]
    chk = [(fl, g) for fl, g in itertools.product(a.floors, a.gates)]
    tasks = [(gp, wi) for gp, wi in pairs if (src / "offer" / f"TT{gp + 1}_w{wi:03d}.npz").exists()]
    rows = Parallel(n_jobs=min(a.workers, len(tasks) * len(chk)), batch_size=1)(
        delayed(check_task)(binary, gp, wi, {u: june[u] for u in june if group[u] == gp},
                            str(src / "offer" / f"TT{gp + 1}_w{wi:03d}.npz"), [c]) for gp, wi in tasks for c in chk)
    c = pd.DataFrame([x for rr in rows for x in rr])
    c.to_csv(RESULTS / f"staged_split_{exp.experiment}{tag}_check.csv", index=False)
    print(c.groupby(["floor", "gate"])[["exact_found_post_hoc", "post_hoc_found_exact", "lost_frac_exact",
                                        "lost_frac_post_hoc"]].agg(["min", "median"]).to_string())


if __name__ == "__main__":
    main()
