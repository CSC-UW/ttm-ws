"""Tagged cells and every other unit of a long sorting at each delivered laser block.

Matches the units of a sorting of the 48 h session (or of a short session,
for a check) to the OPTO2 capture templates by four-wire template cosine on
the same tetrode, so the TT6 cell (OPTO2 unit 22), the TT10 driven cell
(the threshold-4.0 re-sort's cluster 4, id 1004 in ``capture_tt10seed``),
the TT4 and TT11 candidates (OPTO2 units 3 and 39) and every other OPTO2
unit have a counterpart, then measures each unit at every delivered block:
baseline (-50 to 0 ms), onset (1 to 10 ms), in-pulse (10 ms to offset) and
post (0 to 50 and 50 to 100 ms after offset) rates, first-spike latency
and reliability, and the suppression profile on the unit's own 5 ms PSTH
as ``opto_population.py`` defines it on the pooled one (silence latency,
in-pulse fraction from 20 ms to the end, recovery after offset, rebound).
Presence is the spike count in the block and in the hour around it.

The unit templates come from the sorting's own spikes on the preprocessed
binary (up to 1000 spikes per unit from evenly spaced 60 s windows, with
single reads for sparse units), aligned at the four-wire trough like the
OPTO2 templates. Writes ``results/tagged_blocks_<E>_<name>.csv`` (one row
per unit and block), ``results/tagged_match_<E>_<name>.csv`` (one row per
unit with its OPTO2 match and, with ``--id-map``, the seed ids merged
into it), ``figures/tagged_blocks_<name>.png`` (tagged
cells' rates against block time and a presence map of every unit) and
``results/tagged_blocks_<E>_<name>.json`` (the summary printed).

Usage: python subjects/ttm-003/scripts/tagged_cells_blocks.py --subject TTM-003 --experiment TTM-003-48H
           --sorting <session>/sortings_seed42_pcafix/mp_baseline/assembled_reseed_c12 --name baseline_c12
           [--binary <session>/sortings_seed42_pcafix/mp_baseline/binary] [--templates-from DIR] [--extra-templates DIR]
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
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm import config
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
OPTO2_SESSION = "2026-10-16_14-44-22"
TAGGED = {"TT6_cell": ("TT6", 22), "TT10_driven": ("TT10", 1004), "TT4_candidate": ("TT4", 3), "TT11_candidate": ("TT11", 39)}
BIN = 0.005


def load_blocks(session: str) -> pd.DataFrame:
    path = RESULTS / f"stimulation_blocks_{session}.csv"
    if path.exists():
        df = pd.read_csv(path)
        df = df[df.delivered & (df.n_pulses >= 10)]
        df["block"] = df.ttl_train
    else:
        df = pd.read_csv(RESULTS / "stimulation_trains.csv")
        df = df[df.session == session]
        df["onset_h"] = df.onset_s / 3600
        df["block"] = df.ttl_train
    return df.reset_index(drop=True)


def reference_templates(dirs: list[pathlib.Path]) -> dict:
    """{tetrode name: (templates (n, 90, 4), ids, wires)} from capture template files; later dirs add units."""
    out = {}
    for d in dirs:
        for f in sorted(d.glob("templates_TT*.npz")):
            z = np.load(f)
            tt = f.stem.removeprefix("templates_")
            T, ids, wires = z["templates"], [int(u) for u in z["unit_ids"]], [str(w) for w in z["wires"]]
            if tt in out:
                T0, ids0, w0 = out[tt]
                assert w0 == wires, (tt, w0, wires)
                out[tt] = (np.concatenate([T0, T]), ids0 + ids, wires)
            else:
                out[tt] = (T, ids, wires)
    return out


def _window_snippets(binary: str, group: int, a: int, b: int, trains: dict, per_unit: int, pad: int) -> dict:
    rec = si.load(binary)
    x = cp.load_tetrode_traces(rec, group, a, b)
    out = {}
    rng = np.random.default_rng(a)
    for u, fr in trains.items():
        loc = fr[(fr >= a + cp.NBEFORE + pad) & (fr + cp.NAFTER + pad < b)] - a
        if loc.size == 0:
            continue
        if loc.size > per_unit:
            loc = np.sort(rng.choice(loc, per_unit, replace=False))
        out[u] = cp.snippets(x, loc, nbefore=cp.NBEFORE + pad, nafter=cp.NAFTER + pad)
    return out


def unit_templates(binary: str, trains: dict, ugroup: dict, nfr: int, *, n_windows: int = 8, win_s: float = 60.0,
                   max_spikes: int = 1000, workers: int = 16) -> tuple[dict, dict]:
    """Four-wire templates of every unit from evenly spaced windows (single reads for sparse units),
    aligned at the trough. Returns ({unit: (90, 4)}, {unit: n spikes used})."""
    pad = 10
    wf = int(win_s * FS)
    starts = np.linspace(0, nfr - wf, n_windows + 2)[1:-1].astype(int)
    groups = sorted(set(ugroup.values()))
    tasks = [(g, int(s)) for s in starts for g in groups]
    per_unit = max(max_spikes // n_windows, 20)
    res = Parallel(n_jobs=workers)(delayed(_window_snippets)(binary, g, s, s + wf, {u: trains[u] for u in trains if ugroup[u] == g},
                                                              per_unit, pad) for g, s in tasks)
    acc = {u: [] for u in trains}
    for r in res:
        for u, sn in r.items():
            acc[u].append(sn)
    rec = si.load(binary)
    templates, n_used = {}, {}
    rng = np.random.default_rng(0)
    for u in trains:
        sn = np.concatenate(acc[u]) if acc[u] else np.zeros((0, cp.N_SAMP + 2 * pad, 4), np.float32)
        if sn.shape[0] < 50 and trains[u].size:
            fr = trains[u][(trains[u] >= cp.NBEFORE + pad) & (trains[u] + cp.NAFTER + pad < nfr)]
            pick = np.sort(rng.choice(fr, min(200, fr.size), replace=False)) if fr.size else fr
            idx = cp.tetrode_channel_indices(rec, ugroup[u])
            ids = list(np.asarray(rec.channel_ids)[idx])
            extra = [rec.get_traces(start_frame=f - cp.NBEFORE - pad, end_frame=f + cp.NAFTER + pad, channel_ids=ids).astype(np.float32)
                     for f in pick]
            if extra:
                sn = np.concatenate([sn, np.stack(extra)])
        n_used[u] = int(sn.shape[0])
        if sn.shape[0] == 0:
            templates[u] = np.zeros((cp.N_SAMP, 4), np.float32)
            continue
        wide = sn.mean(axis=0)
        k = int(np.argmin(wide[pad:pad + cp.N_SAMP].min(axis=1)))
        shift = int(np.clip(k - cp.NBEFORE, -pad, pad))
        templates[u] = wide[pad + shift:pad + shift + cp.N_SAMP]
    return templates, n_used


def psth_profile(spk_s: np.ndarray, on: np.ndarray, dur: float) -> dict:
    """The opto_population.py profile on one unit's 5 ms PSTH (rates in Hz)."""
    lo, hi = -0.1, dur + 0.3
    edges = np.arange(lo, hi + BIN / 2, BIN)
    i0 = np.searchsorted(spk_s, on + lo); i1 = np.searchsorted(spk_s, on + hi)
    rel = np.concatenate([spk_s[a:b] - e for a, b, e in zip(i0, i1, on)]) if len(on) else np.zeros(0)
    r = np.histogram(rel, edges)[0] / (len(on) * BIN)
    tb = edges[:-1] + BIN / 2
    base = r[tb < 0].mean()
    half = base / 2
    below = (r < half) & (tb >= 0)
    t_sil = float(tb[below][0]) if (below.any() and base > 0) else np.nan
    t_rec = np.nan
    for i in np.where(tb > dur)[0]:
        if i + 3 <= len(r) and np.all(r[i:i + 3] >= half):
            t_rec = float(tb[i]); break
    reb = r[(tb > dur + 0.05) & (tb < dur + 0.25)].max() / base if base > 0 else np.nan
    frac = r[(tb > 0.02) & (tb < dur)].mean() / base if base > 0 else np.nan
    return {"psth_base_hz": float(base), "silence_onset_ms": None if np.isnan(t_sil) else round(t_sil * 1000, 1),
            "in_pulse_frac_20ms_to_end": None if np.isnan(frac) else round(float(frac), 3),
            "silence_after_offset_ms": None if np.isnan(t_rec) else round((t_rec - dur) * 1000, 1),
            "rebound_peak_frac": None if np.isnan(reb) else round(float(reb), 2)}


def rate_in(spk_s: np.ndarray, on: np.ndarray, lo: float, hi: float) -> float:
    i0 = np.searchsorted(spk_s, on + lo); i1 = np.searchsorted(spk_s, on + hi)
    return float((i1 - i0).sum() / (len(on) * (hi - lo)))


def block_rows(u: int, spk_s: np.ndarray, blocks: pd.DataFrame, pulses: pd.DataFrame) -> list[dict]:
    rows = []
    for blk in blocks.itertuples(index=False):
        pp = pulses[pulses.train == blk.block]
        on, off = pp.onset_s.values, pp.offset_s.values
        dur = float(np.median(off - on))
        i0 = np.searchsorted(spk_s, on); i1 = np.searchsorted(spk_s, on + 0.03)
        firsts = np.array([spk_s[a] - e for a, b, e in zip(i0, i1, on) if b > a])
        rel = np.array([(spk_s[a:b] - e) for a, b, e in zip(np.searchsorted(spk_s, on + 0.001), np.searchsorted(spk_s, on + 0.010), on)], dtype=object)
        row = {"unit": u, "block": int(blk.block), "onset_h": round(float(blk.onset_h), 4), "n_pulses": int(len(on)), "pulse_s": round(dur, 3),
               "spikes_in_block": int(np.searchsorted(spk_s, off[-1]) - np.searchsorted(spk_s, on[0])),
               "spikes_hour_around": int(np.searchsorted(spk_s, off[-1] + 1800) - np.searchsorted(spk_s, on[0] - 1800)),
               "base_hz": rate_in(spk_s, on, -0.05, 0.0), "onset_1_10ms_hz": rate_in(spk_s, on, 0.001, 0.010),
               "in_pulse_hz": rate_in(spk_s, on, 0.010, dur), "post_0_50ms_hz": rate_in(spk_s, on + dur, 0.0, 0.05),
               "post_50_100ms_hz": rate_in(spk_s, on + dur, 0.05, 0.10),
               "first_spike_ms_median": round(float(np.median(firsts)) * 1000, 2) if firsts.size else None,
               "reliability_1_10ms": float(np.mean([len(x) > 0 for x in rel])) if len(rel) else 0.0}
        row.update(psth_profile(spk_s, on, dur))
        rows.append(row)
    return rows


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--sorting", required=True, help="folder of a saved sorting with a group property")
    p.add_argument("--name", required=True, help="label for the outputs")
    p.add_argument("--binary", default=None, help="preprocessed binary the sorting's frames refer to (default: mp_baseline/binary)")
    p.add_argument("--templates-from", default=None, help="OPTO2 capture directory with templates_TT*.npz")
    p.add_argument("--extra-templates", default=None, help="a second capture directory whose templates are added (tt10seed)")
    p.add_argument("--id-map", default=None, help="candidate_merges.json whose id_map gives each unit's seed (baseline) ids")
    p.add_argument("--workers", type=int, default=16)
    a = p.parse_args()
    exp = experiment_from_args(a)
    data_root = config.get().data_root
    binary = a.binary or str(exp.sorting_dir("mp_baseline") / "binary")
    tdirs = [pathlib.Path(a.templates_from or data_root / OPTO2_SESSION / "sortings_opto" / "capture"),
             pathlib.Path(a.extra_templates or data_root / OPTO2_SESSION / "sortings_opto" / "capture_tt10seed")]
    rec = si.load(binary)
    nfr = rec.get_num_frames()
    gain = float(rec.get_property("gain_to_uV")[0])
    rec_groups = np.asarray(rec.get_property("group")).astype(int)
    ch_ids = np.asarray(rec.channel_ids)
    sorting = si.load(a.sorting)
    groups = np.asarray(sorting.get_property("group")).astype(int)
    is_mua = np.asarray(sorting.get_property("is_mua")).astype(bool) if "is_mua" in sorting.get_property_keys() else np.zeros(len(groups), bool)
    units = [int(u) for u in sorting.unit_ids]
    ugroup = dict(zip(units, groups.tolist()))
    trains = {u: np.sort(sorting.get_unit_spike_train(u).astype(np.int64)) for u in units}
    print(f"{a.name}: {len(units)} units ({int(is_mua.sum())} multi-unit), {nfr / FS / 3600:.2f} h; templates from {len(units)} units ...", flush=True)
    T, n_used = unit_templates(binary, {u: trains[u] for u in units if not is_mua[units.index(u)]}, ugroup, nfr, workers=a.workers)
    refs = reference_templates(tdirs)

    # match each unit to the OPTO2 templates of its tetrode (wire order matched by name)
    match_rows = []
    tagged_units = {}
    for u in units:
        g = ugroup[u]
        tt = f"TT{g + 1}"
        row = {"unit": u, "tetrode": tt, "is_mua": bool(is_mua[units.index(u)]), "n_spikes": int(trains[u].size),
               "rate_hz": round(trains[u].size / (nfr / FS), 3), "template_trough_uV": None, "n_template_spikes": n_used.get(u, 0),
               "opto2_unit": None, "opto2_cosine": None, "opto2_role": ""}
        if u in T and tt in refs and n_used.get(u, 0) > 0:
            row["template_trough_uV"] = round(float(T[u].min() * gain), 1)
            Tref, ids, wires = refs[tt]
            here = [str(c) for c in ch_ids[cp.tetrode_channel_indices(rec, g)]]
            order = [wires.index(w) for w in here]
            cos = np.array([cp.template_cosine(T[u], Tref[i][:, order]) for i in range(len(ids))])
            j = int(np.argmax(cos))
            row["opto2_unit"], row["opto2_cosine"] = ids[j], round(float(cos[j]), 4)
            for role, (rtt, rid) in TAGGED.items():
                if rtt == tt:
                    k = ids.index(rid) if rid in ids else None
                    if k is not None:
                        row[f"cos_{role}"] = round(float(cos[k]), 4)
        match_rows.append(row)
    match = pd.DataFrame(match_rows)
    if a.id_map:
        idm = json.loads(pathlib.Path(a.id_map).read_text())
        seeds = {}
        for seed, new in idm["id_map"].items():
            seeds.setdefault(int(new), []).append(int(seed))
        for old_mua, new in idm.get("mua_id_map", {}).items():
            seeds.setdefault(int(new), []).append(int(old_mua))
        match["seed_units"] = [" ".join(str(x) for x in sorted(seeds.get(u, []))) for u in match.unit]
    for role, (rtt, rid) in TAGGED.items():
        col = f"cos_{role}"
        if col in match:
            sub = match[(match.tetrode == rtt) & match[col].notna()]
            # The role goes to the unit whose own best OPTO2 match is the tagged unit (mutual best),
            # the highest such cosine winning; a large neighbour with a marginally higher cosine to
            # the tagged template but a better match elsewhere does not take it. Fallback: max cosine.
            mutual = sub[sub.opto2_unit == rid]
            if len(sub):
                pick = mutual if len(mutual) else sub
                best = pick.loc[pick[col].idxmax()]
                tagged_units[role] = {"unit": int(best.unit), "cosine": float(best[col]), "mutual_best": bool(len(mutual)),
                                      "runner_up": None}
                if len(sub) > 1:
                    second = sub.drop(best.name)[col].max()
                    tagged_units[role]["runner_up"] = float(second)
                match.loc[best.name, "opto2_role"] = (match.loc[best.name, "opto2_role"] + " " + role).strip()
    match.to_csv(RESULTS / f"tagged_match_{exp.experiment}_{a.name}.csv", index=False)

    blocks = load_blocks(exp.session)
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    pulses = pulses[~pulses.glitch]
    spk_s = {u: trains[u] / FS for u in units}
    rows = Parallel(n_jobs=a.workers)(delayed(block_rows)(u, spk_s[u], blocks, pulses) for u in units)
    df = pd.DataFrame([r for rr in rows for r in rr])
    df = df.merge(match[["unit", "tetrode", "is_mua", "opto2_unit", "opto2_cosine", "opto2_role"]], on="unit")
    df.round(4).to_csv(RESULTS / f"tagged_blocks_{exp.experiment}_{a.name}.csv", index=False)

    # summary
    summary = {"name": a.name, "sorting": a.sorting, "n_units": len(units), "n_blocks": int(len(blocks)), "tagged": {}}
    for role, d in tagged_units.items():
        sub = df[df.unit == d["unit"]]
        summary["tagged"][role] = d | {
            "blocks_present_10_spikes": int((sub.spikes_in_block >= 10).sum()),
            "median_base_hz": round(float(sub.base_hz.median()), 2), "median_onset_hz": round(float(sub.onset_1_10ms_hz.median()), 2),
            "median_in_pulse_hz": round(float(sub.in_pulse_hz.median()), 2), "median_post_0_50_hz": round(float(sub.post_0_50ms_hz.median()), 2),
            "in_pulse_over_base_median": round(float((sub.in_pulse_hz / sub.base_hz.replace(0, np.nan)).median()), 3),
            "blocks_in_pulse_over_2x_base": int((sub.in_pulse_hz > 2 * sub.base_hz).sum()),
            "median_reliability": round(float(sub.reliability_1_10ms.median()), 3),
            "median_first_spike_ms": None if sub.first_spike_ms_median.isna().all() else round(float(sub.first_spike_ms_median.median()), 2)}
    sua = df[~df.is_mua]
    g = sua.groupby("unit")
    ratio = (g.in_pulse_hz.median() / g.base_hz.median().replace(0, np.nan))
    summary["population"] = {
        "units_present_in_all_blocks_10_spikes": int((g.spikes_in_block.min() >= 10).sum()),
        "median_presence_fraction": round(float((g.spikes_in_block.apply(lambda s: (s >= 10).mean())).median()), 3),
        "units_suppressed_median_ratio_under_0.5": int((ratio < 0.5).sum()), "units_activated_median_ratio_over_2": int((ratio > 2).sum()),
        "units_unchanged": int(((ratio >= 0.5) & (ratio <= 2)).sum()),
        "profile_stability": {"silence_onset_ms_iqr_median": None, "in_pulse_frac_iqr_median": None}}
    sup = sua[sua.unit.isin(ratio[ratio < 0.5].index)]
    if len(sup):
        gg = sup.groupby("unit")
        summary["population"]["profile_stability"] = {
            "silence_onset_ms_iqr_median": round(float(gg.silence_onset_ms.agg(lambda s: s.quantile(0.75) - s.quantile(0.25)).median()), 1),
            "in_pulse_frac_iqr_median": round(float(gg.in_pulse_frac_20ms_to_end.agg(lambda s: s.quantile(0.75) - s.quantile(0.25)).median()), 3),
            "silence_after_offset_ms_iqr_median": round(float(gg.silence_after_offset_ms.agg(lambda s: s.quantile(0.75) - s.quantile(0.25)).median()), 1)}
    (RESULTS / f"tagged_blocks_{exp.experiment}_{a.name}.json").write_text(json.dumps(summary, indent=1, default=str) + "\n")

    # figure
    roles = list(tagged_units)
    fig, axes = plt.subplots(len(roles) + 1, 1, figsize=(16, 3 * len(roles) + 7), gridspec_kw={"height_ratios": [1] * len(roles) + [2.5]})
    axes = np.atleast_1d(axes)
    for ax, role in zip(axes, roles):
        d = tagged_units[role]
        sub = df[df.unit == d["unit"]].sort_values("onset_h")
        ax.plot(sub.onset_h, sub.base_hz, "k.-", label="baseline -50 to 0 ms")
        ax.plot(sub.onset_h, sub.onset_1_10ms_hz, "o-", color="tab:blue", label="onset 1 to 10 ms")
        ax.plot(sub.onset_h, sub.in_pulse_hz, "s-", color="tab:red", label="in pulse 10 ms to offset")
        ax.plot(sub.onset_h, sub.post_0_50ms_hz, "^-", color="tab:gray", label="0 to 50 ms after offset")
        ax.set_ylabel("spikes/s"); ax.set_title(f"{role}: unit {d['unit']} (OPTO2 template cosine {d['cosine']:.3f})", fontsize=9)
        ax.legend(fontsize=7, ncol=4)
    ax = axes[-1]
    order = sorted(units, key=lambda u: (ugroup[u], u))
    pres = df.pivot(index="unit", columns="onset_h", values="spikes_in_block").loc[order]
    im = ax.imshow(np.log10(pres.values + 1), aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(order))); ax.set_yticklabels([f"TT{ugroup[u] + 1} u{u}{'*' if is_mua[units.index(u)] else ''}" for u in order], fontsize=5)
    ax.set_xticks(range(len(pres.columns))); ax.set_xticklabels([f"{h:.1f}" for h in pres.columns], fontsize=7, rotation=90)
    ax.set_xlabel("block onset (session hour)"); ax.set_title("log10(spikes in block + 1) per unit (* = multi-unit)", fontsize=9)
    plt.colorbar(im, ax=ax)
    fig.suptitle(f"{exp.experiment} {a.name}: tagged cells and presence at {len(blocks)} delivered blocks")
    fig.tight_layout(); fig.savefig(exp.figures_dir / f"tagged_blocks_{a.name}.png", dpi=110); plt.close(fig)
    print(json.dumps(summary, indent=1, default=str))
    print(f"wrote results/tagged_blocks_{exp.experiment}_{a.name}.csv, tagged_match_..., figures/tagged_blocks_{a.name}.png")


if __name__ == "__main__":
    main()
