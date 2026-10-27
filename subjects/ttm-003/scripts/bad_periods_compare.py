"""Two sortings of the same recording, one made with bad periods masked: what each puts in flagged time, and clean-time quality.

Reads the bad-period table (``ttm-bad-periods``) and a sorting from each of
two run directories under the session's sortings (default ``mp_baseline``
and ``mp_baseline_masked``). Every spike is scored on the unmasked binary
with ``ttm.purity.score_sorting`` (trough depth, cosine to the unit's
window-local template, residual energy), then each unit's spikes are split by
whether they fall in a bad period of the unit's tetrode. Writes

- ``results/bad_periods_spikes_<E>_<tag>.csv``: per tetrode, flagged seconds and
  each sorting's single-unit and multi-unit spikes in flagged and clean time;
- ``results/bad_periods_units_<E>_<tag>.csv``: per unit of either sorting, its
  clean-time purity numbers (``ttm.purity.purity_table`` over the tetrode's
  valid time) and presence, paired with the other sorting's unit of most
  similar template on the same tetrode (``ttm.purity.match_units``, on
  clean-time spikes);
- ``results/bad_periods_compare_<E>_<tag>.json``: the medians and counts;
- ``figures/bad_periods/counts_<tag>_<tetrode>.png``: clean-time spikes per
  1800 s window of every unit of each tetrode with flagged time, both sortings;
- ``figures/bad_periods/templates_<tag>_<tetrode>.png``: each unit's mean
  waveform in the first windows, from the spikes each sorting gave it.

A unit is present in a window when it has 20 clean spikes or more there;
presence is over the windows from its first such window on in which its
tetrode has at least half the window valid. ``--max-windows N`` restricts
both sortings to the first N windows (for a carry-forward run made with
``mp_baseline.py --max-windows``).

Usage: python subjects/ttm-003/scripts/bad_periods_compare.py --subject TTM-003 --experiment TTM-003-48H
           [--a mp_baseline --b mp_baseline_masked --sorting assembled_prod] [--max-windows N]
           [--tag T] [--workers 48]
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

from ttm import bad_periods as bp
from ttm import capture as cp
from ttm import purity as pu
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
WIN_S = 1800.0
MIN_SPIKES = 20


def load_run(path: pathlib.Path, binary: str, workers: int, nfr: int, cache_tag: str):
    s = si.load(path)
    units = [int(u) for u in s.unit_ids]
    mua = (dict(zip(units, np.asarray(s.get_property("is_mua")).astype(bool)))
           if "is_mua" in s.get_property_keys() else {u: False for u in units})
    cache = path.parent / f"purity_scores_{path.name}{cache_tag}.npz"
    frames, amps, rs, resid, ugroup = pu.score_sorting(binary, s, workers, cache=cache, frame_range=(0, nfr))
    keep = {u: frames[u] < nfr for u in units}
    units = [u for u in units if keep[u].any()]
    return ({u: frames[u][keep[u]] for u in units}, {u: amps[u][keep[u]] for u in units},
            {u: rs[u][keep[u]] for u in units}, {u: resid[u][keep[u]] for u in units},
            {u: ugroup[u] for u in units}, {u: mua[u] for u in units})


def clean_split(frames, ugroup, table):
    return {u: ~bp.bad_mask(table, ugroup[u], f, FS) for u, f in frames.items()}


def window_counts(frames: np.ndarray, n_win: int) -> np.ndarray:
    return np.bincount((frames // int(WIN_S * FS)).astype(int), minlength=n_win)[:n_win]


def valid_per_window(table: pd.DataFrame, g: int, n_win: int, nfr: int) -> np.ndarray:
    out = np.zeros(n_win)
    for w in range(n_win):
        a, b = w * WIN_S, min((w + 1) * WIN_S, nfr / FS)
        sub = table[(table.group == g) & (table.end_s > a) & (table.start_s < b)].copy()
        sub["start_s"] = sub.start_s.clip(a, b) - a
        sub["end_s"] = sub.end_s.clip(a, b) - a
        out[w] = bp.valid_seconds(sub, b - a, [g])[g] / (b - a)
    return out


def presence(counts: np.ndarray, valid_frac: np.ndarray) -> tuple[float, int, int]:
    present = counts >= MIN_SPIKES
    if not present.any():
        return float("nan"), -1, -1
    first = int(np.argmax(present))
    usable = valid_frac[first:] >= 0.5
    return float(present[first:][usable].mean()) if usable.any() else float("nan"), first, int(
        np.flatnonzero(present)[-1])


def unit_table(name, frames, amps, rs, resid, ugroup, mua, clean, table, nfr, n_win):
    valid = bp.valid_seconds(table, nfr / FS, sorted(set(ugroup.values())))
    rows = []
    for g in sorted(set(ugroup.values())):
        members = [u for u in frames if ugroup[u] == g]
        fr = {u: frames[u][clean[u]] for u in members}
        tab = pu.purity_table(fr, {u: amps[u][clean[u]] for u in members}, {u: g for u in members},
                              int(valid[g] * FS), 0.0, rs={u: rs[u][clean[u]] for u in members},
                              resid={u: resid[u][clean[u]] for u in members})
        vf = valid_per_window(table, g, n_win, nfr)
        for _, r in tab.iterrows():
            u = int(r.unit)
            pres, first, last = presence(window_counts(fr[u], n_win), vf)
            rows.append({**r.to_dict(), "sorting": name, "is_mua": bool(mua[u]),
                         "rv_frac": r.rp_violations / r.n_spikes if r.n_spikes else float("nan"),
                         "n_flagged_time": int((~clean[u]).sum()), "presence": pres, "first_window": first,
                         "last_window": last})
    return pd.DataFrame(rows)


def plot_counts(path, tet, units_a, units_b, cnt_a, cnt_b, names, vf, title):
    fig, axs = plt.subplots(2, 1, figsize=(13, 6), sharex=True, sharey=True)
    for ax, units, cnt, name in zip(axs, (units_a, units_b), (cnt_a, cnt_b), names):
        for u in units:
            ax.plot(np.maximum(cnt[u], 0.5), lw=1, marker=".", label=str(u))
        ax.set_yscale("log")
        ax.set_ylabel(f"{name}\nclean spikes / window")
        ax.legend(fontsize=6, ncol=8, loc="upper right")
        bad = np.flatnonzero(vf < 0.999)
        for w in bad:
            ax.axvspan(w - 0.5, w + 0.5, color="r", alpha=0.25 * (1 - vf[w]) + 0.05, lw=0)
    axs[-1].set_xlabel("1800 s window   [red: part of the window flagged on this tetrode]")
    axs[0].set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def window_template(binary, g, frames, w, n=300, seed=0):
    rec = si.load(binary)
    a, b = int(w * WIN_S * FS), int((w + 1) * WIN_S * FS)
    f = frames[(frames >= a + cp.NBEFORE) & (frames < b - cp.NAFTER)]
    if f.size < 10:
        return None
    pick = np.sort(np.random.default_rng(seed).choice(f, min(n, f.size), replace=False))
    ids = list(np.asarray(rec.channel_ids)[cp.tetrode_channel_indices(rec, g)])
    sn = np.stack([rec.get_traces(start_frame=int(x) - cp.NBEFORE, end_frame=int(x) + cp.NAFTER, channel_ids=ids,
                                  return_in_uV=True) for x in pick])
    return sn.mean(axis=0)


def plot_templates(path, binary, g, units_a, units_b, fa, fb, names, windows, title):
    rows = [(names[0], u, fa[u]) for u in units_a] + [(names[1], u, fb[u]) for u in units_b]
    fig, axs = plt.subplots(len(rows), len(windows), figsize=(2.2 * len(windows), 1.3 * len(rows)),
                            squeeze=False, sharey="row")
    for i, (name, u, f) in enumerate(rows):
        for j, w in enumerate(windows):
            ax = axs[i, j]
            t = window_template(binary, g, f, w)
            if t is not None:
                ax.plot(np.concatenate([t[:, c] for c in range(4)]), lw=0.8)
            ax.set_xticks([])
            if j == 0:
                ax.set_ylabel(f"{name}\n{u}", fontsize=7)
            if i == 0:
                ax.set_title(f"window {w}", fontsize=8)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=90)
    plt.close(fig)


def main() -> None:
    line_buffered()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--a", default="mp_baseline")
    ap.add_argument("--a-name", default="baseline")
    ap.add_argument("--b", default="mp_baseline_masked")
    ap.add_argument("--b-name", default="masked")
    ap.add_argument("--sorting", default="assembled_prod")
    ap.add_argument("--max-windows", type=int, default=None)
    ap.add_argument("--tag", default=None, help="output tag (default <sorting>[_w<N>])")
    ap.add_argument("--table", default=None, help="bad-period table (default <session>/bad_periods/bad_periods.csv)")
    ap.add_argument("--workers", type=int, default=48)
    a = ap.parse_args()
    exp = experiment_from_args(a)
    binary = str(exp.sorting_dir("mp_baseline") / "binary")
    table = bp.read_table(a.table or exp.bad_periods_dir / "bad_periods.csv")
    nfr_all = si.load(binary).get_num_frames()
    nfr = nfr_all if a.max_windows is None else min(nfr_all, int(a.max_windows * WIN_S * FS))
    n_win = int(np.ceil(nfr / (WIN_S * FS)))
    cache_tag = "" if a.max_windows is None else f"_w{a.max_windows}"
    tag = a.tag or f"{a.sorting}{cache_tag}"
    table = table[table.start_s < nfr / FS].copy()
    table["end_s"] = table.end_s.clip(upper=nfr / FS)
    names = (a.a_name, a.b_name)
    figdir = exp.figures_dir / "bad_periods"
    figdir.mkdir(parents=True, exist_ok=True)

    runs = {}
    for name, d in zip(names, (a.a, a.b)):
        print(f"{name}: scoring {d}/{a.sorting} over {n_win} windows ...", flush=True)
        fr, am, rs, rz, ug, mua = load_run(exp.sorting_dir(d) / a.sorting, binary, a.workers, nfr, cache_tag)
        runs[name] = dict(frames=fr, amps=am, rs=rs, resid=rz, ugroup=ug, mua=mua,
                          clean=clean_split(fr, ug, table))

    gids = sorted(set().union(*[set(r["ugroup"].values()) for r in runs.values()]))
    tname = {g: f"TT{g + 1}" for g in gids}
    valid = bp.valid_seconds(table, nfr / FS, gids)
    spk = []
    for g in gids:
        row = {"tetrode": tname[g], "flagged_s": round(nfr / FS - valid[g], 1)}
        for name, r in runs.items():
            for kind, sel in (("single", False), ("mua", True)):
                us = [u for u in r["frames"] if r["ugroup"][u] == g and r["mua"][u] == sel]
                row[f"{name}_{kind}_flagged"] = int(sum((~r["clean"][u]).sum() for u in us))
                row[f"{name}_{kind}_clean"] = int(sum(r["clean"][u].sum() for u in us))
        spk.append(row)
    spk = pd.DataFrame(spk)
    spk.to_csv(RESULTS / f"bad_periods_spikes_{exp.experiment}_{tag}.csv", index=False)
    print(spk.to_string(index=False))

    tabs = {name: unit_table(name, r["frames"], r["amps"], r["rs"], r["resid"], r["ugroup"], r["mua"], r["clean"],
                             table, nfr, n_win) for name, r in runs.items()}
    temps = {}
    for name, r in runs.items():
        cf = {u: r["frames"][u][r["clean"][u]] for u in r["frames"]}
        temps[name], _ = pu.unit_templates(binary, cf, r["ugroup"], nfr, workers=min(a.workers, 16))
    sua = {name: [u for u in r["frames"] if not r["mua"][u]] for name, r in runs.items()}
    m = pu.match_units(names, tuple({u: r["frames"][u][r["clean"][u]] for u in r["frames"]} for r in runs.values()),
                       tuple(r["ugroup"] for r in runs.values()), (temps[names[0]], temps[names[1]]),
                       (sua[names[0]], sua[names[1]]))
    cols = ["n_spikes", "rate_hz", "amp_mad_median", "rp_contamination", "rv_frac", "ampsplit_ratio", "cos_median", "cos_p10",
            "resid_median", "amp_cv", "n_flagged_time", "presence", "first_window", "last_window"]
    out = m.copy()
    for name in names:
        t = tabs[name].set_index("unit")
        for c in cols + ["tetrode"]:
            out[f"{c}_{name}"] = out[name].map(lambda u, c=c, t=t: t.loc[int(u), c] if pd.notna(u) else None)
    out.to_csv(RESULTS / f"bad_periods_units_{exp.experiment}_{tag}.csv", index=False)

    summary = {"table": str(a.table or exp.bad_periods_dir / "bad_periods.csv"), "n_windows": n_win,
               "sortings": {n: str(exp.sorting_dir(d) / a.sorting) for n, d in zip(names, (a.a, a.b))},
               "flagged_s_per_tetrode": {tname[g]: round(nfr / FS - valid[g], 1) for g in gids}}
    for name in names:
        t = tabs[name][~tabs[name].is_mua]
        summary[name] = {
            "n_single_units": int(len(t)),
            "single_spikes_in_flagged_time": int(t.n_flagged_time.sum()),
            "units_mostly_in_flagged_time": int((t.n_flagged_time > t.n_spikes).sum()),
            "rp_contamination_median": float(t.rp_contamination.median()),
            "rv_frac_median": float(t.rv_frac.median()),
            "cos_median_median": float(t.cos_median.median()),
            "resid_median_median": float(t.resid_median.median()),
            "presence_median": float(t.presence.median()),
            "units_present_after_window_1": int((t.last_window > 1).sum())}
    (RESULTS / f"bad_periods_compare_{exp.experiment}_{tag}.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k in names}, indent=1))

    flagged_tets = [g for g in gids if nfr / FS - valid[g] >= 60]
    for g in flagged_tets:
        vf = valid_per_window(table, g, n_win, nfr)
        cnts, us = [], []
        for name in names:
            r = runs[name]
            u_g = [u for u in r["frames"] if r["ugroup"][u] == g and not r["mua"][u]]
            us.append(u_g)
            cnts.append({u: window_counts(r["frames"][u][r["clean"][u]], n_win) for u in u_g})
        plot_counts(figdir / f"counts_{tag}_{tname[g]}.png", tname[g], us[0], us[1], cnts[0], cnts[1], names, vf,
                    f"{exp.experiment} {tname[g]}: clean-time spikes per window, {names[0]} and {names[1]}")
        plot_templates(figdir / f"templates_{tag}_{tname[g]}.png", binary, g, us[0], us[1],
                       {u: runs[names[0]]["frames"][u] for u in us[0]},
                       {u: runs[names[1]]["frames"][u] for u in us[1]}, names, list(range(min(n_win, 4))),
                       f"{tname[g]} mean waveform per window (all its spikes, unmasked binary)")
    print(f"figures -> {figdir}", flush=True)


if __name__ == "__main__":
    main()
