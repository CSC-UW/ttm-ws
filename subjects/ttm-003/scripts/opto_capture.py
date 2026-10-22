"""Shape-based capture of one session's spikes from the main sort's templates.

For each tetrode, ``ttm.capture`` detects every trough at ``--thresh-mad``
times the spike-excised noise, assigns it by full-window cosine to the best
same-tetrode template with the amplitude left free, peels accepted templates
off the traces so overlapping spikes get a second look, and enforces each
unit's refractory period. Templates come from the main sort
(``sortings_opto/<variant>/sorting``), or for a named tetrode from one of its
re-sorts (``--seed TT10=resort_TT10_thr4``), or from a template file written
by an earlier run of this script, so one session's templates can be carried
into the other (``--templates-from``). ``--wobble`` runs SpikeInterface's
wobble matcher on the same tetrode with the same templates and cosine gate
as the alternative detector-plus-assigner.

Writes ``<session>/sortings_opto/<out>/`` (``sorting`` with a ``group``
property, ``capture_<TT>.npz``, ``templates_<TT>.npz``, ``wobble_<TT>.npz``),
``results/capture_<experiment>_<out>.json`` (noise, counts, per-unit pulse
responses of the seed and the capture, purity numbers, timing) and
``figures/capture/<out>_<TT>_purity.pdf`` (one purity page per unit).
``--noise-check`` adds the per-block noise-versus-rate measurement that
answers open concern 2 to the JSON and a figure.

Usage: python subjects/ttm-003/scripts/opto_capture.py --subject TTM-003 --experiment TTM-003-OPTO1 --tetrodes TT6 TT10 --noise-check --wobble
"""

from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import spikeinterface as si

from ttm import capture as cp
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
RESORT_ID_BASE = 1000


def pulse_sets(session: str) -> dict:
    """Pulse onsets and widths per pattern set at the strongest power, as ``opto_tag.py`` uses."""
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{session}.csv")
    pulses = pulses[~pulses.glitch]
    trains = pd.read_csv(RESULTS / "stimulation_trains.csv")
    trains = trains[trains.session == session]
    pulses = pulses.assign(power_mW=pulses.train.map(dict(zip(trains.ttl_train, trains.power_mW))),
                           pattern=pulses.train.map(dict(zip(trains.ttl_train, trains.ttl_pattern))))
    out = {}
    for pat in ("A", "B"):
        pp = pulses[pulses.pattern == pat]
        if len(pp):
            pp = pp[pp.power_mW == pp.power_mW.max()]
            out[pat] = {"onset_s": pp.onset_s.values, "offset_s": pp.offset_s.values,
                        "dur": float(np.median(pp.offset_s - pp.onset_s)), "power_mW": float(pp.power_mW.iloc[0])}
    train_bounds = list(zip(trains.onset_s - 1.0, trains.offset_s + 1.0))
    return {"sets": out, "train_bounds": train_bounds, "all_pulses": pulses}


def rate_in(frames: np.ndarray, onsets: np.ndarray, lo: float, hi: float) -> float:
    ts = frames / FS
    i0 = np.searchsorted(ts, onsets + lo)
    i1 = np.searchsorted(ts, onsets + hi)
    return float((i1 - i0).sum() / (len(onsets) * (hi - lo)))


def pulse_response(frames: np.ndarray, ps: dict) -> dict:
    out = {}
    for name, s in ps["sets"].items():
        on, dur = s["onset_s"], s["dur"]
        out[name] = {"base_hz": rate_in(frames, on, -0.05, 0.0), "onset_1_10ms_hz": rate_in(frames, on, 0.001, 0.010),
                     "in_pulse_hz": rate_in(frames, on, 0.010, dur), "post_0_50ms_hz": rate_in(frames, on + dur, 0.0, 0.05),
                     "post_50_200ms_hz": rate_in(frames, on + dur, 0.05, 0.2)}
    return out


def in_pulse_mask(frames: np.ndarray, ps: dict) -> np.ndarray:
    ts = frames / FS
    m = np.zeros(frames.size, bool)
    pp = ps["all_pulses"]
    for a, b in zip(pp.onset_s.values, pp.offset_s.values):
        m |= (ts > a + 0.001) & (ts < b)
    return m


def outside_trains_mask(frames: np.ndarray, ps: dict) -> np.ndarray:
    ts = frames / FS
    m = np.ones(frames.size, bool)
    for a, b in ps["train_bounds"]:
        m &= ~((ts > a) & (ts < b))
    return m


def frac_matched(ref: np.ndarray, other: np.ndarray, tol: int = 15) -> float:
    if ref.size == 0 or other.size == 0:
        return 0.0
    j = np.searchsorted(other, ref)
    lo = other[np.clip(j - 1, 0, other.size - 1)]
    hi = other[np.clip(j, 0, other.size - 1)]
    return float((np.minimum(np.abs(ref - lo), np.abs(ref - hi)) <= tol).mean())


def seed_trains(a, exp, out_dir: pathlib.Path, g: int, main_sorting, main_groups) -> tuple[dict, str]:
    """Spike trains (frames) of the seed units of tetrode ``g``: the main sort, or a re-sort."""
    tt = f"TT{g + 1}"
    for spec in a.seed:
        name, folder = spec.split("=", 1)
        if name == tt:
            s = si.load(str(out_dir / folder / "sorting"))
            return {RESORT_ID_BASE + int(u): s.get_unit_spike_train(u).astype(np.int64) for u in s.unit_ids}, folder
    return {int(u): main_sorting.get_unit_spike_train(u).astype(np.int64)
            for u, gg in zip(main_sorting.unit_ids, main_groups) if gg == g}, "main"


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--variant", default="asis", help="main sort whose preprocessed recording and templates are used")
    p.add_argument("--out", default="capture", help="output variant name under sortings_opto/")
    p.add_argument("--tetrodes", nargs="+", default=["all"])
    p.add_argument("--thresh-mad", type=float, default=3.5)
    p.add_argument("--gate-r", type=float, default=0.6, help="raw cosine the best template must reach")
    p.add_argument("--gate-rc", type=float, default=None, help="noise-corrected cosine it must also reach")
    p.add_argument("--n-pass", type=int, default=3)
    p.add_argument("--margin", type=float, default=0.0, help="best cosine must beat the runner-up by this much")
    p.add_argument("--refr-ms", type=float, default=1.5)
    p.add_argument("--seed", nargs="*", default=[], metavar="TT=FOLDER",
                   help="take a tetrode's templates from this re-sort folder under the variant instead of the main sort")
    p.add_argument("--templates-from", default=None, metavar="DIR",
                   help="a capture output directory whose templates_<TT>.npz replace the seed templates (cross-session)")
    p.add_argument("--wobble", action="store_true", help="also run the wobble matcher with the same templates and gate")
    p.add_argument("--wobble-factor", type=float, default=0.35, help="wobble objective threshold as a fraction of the median template energy")
    p.add_argument("--noise-check", action="store_true", help="per-block noise versus crossing rate (open concern 2)")
    p.add_argument("--block-s", type=float, default=10.0)
    p.add_argument("--analyzer", action="store_true",
                   help="build the SortingAnalyzer of the captured sorting so opto_tag.py and opto_population.py run on --variant <out>")
    p.add_argument("--n-jobs", type=int, default=16)
    p.add_argument("--no-figures", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()
    exp = experiment_from_args(a)
    src = exp.sorting_dir(a.variant, subdir="sortings_opto")
    out = exp.sorting_dir(a.out, subdir="sortings_opto")
    if out.exists() and a.overwrite:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    figdir = exp.figures_dir / "capture"
    figdir.mkdir(parents=True, exist_ok=True)
    rec = si.load(str(src / "preprocessed"))
    main_sorting = si.load(str(src / "sorting"))
    main_groups = np.asarray(main_sorting.get_property("group")).astype(int)
    gain = float(rec.get_property("gain_to_uV")[0])
    ps = pulse_sets(exp.session)
    n_frames = rec.get_num_frames()
    refr_frames = int(round(a.refr_ms * FS / 1000))
    groups = sorted(set(main_groups.tolist())) if a.tetrodes == ["all"] else [int(t.removeprefix("TT")) - 1 for t in a.tetrodes]

    res = {"experiment": exp.experiment, "variant": a.variant, "out": a.out, "params": vars(a) | {"manifest": str(exp.manifest_path)},
           "seconds": float(n_frames / FS), "tetrodes": {}}
    all_trains, all_groups = {}, {}
    for g in groups:
        tt = f"TT{g + 1}"
        t0 = time.perf_counter()
        x = cp.load_tetrode_traces(rec, g)
        wires = [str(c) for c in np.asarray(rec.channel_ids)[cp.tetrode_channel_indices(rec, g)]]
        t_load = time.perf_counter() - t0
        trains, seed_name = seed_trains(a, exp, src, g, main_sorting, main_groups)
        if a.templates_from:
            z = np.load(pathlib.Path(a.templates_from) / f"templates_{tt}.npz")
            T, ids = z["templates"], [int(u) for u in z["unit_ids"]]
            seed_name = f"templates_from:{a.templates_from}"
        else:
            T, ids = cp.templates_from_trains(x, trains)
        np.savez(out / f"templates_{tt}.npz", templates=T, unit_ids=np.array(ids), wires=np.array(wires), seed=seed_name)
        t1 = time.perf_counter()
        C = cp.capture_tetrode(x, T, ids, threshold_mad=a.thresh_mad, gate_r=a.gate_r, gate_rc=a.gate_rc,
                               refr_frames=refr_frames, n_pass=a.n_pass, margin=a.margin)
        t_capture = time.perf_counter() - t1
        np.savez(out / f"capture_{tt}.npz", frames=C.frames, amp_mad=C.amp_mad, label=C.label, best=C.best, r=C.r,
                 rc=C.rc, a=C.a, reason=C.reason.astype(str), pass_index=C.pass_index, noise=C.noise,
                 noise_plain=C.noise_plain, unit_ids=np.array(ids))
        dead_ms = 9 / FS * 1000
        pur = cp.unit_purity(C.trains(), C.amps(), n_frames, dead_ms=dead_ms, rejected=C.rejected_by_refractory(),
                             assign_refr_frames=refr_frames)
        tet = {"wires": wires, "seed": seed_name, "noise_plain_uV": (C.noise_plain * gain).round(3).tolist(),
               "noise_excised_uV": (C.noise * gain).round(3).tolist(), "excised_fraction": round(C.excised_fraction, 4),
               "n_events": int(C.frames.size), "events_per_s": round(C.frames.size / (n_frames / FS), 1),
               "n_accepted": int((C.label >= 0).sum()), "n_pool": int(C.pool.size),
               "n_refractory_rejected": int((C.reason == "refractory").sum()),
               "accepted_by_pass": np.bincount(C.pass_index[C.pass_index >= 0], minlength=a.n_pass).tolist(),
               "seconds": {"load": round(t_load, 1), "capture": round(t_capture, 1)}, "units": {}}
        for i, u in enumerate(ids):
            cap = C.train(i)
            seed = trains.get(u, np.zeros(0, np.int64))
            tet["units"][str(u)] = {
                "template_trough_uV": round(float(T[i].min() * gain), 1),
                "template_wire_profile_uV": (T[i].min(axis=0) * gain).round(1).tolist(),
                "n_seed": int(seed.size), "n_capture": int(cap.size),
                "seed_matched_by_capture": round(frac_matched(seed, cap), 3),
                "seed_matched_by_any_event": round(frac_matched(seed, C.frames), 3),
                "capture_in_seed": round(frac_matched(cap, seed), 3),
                "seed_response": pulse_response(seed, ps), "capture_response": pulse_response(cap, ps),
                "purity": {k: v for k, v in pur[i].items() if k not in ("neighbours", "neighbour_ratios")}
                          | {"neighbours": {str(ids[v]): vv for v, vv in pur[i]["neighbours"].items()},
                             "neighbour_ratios": {str(ids[v]): (None if not np.isfinite(vv) else round(vv, 3))
                                                  for v, vv in pur[i]["neighbour_ratios"].items()}},
                "capture_amp_uV_quantiles": [round(float(v), 1) for v in np.percentile(C.amps()[i] * C.noise.max() * gain, [5, 25, 50, 75, 95])] if cap.size else None,
            }
            all_trains[u] = cap
            all_groups[u] = g
        if a.wobble:
            t2 = time.perf_counter()
            rec_tt = rec.select_channels(list(np.asarray(rec.channel_ids)[cp.tetrode_channel_indices(rec, g)]))
            wf, wl = cp.match_wobble(rec_tt, T, threshold_factor=a.wobble_factor)
            sc = cp.score_events(x, wf, T, C.noise)
            r_assigned = sc["r"][np.arange(wf.size), wl]
            keep = np.isfinite(r_assigned) & (r_assigned >= a.gate_r)
            if a.gate_rc is not None:
                keep &= sc["rc"] >= a.gate_rc
            t_wobble = time.perf_counter() - t2
            np.savez(out / f"wobble_{tt}.npz", frames=wf, label=wl, r=r_assigned, rc=sc["rc"], a=sc["a"], keep=keep)
            wtr = {i: wf[keep & (wl == i)] for i in range(len(ids))}
            wam = {i: cp.snippets(x, wtr[i]).min(axis=(1, 2)) * -1 / C.noise.max() for i in range(len(ids))}
            wpur = cp.unit_purity(wtr, wam, n_frames, dead_ms=0.0)
            tet["wobble"] = {"n_raw": int(wf.size), "n_kept": int(keep.sum()), "seconds": round(t_wobble, 1), "units": {}}
            for i, u in enumerate(ids):
                seed = trains.get(u, np.zeros(0, np.int64))
                tet["wobble"]["units"][str(u)] = {
                    "n": int(wtr[i].size), "seed_matched": round(frac_matched(seed, wtr[i]), 3),
                    "response": pulse_response(wtr[i], ps),
                    "agreement_with_capture": round(frac_matched(C.train(i), wtr[i]), 3),
                    "purity": {k: v for k, v in wpur[i].items() if k not in ("neighbours", "neighbour_ratios")}
                              | {"neighbours": {str(ids[v]): vv for v, vv in wpur[i]["neighbours"].items()}}}
        if a.noise_check:
            nb = cp.noise_by_block(x, block_frames=int(a.block_s * FS), threshold_mad=a.thresh_mad)
            corr, corr_high = {}, {}
            for name in ("plain", "excised"):
                corr[name] = [round(float(np.corrcoef(nb[name][:, c], nb["crossings_per_s"][:, c])[0, 1]), 3)
                              for c in range(4)]
                corr_high[name] = [round(float(np.corrcoef(nb[name][:, c], nb["high_crossings_per_s"][:, c])[0, 1]), 3)
                                   for c in range(4)]
            tet["noise_check"] = {"block_s": a.block_s, "corr_noise_vs_rate": corr, "corr_noise_vs_high_rate": corr_high,
                                  "high_rate_range": [nb["high_crossings_per_s"].min(0).round(1).tolist(), nb["high_crossings_per_s"].max(0).round(1).tolist()],
                                  "block_noise_plain_uV_range": [(nb["plain"] * gain).min(0).round(2).tolist(), (nb["plain"] * gain).max(0).round(2).tolist()],
                                  "block_noise_excised_uV_range": [(nb["excised"] * gain).min(0).round(2).tolist(), (nb["excised"] * gain).max(0).round(2).tolist()],
                                  "block_rate_range": [nb["crossings_per_s"].min(0).round(1).tolist(), nb["crossings_per_s"].max(0).round(1).tolist()]}
            np.savez(out / f"noise_blocks_{tt}.npz", **nb)
            if not a.no_figures:
                fig, axes = plt.subplots(1, 4, figsize=(18, 4))
                for c in range(4):
                    ax = axes[c]
                    ax.scatter(nb["high_crossings_per_s"][:, c], nb["plain"][:, c] * gain, s=6, c="k", label=f"plain r={corr_high['plain'][c]}")
                    ax.scatter(nb["high_crossings_per_s"][:, c], nb["excised"][:, c] * gain, s=6, c="tab:red", label=f"excised r={corr_high['excised'][c]}")
                    ax.set_xlabel("spikes/s (crossings > 7 MAD)"); ax.set_ylabel("noise estimate (uV)"); ax.set_title(wires[c]); ax.legend(fontsize=7)
                fig.suptitle(f"{exp.experiment} {tt}: per-{a.block_s:g} s block noise estimate against crossing rate")
                fig.tight_layout(); fig.savefig(figdir / f"{a.out}_{tt}_noise_vs_rate.png", dpi=100); plt.close(fig)
        res["tetrodes"][tt] = tet
        print(f"{tt}: {C.frames.size} events ({tet['events_per_s']}/s), accepted {tet['n_accepted']}, pool {tet['n_pool']}, "
              f"load {t_load:.0f} s, capture {t_capture:.0f} s" + (f", wobble {tet['wobble']['seconds']} s" if a.wobble else ""))
        for u in ids:
            d = tet["units"][str(u)]
            rb = d["capture_response"].get("B") or next(iter(d["capture_response"].values()))
            sb = d["seed_response"].get("B") or next(iter(d["seed_response"].values()))
            print(f"  unit {u:5d} seed n={d['n_seed']:6d} base {sb['base_hz']:5.1f} onset {sb['onset_1_10ms_hz']:6.1f} in {sb['in_pulse_hz']:6.1f} | "
                  f"capture n={d['n_capture']:6d} base {rb['base_hz']:5.1f} onset {rb['onset_1_10ms_hz']:6.1f} in {rb['in_pulse_hz']:6.1f} post {rb['post_0_50ms_hz']:5.1f} | "
                  f"seed->cap {d['seed_matched_by_capture']:.2f} conflict {d['purity']['refractory_conflict_ratio']:.2f} ampsplit {d['purity']['ampsplit_ratio']:.2f} "
                  f"dup/amb/dist {sum(v == 'duplicate' for v in d['purity']['neighbours'].values())}/{sum(v == 'ambiguous' for v in d['purity']['neighbours'].values())}/{sum(v == 'distinct' for v in d['purity']['neighbours'].values())}")
        if not a.no_figures:
            with PdfPages(figdir / f"{a.out}_{tt}_purity.pdf") as pdf:
                trains_i, amps_i, rs_i = C.trains(), C.amps(), C.rs()
                for i, u in enumerate(ids):
                    if trains_i[i].size == 0:
                        continue
                    fig, axes = plt.subplots(2, 3, figsize=(20, 10))
                    cp.plot_unit_purity(axes, x, gain, i, trains_i, amps_i, rs_i, T, list(range(len(ids))), n_frames=n_frames,
                                        highlight=in_pulse_mask(trains_i[i], ps), title=f"{exp.experiment} {tt} seed {u}:",
                                        dead_ms=dead_ms, n_rejected=C.rejected_by_refractory()[i], assign_refr_frames=refr_frames)
                    fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    # assemble the sorting (unit ids = seed ids; re-sort seeds offset by 1000)
    units = sorted(all_trains)
    sorting = si.NumpySorting.from_unit_dict([{u: np.sort(all_trains[u]) for u in units}], sampling_frequency=FS)
    sorting.set_property("group", np.array([all_groups[u] for u in units]))
    if (out / "sorting").exists():
        shutil.rmtree(out / "sorting")
    sorting.save(folder=str(out / "sorting"))
    if not (out / "preprocessed").exists():
        (out / "preprocessed").symlink_to(src / "preprocessed")
    (RESULTS / f"capture_{exp.experiment}_{a.out}.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(f"wrote {out} and results/capture_{exp.experiment}_{a.out}.json")
    if a.analyzer:
        build_analyzer(sorting, rec, out, n_jobs=a.n_jobs)


def build_analyzer(sorting, rec, out: pathlib.Path, *, n_jobs: int) -> None:
    """The analyzer ``opto_sort.py`` builds, for the captured sorting (group-sparse)."""
    from ttm.sorting import create_group_sparse_analyzer

    t0 = time.perf_counter()
    job_kwargs = dict(n_jobs=n_jobs, chunk_duration="10s", progress_bar=False)
    if (out / "analyzer").exists():
        shutil.rmtree(out / "analyzer")
    analyzer = create_group_sparse_analyzer(sorting, rec, format="binary_folder", folder=str(out / "analyzer"),
                                            n_jobs=n_jobs)
    analyzer.compute({
        "random_spikes": {"max_spikes_per_unit": 1000, "seed": 42},
        "waveforms": {"ms_before": 1.0, "ms_after": 2.0},
        "templates": {"operators": ["average", "std"]},
        "noise_levels": {},
        "spike_amplitudes": {},
        "correlograms": {"window_ms": 100.0, "bin_ms": 1.0},
        "isi_histograms": {"window_ms": 100.0, "bin_ms": 1.0},
        "template_metrics": {},
        "quality_metrics": {"metric_names": ["num_spikes", "firing_rate", "presence_ratio", "snr",
                                             "isi_violation", "rp_violation", "amplitude_cutoff",
                                             "amplitude_median", "sliding_rp_violation"]},
    }, **job_kwargs)
    print(f"analyzer in {time.perf_counter() - t0:.0f}s -> {out / 'analyzer'}")


if __name__ == "__main__":
    main()
