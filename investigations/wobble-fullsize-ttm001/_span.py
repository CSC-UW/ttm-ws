"""Shared pieces of the wobble span test on TTM-001 (paths literal, as run on tononi-2).

The June carry-forward (``windowed_carry_forward_reseed`` in
``investigations/spike-sorting/_mp_common.py``) always starts at window 0
and re-sorts with MountainSort5 at each re-seed. The span test starts at
window ``SPAN[0]`` with the bank June held there and injects at
``RESEED_W`` the units June's re-seed added, so ``carry_span`` repeats the
June loop body (match, optional cosine gate, per-window re-estimation of
units with at least 100 spikes) over a list of windows. The bank and the
re-seed units are rebuilt from June's own spike trains by
``reconstruct_bank`` and ``june_reseed_units``.
"""
from __future__ import annotations

import pathlib
import sys
import time

import numpy as np
import spikeinterface as si

SORTING_DIR = pathlib.Path(__file__).resolve().parents[1] / "spike-sorting"
sys.path.insert(0, str(SORTING_DIR))

from spikeinterface.core.template_tools import get_dense_templates_array  # noqa: E402
from spikeinterface.sortingcomponents.matching import find_spikes_from_templates  # noqa: E402

from _mp_common import (  # noqa: E402
    _window_method_kwargs,
    build_templates_object,
    per_spike_fit,
    templates_from_dense,
)
from ttm.tracking import cosine_from_templates, to_int_numpy_sorting  # noqa: E402

si.set_global_job_kwargs(progress_bar=False)

FS = 30000.0
WIN_S = 1800.0
WLEN = int(WIN_S * FS)
SESSION = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-05-27_09-07-52")
JUNE = SESSION / "sortings_seed42_pcafix/track_eval/mp_long_s2000_d170000"
OUT = SESSION / "sortings_seed42_pcafix/track_eval/wobble_span"
FIG = pathlib.Path("/nvme/neuropixels/tetrode_data/TTM-001/TTM-NOD/figures/wobble_span")
SPAN_START_S = 2000.0          # the binary starts 2000 s into the recording
SPAN = list(range(55, 67))     # 12 windows, 28.1 to 34.1 h of the recording
RESEED_W = 60                  # June re-seeded after matching windows 12, 24, ..., 84
MIN_SPIKES = 100               # re-estimation floor and re-seed admission (June)
RESEED_MIN_SNR, RESEED_ADD_COS = 5.0, 0.8
AMP = {"amplitudes": [0.8, float("inf")]}
LARGE_MAD = 10.0               # template trough on its deepest wire, excised-MAD units

ARMS = {
    "circus": dict(method="circus-omp", method_kwargs=AMP, wobble_factor=None, shape_gate_r=None),
    "w35": dict(method="wobble", method_kwargs=None, wobble_factor=0.35, shape_gate_r=0.60),
    "w25": dict(method="wobble", method_kwargs=None, wobble_factor=0.25, shape_gate_r=0.60),
    "w35a": dict(method="wobble", method_kwargs=None, wobble_factor=0.35, shape_gate_r=0.60, reestimate_min_a=0.8),
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def recording():
    return si.load(JUNE / "binary")


def window_bounds(wi: int, nfr: int) -> tuple[int, int]:
    return wi * WLEN, min((wi + 1) * WLEN, nfr)


def window_counts(sorting, n_windows: int) -> tuple[np.ndarray, np.ndarray]:
    """(unit ids, counts (n_units, n_windows)) of a sorting's spikes per 1800 s window."""
    sv = sorting.to_spike_vector()
    w = sv["sample_index"] // WLEN
    cnt = np.zeros((sorting.get_num_units(), n_windows), np.int64)
    ok = w < n_windows
    np.add.at(cnt, (sv["unit_index"][ok], w[ok]), 1)
    return np.asarray(sorting.unit_ids).astype(np.int64), cnt


def _window_templates(rec, sorting, wi: int, units: list[int], n_jobs: int) -> dict:
    """June's per-window re-estimation for ``units`` (sorted ids, all present in window wi)."""
    a, b = window_bounds(wi, rec.get_num_frames())
    win = rec.frame_slice(a, b)
    win.reset_times()
    groups = dict(zip([int(u) for u in sorting.unit_ids], np.asarray(sorting.get_property("group")).astype(int)))
    samples, labels = [], []
    for u in units:
        tr = sorting.get_unit_spike_train(u).astype(np.int64)
        tr = tr[(tr >= a) & (tr < b)] - a
        samples.append(tr)
        labels.append(np.full(tr.size, u, np.int64))
    ws = si.NumpySorting.from_samples_and_labels([np.concatenate(samples)], [np.concatenate(labels)],
                                                 sampling_frequency=FS, unit_ids=np.asarray(units))
    ws.set_property("group", np.array([groups[u] for u in units]))
    new_t, _ = build_templates_object(ws, win, with_snr=False, n_jobs=n_jobs, ms_before=1.0, ms_after=2.0)
    dense = np.asarray(new_t.get_dense_templates(), dtype=np.float32)
    return {int(u): dense[k] for k, u in enumerate(new_t.unit_ids)}


def reconstruct_bank(rec, sorting, upto_w: int, *, n_jobs: int = 16) -> tuple[dict, dict, dict]:
    """The carry-forward bank after re-estimation in window ``upto_w``: each unit's template from the
    last window <= upto_w in which it had >= MIN_SPIKES spikes, computed as June's loop computed it
    (all units present in that window in one analyzer, ids ascending). Returns
    ({uid: dense (n_samp, n_chan)}, {uid: group}, {uid: source window}); units never present are omitted.
    """
    ids, cnt = window_counts(sorting, upto_w + 1)
    groups = dict(zip(ids.tolist(), np.asarray(sorting.get_property("group")).astype(int).tolist()))
    last = {}
    for k, u in enumerate(ids):
        pres = np.flatnonzero(cnt[k] >= MIN_SPIKES)
        if pres.size:
            last[int(u)] = int(pres[-1])
    bank = {}
    for w in sorted(set(last.values())):
        present = [int(u) for k, u in enumerate(ids) if cnt[k, w] >= MIN_SPIKES]
        tw = _window_templates(rec, sorting, w, present, n_jobs)
        for u, wl in last.items():
            if wl == w:
                bank[u] = tw[u]
    return bank, {u: groups[u] for u in bank}, last


def load_reseed_sort(wi: int):
    """June's MountainSort5 re-sort of window wi (cadence-12 run), aggregated as run_sorter_by_property did."""
    by_group = JUNE / "reseed_sorts_c12" / f"reseed_w{wi}" / "by_group"
    gs = sorted(int(p.name) for p in by_group.iterdir() if p.is_dir())
    from spikeinterface.sorters import read_sorter_folder
    parts = [read_sorter_folder(by_group / str(g)) for g in gs]
    agg = si.aggregate_units(parts)
    agg.set_property("group", np.concatenate([[g] * p.get_num_units() for g, p in zip(gs, parts)]))
    return to_int_numpy_sorting(agg)


def june_reseed_units(rec, bank: dict, ugroup: dict, wi: int, *, next_id: int, n_jobs: int = 16) -> list[dict]:
    """June's re-seed rule at window wi against ``bank`` (the bank after wi's re-estimation): every
    confident cluster (snr >= 5, >= 100 spikes) of the MS5 re-sort whose 4-channel template has
    shift-cosine < 0.8 to every same-tetrode bank template, in cluster order, the bank growing as
    units are added. Returns [{uid, dense, group, birth_cos, ms5_unit}]."""
    a, b = window_bounds(wi, rec.get_num_frames())
    win = rec.frame_slice(a, b)
    win.reset_times()
    rs = load_reseed_sort(wi)
    _, raz = build_templates_object(rs, win, with_snr=True, n_jobs=n_jobs, ms_before=1.0, ms_after=2.0)
    rsnr = raz.get_extension("quality_metrics").get_data()["snr"].to_numpy()
    rdense = get_dense_templates_array(raz, return_in_uV=False)
    rgrp = np.asarray(rs.get_property("group")).astype(int)
    rec_groups = np.asarray(rec.get_property("group")).astype(int)
    cur = dict(bank)
    cg = dict(ugroup)
    added = []
    for j, cu in enumerate([int(x) for x in rs.unit_ids]):
        if rsnr[j] < RESEED_MIN_SNR or len(rs.get_unit_spike_train(cu)) < MIN_SPIKES:
            continue
        g = int(rgrp[j])
        gch = np.flatnonzero(rec_groups == g)
        ct = rdense[j][:, gch]
        best = max((cosine_from_templates(ct, cur[u][:, gch], max_shift_samples=10) for u in cur if cg[u] == g),
                   default=-1.0)
        if best < RESEED_ADD_COS:
            full = np.zeros_like(next(iter(cur.values())))
            full[:, gch] = rdense[j][:, gch]
            cur[next_id] = full
            cg[next_id] = g
            added.append({"uid": next_id, "dense": full, "group": g, "birth_cos": float(best), "ms5_unit": cu,
                          "snr": float(rsnr[j])})
            next_id += 1
    return added


def carry_span(rec, bank: dict, ugroup: dict, windows: list[int], *, method: str, method_kwargs, wobble_factor,
               shape_gate_r, n_jobs: int, reseed_at: dict, raw_dir: pathlib.Path | None = None,
               reestimate_min_a: float | None = None):
    """June's carry-forward loop body over ``windows``, starting from ``bank``; ``reseed_at`` maps a
    window index to units added after that window's re-estimation. With ``raw_dir`` each window's
    ungated matcher output, per-spike (a, r) and bank are saved there. ``reestimate_min_a`` (off in
    June's loop) re-estimates templates only from kept spikes whose amplitude scale to the current
    template is at least this, and counts presence on those spikes. Returns (sorting in absolute
    frames, {uid: per-window counts}, per-window timing rows)."""
    nfr = rec.get_num_frames()
    rec_groups = np.asarray(rec.get_property("group")).astype(int)
    cur = {u: d.copy() for u, d in bank.items()}
    ug = dict(ugroup)
    nbefore = 30
    samples, labels, rows = [], [], []
    counts = {u: [] for u in cur}
    for i, wi in enumerate(windows):
        a, b = window_bounds(wi, nfr)
        win = rec.frame_slice(a, b)
        win.reset_times()
        ids_now = sorted(cur)
        ids_arr = np.asarray(ids_now)
        dense = np.stack([cur[u] for u in ids_now]).astype(np.float32)
        mask = np.stack([rec_groups == ug[u] for u in ids_now])
        templates = templates_from_dense(dense, mask, nbefore, ids_arr, rec.channel_ids)
        mk = _window_method_kwargs(method, templates, method_kwargs=method_kwargs, wobble_factor=wobble_factor)
        t0 = time.time()
        spikes = find_spikes_from_templates(win, templates, method=method, method_kwargs=mk or {},
                                            job_kwargs={"n_jobs": n_jobs, "chunk_duration": "1s",
                                                        "progress_bar": False})
        t_match = time.time() - t0
        n_raw = int(spikes.size)
        t1 = time.time()
        amp, r = per_spike_fit(spikes, templates, win) if (shape_gate_r is not None or raw_dir) else (None, None)
        t_fit = time.time() - t1
        uf = "cluster_index" if "cluster_index" in spikes.dtype.names else "unit_index"
        if raw_dir is not None:
            raw_dir.mkdir(parents=True, exist_ok=True)
            gch = [np.flatnonzero(rec_groups == ug[u]) for u in ids_now]
            np.savez(raw_dir / f"w{wi}.npz", sample_index=spikes["sample_index"].astype(np.int64),
                     cluster_index=spikes[uf].astype(np.int64), a=amp.astype(np.float32), r=r.astype(np.float32),
                     ids=ids_arr, groups=np.array([ug[u] for u in ids_now]),
                     templates4=np.stack([dense[k][:, gch[k]] for k in range(len(ids_now))]),
                     threshold=np.array(mk["parameters"]["threshold"] if method == "wobble" else np.nan))
        if shape_gate_r is not None:
            keep = r >= shape_gate_r
            spikes, amp = spikes[keep], amp[keep]
        s_idx = spikes["sample_index"].astype(np.int64)
        l_idx = spikes[uf].astype(np.int64)
        lab = ids_arr[l_idx]
        samples.append(s_idx + a)
        labels.append(lab)
        cnt = np.bincount(l_idx, minlength=len(ids_now))
        for k, u in enumerate(ids_now):
            counts[u].append(int(cnt[k]))
        t2 = time.time()
        if reestimate_min_a is None:
            r_idx, r_lab, r_cnt = s_idx, lab, cnt
        else:
            full = amp >= reestimate_min_a
            r_idx, r_lab = s_idx[full], lab[full]
            r_cnt = np.bincount(l_idx[full], minlength=len(ids_now))
        present = ids_arr[r_cnt >= MIN_SPIKES]
        if present.size:
            ws = si.NumpySorting.from_samples_and_labels([r_idx], [r_lab], sampling_frequency=FS, unit_ids=ids_arr)
            ws.set_property("group", np.array([ug[u] for u in ids_now]))
            ws = ws.select_units(present)
            new_t, _ = build_templates_object(ws, win, with_snr=False, n_jobs=n_jobs, ms_before=1.0, ms_after=2.0)
            nd = np.asarray(new_t.get_dense_templates(), dtype=np.float32)
            for k, u in enumerate([int(x) for x in new_t.unit_ids]):
                cur[u] = nd[k]
        t_reest = time.time() - t2
        for add in reseed_at.get(wi, []):
            cur[add["uid"]] = add["dense"].copy()
            ug[add["uid"]] = add["group"]
            counts[add["uid"]] = [0] * (i + 1)
        rows.append({"window": wi, "n_units": len(ids_now), "n_raw": n_raw, "n_kept": int(s_idx.size),
                     "match_s": round(t_match, 1), "fit_s": round(t_fit, 1), "reestimate_s": round(t_reest, 1),
                     "threshold": float(mk["parameters"]["threshold"]) if method == "wobble" else None})
        log(f"window {wi}: {len(ids_now)} units, {n_raw:,} matched, {s_idx.size:,} kept; match {t_match:.0f} s, "
            f"fit {t_fit:.0f} s, re-estimate {t_reest:.0f} s")
    all_ids = sorted(cur)
    sorting = si.NumpySorting.from_samples_and_labels([np.concatenate(samples)], [np.concatenate(labels)],
                                                      sampling_frequency=FS, unit_ids=np.asarray(all_ids))
    sorting.set_property("group", np.array([ug[u] for u in all_ids]))
    counts = {u: np.array(c + [0] * (len(windows) - len(c))) for u, c in counts.items()}
    return sorting, counts, rows


def restrict(sorting, a_frame: int, b_frame: int):
    """The sorting's spikes in [a_frame, b_frame), same unit ids and properties (units may be empty)."""
    trains = {}
    for u in sorting.unit_ids:
        tr = sorting.get_unit_spike_train(u).astype(np.int64)
        trains[int(u)] = tr[(tr >= a_frame) & (tr < b_frame)]
    out = si.NumpySorting.from_unit_dict([trains], sampling_frequency=FS)
    for key in sorting.get_property_keys():
        out.set_property(key, np.asarray(sorting.get_property(key)))
    return out


def amp_classes() -> dict:
    """Per stage-1 unit: tetrode, template trough (uV, excised MAD) on the deepest wire of its bank
    template at the span start (birth template for re-seed units), energy over the bank median ``q``,
    wobble's amplitude bar at threshold 0.35 ``(1 + 0.35 / q) / 2``, and large or small at LARGE_MAD."""
    rows = np.load(OUT / "template_amp_w55.npy")
    M = np.median(rows[:, 4])
    out = {}
    for u, g, tr, mad, e in rows:
        q = e / M
        out[int(u)] = {"tetrode": int(g), "template_uv": round(float(tr) * 0.195, 1), "template_mad": round(float(mad), 2),
                       "q": round(float(q), 3), "bar35": round(float((1 + 0.35 / q) / 2), 3),
                       "amp_class": "large" if mad >= LARGE_MAD else "small"}
    return out
