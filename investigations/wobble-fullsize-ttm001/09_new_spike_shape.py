"""Do the spikes wobble adds to a unit have the unit's waveform? Every unit, not a gallery.

In window 61, per stage-1 unit and arm, the spikes are split as in
``04_stage1_compare.py`` (shared, wobble new, June's reassigned to a
neighbour by wobble). For each group with at least 50 spikes in the window,
the mean four-wire snippet (up to 2000 random spikes) is compared with June's
template of the unit in that window (mean of June's spikes): the cosine of
the two waveforms over the full window (best of +-3 sample shifts), the
cosine of their per-wire trough depths (the wire footprint), and the ratio of
their troughs on the template's deepest wire. Writes
``results/new_spike_shape.csv`` and ``results/new_spike_shape.json``.

    uv run python investigations/wobble-fullsize-ttm001/09_new_spike_shape.py
"""
import json
import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

from _span import JUNE, OUT, WLEN, amp_classes, log, recording, restrict
from ttm import capture as cp

HERE = pathlib.Path(__file__).resolve().parent
WINDOW = 61
TOL = 10


def within(q, ref, tol=TOL):
    if q.size == 0 or ref.size == 0:
        return np.zeros(q.size, bool)
    j = np.searchsorted(ref, q)
    lo = ref[np.clip(j - 1, 0, ref.size - 1)]
    hi = ref[np.clip(j, 0, ref.size - 1)]
    return np.minimum(np.abs(q - lo), np.abs(q - hi)) <= tol


def main():
    rec = recording()
    cls = amp_classes()
    a, b = WINDOW * WLEN, (WINDOW + 1) * WLEN
    june = restrict(si.load(JUNE / "assembled_reseed_c12"), a, b)
    gj = dict(zip([int(u) for u in june.unit_ids], np.asarray(june.get_property("group")).astype(int)))
    tj = {u: np.sort(june.get_unit_spike_train(u).astype(np.int64)) for u in gj}
    rng = np.random.default_rng(0)
    traces = {g: cp.load_tetrode_traces(rec, g, a, b) for g in sorted(set(gj.values()))}

    def mean_snip(g, fr):
        loc = fr[(fr >= a + cp.NBEFORE) & (fr + cp.NAFTER < b)] - a
        if loc.size > 2000:
            loc = np.sort(rng.choice(loc, 2000, replace=False))
        return cp.snippets(traces[g], loc).mean(axis=0) if loc.size else None

    rows = []
    for arm in ["circus", "w35", "w25", "w35a"]:
        s = restrict(si.load(OUT / f"span_{arm}" / "sorting"), a, b)
        tw = {int(u): np.sort(s.get_unit_spike_train(u).astype(np.int64)) for u in s.unit_ids}
        for u, g in gj.items():
            C, W = tj[u], tw.get(u, np.zeros(0, np.int64))
            if C.size < 50:
                continue
            oth_c = np.sort(np.concatenate([tj[v] for v in tj if v != u and gj[v] == g] or [np.zeros(0, np.int64)]))
            oth_w = np.sort(np.concatenate([tw[v] for v in tw if v != u and gj.get(v) == g] or [np.zeros(0, np.int64)]))
            T = mean_snip(g, C)
            groups = {"shared": W[within(W, C)], "new": W[~within(W, C) & ~within(W, oth_c)],
                      "june_reassigned": C[~within(C, W) & within(C, oth_w)]}
            k = int(np.argmax(-T.min(axis=0)))
            for name, fr in groups.items():
                if fr.size < 50:
                    continue
                M = mean_snip(g, fr)
                prof_t, prof_m = -T.min(axis=0), -M.min(axis=0)
                rows.append({"arm": arm, "unit": u, "tetrode": g + 1, "amp_class": cls[u]["amp_class"], "group": name,
                             "n_window": int(fr.size), "cosine": round(cp.template_cosine(T, M), 4),
                             "footprint_cosine": round(float(prof_t @ prof_m / (np.linalg.norm(prof_t) * np.linalg.norm(prof_m))), 4),
                             "trough_ratio": round(float(M[:, k].min() / T[:, k].min()), 3)})
        log(f"{arm}: done")
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "results" / "new_spike_shape.csv", index=False)
    out = {}
    for (arm, grp, c), d in df.groupby(["arm", "group", "amp_class"]):
        out.setdefault(arm, {}).setdefault(grp, {})[c] = {
            "units": int(len(d)), "spikes": int(d.n_window.sum()),
            "cosine_median": round(float(d.cosine.median()), 3), "units_cosine_under_0.9": int((d.cosine < 0.9).sum()),
            "footprint_cosine_median": round(float(d.footprint_cosine.median()), 3),
            "trough_ratio_median": round(float(d.trough_ratio.median()), 3)}
    (HERE / "results" / "new_spike_shape.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
