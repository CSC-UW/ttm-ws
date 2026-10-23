"""June's downstream stages over the span, unchanged, for each stage-1 output.

Imports ``97_production_pipeline.py`` and redirects its output directory to
``wobble_span/down_<arm>/``, then runs its stages on the span's windows only:
the CCG-guarded merge (cosine >= 0.90 proposes; CCG duplicate or cosine
>= 0.95 accepts), the residual single-unit pass (unclaimed events >= 10 MAD
offered to the best same-tetrode template at cosine >= 0.8, kept if they do
not break its 1.5 ms refractory period), the per-tetrode multi-unit bucket
(>= 7 MAD, cosine 0.55 to 0.8) and the three-axis scoreboard. June's
residual pass contains no matcher, so "wobble in the residual pass" is the
same pass run on wobble's stage-1 output.

The stage-1 input is ``span_<arm>/sorting`` for the wobble arms. For
``june`` it is June's ``assembled_reseed_c12`` cut to the span, and
``circus`` is the circus re-run of ``02_span_carry.py`` (the reproduction
check). Units without spikes in the span are dropped before the merge. The
detected events are June's ``spike_coverage.npz`` cut to the span. The merge
records are kept (``merges.json``) so units can be followed through it.

    OMP_NUM_THREADS=1 uv run python investigations/wobble-fullsize-ttm001/03_span_downstream.py --arm w35
"""
import argparse
import importlib.util
import json
import shutil
import time

import numpy as np
import spikeinterface as si

from _span import FS, JUNE, OUT, SORTING_DIR, SPAN, WLEN, log, recording, restrict


def load_p97():
    spec = importlib.util.spec_from_file_location("p97", SORTING_DIR / "97_production_pipeline.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def span_peaks(a0, b0):
    z = np.load(JUNE / "spike_coverage.npz", mmap_mode="r")
    ps = z["peak_sample"]
    lo, hi = np.searchsorted(ps, [a0, b0])
    return (np.asarray(ps[lo:hi]), np.asarray(z["peak_group"][lo:hi]), np.asarray(z["amp_mad"][lo:hi]))


def stage1_input(arm, a0, b0):
    if arm == "june":
        return restrict(si.load(JUNE / "assembled_reseed_c12"), a0, b0)
    return si.load(OUT / f"span_{arm}" / "sorting")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", required=True, help="june, circus, w35 or w25")
    ap.add_argument("--n-jobs", type=int, default=16)
    a = ap.parse_args()
    out = OUT / f"down_{a.arm}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    p97 = load_p97()
    p97.OUT = out
    p97.N_JOBS = a.n_jobs
    merges_seen = []
    original_merge = p97.ccg_guarded_merge

    def recording_merge(*args, **kw):
        merged, merges = original_merge(*args, **kw)
        merges_seen.extend(merges)
        return merged, merges

    p97.ccg_guarded_merge = recording_merge

    rec = recording()
    a0, b0 = SPAN[0] * WLEN, (SPAN[-1] + 1) * WLEN
    starts = [w * WLEN for w in SPAN]
    peaks = span_peaks(a0, b0)
    base = stage1_input(a.arm, a0, b0)
    n_spk = {int(u): int(base.get_unit_spike_train(u).size) for u in base.unit_ids}
    keep = [u for u in base.unit_ids if n_spk[int(u)] > 0]
    dropped = [int(u) for u in base.unit_ids if n_spk[int(u)] == 0]
    base = base.select_units(keep)
    base.save(folder=out / "assembled_base")
    log(f"{a.arm}: stage 1 has {len(keep)} units with spikes in windows {SPAN[0]}..{SPAN[-1]} "
        f"({sum(n_spk.values()):,} spikes); dropped without spikes: {dropped}; {len(peaks[0]):,} detected events")
    timing = {}
    t = time.time()
    merged = p97.stage_merge(rec, "assembled_base", "assembled_prod_merge", force=True)
    timing["merge_s"] = round(time.time() - t, 1)
    (out / "merges.json").write_text(json.dumps(
        [{"a": int(x), "b": int(y), "verdict": v, "cosine": round(c, 4)} for x, y, v, c in merges_seen], indent=1))
    t = time.time()
    sua = p97.stage_residual(rec, merged, peaks, starts, "assembled_prod_sua", force=True)
    timing["residual_s"] = round(time.time() - t, 1)
    t = time.time()
    prod = p97.stage_mua(rec, sua, peaks, starts, p97.MUA_THETA_DEFAULT, "assembled_prod", force=True)
    timing["mua_s"] = round(time.time() - t, 1)
    t = time.time()
    p97.stage_score(rec, [("base", base, True, True), ("prod_sua", sua, True, True),
                          ("prod_sua+mua", prod, False, False)], peaks, starts, tag="")
    timing["score_s"] = round(time.time() - t, 1)
    timing.update(arm=a.arm, n_base_units=len(keep), dropped_units=dropped, n_merged=merged.get_num_units(),
                  n_prod=prod.get_num_units(), span_frames=[int(a0), int(b0)], span_hours=(b0 - a0) / FS / 3600)
    (out / "timing.json").write_text(json.dumps(timing, indent=1))
    log(f"{a.arm}: done {timing}")


if __name__ == "__main__":
    main()
