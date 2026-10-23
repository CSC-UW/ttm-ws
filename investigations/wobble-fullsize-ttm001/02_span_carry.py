"""Stage 1 over the span: June's carry-forward with circus-omp or wobble, from June's bank.

Starts at window 55 with the bank from ``01_span_bank.py``, runs the June loop
body (match, cosine gate for wobble, re-estimation of units with >= 100
spikes) through window 66 and adds June's ten re-seed units after window 60.
Arms (``_span.ARMS``): ``circus`` (amplitudes >= 0.8, the June matcher; its
output should equal June's ``assembled_reseed_c12`` on these windows),
``w35`` and ``w25`` (wobble, objective threshold 0.35 or 0.25 times the
window bank's median template energy, then cosine to the assigned template
>= 0.60) and ``w35a`` (``w35`` with templates re-estimated only from spikes at
amplitude scale >= 0.8, the population circus-omp keeps; a diagnostic, not
June's loop). Writes ``wobble_span/span_<arm>/``: the sorting, per-window counts
and timing, and per window the ungated matcher output with each spike's
amplitude scale and cosine (``raw/w<i>.npz``).

    OMP_NUM_THREADS=1 uv run python investigations/wobble-fullsize-ttm001/02_span_carry.py --arm w35
"""
import argparse
import json
import pickle
import shutil
import time

import numpy as np

from _span import ARMS, OUT, RESEED_W, SPAN, carry_span, log, recording


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", choices=sorted(ARMS), required=True)
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--max-windows", type=int, default=None, help="smoke test: first N span windows")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    cfg = ARMS[a.arm]
    windows = SPAN[: a.max_windows] if a.max_windows else SPAN
    out = OUT / f"span_{a.arm}{a.tag}"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    with open(OUT / f"bank_w{SPAN[0]}.pkl", "rb") as f:
        b = pickle.load(f)
    with open(OUT / f"reseed_w{RESEED_W}.pkl", "rb") as f:
        added = pickle.load(f)
    rec = recording()
    log(f"arm {a.arm}: {cfg}; windows {windows[0]}..{windows[-1]}; {len(b['bank'])} units; "
        f"{len(added)} re-seed units after window {RESEED_W}; n_jobs {a.n_jobs}")
    t0 = time.time()
    sorting, counts, rows = carry_span(rec, b["bank"], b["group"], windows, n_jobs=a.n_jobs,
                                       reseed_at={RESEED_W: added}, raw_dir=out / "raw", **cfg)
    wall = time.time() - t0
    sorting.save(folder=out / "sorting")
    ids = sorted(counts)
    np.savez(out / "counts.npz", ids=np.array(ids), windows=np.array(windows), counts=np.stack([counts[u] for u in ids]))
    match_s = sum(r["match_s"] for r in rows)
    fit_s = sum(r["fit_s"] for r in rows)
    reest_s = sum(r["reestimate_s"] for r in rows)
    summary = {"arm": a.arm, "config": {k: v for k, v in cfg.items() if k != "method_kwargs"},
               "method_kwargs": cfg["method_kwargs"], "windows": windows, "n_jobs": a.n_jobs,
               "wall_s": round(wall, 1), "match_s": round(match_s, 1), "fit_s": round(fit_s, 1),
               "reestimate_s": round(reest_s, 1),
               "pipeline_s": round(match_s + reest_s + (fit_s if cfg["shape_gate_r"] is not None else 0.0), 1),
               "n_spikes": int(sum(r["n_kept"] for r in rows)), "n_units": len(ids), "per_window": rows,
               "started": time.strftime("%F %T", time.localtime(t0))}
    (out / "timing.json").write_text(json.dumps(summary, indent=1))
    log(f"arm {a.arm}: {summary['n_spikes']:,} spikes in {len(ids)} units; wall {wall / 60:.1f} min "
        f"(match {match_s / 60:.1f}, fit {fit_s / 60:.1f}, re-estimate {reest_s / 60:.1f})")


if __name__ == "__main__":
    main()
