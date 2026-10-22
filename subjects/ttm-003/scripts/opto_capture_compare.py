"""Compare two sessions' captured units: match by template, compare pulse signatures.

Reads ``results/capture_<experiment>_<out>.json`` and the capture output
directories of two experiments. For every tetrode, units of the first
session are matched to units of the second by the maximum-shift cosine of
their four-wire templates (greedy, best pair first, ``--min-cos``), and the
matched pairs' pulse responses (base, onset, in-pulse and post rates for
the short-pulse set, or whichever set both sessions share) and
``opto_tag.py`` classes (when ``optotag_<experiment>_<out>.csv`` exists)
are tabulated. Writes ``results/capture_xsession_<A>_<B>_<out>.csv`` and
prints the table.

Usage: python subjects/ttm-003/scripts/opto_capture_compare.py --experiments TTM-003-OPTO1 TTM-003-OPTO2 --out capture
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import pandas as pd

from ttm import capture as cp
from ttm.subject import find_subject

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"


def load(exp_name: str, out: str):
    subject = find_subject("TTM-003")
    exp = subject.experiment(exp_name)
    d = exp.sorting_dir(out, subdir="sortings_opto")
    res = json.loads((RESULTS / f"capture_{exp_name}_{out}.json").read_text())
    templates = {}
    for tt in res["tetrodes"]:
        z = np.load(d / f"templates_{tt}.npz")
        templates[tt] = (z["templates"], [int(u) for u in z["unit_ids"]])
    tag_path = RESULTS / f"optotag_{exp_name}_{out}.csv"
    tag = pd.read_csv(tag_path) if tag_path.exists() else None
    return res, templates, tag


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--experiments", nargs=2, required=True)
    p.add_argument("--out", default="capture")
    p.add_argument("--min-cos", type=float, default=0.9)
    p.add_argument("--set", default="B", help="pulse set whose responses are compared")
    a = p.parse_args()
    (ra, ta, tag_a), (rb, tb, tag_b) = load(a.experiments[0], a.out), load(a.experiments[1], a.out)
    rows = []
    for tt in ra["tetrodes"]:
        if tt not in rb["tetrodes"]:
            continue
        TA, ia = ta[tt]
        TB, ib = tb[tt]
        cos = np.array([[cp.template_cosine(TA[i], TB[j], max_shift=5) for j in range(len(ib))] for i in range(len(ia))])
        used_a, used_b = set(), set()
        for flat in np.argsort(cos, axis=None)[::-1]:
            i, j = np.unravel_index(flat, cos.shape)
            if cos[i, j] < a.min_cos:
                break
            if i in used_a or j in used_b:
                continue
            used_a.add(i); used_b.add(j)
            ua, ub = ia[i], ib[j]
            da = ra["tetrodes"][tt]["units"][str(ua)]
            db = rb["tetrodes"][tt]["units"][str(ub)]
            s = a.set if a.set in da["capture_response"] and a.set in db["capture_response"] else next(iter(da["capture_response"]))
            row = {"tetrode": tt, "unit_a": ua, "unit_b": ub, "template_cos": round(float(cos[i, j]), 3), "set": s,
                   "trough_uV_a": da["template_trough_uV"], "trough_uV_b": db["template_trough_uV"],
                   "n_a": da["n_capture"], "n_b": db["n_capture"]}
            for k in ("base_hz", "onset_1_10ms_hz", "in_pulse_hz", "post_0_50ms_hz"):
                row[f"{k}_a"] = round(da["capture_response"][s][k], 2)
                row[f"{k}_b"] = round(db["capture_response"][s][k], 2)
            for name, tag, u in (("a", tag_a, ua), ("b", tag_b, ub)):
                if tag is not None:
                    t = tag[(tag.unit == u) & (tag.set == s)]
                    row[f"class_{name}"] = t["class"].iloc[0] if len(t) else None
            rows.append(row)
        for i, u in enumerate(ia):
            if i not in used_a:
                rows.append({"tetrode": tt, "unit_a": u, "unit_b": None, "template_cos": round(float(cos[i].max()), 3) if cos.size else None})
        for j, u in enumerate(ib):
            if j not in used_b:
                rows.append({"tetrode": tt, "unit_a": None, "unit_b": u, "template_cos": round(float(cos[:, j].max()), 3) if cos.size else None})
    df = pd.DataFrame(rows)
    path = RESULTS / f"capture_xsession_{a.experiments[0]}_{a.experiments[1]}_{a.out}.csv"
    df.to_csv(path, index=False)
    pd.set_option("display.width", 250)
    print(df.to_string(index=False))
    matched = df.dropna(subset=["unit_a", "unit_b"])
    print(f"\n{len(matched)} matched pairs at cosine >= {a.min_cos}; wrote {path}")


if __name__ == "__main__":
    main()
