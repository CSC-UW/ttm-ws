"""Stage A, laser pulses: rules for spikes inside pulses, applied to the permissive per-tetrode fits.

Inside each delivered pulse (onset to offset plus 100 ms) a rule with lower floors replaces the rule used
elsewhere. Judged on window 0 (train 1, 58 pulses):

- driven ground-truth units (two at 9 and 10 MAD, one at 5 MAD): recall of their half-size in-pulse
  spikes and of their first full-size spikes, and the share of their in-pulse spikes that are theirs;
- the other ground-truth units, which are silent from 5 ms after onset to 60 ms after offset: spikes per
  second they collect in that span (all false);
- real units: in-pulse rate against their rate outside pulses (most real units are suppressed by the
  light; a suppressed unit that gains in-pulse spikes has absorbed a neighbour's).

Run: .venv/bin/python investigations/carry-bakeoff/pulse_rules.py
"""
from __future__ import annotations

import json
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402
from admission import load, refractory  # noqa: E402
from score import near  # noqa: E402

OUT_RULE = (4.0, 0.5, 0.6)
PULSE_RULES = [None, (3.5, 0.3, 0.6), (3.0, 0.3, 0.6), (3.0, 0.3, 0.5), (2.5, 0.3, 0.6), (2.5, 0.3, 0.5),
               (3.0, 0.4, 0.6)]


def light_driven(fr, un, on, off, out_s, *, start_ms=15.0, min_ratio=2.0, p_max=1e-4, min_spikes=20):
    """Units that keep firing through the pulse: spikes from ``start_ms`` after onset to offset at least
    ``min_ratio`` times their rate outside pulses (Poisson, p < 1e-4). The onset volley is left out."""
    from scipy.stats import poisson

    out = set()
    st = on + int(start_ms * 1e-3 * bk.FS)
    span = (off - st).sum() / bk.FS
    for u in np.unique(un):
        f = fr[un == u]
        n_in = bk.in_intervals(f, st, off).sum()
        rate = (~bk.in_intervals(f, on, off + int(0.1 * bk.FS))).sum() / out_s
        if n_in >= min_spikes and n_in >= min_ratio * rate * span and poisson.sf(n_in - 1, rate * span) < p_max:
            out.add(int(u))
    return out


def main():
    w = 0
    gt = np.load(bk.BAKE / "gt.npz")
    plan = json.loads((bk.BAKE / "gt_plan.json").read_text())
    driven = [p["gt"] for p in plan if p["kind"] == "driven"]
    z = load("A_pt_open", w)
    on, off = bk.laser_pulses()
    m = (on >= w * bk.WIN) & (on < (w + 1) * bk.WIN)
    on, off = on[m], off[m]
    near_pulse = bk.in_intervals(z["frame"], on, off + int(0.1 * bk.FS))
    silent = bk.in_intervals(z["frame"], on + int(0.005 * bk.FS), off + int(0.06 * bk.FS))
    in_s = ((off + int(0.06 * bk.FS)) - (on + int(0.005 * bk.FS))).sum() / bk.FS
    out_s = bk.WIN / bk.FS - ((off + int(0.1 * bk.FS)) - on).sum() / bk.FS
    gm = (gt["frame"] >= w * bk.WIN) & (gt["frame"] < (w + 1) * bk.WIN)
    base = z["reason"] != 1
    k, r, c = OUT_RULE
    out_keep = base & ~near_pulse & (z["snr"] >= k) & (z["a"] >= r) & (z["cos_pl"] >= c)
    std_keep = refractory(z, out_keep | (base & near_pulse & (z["snr"] >= k) & (z["a"] >= r) & (z["cos_pl"] >= c)))
    resp = light_driven(z["frame"][std_keep], z["unit"][std_keep], on, off, out_s)
    print("units answering the light under the standard rule:", sorted(resp))
    rows = []
    for rule, only_driven in [(x, False) for x in PULSE_RULES] + [(x, True) for x in PULSE_RULES[1:]]:
        kp, rp, cp_ = rule if rule else OUT_RULE
        relax = base & near_pulse & (z["snr"] >= kp) & (z["a"] >= rp) & (z["cos_pl"] >= cp_)
        if only_driven:
            relax &= np.isin(z["unit"], list(resp))
        keep = refractory(z, std_keep | relax)
        row = {"pulse_rule": ("same as outside" if rule is None else f"snr>={kp} a>={rp} cos_pl>={cp_}")
               + (" (driven units only)" if only_driven else "")}
        fr, un = z["frame"][keep], z["unit"][keep]
        for i in driven:
            gi = gm & (gt["unit"] == i)
            mine = np.sort(fr[un == 1000 + i])
            for cls, name in ((4, "half"), (3, "first")):
                g = gt["frame"][gi & (gt["cls"] == cls)]
                row[f"GT{i}_{name}"] = near(g, mine).mean() if g.size else np.nan
            ip = mine[bk.in_intervals(mine, on, off + int(0.1 * bk.FS))]
            row[f"GT{i}_inpulse_prec"] = near(ip, np.sort(gt["frame"][gi])).mean() if ip.size else np.nan
        others = [p["gt"] for p in plan if p["kind"] != "driven"]
        sil = np.isin(un, [1000 + i for i in others]) & bk.in_intervals(fr, on + int(0.005 * bk.FS),
                                                                         off + int(0.06 * bk.FS))
        row["silent_GT_false_hz"] = sil.sum() / in_s / len(others)
        ratios = []
        for u in np.unique(un[un < 1000]):
            f = fr[un == u]
            rin = bk.in_intervals(f, on + int(0.01 * bk.FS), off).sum() / ((off - on - int(0.01 * bk.FS)).sum() / bk.FS)
            rout = (~bk.in_intervals(f, on, off + int(0.1 * bk.FS))).sum() / out_s
            if rout > 0.5:
                ratios.append(rin / rout)
        ratios = np.array(ratios)
        row["real_median_in_out"] = float(np.median(ratios))
        row["real_units_in_out_gt_0.5"] = int((ratios > 0.5).sum())
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(bk.BAKE / "pulse_rules_w00.csv", index=False)
    pd.set_option("display.width", 250)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))


if __name__ == "__main__":
    main()
