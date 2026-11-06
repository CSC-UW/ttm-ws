"""TTM-003-48H after `ttm-sort readmit`: where the fits the replay dropped could have changed a template
update (step 5.5), had the carry run with the 5 ms span. The update needs 100 accepted spikes in the
window, and draws up to 500 of the unit's accepted lead spikes outside the excluded span and bad periods
(lead: no accepted spike of the unit in the 10 ms before). Per unit and window with a dropped fit: whether
the count crosses 100, and how many spikes join or leave the set the update draws from.

Writes template_effect.csv next to this script.

Run from the workspace root: .venv/bin/python investigations/pulse-rule-span/template_effect.py
"""
import json
import pathlib

import numpy as np
import pandas as pd
from ttm import carry as cr

HERE = pathlib.Path(__file__).parent
RUN = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main")
info = json.loads((RUN / "carry_params.json").read_text())
FS = 30000.0
P = cr.CarryParams(**info["params"])
on, off = cr.read_pulses(info["inputs"]["pulses"], FS)
bad = cr.read_bad_periods(info["inputs"]["bad_periods"], FS)
post, lead_f = int(P.pulse_exclude_post_ms * 1e-3 * FS), int(P.lead_ms * 1e-3 * FS)
W = int(P.window_s * FS)


def eligible(f, g, a, n):
    f = np.sort(f)
    keep = (np.r_[np.inf, np.diff(f)] > lead_f) & ~cr.in_intervals(f, on, off + post)
    iv = bad.get(int(g), [])
    if len(iv):
        st, en = np.array(iv).T
        keep &= ~cr.in_intervals(f, st, en)
    keep &= (f - a >= P.nbefore) & (f - a + P.nafter <= n)
    return set(f[keep].tolist())


rows = []
for fpath in sorted((RUN / "carry").glob("spikes_w*.npz")):
    w = int(fpath.stem[-2:])
    z = np.load(fpath)
    old, new = z["reason_carry"], z["reason"]
    changed = np.flatnonzero((old == 0) != (new == 0))
    if not changed.size:
        continue
    bk = cr.Bank.load(RUN / "carry" / f"bankin_w{w:02d}.npz")
    grp = dict(zip(bk.ids.tolist(), bk.group.tolist()))
    for u in np.unique(z["unit"][changed]):
        m = z["unit"] == u
        fo, fn = z["frame"][m & (old == 0)], z["frame"][m & (new == 0)]
        eo, en_ = eligible(fo, grp[int(u)], w * W, W), eligible(fn, grp[int(u)], w * W, W)
        rows.append({"window": w, "unit": int(u), "dropped": int(((old == 0) & (new != 0) & m).sum()),
                     "accepted_old": int(fo.size), "accepted_new": int(fn.size),
                     "crosses_100": (fo.size >= P.rebuild_min_accepted) != (fn.size >= P.rebuild_min_accepted),
                     "eligible_old": len(eo), "eligible_new": len(en_), "joining": len(en_ - eo), "leaving": len(eo - en_)})
df = pd.DataFrame(rows)
df.to_csv(HERE / "template_effect.csv", index=False)
print(df.to_string(index=False))
print("unit-windows:", len(df), "crossing 100:", int(df.crosses_100.sum()), "with a changed draw set:",
      int(((df.joining + df.leaving) > 0).sum()), "spikes joining:", int(df.joining.sum()), "leaving:", int(df.leaving.sum()))
