"""TTM-003-48H carry, dry run of `ttm.carry.readmit_window`: replayed with the parameters the carry ran with
(it must reproduce the saved reasons, except where the pulse rule took fits inside bad periods), and with
the 5 ms span. Nothing is written to the run. After `ttm-sort readmit`, it reads the carry's decisions from
`reason_carry` and its parameters from `carried_with`. Writes readmit_check.csv next to this script.

Run from the workspace root: .venv/bin/python investigations/pulse-rule-span/readmit_check.py
"""
import dataclasses, json, pathlib
import numpy as np, pandas as pd
from joblib import Parallel, delayed
from ttm import carry as cr

RUN = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main")
info = json.loads((RUN / "carry_params.json").read_text())
FS = 30000.0
rec = cr.open_sorting_recording(RUN / "binary", ["CH22"])
N = rec.get_num_frames()
pulses = cr.read_pulses(info["inputs"]["pulses"], FS)
bad = cr.read_bad_periods(info["inputs"]["bad_periods"], FS)
P_old = cr.CarryParams(**cr._params_dict(info.get("carried_with", info["params"])))
assert P_old.pulse_relax_post_ms == 100.0 and P_old.pulse_exclude_post_ms == 100.0
P_new = dataclasses.replace(P_old, pulse_relax_post_ms=5.0)
W = int(P_old.window_s * FS)


def one(w):
    z = np.load(RUN / "carry" / f"spikes_w{w:02d}.npz")
    zz = {k: z[k] for k in ("frame", "unit", "a", "snr", "cos_pl")}
    zz["reason"] = z["reason_carry"] if "reason_carry" in z.files else z["reason"]
    bankin = cr.Bank.load(RUN / "carry" / f"bankin_w{w:02d}.npz")
    a = w * W
    saved = zz["reason"]
    r_old, d_old = cr.readmit_window(zz, bankin, P_old, fs=FS, a=a, n_frames=min(W, N - a), pulses=pulses, bad_periods=bad)
    r_new, d_new = cr.readmit_window(zz, bankin, P_new, fs=FS, a=a, n_frames=min(W, N - a), pulses=pulses, bad_periods=bad)
    grp = dict(zip(bankin.ids.tolist(), bankin.group.tolist()))
    g = np.array([grp[int(u)] for u in zz["unit"]])
    in_bad = np.zeros(saved.size, bool)
    for gg, iv in bad.items():
        m = g == gg
        if m.any() and len(iv):
            st, en = np.array(iv).T
            in_bad[m] = cr.in_intervals(zz["frame"][m], st, en)
    saved_driven = np.load(RUN / "carry" / f"bank_w{w:02d}.npz")["driven"].tolist()
    diff = saved != r_old
    return {"window": w, "fits": int(saved.size), "accepted": int((saved == 0).sum()),
            "replay_differs": int(diff.sum()), "differs_in_bad": int((diff & in_bad).sum()),
            "lost_in_bad": int(((saved == 0) & (r_old != 0) & in_bad).sum()),
            "differs_outside_bad": int((diff & ~in_bad).sum()),
            "sustained_same": d_old == sorted(saved_driven),
            "new_lost": int(((r_old == 0) & (r_new != 0)).sum()), "new_gained": int(((r_old != 0) & (r_new == 0)).sum())}


rows = Parallel(n_jobs=12)(delayed(one)(w) for w in range(len(list((RUN / "carry").glob("spikes_w*.npz")))))
df = pd.DataFrame(rows)
df.to_csv(pathlib.Path(__file__).with_suffix(".csv"), index=False)
print(df[df.replay_differs + df.new_lost + df.new_gained > 0].to_string(index=False))
print(df[["fits", "accepted", "replay_differs", "differs_in_bad", "lost_in_bad", "differs_outside_bad", "new_lost", "new_gained"]].sum())
print("sustained sets reproduced in every window:", bool(df.sustained_same.all()))
