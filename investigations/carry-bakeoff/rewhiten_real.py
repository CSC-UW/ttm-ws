"""Per-window whitening against the first 300 s' whitening on the real 48 h, windows 48 to 59.

Both arms are the per-tetrode whitened carry with blend 0.2 (``R48_W_s45_c50_lnb_b2``); the
re-whitened one (``R48_W_s45_c50_lnb_b2_rw``) starts from the other's bank after window 47. Per
unit and window: accepted spikes in each arm, the share of each arm's spikes the other also
accepted (within 10 samples), refractory conflicts per accepted spike (reason 3), the medians of
snr, peeled cosine and residual of the accepted spikes, the medians of snr and peeled cosine of
the spikes only one arm accepted, and the cosine between the two arms' templates after the
window, in plain space.

Writes ``results/rewhiten_real.csv``; prints per-tetrode sums and medians.

Run: .venv/bin/python investigations/carry-bakeoff/rewhiten_real.py
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402
from score import near  # noqa: E402

HERE = pathlib.Path(__file__).parent
ROOT = bk.BAKE_ROOT / "real48"
FIXED, RW = ROOT / "R48_W_s45_c50_lnb_b2", ROOT / "R48_W_s45_c50_lnb_b2_rw"
GROUPS = np.repeat(np.arange(16), 4)


def plain(bank: bk.Bank, W: np.ndarray) -> bk.Bank:
    return bk.whiten_bank(bank, np.linalg.inv(W.astype(np.float64)), GROUPS)


def med(x: np.ndarray) -> float:
    return float(np.nanmedian(x)) if x.size else np.nan


def main():
    W_fixed = bk.whitening_matrix(False)
    rows = []
    for w in range(48, 60):
        if not (RW / f"spikes_w{w:02d}.npz").exists():
            break
        sp = {k: np.load(d / f"spikes_w{w:02d}.npz") for k, d in (("fixed", FIXED), ("rw", RW))}
        bf = plain(bk.Bank.load(FIXED / f"bank_w{w:02d}.npz"), W_fixed)
        br = plain(bk.Bank.load(RW / f"bank_w{w:02d}.npz"), np.load(RW / f"whiten_w{w:02d}.npy"))
        for k_f, u in enumerate(bf.ids):
            k_r = np.flatnonzero((br.ids == u) & (br.group == bf.group[k_f]))
            if not k_r.size:
                continue
            row = {"window": w, "unit": int(u), "tetrode": f"TT{bf.group[k_f] + 1}"}
            acc, fr = {}, {}
            for arm, s in sp.items():
                m = s["unit"] == u
                acc[arm] = m & (s["reason"] == 0)
                fr[arm] = s["frame"][acc[arm]]
                row[f"n_{arm}"] = int(acc[arm].sum())
                row[f"refr_{arm}"] = float((m & (s["reason"] == 3)).sum() / max(acc[arm].sum(), 1))
                for c in ("snr", "cos_pl", "resid"):
                    row[f"{c}_{arm}"] = med(s[c][acc[arm]])
            for arm, other in (("fixed", "rw"), ("rw", "fixed")):
                both = near(fr[arm], np.sort(fr[other]))
                row[f"shared_{arm}"] = float(both.mean()) if both.size else np.nan
                s = sp[arm]
                only = np.flatnonzero(acc[arm])[~both]
                row[f"only_snr_{arm}"] = med(s["snr"][only])
                row[f"only_cos_pl_{arm}"] = med(s["cos_pl"][only])
            tf, tr = bf.T[k_f].astype(np.float64), br.T[k_r[0]].astype(np.float64)
            row["tmpl_cos"] = float((tf * tr).sum() / np.sqrt((tf * tf).sum() * (tr * tr).sum()))
            row["tmpl_scale"] = float((tr * tf).sum() / (tf * tf).sum())
            rows.append(row)
    df = pd.DataFrame(rows)
    out = HERE / "results" / "rewhiten_real.csv"
    df.to_csv(out, index=False)
    by = df.groupby("tetrode", sort=False).agg(
        n_fixed=("n_fixed", "sum"), n_rw=("n_rw", "sum"), shared_fixed=("shared_fixed", "median"),
        shared_rw=("shared_rw", "median"), refr_fixed=("refr_fixed", "median"), refr_rw=("refr_rw", "median"),
        cos_pl_fixed=("cos_pl_fixed", "median"), cos_pl_rw=("cos_pl_rw", "median"),
        only_snr_fixed=("only_snr_fixed", "median"), only_snr_rw=("only_snr_rw", "median"),
        tmpl_cos_min=("tmpl_cos", "min"), tmpl_scale=("tmpl_scale", "median"))
    pd.set_option("display.width", 200)
    print(f"windows {df.window.min()}-{df.window.max()}, {df.unit.nunique()} units")
    print(by.round(3).to_string())
    print(f"written {out}")


if __name__ == "__main__":
    main()
