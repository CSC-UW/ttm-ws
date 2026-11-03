"""Tables from ceiling.py's per-fit and per-unit outputs.

Run: .venv/bin/python investigations/template-superres/summarize.py w03_w08
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

OUT = ("/tmp/claude-295375638/-Users-gfindlay-ad-wisc-edu-projects-ttm-ws/"
       "e980bd61-79fb-4043-a37f-61959b628a69/scratchpad/superres/")
SNR, COS = 4.5, 0.5


def q(s: pd.Series) -> str:
    return f"{s.median():.4f} [{s.quantile(0.9):.4f}, {s.quantile(0.99):.4f}]"


def main(tag: str) -> None:
    d = pd.read_parquet(OUT + f"fits_{tag}.parquet")
    u = pd.read_csv(OUT + f"units_{tag}.csv")
    d["size"] = np.where(d.trough_mad <= 6, "small", "large")
    d["kind"] = np.where(d.unit >= 1000, "gt", "real")
    d["label"] = np.where(d.kind == "real", "-", np.where(d.hit.fillna(False).astype(bool), "hit", "fp"))
    d["dsnr"] = d.snr_best - d.snr0
    d["dcos"] = d.cos_best - d.cos0
    d["de_rel"] = (d.e_ls0 - d.e_best) / d.e_ls0
    d["de_st_rel"] = (d.e_stored - d.e_best) / d.e_stored
    d["dpeak"] = d.peak_stored - d.peak_best
    acc = d[d.reason == 0]
    print(f"fits scored {len(d):,}; accepted {len(acc):,}; units {d.unit.nunique()} ({d[d.kind=='gt'].unit.nunique()} gt)")
    print("\nAccepted fits: median [p90, p99] of the gain at the best delay")
    rows = []
    for key, x in acc.groupby(["kind", "size", "label"]):
        rows.append({"group": "/".join(key), "n": len(x), "units": x.unit.nunique(),
                     "dsnr": q(x.dsnr), "dsnr_rel": q(x.dsnr / x.snr0), "dcos": q(x.dcos),
                     "resid_drop_vs_ls0": q(x.de_rel), "resid_drop_vs_stored": q(x.de_st_rel),
                     "|d|": f"{x.d_best.abs().median():.3f}"})
    print(pd.DataFrame(rows).to_string(index=False))

    print("\nThreshold crossings (fits scored, reason 0-2): share of scored fits that pass both snr>=4.5 and"
          " cos_pl>=0.5 at the best delay but not at the integer one, against the accepted count")
    rows = []
    for key, x in d.groupby(["kind", "size", "label"]):
        pass0 = (x.snr0 >= SNR) & (x.cos0 >= COS)
        passb = (x.snr_best >= SNR) & (x.cos_best >= COS)
        new = passb & ~pass0
        rows.append({"group": "/".join(key), "n_scored": len(x), "accepted": int((x.reason == 0).sum()),
                     "pass_int": int(pass0.sum()), "new_pass": int(new.sum()),
                     "new/pass_int": f"{new.sum() / max(pass0.sum(), 1):.4f}",
                     "new_by_snr": int((new & (x.snr0 < SNR)).sum()),
                     "new_by_cos": int((new & (x.cos0 < COS)).sum())})
    print(pd.DataFrame(rows).to_string(index=False))

    print("\nGround-truth hits: delay estimate against the injected sub-sample offset")
    h = d[(d.label == "hit") & (d.reason == 0)].copy()
    rows = []
    for (unit, w), x in h.groupby(["unit", "window"]):
        off = x.gt_off - x.gt_off.median()
        rows.append({"unit": unit, "window": w, "trough_mad": round(x.trough_mad.iloc[0], 1), "n": len(x),
                     "sd_true": off.std(), "corr": np.corrcoef(off, x.d_best)[0, 1],
                     "rmse_est": np.sqrt(((x.d_best - x.d_best.median()) - off).pow(2).mean())})
    t = pd.DataFrame(rows)
    print(t.groupby(t.trough_mad <= 6)[["n", "sd_true", "corr", "rmse_est"]].median()
          .rename(index={True: "small", False: "large"}).to_string())

    print("\nLarge-unit residual peaks (accepted, real and gt): stored fit vs best-delay LS fit, noise units")
    big = acc[acc.trough_mad >= 10]
    if len(big):
        for th in (4, 5):
            print(f"  share with residual peak >= {th}: stored {np.mean(big.peak_stored >= th):.4f}, "
                  f"best {np.mean(big.peak_best >= th):.4f} (n {len(big):,}, units {big.unit.nunique()})")
        print(f"  resid energy drop vs stored, trough>=10: {q(big.de_st_rel)}; vs ls0 {q(big.de_rel)}")

    print("\nTemplate blur: energy of the fractionally aligned mean over the integer mean (E1/E0 - 1);"
          " same procedure on template + noise at the integer time (control); oracle on gt")
    u["size"] = np.where(u.trough_mad <= 6, "small", "large")
    u["gain"] = u.E1 / u.E0 - 1
    u["gain_ctrl"] = u.E1_noise_ctrl / u.E0_noise_ctrl - 1
    u["gain_excess"] = u.gain - u.gain_ctrl
    u["trough_gain"] = u.min1 / u.min0 - 1
    if "E2_oracle" in u:
        u["gain_oracle"] = u.E2_oracle / u.E0 - 1
        u["dcos_oracle"] = u.cos2_true - u.cos0_true
        u["dcos_est"] = u.cos1_true - u.cos0_true
    cols = [c for c in ("gain", "gain_ctrl", "gain_excess", "trough_gain", "gain_oracle", "dcos_est",
                        "dcos_oracle") if c in u]
    u["kind"] = np.where(u.unit >= 1000, "gt", "real")
    print(u.groupby(["kind", "size"])[cols].median().to_string())
    print(u.groupby(["kind", "size"]).size().to_string())


if __name__ == "__main__":
    main(sys.argv[1])
