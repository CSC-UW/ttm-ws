"""Channel quality over a long session: dead, saturated or noisy wires and the 400 Hz comb.

Reads 60 s windows at evenly spaced times of the store and, per channel and
window, measures the raw standard deviation, the fraction of samples at the
int16 rails, the 300 to 6000 Hz noise MAD, the 5 MAD crossing rate, and the
power at 400 Hz relative to its neighbours (the comb of
``investigations/acquisition-noise/``). Writes ``results/channel_qc_<E>.csv``
(one row per channel and window), ``results/channel_qc_<E>_summary.csv``
(medians per channel with flags) and ``figures/channel_qc.png``.

Flags: ``dead`` raw std under 20 uV, ``noisy`` raw std over 600 uV,
``saturated`` over 0.1 % of samples at a rail, ``comb`` 400 Hz power more
than 10 times the 380 to 420 Hz flanks.

Usage: python subjects/ttm-003/scripts/channel_qc.py --subject TTM-003 --experiment TTM-003-48H [--n-windows 12]
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal as sps
import spikeinterface as si
from joblib import Parallel, delayed

from ttm.cli._common import add_experiment_args, experiment_from_args

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000
WIN_S = 60
RAIL = 32767


def one_window(store: str, start: int, ch_names: list[str], groups: np.ndarray, gain: float) -> pd.DataFrame:
    rec = si.read_zarr(store)
    raw = rec.get_traces(start_frame=start, end_frame=start + WIN_S * FS, return_in_uV=False)
    x = raw.astype(np.float32) * gain
    sos = sps.butter(3, [300, 6000], btype="bandpass", fs=FS, output="sos")
    xf = sps.sosfiltfilt(sos, x, axis=0)
    mad = np.median(np.abs(xf), axis=0) / 0.6745
    f, pxx = sps.welch(x, fs=FS, nperseg=FS, axis=0)
    comb = pxx[(f >= 398) & (f <= 402)].max(0)
    flank = np.median(pxx[((f >= 380) & (f < 395)) | ((f > 405) & (f <= 420))], axis=0)
    rows = []
    for c, name in enumerate(ch_names):
        below = xf[:, c] < -5 * mad[c]
        rows.append({
            "start_h": start / FS / 3600, "channel": name, "tetrode": f"TT{groups[c] + 1}",
            "raw_std_uV": float(x[:, c].std()),
            "rail_frac": float((np.abs(raw[:, c]) >= RAIL).mean()),
            "mad_uV": float(mad[c]),
            "crossings_5mad_hz": float((below[1:] & ~below[:-1]).sum() / WIN_S),
            "comb_400hz_ratio": float(comb[c] / flank[c]) if flank[c] > 0 else float("nan"),
            "lfp_std_uV": float(sps.sosfiltfilt(sps.butter(4, 300, btype="lowpass", fs=FS, output="sos"), x[:, c]).std()),
        })
    return pd.DataFrame(rows)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--n-windows", type=int, default=12)
    p.add_argument("--workers", type=int, default=12)
    a = p.parse_args()
    exp = experiment_from_args(a)
    rec = si.read_zarr(str(exp.store_path()))
    nfr = rec.get_num_frames()
    groups = np.asarray(rec.get_property("group")).astype(int)
    ch_names = [str(c) for c in rec.get_channel_ids()]
    gain = float(rec.get_channel_gains()[0])
    starts = np.linspace(0, nfr - WIN_S * FS, a.n_windows + 2)[1:-1].astype(int)
    df = pd.concat(Parallel(n_jobs=a.workers)(delayed(one_window)(str(exp.store_path()), int(s), ch_names, groups, gain) for s in starts))
    df.round(5).to_csv(RESULTS / f"channel_qc_{exp.experiment}.csv", index=False)
    med = df.groupby(["channel", "tetrode"], sort=False).median(numeric_only=True).drop(columns="start_h")
    med["rail_frac_max"] = df.groupby(["channel", "tetrode"], sort=False).rail_frac.max()
    med["dead"] = med.raw_std_uV < 20
    med["noisy"] = med.raw_std_uV > 600
    med["saturated"] = med.rail_frac_max > 1e-3
    med["comb"] = med.comb_400hz_ratio > 10
    med["flag"] = med[["dead", "noisy", "saturated", "comb"]].apply(lambda r: ",".join(k for k, v in r.items() if v), axis=1)
    med = med.reset_index()
    med["in_manifest_bad"] = med.channel.isin(exp.badChannels) | med.tetrode.isin(exp.badTetrodes)
    med.round(4).to_csv(RESULTS / f"channel_qc_{exp.experiment}_summary.csv", index=False)

    fig, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True)
    piv = lambda col: df.pivot(index="channel", columns="start_h", values=col).loc[ch_names]
    for ax, col, log in ((axes[0], "raw_std_uV", True), (axes[1], "mad_uV", True), (axes[2], "crossings_5mad_hz", False), (axes[3], "comb_400hz_ratio", True)):
        m = piv(col)
        for k, name in enumerate(ch_names):
            ax.plot(m.columns, m.loc[name].values, lw=0.7, color=plt.cm.tab20(groups[k] % 20), alpha=0.8)
        if log:
            ax.set_yscale("log")
        ax.set_ylabel(col)
    axes[-1].set_xlabel("session hour")
    fig.suptitle(f"{exp.experiment}: per-channel quality in {a.n_windows} windows of {WIN_S} s (colour = tetrode)")
    fig.tight_layout(); fig.savefig(exp.figures_dir / "channel_qc.png", dpi=110); plt.close(fig)
    flagged = med[med.flag != ""]
    print(flagged[["channel", "tetrode", "raw_std_uV", "mad_uV", "crossings_5mad_hz", "rail_frac_max", "comb_400hz_ratio", "flag", "in_manifest_bad"]].round(2).to_string(index=False))
    print(f"manifest bad channels: {exp.badChannels}, bad tetrodes: {exp.badTetrodes}")
    print(f"wrote {RESULTS / f'channel_qc_{exp.experiment}_summary.csv'}")


if __name__ == "__main__":
    main()
