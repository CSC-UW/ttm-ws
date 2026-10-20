"""``ttm-emg``: synthetic EMG (both estimators) for one experiment, plotted against SWA.

Default: a 30 min validation window 10 h into the first segment, plotted only.
``--full`` processes every segment, writes ``synthetic_emg_methods.zarr`` and
adds data-driven 10 min zooms (lowest-EMG, highest-EMG, deprivation midpoint).
``--plot-only`` re-plots from the saved store. SWA comes from
``delta.idelta.zarr`` (``ttm-swa inst``).
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd
import xarray as xr

from ttm import emg as emg_lib
from ttm import experiment as ep
from ttm import power
from ttm.cli._common import (
    add_experiment_args,
    experiment_from_args,
    line_buffered,
)

LIGHT_COLORS = {"on": "gold", "off": "0.35"}


def _overlay(ax, params):
    x0, x1 = ax.get_xlim()
    intervals, labels = ep.get_light_dark_periods(params)
    for (a, b), lab in zip(intervals, labels):
        a, b = max(a / 3600, x0), min(b / 3600, x1)
        if b > a:
            ax.axvspan(a, b, ymin=0.0, ymax=0.04, color=LIGHT_COLORS[lab], zorder=1000, ec="none", clip_on=False)
    dep = ep.get_deprivation_period(params)
    if dep is not None:
        ax.axvspan(dep[0] / 3600, dep[1] / 3600, color="red", alpha=0.10, zorder=0)
    ax.set_xlim(x0, x1)


def _trace(ax, da_, smoothing, color, label, ymax_pct=None, xlim_h=None):
    d = da_.rolling(time=smoothing, center=True, min_periods=1).mean() if smoothing > 1 else da_
    x = d["time"].values / 3600.0
    y = d.values
    ax.plot(x, y, lw=0.6, color=color)
    ax.set_ylabel(label)
    ax.margins(x=0)
    vis = (x >= xlim_h[0]) & (x <= xlim_h[1]) if xlim_h is not None else np.isfinite(y)
    yv = y[vis & np.isfinite(y)]
    if ymax_pct is not None and yv.size:
        top = np.nanpercentile(yv, ymax_pct)
        if np.isfinite(top) and top > 0:
            ax.set_ylim(min(0.0, float(yv.min())), top * 1.08)


def plot_panels(ds, swa, params, out_png, *, xlim_h=None, emg_smooth=1, swa_smooth=1, mode="full", title=""):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g_pct = 99.8 if mode == "full" else 99.0
    swa_pct = 99.5 if mode == "full" else 98.0
    fig, axs = plt.subplots(3, 1, figsize=(14, 8), sharex=True)
    _trace(axs[0], ds["per_window"], emg_smooth, "C0", "EMG\nper_window", xlim_h=xlim_h)
    _trace(axs[1], ds["global"], emg_smooth, "C1", "EMG\nglobal", ymax_pct=g_pct, xlim_h=xlim_h)
    _trace(axs[2], swa, swa_smooth, "C2", "SWA\n(delta uV^2)", ymax_pct=swa_pct, xlim_h=xlim_h)
    axs[0].axhline(1.0, ls=":", c="gray", lw=0.7)
    if xlim_h is not None:
        for ax in axs:
            ax.set_xlim(*xlim_h)
    for ax in axs:
        _overlay(ax, params)
    axs[-1].set_xlabel("session time (h)   [gold=lights-on, dark=lights-off, red=deprivation]")
    axs[0].set_title(title)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_png}")


def _extreme_time(da_, fs, kind, exclude=None, smooth_s=30):
    sm = da_.rolling(time=int(fs * smooth_s), center=True, min_periods=1).mean()
    t = sm["time"].values
    y = sm.values.astype(float).copy()
    if exclude is not None:
        y[(t >= exclude[0]) & (t <= exclude[1])] = np.inf if kind == "min" else -np.inf
    return float(t[np.nanargmin(y) if kind == "min" else np.nanargmax(y)])


def zoom_windows(ds, params):
    dep = ep.get_deprivation_period(params)
    pw, fs = ds["per_window"], emg_lib.EMG_TARGET_SF
    wins = [("NREM_lowEMG", _extreme_time(pw, fs, "min", exclude=dep)),
            ("wake_highEMG", _extreme_time(pw, fs, "max", exclude=dep))]
    if dep is not None:
        wins.append(("deprivation", (dep[0] + dep[1]) / 2))
    return wins


def main(argv: list[str] | None = None) -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--full", action="store_true", help="process every segment and save the store")
    p.add_argument("--plot-only", action="store_true", help="skip generation; re-plot the saved store")
    p.add_argument("--source-compressor", default="blosc-zstd")
    p.add_argument("--tetrodes", nargs="+", default=list(emg_lib.DEFAULT_TETRODES))
    a = p.parse_args(argv)
    exp = experiment_from_args(a)
    params = ep.load_experiment_params(exp.experiment_params_path)
    full = a.full or a.plot_only

    if a.plot_only:
        ds = xr.open_zarr(exp.emg_path)
        print(f"loaded {exp.emg_path} ({ds.sizes['time']} samples)")
    else:
        slice_table = pd.read_csv(exp.slice_table_path)
        ds = emg_lib.synthetic_emg_from_store(
            exp.store_path(a.source_compressor), slice_table, tetrodes=tuple(a.tetrodes),
            window=None if a.full else (10 * 3600, 1800),
        )
        if a.full:
            ds.to_zarr(exp.emg_path, mode="w")
            print(f"saved {exp.emg_path}")

    swa = power.replace_outliers_per_tetrode(power.open_instantaneous_power(exp.idelta_path).compute()).mean("tetrode")
    tag = "full" if full else "window"
    plot_panels(ds, swa, params, exp.figures_dir / f"emg_methods_vs_swa_{tag}.png",
                emg_smooth=(101 if full else 11), swa_smooth=(1250 if full else 250), mode=tag,
                title=f"{exp.subject} {exp.experiment}: EMG (per_window vs global) & SWA "
                      f"({'full session' if full else '30 min window'})")
    if full:
        for label, center_s in zoom_windows(ds, params):
            c_h, half_h = center_s / 3600.0, 5.0 / 60.0
            plot_panels(ds, swa, params, exp.figures_dir / f"emg_methods_vs_swa_zoom_{label}.png",
                        xlim_h=(c_h - half_h, c_h + half_h), emg_smooth=1, swa_smooth=1250, mode="zoom",
                        title=f"{exp.subject} {exp.experiment}: 10-min zoom @ {label} (session {c_h:.2f} h)")
    print("DONE")


if __name__ == "__main__":
    main()
