"""``ttm-swa``: slow-wave activity (delta power) products and figures for one experiment.

Steps (``ttm-swa STEP``):

- ``inst``: 0.5-4 Hz bandpass + Hilbert envelope on the 625 Hz LFP, decimated
  to 125 Hz, averaged per tetrode -> ``delta.idelta.zarr``.
- ``stft``: 4 s DPSS STFT band powers (delta, vlad, theta, sigma, gamma) per
  tetrode -> ``stft_bandpowers.zarr``.
- ``plot``: mean and per-tetrode SWA traces with the light/dark cycle and
  deprivation window -> ``figures/swa_{inst,stft}_{mean,per_tetrode}.png``.
- ``all``: the three in order.

Outlier transients are replaced at plot time only; the stores stay raw.
"""

from __future__ import annotations

import argparse

from ttm import experiment as ep
from ttm import power
from ttm.cli._common import (
    add_experiment_args,
    experiment_from_args,
    line_buffered,
)
from ttm.lfp import open_lfps_dataarray
from ttm.subject import Experiment

STFT_SMOOTHING, STFT_Q = 14, 10  # findlay2025a fig_2a
YMAX_PCT = 99.5


def step_inst(exp: Experiment) -> None:
    from dask.diagnostics import ProgressBar

    low, high = power.BANDS["delta"]
    lfp = open_lfps_dataarray(exp.lfp_path)
    ipwr = power.extract_instantaneous_bandpower(lfp, low, high)
    print(f"instantaneous delta -> {ipwr.dims} {ipwr.shape} @ {ipwr.fs} Hz")
    with ProgressBar():
        power.save_instantaneous_power(ipwr, exp.idelta_path)
    print(f"Wrote {exp.idelta_path}")


def step_stft(exp: Experiment) -> None:
    lfp = open_lfps_dataarray(exp.lfp_path)
    ds = power.compute_stft_bandpowers(lfp)
    print(f"STFT band powers: {list(ds.data_vars)} {ds['delta'].dims}{ds['delta'].shape} "
          f"@ {ds.attrs['fs']} Hz (n_fft={ds.attrs['stft_n_fft']})")
    power.save_stft_bandpowers(ds, exp.stft_bandpowers_path)
    print(f"Wrote {exp.stft_bandpowers_path}")


def step_plot(exp: Experiment) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from ttm import plotting

    params = ep.load_experiment_params(exp.experiment_params_path)
    inst = power.replace_outliers_per_tetrode(power.open_instantaneous_power(exp.idelta_path).compute())
    stft = power.replace_outliers_per_tetrode(power.open_stft_bandpowers(exp.stft_bandpowers_path)["delta"].compute())

    def dt(da):
        return float(np.median(np.diff(da["time"].values)))

    scale = dt(stft) / dt(inst)
    inst_smoothing, inst_q = round(STFT_SMOOTHING * scale), round(STFT_Q * scale)
    exp.figures_dir.mkdir(parents=True, exist_ok=True)
    products = [
        ("inst", "Instantaneous delta power", inst, inst_smoothing, inst_q),
        ("stft", "STFT delta band power", stft, STFT_SMOOTHING, STFT_Q),
    ]
    for key, label, da, smoothing, q in products:
        title = f"{exp.subject} {exp.experiment} - SWA ({label})"
        fig, ax = plt.subplots(figsize=(14, 3))
        plotting.plot_swa_overview(da.mean("tetrode"), params, ax, smoothing=smoothing, q=q,
                                   title=title, ymax_percentile=YMAX_PCT)
        out = exp.figures_dir / f"swa_{key}_mean.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out}")
        fig = plotting.plot_swa_small_multiples(da, params, smoothing=smoothing, q=q,
                                                title=title, ymax_percentile=YMAX_PCT)
        out = exp.figures_dir / f"swa_{key}_per_tetrode.png"
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out}")


STEPS = {"inst": step_inst, "stft": step_stft, "plot": step_plot}


def main(argv: list[str] | None = None) -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("step", choices=[*STEPS, "all"])
    add_experiment_args(p)
    a = p.parse_args(argv)
    exp = experiment_from_args(a)
    for name in (list(STEPS) if a.step == "all" else [a.step]):
        print(f"== {name}")
        STEPS[name](exp)


if __name__ == "__main__":
    main()
