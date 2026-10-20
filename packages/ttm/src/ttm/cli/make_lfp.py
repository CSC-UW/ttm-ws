"""``ttm-make-lfp``: 30 kHz store -> 625 Hz LFP store, optionally -> 125 Hz tetrode-lead store.

``make_lfp`` and ``make_subsampled_lfp`` in :mod:`ttm.lfp` do the
work; ``--check-plot`` writes anti-aliasing verification figures to
``<experiment>/figures/``.
"""

from __future__ import annotations

import argparse

import numpy as np

from ttm.cli._common import (
    add_experiment_args,
    experiment_from_args,
    line_buffered,
)
from ttm.lfp import make_lfp, make_subsampled_lfp, open_lfps_dataarray

CHECK_T0_S = 3600.0


def _window(da, channel, t0, dur):
    t = da["time"].values
    i0 = int(np.searchsorted(t, t[0] + t0))
    i1 = int(np.searchsorted(t, t[0] + t0 + dur))
    return da.sel(channel=channel).values[i0:i1]


def main(argv: list[str] | None = None) -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--source-compressor", default="blosc-zstd", help="Which 30 kHz store to read")
    p.add_argument("--subsample", type=int, default=125, help="Tetrode-lead sub-LFP rate (Hz); 0 to skip")
    p.add_argument("--n-jobs", type=int, default=16)
    p.add_argument("--check-plot", action="store_true")
    a = p.parse_args(argv)
    exp = experiment_from_args(a)

    src = exp.store_path(a.source_compressor)
    if exp.lfp_path.exists():
        print(f"exists, skipping: {exp.lfp_path}")
    else:
        make_lfp(src, exp.lfp_path, n_jobs=a.n_jobs)
    lfp = open_lfps_dataarray(exp.lfp_path)
    print(f"LFP {lfp.dims} {lfp.shape} @ {lfp.attrs['fs']} Hz")

    sub = None
    if a.subsample:
        sub_path = exp.subsampled_lfp_path(a.subsample)
        if sub_path.exists():
            print(f"exists, skipping: {sub_path}")
        else:
            make_subsampled_lfp(exp.lfp_path, sub_path, resample_rate=a.subsample, n_jobs=a.n_jobs)
        sub = open_lfps_dataarray(sub_path)
        print(f"sub-LFP {sub.dims} {sub.shape} @ {sub.attrs['fs']} Hz, tetrodes "
              f"{list(map(str, sub['tetrode'].values)) if 'tetrode' in sub.coords else 'n/a'}")

    if a.check_plot:
        import spikeinterface as si

        from ttm.plotting import plot_decimation_check

        exp.figures_dir.mkdir(parents=True, exist_ok=True)
        rec = si.read_zarr(str(src))
        fs = rec.sampling_frequency
        ch = rec.get_channel_ids()[0]
        raw = rec.get_traces(start_frame=int(CHECK_T0_S * fs), end_frame=int((CHECK_T0_S + 60) * fs),
                             channel_ids=[ch], return_in_uV=True)[:, 0]
        lead = str(lfp["channel"].values[0])
        fig = plot_decimation_check(raw, fs, _window(lfp, lead, CHECK_T0_S, 60), float(lfp.attrs["fs"]),
                                    title=f"{exp.session} raw -> LFP ({lead}, 60 s @ {CHECK_T0_S:.0f} s)")
        out = exp.figures_dir / "lfp_check.png"
        fig.savefig(out, dpi=110)
        print(f"saved {out}")
        if sub is not None:
            lead = str(sub["channel"].values[0])
            fig = plot_decimation_check(
                _window(lfp, lead, CHECK_T0_S, 120), float(lfp.attrs["fs"]),
                _window(sub, lead, CHECK_T0_S, 120), float(sub.attrs["fs"]),
                title=f"{exp.session} LFP -> sub-LFP ({lead}, 120 s @ {CHECK_T0_S:.0f} s)",
            )
            out = exp.figures_dir / "lfp_subsample_check.png"
            fig.savefig(out, dpi=110)
            print(f"saved {out}")
    print("DONE")


if __name__ == "__main__":
    main()
