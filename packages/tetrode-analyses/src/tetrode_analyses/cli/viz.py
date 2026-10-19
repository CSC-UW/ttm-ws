"""The loupe visualization hand-off: build a bundle on the server, view it on a laptop.

``ttm-viz-bundle`` (server) assembles ``<session>/viz_bundle/``: the EMG store
copied verbatim, the 125 Hz tetrode-lead LFP re-exported as a plain xarray
zarr, the sorting flattened to ``spikes.parquet`` in session seconds, the loupe
``state_definitions.json`` and a ``manifest.json``. ``ttm-fetch-bundle viz``
copies it; ``ttm-loupe --data-dir DIR`` views it with no SpikeInterface, no
source recording and no server mounts (install extra ``viz``).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import pathlib
import re
import shutil
from importlib import resources

import numpy as np
import polars as pl

from tetrode_analyses.cli._common import (
    add_experiment_args,
    experiment_from_args,
    line_buffered,
)
from tetrode_analyses.viz import _tetrode_sort_key, tetrode_color_map

EMG_FILENAME = "synthetic_emg_methods.zarr"
LFP_FILENAME = "lfp.125hz.zarr"
SPIKES_FILENAME = "spikes.parquet"
STATE_DEFINITIONS_FILENAME = "state_definitions.json"
DEFAULT_SORTINGS_SUBDIR = "sortings_seed42_pcafix"
DEFAULT_SORTING = "blosc-43200s-train3600s"
DEFAULT_SOURCE_COMPRESSOR = "blosc-zstd"
EMG_METHODS = ("per_window", "global")
EMG_COLORS = {"per_window": "#1f77b4", "global": "#ff7f0e"}
SEPARATOR_PARAMS = {"gap": 0.6, "color": "#888888", "width": 1.0}


def build_spike_dataframe(aggregated, source_recording) -> pl.DataFrame:
    """Flatten an aggregated sorting to (time, unit_id, tetrode) in one pass.

    One ``sample_index_to_time`` call over the whole spike vector decompresses
    the recording's lazy time vector once, not once per unit.
    """
    sv = aggregated.to_spike_vector()
    frames = np.asarray(sv["sample_index"])
    unit_idx = np.asarray(sv["unit_index"])
    times = np.asarray(source_recording.sample_index_to_time(frames, segment_index=0))
    unit_ids = np.asarray(aggregated.unit_ids)
    groups = np.asarray(aggregated.get_property("group"))
    return pl.DataFrame({
        "time": times.astype(np.float64),
        "unit_id": unit_ids[unit_idx].astype(np.int32),
        "tt_num": (groups[unit_idx] + 1).astype(np.int16),
    }).with_columns(("TT" + pl.col("tt_num").cast(pl.String)).alias("tetrode")).drop("tt_num")


def tetrode_separators(spikes: pl.DataFrame) -> tuple[list[str], list[int] | None]:
    """Ordered tetrode labels + the last unit_id of each tetrode but the final one.

    Returns ``None`` boundaries when unit_ids are not tetrode-contiguous.
    """
    per_unit = spikes.select("unit_id", "tetrode").unique().sort("unit_id")
    tetrodes = sorted(per_unit["tetrode"].unique().to_list(), key=_tetrode_sort_key)
    spans = {tt: (int(per_unit.filter(pl.col("tetrode") == tt)["unit_id"].min()),
                  int(per_unit.filter(pl.col("tetrode") == tt)["unit_id"].max())) for tt in tetrodes}
    ordered = [spans[tt] for tt in tetrodes]
    contiguous = all(ordered[i][1] < ordered[i + 1][0] for i in range(len(ordered) - 1))
    if not contiguous:
        return tetrodes, None
    return tetrodes, [spans[tt][1] for tt in tetrodes[:-1]]


def _copytree(src: pathlib.Path, dst: pathlib.Path, *, overwrite: bool) -> None:
    if dst.exists():
        if not overwrite:
            print(f"  exists, skipping: {dst}")
            return
        shutil.rmtree(dst)
    print(f"  copy {src} -> {dst}")
    shutil.copytree(src, dst)


def main_build(argv: list[str] | None = None) -> None:
    line_buffered()
    import spikeinterface as si

    from tetrode_analyses.lfp import make_subsampled_lfp, open_lfps_dataarray

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--sorting", default=DEFAULT_SORTING)
    p.add_argument("--sortings-subdir", default=DEFAULT_SORTINGS_SUBDIR)
    p.add_argument("--source-compressor", default=DEFAULT_SOURCE_COMPRESSOR)
    p.add_argument("--bundle", type=pathlib.Path, default=None, help="Output dir (default <session>/viz_bundle)")
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args(argv)
    exp = experiment_from_args(a)

    sub_lfp = exp.subsampled_lfp_path(125)
    aggregated_dir = exp.sorting_dir(a.sorting, a.sortings_subdir) / "aggregated"
    bundle = a.bundle or exp.viz_bundle_dir
    bundle.mkdir(parents=True, exist_ok=True)
    print(f"Bundle -> {bundle}")

    print("[1/5] EMG")
    _copytree(exp.emg_path, bundle / EMG_FILENAME, overwrite=a.overwrite)

    print("[2/5] 125 Hz LFP")
    if not sub_lfp.exists():
        print(f"  building 125 Hz sub-LFP (missing): {sub_lfp}")
        make_subsampled_lfp(exp.lfp_path, sub_lfp, resample_rate=125, n_jobs=16)
    bundle_lfp = bundle / LFP_FILENAME
    da = open_lfps_dataarray(sub_lfp)
    if bundle_lfp.exists() and not a.overwrite:
        print(f"  exists, skipping: {bundle_lfp}")
    else:
        if bundle_lfp.exists():
            shutil.rmtree(bundle_lfp)
        print(f"  re-export {da.shape} @ {float(da.attrs['fs']):.1f} Hz -> {bundle_lfp}")
        da.to_dataset(name="lfp").to_zarr(bundle_lfp, mode="w")
    tetrodes_lfp = [str(t) for t in np.asarray(da["tetrode"].values)]

    print("[3/5] spikes")
    spikes_path = bundle / SPIKES_FILENAME
    if spikes_path.exists() and not a.overwrite:
        print(f"  exists, skipping: {spikes_path}")
        spikes = pl.read_parquet(spikes_path)
    else:
        agg = si.load(str(aggregated_dir))
        spikes = build_spike_dataframe(agg, si.read_zarr(str(exp.store_path(a.source_compressor))))
        spikes.write_parquet(spikes_path, compression="zstd")
    tetrode_order, boundaries = tetrode_separators(spikes)
    if boundaries is None:
        raise ValueError("unit_ids are not tetrode-contiguous; separators would be wrong")
    n_units = int(spikes["unit_id"].n_unique())
    print(f"  {spikes.height} spikes, {n_units} units, {len(tetrode_order)} tetrodes -> {spikes_path}")

    print("[4/5] state definitions")
    (bundle / STATE_DEFINITIONS_FILENAME).write_text(
        (resources.files("tetrode_analyses.assets") / STATE_DEFINITIONS_FILENAME).read_text()
    )

    print("[5/5] manifest")
    manifest = {
        "subject": exp.subject, "experiment": exp.experiment, "session": exp.session,
        "sorting": a.sorting, "sortings_subdir": a.sortings_subdir,
        "source_compressor": a.source_compressor,
        "created": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "lfp": {"fs": float(da.attrs["fs"]), "n_samples": int(da.sizes["time"]),
                "n_channels": int(da.sizes["channel"]), "tetrodes": tetrodes_lfp},
        "spikes": {"n_spikes": int(spikes.height), "n_units": n_units,
                   "tetrode_order": tetrode_order, "separator_boundaries": boundaries},
        "emg_methods": list(EMG_METHODS),
        "state_definitions": STATE_DEFINITIONS_FILENAME,
    }
    (bundle / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    print(f"\nDONE. Fetch on your laptop with:\n  ttm-fetch-bundle viz --subject {exp.subject} "
          f"--experiment {exp.experiment} --dest <drive>/{exp.session}_viz")


def main_loupe(argv: list[str] | None = None) -> None:
    import loupe as lp
    import xarray as xr

    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", type=pathlib.Path, required=True, help="Bundle folder from ttm-fetch-bundle viz")
    p.add_argument("--window-len", type=float, default=10.0, help="Initial window length (s)")
    a = p.parse_args(argv)
    d = a.data_dir
    for name in (SPIKES_FILENAME, LFP_FILENAME, EMG_FILENAME, STATE_DEFINITIONS_FILENAME):
        if not (d / name).exists():
            raise FileNotFoundError(f"{name} not found in {d}; build with ttm-viz-bundle and fetch with ttm-fetch-bundle viz")

    emg = xr.open_zarr(d / EMG_FILENAME)[list(EMG_METHODS)].to_dataarray(dim="method").load()
    da = xr.open_zarr(d / LFP_FILENAME)["lfp"].load()
    tetrodes = [str(t) for t in np.asarray(da["tetrode"].values)]
    da = da.assign_coords(tt_num=("channel", np.array([int(re.findall(r"\d+", t)[0]) for t in tetrodes])))
    spikes = pl.read_parquet(d / SPIKES_FILENAME)
    tt_order, separators = tetrode_separators(spikes)
    palette = tetrode_color_map(tt_order) | tetrode_color_map(tetrodes)
    print(f"EMG {dict(emg.sizes)} | sub-LFP {da.shape} @ {float(da.attrs['fs']):.1f} Hz ({len(tetrodes)} leads) | "
          f"{spikes.height} spikes from {spikes['unit_id'].n_unique()} units across {len(tt_order)} tetrodes")
    if separators is None:
        print("  NOTE: unit_ids not tetrode-contiguous; falling back to split_by='tetrode'.")
        raster = lp.RasterConfig(spikes, time_col="time", order_by="unit_id", split_by="tetrode",
                                 palette=palette, array_name="units")
    else:
        raster = lp.RasterConfig(spikes, time_col="time", order_by="unit_id", split_by=None, hue="tetrode",
                                 palette=palette, horizontal_separators=separators,
                                 separator_params=SEPARATOR_PARAMS, array_name="units")
    lp.view(
        [
            lp.TraceConfig(emg, mode="stacked-subplots", hue="method", palette=EMG_COLORS,
                           array_name="EMG", add_bottom_spine=True),
            lp.TraceConfig(da, mode="dense", order_by="tt_num", descending=True, hue="tetrode",
                           palette=palette, array_name="LFP"),
            raster,
        ],
        window_len=a.window_len,
        state_definitions=str(d / STATE_DEFINITIONS_FILENAME),
    )


if __name__ == "__main__":
    main_build()
