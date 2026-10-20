"""``ttm-convert``: Open Ephys acquisition -> tetrode-aligned 30 kHz Zarr store(s).

Reads the raw session from the ``archive`` remote, concatenates its experiments
on the real sync clock (:func:`ttm.spikeinterface.get_recording`),
writes ``<session>/slice_table.csv`` and one store per ``--compressor``
(:func:`ttm.spikeinterface.convert_recording`), then appends
timings to ``<session>/conversion_results.json``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time

from ttm.cli._common import (
    add_experiment_args,
    experiment_from_args,
    line_buffered,
)
from ttm.spikeinterface import convert_recording, get_recording


def main(argv: list[str] | None = None) -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument(
        "--compressor", action="append", choices=["blosc-zstd", "wavpack"], default=None,
        help="Store(s) to write; repeatable. Default: blosc-zstd (lossless).",
    )
    p.add_argument("--bps", type=float, default=2.25, help="WavPack bits/sample; 0 = lossless")
    p.add_argument("--oe-experiment", type=int, default=None, help="Single OE experiment index; default all")
    p.add_argument("--time-chunk-s", type=float, default=30.0)
    p.add_argument("--inter-tetrode-um", type=float, default=300.0)
    p.add_argument("--n-jobs", type=int, default=32)
    p.add_argument("--overwrite", action="store_true", help="Replace an existing store")
    a = p.parse_args(argv)
    exp = experiment_from_args(a)
    compressors = a.compressor or ["blosc-zstd"]

    print(f"[{time.strftime('%T')}] get_recording({exp.acquisition_dir}, stream {exp.oe_stream_index})")
    t0 = time.perf_counter()
    rec, slice_table = get_recording(
        exp.acquisition_dir,
        oe_experiment_index=a.oe_experiment,
        oe_stream_index=exp.oe_stream_index,
        inter_tetrode_um=a.inter_tetrode_um,
    )
    n_samples = rec.get_num_frames()
    print(
        f"  loaded in {time.perf_counter() - t0:.0f}s | {rec.get_num_channels()} ch "
        f"{n_samples / rec.sampling_frequency / 3600:.2f} h | has_time_vector={rec.has_time_vector(0)}"
    )
    exp.session_dir.mkdir(parents=True, exist_ok=True)
    slice_table.to_csv(exp.slice_table_path, index=False)
    print(f"  wrote {exp.slice_table_path}")

    src_bytes = n_samples * rec.get_num_channels() * rec.get_dtype().itemsize
    results_path = exp.session_dir / "conversion_results.json"
    results = json.loads(results_path.read_text()) if results_path.exists() else []
    for compressor in compressors:
        bps = a.bps if compressor == "wavpack" else None
        dst = exp.store_path(compressor, bps)
        if dst.exists():
            if not a.overwrite:
                print(f"  exists, skipping (pass --overwrite): {dst}")
                continue
            shutil.rmtree(dst)
        print(f"[{time.strftime('%T')}] convert {compressor} -> {dst}")
        t0 = time.perf_counter()
        convert_recording(
            rec, dst, compressor=compressor, bps=a.bps,
            time_chunk_s=a.time_chunk_s, n_jobs=a.n_jobs,
        )
        secs = time.perf_counter() - t0
        disk = sum(f.stat().st_size for f in dst.rglob("*") if f.is_file())
        row = {
            "store": dst.name, "compressor": compressor, "bps": bps, "seconds": round(secs, 1),
            "out_GB": round(disk / 1e9, 1), "ratio": round(src_bytes / disk, 3),
            "source_read_MBps": round(src_bytes / 1e6 / secs, 1),
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        }
        results.append(row)
        results_path.write_text(json.dumps(results, indent=2) + "\n")
        print(f"[{time.strftime('%T')}] DONE {json.dumps(row)}")
    print("ALL DONE")


if __name__ == "__main__":
    main()
