"""Snapshot and compare the derived products of one experiment.

``snapshot`` writes a JSON summary of the products that the pipeline
entry points produce for an experiment: ``experiment_params.json``, the
production sorting's unit set, and the SWA / EMG Zarr stores (shape,
attributes, per-tetrode statistics, and a hash of the raw values).
``compare`` diffs two such summaries and exits non-zero on any
difference. Both are cheap enough to run before and after a change to
the library or a re-run of an entry point.

    uv run python scripts/snapshot_derived.py snapshot \
        --subject TTM-001 --experiment TTM-NOD -o docs/baselines/x.json
    uv run python scripts/snapshot_derived.py compare a.json b.json
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import pathlib
import subprocess
import sys
from typing import Any

import numpy as np
import xarray as xr
from tetrode_analyses.cli._common import add_experiment_args, experiment_from_args
from tetrode_analyses.subject import Experiment

DEFAULT_SORTING = "track_eval/mp_long_s2000_d170000/assembled_prod"


def _sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _array_hash(a: np.ndarray) -> str:
    return _sha256(np.ascontiguousarray(a).tobytes())


def _jsonable(v: Any) -> Any:
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    return v


def _stats(a: np.ndarray) -> dict[str, Any]:
    a = np.asarray(a, dtype=np.float64)
    finite = np.isfinite(a)
    out: dict[str, Any] = {"n": int(a.size), "n_nonfinite": int((~finite).sum())}
    if finite.any():
        f = a[finite]
        out.update(
            mean=float(f.mean()),
            std=float(f.std()),
            min=float(f.min()),
            max=float(f.max()),
        )
    return out


def _dataarray_summary(da: xr.DataArray, per: str = "tetrode") -> dict[str, Any]:
    values = da.values
    out: dict[str, Any] = {
        "dims": list(da.dims),
        "shape": list(da.shape),
        "dtype": str(da.dtype),
        "attrs": _jsonable(dict(da.attrs)),
        "sha256": _array_hash(values),
        "all": _stats(values),
    }
    if "time" in da.coords:
        t = np.asarray(da["time"].values)
        out["time"] = {
            "n": int(t.size),
            "first": _jsonable(t[0]),
            "last": _jsonable(t[-1]),
        }
    if per in da.dims:
        axis = da.dims.index(per)
        out[per] = {
            str(label): _stats(np.take(values, i, axis=axis))
            for i, label in enumerate(np.asarray(da[per].values))
        }
    return out


def _dataset_summary(ds: xr.Dataset) -> dict[str, Any]:
    return {
        "attrs": _jsonable(dict(ds.attrs)),
        "data_vars": {
            name: _dataarray_summary(ds[name]) for name in sorted(ds.data_vars)
        },
    }


def _zarr_summary(
    path: pathlib.Path, variable: str | None = None
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    ds = xr.open_zarr(str(path))
    if variable is not None:
        return _dataarray_summary(ds[variable])
    return _dataset_summary(ds)


def _params_summary(path: pathlib.Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    content = json.loads(path.read_text())
    canonical = json.dumps(content, sort_keys=True, separators=(",", ":")).encode()
    return {"sha256": _sha256(canonical), "content": content}


def _sorting_summary(path: pathlib.Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    import spikeinterface as si

    s = si.load(path)
    unit_ids = np.asarray(s.unit_ids)
    counts = s.count_num_spikes_per_unit(outputs="dict")
    per_unit = {str(u): int(counts[u]) for u in unit_ids}
    out: dict[str, Any] = {
        "num_units": int(s.get_num_units()),
        "num_spikes": int(sum(per_unit.values())),
        "unit_ids_sha256": _sha256(json.dumps([str(u) for u in unit_ids]).encode()),
        "spike_counts_sha256": _sha256(json.dumps(per_unit, sort_keys=True).encode()),
        "spike_counts": per_unit,
    }
    if "is_mua" in s.get_property_keys():
        im = np.asarray(s.get_property("is_mua"), dtype=bool)
        out["n_sua"] = int((~im).sum())
        out["n_mua"] = int(im.sum())
    if "group" in s.get_property_keys():
        g = np.asarray(s.get_property("group"))
        out["units_per_group"] = {
            str(k): int(v) for k, v in zip(*np.unique(g, return_counts=True))
        }
    return out


def _git_head(path: pathlib.Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
        ).strip()
    except subprocess.CalledProcessError, OSError:
        return None


def snapshot(exp: Experiment, sorting: str) -> dict[str, Any]:
    workspace = pathlib.Path(__file__).resolve().parents[1]
    return {
        "meta": {
            "subject": exp.subject,
            "experiment": exp.experiment,
            "session": exp.session,
            "taken": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
            "workspace_commit": _git_head(workspace),
            "sorting": sorting,
        },
        "experiment_params": _params_summary(exp.experiment_params_path),
        "sorting": _sorting_summary(exp.sorting_dir(sorting)),
        "idelta": _zarr_summary(exp.idelta_path, variable="pwr"),
        "stft_bandpowers": _zarr_summary(exp.stft_bandpowers_path),
        "synthetic_emg": _zarr_summary(exp.emg_path),
    }


IGNORED = {("meta", "taken"), ("meta", "workspace_commit")}


def _flatten(d: Any, prefix: tuple[str, ...] = ()) -> dict[tuple[str, ...], Any]:
    if isinstance(d, dict):
        out: dict[tuple[str, ...], Any] = {}
        for k, v in d.items():
            out.update(_flatten(v, prefix + (str(k),)))
        return out
    return {prefix: d}


def compare(
    a: dict[str, Any], b: dict[str, Any], rtol: float = 0.0
) -> tuple[list[str], list[str]]:
    """Return (failures, notes). With ``rtol`` > 0, floats within that relative
    tolerance match and differing value hashes are notes rather than failures."""
    fa, fb = _flatten(a), _flatten(b)
    failures, notes = [], []
    for key in sorted(set(fa) | set(fb)):
        if key[:2] in IGNORED:
            continue
        name = ".".join(key)
        if key not in fa:
            failures.append(f"+ {name} = {fb[key]!r}")
        elif key not in fb:
            failures.append(f"- {name} = {fa[key]!r}")
        elif fa[key] != fb[key]:
            x, y = fa[key], fb[key]
            if rtol > 0 and key[-1] == "sha256" and key[0] != "experiment_params":
                notes.append(f"~ {name}: values not bit-identical")
            elif (
                rtol > 0
                and isinstance(x, float)
                and isinstance(y, float)
                and math.isclose(x, y, rel_tol=rtol)
            ):
                notes.append(f"~ {name}: {x!r} -> {y!r} (within rtol={rtol:g})")
            else:
                failures.append(f"! {name}: {x!r} -> {y!r}")
    return failures, notes


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    ps = sub.add_parser(
        "snapshot", help="write a JSON summary of an experiment's derived products"
    )
    add_experiment_args(ps)
    ps.add_argument(
        "--sorting",
        default=DEFAULT_SORTING,
        help="sorting path under sortings_seed42_pcafix/ (default: %(default)s)",
    )
    ps.add_argument("-o", "--output", type=pathlib.Path, required=True)
    pc = sub.add_parser("compare", help="diff two summaries; exit 1 on any difference")
    pc.add_argument("a", type=pathlib.Path)
    pc.add_argument("b", type=pathlib.Path)
    pc.add_argument(
        "--rtol",
        type=float,
        default=0.0,
        help="relative tolerance for float statistics (default: exact)",
    )
    args = p.parse_args(argv)

    if args.cmd == "snapshot":
        exp = experiment_from_args(args)
        summary = snapshot(exp, args.sorting)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n")
        present = [k for k, v in summary.items() if k != "meta" and v is not None]
        print(f"Wrote {args.output} ({', '.join(present)})")
        return

    failures, notes = compare(
        json.loads(args.a.read_text()), json.loads(args.b.read_text()), rtol=args.rtol
    )
    if notes or failures:
        print("\n".join(notes + failures))
    if failures:
        sys.exit(1)
    print(
        f"{args.a} and {args.b} match"
        + (f" (rtol={args.rtol:g}, {len(notes)} tolerated)" if notes else "")
    )


if __name__ == "__main__":
    main()
