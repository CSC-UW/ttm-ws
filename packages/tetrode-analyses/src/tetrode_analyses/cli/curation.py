"""Curation in spikeinterface-gui, on the server or over a downloaded bundle.

``ttm-curation-manifest`` (server) writes
``<session>/curation_bundles/<sorting>.manifest.json`` listing the analyzer and,
unless ``--no-traces``, the recording the sort was produced from, each with the
relative path it must land at so the analyzer's relative recording link keeps
resolving. ``ttm-fetch-bundle curation`` copies those items to a laptop.

``ttm-curate`` launches the GUI: on the server with ``--analyzer`` (default: the
chunk-tracked ``tracked_48h/analyzer_clustered.zarr`` of the selected
experiment), on a laptop with ``--bundle DIR``. ``--mode web`` serves the Panel
backend on ``localhost:<port>`` without opening a browser (reach it over an SSH
tunnel); ``--mode desktop`` needs a display. ``--style grahams_curation`` applies
the layout, per-view settings and unit-list columns shipped in
``tetrode_analyses/assets``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import pathlib
import subprocess
import sys
from importlib import resources

from tetrode_analyses.cli._common import add_experiment_args, experiment_from_args

DEFAULT_SORTINGS_SUBDIR = "sortings_seed42_pcafix"
DEFAULT_SORTING = "tracked_48h"
DEFAULT_ANALYZER_NAME = "analyzer_clustered.zarr"
RECORDING_COMPRESSOR = "blosc-zstd"
BUNDLE_MANIFEST_FILENAME = "curation_manifest.json"
DEFAULT_PORT = 8000

STYLES = {
    "grahams_curation": {
        "layout_file": "grahams_curation_layout.json",
        "settings_file": "grahams_curation_settings.json",
        "displayed_unit_properties": [
            "group", "tier", "n_chunks", "track_hours", "n_windows", "identity_min_cos",
            "is_reseeded", "unitrefine_label", "unitrefine_neural_prob", "firing_rate",
            "rp_contamination", "sliding_rp_violation",
        ],
    },
}


def _asset(name: str) -> str:
    return (resources.files("tetrode_analyses.assets") / name).read_text()


def style_kwargs(style_name: str | None) -> dict:
    """``run_mainwindow`` kwargs for a named style (empty for ``None``)."""
    if style_name is None:
        return {}
    style = STYLES[style_name]
    kw = {
        "layout": json.loads(_asset(style["layout_file"])),
        "user_settings": json.loads(_asset(style["settings_file"])),
    }
    if style.get("displayed_unit_properties"):
        kw["displayed_unit_properties"] = style["displayed_unit_properties"]
    return kw


def du_size(path: pathlib.Path, timeout: float = 120.0) -> str | None:
    try:
        out = subprocess.run(["du", "-sh", str(path)], capture_output=True, text=True,
                             timeout=timeout, check=True)
        return out.stdout.split("\t", 1)[0].strip() or None
    except (subprocess.SubprocessError, OSError):
        return None


def main_manifest(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--sorting", default=DEFAULT_SORTING)
    p.add_argument("--sortings-subdir", default=DEFAULT_SORTINGS_SUBDIR)
    p.add_argument("--analyzer-name", default=DEFAULT_ANALYZER_NAME)
    p.add_argument("--no-traces", action="store_true", help="Omit the recording (no raw-traces view)")
    p.add_argument("--no-size", action="store_true", help="Skip du size lookup")
    p.add_argument("--manifest-out", type=pathlib.Path, default=None)
    a = p.parse_args(argv)
    exp = experiment_from_args(a)
    with_traces = not a.no_traces

    analyzer_path = exp.sorting_dir(a.sorting, a.sortings_subdir) / a.analyzer_name
    if not analyzer_path.exists():
        raise FileNotFoundError(f"analyzer not found: {analyzer_path}")
    recording_path = exp.store_path(RECORDING_COMPRESSOR)
    if with_traces and not recording_path.exists():
        raise FileNotFoundError(f"recording not found: {recording_path} (use --no-traces)")

    def item(role, path, dest_relpath, required):
        return {"role": role, "server_path": str(path), "dest_relpath": dest_relpath,
                "size": None if a.no_size else du_size(path), "required": required}

    items = [item("analyzer", analyzer_path, str(analyzer_path.relative_to(exp.session_dir)), True)]
    if with_traces:
        items.append(item("recording", recording_path, recording_path.name, False))
    manifest = {
        "kind": "curation_manifest",
        "subject": exp.subject, "experiment": exp.experiment, "session": exp.session,
        "sorting": a.sorting, "analyzer_name": a.analyzer_name, "sortings_subdir": a.sortings_subdir,
        "recording_compressor": RECORDING_COMPRESSOR, "with_traces": with_traces,
        "created": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "items": items,
    }
    out = a.manifest_out or (exp.curation_bundles_dir / f"{a.sorting}.manifest.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    print(f"\nWrote manifest -> {out}")
    print(f"Fetch it on your laptop with:\n  ttm-fetch-bundle curation --subject {exp.subject} "
          f"--experiment {exp.experiment} --sorting {a.sorting} --dest <drive>/{exp.session}")


def _analyzer_from_bundle(bundle: pathlib.Path, no_traces: bool):
    manifest_path = bundle / BUNDLE_MANIFEST_FILENAME
    if not manifest_path.exists():
        raise FileNotFoundError(f"{BUNDLE_MANIFEST_FILENAME} not in {bundle}; run ttm-fetch-bundle curation first")
    manifest = json.loads(manifest_path.read_text())
    rel = next(it["dest_relpath"] for it in manifest["items"] if it["role"] == "analyzer")
    want_traces = not no_traces and manifest.get("with_traces", True)
    return bundle / rel, want_traces, manifest["sorting"]


def main_curate(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--bundle", type=pathlib.Path, help="Folder populated by ttm-fetch-bundle curation (laptop)")
    p.add_argument("--analyzer", type=pathlib.Path, help="Analyzer to curate (server); default tracked_48h/analyzer_clustered.zarr")
    p.add_argument("--sorting", default=DEFAULT_SORTING)
    p.add_argument("--sortings-subdir", default=DEFAULT_SORTINGS_SUBDIR)
    p.add_argument("--analyzer-name", default=DEFAULT_ANALYZER_NAME)
    p.add_argument("--mode", choices=["web", "desktop"], default="web")
    p.add_argument("--port", type=int, default=int(os.environ.get("SIGUI_PORT", DEFAULT_PORT)))
    p.add_argument("--style", choices=sorted(STYLES), default=None)
    p.add_argument("--no-traces", action="store_true", help="Curate without the raw-traces view")
    a = p.parse_args(argv)

    if a.bundle:
        analyzer_path, want_traces, sorting = _analyzer_from_bundle(a.bundle, a.no_traces)
    elif a.analyzer:
        analyzer_path, want_traces, sorting = a.analyzer, not a.no_traces, a.analyzer.parent.name
    else:
        exp = experiment_from_args(a)
        analyzer_path = exp.sorting_dir(a.sorting, a.sortings_subdir) / a.analyzer_name
        want_traces, sorting = not a.no_traces, a.sorting
    if not analyzer_path.exists():
        raise FileNotFoundError(f"analyzer not found: {analyzer_path}")

    sys.argv = sys.argv[:1]  # keep Qt from reading --style as its own option
    if a.mode == "web":
        os.environ.setdefault("BROWSER", "echo")

    import spikeinterface as si
    import spikeinterface_gui as sg

    analyzer = si.load_sorting_analyzer(analyzer_path, load_extensions=False)
    if want_traces and not analyzer.has_recording():
        raise FileNotFoundError(
            "Traces requested but the recording did not resolve next to the analyzer "
            f"({RECORDING_COMPRESSOR} store expected at the session/bundle root). "
            "Pass --no-traces to curate without the raw-traces view."
        )
    print(f"Curation: {sorting} | {analyzer.sorting.get_num_units()} units | "
          f"traces={'on' if want_traces else 'off'} | mode={a.mode}", flush=True)
    skw = style_kwargs(a.style)
    if a.mode == "web":
        print(f"\nServing curation GUI at  http://localhost:{a.port}/   (Ctrl-C to stop)")
        print(f"Remote? tunnel first:     ssh -N -L {a.port}:localhost:{a.port} <you>@<server>\n", flush=True)
        sg.run_mainwindow(analyzer, mode="web", curation=True, with_traces=want_traces,
                          address="localhost", port=a.port,
                          panel_start_server_kwargs={"show": False}, **skw)
    else:
        print("Launching desktop (Qt) curation GUI.", flush=True)
        sg.run_mainwindow(analyzer, mode="desktop", curation=True, with_traces=want_traces, **skw)


if __name__ == "__main__":
    main_curate()
