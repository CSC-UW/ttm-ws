"""``ttm-fetch-bundle``: copy a server-built bundle to this machine with rclone.

The server is the ``server`` remote of the per-machine config, an rclone
``remote:path`` pointing at the server's data root (for example
``tononi-2:/nvme/neuropixels/tetrode_data``). ``viz`` copies
``<session>/viz_bundle/`` into ``--dest``; ``curation`` first copies the
curation manifest, then each item it lists to its ``dest_relpath`` under
``--dest`` so the analyzer's relative recording link resolves.
``--no-check-dest`` is the default for a first one-way copy; ``--check``
re-enables checksumming for resumable re-runs.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import subprocess
import sys
import time

from tetrode_analyses import config
from tetrode_analyses.cli._common import add_experiment_args, experiment_from_args
from tetrode_analyses.cli.curation import BUNDLE_MANIFEST_FILENAME, DEFAULT_SORTING


def rclone(args: list[str]) -> None:
    result = subprocess.run(["rclone", *args], check=False)
    if result.returncode != 0:
        sys.exit(f"rclone failed (exit {result.returncode}): rclone {' '.join(args)}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=["viz", "curation"])
    add_experiment_args(p)
    p.add_argument("--dest", type=pathlib.Path, required=True)
    p.add_argument("--remote", default="server", help="Name of the rclone remote in [remotes]")
    p.add_argument("--sorting", default=DEFAULT_SORTING, help="curation: which manifest to fetch")
    p.add_argument("--no-traces", action="store_true", help="curation: skip the recording item")
    p.add_argument("--transfers", type=int, default=16)
    p.add_argument("--checkers", type=int, default=16)
    p.add_argument("--check", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args(argv)
    if not shutil.which("rclone"):
        sys.exit("rclone not found on PATH. Install it: https://rclone.org/install/")
    exp = experiment_from_args(a)
    remote = config.get().remote(a.remote).rstrip("/")
    session_remote = f"{remote}/{exp.session}"
    a.dest.mkdir(parents=True, exist_ok=True)

    flags = [f"--transfers={a.transfers}", f"--checkers={a.checkers}", "--progress"]
    if not a.check:
        flags.append("--no-check-dest")
    if a.dry_run:
        flags.append("--dry-run")
    t0 = time.perf_counter()

    if a.kind == "viz":
        src = f"{session_remote}/viz_bundle"
        print(f"Downloading {src} -> {a.dest}")
        rclone(["copy", src, str(a.dest), *flags])
        print(f"Done in {time.perf_counter() - t0:.1f}s\nView with:\n  ttm-loupe --data-dir {a.dest}")
        return

    local_manifest = a.dest / BUNDLE_MANIFEST_FILENAME
    src_manifest = f"{session_remote}/curation_bundles/{a.sorting}.manifest.json"
    print(f"Fetching manifest {src_manifest}")
    rclone(["copyto", src_manifest, str(local_manifest)])
    manifest = json.loads(local_manifest.read_text())
    skip_recording = a.no_traces or not manifest.get("with_traces", True)
    to_copy = [it for it in manifest["items"] if not (it["role"] == "recording" and skip_recording)]
    print(f"\nBundle: {manifest['subject']}/{manifest['experiment']} sort={manifest['sorting']} "
          f"(with_traces={not skip_recording})")
    for it in to_copy:
        print(f"  - {it['role']:9s} {it.get('size') or '?':>6}  -> {it['dest_relpath']}")
        server_rel = pathlib.Path(it["server_path"]).relative_to(manifest_session_root(manifest, it))
        src = f"{session_remote}/{server_rel}"
        dst = a.dest / it["dest_relpath"]
        dst.parent.mkdir(parents=True, exist_ok=True)
        print(f"\n=== {it['role']} -> {dst} ===")
        rclone(["copy", src, str(dst), *flags])
    print(f"\nDone in {time.perf_counter() - t0:.1f}s\nCurate with:\n  ttm-curate --bundle {a.dest}"
          + (" --no-traces" if skip_recording else ""))


def manifest_session_root(manifest: dict, item: dict) -> pathlib.Path:
    """The server's ``<session>`` directory, recovered from an item's absolute server path."""
    path = pathlib.Path(item["server_path"])
    for parent in [path, *path.parents]:
        if parent.name == manifest["session"]:
            return parent
    raise ValueError(f"{item['server_path']} is not under a {manifest['session']!r} directory")


if __name__ == "__main__":
    main()
