from __future__ import annotations

import argparse
import sys

from ttm.subject import Experiment, resolve_experiment


def add_experiment_args(parser: argparse.ArgumentParser) -> None:
    g = parser.add_argument_group("experiment selection (manifest > subject/experiment > $TTM_SUBJECT/$TTM_EXPERIMENT)")
    g.add_argument("--manifest", help="Path to subjects/<id>/experiments/<NAME>.toml")
    g.add_argument("--subject", help="Subject id, e.g. TTM-001")
    g.add_argument("--experiment", help="Experiment name, e.g. TTM-NOD")


def experiment_from_args(args: argparse.Namespace) -> Experiment:
    exp = resolve_experiment(
        manifest=args.manifest, subject=args.subject, experiment=args.experiment
    )
    print(f"{exp.subject}/{exp.experiment}: session {exp.session} ({exp.manifest_path})", flush=True)
    return exp


def line_buffered() -> None:
    sys.stdout.reconfigure(line_buffering=True)
