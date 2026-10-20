"""``ttm-experiment-params``: write ``<subject>/<experiment>/experiment_params.json``.

The manifest supplies everything except the session anchor, which is read from
the session's ``slice_table.csv`` when it exists (written by ``ttm-convert``)
and otherwise from the raw acquisition's first ``sync_messages.txt``.
"""

from __future__ import annotations

import argparse

from ttm import experiment as ep
from ttm.cli._common import add_experiment_args, experiment_from_args


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--from-acquisition", action="store_true",
                   help="Read the anchor from the raw acquisition even if slice_table.csv exists")
    a = p.parse_args(argv)
    exp = experiment_from_args(a)

    source = exp.slice_table_path
    if a.from_acquisition or not source.exists():
        source = exp.acquisition_dir
    t0_unix = ep.read_session_t0_unix(source)
    params = ep.ExperimentParams.from_experiment(exp, t0_unix)
    out = ep.save_experiment_params(params, exp.experiment_params_path)

    print(f"Wrote {out}")
    print(f"Session anchor: {params.t0_datetime.isoformat()} (t0_unix={t0_unix}, from {source})")
    intervals, labels = ep.get_light_dark_periods(params)
    print("Light/dark periods (session seconds):")
    for (s, e), lab in zip(intervals, labels):
        print(f"  {lab:>3}  {s:>10.0f} -> {e:>10.0f}")
    print(f"Deprivation (session seconds): {ep.get_deprivation_period(params)}")


if __name__ == "__main__":
    main()
