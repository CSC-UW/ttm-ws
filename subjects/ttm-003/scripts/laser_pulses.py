"""The laser pulses of a session in the form ``ttm-sort --pulses`` reads.

From ``results/stimulation_pulses_<session>.csv`` and ``results/stimulation_blocks_<session>.csv``
(``stim_trains.py``), the pulses of delivered trains that are not glitches, as the carry-forward
bake-off selected them (``investigations/carry-bakeoff/_bake.py``, ``laser_pulses``). Writes
``results/laser_pulses_<session>.csv`` with ``onset_s`` and ``offset_s`` in session seconds.

Usage: python subjects/ttm-003/scripts/laser_pulses.py --subject TTM-003 --experiment TTM-003-48H
"""

from __future__ import annotations

import argparse
import pathlib

import pandas as pd

from ttm.cli._common import add_experiment_args, experiment_from_args

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    exp = experiment_from_args(ap.parse_args())
    blocks = pd.read_csv(RESULTS / f"stimulation_blocks_{exp.session}.csv")
    delivered = set(blocks.loc[blocks.delivered, "ttl_train"])
    p = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    p = p[~p.glitch & p.train.isin(delivered)].sort_values("onset_s", kind="stable")
    out = RESULTS / f"laser_pulses_{exp.session}.csv"
    p[["onset_s", "offset_s"]].to_csv(out, index=False)
    print(f"{len(p)} pulses of {len(delivered)} delivered trains -> {out}")


if __name__ == "__main__":
    main()
