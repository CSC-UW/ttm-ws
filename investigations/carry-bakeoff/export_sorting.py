"""Save an arm's accepted spikes as a SpikeInterface sorting with a ``group`` property per unit.

The sorting goes to ``<run dir>/<name>`` (default: ``<arm dir>/sorting``); with ``--run-dir`` it can be
placed under the session's sortings beside the June runs, where ``bad_periods_compare.py`` and the
purity scripts read it.

Run: .venv/bin/python investigations/carry-bakeoff/export_sorting.py ARM [--run-dir DIR --name NAME]
"""
from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm")
    ap.add_argument("--run-dir", type=pathlib.Path, default=None)
    ap.add_argument("--name", default="sorting")
    ap.add_argument("--root", type=pathlib.Path, default=bk.BAKE / "arms")
    a = ap.parse_args()
    d = a.root / a.arm
    files = sorted(d.glob("spikes_w*.npz"))
    fr, un = [], []
    for f in files:
        z = np.load(f)
        acc = z["reason"] == 0
        fr.append(z["frame"][acc])
        un.append(z["unit"][acc])
    fr, un = np.concatenate(fr), np.concatenate(un)
    last = bk.Bank.load(sorted(d.glob("bank_w*.npz"))[-1])
    o = np.argsort(fr, kind="stable")
    s = si.NumpySorting.from_samples_and_labels([fr[o]], [un[o]], sampling_frequency=bk.FS, unit_ids=last.ids)
    s.set_property("group", last.group)
    out = (a.run_dir or d) / a.name
    s.save(folder=out, overwrite=True)
    print(f"{len(files)} windows, {fr.size:,} spikes, {len(last.ids)} units -> {out}")


if __name__ == "__main__":
    main()
