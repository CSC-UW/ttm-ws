"""Which unit carries each light-driven cell's response, block by block, in the 48 h sortings.

Reads the per-block tables of ``subjects/ttm-003/scripts/tagged_cells_blocks.py`` (tags
``masked_c12``, ``carry48_joint``, ``carry48_pt``, ``carry48_joint_rules``, ``carry48_pt_rc``) and prints, per tetrode with a tagged cell, the
spikes per second 1 to 10 ms after onset of every unit that reaches 50 Hz in some block. Writes
``results/tagged_units_onset.csv``.

``--waveforms`` adds the window-0 check on TT4: the mean of each unit's spikes 1 to 10 ms after
the onsets against the mean of its own and the other units' spikes away from the pulses (cosine,
plain traces).

Run: .venv/bin/python investigations/carry-bakeoff/tagged_units.py [--waveforms]
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import _bake as bk  # noqa: E402

HERE = pathlib.Path(__file__).parent
RESULTS = HERE.parents[1] / "subjects" / "ttm-003" / "results"
TAGS = {"production": "masked_c12", "joint_whitened": "carry48_joint", "per_tetrode": "carry48_pt",
        "joint_rules": "carry48_joint_rules", "per_tetrode_rc": "carry48_pt_rc"}
SORTS = {"production": "mp_baseline_masked", "joint_whitened": "carry_R48_joint", "per_tetrode": "carry_R48_pt",
         "joint_rules": "carry_R48_joint_rules", "per_tetrode_rc": "carry_R48_pt_rc"}
TETRODES = ("TT4", "TT6", "TT10", "TT11")


def onset_tables() -> pd.DataFrame:
    parts = []
    for name, tag in TAGS.items():
        d = pd.read_csv(RESULTS / f"tagged_blocks_TTM-003-48H_{tag}.csv")
        d = d[d.tetrode.isin(TETRODES)]
        hot = d.groupby("unit").onset_1_10ms_hz.max()
        parts.append(d[d.unit.isin(hot.index[hot >= 50])].assign(sorting=name))
    df = pd.concat(parts)[["sorting", "tetrode", "unit", "block", "onset_h", "onset_1_10ms_hz", "base_hz",
                           "in_pulse_hz", "reliability_1_10ms"]]
    df.to_csv(HERE / "results" / "tagged_units_onset.csv", index=False)
    pd.set_option("display.width", 250)
    for tet in TETRODES:
        for name in TAGS:
            d = df[(df.tetrode == tet) & (df.sorting == name)]
            print(f"== {tet} {name}: onset 1-10 ms, Hz")
            print(d.pivot_table(index="unit", columns="block", values="onset_1_10ms_hz").round(0).astype("Int64")
                  .to_string())
    return df


def waveforms(g: int = 3, units=(5, 6, 7)) -> None:
    import spikeinterface as si

    from ttm import capture as cp

    rec = si.load(bk.BINARY)
    on, off = bk.laser_pulses()
    a, b = 0, bk.WIN
    m = (on > a) & (on < b)
    on, off = on[m], off[m]
    x = cp.load_tetrode_traces(rec, g, a, b)
    spont, evok = {}, {}
    for name, p in SORTS.items():
        s = si.load(bk.SORTINGS / p / "assembled_reseed_c12")
        for u in units:
            f = np.sort(s.get_unit_spike_train(u, start_frame=a + bk.NB, end_frame=b - bk.NA))
            ev = bk.in_intervals(f, on + 30, on + 300)
            quiet = ~bk.in_intervals(f, on - 3000, off + 3000)
            spont[(name, u)] = cp.snippets(x, f[quiet][:2000] - a).mean(0)
            print(f"{name} unit {u}: {f.size} spikes in window 0, {ev.sum()} 1-10 ms after the {on.size} onsets")
            if ev.sum() >= 10:
                evok[(name, u)] = cp.snippets(x, f[ev] - a).mean(0)

    def cos(A, B):
        return float((A * B).sum() / np.sqrt((A * A).sum() * (B * B).sum()))

    keys = list(spont)
    print("cosine of the evoked means (rows) to the spontaneous means (columns)")
    print(pd.DataFrame([[cos(E, spont[k]) for k in keys] for E in evok.values()],
                       index=[f"{n} {u}" for n, u in evok], columns=[f"{n} {u}" for n, u in keys]).round(3)
          .to_string())


if __name__ == "__main__":
    onset_tables()
    if "--waveforms" in sys.argv:
        waveforms()
