"""Template-shaped events inside laser pulses against sham times, by size, for units of ``opto_units.py``.

A light-driven cell firing fast through a pulse can fire spikes too small for the sort. For each requested
unit of an ``opto_units.py`` stretch, this takes the unit's mean waveform on its tetrode (1 ms before to 2 ms
after, from 1,000 of its spikes outside pulses and the 100 ms after them, seed 0), slides it over the run's
preprocessed binary from 15 ms after each pulse onset of the stretch to the pulse's median length, and keeps
the local maxima (at least 1 ms apart) of the projection's amplitude (in units of the template) with a cosine
of 0.5 or more to the template. The same at sham times 500 and 250 ms before each onset gives the background.
An excess of small events inside pulses that a non-driven unit's template does not also show would be small
spikes of the driven cell.

Writes ``results/opto_small_spikes_<E>.csv``: per stretch, unit and window kind, events per second in
amplitude bins.

Usage: python subjects/ttm-003/scripts/opto_small_spikes.py --subject TTM-003 --experiment TTM-003-48H-2
           --units day1:6 sd:7 day1:9 day1:12 day1:50
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np
import pandas as pd
import scipy.signal as sps
import spikeinterface as si
from ttm import carry as cr
from ttm.cli._common import add_experiment_args, experiment_from_args
from ttm.sorting import bad_channel_ids

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
BINS = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.5, 3.0]
NB, NA = 30, 60


def events(rec, ids, T, times, a_s, b_s, fs):
    TT = float((T ** 2).sum())
    Tn = np.sqrt(TT)
    out = []
    for t0 in times:
        f0, f1 = round((t0 + a_s) * fs), round((t0 + b_s) * fs)
        x = rec.get_traces(start_frame=f0 - NB, end_frame=f1 + NA, channel_ids=ids, return_in_uV=True).astype(np.float64)
        proj = sum(np.correlate(x[:, k], T[:, k], mode="valid") for k in range(len(ids))) / TT
        energy = np.sqrt(np.convolve((x ** 2).sum(1), np.ones(NB + NA), mode="valid"))
        cos = proj * TT / (energy * Tn + 1e-9)
        pk, _ = sps.find_peaks(proj, height=BINS[0], distance=30)
        pk = pk[pk < f1 - f0]
        out.append(np.c_[proj[pk], cos[pk]])
    return np.concatenate(out) if out else np.zeros((0, 2))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(ap)
    ap.add_argument("--units", nargs="+", required=True, help="stretch:unit pairs, the stretch an opto_units.py --name")
    ap.add_argument("--run", default="main")
    a = ap.parse_args()
    exp = experiment_from_args(a)
    store = si.read_zarr(str(exp.store_path()))
    drop = bad_channel_ids(store, exp.badChannels, exp.badTetrodes)
    rec = cr.open_sorting_recording(exp.sorting_dir(a.run, subdir="sortings_carry") / "binary", drop)
    fs = rec.get_sampling_frequency()
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    blocks = pd.read_csv(RESULTS / f"stimulation_blocks_{exp.session}.csv")
    pulses = pulses[~pulses.glitch & pulses.train.isin(blocks.loc[blocks.delivered, "ttl_train"])]
    rows = []
    for spec in a.units:
        name, unit = spec.split(":")
        unit = int(unit)
        info = json.loads((RESULTS / f"opto_units_{exp.experiment}_{name}.json").read_text())
        h0, h1 = info["start_h"], info["end_h"]
        fa = int(h0 * 3600 * fs)
        sorting = cr.sort_window(None, exp.experiment_dir / "opto_units" / name, cr.CarryParams())
        g = int(np.asarray(sorting.get_property("group"))[list(sorting.unit_ids).index(unit)])
        ids = [rec.channel_ids[k] for k in cr.cp.tetrode_channel_indices(rec, g)]
        pp = pulses[(pulses.onset_s >= h0 * 3600 + 1) & (pulses.onset_s < h1 * 3600 - 1)].sort_values("onset_s")
        on, off = pp.onset_s.to_numpy(), pp.offset_s.to_numpy()
        fr = sorting.get_unit_spike_train(unit) + fa
        t = fr / fs
        j = np.searchsorted(on, t, "right") - 1
        near = (j >= 0) & (t < off[np.maximum(j, 0)] + 0.1)
        far = np.flatnonzero(~near & (fr >= NB) & (fr + NA <= rec.get_num_frames()))
        pick = np.random.default_rng(0).choice(far, min(1000, far.size), replace=False)
        T = np.stack([rec.get_traces(start_frame=int(f) - NB, end_frame=int(f) + NA, channel_ids=ids, return_in_uV=True)
                      for f in fr[pick]]).astype(np.float64).mean(0)
        b_s = float(np.median(off - on))
        dur = on.size * (b_s - 0.015)
        for kind, times in (("pulse", on), ("sham_500ms", on - 0.5), ("sham_250ms", on - 0.25)):
            E = events(rec, ids, T, times, 0.015, b_s, fs)
            E = E[E[:, 1] >= 0.5]
            h = np.histogram(E[:, 0], BINS)[0] / dur
            rows.append({"stretch": name, "unit": unit, "tetrode": f"TT{g + 1}", "kind": kind, "pulses": int(on.size),
                         "window_ms": f"15-{b_s * 1000:.0f}",
                         **{f"a_{lo}_{hi}": round(float(v), 2) for lo, hi, v in zip(BINS[:-1], BINS[1:], h)}})
        print(pd.DataFrame(rows[-3:]).to_string(index=False))
    pd.DataFrame(rows).to_csv(RESULTS / f"opto_small_spikes_{exp.experiment}.csv", index=False)


if __name__ == "__main__":
    main()
