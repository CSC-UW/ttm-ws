"""Are the October 16 tagged cells present at each 48 h block? Capture with the OPTO2 templates.

Independent of any 48 h sorting: for each delivered block and each named
tetrode, the baseline binary's traces from 60 s before the train to 60 s
after are captured with ``ttm.capture`` using the OPTO2 capture templates
of that tetrode (``capture/templates_<TT>.npz`` plus
``capture_tt10seed`` for TT10), and every template's baseline (-50 to
0 ms), onset (1 to 10 ms), in-pulse (10 ms to offset) and post (0 to
50 ms) rates, reliability and spike count in the window are tabulated.
A template whose cell is present shows its October 16 signature; a
template whose cell is gone collects only hash. Writes
``results/tagged_xsession_<E>.csv``, a JSON summary of the tagged
templates and ``figures/tagged_xsession.png``.

Usage: python subjects/ttm-003/scripts/tagged_cells_xsession.py --subject TTM-003 --experiment TTM-003-48H
           [--tetrodes TT6 TT10 TT4 TT11] [--margin-s 60] [--workers 24]
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spikeinterface as si
from joblib import Parallel, delayed

from ttm import capture as cp
from ttm import config
from ttm.cli._common import add_experiment_args, experiment_from_args, line_buffered

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
RESULTS = SUBJECT_DIR / "results"
FS = 30000.0
OPTO2_SESSION = "2026-10-16_14-44-22"
TAGGED = {"TT6_cell": ("TT6", 22), "TT10_driven": ("TT10", 1004), "TT4_candidate": ("TT4", 3), "TT11_candidate": ("TT11", 39)}


def load_blocks(session: str) -> pd.DataFrame:
    df = pd.read_csv(RESULTS / f"stimulation_blocks_{session}.csv")
    df = df[df.delivered & (df.n_pulses >= 10)]
    return df.reset_index(drop=True)


def templates_for(tt: str, dirs: list[pathlib.Path], wires_here: list[str]):
    T, ids = [], []
    for d in dirs:
        f = d / f"templates_{tt}.npz"
        if f.exists():
            z = np.load(f)
            order = [[str(w) for w in z["wires"]].index(w) for w in wires_here]
            T.append(z["templates"][:, :, order]); ids += [int(u) for u in z["unit_ids"]]
    return (np.concatenate(T), ids) if T else (np.zeros((0, cp.N_SAMP, 4), np.float32), [])


def rate_in(ts: np.ndarray, on: np.ndarray, lo: float, hi: float) -> float:
    i0 = np.searchsorted(ts, on + lo); i1 = np.searchsorted(ts, on + hi)
    return float((i1 - i0).sum() / (len(on) * (hi - lo)))


def one(binary: str, group: int, blk, pulses: pd.DataFrame, T: np.ndarray, ids: list, margin_s: float,
        thresh_mad: float, gate_r: float) -> list[dict]:
    rec = si.load(binary)
    gain = float(rec.get_property("gain_to_uV")[0])
    a = max(int((blk.onset_s - margin_s) * FS), 0)
    b = min(int((blk.offset_s + margin_s) * FS), rec.get_num_frames())
    x = cp.load_tetrode_traces(rec, group, a, b)
    C = cp.capture_tetrode(x, T, ids, threshold_mad=thresh_mad, gate_r=gate_r)
    pp = pulses[pulses.train == blk.ttl_train]
    on = pp.onset_s.values - a / FS
    dur = float(np.median(pp.offset_s.values - pp.onset_s.values))
    rows = []
    for i, u in enumerate(ids):
        fr = C.train(i)
        ts = fr / FS
        rel = [((ts >= e + 0.001) & (ts < e + 0.010)).any() for e in on]
        rows.append({"block": int(blk.ttl_train), "onset_h": round(float(blk.onset_h), 4), "tetrode": f"TT{group + 1}",
                     "opto2_unit": u, "template_trough_uV": round(float(T[i].min() * gain), 1),
                     "n_window": int(fr.size), "window_s": round((b - a) / FS, 1),
                     "rate_window_hz": round(fr.size / ((b - a) / FS), 3),
                     "base_hz": rate_in(ts, on, -0.05, 0.0), "onset_1_10ms_hz": rate_in(ts, on, 0.001, 0.010),
                     "in_pulse_hz": rate_in(ts, on, 0.010, dur), "post_0_50ms_hz": rate_in(ts, on + dur, 0.0, 0.05),
                     "reliability_1_10ms": float(np.mean(rel)) if len(rel) else 0.0,
                     "amp_mad_median": round(float(np.median(C.amps()[i])), 2) if fr.size else None,
                     "r_median": round(float(np.median(C.rs()[i])), 3) if fr.size else None,
                     "noise_uV_max": round(float(C.noise.max() * gain), 2)})
    return rows


def main() -> None:
    line_buffered()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_experiment_args(p)
    p.add_argument("--tetrodes", nargs="+", default=["TT6", "TT10", "TT4", "TT11"])
    p.add_argument("--binary", default=None)
    p.add_argument("--margin-s", type=float, default=60.0)
    p.add_argument("--thresh-mad", type=float, default=3.5)
    p.add_argument("--gate-r", type=float, default=0.6)
    p.add_argument("--workers", type=int, default=24)
    a = p.parse_args()
    exp = experiment_from_args(a)
    data_root = config.get().data_root
    binary = a.binary or str(exp.sorting_dir("mp_baseline") / "binary")
    dirs = [data_root / OPTO2_SESSION / "sortings_opto" / "capture", data_root / OPTO2_SESSION / "sortings_opto" / "capture_tt10seed"]
    rec = si.load(binary)
    ch_ids = np.asarray(rec.channel_ids)
    blocks = load_blocks(exp.session)
    pulses = pd.read_csv(RESULTS / f"stimulation_pulses_{exp.session}.csv")
    pulses = pulses[~pulses.glitch]
    tasks = []
    for tt in a.tetrodes:
        g = int(tt.removeprefix("TT")) - 1
        wires = [str(c) for c in ch_ids[cp.tetrode_channel_indices(rec, g)]]
        T, ids = templates_for(tt, dirs, wires)
        for blk in blocks.itertuples(index=False):
            tasks.append((g, blk, T, ids))
    print(f"{len(tasks)} tasks ({len(blocks)} blocks x {len(a.tetrodes)} tetrodes)", flush=True)
    rows = Parallel(n_jobs=a.workers)(delayed(one)(binary, g, blk, pulses, T, ids, a.margin_s, a.thresh_mad, a.gate_r)
                                      for g, blk, T, ids in tasks)
    df = pd.DataFrame([r for rr in rows for r in rr])
    df.round(4).to_csv(RESULTS / f"tagged_xsession_{exp.experiment}.csv", index=False)
    summ = {}
    for role, (tt, rid) in TAGGED.items():
        s = df[(df.tetrode == tt) & (df.opto2_unit == rid)]
        if not len(s):
            continue
        summ[role] = {"opto2_unit": rid, "blocks": int(len(s)), "median_base_hz": round(float(s.base_hz.median()), 2),
                      "median_onset_hz": round(float(s.onset_1_10ms_hz.median()), 2), "median_in_pulse_hz": round(float(s.in_pulse_hz.median()), 2),
                      "median_post_0_50_hz": round(float(s.post_0_50ms_hz.median()), 2), "median_reliability": round(float(s.reliability_1_10ms.median()), 3),
                      "blocks_in_pulse_over_2x_base": int((s.in_pulse_hz > 2 * s.base_hz).sum()),
                      "blocks_onset_over_3x_base": int((s.onset_1_10ms_hz > 3 * s.base_hz).sum()),
                      "median_rate_window_hz": round(float(s.rate_window_hz.median()), 2)}
    (RESULTS / f"tagged_xsession_{exp.experiment}.json").write_text(json.dumps(summ, indent=1) + "\n")
    print(json.dumps(summ, indent=1))
    fig, axes = plt.subplots(len(a.tetrodes), 1, figsize=(16, 3.2 * len(a.tetrodes)), sharex=True)
    for ax, tt in zip(np.atleast_1d(axes), a.tetrodes):
        s = df[df.tetrode == tt]
        for u, su in s.groupby("opto2_unit"):
            su = su.sort_values("onset_h")
            role = [k for k, (t_, r_) in TAGGED.items() if t_ == tt and r_ == u]
            ax.plot(su.onset_h, su.in_pulse_hz, "s-" if role else ".-", lw=1.5 if role else 0.6, label=f"u{u}{' ' + role[0] if role else ''}")
        ax.set_ylabel("in-pulse spikes/s"); ax.set_title(f"{tt}: in-pulse rate per OPTO2 template at each block", fontsize=9); ax.legend(fontsize=6, ncol=6)
    np.atleast_1d(axes)[-1].set_xlabel("block onset (session hour)")
    fig.tight_layout(); fig.savefig(exp.figures_dir / "tagged_xsession.png", dpi=110); plt.close(fig)
    print(f"wrote results/tagged_xsession_{exp.experiment}.csv and figures/tagged_xsession.png")


if __name__ == "__main__":
    main()
