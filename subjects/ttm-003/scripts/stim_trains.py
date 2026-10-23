"""Join the laser TTL record with the hand stimulation log.

Reads line 1 of the acquisition-board TTL events of every session named in
``subjects/ttm-003/stimulation_log.toml``, groups pulses into trains (gaps
> 2 s), drops the sub-millisecond toggles that precede each train, and
joins the trains in order with the hand log. Writes
``subjects/ttm-003/results/stimulation_trains.csv`` (one row per train) and
``results/stimulation_pulses_<session>.csv`` (one row per pulse, session
seconds relative to the first acquisition-board sample).

A session listed without ``trains`` has no hand log; its trains are written
unjoined to ``results/stimulation_blocks_<session>.csv`` with the onset in
session hours and local time, the pattern matched from the TTL alone,
``anomalous`` set where a pulse width lies outside every pattern, and, when
the entry carries an ``expected`` protocol, its power and ``as_expected`` set
where the train's pattern and pulse count match it. ``delivered`` is False
for trains inside the entry's ``not_delivered`` windows (session hours).

Usage: python subjects/ttm-003/scripts/stim_trains.py
"""

from __future__ import annotations

import datetime as dt
import pathlib
import tomllib
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd

from ttm.experiment import read_session_t0_unix
from ttm.subject import find_subject

SUBJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
LOG = SUBJECT_DIR / "stimulation_log.toml"
RESULTS = SUBJECT_DIR / "results"
LASER_LINE = 1
TRAIN_GAP_S = 2.0
MIN_PULSE_MS = 1.0
STREAM = "Acquisition_Board-101.acquisition_board"


def read_laser_pulses(acquisition_dir: pathlib.Path) -> pd.DataFrame:
    rec = next(acquisition_dir.glob("Record Node */experiment1/recording1"))
    ttl = rec / "events" / STREAM / "TTL"
    states = np.load(ttl / "states.npy")
    samples = np.load(ttl / "sample_numbers.npy")
    first = np.load(rec / "continuous" / STREAM / "sample_numbers.npy", mmap_mode="r")[0]
    keep = np.abs(states) == LASER_LINE
    t = (samples[keep] - first) / 30000.0
    st = states[keep]
    on, off = t[st > 0], t[st < 0]
    if len(on) != len(off) or np.any(off < on):
        raise ValueError(f"{acquisition_dir}: unpaired laser TTL edges")
    df = pd.DataFrame({"onset_s": on, "offset_s": off})
    df["duration_ms"] = (df.offset_s - df.onset_s) * 1000
    df["glitch"] = df.duration_ms < MIN_PULSE_MS
    real = df[~df.glitch]
    train = np.r_[0, np.cumsum(np.diff(real.onset_s.values) > TRAIN_GAP_S)]
    df["train"] = -1
    df.loc[~df.glitch, "train"] = train + 1
    df["pulse_in_train"] = -1
    df.loc[~df.glitch, "pulse_in_train"] = df[~df.glitch].groupby("train").cumcount() + 1
    return df


def summarize_trains(pulses: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for k, g in pulses[~pulses.glitch].groupby("train"):
        ipi = np.diff(g.onset_s.values)
        rows.append({
            "ttl_train": k,
            "onset_s": g.onset_s.iloc[0],
            "offset_s": g.offset_s.iloc[-1],
            "duration_s": round(float(g.offset_s.iloc[-1] - g.onset_s.iloc[0]), 2),
            "n_pulses": len(g),
            "pulse_ms_median": round(float(g.duration_ms.median()), 2),
            "pulse_ms_min": round(float(g.duration_ms.min()), 2),
            "pulse_ms_max": round(float(g.duration_ms.max()), 2),
            "ipi_ms_median": round(float(np.median(ipi)) * 1000, 2),
            "rate_hz": round(1 / float(np.median(ipi)), 3),
            "n_glitch_edges_before": int(((pulses.glitch) & (pulses.onset_s < g.onset_s.iloc[0])
                                          & (pulses.onset_s > g.onset_s.iloc[0] - TRAIN_GAP_S)).sum()),
        })
    return pd.DataFrame(rows)


def ttl_pattern(row: pd.Series, patterns: dict) -> str:
    for name, p in patterns.items():
        if abs(row.pulse_ms_median - p["pulse_ms"]) < 25 and abs(row.rate_hz - p["rate_hz"]) < 0.3:
            return name
    return "?"


def unlogged_blocks(trains: pd.DataFrame, patterns: dict, t0_unix: float, tz: str,
                    expected: dict | None = None, not_delivered: list[dict] = ()) -> pd.DataFrame:
    """Per-train table for a session with no hand log."""
    blocks = trains.copy()
    blocks["onset_h"] = (blocks.onset_s / 3600).round(4)
    t0 = dt.datetime.fromtimestamp(t0_unix, ZoneInfo(tz))
    blocks["onset_local"] = [(t0 + dt.timedelta(seconds=float(s))).strftime("%Y-%m-%d %H:%M:%S") for s in blocks.onset_s]
    blocks["ttl_pattern"] = blocks.apply(ttl_pattern, axis=1, patterns=patterns)
    widths = [p["pulse_ms"] for p in patterns.values()]
    blocks["anomalous"] = [not any(abs(w - pw) < 25 for pw in widths) for w in blocks.pulse_ms_max]
    if expected:
        blocks["expected_pattern"] = expected["pattern"]
        blocks["laser_level"] = expected["laser_level"]
        blocks["power_mW"] = expected["power_mW"]
        blocks["as_expected"] = (blocks.ttl_pattern == expected["pattern"]) & (blocks.n_pulses == expected["pulses_per_train"])
    blocks["delivered"] = True
    for w in not_delivered:
        blocks.loc[(blocks.onset_h >= w["from_h"]) & (blocks.onset_h < w["to_h"]), "delivered"] = False
    return blocks


def main() -> None:
    log = tomllib.loads(LOG.read_text())
    subject = find_subject("TTM-003")
    RESULTS.mkdir(exist_ok=True)
    out = []
    for session, entry in log["sessions"].items():
        exp = subject.experiment(entry["experiment"])
        assert exp.session == session, (exp.session, session)
        pulses = read_laser_pulses(exp.acquisition_dir)
        pulses.round(5).to_csv(RESULTS / f"stimulation_pulses_{session}.csv", index=False)
        trains = summarize_trains(pulses)
        if "trains" not in entry:
            blocks = unlogged_blocks(trains, log["patterns"], read_session_t0_unix(exp.acquisition_dir), exp.timezone,
                                     expected=entry.get("expected"), not_delivered=entry.get("not_delivered", ()))
            blocks.insert(0, "experiment", entry["experiment"])
            blocks.insert(0, "session", session)
            blocks.round(5).to_csv(RESULTS / f"stimulation_blocks_{session}.csv", index=False)
            show = ["ttl_train", "onset_h", "onset_local", "n_pulses", "pulse_ms_median", "pulse_ms_max",
                    "rate_hz", "duration_s", "ttl_pattern", "anomalous", "delivered"] + (["as_expected"] if "as_expected" in blocks else [])
            print(blocks[show].to_string(index=False))
            continue
        hand = pd.DataFrame(entry["trains"]).rename(columns={"pattern": "logged_pattern"})
        if len(hand) != len(trains):
            raise ValueError(f"{session}: {len(hand)} logged trains, {len(trains)} in the TTL")
        merged = pd.concat([hand, trains], axis=1)
        merged.insert(0, "experiment", entry["experiment"])
        merged.insert(0, "session", session)
        merged["ttl_pattern"] = merged.apply(ttl_pattern, axis=1, patterns=log["patterns"])
        merged["onset_min"] = (merged.onset_s / 60).round(2)
        out.append(merged)
    df = pd.concat(out, ignore_index=True)
    df.to_csv(RESULTS / "stimulation_trains.csv", index=False)
    cols = ["session", "n", "logged_pattern", "ttl_pattern", "laser_level", "power_mW", "approx_time",
            "onset_min", "n_pulses", "pulse_ms_median", "rate_hz", "state", "neural_response"]
    print(df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
