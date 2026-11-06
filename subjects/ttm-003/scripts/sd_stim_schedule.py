"""Pulse and sham schedule for stimulating only the first and last hour of a sleep deprivation.

Each window gets ``E = R * k / (k - 1)`` events, R light pulses and R / (k - 1) shams, with
``k = --sham-every``. Event times are E sorted uniform draws on the window's free time, the j-th
(from 0) shifted by ``j * --min-interval-s``: every interval is at least the minimum, the excess over it
is close to exponential, so the moment of the next event cannot be anticipated, and the count is
exact. A uniformly drawn set of exactly E / k events are shams, so each event, pulse or sham, follows a
sham with the same probability and the time since the previous light pulse has one distribution for
both. With ``--independent-windows`` each window gets its own draw; by default both
windows use the same schedule relative to their start.

Writes ``<out>.csv`` (one row per event: window, index, type, offset from the window start in ms,
clock time, interval from the previous event in ms) and ``<out>.json`` (the parameters, seed and
checks). Times are rounded to 1 ms after the draw; the checks are on the rounded times.

Usage: python subjects/ttm-003/scripts/sd_stim_schedule.py --pulses-per-hour 480 --sham-every 2 --seed 1
           --date 2026-12-01 --windows 09:00-10:00 13:00-14:00 --out sched_480
"""

from __future__ import annotations

import argparse
import datetime as dt
import json

import numpy as np
import pandas as pd


def draw(n_events: int, span_ms: int, margin_ms: int, min_ms: int, k: int, rng) -> pd.DataFrame:
    free = span_ms - 2 * margin_ms - (n_events - 1) * min_ms
    if free <= 0:
        raise SystemExit(f"{n_events} events at {min_ms} ms minimum do not fit in {span_ms} ms")
    u = np.sort(rng.uniform(0, free, n_events))
    t = np.round(margin_ms + u + np.arange(n_events) * min_ms).astype(np.int64)
    kind = np.full(n_events, "pulse", dtype=object)
    kind[rng.choice(n_events, n_events // k, replace=False)] = "sham"
    return pd.DataFrame({"offset_ms": t, "type": kind})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pulses-per-hour", type=int, required=True, help="light pulses per hour of window")
    ap.add_argument("--sham-every", type=int, default=2, help="one event in this many is a sham")
    ap.add_argument("--min-interval-s", type=float, default=1.0, help="shortest time between events, onset to onset")
    ap.add_argument("--margin-s", type=float, default=5.0, help="no event this close to a window's start or end")
    ap.add_argument("--windows", nargs="+", default=["09:00-10:00", "13:00-14:00"], help="HH:MM-HH:MM, local")
    ap.add_argument("--date", required=True, help="YYYY-MM-DD of the sleep deprivation")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--independent-windows", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    k = a.sham_every
    rng = np.random.default_rng(a.seed)
    rows, checks, base = [], [], None
    for w, spec in enumerate(a.windows):
        s, e = (dt.datetime.fromisoformat(f"{a.date}T{x}") for x in spec.split("-"))
        span_ms = round((e - s).total_seconds() * 1000)
        n_p = round(a.pulses_per_hour * span_ms / 3_600_000)
        if n_p % (k - 1):
            raise SystemExit(f"{n_p} pulses per window is not a multiple of {k - 1}")
        n_e = n_p * k // (k - 1)
        if base is None or a.independent_windows:
            d = draw(n_e, span_ms, round(a.margin_s * 1000), round(a.min_interval_s * 1000), k, rng)
            base = d
        d = base.copy()
        d.insert(0, "window", w + 1)
        d.insert(1, "index", np.arange(len(d)))
        d["clock"] = [(s + dt.timedelta(milliseconds=int(x))).isoformat(timespec="milliseconds") for x in d.offset_ms]
        d["interval_ms"] = d.offset_ms.diff().astype("Int64")
        iv = d.interval_ms.dropna()
        light = d.offset_ms.where(d.type == "pulse").ffill().shift()
        since = (d.offset_ms - light) / 1000
        checks.append({"window": w + 1, "start": s.isoformat(), "end": e.isoformat(), "pulses": int((d.type == "pulse").sum()),
                       "shams": int((d.type == "sham").sum()), "first_ms": int(d.offset_ms.iloc[0]),
                       "last_ms": int(d.offset_ms.iloc[-1]), "interval_ms_min": int(iv.min()),
                       "interval_ms_mean": float(iv.mean()), "interval_ms_max": int(iv.max()),
                       "since_pulse_s_median_pulse": float(since[d.type == "pulse"].median()),
                       "since_pulse_s_median_sham": float(since[d.type == "sham"].median()),
                       "min_interval_ok": bool(iv.min() >= round(a.min_interval_s * 1000)),
                       "margins_ok": bool(d.offset_ms.iloc[0] >= a.margin_s * 1000
                                          and d.offset_ms.iloc[-1] <= span_ms - a.margin_s * 1000)})
        rows.append(d)
    out = pd.concat(rows)
    out.to_csv(f"{a.out}.csv", index=False)
    with open(f"{a.out}.json", "w") as f:
        json.dump({"params": vars(a), "checks": checks}, f, indent=1)
    print(pd.DataFrame(checks).to_string(index=False))


if __name__ == "__main__":
    main()
