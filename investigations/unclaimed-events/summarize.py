"""Tables of what the fold and multi-unit steps would do, from replay.py's per-window files.

    .venv/bin/python investigations/unclaimed-events/summarize.py <replay out dir>

results/tetrodes.csv, per tetrode: events of 7 and of 10 MAD or more and the share no unit claimed; the
fold's proposals, refractory rejections and folded spikes; the multi-unit events and their rate over the
tetrode's time outside bad intervals, next to the sorted spikes' rate; how the unclaimed events of 7 MAD
or more split by best cosine (under 0.55, 0.55 to 0.8, 0.8 or more) next to the units' own template
spikes of 7 MAD or more; the share of claimed, folded and multi-unit events with an event on 3 or more
other tetrodes within 0.5 ms; the share inside laser pulses (onset to offset plus 100 ms).

results/units.csv, per unit: sorted spikes; folded spikes and their share; the share of its template
spikes (up to 500 per window) whose own cosine is 0.8 or more and whose best template is its own (the
fold's rule applied to spikes the sort gave it), and whose best cosine is 0.55 to 0.8 (the multi-unit
band); the median scale of folded events and of its own spikes relative to its template, and the 99th
percentile of the latter; the share of folded events where the matcher fitted this unit within 0.5 ms
and why that fit was rejected; folded pairs closer than 1.5 ms; and, with 200 or more folded events, the
cross-correlogram of folded events with the unit's sorted spikes next to the unit's autocorrelogram and
the largest cross-correlogram of the folded events with a same-tetrode neighbour's spikes (mean count
per lag over 0.5 to 1.5 ms, 1.5 to 5 ms and 5 to 15 ms, over the mean per lag 20 to 50 ms away); the
multi-unit events whose best template is the unit's and, with 200 or more, their cross-correlogram with
its sorted spikes in the same terms.

results/own_cosine_by_amplitude.csv: the share of the units' own template spikes on TT5 to TT12 with
cosine 0.8 or more to their own template, by amplitude band.

results/same_tetrode_pairs.csv: every pair of units on one tetrode of TT5 to TT12, the cross-correlogram
in the terms of units.csv.

results/summary.json: totals and quantiles over all tetrodes (fold counts, also with the best of shifts
up to 3 samples; scale of folded and multi-unit events; the folded events' margin between best and
second-best cosine; multi-unit amplitude; the share of time inside laser pulses), and per TT5 to TT12
unit with 20 or more folded events, the share of its own spikes at least as large as its median folded
event.

results/fit_reasons.csv: for folded and multi-unit events, the share whose highest-snr matcher fit on the
tetrode within 0.5 ms, and whose highest-snr fit of the best-matching unit there, was rejected for each
reason (the floor split into snr under 4.5 and amplitude under the relative floor), or that had no such
fit.
"""

import itertools
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from ttm.carry import read_bad_periods

from replay import BAD, FS, PULSE_POST, PULSES, REFR, RUN, WIN

OUT = pathlib.Path(__file__).parent / "results"
REASONS = ("no fit", "snr floor", "relative floor", "shape", "refractory", "bad period")
SNR_FLOOR = 4.5
LAG_BINS = {"0p5_1p5ms": (16, 45), "1p5_5ms": (46, 150), "5_15ms": (151, 450)}
BASE = (600, 1500)
MIN_FOLDED = 200
LIVE = range(4, 12)
AMP_BANDS = [0, 5.5, 7, 10, 15, np.inf]


def load(d):
    ev, host, counts, wins = [], [], None, []
    for f in sorted(pathlib.Path(d).glob("w*.npz")):
        z = np.load(f)
        ev.append(pd.DataFrame({k[3:]: z[k] for k in z.files if k.startswith("ev_")}))
        host.append(pd.DataFrame({k[5:]: z[k] for k in z.files if k.startswith("host_")}))
        counts = z["counts"] if counts is None else counts + z["counts"]
        wins.append(int(f.stem[1:]))
    return pd.concat(ev, ignore_index=True), pd.concat(host, ignore_index=True), counts, wins


def lag_counts(ref, q, lo, hi):
    """Number of (q, ref) pairs with lo <= |ref - q| <= hi, for sorted ``ref``."""
    n = 0
    for s in (1, -1):
        a = np.searchsorted(ref, q + s * lo, "left" if s > 0 else "right")
        b = np.searchsorted(ref, q + s * hi, "right" if s > 0 else "left")
        n += np.abs(b - a).sum()
    return n


def profile(ref, q):
    """Mean pair count per lag in each LAG_BINS range over the mean per lag in BASE."""
    base = lag_counts(ref, q, *BASE) / (2 * (BASE[1] - BASE[0] + 1))
    return {k: lag_counts(ref, q, lo, hi) / (2 * (hi - lo + 1)) / base for k, (lo, hi) in LAG_BINS.items()}


def cos_split(c):
    return {"lt_0p55": float((c < 0.55).mean()), "0p55_0p8": float(((c >= 0.55) & (c < 0.8)).mean()),
            "ge_0p8": float((c >= 0.8).mean())} if len(c) else {}


def reason(code, snr):
    """Rejection reason per fit (carry.REASON_*), with the floor split by whether snr cleared SNR_FLOOR."""
    out = np.full(len(code), "accepted", object)
    for c, name in ((-1, "no fit"), (2, "shape"), (3, "refractory"), (5, "bad period")):
        out[code == c] = name
    out[(code == 1) & ~(snr >= SNR_FLOOR)] = "snr floor"
    out[(code == 1) & (snr >= SNR_FLOOR)] = "relative floor"
    return out


def share(m):
    return float(np.mean(m)) if len(m) else np.nan


def main():
    d = sys.argv[1]
    ev, host, counts, wins = load(d)
    s = si.load(f"{RUN}/sorting")
    nfr = si.load(f"{RUN}/binary").get_num_frames()
    ug = dict(zip(np.asarray(s.unit_ids).astype(int).tolist(), np.asarray(s.get_property("group")).astype(int).tolist()))
    spans = [(w * WIN, min((w + 1) * WIN, nfr)) for w in wins]
    trains = {}
    for u in ug:
        t = np.sort(s.get_unit_spike_train(u)).astype(np.int64)
        trains[u] = np.concatenate([t[np.searchsorted(t, a):np.searchsorted(t, b)] for a, b in spans])
    bad = read_bad_periods(BAD, FS)

    def observed_s(g):
        t = sum(b - a for a, b in spans)
        for a, b in bad.get(g, []):
            for lo, hi in spans:
                t -= max(0, min(b, hi) - max(a, lo))
        return t / FS

    rows = []
    for g in sorted(set(ug.values())):
        c = counts[g]
        e = ev[ev["group"] == g]
        h = host[host["uid"].map(ug) == g]
        n7, n10 = c[:, 1:].sum(), c[:, 2].sum()
        folded, mua = e["folded"].to_numpy(), e["mua"].to_numpy()
        row = {"tetrode": g + 1, "units": sum(v == g for v in ug.values()),
               "events_7mad": int(n7), "unclaimed_7mad": float(c[0, 1:].sum() / n7) if n7 else np.nan,
               "events_10mad": int(n10), "unclaimed_10mad": float(c[0, 2].sum() / n10) if n10 else np.nan,
               "fold_proposed": int(e["propose"].sum()), "fold_refractory": int((e["propose"] & ~e["refr_ok"]).sum()),
               "folded": int(folded.sum()), "multi_unit": int(mua.sum()),
               "no_template": int((e["best_uid"] < 0).sum()),
               "multi_unit_hz": float(mua.sum() / observed_s(g)),
               "sorted_hz": float(sum(len(trains[u]) for u in ug if ug[u] == g) / observed_s(g))}
        row.update({f"unclaimed_{k}": v for k, v in cos_split(e.loc[e["best_uid"] >= 0, "best_cos"]).items()})
        row.update({f"own7_{k}": v for k, v in cos_split(h.loc[h["amp_mad"] >= 7, "best_cos"]).items()})
        row.update({"claimed_3_other_tetrodes": float(c[1, 1:, 3:].sum() / c[1, 1:].sum()) if c[1, 1:].sum() else np.nan,
                    "folded_3_other_tetrodes": share(e.loc[folded, "others"] >= 3),
                    "multi_unit_3_other_tetrodes": share(e.loc[mua, "others"] >= 3),
                    "folded_in_pulse": share(e.loc[folded, "in_pulse"]),
                    "multi_unit_in_pulse": share(e.loc[mua, "in_pulse"]),
                    "multi_unit_amp_mad_median": float(e.loc[mua, "amp_mad"].median()) if mua.any() else np.nan})
        rows.append(row)
    tet = pd.DataFrame(rows)
    tet.to_csv(OUT / "tetrodes.csv", index=False)

    urows = []
    for u, g in ug.items():
        f = ev[ev["folded"] & (ev["best_uid"] == u)]
        h = host[host["uid"] == u]
        fr = np.sort(f["frame"].to_numpy())
        row = {"unit": u, "tetrode": g + 1, "spikes": len(trains[u]), "folded": len(fr),
               "folded_share": len(fr) / len(trains[u]) if len(trains[u]) else np.nan,
               "own_rule_share": share((h["own_cos"] >= 0.8) & (h["best_uid"] == u)),
               "own_multi_unit_band": share((h["best_cos"] >= 0.55) & (h["best_cos"] < 0.8)),
               "folded_alpha_median": float(f["alpha"].median()) if len(f) else np.nan,
               "own_alpha_median": float(h["alpha"].median()) if len(h) else np.nan,
               "own_alpha_p99": float(h["alpha"].quantile(0.99)) if len(h) else np.nan,
               "folded_amp_mad_median": float(f["amp_mad"].median()) if len(f) else np.nan,
               "own_amp_mad_median": float(h["amp_mad"].median()) if len(h) else np.nan,
               "folded_pairs_within_1p5ms": int((np.diff(fr) <= REFR).sum()) if len(fr) > 1 else 0}
        r = reason(f["hfit_reason"].to_numpy(), f["hfit_snr"].to_numpy())
        for name in REASONS:
            row[f"own_fit_{name.replace(' ', '_')}"] = share(r == name) if len(f) else np.nan
        mu = np.sort(ev.loc[ev["mua"] & (ev["best_uid"] == u), "frame"].to_numpy())
        row["multi_unit_best"] = len(mu)
        if len(mu) >= MIN_FOLDED:
            row.update({f"multi_unit_{k}": v for k, v in profile(trains[u], mu).items()})
        if len(fr) >= MIN_FOLDED or len(mu) >= MIN_FOLDED:
            row.update({f"own_{k}": v for k, v in profile(trains[u], trains[u]).items()})
        if len(fr) >= MIN_FOLDED:
            row.update({f"folded_{k}": v for k, v in profile(trains[u], fr).items()})
            nb = [profile(trains[v], fr) for v in ug if ug[v] == g and v != u and len(trains[v])]
            for k in LAG_BINS:
                row[f"folded_neighbour_{k}_max"] = max((p[k] for p in nb), default=np.nan)
        urows.append(row)
    units = pd.DataFrame(urows).sort_values(["tetrode", "unit"])
    units.to_csv(OUT / "units.csv", index=False)

    frows = []
    for name, m in (("folded", ev["folded"]), ("multi-unit", ev["mua"])):
        e = ev[m]
        for prefix, whose in (("fit", "strongest fit"), ("hfit", "best-matching unit's fit")):
            r = reason(e[f"{prefix}_reason"].to_numpy(), e[f"{prefix}_snr"].to_numpy())
            for rn in REASONS:
                frows.append({"events": name, "fit": whose, "reason": rn, "share": share(r == rn)})
    fits = pd.DataFrame(frows)
    fits.to_csv(OUT / "fit_reasons.csv", index=False)

    hl = host[host["uid"].map(ug).isin(LIVE)]
    band = pd.cut(hl["amp_mad"], AMP_BANDS, right=False)
    amp = (hl["own_cos"] >= 0.8).groupby(band, observed=True).agg(["mean", "size"])
    amp = amp.rename(columns={"mean": "own_cos_0p8_share", "size": "spikes"}).rename_axis("amp_mad").reset_index()
    amp["amp_mad"] = amp["amp_mad"].astype(str)
    amp.to_csv(OUT / "own_cosine_by_amplitude.csv", index=False)

    prow = []
    for g in LIVE:
        for a, b in itertools.combinations([u for u in ug if ug[u] == g], 2):
            prow.append({"unit_a": a, "unit_b": b, "tetrode": g + 1, **profile(trains[a], trains[b])})
    pairs = pd.DataFrame(prow)
    pairs.to_csv(OUT / "same_tetrode_pairs.csv", index=False)

    f = ev[ev["folded"]]
    m = ev[ev["mua"]]
    pr = pd.read_csv(PULSES)
    fl = units[units["tetrode"].isin([g + 1 for g in LIVE]) & (units["folded"] >= 20)]
    at_least = {int(r.unit): share(host.loc[host["uid"] == r.unit, "alpha"] >= r.folded_alpha_median)
                for r in fl.itertuples()}
    q = lambda x, ps: {str(p): float(x.quantile(p)) for p in ps}
    summary = {
        "fold_proposed": int(ev["propose"].sum()),
        "fold_refractory": int((ev["propose"] & ~ev["refr_ok"]).sum()),
        "folded": len(f), "folded_tt5_12": int(f["group"].isin(LIVE).sum()),
        "fold_proposed_best_shift": int(((ev["amp_mad"] >= 10) & (ev["best_shift_cos"] >= 0.8) & (ev["best_uid"] >= 0)).sum()),
        "folded_other_unit_best_shift": share(f["best_uid"] != f["best_shift_uid"]),
        "folded_alpha": q(f["alpha"], (0.1, 0.5, 0.9)),
        "folded_cosine_margin_median": float((f["best_cos"] - f["second_cos"]).median()),
        "multi_unit": len(m), "multi_unit_tt5_12": int(m["group"].isin(LIVE).sum()),
        "multi_unit_over_sorted": len(m) / sum(len(t) for t in trains.values()),
        "multi_unit_amp_mad": q(m["amp_mad"], (0.5, 0.9)),
        "multi_unit_alpha": q(m["alpha"], (0.1, 0.5, 0.9)),
        "pulse_time_share": float((pr["offset_s"] + PULSE_POST / FS - pr["onset_s"]).sum() * FS / sum(b - a for a, b in spans)),
        "own_at_least_folded_median": at_least,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 60)
    print(f"windows {len(wins)}")
    print(tet.round(3).T.to_string())
    print(units.round(3).to_string(index=False))
    print(fits.round(3).to_string(index=False))
    print(amp.round(3).to_string(index=False))
    print(pairs.describe().round(2).to_string())
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
