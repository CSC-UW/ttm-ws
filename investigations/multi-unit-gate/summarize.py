"""Per-tetrode cosine cutoffs for the multi-unit gate on the final TTM-003-48H sorting, from gate.py's
per-window files.

    .venv/bin/python investigations/multi-unit-gate/summarize.py <gate out dir> [--sub]

With --sub, every score is the one after subtracting the neighbouring sorted spikes (gate.py), and each
file below gets the suffix _sub. Kept counts are then also given for unclaimed events whose residual
amplitude is 7 or more (kept_res7_*, hz_res7_*).

The cutoff losing x % on a tetrode is the x-th percentile of its reference events' best cosine, so that
a share x of the events exactly one sorted spike claims scores below it. The amplitude-matched cutoff
weights the reference events so that their amplitudes (0.25 bins from 7 to 15, then 15 or more) are
distributed like the unclaimed events'. Kept counts and rates drop an event within 0.3 ms of the previous
kept one of the tetrode, as 97_production_pipeline.py does. Writes to results/:

tetrodes.csv: per tetrode with templates, the events of 7 or more (excluded at window edges or in bad
    intervals, claimed by two or more sorted spikes, reference, unclaimed, unclaimed with no template),
    the reference share below 0.55, the cutoffs losing 1, 5 and 10 % (plain and amplitude-matched, and
    the range of the 5 % cutoff over eight blocks of 12 windows), the share and rate of unclaimed events
    each cutoff keeps, the share of background stretches each lets through, the share of events with an
    event on 3 or more other tetrodes, the sorted rate, and unclaimed_7mad_replay and band_replay from
    investigations/unclaimed-events/results/tetrodes.csv next to the same counts here.
by_amplitude.csv: per tetrode and amplitude band, the reference cutoffs losing 1, 5 and 10 % and the
    share of unclaimed events kept at the tetrode's 5 % cutoff and at 0.55.
banded.csv: per tetrode and share lost (5 and 10 %), the unclaimed events kept by the tetrode's one
    cutoff and by a cutoff per amplitude band (the bands of by_amplitude.csv; the tetrode's cutoff in a band
    with fewer than 100 reference events), those kept by only one of the two, and the share of reference
    events of 7 to 8 and of 12 or more that the one cutoff loses.
cosine_hist.csv: per tetrode, counts of best cosine in 0.02 bins for reference, unclaimed and background.
summary.json: totals over TT5 to TT12.
"""

import json
import pathlib
import sys

import numpy as np
import pandas as pd
from gate import AFTER, BEFORE, FLOOR, FS, Setup

OUT = pathlib.Path(__file__).parent / "results"
REPLAY = pathlib.Path(__file__).parents[1] / "unclaimed-events/results/tetrodes.csv"
LOSS = (1, 5, 10)
AMP_EDGES = np.r_[np.arange(7.0, 15.0, 0.25), 15.0, np.inf]
AMP_BANDS = ((7, 8), (8, 9), (9, 10), (10, 12), (12, 15), (15, np.inf))
HIST = np.linspace(-1, 1, 101)
COINC = 9
N_BLOCKS = 8
LIVE = range(4, 12)
SUB = "--sub" in sys.argv
SFX = "_sub" if SUB else ""


def wpct(v, w, q):
    """Weighted q-th percentile (0-100) of ``v``: the smallest value whose cumulative weight reaches q %."""
    o = np.argsort(v)
    cw = np.cumsum(w[o])
    return float(v[o][np.searchsorted(cw, q / 100 * cw[-1])])


def kept(frames):
    """Number of frames left after dropping each within COINC of the previous kept one."""
    f = np.sort(frames)
    n, last = 0, -np.inf
    for v in f:
        if v - last > COINC:
            n, last = n + 1, v
    return n


def kept_fast(frames):
    """kept() without the loop when no two frames are within COINC of each other; else kept()."""
    f = np.sort(frames)
    return f.size if f.size < 2 or np.diff(f).min() > COINC else kept(f)


def main():
    d = pathlib.Path(sys.argv[1])
    st = Setup()
    files = sorted(d.glob("w*.npz"))
    ev, bg, counts, spans = {}, {"group": [], "cos": []}, {}, []
    for wi, f in enumerate(files):
        z = np.load(f)
        spans.append((int(z["a"]) + BEFORE, int(z["b"]) - AFTER))
        n = z["ev_frame"].size
        for k in z.files:
            if k.startswith("ev_"):
                ev.setdefault(k[3:], []).append(z[k])
        ev.setdefault("window", []).append(np.full(n, wi, np.int16))
        bg["group"].append(z["bg_group"])
        bg["cos"].append(z["bg_sub_cos"] if SUB else z["bg_cos"])
        for k in ("n_events", "n_excluded", "n_overlap"):
            counts[k] = counts.get(k, 0) + z[k]
    ev = {k: np.concatenate(v) for k, v in ev.items()}
    bg = {k: np.concatenate(v) for k, v in bg.items()}
    assert len(files) == 96
    if SUB:
        ev["best_cos"], ev["best_uid"], ev["alpha"] = ev["sub_cos"], ev["sub_uid"], ev["sub_alpha"]
    res7 = ev["res_amp"] >= FLOOR if "res_amp" in ev else np.ones(ev["frame"].size, bool)

    def observed_s(g):
        tot = sum(b - a for a, b in spans)
        if g in st.bad:
            for s0, e0 in zip(*st.bad[g]):
                tot -= sum(max(0, min(e0, b) - max(s0, a)) for a, b in spans)
        return tot / FS

    def sorted_spikes(g):
        if g not in st.pool_f:
            return 0
        p = st.pool_f[g]
        inside = np.zeros(p.size, bool)
        for a, b in spans:
            inside[np.searchsorted(p, a):np.searchsorted(p, b)] = True
        return int((inside & ~st.in_bad(p, g)).sum())

    replay = pd.read_csv(REPLAY).set_index("tetrode")
    rows, amp_rows, hist_rows, band_rows = [], [], [], []
    for g in sorted(np.unique(ev["group"][ev["best_uid"] >= 0]).tolist()):
        on = ev["group"] == g
        ref = on & (ev["claims"] == 1) & (ev["best_uid"] >= 0)
        un_all = on & (ev["claims"] == 0)
        un = un_all & (ev["best_uid"] >= 0)
        rc, uc = ev["best_cos"][ref], ev["best_cos"][un]
        bgc = bg["cos"][bg["group"] == g]
        ha, _ = np.histogram(ev["amp"][un], AMP_EDGES)
        hr, _ = np.histogram(ev["amp"][ref], AMP_EDGES)
        ratio = np.where(hr > 0, ha / np.maximum(hr, 1), 0.0)
        w = ratio[np.clip(np.searchsorted(AMP_EDGES, ev["amp"][ref], "right") - 1, 0, len(ratio) - 1)]
        cuts = {q: float(np.percentile(rc, q)) for q in LOSS}
        obs = observed_s(g)
        row = {"tetrode": g + 1, "units": sum(1 for u in st.ugroup.values() if u == g),
               "events_7": int(counts["n_events"][g]), "excluded_7": int(counts["n_excluded"][g]),
               "claimed_2_or_more": int(counts["n_overlap"][g]), "reference": int(ref.sum()),
               "unclaimed": int(un_all.sum()), "unclaimed_no_template": int((un_all & ~un).sum()),
               "reference_below_0p55": float(np.mean(rc < 0.55)),
               "observed_h": obs / 3600, "sorted_hz": sorted_spikes(g) / obs}
        for q in LOSS:
            row[f"cut_{q}"] = cuts[q]
            row[f"cut_matched_{q}"] = wpct(rc, w, q)
        win = ev["window"][ref]
        block = [float(np.percentile(rc[(win // 12) == i], 5)) for i in range(N_BLOCKS)
                 if np.any((win // 12) == i)]
        row["cut_5_block_min"], row["cut_5_block_max"] = min(block), max(block)
        for name, c in [("0p55", 0.55)] + [(f"cut_{q}", cuts[q]) for q in LOSS]:
            k = un & (ev["best_cos"] >= c)
            row[f"kept_{name}"] = float(k.sum() / max(un.sum(), 1))
            row[f"hz_{name}"] = kept_fast(ev["frame"][k]) / obs
            row[f"background_{name}"] = float(np.mean(bgc >= c))
            if SUB:
                row[f"kept_res7_{name}"] = float((k & res7).sum() / max(un.sum(), 1))
                row[f"hz_res7_{name}"] = kept_fast(ev["frame"][k & res7]) / obs
        if "res_amp" in ev:
            row["reference_res_below_7"] = float(np.mean(~res7[ref]))
            row["unclaimed_res_below_7"] = float(np.mean(~res7[un]))
        o3 = ev["others"] >= 3
        row["others3_reference"] = float(np.mean(o3[ref]))
        row["others3_below_cut_5"] = float(np.mean(o3[un & (ev["best_cos"] < cuts[5])]))
        row["others3_above_cut_5"] = float(np.mean(o3[un & (ev["best_cos"] >= cuts[5])]))
        band = un & (ev["best_cos"] >= 0.55) & (ev["best_cos"] < 0.8)
        row["band_0p55_0p8"] = kept_fast(ev["frame"][band])
        if g + 1 in replay.index:
            r = replay.loc[g + 1]
            row["unclaimed_7mad_replay"] = round(float(r["events_7mad"] * r["unclaimed_7mad"]))
            row["band_replay"] = int(r["multi_unit"])
        rows.append(row)
        for lo, hi in AMP_BANDS:
            ra = ref & (ev["amp"] >= lo) & (ev["amp"] < hi)
            ua = un & (ev["amp"] >= lo) & (ev["amp"] < hi)
            r = {"tetrode": g + 1, "amp_lo": lo, "amp_hi": hi, "reference": int(ra.sum()), "unclaimed": int(ua.sum())}
            for q in LOSS:
                r[f"cut_{q}"] = float(np.percentile(ev["best_cos"][ra], q)) if ra.sum() >= 100 else np.nan
            r["kept_cut_5"] = float(np.mean(ev["best_cos"][ua] >= cuts[5])) if ua.any() else np.nan
            r["kept_0p55"] = float(np.mean(ev["best_cos"][ua] >= 0.55)) if ua.any() else np.nan
            amp_rows.append(r)
        amp_on, cos_on, ref_on, un_on = ev["amp"][on], ev["best_cos"][on], ref[on], un[on]
        for q in (5, 10):
            bc = np.full(amp_on.size, cuts[q])
            for lo, hi in AMP_BANDS:
                ib = (amp_on >= lo) & (amp_on < hi)
                if (ib & ref_on).sum() >= 100:
                    bc[ib] = np.percentile(cos_on[ib & ref_on], q)
            one, banded = un_on & (cos_on >= cuts[q]), un_on & (cos_on >= bc)
            band_rows.append({"tetrode": g + 1, "loss": q, "unclaimed": int(un_on.sum()), "kept_one": int(one.sum()),
                              "kept_banded": int(banded.sum()), "one_only": int((one & ~banded).sum()),
                              "banded_only": int((banded & ~one).sum()),
                              "reference_lost_7_8": float(np.mean(cos_on[ref_on & (amp_on < 8)] < cuts[q])),
                              "reference_lost_12_up": float(np.mean(cos_on[ref_on & (amp_on >= 12)] < cuts[q]))})
        for kind, v in (("reference", rc), ("unclaimed", uc), ("background", bgc)):
            h, _ = np.histogram(v, HIST)
            hist_rows += [{"tetrode": g + 1, "kind": kind, "bin_lo": round(float(lo), 2), "count": int(c)}
                          for lo, c in zip(HIST[:-1], h)]

    t = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    t.to_csv(OUT / f"tetrodes{SFX}.csv", index=False)
    pd.DataFrame(amp_rows).to_csv(OUT / f"by_amplitude{SFX}.csv", index=False)
    pd.DataFrame(hist_rows).to_csv(OUT / f"cosine_hist{SFX}.csv", index=False)
    pd.DataFrame(band_rows).to_csv(OUT / f"banded{SFX}.csv", index=False)
    live = t[t["tetrode"].isin([g + 1 for g in LIVE])]
    summary = {"windows": len(files),
               "tt5_12": {"unclaimed": int(live["unclaimed"].sum()), "reference": int(live["reference"].sum()),
                          "claimed_2_or_more": int(live["claimed_2_or_more"].sum()),
                          "unclaimed_7mad_replay": int(live["unclaimed_7mad_replay"].sum()),
                          "band_0p55_0p8": int(live["band_0p55_0p8"].sum()), "band_replay": int(live["band_replay"].sum()),
                          **{f"multi_unit_{n}": float((live[f"hz_{n}"] * live["observed_h"] * 3600).sum())
                             for n in ("0p55", "cut_1", "cut_5", "cut_10")},
                          **({f"multi_unit_res7_{n}": float((live[f"hz_res7_{n}"] * live["observed_h"] * 3600).sum())
                              for n in ("0p55", "cut_1", "cut_5", "cut_10")} if SUB else {}),
                          "sorted_spikes": float((live["sorted_hz"] * live["observed_h"] * 3600).sum())}}
    (OUT / f"summary{SFX}.json").write_text(json.dumps(summary, indent=1))
    pd.set_option("display.width", 250)
    print(t.round(3).T.to_string())
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
