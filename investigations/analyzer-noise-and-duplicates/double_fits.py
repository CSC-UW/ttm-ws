"""Whether the zero-lag coincidences between same-tetrode units of the final TTM-003-48H sorting are
one spike fitted twice or two spikes, from the recording itself.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/double_fits.py

For a pair (a, b), a coincident event is a spike of a with a spike of b within 3 samples. Its
waveform on the tetrode (1 ms before to 2 ms after a's spike, µV, the int16 binary) is compared
with three predictions made from the units' mean waveforms over spikes with no partner within
1.5 ms: a's alone, b's alone (shifted by the lag), and their sum. One spike fitted twice is
explained by one template; two spikes by the sum. Spikes of a without a partner are scored the same
way as a check. Also counts each unit's own spikes closer than 3 and 45 samples (0.1 and 1.5 ms).
"""

import pathlib

import numpy as np
import pandas as pd
import spikeinterface as si

RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
OUT = pathlib.Path(__file__).parent / "results"
NEAR, ISO = 3, 45
BEFORE, AFTER = 30, 60
N_EVENTS = 1000
LIVE = range(5, 13)


def snippets(rec, frames, channels):
    out = np.empty((len(frames), BEFORE + AFTER, len(channels)), np.float32)
    for k, f in enumerate(frames):
        out[k] = rec.get_traces(start_frame=int(f) - BEFORE, end_frame=int(f) + AFTER, channel_ids=channels,
                                return_in_uV=True)
    return out


def shift(t, lag):
    """Template ``t`` (samples, channels) moved later by ``lag`` samples, zero-filled."""
    s = np.zeros_like(t)
    if lag >= 0:
        s[lag:] = t[:len(t) - lag]
    else:
        s[:lag] = t[-lag:]
    return s


def partners(ta, tb, lo, hi):
    """For each frame of ``ta``, the offset to the nearest frame of ``tb`` within [lo, hi], or none."""
    i = np.clip(np.searchsorted(tb, ta), 1, len(tb) - 1)
    d = np.where(np.abs(tb[i] - ta) < np.abs(tb[i - 1] - ta), tb[i] - ta, tb[i - 1] - ta)
    return np.where((d >= lo) & (d <= hi), d, np.iinfo(np.int64).min)


def main():
    rng = np.random.default_rng(0)
    rec = si.load(f"{RUN}/binary")
    s = si.load(f"{RUN}/sorting")
    uids = list(s.unit_ids)
    group = dict(zip(uids, np.asarray(s.get_property("group"))))
    rec_group = np.asarray(rec.get_property("group"))
    n_frames = rec.get_num_frames()
    tr = {u: np.sort(s.get_unit_spike_train(u)).astype(np.int64) for u in uids}
    tr = {u: t[(t >= BEFORE) & (t < n_frames - AFTER)] for u, t in tr.items()}

    isi = pd.DataFrame([{"unit": u, "tetrode": int(group[u]) + 1, "spikes": len(t),
                         "isi_le_3": int((np.diff(t) <= NEAR).sum()), "isi_le_45": int((np.diff(t) <= ISO).sum())}
                        for u, t in tr.items()])
    isi.to_csv(OUT / "within_unit_isi.csv", index=False)
    print(f"within-unit intervals <= 3 samples: {isi.isi_le_3.sum()}, <= 45 samples: {isi.isi_le_45.sum()}")

    z = pd.read_csv(OUT / "zero_lag_pairs.csv")
    z = z[(z["pair"] == "same tetrode") & z["tetrode_a"].isin(LIVE)]
    test = z[z["excess_of_smaller_unit"] > 0.01]
    ctrl = z[z["excess_of_smaller_unit"] <= 0.01].assign(d=lambda x: (x.near_over_expected - 1).abs()).nsmallest(3, "d")
    rows = []
    for kind, pairs in (("excess", test), ("no excess", ctrl)):
        for a, b in zip(pairs["unit_a"], pairs["unit_b"]):
            ch = np.asarray(rec.channel_ids)[rec_group == group[a]]
            lag_ab = partners(tr[a], tr[b], -NEAR, NEAR)
            iso_a = partners(tr[a], tr[b], -ISO, ISO) == np.iinfo(np.int64).min
            iso_b = partners(tr[b], tr[a], -ISO, ISO) == np.iinfo(np.int64).min
            Ta = snippets(rec, rng.choice(tr[a][iso_a], min(N_EVENTS, iso_a.sum()), replace=False), ch).mean(0)
            Tb = snippets(rec, rng.choice(tr[b][iso_b], min(N_EVENTS, iso_b.sum()), replace=False), ch).mean(0)
            co = np.flatnonzero(lag_ab != np.iinfo(np.int64).min)
            pick = rng.choice(co, min(N_EVENTS, len(co)), replace=False)
            alone = rng.choice(np.flatnonzero(iso_a), min(N_EVENTS, int(iso_a.sum())), replace=False)
            main = int(np.argmax(np.abs(Ta).max(0)))
            peak = int(np.argmax(np.abs(Ta[:, main])))
            for events, label in ((pick, "coincident"), (alone, "a alone")):
                X = snippets(rec, tr[a][events], ch)
                lags = lag_ab[events] if label == "coincident" else np.zeros(len(events), np.int64)
                best = []
                for x, lag in zip(X, lags):
                    Tbs = shift(Tb, int(lag))
                    err = [((x - Ta) ** 2).sum(), ((x - Tbs) ** 2).sum(), ((x - Ta - Tbs) ** 2).sum()]
                    best.append(int(np.argmin(err)))
                best = np.bincount(best, minlength=3) / len(best)
                rows.append({"pairs": kind, "unit_a": a, "unit_b": b, "tetrode": int(group[a]) + 1,
                             "cosine": float(z.loc[(z.unit_a == a) & (z.unit_b == b), "cosine"].iloc[0]),
                             "events": label, "n": len(events), "best_a": best[0], "best_b": best[1], "best_sum": best[2],
                             "peak_over_a_template": float(np.median(X[:, peak, main] / Ta[peak, main])),
                             "b_template_at_a_peak": float(Tb[peak, main] / Ta[peak, main])})
            print(rows[-2], rows[-1], sep="\n", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "double_fits.csv", index=False)
    pd.set_option("display.width", 250)
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
