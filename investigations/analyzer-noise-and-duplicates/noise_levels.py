"""How many random chunks SpikeInterface's noise levels need on the final TTM-003-48H binary, with and
without each tetrode's bad intervals.

    .venv/bin/python investigations/analyzer-noise-and-duplicates/noise_levels.py <scratch dir>

`get_noise_levels` (method "mad") takes the MAD of each of 20 random 0.5 s chunks and averages them.
One pass over the int16 binary computes the MAD of every 0.5 s chunk of a grid on every channel
(<scratch dir>/chunk_mad.npz), so the mean over all chunks, with or without the chunks that overlap
a channel's bad intervals, is exact, and draws of any size are subsets of that table. The pass also
keeps a histogram of every channel's samples for the MAD of the pooled recording.
"""

import pathlib
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import spikeinterface as si
from spikeinterface.core import get_noise_levels

SESSION = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12"
BINARY = f"{SESSION}/sortings_carry/main/binary"
BAD = f"{SESSION}/bad_periods/bad_periods.csv"
OUT = pathlib.Path(__file__).parent / "results"
CHUNK = 15000
READ = 20 * CHUNK
K_LIST = (20, 100, 500, 2000, 10000)
N_SEEDS = 100
SI_SEEDS = (0, 1, 2, 3, 4)
MAD_TO_SD = 0.6744897501960817


def bad_chunks(n_chunks, groups):
    """(chunk, channel) True where the 0.5 s chunk overlaps a bad interval of the channel's tetrode."""
    t = pd.read_csv(BAD)
    out = np.zeros((n_chunks, len(groups)), bool)
    for g, a, b in zip(t["group"], t["start_s"], t["end_s"]):
        lo, hi = int(a * 30000) // CHUNK, int(np.ceil(b * 30000 / CHUNK))
        out[lo:min(hi, n_chunks), groups == g] = True
    return out


def block(args):
    first, last, bad = args
    rec = si.load(BINARY)
    n_ch = rec.get_num_channels()
    mad = np.empty((last - first, n_ch), np.float32)
    h_all = np.zeros((n_ch, 65536), np.int64)
    h_clean = np.zeros((n_ch, 65536), np.int64)
    for c0 in range(first, last, READ // CHUNK):
        c1 = min(last, c0 + READ // CHUNK)
        x = rec.get_traces(start_frame=c0 * CHUNK, end_frame=c1 * CHUNK, return_in_uV=False)
        x = x.reshape(c1 - c0, CHUNK, n_ch)
        med = np.median(x, axis=1, keepdims=True)
        mad[c0 - first:c1 - first] = np.median(np.abs(x - med), axis=1) / MAD_TO_SD
        for ch in range(n_ch):
            v = x[:, :, ch].astype(np.int32) + 32768
            h_all[ch] += np.bincount(v.ravel(), minlength=65536)
            good = ~bad[c0 - first:c1 - first, ch]
            if good.any():
                h_clean[ch] += np.bincount(v[good].ravel(), minlength=65536)
    return first, mad, h_all, h_clean


def pooled_mad(h):
    """MAD / 0.6745 of the samples whose int16 histogram (offset 32768) is ``h``."""
    c = np.cumsum(h)
    med = np.searchsorted(c, c[-1] / 2) - 32768
    dev = np.abs(np.arange(65536) - 32768 - med)
    order = np.argsort(dev, kind="stable")
    cd = np.cumsum(h[order])
    return dev[order][np.searchsorted(cd, cd[-1] / 2)] / MAD_TO_SD


def main():
    scratch = pathlib.Path(sys.argv[1])
    scratch.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(exist_ok=True)
    rec = si.load(BINARY)
    groups = np.asarray(rec.get_property("group")).astype(int)
    n_chunks = rec.get_num_frames() // CHUNK
    bad = bad_chunks(n_chunks, groups)
    table = scratch / "chunk_mad.npz"
    if table.exists():
        z = np.load(table)
        mad, h_all, h_clean = z["mad"], z["h_all"], z["h_clean"]
    else:
        t0 = time.time()
        edges = np.linspace(0, n_chunks, 49).astype(int)
        mad = np.empty((n_chunks, rec.get_num_channels()), np.float32)
        h_all = np.zeros((rec.get_num_channels(), 65536), np.int64)
        h_clean = np.zeros_like(h_all)
        with ProcessPoolExecutor(48) as ex:
            for first, m, ha, hc in ex.map(block, [(a, b, bad[a:b]) for a, b in zip(edges[:-1], edges[1:])]):
                mad[first:first + len(m)] = m
                h_all += ha
                h_clean += hc
        np.savez(table, mad=mad, h_all=h_all, h_clean=h_clean)
        print(f"pass over the binary in {(time.time() - t0) / 60:.1f} min", flush=True)

    ch = np.asarray(rec.channel_ids)
    exact_all = mad.mean(0)
    exact_clean = np.array([mad[~bad[:, c], c].mean() for c in range(len(ch))])
    per_channel = pd.DataFrame({
        "channel": ch, "tetrode": groups + 1, "bad_chunk_fraction": bad.mean(0),
        "mean_chunk_mad_all": exact_all, "mean_chunk_mad_clean": exact_clean,
        "pooled_mad_all": [pooled_mad(h) for h in h_all], "pooled_mad_clean": [pooled_mad(h) for h in h_clean]})
    per_channel["all_over_clean"] = per_channel["mean_chunk_mad_all"] / per_channel["mean_chunk_mad_clean"]
    per_channel["pooled_over_mean_clean"] = per_channel["pooled_mad_clean"] / per_channel["mean_chunk_mad_clean"]
    per_channel.to_csv(OUT / "noise_levels_channels.csv", index=False)

    rows = []
    rng = np.random.default_rng(0)
    for k in K_LIST:
        err_all, err_clean = [], []
        for _ in range(N_SEEDS):
            pick = rng.choice(n_chunks, size=k, replace=False)
            m, b = mad[pick], bad[pick]
            err_all.append(m.mean(0) / exact_clean - 1)
            err_clean.append(np.array([m[~b[:, c], c].mean() for c in range(len(ch))]) / exact_clean - 1)
        for name, e in (("all chunks", np.array(err_all)), ("clean chunks", np.array(err_clean))):
            rows.append({"chunks": k, "seconds": k * CHUNK / 30000, "draw": name,
                         "rel_sd_median_channel": float(np.median(e.std(0))),
                         "rel_sd_max_channel": float(e.std(0).max()),
                         "abs_err_p95_median_channel": float(np.median(np.percentile(np.abs(e), 95, axis=0))),
                         "abs_err_p95_max_channel": float(np.percentile(np.abs(e), 95, axis=0).max()),
                         "bias_max_channel": float(np.abs(e.mean(0)).max())})
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "noise_levels_draws.csv", index=False)

    si_rows = []
    for seed in SI_SEEDS:
        t0 = time.time()
        nl = get_noise_levels(rec, return_in_uV=False, method="mad", force_recompute=True,
                              random_slices_kwargs={"num_chunks_per_segment": 20, "chunk_duration": "500ms",
                                                    "seed": seed}, n_jobs=8)
        si_rows.append({"seed": seed, "seconds": time.time() - t0,
                        **{c: float(v) for c, v in zip(ch, nl / exact_clean - 1)}})
    si_df = pd.DataFrame(si_rows)
    si_df.to_csv(OUT / "noise_levels_spikeinterface_20.csv", index=False)
    si_err = si_df[list(ch)].to_numpy()
    pd.set_option("display.width", 250)
    print(per_channel.round(4).to_string(index=False))
    print(s.round(4).to_string(index=False))
    print(f"get_noise_levels, 20 chunks, {len(SI_SEEDS)} seeds: median channel SD of the relative error "
          f"{np.median(si_err.std(0)):.4f}, largest |error| {np.abs(si_err).max():.4f}")


if __name__ == "__main__":
    main()
