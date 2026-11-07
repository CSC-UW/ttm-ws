"""Share of the refractory-rejected fits (every 8th window of TTM-003-48H) with an event of their tetrode's
multi-unit unit within 10 samples, by lag to the survivor."""
import numpy as np
import spikeinterface as si
RUN = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main"
s = si.load(f"{RUN}/sorting_mua")
fs = s.get_sampling_frequency()
RF, TOL = int(1.5e-3 * fs), 10
mua_ids = [u for u in s.unit_ids if s.get_unit_property(u, "is_mua")]
mua = {int(s.get_unit_property(u, "group")): np.sort(s.get_unit_spike_train(u)) for u in mua_ids}
hit_all, lag_all = [], []
for w in range(0, 96, 8):
    z = np.load(f"{RUN}/carry/spikes_w{w:02d}.npz")
    bk = np.load(f"{RUN}/carry/bankin_w{w:02d}.npz")
    gmap = dict(zip(bk["ids"].tolist(), bk["group"].tolist()))
    keep = np.isin(z["reason"], (0, 3))
    f, u, reason = z["frame"][keep], z["unit"][keep], z["reason"][keep]
    o = np.lexsort((f, u))
    f, u, reason = f[o], u[o], reason[o]
    close = np.r_[False, (np.diff(f) < RF) & (u[1:] == u[:-1])]
    chain = np.cumsum(~close) - 1
    acc = reason == 0
    sf = np.full(chain.max() + 1, -1, np.int64)
    sf[chain[acc]] = f[acc]
    rej = np.flatnonzero((reason == 3) & (sf[chain] >= 0))
    g = np.array([gmap[x] for x in u[rej].tolist()])
    hit = np.zeros(rej.size, bool)
    for gg in np.unique(g):
        if gg not in mua:
            continue
        m = g == gg
        M = mua[gg]
        hit[m] = (np.searchsorted(M, f[rej][m] + TOL, "right") - np.searchsorted(M, f[rej][m] - TOL)) > 0
    hit_all.append(hit)
    lag_all.append(np.abs(f[rej] - sf[chain[rej]]))
hit, lag = np.concatenate(hit_all), np.concatenate(lag_all) / fs * 1e3
print(f"{hit.size:,} rejected fits; with a multi-unit event within {TOL} samples: {hit.mean():.3f}")
for lo, hi in ((0, 0.5), (0.5, 1.0), (1.0, 1.5)):
    m = (lag >= lo) & (lag < hi)
    print(f"  lag {lo}-{hi} ms: {hit[m].mean():.3f} of {m.sum():,}")
