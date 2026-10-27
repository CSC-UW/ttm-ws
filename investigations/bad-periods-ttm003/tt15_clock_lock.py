"""TT15's units and the 300-sample (100 Hz) component of its traces.

Per unit with 1000 spikes or more in the first 14 windows of the baseline and
of the masked run: the share of spikes in each tenth of a 300-sample cycle
counted from the first sample. Then the 300-sample cycle mean of each channel
over 60 s at four times (peak to peak over the channel's MAD), the CH22 cycle
waveform's lag against its first-hour shape at six times, and the tenth of the
cycle that holds most spikes of baseline units 43 and 74 in each 6 h block of
the 48 h `assembled_prod`.

Run: .venv/bin/python investigations/bad-periods-ttm003/tt15_clock_lock.py
"""
import numpy as np
import spikeinterface as si

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
FS = 30000
P = 300
N = int(14 * 1800 * FS)

rec = si.load(f"{S}/mp_baseline/binary")
ch = list(rec.channel_ids)
tt = rec.get_property("tetrode")
gt = dict(zip(rec.get_property("group"), tt))

for name, d in [("baseline", "mp_baseline"), ("masked", "mp_baseline_masked_w14")]:
    s = si.load(f"{S}/{d}/assembled_reseed_c12")
    print(f"{name}, first 14 windows: share of spikes per tenth of the 300-sample cycle")
    for u, g in zip(s.unit_ids, s.get_property("group")):
        t = s.get_unit_spike_train(u)
        t = t[t < N]
        if t.size < 1000:
            continue
        h = np.bincount((t % P) * 10 // P, minlength=10) / t.size
        isi = np.diff(t) / FS * 1000
        on10 = np.mean((np.abs(isi / 10 - np.round(isi / 10)) * 10 < 0.1) & (isi > 5))
        mark = "  <-- locked" if h.max() > 0.2 else ""
        print(f"  {gt[g]:5s} unit {int(u):3d} n {t.size:8d}  max tenth {h.max():.3f}  empty tenths {(h < 0.01).sum()}"
              f"  ISI within 0.1 ms of k*10 ms {on10:.3f}{mark}")

print("\n300-sample cycle mean over 60 s, peak to peak / MAD")
for t0 in (1000, 30000, 100000, 160000):
    x = rec.get_traces(start_frame=t0 * FS, end_frame=(t0 + 60) * FS).astype(np.float32)
    mad = np.median(np.abs(x - np.median(x, 0)), 0) / 0.6745
    f = x[:(len(x) // P) * P].reshape(-1, P, x.shape[1]).mean(0)
    r = np.ptp(f, 0) / mad
    top = np.argsort(-r)[:6]
    print(f"  t={t0} s: " + ", ".join(f"{ch[i]} ({tt[i]}) {r[i]:.2f}" for i in top)
          + f" | median channel {np.median(r):.3f}")

print("\nCH22 cycle waveform against its shape at 0.3 h")
w = {}
for t0 in (1000, 20000, 60000, 100000, 140000, 170000):
    x = rec.get_traces(start_frame=t0 * FS, end_frame=(t0 + 60) * FS, channel_ids=["CH22"]).astype(np.float32)[:, 0]
    w[t0] = x[:(len(x) // P) * P].reshape(-1, P).mean(0)
ref = w[1000] - w[1000].mean()
for t0, v in w.items():
    v = v - v.mean()
    xc = np.array([np.dot(ref, np.roll(v, -k)) for k in range(P)])
    k = int(np.argmax(xc))
    print(f"  {t0 / 3600:5.1f} h: lag {k - P if k > P // 2 else k:+d} samples, "
          f"corr {xc.max() / np.linalg.norm(ref) / np.linalg.norm(v):.3f}, trough at sample {int(np.argmin(v))}")

s = si.load(f"{S}/mp_baseline/assembled_prod")
print("\nassembled_prod, tenth of the cycle with most spikes per 6 h block")
for u in (43, 74):
    t = s.get_unit_spike_train(u)
    blk = t // int(6 * 3600 * FS)
    print(f"  unit {u}: " + " ".join(
        str(int(np.argmax(np.bincount((t[blk == b] % P) * 10 // P, minlength=10)))) if (blk == b).any() else "-"
        for b in range(8)))
