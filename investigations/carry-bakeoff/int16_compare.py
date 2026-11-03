import sys
import numpy as np
import pandas as pd

A = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix/carry_bakeoff/aligned48/A48_W_s45_c50_lnb_b2"
B = sys.argv[2]
windows = range(int(sys.argv[1]))
OUT = sys.argv[3] if len(sys.argv) > 3 else None


def key(z, m):
    return z["frame"][m].astype(np.int64) * 100000 + z["unit"][m]


rows, units = [], []
for w in windows:
    a = np.load(f"{A}/spikes_w{w:02d}.npz")
    b = np.load(f"{B}/spikes_w{w:02d}.npz")
    ka, kb = key(a, slice(None)), key(b, slice(None))
    common, ia, ib = np.intersect1d(ka, kb, return_indices=True)
    acc_a, acc_b = a["reason"] == 0, b["reason"] == 0
    kaa, kba = key(a, acc_a), key(b, acc_b)
    both = np.intersect1d(kaa, kba)
    # accepted spikes of one side with a same-unit accepted spike within 1 sample on the other
    def near(k1, k2):
        u1, f1 = k1 % 100000, k1 // 100000
        u2, f2 = k2 % 100000, k2 // 100000
        hit = 0
        for u in np.unique(u1):
            x, y = np.sort(f1[u1 == u]), np.sort(f2[u2 == u])
            if not len(y):
                continue
            j = np.clip(np.searchsorted(y, x), 1, len(y) - 1)
            d = np.minimum(np.abs(y[j] - x), np.abs(y[j - 1] - x))
            hit += int((d <= 1).sum())
        return hit
    da = a["a"][ia] - b["a"][ib]
    dsnr = a["snr"][ia] - b["snr"][ib]
    same_reason = float((a["reason"][ia] == b["reason"][ib]).mean())
    rows.append(dict(window=w, fits_f32=len(ka), fits_i16=len(kb), fits_same=len(common),
                     acc_f32=int(acc_a.sum()), acc_i16=int(acc_b.sum()), acc_same=len(both),
                     acc_f32_within1=near(kaa, kba), acc_i16_within1=near(kba, kaa),
                     reason_agree=round(same_reason, 6),
                     amp_absdiff_p99=float(np.percentile(np.abs(da), 99)), amp_absdiff_max=float(np.abs(da).max()),
                     snr_absdiff_p99=float(np.nanpercentile(np.abs(dsnr), 99)), snr_nan_same=bool((np.isnan(a["snr"][ia]) == np.isnan(b["snr"][ib])).all())))
    ua = pd.Series(a["unit"][acc_a]).value_counts()
    ub = pd.Series(b["unit"][acc_b]).value_counts()
    u = pd.concat([ua.rename("f32"), ub.rename("i16")], axis=1).fillna(0).astype(int)
    u["window"] = w
    units.append(u)
    ba, bb = np.load(f"{A}/bank_w{w:02d}.npz"), np.load(f"{B}/bank_w{w:02d}.npz")
    assert (ba["ids"] == bb["ids"]).all(), "bank ids differ"
    Ta, Tb = ba["T"].reshape(len(ba["ids"]), -1), bb["T"].reshape(len(bb["ids"]), -1)
    cos = (Ta * Tb).sum(1) / np.linalg.norm(Ta, axis=1) / np.linalg.norm(Tb, axis=1)
    rel = np.abs(Ta - Tb).max(1) / np.abs(Ta).max(1)
    rows[-1].update(bank_cos_min=float(cos.min()), bank_maxdiff_rel=float(rel.max()))

d = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(d.to_string(index=False))
u = pd.concat([x.rename_axis("unit").reset_index() for x in units], ignore_index=True)
u["diff"] = u["i16"] - u["f32"]
u["rel"] = u["diff"] / u["f32"].clip(lower=1)
print("\nper unit and window, accepted spikes: largest relative differences")
print(u.reindex(u["rel"].abs().sort_values(ascending=False).index).head(12).to_string())
print(f"\nunit-windows with any difference: {(u['diff'] != 0).sum()} of {len(u)}; "
      f"total accepted f32 {u['f32'].sum()} i16 {u['i16'].sum()}")
if OUT:
    d.to_csv(f"{OUT}_windows.csv", index=False)
    u.to_csv(f"{OUT}_units.csv", index=False)
