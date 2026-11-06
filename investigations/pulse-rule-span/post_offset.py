"""TTM-003-48H carry (run main, carried with the pulse rule's span to 100 ms after the offset): for each
window's sustained responders, the accepted spikes by time from the pulse, split into those the normal
rules accept and those only the pulse rule accepts, with their amplitude (relative to the unit's
median of rule 5.3.4) and cosine; and the rate of both in 2 ms bins around the offset.

It reads the carry's own decisions: ``reason_carry`` where ``ttm-sort readmit`` has replayed the window
(the replay keeps them there), else ``reason``. Writes post_offset_segments.csv and offset_profile.csv
next to this script.

Run from the workspace root: .venv/bin/python investigations/pulse-rule-span/post_offset.py
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).parent
D = pathlib.Path("/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_carry/main/carry")
FS = 30000.0
P = pd.read_csv("subjects/ttm-003/results/laser_pulses_2026-10-20_09-02-12.csv").sort_values("onset_s")
ON = (P.onset_s.to_numpy() * FS).astype(np.int64)
OFF = (P.offset_s.to_numpy() * FS).astype(np.int64)
POST, WIN = int(0.1 * FS), int(1800 * FS)
SEGS = [("onset 0-15 ms", "on", 0, 15), ("15 ms-offset", "mid", None, None), ("off +0-10", "off", 0, 10),
        ("off +10-20", "off", 10, 20), ("off +20-50", "off", 20, 50), ("off +50-100", "off", 50, 100)]
EDGES = np.arange(-20, 31, 2)


def classify(U, fr, a, snr, cos, rs, near):
    """Unit ``U``'s normal-rule pass (5.3.1 to 5.3.4, refractory aside) and its median of rule 5.3.4."""
    S = U & (snr >= 4.5) & (rs != 5) & (cos >= 0.5)
    idx = np.flatnonzero(S & ~near)
    idx = idx[np.argsort(fr[idx])]
    lead = np.r_[True, np.diff(fr[idx]) > int(0.010 * FS)]
    ref = float(np.median(a[idx[lead]])) if lead.sum() >= 30 else 1.0
    return S & (a >= 0.5 * ref), ref


rows, prof = [], []
for f in sorted(D.glob("bank_w*.npz")):
    w = int(f.stem[-2:])
    bk = np.load(f)
    drv = bk["driven"]
    if drv.size == 0:
        continue
    grp = dict(zip(bk["ids"], bk["group"]))
    a0, a1 = w * WIN, (w + 1) * WIN
    m = (ON >= a0) & (ON < a1)
    won, woff = ON[m], OFF[m]
    z = np.load(D / f"spikes_w{w:02d}.npz")
    sel = np.isin(z["unit"], drv)
    rkey = "reason_carry" if "reason_carry" in z.files else "reason"
    fr, u, amp, snr, cos, rs = (z[k][sel] for k in ("frame", "unit", "a", "snr", "cos_pl", rkey))
    k = np.searchsorted(won, fr, side="right") - 1
    kk = np.clip(k, 0, None)
    t_on = (fr - won[kk]) / FS * 1e3
    t_off = (fr - woff[kk]) / FS * 1e3
    near = (k >= 0) & (fr < woff[kk] + POST)
    nxt = won[np.clip(np.searchsorted(won, fr), 0, won.size - 1)]
    far = ~near & (np.minimum(np.abs(fr - nxt), np.abs(fr - won[kk])) / FS > 2.0)
    st = np.r_[won - int(2 * FS), won]
    en = np.r_[won + int(2 * FS), woff + POST]
    o = np.argsort(st)
    st, en = st[o], np.maximum.accumulate(en[o])
    new = np.r_[True, st[1:] > en[:-1]]
    ends = np.r_[en[np.flatnonzero(new)[1:] - 1], en[-1]]
    far_s = (a1 - a0) / FS - (np.minimum(ends, a1) - np.maximum(st[new], a0)).clip(0).sum() / FS
    j = np.clip(np.searchsorted(woff, fr), 0, woff.size - 1)
    jm = np.clip(j - 1, 0, None)
    t_near_off = np.where(np.abs(fr - woff[jm]) < np.abs(fr - woff[j]), fr - woff[jm], fr - woff[j]) / FS * 1e3
    for uu in drv:
        U = u == uu
        normal, ref = classify(U, fr, amp, snr, cos, rs, near)
        acc = U & (rs == 0)
        rel_only = acc & near & ~normal
        bg = U & (rs != 5) & (snr >= 2.5) & (amp >= 0.3 * ref) & (cos >= 0.5) & ~normal & far
        base = {"window": w, "unit": int(uu), "tt": int(grp[uu]) + 1, "pulses": int(won.size),
                "rate_far_normal": (normal & acc & far).sum() / far_s,
                "rate_far_relax_eligible": bg.sum() / far_s,
                "cos_far_normal": float(np.median(cos[normal & acc & far])),
                "cos_far_relax_eligible": float(np.median(cos[bg])) if bg.any() else np.nan}
        for name, kind, lo, hi in SEGS:
            if kind == "on":
                s, dur = near & (t_on >= lo) & (t_on < hi), (hi - lo) / 1e3
            elif kind == "mid":
                s, dur = near & (t_on >= 15) & (t_off < 0), float(np.mean((woff - won) / FS)) - 0.015
            else:
                s, dur = near & (t_off >= lo) & (t_off < hi), (hi - lo) / 1e3
            nn, nr = acc & s & normal, rel_only & s
            rows.append({**base, "seg": name, "n_normal": int(nn.sum()), "n_relax_only": int(nr.sum()),
                         "rate_normal": nn.sum() / (won.size * dur), "rate_relax_only": nr.sum() / (won.size * dur),
                         "amp_rel_normal": float(np.median(amp[nn]) / ref) if nn.any() else np.nan,
                         "amp_rel_relax": float(np.median(amp[nr]) / ref) if nr.any() else np.nan,
                         "cos_normal": float(np.median(cos[nn])) if nn.any() else np.nan,
                         "cos_relax": float(np.median(cos[nr])) if nr.any() else np.nan})
        h_n = np.histogram(t_near_off[acc & normal], EDGES)[0]
        h_r = np.histogram(t_near_off[acc & ~normal], EDGES)[0]
        for e, a_, b_ in zip(EDGES[:-1], h_n, h_r):
            prof.append({"window": w, "unit": int(uu), "pulses": int(won.size), "ms_from_offset": int(e),
                         "n_normal": int(a_), "n_relax_only": int(b_)})
    print(w, drv.tolist(), flush=True)
pd.DataFrame(rows).to_csv(HERE / "post_offset_segments.csv", index=False)
pd.DataFrame(prof).to_csv(HERE / "offset_profile.csv", index=False)
