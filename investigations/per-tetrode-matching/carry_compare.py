"""Joint against per-tetrode circus-omp over the first 14 windows of the carry-forward, masked and unmasked.

Runs are `mp_baseline.py --max-windows 14` carries (`mp_baseline` is the 48 h
run, cut to its first 14 windows). Single units only, clean time only (outside
the unit's tetrode's rows of the bad-period table). Per run and tetrode, then
per class of tetrodes: units with 20 clean spikes or more, clean spikes, the
share of the baseline's detected events at 10 MAD or more that a unit of the
same tetrode claims within 0.5 ms, the share of spikes with an interval under
1.5 ms to the unit's previous spike, the share with a spike of the same run on
another tetrode within 0.2 ms, and the median cosine of each spike to its
unit's window-local template (`purity_scores_*_w14.npz`, written by
`bad_periods_compare.py --max-windows 14`).

Run: .venv/bin/python investigations/per-tetrode-matching/carry_compare.py \
         joint=mp_baseline pertet=mp_baseline_pertet_w14 \
         joint_masked=mp_baseline_masked_w14 pertet_masked=mp_baseline_masked_pertet_w14
"""
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, "/Users/gfindlay@ad.wisc.edu/projects/ttm-ws/investigations/spike-sorting")
from _wobble_eval import TOL_MS, _within_tol
from ttm import bad_periods as bp

S = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix"
TABLE = "/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/bad_periods/bad_periods.csv"
SORTING = "assembled_reseed_c12"
FS = 30000.0
NFR = int(14 * 1800 * FS)
RP = int(1.5e-3 * FS)
COINC = int(0.2e-3 * FS)
CLASSES = {"TT4": [3], "TT5-TT12": list(range(4, 12)), "few-unit": [0, 1, 2, 12, 13, 14, 15]}
runs = dict(r.split("=") for r in sys.argv[1:])
table = bp.read_table(TABLE)
tol = int(TOL_MS * 1e-3 * FS)

z = np.load(f"{S}/mp_baseline/spike_coverage.npz")
m = (z["amp_mad"] >= 10) & (z["peak_sample"] < NFR)
ps, pg = z["peak_sample"][m], z["peak_group"][m].astype(int)
ok = np.ones(len(ps), bool)
for g in range(16):
    k = pg == g
    ok[k] = ~bp.bad_mask(table, g, ps[k], FS)
ps, pg = ps[ok], pg[ok]

rows = {}
for name, d in runs.items():
    s = si.load(f"{S}/{d}/{SORTING}")
    sc = np.load(f"{S}/{d}/purity_scores_{SORTING}_w14.npz")
    mua = (np.asarray(s.get_property("is_mua")).astype(bool) if "is_mua" in s.get_property_keys()
           else np.zeros(len(s.unit_ids), bool))
    grp = np.asarray(s.get_property("group")).astype(int)
    fr, ug, cos, short = [], [], [], []
    n_units = np.zeros(16, int)
    for u, g, mu in zip(s.unit_ids, grp, mua):
        if mu:
            continue
        f, r = sc[f"frames_{u}"], sc[f"rs_{u}"]
        k = f < NFR
        f, r = f[k], r[k]
        o = np.argsort(f)
        f, r = f[o], r[o]
        sh = np.r_[False, np.diff(f) < RP]
        c = ~bp.bad_mask(table, int(g), f, FS)
        n_units[g] += c.sum() >= 20
        fr.append(f[c])
        ug.append(np.full(c.sum(), g))
        cos.append(r[c])
        short.append(sh[c])
    fr, ug, cos, short = map(np.concatenate, (fr, ug, cos, short))
    by_tet = {g: np.sort(fr[ug == g]) for g in range(16)}
    claimed = _within_tol(ps, pg, by_tet, tol)
    other = np.zeros(len(fr), bool)
    for g in range(16):
        k = ug == g
        st = np.sort(fr[~k])
        if k.any() and len(st):
            j = np.searchsorted(st, fr[k])
            d_ = np.minimum(np.abs(fr[k] - st[np.clip(j - 1, 0, len(st) - 1)]), np.abs(st[np.clip(j, 0, len(st) - 1)] - fr[k]))
            other[k] = d_ <= COINC
    rows[name] = (n_units, ug, claimed, short, other, cos)


def line(name, gs):
    n_units, ug, claimed, short, other, cos = rows[name]
    k, e = np.isin(ug, gs), np.isin(pg, gs)
    cl = f"{claimed[e].mean() * 100:7.1f}" if e.any() else "      -"
    if not k.any():
        return f"  {name:14s} {n_units[gs].sum():5d} {0:9d} {cl}"
    return (f"  {name:14s} {n_units[gs].sum():5d} {k.sum():9d} {cl} {short[k].mean() * 100:8.2f} "
            f"{other[k].mean() * 100:8.1f} {np.nanmedian(cos[k]):7.3f}")


head = f"  {'run':14s} {'units':>5s} {'spikes':>9s} {'claim%':>7s} {'isi<1.5%':>8s} {'xtet%':>8s} {'cos_med':>7s}"
print(f"windows 0 to 13, clean time, single units; {len(ps)} clean events at 10 MAD or more")
for cname, gs in CLASSES.items():
    print(f"\n{cname}\n{head}")
    for name in runs:
        print(line(name, gs))
for g in range(16):
    print(f"\nTT{g + 1} ({(pg == g).sum()} events)\n{head}")
    for name in runs:
        print(line(name, [g]))
