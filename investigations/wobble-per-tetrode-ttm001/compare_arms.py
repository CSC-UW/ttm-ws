"""circus-omp and wobble, each on all 64 channels at once and per tetrode, on a clean stretch of TTM-001.

Bank, settings and matchers as in `_stretch.py`. Arms: circus-omp joint and
per tetrode; wobble joint and per tetrode, before and after the cosine gate.
The stretch lies inside window 55; the bad-period table
(`bad_periods_crop_s2000_d170000`, times in the binary's frame) is checked
for flagged time in it.

Per spike: its tetrode, amplitude scale and cosine to its unit's template on
the unit's four wires (`_mp_common.per_spike_fit`), trough in excised-MAD
units (deepest sample over the snippet and the four wires, as `ttm.purity`
scores it), whether its unit's previous spike is within 1.5 ms, and whether
a spike of the same arm on another tetrode is within 0.2 ms. The chance level
of the last is, per tetrode, the share of the stretch's samples within 0.2 ms
of a spike of that arm on another tetrode; for a set of spikes it is that
share averaged over the set's spikes. Claimed: the share of the detected
events of `spike_coverage.npz` (negative peaks at 5.5 MAD per tetrode over
the 47 h binary) with a spike of the same tetrode within 0.5 ms, at 10 MAD
and up and from 5.5 to 10 MAD.

Sections: per arm, pooled and by unit class (template trough at the span
start under or at least 10 MAD, `template_amp_w55.npy`); per tetrode; per
tetrode against joint of the same matcher, spike for spike (same unit within
3 samples); wobble against circus-omp spike for spike, joint pair and per
tetrode pair; per unit; all arms on white noise at the clean noise level
(`mask_bad_periods` with every tetrode flagged over the stretch). Spike arrays
are saved under `track_eval/wobble_pertet/`.

Run: OPENBLAS_NUM_THREADS=16 .venv/bin/python investigations/wobble-per-tetrode-ttm001/compare_arms.py T0 T1
"""
import os
import sys

import numpy as np
import pandas as pd
import spikeinterface as si
from _stretch import (
    BINARY,
    COVERAGE,
    FS,
    OUT,
    TABLE,
    TEMPLATE_AMP,
    bank_templates,
    gate,
    match,
    method_kwargs,
    per_spike_fit,
)
from _wobble_eval import _within_tol
from ttm import bad_periods as bp
from ttm import capture as cp

RP = int(1.5e-3 * FS)
COINC = int(0.2e-3 * FS)
TOL = int(0.5e-3 * FS)
SAME = 3
T0, T1 = float(sys.argv[1]), float(sys.argv[2])
f0, f1 = int(T0 * FS), int(T1 * FS)
n = f1 - f0
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)

full = si.load(BINARY)
tmpl, ug = bank_templates(full)
ids = np.asarray(tmpl.unit_ids).astype(int)
amp_rows = np.load(TEMPLATE_AMP)
tmad = dict(zip(amp_rows[:, 0].astype(int).tolist(), amp_rows[:, 3].tolist()))
small = np.array([tmad[u] < 10 for u in ids])
table = bp.read_table(TABLE)
flagged = table[(table.start_s < T1) & (table.end_s > T0)]
thr = method_kwargs("wobble", tmpl)["parameters"]["threshold"]
print(f"{T0:.0f} to {T1:.0f} s of the binary (session {T0 + 2000:.0f} to {T1 + 2000:.0f} s), window 55; "
      f"{len(ids)} units, {int(small.sum())} with template trough under 10 MAD; flagged tetrodes here: "
      f"{sorted(set(flagged.tetrode)) or 'none'}; wobble threshold {thr:.0f} (0.35 x bank median energy), "
      f"the same in both wobble arms")
print("units per tetrode (small): " + ", ".join(f"TT{g + 1} {(ug == g).sum()} ({(small & (ug == g)).sum()})"
                                                 for g in range(16)))

seg = full.frame_slice(f0, f1)
seg.reset_times()
sp, wall = {}, {}
for matcher, tag in (("circus-omp", "omp"), ("wobble", "wob")):
    for per, how in ((False, "joint"), (True, "per-tet")):
        s, dt = match(seg, tmpl, matcher, per)
        if matcher == "wobble":
            sp[f"wob {how} ungated"], wall[f"wob {how} ungated"] = s, dt
            g_, _, _, dg = gate(s, tmpl, seg)
            sp[f"wob {how}"], wall[f"wob {how}"] = g_, dt + dg
        else:
            sp[f"omp {how}"], wall[f"omp {how}"] = s, dt
ARMS = ["omp joint", "omp per-tet", "wob joint ungated", "wob per-tet ungated", "wob joint", "wob per-tet"]
os.makedirs(OUT, exist_ok=True)
np.savez(f"{OUT}/spikes_{T0:.0f}_{T1:.0f}.npz", unit_ids=ids, groups=ug,
         **{a.replace(" ", "_"): sp[a] for a in ARMS})

traces = {g: cp.load_tetrode_traces(seg, g) for g in range(16)}
noise = {g: cp.noise_mad_excised(traces[g])[0] for g in range(16)}


def coincidence(st, tet):
    out, chance = np.zeros(len(st), bool), np.zeros(16)
    for g in range(16):
        other = np.sort(st[tet != g])
        m = tet == g
        if other.size == 0:
            continue
        if m.any():
            j = np.searchsorted(other, st[m])
            d = np.minimum(np.abs(st[m] - other[np.clip(j - 1, 0, len(other) - 1)]),
                           np.abs(other[np.clip(j, 0, len(other) - 1)] - st[m]))
            out[m] = d <= COINC
        cov = (np.bincount(np.clip(other - COINC, 0, n), minlength=n + 1)
               - np.bincount(np.clip(other + COINC + 1, 0, n), minlength=n + 1))
        chance[g] = (np.cumsum(cov[:-1]) > 0).mean()
    return out, chance


def features(s):
    st, ci = s["sample_index"].astype(np.int64), s["cluster_index"].astype(np.int64)
    tet = ug[ci]
    a, r = per_spike_fit(s, tmpl, seg)
    trough = np.full(len(s), np.nan)
    for g in range(16):
        k = np.flatnonzero((tet == g) & (st >= cp.NBEFORE) & (st + cp.NAFTER <= n))
        if k.size:
            trough[k] = -(cp.snippets(traces[g], st[k]) / noise[g][None, None, :]).min(axis=(1, 2))
    short = np.zeros(len(s), bool)
    for u in np.unique(ci):
        idx = np.flatnonzero(ci == u)
        o = idx[np.argsort(st[idx], kind="stable")]
        short[o[1:]] = np.diff(st[o]) < RP
    xtet, ch = coincidence(st, tet)
    return {"st": st, "ci": ci, "tet": tet, "small": small[ci], "a": a, "r": r, "trough": trough, "short": short,
            "xtet": xtet, "chance": ch[tet]}


F = {a: features(sp[a]) for a in ARMS}
z = np.load(COVERAGE)
k = (z["peak_sample"] >= f0) & (z["peak_sample"] < f1)
ps, pg, pa = z["peak_sample"][k] - f0, z["peak_group"][k].astype(int), z["amp_mad"][k]


def claimed(f):
    return _within_tol(ps, pg, {g: np.sort(f["st"][f["tet"] == g]) for g in range(16)}, TOL)


C = {a: claimed(F[a]) for a in ARMS}
big, mid = pa >= 10, pa < 10


def stats(f, sel):
    if not sel.any():
        return {"spikes": 0}
    return {"spikes": int(sel.sum()), "isi<1.5ms%": round(f["short"][sel].mean() * 100, 2),
            "cos_med": round(float(np.nanmedian(f["r"][sel])), 3),
            "cos<0.6%": round(float(np.mean(f["r"][sel] < 0.6) * 100), 1),
            "trough_mad_med": round(float(np.nanmedian(f["trough"][sel])), 2),
            "amp_scale_med": round(float(np.nanmedian(f["a"][sel])), 3),
            "xtet%": round(f["xtet"][sel].mean() * 100, 1), "chance%": round(f["chance"][sel].mean() * 100, 1)}


print(f"\nper arm, pooled ({int(big.sum())} detected events at 10 MAD and up, {int(mid.sum())} from 5.5 to 10 MAD)")
rows = []
for a in ARMS:
    f = F[a]
    rows.append({"arm": a, **stats(f, np.ones(len(f["st"]), bool)),
                 "claim10%": round(C[a][big].mean() * 100, 1), "claim5.5-10%": round(C[a][mid].mean() * 100, 1),
                 "wall_s": round(wall[a], 1)})
print(pd.DataFrame(rows).to_string(index=False))
for cls, want in (("small units (template under 10 MAD)", True), ("large units", False)):
    print(f"\nper arm, {cls}")
    print(pd.DataFrame([{"arm": a, **stats(F[a], F[a]["small"] == want)} for a in ARMS]).to_string(index=False))


def per_tet(fn, label):
    print(f"\n{label} per tetrode")
    print("tetrode " + "".join(f"{a:>21s}" for a in ARMS))
    for g in range(16):
        print(f"  TT{g + 1:<4d}" + "".join(f"{fn(a, g):>21s}" for a in ARMS))


per_tet(lambda a, g: f"{(F[a]['tet'] == g).sum():d}", "spikes")
per_tet(lambda a, g: f"{F[a]['short'][F[a]['tet'] == g].mean() * 100:.2f}", "intervals under 1.5 ms, %")
per_tet(lambda a, g: f"{np.nanmedian(F[a]['r'][F[a]['tet'] == g]):.3f}", "median cosine")
per_tet(lambda a, g: f"{F[a]['xtet'][F[a]['tet'] == g].mean() * 100:.1f} ({F[a]['chance'][F[a]['tet'] == g][:1].mean() * 100:.1f})",
        "spikes with a spike on another tetrode within 0.2 ms, % (chance)")
per_tet(lambda a, g: f"{C[a][big & (pg == g)].mean() * 100:.1f}" if (big & (pg == g)).any() else "-",
        "claimed events at 10 MAD and up, %")


def in_other(x, y):
    hit = np.zeros(len(x["st"]), bool)
    for u in np.unique(x["ci"]):
        mx = x["ci"] == u
        sy = np.sort(y["st"][y["ci"] == u])
        if sy.size:
            sx = x["st"][mx]
            j = np.searchsorted(sy, sx)
            d = np.minimum(np.abs(sx - sy[np.clip(j - 1, 0, len(sy) - 1)]), np.abs(sy[np.clip(j, 0, len(sy) - 1)] - sx))
            hit[mx] = d <= SAME
    return hit


def pair(xa, ya, xname, yname):
    x, y = F[xa], F[ya]
    xb, yb = in_other(x, y), in_other(y, x)
    print(f"\n{ya} against {xa}: {len(x['st'])} and {len(y['st'])} spikes; in both (same unit, 3 samples) "
          f"{int(yb.sum())}, {xname} only {int((~xb).sum())}, {yname} only {int((~yb).sum())}")
    rows = []
    for cls, cm in (("all", None), ("small units", True), ("large units", False)):
        for name, f, sel in ((f"both ({yname})", y, yb), (f"{xname} only", x, ~xb), (f"{yname} only", y, ~yb)):
            if cm is not None:
                sel = sel & (f["small"] == cm)
            rows.append({"units": cls, "category": name, **stats(f, sel)})
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"  per tetrode: {xname} / {yname} / both / {yname} only with median cosine, trough (MAD), "
          f"intervals under 1.5 ms %, on another tetrode within 0.2 ms % (chance)")
    for g in range(16):
        mx, my = x["tet"] == g, y["tet"] == g
        o = my & ~yb
        extra = (f"{np.nanmedian(y['r'][o]):.3f} {np.nanmedian(y['trough'][o]):5.1f} {y['short'][o].mean() * 100:5.1f} "
                 f"{y['xtet'][o].mean() * 100:5.1f} ({y['chance'][o].mean() * 100:4.1f})") if o.any() else "-"
        print(f"  TT{g + 1:<4d} {mx.sum():7d} {my.sum():7d} {(my & yb).sum():7d} {o.sum():7d}  {extra}")


print("\n\nper tetrode against joint, same matcher")
pair("omp joint", "omp per-tet", "joint", "per-tet")
pair("wob joint ungated", "wob per-tet ungated", "joint", "per-tet")
pair("wob joint", "wob per-tet", "joint", "per-tet")
print("\n\nwobble against circus-omp")
pair("omp joint", "wob joint", "omp", "wob")
pair("omp per-tet", "wob per-tet", "omp", "wob")
pair("omp joint", "wob joint ungated", "omp", "wob")
pair("omp per-tet", "wob per-tet ungated", "omp", "wob")

print("\n\nper unit: spikes / intervals under 1.5 ms")
print("unit tetrode tmpl_mad " + "".join(f"{a:>21s}" for a in ARMS))
for k_, u in enumerate(ids):
    cells = "".join(f"{(F[a]['ci'] == k_).sum():>13d} {F[a]['short'][F[a]['ci'] == k_].sum():>7d}" for a in ARMS)
    print(f"{u:4d}   TT{ug[k_] + 1:<4d} {tmad[u]:7.1f} {cells}")

flag_all = pd.DataFrame({"group": range(16), "tetrode": [f"TT{g + 1}" for g in range(16)],
                         "start_s": T0 - 5.0, "end_s": T1 + 5.0, "kind": "null", "channels": "",
                         "max_loud": 0.0, "max_corr_excess": 0.0, "flagged_s": T1 - T0})
groups = np.asarray(full.get_property("group")).astype(int)
null = bp.mask_bad_periods(full, flag_all, noise_levels=bp.clean_noise_levels(
    full, [(int(a_ * FS), int(b_ * FS), np.flatnonzero(groups == g).tolist())
           for g, a_, b_ in table[["group", "start_s", "end_s"]].itertuples(index=False)]))
nseg = null.frame_slice(f0, f1)
nseg.reset_times()
print(f"\nwhite noise at the clean noise level, {T1 - T0:.0f} s")
for matcher, tag in (("circus-omp", "omp"), ("wobble", "wob")):
    for per, how in ((False, "joint"), (True, "per-tet")):
        s, dt = match(nseg, tmpl, matcher, per)
        out = {f"{tag} {how}": s}
        if matcher == "wobble":
            out = {f"{tag} {how} ungated": s, f"{tag} {how}": gate(s, tmpl, nseg)[0]}
        for name, s_ in out.items():
            cnt = np.bincount(ug[s_["cluster_index"]], minlength=16) if len(s_) else np.zeros(16, int)
            print(f"  {name:20s} {len(s_):7d} spikes; per tetrode {cnt.tolist()}")
