"""Run one arm of the carry-forward bake-off over windows of the hybrid binary.

Every arm starts from ``bank_start.npz``: mp_baseline's seed bank (45 units,
MountainSort5 on window 0, snr 5 or more, 100 spikes or more, deduplicated at
cosine 0.95) plus one template per ground-truth unit (ids 1000 and up), the
mean of up to 500 of its injected spikes in window 0 of the hybrid binary.
``--make-bank`` writes that file. Output per window goes to
``<_bake.BAKE>/arms/<name>/``.

Run: .venv/bin/python investigations/carry-bakeoff/run_arm.py NAME [--windows 0-11] [--n-jobs 16] [--workers 8]
     .venv/bin/python investigations/carry-bakeoff/run_arm.py --make-bank
"""
from __future__ import annotations

import argparse
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
import _bake as bk  # noqa: E402
from _bake import Arm, Recipe  # noqa: E402
from ttm import capture as cp  # noqa: E402

JUNE = Recipe()
FROZEN = Recipe(blend=0.0)
RECIPES = {
    "june": JUNE,                                                     # all accepted, raw snippets, mean
    "pl": Recipe(snippet="peeled"),                                   # other fitted spikes subtracted
    "ln": Recipe(select="lead_nostim", snippet="peeled"),             # no same-unit spike in the 10 ms before, no laser
    "ln_med": Recipe(select="lead_nostim", snippet="peeled", estimator="median"),
    "ln_trim": Recipe(select="lead_nostim", snippet="peeled", estimator="trimmed"),
    "ln_b5": Recipe(select="lead_nostim", snippet="peeled", blend=0.5),
    "full": Recipe(select="full"),                                    # a >= 0.8 only (Track G's w35a)
    "iso": Recipe(select="isolated"),                                 # no other unit's spike within 1.5 ms
    "frozen": FROZEN,
    # bad periods of the tetrode excluded from re-estimation ("match but don't learn")
    "lnb": Recipe(select="lead_nostim_nobad", snippet="peeled"),
    # and the scale corrected for the truncation the admission floors make (truncated-normal mean of snr)
    "lnbt": Recipe(select="lead_nostim_nobad", snippet="peeled", scale="trunc"),
    "lnbt_b5": Recipe(select="lead_nostim_nobad", snippet="peeled", scale="trunc", blend=0.5),
    "t": Recipe(scale="trunc"),
}
ARMS = {a.name: a for a in [
    # Stage A, windows 0 and 1 on the starting bank, templates frozen: what each matcher admits.
    # June: circus-omp on all channels, a >= 0.8; masked as in production; per tetrode at June's
    # floor; per tetrode with a permissive floor (a >= 0.3 and snr >= 2.5) for rules applied afterwards
    Arm("A_joint", per_tetrode=False, recipe=FROZEN),
    Arm("A_joint_masked", per_tetrode=False, masked=True, recipe=FROZEN),
    Arm("A_pt", recipe=FROZEN),
    Arm("A_pt_open", amp_rel=0.3, snr_floor=2.5, recipe=FROZEN),
    Arm("A_pt_open_white", amp_rel=0.3, snr_floor=2.5, recipe=FROZEN, whiten=True),
    # Stage B, carries over windows 0 to 11
    Arm("joint_june", per_tetrode=False),
    Arm("joint_june_masked", per_tetrode=False, masked=True),
    Arm("pt_june"),
] + [
    # admission from stage A (snr floor, relative floor 0.5, peeled cosine 0.6, refractory 1.5 ms) x recipe
    Arm(f"pt_{adm}_{rname}", amp_rel=0.5, snr_floor=k, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=r)
    for adm, k in (("s4", 4.0), ("s45", 4.5)) for rname, r in RECIPES.items()
    if (adm == "s4" and rname not in ("t",)) or rname in ("june", "ln", "ln_b5", "lnbt")
] + [
    Arm("pt_june_ln", recipe=RECIPES["ln"]),
    # the snr floor against the projection noise of the window's residual (busy tetrodes' own spikes removed)
    Arm("pt_s4r_ln", amp_rel=0.5, snr_floor=4.0, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES["ln"],
        sigma_from="resid"),
    Arm("pt_s45r_ln", amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES["ln"],
        sigma_from="resid"),
    Arm("joint_june_masked_ln", per_tetrode=False, masked=True, recipe=RECIPES["ln"]),
    Arm("joint_june_masked_t", per_tetrode=False, masked=True, recipe=RECIPES["t"]),
    Arm("joint_june_masked_lnbt", per_tetrode=False, masked=True, recipe=RECIPES["lnbt"]),
    Arm("pt_s4r_lnbt", amp_rel=0.5, snr_floor=4.0, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES["lnbt"],
        sigma_from="resid"),
] + [
    # Stage D, the candidates: permissive matcher (a >= 0.3, snr >= 2.5), admission afterwards (snr floor,
    # a >= 0.5, peeled cosine 0.6, refractory 1.5 ms, fits in the tetrode's own bad periods dropped), the
    # relaxed in-pulse rule for units that fire through the pulses, and re-estimation from lead spikes
    # outside pulses and bad periods with the scale corrected for truncation
    Arm(f"D_{adm}_{rname}", amp_rel=0.5, snr_floor=k, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES[rname],
        match_floor=(0.3, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True)
    for adm, k in (("s4", 4.0), ("s45", 4.5)) for rname in ("lnbt", "lnbt_b5", "frozen")
]}


def make_bank() -> None:
    rec = si.load(bk.HYBRID)
    seed = bk.Bank.load(bk.BAKE_ROOT / "bank0.npz")
    gt = np.load(bk.BAKE / "gt.npz")
    rng = np.random.default_rng(0)
    ids, group, T = list(seed.ids), list(seed.group), list(seed.T)
    for i, g in enumerate(gt["group"]):
        fr = gt["frame"][(gt["unit"] == i) & (gt["frame"] < bk.WIN)]
        fr = fr[(fr > bk.NB) & (fr < bk.WIN - bk.NA)]
        fr = np.sort(rng.choice(fr, min(500, fr.size), replace=False))
        x = cp.load_tetrode_traces(rec, int(g), 0, bk.WIN)
        T.append(cp.snippets(x, fr).mean(axis=0))
        ids.append(1000 + i)
        group.append(int(g))
        print(f"GT {i} TT{g + 1}: {fr.size} spikes, cosine to truth "
              f"{cp.template_cosine(T[-1], gt['T_win'][i, 0] * gt['size_win'][i, 0], max_shift=0):.3f}")
    bk.Bank(np.array(ids, np.int64), np.array(group), np.stack(T).astype(np.float32)).save(bk.BAKE / "bank_start.npz")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name", nargs="?")
    ap.add_argument("--make-bank", action="store_true")
    ap.add_argument("--windows", default="0-11")
    ap.add_argument("--n-jobs", type=int, default=16)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    if a.make_bank:
        make_bank()
        return
    lo, hi = (int(v) for v in a.windows.split("-"))
    arm = ARMS[a.name]
    bank = bk.Bank.load(bk.BAKE / "bank_start.npz")
    if arm.whiten:
        groups = np.asarray(si.load(bk.HYBRID if arm.hybrid else bk.BINARY).get_property("group")).astype(int)
        bank = bk.whiten_bank(bank, bk.whitening_matrix(arm.hybrid), groups)
    bk.log(f"{arm.name}: {arm}")
    bk.run_carry(arm, bank, bk.BAKE / "arms" / arm.name, windows=range(lo, hi + 1), n_jobs=a.n_jobs,
                 workers=a.workers)


if __name__ == "__main__":
    main()
