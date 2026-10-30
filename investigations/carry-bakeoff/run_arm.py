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
import pathlib
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
    "lnb_b5": Recipe(select="lead_nostim_nobad", snippet="peeled", blend=0.5),
    "lnb_b2": Recipe(select="lead_nostim_nobad", snippet="peeled", blend=0.2),
    # only confident spikes: snr 1.5 above the admission floor and peeled cosine 0.7 or more
    "conf": Recipe(select="lead_nostim_nobad", snippet="peeled", conf_snr=1.5, conf_cos=0.7),
    "conf_b5": Recipe(select="lead_nostim_nobad", snippet="peeled", conf_snr=1.5, conf_cos=0.7, blend=0.5),
    # and the spikes kept in a band symmetric about the mode of their snr, so the admission floors' cut of
    # the amplitude distribution does not raise the template (a truncated-normal fit, scale="trunc",
    # shrank real units' templates to the 0.5 clip: their snr is not one truncated normal)
    "lnbs": Recipe(select="lead_nostim_nobad", snippet="peeled", scale="sym"),
    "lnbs_b5": Recipe(select="lead_nostim_nobad", snippet="peeled", scale="sym", blend=0.5),
    "s": Recipe(scale="sym"),
    # rebuilt templates shifted back so their reference extremum stays at sample 30
    "lnb_b2_rc": Recipe(select="lead_nostim_nobad", snippet="peeled", blend=0.2, recenter=True),
    "june_rc": Recipe(recenter=True),
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
    if (adm == "s4" and rname not in ("s",)) or rname in ("june", "ln", "ln_b5", "lnbs")
] + [
    Arm("pt_june_ln", recipe=RECIPES["ln"]),
    # the snr floor against the projection noise of the window's residual (busy tetrodes' own spikes removed)
    Arm("pt_s4r_ln", amp_rel=0.5, snr_floor=4.0, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES["ln"],
        sigma_from="resid"),
    Arm("pt_s45r_ln", amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES["ln"],
        sigma_from="resid"),
    Arm("joint_june_masked_ln", per_tetrode=False, masked=True, recipe=RECIPES["ln"]),
    Arm("joint_june_masked_s", per_tetrode=False, masked=True, recipe=RECIPES["s"]),
    Arm("joint_june_masked_lnbs", per_tetrode=False, masked=True, recipe=RECIPES["lnbs"]),
    Arm("pt_s4r_lnbs", amp_rel=0.5, snr_floor=4.0, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES["lnbs"],
        sigma_from="resid"),
] + [
    # Stage D, the candidates: permissive matcher (a >= 0.2, snr >= 2.5), admission afterwards (snr floor,
    # peeled cosine 0.6, a >= 0.5 of the unit's median fit amplitude in the window, refractory 1.5 ms, fits
    # in the tetrode's own bad periods dropped), the relaxed in-pulse rule for units that fire through the
    # pulses, and re-estimation from lead spikes outside pulses and bad periods (peeled, mean)
    Arm(f"D_{adm}_{rname}", amp_rel=0.5, snr_floor=k, gate=("cos_pl", 0.6), refractory_ms=1.5, recipe=RECIPES[rname],
        match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True, amp_ref="median")
    for adm, k in (("s4", 4.0), ("s45", 4.5)) for rname in ("lnb", "lnb_b5", "lnb_b2", "june", "frozen", "conf")
    if adm == "s4" or rname not in ("lnb_b2", "conf")
] + [
    # Stage W: the same on traces whitened per tetrode (4x4 ZCA from spike-excised noise of the first
    # 300 s), where stage A put the peeled-cosine gate at 0.5
    Arm(f"W_s45_c{int(c * 100)}_{rname}", amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", c), refractory_ms=1.5,
        recipe=RECIPES[rname], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True, amp_ref="median",
        whiten=True)
    for c, rname in ((0.5, "lnb"), (0.5, "lnb_b5"), (0.5, "lnb_b2"), (0.5, "june"), (0.5, "frozen"), (0.55, "lnb"),
                     (0.5, "conf"), (0.5, "conf_b5"))
] + [
    Arm("W_joint_masked_june", per_tetrode=False, masked=True, whiten=True),
] + [
    # The 48 h of the original binary, unmasked, with June's re-seeding every 12 windows from
    # mp_baseline's MountainSort5 re-seed sorts (run with --bank <BAKE_ROOT>/bank0.npz --windows 0-95)
    Arm("R48_W_s45_c50_lnb_b2", hybrid=False, amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.5), refractory_ms=1.5,
        recipe=RECIPES["lnb_b2"], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True, amp_ref="median",
        whiten=True, reseed_every=12, reseed_sorts=str(bk.SORTINGS / "mp_baseline" / "reseed_sorts")),
    Arm("R48_W_joint_masked_june", per_tetrode=False, hybrid=False, masked=True, whiten=True, reseed_every=12,
        reseed_sorts=str(bk.SORTINGS / "mp_baseline_masked" / "reseed_sorts")),
] + [
    # W_s45_c50_lnb_b2 again, after the driven test moved to the relaxed rule
    Arm("W_s45_c50_lnb_b2_dt2", amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.5), refractory_ms=1.5,
        recipe=RECIPES["lnb_b2"], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True,
        amp_ref="median", whiten=True),
    # joint matching with the per-tetrode arm's floors, admission and recipe: matcher against recipe
    Arm("W_joint_masked_s45_c50_lnb_b2", per_tetrode=False, masked=True, amp_rel=0.5, snr_floor=4.5,
        gate=("cos_pl", 0.5), refractory_ms=1.5, recipe=RECIPES["lnb_b2"], match_floor=(0.2, 2.5),
        pulse_rule=(2.5, 0.3, 0.5), drop_bad=True, amp_ref="median", whiten=True),
    Arm("R48_W_joint_masked_s45_c50_lnb_b2", per_tetrode=False, hybrid=False, masked=True, amp_rel=0.5,
        snr_floor=4.5, gate=("cos_pl", 0.5), refractory_ms=1.5, recipe=RECIPES["lnb_b2"], match_floor=(0.2, 2.5),
        pulse_rule=(2.5, 0.3, 0.5), drop_bad=True, amp_ref="median", whiten=True, reseed_every=12,
        reseed_sorts=str(bk.SORTINGS / "mp_baseline_masked" / "reseed_sorts")),
    # A whitening matrix per window instead of the first 300 s'. The hybrid arm runs windows 0-11; the
    # real one runs windows 48-59 from R48's bank after window 47 (--bank of that bank unwhitened).
    Arm("W_s45_c50_lnb_b2_rw", amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.5), refractory_ms=1.5,
        recipe=RECIPES["lnb_b2"], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True,
        amp_ref="median", whiten=True, rewhiten=True),
    Arm("R48_W_s45_c50_lnb_b2_rw", hybrid=False, amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.5),
        refractory_ms=1.5, recipe=RECIPES["lnb_b2"], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5),
        drop_bad=True, amp_ref="median", whiten=True, rewhiten=True, reseed_every=12,
        reseed_sorts=str(bk.SORTINGS / "mp_baseline" / "reseed_sorts")),
] + [
    # re-centring the rebuilt templates
    Arm("W_s45_c50_lnb_b2_rc", amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.5), refractory_ms=1.5,
        recipe=RECIPES["lnb_b2_rc"], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5), drop_bad=True,
        amp_ref="median", whiten=True),
    Arm("W_joint_masked_june_rc", per_tetrode=False, masked=True, whiten=True, recipe=RECIPES["june_rc"]),
    Arm("R48_W_s45_c50_lnb_b2_rc", hybrid=False, amp_rel=0.5, snr_floor=4.5, gate=("cos_pl", 0.5),
        refractory_ms=1.5, recipe=RECIPES["lnb_b2_rc"], match_floor=(0.2, 2.5), pulse_rule=(2.5, 0.3, 0.5),
        drop_bad=True, amp_ref="median", whiten=True, reseed_every=12,
        reseed_sorts=str(bk.SORTINGS / "mp_baseline" / "reseed_sorts")),
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
    ap.add_argument("--bank", default=None, help="starting bank (default <BAKE>/bank_start.npz)")
    ap.add_argument("--root", default=None, help="output root (default <BAKE>/arms)")
    a = ap.parse_args()
    if a.make_bank:
        make_bank()
        return
    lo, hi = (int(v) for v in a.windows.split("-"))
    arm = ARMS[a.name]
    bank = bk.Bank.load(a.bank or bk.BAKE / "bank_start.npz")
    if arm.whiten:
        groups = np.asarray(si.load(bk.HYBRID if arm.hybrid else bk.BINARY).get_property("group")).astype(int)
        bank = bk.whiten_bank(bank, bk.whitening_matrix(arm.hybrid), groups)
    bk.log(f"{arm.name}: {arm}")
    root = pathlib.Path(a.root) if a.root else bk.BAKE / "arms"
    bk.run_carry(arm, bank, root / arm.name, windows=range(lo, hi + 1), n_jobs=a.n_jobs,
                 workers=a.workers)


if __name__ == "__main__":
    main()
