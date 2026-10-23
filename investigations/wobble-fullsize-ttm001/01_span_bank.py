"""Rebuild June's carry-forward bank at the start of the span and June's re-seed units at its boundary.

The bank entering window 55 is each unit's template from the last window up to
54 in which June's ``assembled_reseed_c12`` gave it at least 100 spikes,
computed as June's loop computed it. The re-seed units are June's rule
(snr >= 5, >= 100 spikes, shift-cosine < 0.8 to the live same-tetrode bank)
applied to June's MountainSort5 re-sort of window 60 against the bank after
window 60's re-estimation, rebuilt the same way. The check is that the rule
gives June's births at window 60 (ids 80 to 89) with June's at-birth cosines
(``reseed_cadence_sweep.npz``).

It also writes ``template_amp_w55.npy``: per unit (id, tetrode, trough of the
template on its deepest wire in raw units, the same in excised-MAD units of
the first 5 min of window 55, template energy), the source of the amplitude
classes (``_span.amp_classes``). ``--amps-only`` rewrites only that file.

    uv run python investigations/wobble-fullsize-ttm001/01_span_bank.py [--amps-only]
"""
import argparse
import json
import pickle

import numpy as np
import spikeinterface as si

from _span import JUNE, OUT, RESEED_W, SPAN, WLEN, june_reseed_units, log, reconstruct_bank, recording
from ttm import capture as cp


def template_amplitudes(rec):
    with open(OUT / f"bank_w{SPAN[0]}.pkl", "rb") as f:
        b = pickle.load(f)
    with open(OUT / f"reseed_w{RESEED_W}.pkl", "rb") as f:
        added = pickle.load(f)
    groups = np.asarray(rec.get_property("group")).astype(int)
    a = SPAN[0] * WLEN
    noise = {g: cp.noise_mad_excised(cp.load_tetrode_traces(rec, g, a, a + 300 * 30000))[0] for g in range(16)}
    bank, ug = dict(b["bank"]), dict(b["group"])
    for d in added:
        bank[d["uid"]], ug[d["uid"]] = d["dense"], d["group"]
    rows = []
    for u in sorted(bank):
        g = ug[u]
        t = bank[u][:, groups == g]
        tr = -t.min(axis=0)
        c = int(np.argmax(tr))
        rows.append((u, g, tr[c], tr[c] / noise[g][c], float((t.astype(float) ** 2).sum())))
    np.save(OUT / "template_amp_w55.npy", np.array(rows))
    log(f"template amplitudes for {len(rows)} units -> template_amp_w55.npy")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--amps-only", action="store_true")
    OUT.mkdir(parents=True, exist_ok=True)
    rec = recording()
    if ap.parse_args().amps_only:
        template_amplitudes(rec)
        return
    asm = si.load(JUNE / "assembled_reseed_c12")
    log(f"June assembly: {asm.get_num_units()} units")
    bank, ug, src = reconstruct_bank(rec, asm, SPAN[0] - 1)
    log(f"bank entering window {SPAN[0]}: {len(bank)} units; source windows {sorted(set(src.values()))}")
    with open(OUT / f"bank_w{SPAN[0]}.pkl", "wb") as f:
        pickle.dump({"bank": bank, "group": ug, "source_window": src}, f)

    bank60, ug60, _ = reconstruct_bank(rec, asm, RESEED_W)
    added = june_reseed_units(rec, bank60, ug60, RESEED_W, next_id=max(bank60) + 1)
    sweep = np.load(JUNE / "reseed_cadence_sweep.npz")
    jb = {int(u): float(c) for u, c, w in zip(sweep["born_ids_c12"], sweep["birth_cos_c12"], sweep["birth_win_c12"])
          if int(w) == RESEED_W}
    rows = [{"uid": d["uid"], "group": d["group"], "ms5_unit": d["ms5_unit"], "snr": round(d["snr"], 2),
             "birth_cos": round(d["birth_cos"], 4), "june_birth_cos": round(jb.get(d["uid"], float("nan")), 4)}
            for d in added]
    for r in rows:
        log(f"re-seed unit {r}")
    ok = sorted(jb) == [d["uid"] for d in added] and all(
        abs(d["birth_cos"] - jb[d["uid"]]) < 1e-3 for d in added)
    log(f"re-seed at window {RESEED_W}: {len(added)} units; June had {len(jb)}; identical: {ok}")
    with open(OUT / f"reseed_w{RESEED_W}.pkl", "wb") as f:
        pickle.dump(added, f)
    (OUT / "bank_check.json").write_text(json.dumps(
        {"span": SPAN, "bank_units": len(bank), "reseed_window": RESEED_W, "reseed_units": rows,
         "june_births_at_window": {str(k): v for k, v in jb.items()}, "reseed_identical": ok}, indent=1))
    template_amplitudes(rec)


if __name__ == "__main__":
    main()
