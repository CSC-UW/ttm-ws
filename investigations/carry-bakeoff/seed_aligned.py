"""The seed bank re-sorted on the aligned-reference binary, as mp_baseline built its own.

MountainSort5 on the first 1800 s of ``mp_aligned/binary_w00_11`` (``ttm.tracking.sort_chunk``, as
``subjects/ttm-003/scripts/mp_baseline.py`` sorted window 0 of ``mp_baseline/binary``), then units
with snr >= 5 and >= 100 spikes, deduplicated at cosine 0.95, templates from ``_mp_common``, cut to
each unit's tetrode and 30 samples before to 60 after the spike. Writes
``v2_aligned/bank0_aligned.npz`` and the seed sorting beside the binary (``mp_aligned/seed_sort``,
``mp_aligned/seed_bank_sorting``), and prints the units per tetrode against ``bank0.npz``.

Run: BAKE_REFERENCE=aligned .venv/bin/python investigations/carry-bakeoff/seed_aligned.py
"""
from __future__ import annotations

import pathlib
import shutil
import sys

import numpy as np
import spikeinterface as si

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "spike-sorting"))
import _bake as bk  # noqa: E402
from _mp_common import build_templates_object, dedup_sorting  # noqa: E402
from ttm import capture as cp  # noqa: E402
from ttm.tracking import sort_chunk, to_int_numpy_sorting  # noqa: E402

OUT = bk.SORTINGS / "mp_aligned"


def main():
    assert bk.REFERENCE == "aligned", "run with BAKE_REFERENCE=aligned"
    rec = si.load(bk.BINARY)
    win0 = rec.frame_slice(0, bk.WIN)
    win0.reset_times()
    shutil.rmtree(OUT / "seed_sort", ignore_errors=True)
    ref0 = to_int_numpy_sorting(sort_chunk(win0, OUT / "seed_sort"))
    _, az = build_templates_object(ref0, win0, with_snr=True, n_jobs=32)
    snr = az.get_extension("quality_metrics").get_data()["snr"].to_numpy()
    nsp = np.array([len(ref0.get_unit_spike_train(u)) for u in ref0.unit_ids])
    conf = np.asarray(ref0.unit_ids)[(snr >= 5) & (nsp >= 100)]
    seed = dedup_sorting(ref0.select_units(conf), win0, cosine_min=0.95)
    tm, _ = build_templates_object(seed, win0, with_snr=False, n_jobs=32)
    print(f"{ref0.get_num_units()} seed units -> {len(conf)} confident -> {seed.get_num_units()} after dedup")
    seed.save(folder=OUT / "seed_bank_sorting", overwrite=True)
    dense, nb = tm.get_dense_templates(), tm.nbefore
    groups = np.asarray(seed.get_property("group")).astype(int)
    T = np.stack([dense[i, nb - bk.NB:nb + bk.NA][:, cp.tetrode_channel_indices(rec, g)] for i, g in enumerate(groups)])
    bank = bk.Bank(np.arange(len(groups), dtype=np.int64), groups, T.astype(np.float32))
    bank.save(bk.BAKE / "bank0_aligned.npz")
    old = bk.Bank.load(bk.BAKE_ROOT / "bank0.npz")
    print("units per tetrode, global-median seed bank / aligned:")
    print("  " + "  ".join(f"TT{g + 1} {(old.group == g).sum()}/{(groups == g).sum()}" for g in range(16)))


if __name__ == "__main__":
    main()
