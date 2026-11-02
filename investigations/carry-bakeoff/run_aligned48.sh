#!/usr/bin/env bash
# The aligned-reference 48 h: binary, whitening, re-seed sorts, then both carries and their exports.
# The per-tetrode carry runs windows 0-11 while the re-seed sorts run, and continues once they exist.
# Run from the workspace root: nohup bash investigations/carry-bakeoff/run_aligned48.sh > LOG 2>&1 &
set -euo pipefail
export BAKE_REFERENCE=aligned48
PY=.venv/bin/python
CB=investigations/carry-bakeoff
S=/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix
ROOT=$S/carry_bakeoff/aligned48
BANK=$S/carry_bakeoff/v2_aligned/bank0_aligned.npz
mkdir -p "$ROOT/logs"
stamp() { date '+%F %T'; }

if [ ! -f "$S/mp_aligned/binary/binary.json" ]; then
  echo "$(stamp) binary"
  $PY investigations/common-reference/aligned_binary.py --binary --full --jobs 64
fi
echo "$(stamp) whitening"
$PY -c "import sys; sys.path.insert(0, '$CB'); import _bake as bk; bk.BAKE.mkdir(parents=True, exist_ok=True); bk.whitening_matrix(False)"

echo "$(stamp) re-seed sorts and per-tetrode windows 0-11"
$PY $CB/reseed_aligned.py --workers 14 > "$ROOT/logs/reseed_aligned.log" 2>&1 &
P_SORT=$!
$PY $CB/run_arm.py A48_W_s45_c50_lnb_b2 --bank "$BANK" --windows 0-11 --root "$ROOT" \
  > "$ROOT/logs/A48_W_s45_c50_lnb_b2_w00_11.log" 2>&1 &
P_PT=$!
wait $P_SORT
echo "$(stamp) re-seed sorts done"
$PY $CB/run_arm.py A48_W_joint_masked_s45_c50_lnb_b2 --bank "$BANK" --windows 0-95 --root "$ROOT" \
  > "$ROOT/logs/A48_W_joint_masked_s45_c50_lnb_b2.log" 2>&1 &
P_JT=$!
wait $P_PT
echo "$(stamp) per-tetrode windows 0-11 done"
$PY $CB/run_arm.py A48_W_s45_c50_lnb_b2 --bank "$BANK" --windows 0-95 --root "$ROOT" \
  > "$ROOT/logs/A48_W_s45_c50_lnb_b2.log" 2>&1
echo "$(stamp) per-tetrode done"
wait $P_JT
echo "$(stamp) joint done"
$PY $CB/export_sorting.py A48_W_s45_c50_lnb_b2 --root "$ROOT" --run-dir "$S/carry_A48_pt" --name assembled_reseed_c12
$PY $CB/export_sorting.py A48_W_joint_masked_s45_c50_lnb_b2 --root "$ROOT" --run-dir "$S/carry_A48_joint_rules" \
  --name assembled_reseed_c12
echo "$(stamp) exported"
