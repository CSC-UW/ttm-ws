#!/usr/bin/env bash
# After run_aligned48.sh: the bad-period and coverage comparison of the two aligned 48 h carries,
# and the tagged-cell blocks of each, all scored on the aligned binary they were sorted from.
# Run from the workspace root: nohup bash investigations/carry-bakeoff/compare_aligned48.sh > LOG 2>&1 &
set -euo pipefail
PY=.venv/bin/python
S=/nvme/neuropixels/tetrode_data/2026-10-20_09-02-12/sortings_seed42_pcafix
L=$S/carry_bakeoff/aligned48/logs
EXP="--subject TTM-003 --experiment TTM-003-48H"
stamp() { date '+%F %T'; }
until [ -f "$S/mp_aligned/spike_coverage.npz" ]; do sleep 60; done
echo "$(stamp) bad periods and coverage"
$PY subjects/ttm-003/scripts/bad_periods_compare.py $EXP --a carry_A48_joint_rules --a-name joint_rules \
  --b carry_A48_pt --b-name per_tetrode --sorting assembled_reseed_c12 --binary-run mp_aligned \
  --coverage "$S/mp_aligned/spike_coverage.npz" --tag carry48_aligned_jr_pt > "$L/compare_aligned_jr_pt.log" 2>&1
for arm in pt joint_rules; do
  echo "$(stamp) tagged blocks $arm"
  $PY subjects/ttm-003/scripts/tagged_cells_blocks.py $EXP --sorting "$S/carry_A48_$arm/assembled_reseed_c12" \
    --name "carry48_aligned_$arm" --binary "$S/mp_aligned/binary" > "$L/tagged_aligned_$arm.log" 2>&1
done
echo "$(stamp) done"
