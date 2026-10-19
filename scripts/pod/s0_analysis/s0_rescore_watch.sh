#!/bin/bash
# Fills in restraint-free rescoring for any run that has been scored but not yet
# rescored. Independent of s0_watch.sh so neither has to be restarted; both are
# idempotent and skip work already on disk.
set -u
BENCH=/workspace/conotoxin_bench
AN=$BENCH/s0_analysis
LI=$BENCH/structures/lightdock_inputs
LOG=$AN/rescore_watch.log
mkdir -p "$AN/rescore"
echo "=== rescore watcher START $(date -u +%FT%TZ) pid=$$" >> "$LOG"
while :; do
  pending=0
  for j in "$AN"/rmsd/*.json; do
    [ -f "$j" ] || continue
    r=$(basename "$j" .json)
    [ -f "$AN/rescore/$r.npz" ] && continue
    # only score complete runs
    python3 -c "import json,sys; sys.exit(0 if json.load(open('$j')).get('complete') else 1)" || { pending=1; continue; }
    case "$r" in
      A1_*) nat=$LI/9NX0_ligand_I.pdb ;;
      A2_*) nat=$LI/7EGR_ligand_P.pdb ;;
      *) continue ;;
    esac
    echo "--- rescoring $r $(date -u +%FT%TZ)" >> "$LOG"
    nice -n 19 python3 "$AN/bin/s0_rescore.py" "$BENCH/s0_gate/$r" "$nat" \
        "$AN/rmsd/$r.npz" 100 "$AN/rescore/$r.npz" >> "$LOG" 2>&1 \
      || echo "!!! rescore FAILED for $r" >> "$LOG"
  done
  n=$(ls "$AN"/rescore/*.npz 2>/dev/null | wc -l)
  [ "$n" -ge 12 ] && { echo "=== rescore watcher DONE 12/12 $(date -u +%FT%TZ)" >> "$LOG"; break; }
  sleep 300
done
