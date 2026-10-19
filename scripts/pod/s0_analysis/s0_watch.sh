#!/bin/bash
# Watch the two S0 chains and score each run's ligand-RMSD the moment it finishes.
#
# Completion is taken from the chain queue logs ("--- <run> rc=<n> END"), not from
# counting gso files, so a half-written gso_100.out can never be scored.
# One nice'd python process at a time; the docking itself owns 40 of 48 cores.
set -u
BENCH=/workspace/conotoxin_bench
S0=$BENCH/s0_gate
AN=$BENCH/s0_analysis
BIN=$AN/bin
LI=$BENCH/structures/lightdock_inputs
RMSD=$AN/rmsd
LOG=$AN/watch.log
mkdir -p "$RMSD" "$AN/done"

native_for () {
  case "$1" in
    A1_*) echo "$LI/9NX0_ligand_I.pdb" ;;
    A2_*) echo "$LI/7EGR_ligand_P.pdb" ;;
    *) echo "" ;;
  esac
}

RUNS="A1_9NX0_BLIND_rep1 A1_9NX0_BLIND_rep2 A1_9NX0_XDON_rep1 A1_9NX0_XDON_rep2 \
A1_9NX0_NATIVE_rep1 A1_9NX0_NATIVE_rep2 A2_7EGR_BLIND_rep1 A2_7EGR_BLIND_rep2 \
A2_7EGR_XDON_rep1 A2_7EGR_XDON_rep2 A2_7EGR_NATIVE_rep1 A2_7EGR_NATIVE_rep2"

echo "=== watcher START $(date -u +%FT%TZ) pid=$$" >> "$LOG"
while :; do
  newly=0
  for r in $RUNS; do
    [ -f "$AN/done/$r" ] && continue
    # chain queue log records the run's exit
    if ! grep -qE "^--- $r rc=[0-9]+ END" "$S0"/queue_v3_A1.log "$S0"/queue_v3_A2.log 2>/dev/null; then
      continue
    fi
    rc=$(grep -hoE "^--- $r rc=[0-9]+" "$S0"/queue_v3_A*.log 2>/dev/null | tail -1 | sed 's/.*rc=//')
    nat=$(native_for "$r")
    echo "--- $r finished rc=$rc, scoring $(date -u +%FT%TZ)" >> "$LOG"
    if nice -n 19 python3 "$BIN/s0_rmsd_run.py" "$S0/$r" "$nat" 100 "$RMSD" >> "$LOG" 2>&1; then
      echo "$rc" > "$AN/done/$r"
      newly=1
    else
      echo "!!! $r scoring FAILED (will retry next cycle)" >> "$LOG"
    fi
  done
  n=$(ls "$AN/done" 2>/dev/null | wc -l)
  if [ "$newly" = "1" ]; then
    # s0_ladder.py also emits the score-vs-RMSD table; it is re-executed fresh
    # each cycle, so improvements to it land without restarting this watcher.
    nice -n 19 python3 "$BIN/s0_ladder.py" "$RMSD" "$AN/ladder" >> "$LOG" 2>&1 \
      || echo "!!! ladder failed" >> "$LOG"
    echo "=== scored $n/12 $(date -u +%FT%TZ)" >> "$LOG"
  fi
  [ "$n" -ge 12 ] && { echo "=== watcher DONE all 12 $(date -u +%FT%TZ)" >> "$LOG"; break; }
  sleep 300
done
