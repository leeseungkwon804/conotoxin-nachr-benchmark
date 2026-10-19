#!/bin/bash
# Emits a line only when the scored-run count changes, or on a terminal condition.
# Must run under bash: zsh does not word-split an unquoted "$SSH" expansion, which
# silently broke an earlier version. Each probe uses $(... | wc -l) so it emits
# exactly one line -- `pgrep -c` both prints and exits non-zero, which shifted
# every field when combined with `|| echo 0`.
AN=/workspace/conotoxin_bench/s0_analysis
Q=/workspace/conotoxin_bench/s0_gate/queue_v3_A1.log
prev=-1
while :; do
  st=$(ssh -i /Users/dkyoo/.ssh/id_ed25519 -o StrictHostKeyChecking=no -o BatchMode=yes \
        -o ConnectTimeout=20 -p 14087 root@213.173.107.137 \
        "ls $AN/done 2>/dev/null | wc -l; grep -c 'chain A1 DONE' $Q 2>/dev/null || true; pgrep -f 'lightdock3.py setup.json 100' | wc -l; pgrep -f s0_watch.sh | wc -l" 2>/dev/null)
  if [ -z "$st" ]; then sleep 120; continue; fi
  n=$(sed -n 1p <<<"$st"); a1=$(sed -n 2p <<<"$st")
  ld=$(sed -n 3p <<<"$st"); w=$(sed -n 4p <<<"$st")
  if [ "${n:-0}" != "$prev" ]; then
    echo "$(date -u +%H:%M:%SZ) S0 scored ${n}/12   A1-done=${a1}  lightdock-procs=${ld}  watcher=${w}"
    prev=$n
  fi
  if [ "${n:-0}" -ge 12 ]; then
    echo "ALL 12 RUNS SCORED - finalise ladder and report"; exit 0
  fi
  if [ "${a1:-0}" -ge 1 ] && [ "${ld:-0}" -eq 0 ]; then
    echo "$(date -u +%H:%M:%SZ) A1 chain finished; ${n}/12 scored, waiting on watcher (watcher procs=${w})"
  fi
  if [ "${ld:-0}" -eq 0 ] && [ "${a1:-0}" -eq 0 ]; then
    echo "ALERT no lightdock process and no A1 DONE marker - chain died at ${n}/12"; exit 3
  fi
  if [ "${w:-0}" -eq 0 ]; then
    echo "ALERT watcher process gone at ${n}/12 - completed runs will not be scored"; exit 4
  fi
  sleep 300
done
