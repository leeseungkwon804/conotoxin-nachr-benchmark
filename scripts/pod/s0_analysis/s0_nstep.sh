#!/bin/bash
# Decisive one-step experiment.
#
# Builds a 2-swarm copy of an existing S0 run and runs lightdock3.py for exactly
# ONE step.  With one step the GSO does:
#     update_luciferin()   -> evaluates the objective at P_0 (the initial positions)
#     movement_phase()     -> moves P_0 -> P_1
#     save(1)              -> writes gso_1.out
# P_0 is known exactly from init/initial_positions_0.dat, so we can ask which
# scoring variant, evaluated at WHICH position, reproduces the Scoring column
# in gso_1.out to floating point.  That settles both questions at once:
#   Q1 does f_R (receptor restraints) enter the score?
#   Q2 is the Scoring column paired with P_0 (pre-move) or P_1 (post-move)?
set -eu
SRC="$1"; STEPS="${3:-1}"
DST="$2"                       # scratch dir
rm -rf "$DST"; mkdir -p "$DST/init"
for f in receptor.pdb ligand.pdb lightdock_receptor.pdb lightdock_ligand.pdb \
         lightdock_receptor_mask.npy lightdock_ligand_mask.npy setup.json; do
  [ -f "$SRC/$f" ] && cp "$SRC/$f" "$DST/$f"
done
[ -f "$SRC/restraints.list" ] && cp "$SRC/restraints.list" "$DST/restraints.list"
cp "$SRC/init/initial_positions_0.dat" "$DST/init/"
cp "$SRC/init/initial_positions_1.dat" "$DST/init/"
[ -f "$SRC/init/swarm_centers.pdb" ] && cp "$SRC/init/swarm_centers.pdb" "$DST/init/"
python3 - "$DST/setup.json" <<'EOF'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
d["swarms"] = 2
json.dump(d, open(p, "w"), indent=4)
print("  setup.json swarms ->", d["swarms"], "| restraints:", d["restraints"])
EOF
cd "$DST"
# lightdock3.py writes into swarm_<i>/ but does not create them -- that is
# lightdock3_setup.py's job, and we deliberately did not re-run setup.
mkdir -p swarm_0 swarm_1
nice -n 19 /usr/local/bin/lightdock3.py setup.json ${STEPS:-1} -c 1 -s vdw -sg 324324 > sim.log 2>&1
echo "  lightdock3 rc=$? ; wrote: $(ls swarm_0/ 2>/dev/null | tr '\n' ' ')"
grep -iE "error|traceback" sim.log | head -3 || true
