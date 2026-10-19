#!/usr/bin/env python3
"""
Pin the pose reconstruction against LightDock's own lgd_generate_conformations.py.

Nothing about the quaternion convention is assumed.  A gso_<step>.out from a
live run is copied into a scratch directory, lgd_generate_conformations.py is
asked to write real PDB poses there, and each candidate convention is scored by
max per-atom deviation from those PDBs.  A convention only passes at <= 1e-2 A.

Run:  s0_calibrate.py <run_dir> <gso_file> <n_glowworms> <scratch_dir>
"""
import json
import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import (QUAT_CONVENTIONS, parse_pdb_atoms, poses, read_gso,
                    receptor_shift, xyz_array)

run_dir = os.path.abspath(sys.argv[1])
gso_file = os.path.abspath(sys.argv[2])
n_glow = int(sys.argv[3])
# absolute: lgd_generate_conformations runs with cwd=scratch, so a relative
# scratch path would be resolved twice.
scratch = os.path.abspath(sys.argv[4])

os.makedirs(scratch, exist_ok=True)
for f in ("lightdock_receptor.pdb", "lightdock_ligand.pdb",
          "receptor.pdb", "ligand.pdb", "setup.json"):
    shutil.copy(os.path.join(run_dir, f), os.path.join(scratch, f))
local_gso = os.path.join(scratch, os.path.basename(gso_file))
shutil.copy(gso_file, local_gso)

# get_lightdock_structures() prepends "lightdock_" to the basename it is given,
# so hand it the plain names.
cmd = ["/usr/local/bin/lgd_generate_conformations.py",
       os.path.join(scratch, "receptor.pdb"),
       os.path.join(scratch, "ligand.pdb"),
       local_gso, str(n_glow),
       "--setup", os.path.join(scratch, "setup.json")]
proc = subprocess.run(cmd, capture_output=True, text=True, cwd=scratch)
print("lgd_generate_conformations rc =", proc.returncode)
if proc.returncode != 0:
    print("STDERR:", proc.stderr[-3000:])
    sys.exit(1)

lig = parse_pdb_atoms(os.path.join(run_dir, "lightdock_ligand.pdb"))
lig_xyz = xyz_array(lig)
n_lig = len(lig)
lig_names = [(a["chain"], a["resid"], a["name"]) for a in lig]

trans, quat, score = read_gso(local_gso)
print(f"gso lines = {len(trans)}, ligand atoms = {n_lig}")

# The generated lightdock_<i>.pdb holds receptor then ligand; take the tail.
truth = []
used = []
for i in range(n_glow):
    p = os.path.join(scratch, f"lightdock_{i}.pdb")
    if not os.path.exists(p):
        continue
    atoms = parse_pdb_atoms(p)
    tail = atoms[-n_lig:]
    if [(a["chain"], a["resid"], a["name"]) for a in tail] != lig_names:
        raise RuntimeError(f"{p}: ligand block tail does not match lightdock_ligand.pdb")
    truth.append(xyz_array(tail))
    used.append(i)
truth = np.asarray(truth)
print(f"reference poses read = {len(truth)} (glowworms {used[:5]}...)")

shift, n_shift_atoms, spread = receptor_shift(run_dir)
print(f"receptor shift = {np.round(shift, 4)}  (from {n_shift_atoms} atoms, "
      f"component spread {spread:.2e} A)")

results = {}
for conv in QUAT_CONVENTIONS:
    P = poses(lig_xyz, trans[used], quat[used], conv, shift=None)
    dev = np.abs(P - truth)
    results[conv["name"]] = {
        "max_atom_dev_A": float(np.sqrt((( P - truth) ** 2).sum(-1)).max()),
        "max_coord_dev_A": float(dev.max()),
        "mean_atom_dev_A": float(np.sqrt(((P - truth) ** 2).sum(-1)).mean()),
    }
    print(f"  {conv['name']:>8}  max per-atom dev {results[conv['name']]['max_atom_dev_A']:12.6f} A")

winners = [k for k, v in results.items() if v["max_atom_dev_A"] <= 1e-2]
print("PASS:", winners)
if len(winners) != 1:
    print("FAIL: expected exactly one convention to reproduce lgd output")
    json.dump({"results": results, "winner": None}, open(os.path.join(scratch, "calibration.json"), "w"), indent=2)
    sys.exit(2)

winner = winners[0]
conv = [c for c in QUAT_CONVENTIONS if c["name"] == winner][0]

# Independent sanity check: LightDock poses are receptor-centred, so adding the
# receptor shift must land the pose in the deposited-PDB frame.  Verify by
# confirming the shifted pose sits inside the deposited receptor's bounding box.
rec_raw = xyz_array(parse_pdb_atoms(os.path.join(run_dir, "receptor.pdb")))
P_dep = poses(lig_xyz, trans[used], quat[used], conv, shift=shift)
lo, hi = rec_raw.min(0) - 40, rec_raw.max(0) + 40
inside = int(np.all((P_dep.reshape(-1, 3) >= lo) & (P_dep.reshape(-1, 3) <= hi), axis=1).sum())
print(f"frame check: {inside}/{P_dep.reshape(-1,3).shape[0]} pose atoms within "
      f"receptor bbox+40 A in deposited frame")

json.dump({"results": results, "winner": winner,
           "receptor_shift": shift.tolist(),
           "receptor_shift_spread_A": spread,
           "n_reference_poses": len(truth),
           "frame_check_inside": inside,
           "frame_check_total": int(P_dep.reshape(-1, 3).shape[0])},
          open(os.path.join(scratch, "calibration.json"), "w"), indent=2)
print("WINNER =", winner)
