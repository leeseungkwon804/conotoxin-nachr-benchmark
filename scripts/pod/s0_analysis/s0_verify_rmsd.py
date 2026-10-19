#!/usr/bin/env python3
"""
Independent check of the analytic RMSD, through real PDB files.

For a chosen swarm, lgd_generate_conformations.py writes actual pose PDBs.  The
RMSD is then recomputed from those files by a code path that shares nothing with
s0_rmsd_run.py except the receptor shift: atoms are matched by (residue-order,
atom-name) read straight off the generated PDB text.  The two numbers must agree
to PDB write precision.

Usage: s0_verify_rmsd.py <run_dir> <native_ligand_pdb> <swarm_id> <step> <scratch_dir>
"""
import os
import shutil
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import (QUAT_CONVENTIONS, build_correspondence, parse_pdb_atoms,
                    read_gso, receptor_shift, residue_blocks, rmsd_to_native,
                    xyz_array)

run_dir = os.path.abspath(sys.argv[1])
native_pdb = os.path.abspath(sys.argv[2])
swarm_id = int(sys.argv[3])
step = int(sys.argv[4])
scratch = os.path.abspath(sys.argv[5])
conv = [c for c in QUAT_CONVENTIONS if c["name"] == "wxyz"][0]

os.makedirs(scratch, exist_ok=True)
for f in ("lightdock_receptor.pdb", "lightdock_ligand.pdb",
          "receptor.pdb", "ligand.pdb", "setup.json"):
    shutil.copy(os.path.join(run_dir, f), os.path.join(scratch, f))
gso = os.path.join(run_dir, f"swarm_{swarm_id}", f"gso_{step}.out")
local_gso = os.path.join(scratch, f"gso_{step}.out")
shutil.copy(gso, local_gso)

trans, quat, score = read_gso(local_gso)
n_glow = len(trans)

# ---- path A: analytic
lig_all = parse_pdb_atoms(os.path.join(run_dir, "lightdock_ligand.pdb"))
pose_idx, native_xyz, corr = build_correspondence(
    os.path.join(run_dir, "lightdock_ligand.pdb"), native_pdb)
shift, _, _ = receptor_shift(run_dir)
rms_analytic = rmsd_to_native(xyz_array(lig_all)[pose_idx], trans, quat,
                              conv, shift, native_xyz)

# ---- path B: through real PDB files
proc = subprocess.run(
    ["/usr/local/bin/lgd_generate_conformations.py",
     os.path.join(scratch, "receptor.pdb"), os.path.join(scratch, "ligand.pdb"),
     local_gso, str(n_glow), "--setup", os.path.join(scratch, "setup.json")],
    capture_output=True, text=True, cwd=scratch)
if proc.returncode != 0:
    print("lgd_generate_conformations failed:", proc.stderr[-2000:])
    sys.exit(1)

nat = parse_pdb_atoms(native_pdb)
n_order, n_blocks = residue_blocks(nat)
nat_seq = [n_blocks[k]["resname"] for k in n_order]

n_lig = len(lig_all)
rms_pdb = np.full(n_glow, np.nan)
for i in range(n_glow):
    p = os.path.join(scratch, f"lightdock_{i}.pdb")
    if not os.path.exists(p):
        continue
    tail = parse_pdb_atoms(p)[-n_lig:]
    d_order, d_blocks = residue_blocks(tail)
    dock_seq = [d_blocks[k]["resname"] for k in d_order]
    # native may carry one extra trailing residue (an amidation cap); pair by
    # order over the common prefix and require the residue names to agree.
    pairs = [(a, b) for a, b in zip(range(len(d_order)), range(len(n_order)))
             if dock_seq[a] == nat_seq[b]]
    a_xyz, b_xyz = [], []
    for di, ni in pairs:
        dm = {tail[k]["name"]: tail[k]["xyz"] for k in d_blocks[d_order[di]]["idx"]}
        nm = {nat[k]["name"]: nat[k]["xyz"] for k in n_blocks[n_order[ni]]["idx"]}
        for name in sorted(set(dm) & set(nm)):
            a_xyz.append(dm[name])
            b_xyz.append(nm[name])
    # the generated pose is receptor-centred; put it back in the deposited frame
    A = np.asarray(a_xyz) + shift
    B = np.asarray(b_xyz)
    rms_pdb[i] = float(np.sqrt(((A - B) ** 2).sum(1).mean()))
    if i == 0:
        print(f"path B atom pairs: {len(a_xyz)}  (path A: {corr['matched_atoms']})")

ok = ~np.isnan(rms_pdb)
diff = np.abs(rms_analytic[ok] - rms_pdb[ok])
print(f"swarm {swarm_id}: compared {ok.sum()} poses")
print(f"  analytic best {rms_analytic[ok].min():.4f} A   pdb-path best {rms_pdb[ok].min():.4f} A")
print(f"  max |analytic - pdb| = {diff.max():.6f} A   mean {diff.mean():.6f} A")
print("AGREE" if diff.max() < 0.01 else "DISAGREE")
