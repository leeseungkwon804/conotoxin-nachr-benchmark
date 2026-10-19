"""Independently recompute native 5.0 A contact sets and their overlap with each
restraint set actually used, so the dose-response table cites measured coverage."""
import glob, json, os, sys
import numpy as np
sys.path.insert(0, "/workspace/conotoxin_bench/s0_analysis/bin")
from s0_lib import parse_pdb_atoms

LI = "/workspace/conotoxin_bench/structures/lightdock_inputs"
CUT = 5.0

def native_set(pair_pdb, lig_pdb):
    rec = parse_pdb_atoms(pair_pdb); lig = parse_pdb_atoms(lig_pdb)
    L = np.asarray([a["xyz"] for a in lig])
    R = np.asarray([a["xyz"] for a in rec])
    d = np.sqrt(((L[:, None, :] - R[None, :, :]) ** 2).sum(-1))
    hit = (d < CUT).any(0)
    return {(rec[i]["chain"], rec[i]["resname"], rec[i]["resid"].strip())
            for i in np.where(hit)[0]}

def restraint_set(path):
    out = set()
    for line in open(path):
        if line.startswith("R "):
            c, rn, num = line.split()[1].split(".")
            out.add((c.upper(), rn[:3].upper(), num))
    return out

for tag, pair, lig, rfiles in [
    ("A1 9NX0 (pair D,E)", f"{LI}/9NX0_receptorPAIR_DE.pdb", f"{LI}/9NX0_ligand_I.pdb",
     ["/workspace/conotoxin_bench/s0_gate/A1_9NX0_XDON_rep1/restraints.list",
      "/workspace/conotoxin_bench/s0_gate/A1_9NX0_NATIVE_rep1/restraints.list"]),
    ("A2 7EGR (pair H,G)", f"{LI}/7EGR_receptorPAIR_HG.pdb", f"{LI}/7EGR_ligand_P.pdb",
     ["/workspace/conotoxin_bench/s0_gate/A2_7EGR_XDON_rep1/restraints.list",
      "/workspace/conotoxin_bench/s0_gate/A2_7EGR_NATIVE_rep1/restraints.list"])]:
    nat = native_set(pair, lig)
    print(f"{tag}: native 5.0 A contact residues = {len(nat)}")
    for rf in rfiles:
        if not os.path.exists(rf):
            print(f"    {rf}: missing"); continue
        rs = restraint_set(rf)
        ov = rs & nat
        rung = rf.split("/")[-2]
        print(f"    {rung:<24} R={len(rs):>2}  overlap={len(ov):>2}  "
              f"recall(of native)={len(ov)/len(nat):.3f}  "
              f"precision(of restraints)={len(ov)/len(rs):.3f}")
