#!/usr/bin/env python3
"""
The hard floor a rigid-body docking of an APO ligand can ever reach.

S0 docks an apo NMR conopeptide as a rigid body.  If the apo conformer differs
internally from the bound conformer, then no rotation+translation can superimpose
them, and the best achievable ligand-RMSD is bounded below by the optimally
superposed (Kabsch) RMSD between the two conformers.  That number is a property
of the two structures alone -- no docking involved -- and it must be known before
a "< 3 A" gate can be interpreted: if the floor is already near 3 A, the gate is
testing the ligand's conformational change, not the docking.

Usage: s0_intrinsic_floor.py <run_dir> <native_ligand_pdb> [<run_dir> <native> ...]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import build_correspondence, parse_pdb_atoms, xyz_array


def kabsch_rmsd(P, Q):
    """RMSD after optimal rigid superposition of P onto Q."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Pr = Pc @ R
    return float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean())), R


out = []
args = sys.argv[1:]
for i in range(0, len(args), 2):
    run_dir, native_pdb = os.path.abspath(args[i]), args[i + 1]
    name = os.path.basename(run_dir)
    lig_pdb = os.path.join(run_dir, "lightdock_ligand.pdb")
    idx, nat_xyz, corr = build_correspondence(lig_pdb, native_pdb)
    dock_xyz = xyz_array(parse_pdb_atoms(lig_pdb))[idx]

    all_rmsd, _ = kabsch_rmsd(dock_xyz, nat_xyz)
    # backbone-only view: side chains of Arg/Trp are floppy and inflate the number
    names = []
    dl = parse_pdb_atoms(lig_pdb)
    for k in idx:
        names.append(dl[k]["name"])
    bb = np.array([n in ("N", "CA", "C", "O") for n in names])
    ca = np.array([n == "CA" for n in names])
    bb_rmsd, _ = kabsch_rmsd(dock_xyz[bb], nat_xyz[bb])
    ca_rmsd, _ = kabsch_rmsd(dock_xyz[ca], nat_xyz[ca])
    # per-atom deviation after all-atom superposition, to see where it lives
    Pc = dock_xyz - dock_xyz.mean(0)
    Qc = nat_xyz - nat_xyz.mean(0)
    _, R = kabsch_rmsd(dock_xyz, nat_xyz)
    dev = np.linalg.norm(Pc @ R - Qc, axis=1)
    worst = sorted(zip(dev, names, range(len(dev))), reverse=True)[:6]

    rec = {"run": name, "native": os.path.basename(native_pdb),
           "matched_atoms": int(len(idx)),
           "kabsch_rmsd_allatom_A": round(all_rmsd, 3),
           "kabsch_rmsd_backbone_A": round(bb_rmsd, 3),
           "kabsch_rmsd_CA_A": round(ca_rmsd, 3),
           "n_backbone": int(bb.sum()), "n_CA": int(ca.sum()),
           "max_atom_dev_after_superposition_A": round(float(dev.max()), 3),
           "worst_atoms": [(round(float(d), 2), n) for d, n, _ in worst]}
    out.append(rec)
    print(f"{name}  vs {rec['native']}")
    print(f"  matched atoms {rec['matched_atoms']}  "
          f"| Kabsch RMSD all-atom {all_rmsd:.3f} A, backbone(N,CA,C,O) {bb_rmsd:.3f} A "
          f"({int(bb.sum())} atoms), CA-only {ca_rmsd:.3f} A ({int(ca.sum())} atoms)")
    print(f"  worst atoms after superposition: {rec['worst_atoms']}")

dest = os.environ.get("FLOOR_OUT")
if dest:
    json.dump(out, open(dest, "w"), indent=2)
    print("wrote", dest)
