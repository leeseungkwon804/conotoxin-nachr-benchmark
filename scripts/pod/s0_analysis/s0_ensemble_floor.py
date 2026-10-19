#!/usr/bin/env python3
"""
Intrinsic rigid-docking floor across a whole apo NMR ensemble.

S0 docked model 1 of a 20-model ensemble.  Model 1 is an arbitrary choice, so the
floor it imposes may be avoidable.  This computes the Kabsch floor for every
model, on three atom bases, and also for any bound copy of the same toxin taken
from a different interface of the same deposited structure -- that last number is
the experimental reproducibility floor, i.e. how much the same toxin differs
between copies within one entry, which bounds how much of the apo-bound
difference is real conformational change rather than modelling noise.

Usage:
  s0_ensemble_floor.py <native_ligand_pdb> <glob_of_apo_models> [<extra_pdb> ...]
"""
import glob as globmod
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import build_correspondence, parse_pdb_atoms, xyz_array


def kabsch_rmsd(P, Q):
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(V @ Wt))])
    return float(np.sqrt((((Pc @ (V @ D @ Wt)) - Qc) ** 2).sum(1).mean()))


native = sys.argv[1]
pattern = sys.argv[2]
extras = sys.argv[3:]
files = sorted(globmod.glob(pattern)) + list(extras)
files = [f for f in files if "ensemble" not in os.path.basename(f)]

rows = []
for f in files:
    idx, nat_xyz, corr = build_correspondence(f, native)
    atoms = parse_pdb_atoms(f)
    names = [atoms[k]["name"] for k in idx]
    resn = [atoms[k]["resname"] for k in idx]
    # match the docking prep, which stripped the C-terminal NH2 cap
    keep = np.array([r != "NH2" for r in resn])
    P = xyz_array(atoms)[idx][keep]
    Q = nat_xyz[keep]
    nm = [n for n, k in zip(names, keep) if k]
    bb = np.array([n in ("N", "CA", "C", "O") for n in nm])
    ca = np.array([n == "CA" for n in nm])
    rows.append({
        "file": os.path.basename(f),
        "atoms": int(keep.sum()),
        "allatom": round(kabsch_rmsd(P, Q), 3),
        "backbone": round(kabsch_rmsd(P[bb], Q[bb]), 3),
        "ca": round(kabsch_rmsd(P[ca], Q[ca]), 3),
    })

print(f"native reference: {os.path.basename(native)}")
print(f"{'file':<28} {'atoms':>5} {'allatom':>8} {'backbone':>9} {'CA':>7}")
for r in rows:
    print(f"{r['file']:<28} {r['atoms']:>5} {r['allatom']:>8.3f} "
          f"{r['backbone']:>9.3f} {r['ca']:>7.3f}")

model_rows = [r for r in rows if "_m" in r["file"] and "bound" not in r["file"]]
if model_rows:
    print("\nacross ensemble models:")
    for basis in ("allatom", "backbone", "ca"):
        v = np.array([r[basis] for r in model_rows])
        best = model_rows[int(np.argmin(v))]["file"]
        print(f"  {basis:<9} min {v.min():.3f} ({best})  "
              f"median {np.median(v):.3f}  max {v.max():.3f}")
bound_rows = [r for r in rows if "bound" in r["file"]]
for r in bound_rows:
    print(f"\nbound copy from another interface of the same entry: {r['file']}")
    print(f"  allatom {r['allatom']:.3f}  backbone {r['backbone']:.3f}  CA {r['ca']:.3f}"
          "   <- experimental reproducibility floor")

dest = os.environ.get("ENSEMBLE_OUT")
if dest:
    prev = json.load(open(dest)) if os.path.exists(dest) else {}
    prev[os.path.basename(native)] = rows
    json.dump(prev, open(dest, "w"), indent=2)
    print("wrote", dest)
