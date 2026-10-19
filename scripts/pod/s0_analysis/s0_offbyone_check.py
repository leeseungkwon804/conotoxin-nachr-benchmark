#!/usr/bin/env python3
"""
Exact test of the gso_<step>.out score/coordinate pairing.

A run of 11 steps saves at step 10 (10 % 10 == 0) and step 11 (11 >= steps), so
two CONSECUTIVE states are on disk.  If the Scoring column belongs to the
position held at the start of the step -- i.e. the pre-move position -- then

    Scoring in gso_11.out  ==  objective evaluated at the coordinates in gso_10.out

and it should NOT equal the objective at gso_11.out's own coordinates.  Both are
checkable directly; no unknown quantities.

Usage: s0_offbyone_check.py <nstep_dir> <swarm_id> <lo_step> <hi_step>
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import QUAT_CONVENTIONS, parse_pdb_atoms, read_gso, rot_matrices, xyz_array

from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.prep.simulation import parse_restraints_file
from lightdock.scoring.vdw.driver import DefinedModelAdapter, DefinedScoringFunction

d = os.path.abspath(sys.argv[1])
sid = int(sys.argv[2])
lo, hi = int(sys.argv[3]), int(sys.argv[4])
conv = [c for c in QUAT_CONVENTIONS if c["name"] == "wxyz"][0]

rp = os.path.join(d, "lightdock_receptor.pdb")
lp = os.path.join(d, "lightdock_ligand.pdb")
a, r, c = parse_complex_from_file(rp)
rec = Complex(c, a, structure_file_name=rp)
a, r, c = parse_complex_from_file(lp)
lig = Complex(c, a, structure_file_name=lp)
rec_keys = None
if os.path.exists(os.path.join(d, "restraints.list")):
    p = parse_restraints_file(os.path.join(d, "restraints.list"))
    rec_keys = p["receptor"]["active"] + p["receptor"]["passive"]
ad = DefinedModelAdapter(rec, lig, rec_keys, None)
sf = DefinedScoringFunction()
print(f"receptor restraints wired into scorer: {len(ad.receptor_model.restraints or {})}")

lig_base = np.asarray(ad.ligand_model.coordinates[0].coordinates, dtype=np.float64).copy()
rc = ad.receptor_model.coordinates[0]
lc = ad.ligand_model.coordinates[0]

t_lo, q_lo, s_lo = read_gso(os.path.join(d, f"swarm_{sid}", f"gso_{lo}.out"))
t_hi, q_hi, s_hi = read_gso(os.path.join(d, f"swarm_{sid}", f"gso_{hi}.out"))
n = min(len(t_lo), len(t_hi))
moved = ~np.all(np.isclose(t_lo[:n], t_hi[:n], atol=1e-9), axis=1)
print(f"glowworms {n}; translation changed between step {lo} and {hi}: {int(moved.sum())}")


def energies(tt, qq):
    R = rot_matrices(qq[:n], conv)
    P = np.einsum("gij,aj->gai", R, lig_base) + tt[:n][:, None, :]
    out = np.empty(n)
    for k in range(n):
        lc.coordinates[:] = P[k]
        out[k] = sf(ad.receptor_model, rc, ad.ligand_model, lc)
    return out


e_lo = energies(t_lo, q_lo)      # objective at P_lo
e_hi = energies(t_hi, q_hi)      # objective at P_hi

res = {"dir": d, "swarm": sid, "lo_step": lo, "hi_step": hi, "n": int(n),
       "n_moved": int(moved.sum()), "tests": {}}
for label, recomputed, reference in (
        (f"E(P_{lo}) vs Scoring in gso_{hi}.out  [pre-move hypothesis]", e_lo, s_hi[:n]),
        (f"E(P_{hi}) vs Scoring in gso_{hi}.out  [post-move hypothesis]", e_hi, s_hi[:n]),
        (f"E(P_{lo}) vs Scoring in gso_{lo}.out  [same-file control]", e_lo, s_lo[:n])):
    diff = np.abs(recomputed - reference)
    res["tests"][label] = {"median_abs_diff": float(np.median(diff)),
                           "max_abs_diff": float(diff.max()),
                           "n_within_0.05": int((diff < 0.05).sum())}
    t = res["tests"][label]
    print(f"  {label}")
    print(f"      median|diff| {t['median_abs_diff']:.6g}   max {t['max_abs_diff']:.6g}   "
          f"within 0.05: {t['n_within_0.05']}/{n}")

pre = res["tests"][f"E(P_{lo}) vs Scoring in gso_{hi}.out  [pre-move hypothesis]"]
post = res["tests"][f"E(P_{hi}) vs Scoring in gso_{hi}.out  [post-move hypothesis]"]
res["verdict"] = ("PRE-MOVE: gso_N.out pairs post-move coordinates with the score of the "
                  "PREVIOUS position" if pre["median_abs_diff"] < post["median_abs_diff"]
                  else "POST-MOVE: the Scoring column matches its own coordinates")
print("VERDICT:", res["verdict"])
dest = os.environ.get("OFFBYONE_OUT")
if dest:
    json.dump(res, open(dest, "w"), indent=2)
    print("wrote", dest)
