#!/usr/bin/env python3
"""
Reads the one-step scratch run and reports which (scoring variant, position)
combination reproduces gso_1.out's Scoring column exactly.

Usage: s0_onestep_check.py <onestep_dir> <swarm_id>
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
sid = int(sys.argv[2]) if len(sys.argv) > 2 else 0
conv = [c for c in QUAT_CONVENTIONS if c["name"] == "wxyz"][0]

# P_0 : initial positions, same 7-column layout, no header lines starting with '#'
P0 = []
for line in open(os.path.join(d, "init", f"initial_positions_{sid}.dat")):
    if line.startswith("#") or not line.strip():
        continue
    v = [float(x) for x in line.split()]
    P0.append(v[:7])
P0 = np.asarray(P0)
t0, q0 = P0[:, 0:3], P0[:, 3:7]

# P_1 and the Scoring column that lightdock wrote next to it
t1, q1, sco = read_gso(os.path.join(d, f"swarm_{sid}", "gso_1.out"))
print(f"glowworms: init {len(t0)}, gso_1 {len(t1)}")
n = min(len(t0), len(t1))
moved = ~np.all(np.isclose(t0[:n], t1[:n], atol=1e-9), axis=1)
print(f"translation changed on step 1 for {int(moved.sum())}/{n} glowworms")


def build(rec_restr, lig_restr):
    rp = os.path.join(d, "lightdock_receptor.pdb")
    lp = os.path.join(d, "lightdock_ligand.pdb")
    a, r, c = parse_complex_from_file(rp)
    rec = Complex(c, a, structure_file_name=rp)
    a, r, c = parse_complex_from_file(lp)
    lig = Complex(c, a, structure_file_name=lp)
    return DefinedModelAdapter(rec, lig, rec_restr, lig_restr), DefinedScoringFunction()


rec_keys = lig_keys = None
rp = os.path.join(d, "restraints.list")
if os.path.exists(rp):
    parsed = parse_restraints_file(rp)
    rec_keys = parsed["receptor"]["active"] + parsed["receptor"]["passive"]
    lig_keys = parsed["ligand"]["active"] + parsed["ligand"]["passive"]

variants = {"A_restraint_free": build(None, None),
            "B_restraints_wired": build(rec_keys, lig_keys)}
for tag, (ad, _) in variants.items():
    print(f"{tag}: parsed receptor restraints = "
          f"{len(ad.receptor_model.restraints or {})}, "
          f"ligand = {len(ad.ligand_model.restraints or {})}")

lig_base = np.asarray(list(variants.values())[0][0]
                      .ligand_model.coordinates[0].coordinates, dtype=np.float64).copy()
res = {"dir": d, "swarm": sid, "n_glowworms": int(n),
       "n_moved_on_step1": int(moved.sum()), "combinations": {}}

for pos_tag, (tt, qq) in (("P0_premove", (t0, q0)), ("P1_postmove", (t1, q1))):
    R = rot_matrices(qq[:n], conv)
    P = np.einsum("gij,aj->gai", R, lig_base) + tt[:n][:, None, :]
    for tag, (ad, sf) in variants.items():
        rc = ad.receptor_model.coordinates[0]
        lc = ad.ligand_model.coordinates[0]
        vals = np.empty(n)
        for k in range(n):
            lc.coordinates[:] = P[k]
            vals[k] = sf(ad.receptor_model, rc, ad.ligand_model, lc)
        diff = np.abs(vals - sco[:n])
        key = f"{tag} @ {pos_tag}"
        res["combinations"][key] = {
            "n_exact_1e-9": int((diff < 1e-9).sum()),
            "n_exact_1e-6": int((diff < 1e-6).sum()),
            "median_abs_diff": float(np.median(diff)),
            "max_abs_diff": float(diff.max()),
        }
        c = res["combinations"][key]
        print(f"  {key:<34} exact(<1e-9)={c['n_exact_1e-9']:>4}/{n}  "
              f"exact(<1e-6)={c['n_exact_1e-6']:>4}/{n}  "
              f"median|diff|={c['median_abs_diff']:.6g}  max={c['max_abs_diff']:.6g}")

win = max(res["combinations"], key=lambda k: res["combinations"][k]["n_exact_1e-6"])
res["winner"] = win
res["winner_exact_frac"] = (res["combinations"][win]["n_exact_1e-6"] / n)
print(f"WINNER: {win}  ({res['combinations'][win]['n_exact_1e-6']}/{n} exact)")
dest = os.environ.get("ONESTEP_OUT")
if dest:
    prev = json.load(open(dest)) if os.path.exists(dest) else []
    json.dump([p for p in prev if p["dir"] != d] + [res], open(dest, "w"), indent=2)
    print("wrote", dest)
