#!/usr/bin/env python3
"""
Two questions about LightDock 0.9.4's gso_<step>.out, settled by measurement.

Q1  Do receptor restraints reach the SCORING function during docking, or only
    swarm placement?  Source reading says only placement: the docking-phase
    CommandLineParser defines no `receptor_restraints`, and
    docking_multiprocessing.set_scoring_function() wraps
    `parser.args.receptor_restraints["active"]` in a bare `except: pass`.
    If that is right, restraints_satisfied() sees an empty dict and returns 0.0,
    so score = energy * (1 + 0 + 0) for every run and rung-to-rung score
    comparison is legitimate.

Q2  Is the `Scoring` column paired with the coordinates written beside it?
    The GSO step is  update_luciferin() -> movement_phase() -> save(), and
    Glowworm.compute_luciferin() re-evaluates only `if self.moved or step == 0`.
    So `scoring` should be the objective at the position held at the START of the
    step, while the saved coordinates are the position AFTER that step's move.

Test.  Recompute the vdW score at the WRITTEN coordinates two ways --
restraint-free, and with the run's own restraints wired in the way
lightdock3_setup.py parses them -- and look for an exact-match spike.
Glowworms that did not move on the final step have P_N == P_{N-1}, so for them
the correctly-specified scoring function must reproduce the file's value to
floating-point.  Whichever variant produces that spike is the one the simulation
actually used; the size of the spike is the non-mover fraction.

Usage: s0_probe_scoring.py <run_dir> <step> [n_swarms]
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

run_dir = os.path.abspath(sys.argv[1])
step = int(sys.argv[2])
n_swarms = int(sys.argv[3]) if len(sys.argv) > 3 else 6
conv = [c for c in QUAT_CONVENTIONS if c["name"] == "wxyz"][0]
run_name = os.path.basename(run_dir)

rec_pdb = os.path.join(run_dir, "lightdock_receptor.pdb")
lig_pdb = os.path.join(run_dir, "lightdock_ligand.pdb")


def build(rec_restr, lig_restr):
    a, r, c = parse_complex_from_file(rec_pdb)
    rec = Complex(c, a, structure_file_name=rec_pdb)
    a, r, c = parse_complex_from_file(lig_pdb)
    lig = Complex(c, a, structure_file_name=lig_pdb)
    ad = DefinedModelAdapter(rec, lig, rec_restr, lig_restr)
    return ad, DefinedScoringFunction()


# variant A: exactly what the docking phase ends up doing (restraints = None)
adA, sfA = build(None, None)
# variant B: restraints wired in, as strings in the form the vdW adapter matches
restr_path = os.path.join(run_dir, "restraints.list")
rec_keys = lig_keys = None
if os.path.exists(restr_path):
    parsed = parse_restraints_file(restr_path)
    rec_keys = parsed["receptor"]["active"] + parsed["receptor"]["passive"]
    lig_keys = parsed["ligand"]["active"] + parsed["ligand"]["passive"]
adB, sfB = build(rec_keys, lig_keys)

print(f"{run_name}")
print(f"  restraints.list R-active+passive = {len(rec_keys) if rec_keys else 0}, "
      f"L = {len(lig_keys) if lig_keys else 0}")
print(f"  variant A parsed_restraints: receptor {len(adA.receptor_model.restraints or {})}, "
      f"ligand {len(adA.ligand_model.restraints or {})}")
print(f"  variant B parsed_restraints: receptor {len(adB.receptor_model.restraints or {})}, "
      f"ligand {len(adB.ligand_model.restraints or {})}")
if rec_keys:
    print(f"  variant B sample keys matched: "
          f"{sorted(list((adB.receptor_model.restraints or {}).keys()))[:4]}")
    print(f"  restraints.list sample: {rec_keys[:4]}")

lig_base = np.asarray(adA.ligand_model.coordinates[0].coordinates, dtype=np.float64).copy()
file_xyz = xyz_array(parse_pdb_atoms(lig_pdb))
assert np.allclose(file_xyz, lig_base, atol=1e-3), "adapter ligand coords != file"

res = {"run": run_name, "step": step, "variants": {}}
for tag, ad, sf in (("A_restraint_free", adA, sfA), ("B_restraints_wired", adB, sfB)):
    rc = ad.receptor_model.coordinates[0]
    lc = ad.ligand_model.coordinates[0]
    diffs, rels, vals, refs = [], [], [], []
    for sid in range(n_swarms):
        f = os.path.join(run_dir, f"swarm_{sid}", f"gso_{step}.out")
        if not os.path.exists(f):
            continue
        trans, quat, score = read_gso(f)
        R = rot_matrices(quat, conv)
        P = np.einsum("gij,aj->gai", R, lig_base) + trans[:, None, :]
        for k in range(len(trans)):
            lc.coordinates[:] = P[k]
            v = sf(ad.receptor_model, rc, ad.ligand_model, lc)
            vals.append(v)
            refs.append(score[k])
            diffs.append(abs(v - score[k]))
            rels.append(abs(v - score[k]) / max(1e-9, abs(score[k])))
    diffs = np.asarray(diffs)
    vals = np.asarray(vals)
    refs = np.asarray(refs)
    exact = int((diffs < 1e-6).sum())
    res["variants"][tag] = {
        "n_poses": int(len(diffs)),
        "n_exact_match": exact,
        "frac_exact_match": float(exact / len(diffs)) if len(diffs) else None,
        "median_abs_diff": float(np.median(diffs)),
        "max_abs_diff": float(diffs.max()),
        "ratio_recomputed_over_file_median": float(np.median(vals / np.where(
            np.abs(refs) < 1e-9, np.nan, refs))),
    }
    v = res["variants"][tag]
    print(f"  [{tag}] n={v['n_poses']}  exact(<1e-6)={v['n_exact_match']} "
          f"({v['frac_exact_match']:.3f})  median|diff|={v['median_abs_diff']:.6g}  "
          f"max|diff|={v['max_abs_diff']:.6g}  median(recomp/file)="
          f"{v['ratio_recomputed_over_file_median']:.4f}")

a = res["variants"]["A_restraint_free"]
b = res["variants"]["B_restraints_wired"]
if a["n_exact_match"] > 0 and a["n_exact_match"] >= b["n_exact_match"]:
    verdict = ("Q1: restraints do NOT reach the scoring function -- the "
               "restraint-free variant reproduces the file exactly for non-movers. "
               "score = energy x 1 in every rung.")
elif b["n_exact_match"] > a["n_exact_match"]:
    verdict = ("Q1: restraints DO reach the scoring function -- only the wired "
               "variant reproduces the file. score = energy x (1 + f_R).")
else:
    verdict = "Q1: inconclusive -- neither variant reproduces the file exactly."
res["verdict_Q1"] = verdict
res["verdict_Q2"] = (
    f"Q2: {a['n_poses'] - a['n_exact_match']} of {a['n_poses']} poses differ, "
    f"consistent with the Scoring column belonging to the pre-move position; "
    f"the exact-match fraction {a['frac_exact_match']:.3f} is the non-mover fraction.")
print(" ", verdict)
print(" ", res["verdict_Q2"])
dest = os.environ.get("PROBE_OUT")
if dest:
    prev = json.load(open(dest)) if os.path.exists(dest) else []
    json.dump([p for p in prev if p["run"] != run_name] + [res], open(dest, "w"), indent=2)
    print("  wrote", dest)
