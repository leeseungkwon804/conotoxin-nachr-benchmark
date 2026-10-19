#!/usr/bin/env python3
"""
Re-score poses with a RESTRAINT-FREE vdW function, so scores are comparable
across rungs.

Why this is necessary.  lightdock/scoring/vdw/driver.py returns

    score = energy * (1 + f_R + f_L),   energy = -vdw_energy

and ScoringFunction.restraints_satisfied() returns 0.0 for an *empty* restraint
set.  Every S0 run has L = 0, so f_L = 0 throughout -- but f_R is identically 0
in the BLIND rung (no receptor restraints declared) while it ranges over [0, 1]
in XDON and NATIVE.  The achievable multiplier is therefore x1 for BLIND and up
to x2 for the restrained rungs.  This is the same class of defect as the x3-vs-x2
ligand-restraint asymmetry recorded in AGENT.md section 2 item 1, reappearing
between BLIND and restrained rungs.  Raw `Scoring` columns from gso_*.out must
not be compared across rungs.

Instantiating the adapter with no restraints (exactly what lgd_calculate_scoring.py
does) forces f_R = f_L = 0 for every run, so what comes out is the bare
truncated vdW energy on a single common scale.

Usage:
  s0_rescore.py <run_dir> <native_ligand_pdb> <rmsd_npz> <step> <out_npz>
                [--mode all|subset] [--seed 20260730]

subset mode (default) rescores: the top 500 poses by native score, the 500
lowest-RMSD poses, and a seeded random 2000 -- enough for cross-rung comparison
of the poses that matter plus an unbiased correlation sample, without a
million-pose sweep.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import (QUAT_CONVENTIONS, parse_pdb_atoms, read_gso,
                    receptor_shift, rot_matrices, xyz_array)

from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from lightdock.scoring.vdw.driver import DefinedModelAdapter, DefinedScoringFunction

run_dir = os.path.abspath(sys.argv[1])
native_pdb = sys.argv[2]
rmsd_npz = sys.argv[3]
step = int(sys.argv[4])
out_npz = sys.argv[5]
mode = sys.argv[sys.argv.index("--mode") + 1] if "--mode" in sys.argv else "subset"
seed = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20260730
conv = [c for c in QUAT_CONVENTIONS if c["name"] == "wxyz"][0]
run_name = os.path.basename(run_dir)

# --- restraint-free scoring function on the lightdock-frame structures
rec_pdb = os.path.join(run_dir, "lightdock_receptor.pdb")
lig_pdb = os.path.join(run_dir, "lightdock_ligand.pdb")
a, r, c = parse_complex_from_file(rec_pdb)
receptor = Complex(c, a, structure_file_name=rec_pdb)
a, r, c = parse_complex_from_file(lig_pdb)
ligand = Complex(c, a, structure_file_name=lig_pdb)
adapter = DefinedModelAdapter(receptor, ligand)          # no restraints -> f_R = f_L = 0
sf = DefinedScoringFunction()
assert not adapter.receptor_model.restraints, "receptor restraints leaked into rescoring"
assert not adapter.ligand_model.restraints, "ligand restraints leaked into rescoring"

rec_coords = adapter.receptor_model.coordinates[0]
lig_coords0 = adapter.ligand_model.coordinates[0]
lig_base = np.asarray(lig_coords0.coordinates, dtype=np.float64).copy()
n_lig = len(lig_base)

# The adapter reorders/filters nothing for vdw, but verify the atom count matches
# the file we build poses from, otherwise the transform would be applied to the
# wrong rows.
lig_file_xyz = xyz_array(parse_pdb_atoms(lig_pdb))
if len(lig_file_xyz) != n_lig or not np.allclose(lig_file_xyz, lig_base, atol=1e-3):
    raise RuntimeError(f"{run_name}: adapter ligand coordinates do not match "
                       f"{lig_pdb} ({n_lig} vs {len(lig_file_xyz)})")

# --- which poses
d = np.load(rmsd_npz)
sw, gw, rms, sco = d["swarm"], d["glowworm"], d["rmsd"], d["scoring"]
n_all = len(rms)
rand_sel = None
if mode == "all":
    sel = np.arange(n_all)
    rand_sel = sel
else:
    rng = np.random.default_rng(seed)
    top_score = np.argsort(-sco)[:500]
    top_rmsd = np.argsort(rms)[:500]
    # the random block is the only unbiased part; percentiles and rho use it alone
    rand_sel = rng.choice(n_all, size=min(2000, n_all), replace=False)
    sel = np.unique(np.concatenate([top_score, top_rmsd, rand_sel]))
print(f"{run_name}: rescoring {len(sel)} of {n_all} poses (mode={mode})")

# --- group by swarm so each gso file is read once
order = np.argsort(sw[sel], kind="mergesort")
sel = sel[order]
new_score = np.full(len(sel), np.nan)
t0 = time.time()
i = 0
while i < len(sel):
    sid = int(sw[sel[i]])
    j = i
    while j + 1 < len(sel) and int(sw[sel[j + 1]]) == sid:
        j += 1
    trans, quat, _ = read_gso(os.path.join(run_dir, f"swarm_{sid}", f"gso_{step}.out"))
    idx = gw[sel[i:j + 1]].astype(int)
    R = rot_matrices(quat[idx], conv)
    P = np.einsum("gij,aj->gai", R, lig_base) + trans[idx][:, None, :]
    for k in range(len(idx)):
        lig_coords0.coordinates[:] = P[k]
        new_score[i + k] = sf(adapter.receptor_model, rec_coords,
                              adapter.ligand_model, lig_coords0)
    i = j + 1
dt = time.time() - t0
print(f"  {len(sel)} rescored in {dt:.1f}s ({dt / max(1, len(sel)) * 1000:.2f} ms/pose)")

restr = os.path.join(run_dir, "restraints.list")
n_R = 0
if os.path.exists(restr):
    n_R = sum(1 for line in open(restr) if line.startswith("R "))

# Expected disagreement with LightDock's own Scoring column, for reference: the
# column is the objective at the PREVIOUS position (FINDINGS.md F3), and for a
# restrained run it also carries the (1 + f_R) factor this rescoring drops (F4).
# So a nonzero difference here is the documented behaviour, not an error.
ref = sco[sel].astype(np.float64)
diff = np.abs(new_score - ref)
ratio = new_score / np.where(np.abs(ref) < 1e-9, np.nan, ref)
check = {"median_abs_diff_vs_file": float(np.nanmedian(diff)),
         "max_abs_diff_vs_file": float(np.nanmax(diff)),
         "median_ratio_recomputed_over_file": float(np.nanmedian(ratio)),
         "note": ("differs from the file by design: file column is pre-move (F3) and "
                  "carries (1 + f_R) when restraints are present (F4)")}
print(f"  vs LightDock's Scoring column: median|diff| {check['median_abs_diff_vs_file']:.4g}, "
      f"median ratio {check['median_ratio_recomputed_over_file']:.4f} "
      f"(expected <1 for R>0, ~1 for R=0)")

# --- rung-comparable metrics from the restraint-free score
r_sub = rms[sel].astype(np.float64)
rand_mask = np.isin(sel, rand_sel) if mode != "all" else np.ones(len(sel), bool)
best_i = int(np.argmin(r_sub))
top_i = int(np.argmax(new_score))
# percentile of the best-RMSD pose, estimated against the unbiased random sample
rnd_scores = new_score[rand_mask]
pctile = float((rnd_scores <= new_score[best_i]).mean() * 100.0)
stale_rank_of_corrected_top = int((sco[sel] > sco[sel][top_i]).sum()) + 1
comparable = {
    "best_rmsd_in_subset": float(r_sub[best_i]),
    "restraintfree_score_of_best_rmsd_pose": float(new_score[best_i]),
    "best_rmsd_pose_score_percentile_est": pctile,
    "percentile_estimated_from_n_random": int(rand_mask.sum()),
    "rmsd_of_top1_by_restraintfree_score": float(r_sub[top_i]),
    "restraintfree_score_max_in_subset": float(new_score.max()),
    "restraintfree_score_median_random": float(np.median(rnd_scores)),
    # if the corrected argmax sits deep inside the stale-ranked set, the subset
    # almost certainly contains the true global argmax
    "stale_rank_of_corrected_top1": stale_rank_of_corrected_top,
    "spearman_restraintfree_score_vs_rmsd_random": None,
}
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from s0_lib import spearman
    comparable["spearman_restraintfree_score_vs_rmsd_random"] = spearman(
        new_score[rand_mask], r_sub[rand_mask])
except Exception as exc:
    comparable["spearman_error"] = repr(exc)

np.savez_compressed(out_npz, sel=sel, swarm=sw[sel], glowworm=gw[sel],
                    rmsd=rms[sel], score_native=sco[sel],
                    score_restraintfree=new_score.astype(np.float32),
                    is_random_sample=rand_mask)
meta = {"run": run_name, "mode": mode, "seed": seed, "n_rescored": int(len(sel)),
        "n_all": int(n_all), "restraints_R": n_R,
        "ms_per_pose": dt / max(1, len(sel)) * 1000,
        "vs_file_column": check, "rung_comparable": comparable}
json.dump(meta, open(out_npz.replace(".npz", ".meta.json"), "w"), indent=2)
print(f"  [rung-comparable] top-1 by restraint-free score -> RMSD "
      f"{comparable['rmsd_of_top1_by_restraintfree_score']:.3f} A  |  "
      f"best-RMSD pose percentile {pctile:.2f}  |  "
      f"rho {comparable['spearman_restraintfree_score_vs_rmsd_random']}  |  "
      f"stale rank of corrected top1 = {stale_rank_of_corrected_top}")
print("  wrote", out_npz)
