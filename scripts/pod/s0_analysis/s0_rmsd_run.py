#!/usr/bin/env python3
"""
Ligand-RMSD of every LightDock pose in one S0 run against the native bound pose.

Usage:
  s0_rmsd_run.py <run_dir> <native_ligand_pdb> <step> <out_dir> [--conv wxyz]

Writes into <out_dir>/<run_name>.{json,npz}:
  json  -- correspondence report, completion, and every summary statistic
  npz   -- swarm_id / glowworm_id / rmsd / scoring for all poses read

Metrics, and why each one is here:
  best_rmsd                    min over all poses.  Sampling ceiling: "could the
                               search reach the native pose at all".  This is the
                               kill-gate number the brief asks for.
  rmsd_of_top1_by_score        the pose you would actually hand to a biologist.
  best_rmsd_in_topN_by_score   how deep into the score-ranked list you must go.
  per_swarm_best_rmsd          one value per swarm; the unit for the unequal-N
                               bootstrap between anchors (390 vs 134 swarms).
"""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import (QUAT_CONVENTIONS, build_correspondence, parse_pdb_atoms,
                    poses, read_gso, receptor_shift, rmsd_to_native, xyz_array)

run_dir = os.path.abspath(sys.argv[1])
native_pdb = sys.argv[2]
step = int(sys.argv[3])
out_dir = sys.argv[4]
conv_name = "wxyz"
if "--conv" in sys.argv:
    conv_name = sys.argv[sys.argv.index("--conv") + 1]
conv = [c for c in QUAT_CONVENTIONS if c["name"] == conv_name][0]

run_name = os.path.basename(run_dir)
os.makedirs(out_dir, exist_ok=True)

lig_all = parse_pdb_atoms(os.path.join(run_dir, "lightdock_ligand.pdb"))
lig_xyz_all = xyz_array(lig_all)
pose_idx, native_xyz, corr = build_correspondence(
    os.path.join(run_dir, "lightdock_ligand.pdb"), native_pdb)
lig_sub = lig_xyz_all[pose_idx]
shift, n_shift_atoms, spread = receptor_shift(run_dir)

# Atom subsets.  All-atom RMSD is dominated by long flexible side chains (Arg,
# Trp, Tyr) whose apo NMR rotamers bear no relation to the bound ones, so the
# same run is also scored on backbone and CA only.  Which basis is meaningful is
# decided by the intrinsic floor below, not by preference.
sub_names = [lig_all[k]["name"] for k in pose_idx]
MASKS = {
    "allatom": np.ones(len(sub_names), dtype=bool),
    "backbone": np.array([n in ("N", "CA", "C", "O") for n in sub_names]),
    "ca": np.array([n == "CA" for n in sub_names]),
}


def kabsch_rmsd(P, Q):
    """RMSD after optimal rigid superposition -- the floor for rigid docking."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(V @ Wt))])
    return float(np.sqrt((((Pc @ (V @ D @ Wt)) - Qc) ** 2).sum(1).mean()))


# A rigid body cannot beat its own optimally-superposed deviation from the target
# conformer.  Recorded per basis so the gate threshold can be read against it.
intrinsic_floor = {k: round(kabsch_rmsd(lig_sub[m], native_xyz[m]), 3)
                   for k, m in MASKS.items() if m.sum() >= 3}

setup = json.load(open(os.path.join(run_dir, "setup.json")))
if setup.get("use_anm"):
    raise RuntimeError(f"{run_name}: use_anm is true -- pose reconstruction "
                       "would need normal-mode extents, refusing")

restr = os.path.join(run_dir, "restraints.list")
n_R = n_L = 0
if os.path.exists(restr):
    for line in open(restr):
        if line.startswith("R "):
            n_R += 1
        elif line.startswith("L "):
            n_L += 1

swarm_dirs = sorted(glob.glob(os.path.join(run_dir, "swarm_*")),
                    key=lambda p: int(p.rsplit("_", 1)[1]))
n_swarms = len(swarm_dirs)

# --- sampling-floor bookkeeping -------------------------------------------
# The translation in gso_*.out is the ligand centroid (lightdock_ligand.pdb is
# origin-centred), so |t + shift - native_centroid| is the pose's centroid
# displacement.  RMSD can never be smaller than that displacement, which makes
# it a hard floor: if the floor is already above 3 A the gate result is a
# swarm-placement artefact, not a statement about scoring or physics.
native_centroid = native_xyz.mean(0)
centres_pdb = os.path.join(run_dir, "init", "swarm_centers.pdb")
swarm_centres = None
if os.path.exists(centres_pdb):
    swarm_centres = xyz_array(parse_pdb_atoms(centres_pdb)) + shift
swarm_radius = float(setup.get("swarm_radius", 10.0))

sw_id, gw_id, rms, sco = [], [], [], []
per_swarm_best, per_swarm_ids, per_swarm_top1_rmsd = [], [], []
cdisp_all, escaped, n_trans = [], 0, 0
basis_rms, basis_swarm_best = {}, {}
missing = []
for d in swarm_dirs:
    sid = int(d.rsplit("_", 1)[1])
    f = os.path.join(d, f"gso_{step}.out")
    if not os.path.exists(f):
        missing.append(sid)
        continue
    try:
        trans, quat, score = read_gso(f)
    except RuntimeError as exc:            # truncated file from a live run
        missing.append(sid)
        print(f"  skip swarm_{sid}: {exc}", file=sys.stderr)
        continue
    P = poses(lig_sub, trans, quat, conv, shift)
    dsq = ((P - native_xyz[None, :, :]) ** 2).sum(-1)
    r = np.sqrt(dsq.mean(-1))
    for key, m in MASKS.items():
        if key == "allatom":
            continue
        rr = np.sqrt(dsq[:, m].mean(-1))
        basis_rms.setdefault(key, []).append(rr.astype(np.float32))
        basis_swarm_best.setdefault(key, []).append(float(rr.min()))
    cent = trans + shift
    cdisp_all.append(np.linalg.norm(cent - native_centroid, axis=1).astype(np.float32))
    if swarm_centres is not None and sid < len(swarm_centres):
        escaped += int((np.linalg.norm(cent - swarm_centres[sid], axis=1)
                        > swarm_radius).sum())
        n_trans += len(cent)
    sw_id.append(np.full(len(r), sid, dtype=np.int32))
    gw_id.append(np.arange(len(r), dtype=np.int32))
    rms.append(r.astype(np.float32))
    sco.append(score.astype(np.float32))
    per_swarm_best.append(float(r.min()))
    per_swarm_ids.append(sid)
    per_swarm_top1_rmsd.append(float(r[int(np.argmax(score))]))

n_done = len(per_swarm_ids)
summary = {
    "run": run_name,
    "run_dir": run_dir,
    "step": step,
    "quat_convention": conv_name,
    "receptor_shift": shift.tolist(),
    "receptor_shift_component_spread_A": spread,
    "receptor_atoms": n_shift_atoms,
    "restraints_R": n_R,
    "restraints_L": n_L,
    "swarms_total": n_swarms,
    "swarms_scored": n_done,
    "swarms_missing_gso": missing if len(missing) <= 40 else
                          missing[:40] + [f"...+{len(missing)-40} more"],
    "complete": n_done == n_swarms and n_swarms > 0,
    "correspondence": corr,
}

if n_done == 0:
    summary["note"] = f"no gso_{step}.out yet"
    json.dump(summary, open(os.path.join(out_dir, run_name + ".json"), "w"), indent=2)
    print(json.dumps({k: summary[k] for k in
                      ("run", "swarms_total", "swarms_scored", "note")}, indent=2))
    sys.exit(0)

sw_id = np.concatenate(sw_id)
gw_id = np.concatenate(gw_id)
rms = np.concatenate(rms)
sco = np.concatenate(sco)
per_swarm_best = np.asarray(per_swarm_best)
per_swarm_top1_rmsd = np.asarray(per_swarm_top1_rmsd)

order_score = np.argsort(-sco)          # LightDock scoring: higher is better
best_i = int(np.argmin(rms))
top1_i = int(order_score[0])

summary.update({
    "poses_scored": int(len(rms)),
    "glowworms_per_swarm": int(setup["glowworms"]),
    "best_rmsd": float(rms[best_i]),
    "best_rmsd_swarm": int(sw_id[best_i]),
    "best_rmsd_glowworm": int(gw_id[best_i]),
    "best_rmsd_scoring": float(sco[best_i]),
    "best_rmsd_score_rank": int(np.where(order_score == best_i)[0][0]) + 1,
    "rmsd_of_top1_by_score": float(rms[top1_i]),
    "top1_by_score_scoring": float(sco[top1_i]),
    "top1_by_score_swarm": int(sw_id[top1_i]),
    "best_rmsd_in_top_by_score": {
        str(n): float(rms[order_score[:n]].min())
        for n in (1, 10, 100, 200, 400, 1000) if n <= len(rms)
    },
    "rmsd_percentiles_all_poses": {
        str(p): float(np.percentile(rms, p)) for p in (0.1, 1, 5, 25, 50, 75, 100)
    },
    "per_swarm_best_rmsd": {
        "n": int(len(per_swarm_best)),
        "min": float(per_swarm_best.min()),
        "p5": float(np.percentile(per_swarm_best, 5)),
        "p25": float(np.percentile(per_swarm_best, 25)),
        "median": float(np.median(per_swarm_best)),
        "p75": float(np.percentile(per_swarm_best, 75)),
        "max": float(per_swarm_best.max()),
        "n_under_2A": int((per_swarm_best < 2.0).sum()),
        "n_under_3A": int((per_swarm_best < 3.0).sum()),
        "n_under_5A": int((per_swarm_best < 5.0).sum()),
        "n_under_10A": int((per_swarm_best < 10.0).sum()),
    },
    "per_swarm_top1_by_score_rmsd_median": float(np.median(per_swarm_top1_rmsd)),
    "scoring": {
        "max": float(sco.max()), "median": float(np.median(sco)),
        "min": float(sco.min()),
    },
    "n_poses_under_2A": int((rms < 2.0).sum()),
    "n_poses_under_3A": int((rms < 3.0).sum()),
    "n_poses_under_5A": int((rms < 5.0).sum()),
})

summary["atom_basis"] = {k: int(m.sum()) for k, m in MASKS.items()}
summary["intrinsic_floor_A"] = intrinsic_floor
summary["floor_note"] = (
    "intrinsic_floor_A is the Kabsch (optimally superposed) RMSD between the apo "
    "docked conformer and the native bound conformer on that atom basis. No rigid-body "
    "pose can score below it, so a gate threshold at or under the floor is unreachable "
    "by construction.")
basis_summary = {}
for key in ("backbone", "ca"):
    if key not in basis_rms:
        continue
    rr = np.concatenate(basis_rms[key])
    psb_k = np.asarray(basis_swarm_best[key])
    bi = int(np.argmin(rr))
    basis_summary[key] = {
        "atoms": int(MASKS[key].sum()),
        "intrinsic_floor_A": intrinsic_floor.get(key),
        "best_rmsd": float(rr[bi]),
        "best_rmsd_excess_over_floor_A": (float(rr[bi]) - intrinsic_floor[key])
        if key in intrinsic_floor else None,
        "n_poses_under_2A": int((rr < 2.0).sum()),
        "n_poses_under_3A": int((rr < 3.0).sum()),
        "per_swarm_best_median": float(np.median(psb_k)),
        "per_swarm_best_min": float(psb_k.min()),
        "n_swarms_under_2A": int((psb_k < 2.0).sum()),
        "n_swarms_under_3A": int((psb_k < 3.0).sum()),
    }
summary["by_atom_basis"] = basis_summary
summary["best_rmsd_excess_over_floor_A"] = (
    summary["best_rmsd"] - intrinsic_floor["allatom"]
    if "allatom" in intrinsic_floor else None)

cdisp = np.concatenate(cdisp_all)
sampling = {
    "swarm_radius_A": swarm_radius,
    "min_pose_centroid_displacement_A": float(cdisp.min()),
    "median_pose_centroid_displacement_A": float(np.median(cdisp)),
    "n_poses_centroid_within_3A_of_native": int((cdisp < 3.0).sum()),
    "glowworms_escaping_own_swarm_sphere": escaped,
    "glowworms_with_centre_known": n_trans,
    "escape_fraction": (escaped / n_trans) if n_trans else None,
}
if swarm_centres is not None:
    dc = np.linalg.norm(swarm_centres - native_centroid, axis=1)
    da = np.sqrt(((swarm_centres[:, None, :] - native_xyz[None, :, :]) ** 2)
                 .sum(-1)).min(1)
    sampling.update({
        "nearest_swarm_centre_to_native_centroid_A": float(dc.min()),
        "nearest_swarm_centre_to_native_atom_A": float(da.min()),
        "n_swarm_centres_within_radius_of_native_centroid": int((dc <= swarm_radius).sum()),
        "native_centroid_inside_some_swarm_sphere": bool(dc.min() <= swarm_radius),
    })
summary["sampling_floor"] = sampling

extra = {}
for key in ("backbone", "ca"):
    if key in basis_rms:
        extra[f"rmsd_{key}"] = np.concatenate(basis_rms[key])
        extra[f"per_swarm_best_{key}"] = np.asarray(basis_swarm_best[key], dtype=np.float32)
np.savez_compressed(os.path.join(out_dir, run_name + ".npz"),
                    swarm=sw_id, glowworm=gw_id, rmsd=rms, scoring=sco,
                    centroid_disp=cdisp,
                    per_swarm_best=per_swarm_best.astype(np.float32),
                    per_swarm_ids=np.asarray(per_swarm_ids, dtype=np.int32),
                    per_swarm_top1_rmsd=per_swarm_top1_rmsd.astype(np.float32),
                    **extra)
json.dump(summary, open(os.path.join(out_dir, run_name + ".json"), "w"), indent=2)

print(f"{run_name}: swarms {n_done}/{n_swarms}  poses {len(rms)}  "
      f"matched_atoms {corr['matched_atoms']}  "
      f"BEST {summary['best_rmsd']:.3f} A  "
      f"top1-by-score {summary['rmsd_of_top1_by_score']:.3f} A  "
      f"swarms<3A {summary['per_swarm_best_rmsd']['n_under_3A']}  "
      f"floor(min centroid disp) {sampling['min_pose_centroid_displacement_A']:.3f} A  "
      f"escape {sampling['escape_fraction']}")
print(f"  intrinsic floor (Kabsch apo-vs-bound): " +
      ", ".join(f"{k} {v:.3f} A" for k, v in intrinsic_floor.items()) +
      f"  |  best-minus-floor allatom "
      f"{summary['best_rmsd_excess_over_floor_A']:.3f} A")
for k, v in basis_summary.items():
    print(f"  [{k}] best {v['best_rmsd']:.3f} A (floor {v['intrinsic_floor_A']:.3f}, "
          f"excess {v['best_rmsd_excess_over_floor_A']:.3f})  "
          f"swarms<3A {v['n_swarms_under_3A']}  poses<3A {v['n_poses_under_3A']}")
