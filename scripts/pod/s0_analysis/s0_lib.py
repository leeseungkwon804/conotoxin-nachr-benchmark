#!/usr/bin/env python3
"""
s0_lib -- shared machinery for the S0 3-rung ladder ligand-RMSD analysis.

Geometry, stated explicitly because getting it wrong is the classic trap here:

  lightdock3_setup.py writes lightdock_receptor.pdb / lightdock_ligand.pdb with
  each body translated so its own centroid sits at the origin.  The receptor is
  rigid for the whole simulation (use_anm is false in every S0 run, verified in
  setup.json), so the only transform recorded per glowworm in gso_<step>.out is
  the ligand's (translation t, quaternion q) in the *receptor-centred* frame:

      pose_lightdock = R(q) @ lightdock_ligand_xyz + t

  To compare against the native bound ligand -- which lives in the deposited
  PDB frame -- undo the receptor centring:

      pose_deposited = R(q) @ lightdock_ligand_xyz + t + rec_shift
      rec_shift = receptor.pdb_xyz - lightdock_receptor.pdb_xyz   (a constant)

  No superposition.  The receptor never moves, so an interface-frame RMSD is
  exactly the question the gate asks ("did the peptide land where it belongs"),
  and superposing would hide a translated pose.

The quaternion convention is NOT hardcoded -- calibrate() determines it by
reproducing lgd_generate_conformations.py output, and the resolved convention is
stored in the run's calibration record.
"""

import os
import numpy as np

# ---------------------------------------------------------------- PDB parsing

# (chain, resseq_with_icode, resname, atomname, altloc, xyz)
def parse_pdb_atoms(path, keep_altloc_blank_or_A=True):
    atoms = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            altloc = line[16]
            if keep_altloc_blank_or_A and altloc not in (" ", "A"):
                continue
            atoms.append({
                "chain": line[21],
                "resid": line[22:27].rstrip(),          # resseq + icode
                "resname": line[17:20].strip(),
                "name": line[12:16].strip(),
                "altloc": altloc,
                "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            })
    return atoms


def residue_blocks(atoms):
    """Group atoms into residues, preserving file order."""
    order, blocks = [], {}
    for i, a in enumerate(atoms):
        key = (a["chain"], a["resid"])
        if key not in blocks:
            blocks[key] = {"resname": a["resname"], "idx": []}
            order.append(key)
        blocks[key]["idx"].append(i)
    return order, blocks


def xyz_array(atoms):
    return np.asarray([a["xyz"] for a in atoms], dtype=np.float64)


# ------------------------------------------------------- receptor frame shift

def receptor_shift(run_dir):
    """Constant vector taking lightdock_receptor.pdb coords back to receptor.pdb."""
    raw = parse_pdb_atoms(os.path.join(run_dir, "receptor.pdb"))
    ld = parse_pdb_atoms(os.path.join(run_dir, "lightdock_receptor.pdb"))
    if len(raw) != len(ld):
        # --noh / --noxt dropped atoms: match on identity instead of position.
        key = lambda a: (a["chain"], a["resid"], a["name"])
        raw_map = {key(a): a["xyz"] for a in raw}
        pairs = [(raw_map[key(a)], a["xyz"]) for a in ld if key(a) in raw_map]
        if not pairs:
            raise RuntimeError(f"{run_dir}: no receptor atom correspondence")
        d = np.asarray([p[0] for p in pairs]) - np.asarray([p[1] for p in pairs])
    else:
        d = xyz_array(raw) - xyz_array(ld)
    spread = float(np.max(d.max(0) - d.min(0)))
    if spread > 1e-3:
        raise RuntimeError(
            f"{run_dir}: receptor shift is not a pure translation "
            f"(component spread {spread:.4f} A) -- receptor was moved, abort")
    return d.mean(0), len(d), spread


# -------------------------------------------------------- residue-level align

_GAP = -4.0

def align_resname_sequences(seq_a, seq_b):
    """
    Needleman-Wunsch on 3-letter residue names (identity +2 / mismatch -1).
    Returns list of (i, j) index pairs with identical residue names only.
    Used instead of naive resnum arithmetic so an amidation cap, a numbering
    offset, or a missing residue cannot silently shift the correspondence.
    """
    n, m = len(seq_a), len(seq_b)
    F = np.zeros((n + 1, m + 1))
    F[:, 0] = np.arange(n + 1) * _GAP
    F[0, :] = np.arange(m + 1) * _GAP
    ptr = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left
    ptr[:, 0] = 1
    ptr[0, :] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = 2.0 if seq_a[i - 1] == seq_b[j - 1] else -1.0
            cand = (F[i - 1, j - 1] + s, F[i - 1, j] + _GAP, F[i, j - 1] + _GAP)
            k = int(np.argmax(cand))
            F[i, j] = cand[k]
            ptr[i, j] = k
    pairs, i, j = [], n, m
    while i > 0 or j > 0:
        k = ptr[i, j]
        if k == 0:
            i, j = i - 1, j - 1
            if seq_a[i] == seq_b[j]:
                pairs.append((i, j))
        elif k == 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def build_correspondence(docked_lig_pdb, native_lig_pdb):
    """
    Map atoms of lightdock_ligand.pdb (pose order) onto the native bound ligand.

    Returns (pose_index, native_xyz, report).  pose_index[k] is the row in the
    lightdock_ligand coordinate array matching native_xyz[k].
    """
    dock = parse_pdb_atoms(docked_lig_pdb)
    nat = parse_pdb_atoms(native_lig_pdb)
    d_order, d_blocks = residue_blocks(dock)
    n_order, n_blocks = residue_blocks(nat)
    d_seq = [d_blocks[k]["resname"] for k in d_order]
    n_seq = [n_blocks[k]["resname"] for k in n_order]

    res_pairs = align_resname_sequences(d_seq, n_seq)

    pose_index, native_xyz, matched_atoms = [], [], []
    dock_only, native_only = [], []
    paired_d = {i for i, _ in res_pairs}
    paired_n = {j for _, j in res_pairs}

    for di, ni in res_pairs:
        dk, nk = d_order[di], n_order[ni]
        dmap = {dock[i]["name"]: i for i in d_blocks[dk]["idx"]}
        nmap = {nat[i]["name"]: i for i in n_blocks[nk]["idx"]}
        for name in dmap:
            if name in nmap:
                pose_index.append(dmap[name])
                native_xyz.append(nat[nmap[name]]["xyz"])
                matched_atoms.append((d_blocks[dk]["resname"], dk[1], nk[1], name))
            else:
                dock_only.append((d_blocks[dk]["resname"], dk[1], name))
        for name in nmap:
            if name not in dmap:
                native_only.append((n_blocks[nk]["resname"], nk[1], name))

    for i, k in enumerate(d_order):
        if i not in paired_d:
            for a in d_blocks[k]["idx"]:
                dock_only.append((d_blocks[k]["resname"], k[1], dock[a]["name"]))
    for j, k in enumerate(n_order):
        if j not in paired_n:
            for a in n_blocks[k]["idx"]:
                native_only.append((n_blocks[k]["resname"], k[1], nat[a]["name"]))

    report = {
        "docked_file": docked_lig_pdb,
        "native_file": native_lig_pdb,
        "docked_atoms": len(dock),
        "native_atoms": len(nat),
        "docked_residues": len(d_order),
        "native_residues": len(n_order),
        "residue_pairs": len(res_pairs),
        "matched_atoms": len(pose_index),
        "docked_only_atoms": sorted(dock_only),
        "native_only_atoms": sorted(native_only),
        "docked_seq": d_seq,
        "native_seq": n_seq,
        "residue_pair_map": [(d_order[i][1], n_order[j][1], d_seq[i])
                             for i, j in res_pairs],
        "backbone_matched": sum(1 for m in matched_atoms
                                if m[3] in ("N", "CA", "C", "O")),
        "ca_matched": sum(1 for m in matched_atoms if m[3] == "CA"),
    }
    return (np.asarray(pose_index, dtype=np.int64),
            np.asarray(native_xyz, dtype=np.float64).reshape(-1, 3),
            report)


# ------------------------------------------------------------ gso_<step>.out

def read_gso(path):
    """
    Parse a LightDock gso_<step>.out file.

    Layout (0.9.4):
      #Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring
      (tx, ty, tz, q0, q1, q2, q3)   recid  ligid  luciferin  nn  vision  scoring
    ANM extents would follow the quaternion inside the parentheses; every S0 run
    has use_anm false, so a longer tuple is an error rather than something to
    silently ignore.
    """
    trans, quat, score = [], [], []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("("):
                continue
            close = line.index(")")
            coord = [float(v) for v in line[1:close].split(",")]
            if len(coord) != 7:
                raise RuntimeError(
                    f"{path}: {len(coord)} coordinate fields, expected 7 "
                    "(ANM extents present? use_anm should be false)")
            trans.append(coord[0:3])
            quat.append(coord[3:7])
            score.append(float(line[close + 1:].split()[-1]))
    if not trans:
        raise RuntimeError(f"{path}: no coordinate lines")
    return (np.asarray(trans, dtype=np.float64),
            np.asarray(quat, dtype=np.float64),
            np.asarray(score, dtype=np.float64))


# --------------------------------------------------------------- quaternions

# Candidate conventions.  `order` says how the four numbers in the file map to
# (w, x, y, z); `transpose` flips between v' = q v q* and its inverse.
QUAT_CONVENTIONS = [
    {"name": "wxyz",       "order": (0, 1, 2, 3), "transpose": False},
    {"name": "wxyz_T",     "order": (0, 1, 2, 3), "transpose": True},
    {"name": "xyzw",       "order": (3, 0, 1, 2), "transpose": False},
    {"name": "xyzw_T",     "order": (3, 0, 1, 2), "transpose": True},
]


def rot_matrices(quat, conv):
    """Stack of rotation matrices, shape (n, 3, 3)."""
    o = conv["order"]
    q = quat[:, list(o)]
    q = q / np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    R = np.empty((len(q), 3, 3), dtype=np.float64)
    R[:, 0, 0] = 1 - 2 * (y * y + z * z)
    R[:, 0, 1] = 2 * (x * y - w * z)
    R[:, 0, 2] = 2 * (x * z + w * y)
    R[:, 1, 0] = 2 * (x * y + w * z)
    R[:, 1, 1] = 1 - 2 * (x * x + z * z)
    R[:, 1, 2] = 2 * (y * z - w * x)
    R[:, 2, 0] = 2 * (x * z - w * y)
    R[:, 2, 1] = 2 * (y * z + w * x)
    R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    if conv["transpose"]:
        R = np.transpose(R, (0, 2, 1))
    return R


def poses(lig_xyz, trans, quat, conv, shift=None):
    """
    Apply every (q, t) to lig_xyz.  Returns (n_glowworms, n_atoms, 3).
    lig_xyz must be the lightdock_ligand.pdb coordinates (already origin-centred).
    """
    R = rot_matrices(quat, conv)
    out = np.einsum("gij,aj->gai", R, lig_xyz) + trans[:, None, :]
    if shift is not None:
        out = out + np.asarray(shift, dtype=np.float64)[None, None, :]
    return out


# ------------------------------------------------------------- rank statistics

def rank_avg(v):
    """Average ranks, ties shared. 1-based."""
    v = np.asarray(v, dtype=np.float64)
    o = np.argsort(v, kind="mergesort")
    r = np.empty(len(v), dtype=np.float64)
    sv = v[o]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        r[o[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return r


def spearman(x, y):
    """Spearman rho. scipy is not installed on the pod."""
    rx, ry = rank_avg(x), rank_avg(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


def score_vs_rmsd_row(summary, npz):
    """
    The score-versus-RMSD block: can the scoring function find the pose the
    sampler found?

    Sign convention: LightDock scoring is maximised (write_ranking_to_file sorts
    reverse=True), so a scoring function that discriminates correctly yields
    NEGATIVE rho against RMSD.  Positive rho would mean actively anti-correlated.
    """
    rms = npz["rmsd"].astype(np.float64)
    sco = npz["scoring"].astype(np.float64)
    sw = npz["swarm"]
    psb = npz["per_swarm_best"].astype(np.float64)
    ids = npz["per_swarm_ids"]
    n = len(rms)
    rk = summary["best_rmsd_score_rank"]

    # per-swarm best score, aligned with per_swarm_best
    order = np.argsort(sw, kind="mergesort")
    sw_s, sco_s = sw[order], sco[order]
    bounds = np.searchsorted(sw_s, ids, side="left")
    ends = np.searchsorted(sw_s, ids, side="right")
    best_sco = np.array([sco_s[a:b].max() if b > a else np.nan
                         for a, b in zip(bounds, ends)])

    rhos = []
    for a, b in zip(bounds, ends):
        if b - a > 5:
            idx = order[a:b]
            rhos.append(spearman(sco[idx], rms[idx]))
    rhos = np.asarray([r for r in rhos if np.isfinite(r)])

    return {
        "run": summary["run"],
        "restraints_R": summary.get("restraints_R"),
        "poses": int(n),
        "best_rmsd": summary["best_rmsd"],
        "best_pose_score_rank": rk,
        # 100 = best-scoring pose in the run
        "best_pose_score_percentile": float((n - rk + 1) / n * 100.0),
        "rmsd_of_top1_by_score": summary["rmsd_of_top1_by_score"],
        "best_rmsd_in_top10_by_score": summary["best_rmsd_in_top_by_score"].get("10"),
        "best_rmsd_in_top100_by_score": summary["best_rmsd_in_top_by_score"].get("100"),
        "spearman_pose_level": spearman(sco, rms),
        "spearman_swarm_level": spearman(best_sco, psb),
        "n_swarms": int(len(ids)),
        "median_within_swarm_spearman": float(np.median(rhos)) if len(rhos) else None,
        "frac_swarms_spearman_negative": float((rhos < 0).mean()) if len(rhos) else None,
        "scoring_max": summary["scoring"]["max"],
        "scoring_median": summary["scoring"]["median"],
    }


def rmsd_to_native(lig_xyz_sub, trans, quat, conv, shift, native_xyz):
    """
    Interface-frame ligand RMSD, no superposition, vectorised over glowworms.
    lig_xyz_sub / native_xyz are the matched-atom subsets, same order.
    """
    P = poses(lig_xyz_sub, trans, quat, conv, shift)
    d = P - native_xyz[None, :, :]
    return np.sqrt((d * d).sum(-1).mean(-1))
