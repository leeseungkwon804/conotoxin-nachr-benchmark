#!/usr/bin/env python3
"""Compact per-run table, including the score-vs-RMSD columns HQ asked for."""
import glob
import json
import os
import sys

import numpy as np

rmsd_dir = sys.argv[1] if len(sys.argv) > 1 else "rmsd"


def spearman(x, y):
    """Spearman rho with average ranks for ties; no scipy on this pod."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)

    def rank(v):
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

    rx, ry = rank(x), rank(y)
    rx -= rx.mean()
    ry -= ry.mean()
    d = np.sqrt((rx * rx).sum() * (ry * ry).sum())
    return float((rx * ry).sum() / d) if d > 0 else float("nan")


hdr = ["run", "R", "swarms", "best_RMSD", "rank_of_best", "best_score_pctile",
       "top1score_RMSD", "top10", "top100", "rho_pose", "rho_swarm",
       "frac_swarm_rho_neg", "psb_p5", "psb_med", "n<5A", "n<10A",
       "floor", "near_ctr", "score_max", "score_med"]
print("\t".join(hdr))
out = []
for f in sorted(glob.glob(os.path.join(rmsd_dir, "*.json"))):
    s = json.load(open(f))
    if "best_rmsd" not in s:
        continue
    npz = f[:-5] + ".npz"
    d = np.load(npz)
    rms, sco, sw = d["rmsd"].astype(np.float64), d["scoring"].astype(np.float64), d["swarm"]
    psb, psbid = d["per_swarm_best"].astype(np.float64), d["per_swarm_ids"]
    n = len(rms)
    rk = s["best_rmsd_score_rank"]
    pct = (n - rk + 1) / n * 100.0
    # pose-level: negative rho = scoring discriminates correctly (higher score, lower RMSD)
    rho_pose = spearman(sco, rms)
    # swarm-level: each swarm's best score vs that swarm's best RMSD (independent units)
    best_sco = np.array([sco[sw == i].max() for i in psbid])
    rho_swarm = spearman(best_sco, psb)
    # per-swarm internal rho, then how often it points the right way
    rhos = []
    for i in psbid:
        m = sw == i
        if m.sum() > 5:
            rhos.append(spearman(sco[m], rms[m]))
    rhos = np.asarray([r for r in rhos if np.isfinite(r)])
    frac_neg = float((rhos < 0).mean()) if len(rhos) else float("nan")
    sf = s["sampling_floor"]
    row = dict(run=s["run"], R=s["restraints_R"],
               swarms=f"{s['swarms_scored']}/{s['swarms_total']}",
               best_RMSD=round(s["best_rmsd"], 3), rank_of_best=rk,
               best_score_pctile=round(pct, 2),
               top1score_RMSD=round(s["rmsd_of_top1_by_score"], 2),
               top10=round(s["best_rmsd_in_top_by_score"]["10"], 2),
               top100=round(s["best_rmsd_in_top_by_score"]["100"], 2),
               rho_pose=round(rho_pose, 4), rho_swarm=round(rho_swarm, 4),
               frac_swarm_rho_neg=round(frac_neg, 3),
               psb_p5=round(s["per_swarm_best_rmsd"]["p5"], 2),
               psb_med=round(s["per_swarm_best_rmsd"]["median"], 2),
               n_under5=s["per_swarm_best_rmsd"]["n_under_5A"],
               n_under10=s["per_swarm_best_rmsd"]["n_under_10A"],
               floor=round(sf["min_pose_centroid_displacement_A"], 2),
               near_ctr=round(sf["nearest_swarm_centre_to_native_centroid_A"], 2),
               score_max=round(s["scoring"]["max"], 2),
               score_med=round(s["scoring"]["median"], 2))
    out.append(row)
    print("\t".join(str(row[k]) for k in
                    ["run", "R", "swarms", "best_RMSD", "rank_of_best",
                     "best_score_pctile", "top1score_RMSD", "top10", "top100",
                     "rho_pose", "rho_swarm", "frac_swarm_rho_neg", "psb_p5",
                     "psb_med", "n_under5", "n_under10", "floor", "near_ctr",
                     "score_max", "score_med"]))
json.dump(out, open(os.path.join(rmsd_dir, "..", "score_vs_rmsd_table.json"), "w"), indent=2)
