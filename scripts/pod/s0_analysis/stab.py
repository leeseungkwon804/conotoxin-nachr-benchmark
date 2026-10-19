"""Rep-to-rep stability of each candidate metric, to decide which are usable."""
import json, glob
import numpy as np
runs = {}
for f in sorted(glob.glob("rmsd/*.json")):
    s = json.load(open(f))
    if "best_rmsd" not in s: continue
    d = np.load(f[:-5] + ".npz")
    cd = d["centroid_disp"]; sco = d["scoring"]
    top = int(np.argmax(sco))
    runs[s["run"]] = dict(
        best=s["best_rmsd"],
        psb_med=s["per_swarm_best_rmsd"]["median"],
        psb_p5=s["per_swarm_best_rmsd"]["p5"],
        n10=s["per_swarm_best_rmsd"]["n_under_10A"],
        top1=s["rmsd_of_top1_by_score"],
        top1_cdisp=float(cd[top]),
        max_cdisp=float(cd.max()))
groups = {}
for r in runs:
    groups.setdefault(r.rsplit("_rep", 1)[0], []).append(r)
print("%-24s %-16s %8s %8s %8s %8s %8s" % ("rung","metric","rep1","rep2","|diff|","","" ))
for g, rs in sorted(groups.items()):
    if len(rs) < 2: continue
    a, b = runs[sorted(rs)[0]], runs[sorted(rs)[1]]
    for m, lab in (("best","best_rmsd"),("psb_med","per-swarm-best med"),
                   ("psb_p5","per-swarm-best p5"),("n10","n swarms <10A"),
                   ("top1","top1-by-score RMSD")):
        print("%-24s %-16s %8.2f %8.2f %8.2f" % (g, lab, a[m], b[m], abs(a[m]-b[m])))
    print("%-24s %-16s %8.1f %8.1f" % (g, "top1 cdisp", a["top1_cdisp"], b["top1_cdisp"]))
    print("%-24s %-16s %8.1f %8.1f" % (g, "max cdisp (extent)", a["max_cdisp"], b["max_cdisp"]))
    print()
