#!/usr/bin/env python3
"""
Assemble the S0 3-rung ladder and answer the gate questions.

Usage: s0_ladder.py <rmsd_dir> <out_prefix> [--boot 20000] [--seed 20260730]

Reads every <run>.json / <run>.npz written by s0_rmsd_run.py and emits
<out_prefix>.json plus <out_prefix>.md.

Two comparisons, deliberately kept apart:

  WITHIN an anchor  BLIND / XDON / NATIVE were set up at matched swarm counts
                    (A1 390/390/392, A2 134/134/132), so raw best-RMSD is
                    directly comparable.  Reported as-is, plus a rank test on
                    the per-swarm-best distributions.

  BETWEEN anchors   A1 has ~2.9x the swarms of A2, and best-RMSD is a minimum,
                    so raw values are unequal-N order statistics.  Every run is
                    subsampled without replacement to the common minimum swarm
                    count; the min is retaken inside each subsample.  This gives
                    a distribution of "best RMSD you would have seen with N
                    swarms" and is the only across-anchor number quoted.
"""
import glob
import json
import math
import os
import re
import sys

import numpy as np

rmsd_dir = sys.argv[1]
out_prefix = sys.argv[2]
B = int(sys.argv[sys.argv.index("--boot") + 1]) if "--boot" in sys.argv else 20000
SEED = int(sys.argv[sys.argv.index("--seed") + 1]) if "--seed" in sys.argv else 20260730

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import score_vs_rmsd_row, spearman

RUN_RE = re.compile(r"^(A1|A2)_(\w+?)_(BLIND|XDON|NATIVE)_rep(\d+)$")
RUNG_ORDER = {"BLIND": 0, "XDON": 1, "NATIVE": 2}

runs = []
for jf in sorted(glob.glob(os.path.join(rmsd_dir, "*.json"))):
    s = json.load(open(jf))
    m = RUN_RE.match(s["run"])
    if not m:
        continue
    npz = jf[:-5] + ".npz"
    rec = {"summary": s, "anchor": m.group(1), "pdb": m.group(2),
           "rung": m.group(3), "rep": int(m.group(4)),
           "per_swarm_best": None, "score_row": None}
    if os.path.exists(npz) and "best_rmsd" in s:
        d = np.load(npz)
        rec["per_swarm_best"] = d["per_swarm_best"].astype(np.float64)
        try:
            rec["score_row"] = score_vs_rmsd_row(s, d)
        except Exception as exc:                       # never lose the ladder to this
            rec["score_row"] = {"run": s["run"], "error": repr(exc)}
    runs.append(rec)

runs.sort(key=lambda r: (r["anchor"], RUNG_ORDER[r["rung"]], r["rep"]))
if not runs:
    print("no run summaries found in", rmsd_dir)
    sys.exit(1)


def mannwhitney_u(a, b):
    """Two-sided Mann-Whitney U with a normal approximation and tie correction."""
    a, b = np.asarray(a), np.asarray(b)
    n1, n2 = len(a), len(b)
    allv = np.concatenate([a, b])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv))
    sv = allv[order]
    i = 0
    tie_term = 0.0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        t = j - i + 1
        if t > 1:
            tie_term += t ** 3 - t
        i = j + 1
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    mu = n1 * n2 / 2.0
    n = n1 + n2
    sd = np.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1))))
    if sd == 0:
        return float(min(u1, u2)), 1.0
    z = (min(u1, u2) - mu) / sd
    # two-sided normal tail; scipy is not installed on the pod
    p = 1.0 + math.erf(-abs(z) / math.sqrt(2.0))
    return float(min(u1, u2)), float(min(1.0, p))


def boot_min(psb, n_sub, B, rng):
    """Distribution of min(per-swarm-best) over subsamples of size n_sub."""
    N = len(psb)
    if n_sub >= N:
        return np.full(1, float(psb.min())), True   # exact, no subsampling possible
    idx = rng.random((B, N)).argsort(1)[:, :n_sub]
    return psb[idx].min(1), False


complete = [r for r in runs if r["summary"].get("complete")]
n_common = None
psb_runs = [r for r in runs if r["per_swarm_best"] is not None]
if psb_runs:
    n_common = min(len(r["per_swarm_best"]) for r in psb_runs)

rng = np.random.default_rng(SEED)
rows = []
for r in runs:
    s = r["summary"]
    row = {
        "run": s["run"], "anchor": r["anchor"], "rung": r["rung"], "rep": r["rep"],
        "restraints_R": s.get("restraints_R"), "restraints_L": s.get("restraints_L"),
        "swarms_total": s.get("swarms_total"), "swarms_scored": s.get("swarms_scored"),
        "complete": bool(s.get("complete")),
        "matched_atoms": s.get("correspondence", {}).get("matched_atoms"),
        "best_rmsd": s.get("best_rmsd"),
        "best_rmsd_score_rank": s.get("best_rmsd_score_rank"),
        "rmsd_of_top1_by_score": s.get("rmsd_of_top1_by_score"),
        "best_in_top10_by_score": (s.get("best_rmsd_in_top_by_score") or {}).get("10"),
        "best_in_top100_by_score": (s.get("best_rmsd_in_top_by_score") or {}).get("100"),
        "per_swarm_best_median": (s.get("per_swarm_best_rmsd") or {}).get("median"),
        "per_swarm_best_p5": (s.get("per_swarm_best_rmsd") or {}).get("p5"),
        "swarms_under_2A": (s.get("per_swarm_best_rmsd") or {}).get("n_under_2A"),
        "swarms_under_3A": (s.get("per_swarm_best_rmsd") or {}).get("n_under_3A"),
        "swarms_under_5A": (s.get("per_swarm_best_rmsd") or {}).get("n_under_5A"),
        "swarms_under_10A": (s.get("per_swarm_best_rmsd") or {}).get("n_under_10A"),
        "scoring_max": (s.get("scoring") or {}).get("max"),
        "intrinsic_floor_allatom_A": (s.get("intrinsic_floor_A") or {}).get("allatom"),
        "intrinsic_floor_backbone_A": (s.get("intrinsic_floor_A") or {}).get("backbone"),
        "intrinsic_floor_ca_A": (s.get("intrinsic_floor_A") or {}).get("ca"),
        "best_rmsd_excess_over_floor_A": s.get("best_rmsd_excess_over_floor_A"),
        "best_rmsd_backbone": (s.get("by_atom_basis", {}).get("backbone") or {})
            .get("best_rmsd"),
        "best_backbone_excess_A": (s.get("by_atom_basis", {}).get("backbone") or {})
            .get("best_rmsd_excess_over_floor_A"),
        "best_rmsd_ca": (s.get("by_atom_basis", {}).get("ca") or {}).get("best_rmsd"),
        "best_ca_excess_A": (s.get("by_atom_basis", {}).get("ca") or {})
            .get("best_rmsd_excess_over_floor_A"),
        "swarms_under_3A_backbone": (s.get("by_atom_basis", {}).get("backbone") or {})
            .get("n_swarms_under_3A"),
        "sampling_floor_A": (s.get("sampling_floor") or {})
            .get("min_pose_centroid_displacement_A"),
        "escape_fraction": (s.get("sampling_floor") or {}).get("escape_fraction"),
        "nearest_swarm_centre_A": (s.get("sampling_floor") or {})
            .get("nearest_swarm_centre_to_native_centroid_A"),
    }
    if r["per_swarm_best"] is not None and n_common:
        d, exact = boot_min(r["per_swarm_best"], n_common, B, rng)
        row["bootN"] = n_common
        row["boot_best_median"] = float(np.median(d))
        row["boot_best_p5"] = float(np.percentile(d, 5))
        row["boot_best_p95"] = float(np.percentile(d, 95))
        row["boot_exact_no_subsampling"] = bool(exact)
        row["boot_frac_under_3A"] = float((d < 3.0).mean())
    rows.append(row)

# ---- within-anchor ladder
ladder = {}
for anchor in sorted({r["anchor"] for r in runs}):
    ar = [r for r in runs if r["anchor"] == anchor]
    entry = {"rungs": {}, "pdb": ar[0]["pdb"]}
    for rung in ("BLIND", "XDON", "NATIVE"):
        rr = [r for r in ar if r["rung"] == rung]
        if not rr:
            continue
        bests = [r["summary"]["best_rmsd"] for r in rr if "best_rmsd" in r["summary"]]
        pooled = np.concatenate([r["per_swarm_best"] for r in rr
                                 if r["per_swarm_best"] is not None]) \
            if any(r["per_swarm_best"] is not None for r in rr) else None
        entry["rungs"][rung] = {
            "reps": [r["rep"] for r in rr],
            "all_complete": all(r["summary"].get("complete") for r in rr),
            "swarms_scored": [r["summary"].get("swarms_scored") for r in rr],
            "swarms_total": [r["summary"].get("swarms_total") for r in rr],
            "best_rmsd_per_rep": bests,
            "best_rmsd_pooled": float(min(bests)) if bests else None,
            "top1_by_score_rmsd_per_rep": [r["summary"].get("rmsd_of_top1_by_score")
                                           for r in rr],
            "restraints_R": rr[0]["summary"].get("restraints_R"),
            "pooled_per_swarm_best_median": float(np.median(pooled))
            if pooled is not None else None,
            "pooled_per_swarm_best_p5": float(np.percentile(pooled, 5))
            if pooled is not None else None,
            "pooled_swarms_under_3A": int((pooled < 3.0).sum())
            if pooled is not None else None,
            "pooled_swarms_under_5A": int((pooled < 5.0).sum())
            if pooled is not None else None,
            "pooled_swarms_under_10A": int((pooled < 10.0).sum())
            if pooled is not None else None,
            "pooled_n_swarms": int(len(pooled)) if pooled is not None else None,
            "sampling_floor_A": [r["summary"].get("sampling_floor", {})
                                 .get("min_pose_centroid_displacement_A") for r in rr],
        }
    # rank tests on pooled per-swarm-best
    tests = {}
    for a, b in (("BLIND", "XDON"), ("XDON", "NATIVE"), ("BLIND", "NATIVE")):
        pa = [r["per_swarm_best"] for r in ar if r["rung"] == a and r["per_swarm_best"] is not None]
        pb = [r["per_swarm_best"] for r in ar if r["rung"] == b and r["per_swarm_best"] is not None]
        if pa and pb:
            A, Bv = np.concatenate(pa), np.concatenate(pb)
            u, p = mannwhitney_u(A, Bv)
            tests[f"{a}_vs_{b}"] = {
                "n": [int(len(A)), int(len(Bv))],
                "median": [float(np.median(A)), float(np.median(Bv))],
                "U": u, "p_two_sided": p,
                "better": b if np.median(Bv) < np.median(A) else a,
            }
    entry["rank_tests_per_swarm_best"] = tests
    ladder[anchor] = entry

# ---- verdict
GATE = 3.0
verdict = {"gate_threshold_A": GATE, "questions": {}}
any_sub3 = [r["run"] for r in rows if r["best_rmsd"] is not None and r["best_rmsd"] < GATE]
verdict["questions"]["a_sub3A_reached"] = {
    "runs_under_3A": any_sub3,
    "answer": bool(any_sub3),
}
a2 = [r for r in rows if r["anchor"] == "A2" and r["best_rmsd"] is not None]
verdict["questions"]["a_sub3A_full_interface_A2"] = {
    "runs_under_3A": [r["run"] for r in a2 if r["best_rmsd"] < GATE],
    "best_rmsd_by_rung": {rung: min([r["best_rmsd"] for r in a2 if r["rung"] == rung],
                                    default=None)
                          for rung in ("BLIND", "XDON", "NATIVE")},
}
mono = {}
for anchor, e in ladder.items():
    b = {k: v["best_rmsd_pooled"] for k, v in e["rungs"].items()}
    med = {k: v["pooled_per_swarm_best_median"] for k, v in e["rungs"].items()}
    p5 = {k: v["pooled_per_swarm_best_p5"] for k, v in e["rungs"].items()}
    if all(b.get(k) is not None for k in ("BLIND", "XDON", "NATIVE")):
        # the min over ~390 swarms is a single order statistic and therefore
        # noisy; the median and p5 of per-swarm-best carry the same ordering
        # question with far more data behind it, so both are reported.
        mono[anchor] = {
            "per_swarm_best_median": med,
            "per_swarm_best_p5": p5,
            "monotonic_by_median": (med["NATIVE"] < med["XDON"] < med["BLIND"])
            if all(v is not None for v in med.values()) else None,
            "monotonic_by_p5": (p5["NATIVE"] < p5["XDON"] < p5["BLIND"])
            if all(v is not None for v in p5.values()) else None,
            "best_rmsd": b,
            "NATIVE_better_than_XDON": b["NATIVE"] < b["XDON"],
            "XDON_better_than_BLIND": b["XDON"] < b["BLIND"],
            "monotonic_in_disclosure": b["NATIVE"] < b["XDON"] < b["BLIND"],
            "gap_BLIND_minus_NATIVE_A": b["BLIND"] - b["NATIVE"],
            "gap_BLIND_minus_XDON_A": b["BLIND"] - b["XDON"],
            "gap_XDON_minus_NATIVE_A": b["XDON"] - b["NATIVE"],
        }
verdict["questions"]["b_monotonic_in_information"] = mono
verdict["questions"]["c_gap_magnitude"] = {
    anchor: {"gap_BLIND_minus_NATIVE_A": v["gap_BLIND_minus_NATIVE_A"],
             "interpretation_note":
                 "how much of the restrained-docking accuracy is attributable to "
                 "disclosed interface information, at this anchor"}
    for anchor, v in mono.items()}
verdict["all_12_runs_complete"] = len(complete) == 12 and len(runs) == 12
verdict["runs_seen"] = len(runs)
verdict["runs_complete"] = len(complete)
verdict["incomplete_runs"] = [{"run": r["summary"]["run"],
                               "scored": r["summary"].get("swarms_scored"),
                               "total": r["summary"].get("swarms_total")}
                              for r in runs if not r["summary"].get("complete")]

score_rows = [r["score_row"] for r in runs if r["score_row"]]
out = {"rmsd_dir": os.path.abspath(rmsd_dir), "bootstrap_B": B, "bootstrap_seed": SEED,
       "bootstrap_common_N": n_common, "rows": rows, "ladder": ladder,
       "score_vs_rmsd": score_rows, "verdict": verdict}
json.dump(out, open(out_prefix + ".json", "w"), indent=2)
json.dump(score_rows, open(os.path.join(os.path.dirname(out_prefix) or ".",
                                        "score_vs_rmsd_table.json"), "w"), indent=2)

# ---- markdown
def f(x, nd=2):
    return "—" if x is None else f"{x:.{nd}f}"


L = []
L.append("# S0 3-rung ladder — ligand-RMSD vs native bound pose\n")
L.append(f"- runs seen: **{len(runs)}/12**, complete: **{len(complete)}/12**")
L.append(f"- RMSD: interface-frame, no superposition, matched heavy atoms only")
L.append(f"- across-anchor bootstrap: B={B}, common N={n_common} swarms, seed={SEED}\n")
L.append("## Per-run\n")
L.append("| run | R restr | swarms | atoms | **best RMSD (Å)** | rank of best by score | top-1 by score (Å) | best in top-10 (Å) | swarms <3Å | <5Å | <10Å | per-swarm-best p5 / median (Å) | sampling floor (Å) |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    sw = f"{r['swarms_scored']}/{r['swarms_total']}" + ("" if r["complete"] else " ⚠")
    L.append(f"| {r['run']} | {r['restraints_R']} | {sw} | {r['matched_atoms']} | "
             f"**{f(r['best_rmsd'])}** | {r['best_rmsd_score_rank'] or '—'} | "
             f"{f(r['rmsd_of_top1_by_score'])} | {f(r['best_in_top10_by_score'])} | "
             f"{r['swarms_under_3A']} | {r['swarms_under_5A']} | {r['swarms_under_10A']} | "
             f"{f(r['per_swarm_best_p5'])} / {f(r['per_swarm_best_median'])} | "
             f"{f(r['sampling_floor_A'])} |")

L.append("\n## Within-anchor ladder (swarm counts matched by design)\n")
for anchor, e in ladder.items():
    L.append(f"### {anchor} — {e['pdb']}\n")
    L.append("| rung | R restr | swarms | best RMSD (Å) rep1 / rep2 | pooled best (Å) | top-1 by score (Å) | pooled per-swarm-best p5 / median (Å) | pooled <3Å / <5Å / <10Å of n |")
    L.append("|---|---|---|---|---|---|---|---|")
    for rung in ("BLIND", "XDON", "NATIVE"):
        v = e["rungs"].get(rung)
        if not v:
            continue
        per = " / ".join(f(x) for x in v["best_rmsd_per_rep"])
        t1 = " / ".join(f(x) for x in v["top1_by_score_rmsd_per_rep"])
        sw = " / ".join(f"{a}/{b}" for a, b in zip(v["swarms_scored"], v["swarms_total"]))
        L.append(f"| {rung} | {v['restraints_R']} | {sw} | {per} | "
                 f"**{f(v['best_rmsd_pooled'])}** | {t1} | "
                 f"{f(v['pooled_per_swarm_best_p5'])} / {f(v['pooled_per_swarm_best_median'])} | "
                 f"{v['pooled_swarms_under_3A']} / {v['pooled_swarms_under_5A']} / "
                 f"{v['pooled_swarms_under_10A']} of {v['pooled_n_swarms']} |")
    if e["rank_tests_per_swarm_best"]:
        L.append("\nMann–Whitney on per-swarm-best RMSD:\n")
        L.append("| comparison | n | medians (Å) | p (two-sided) | better |")
        L.append("|---|---|---|---|---|")
        for k, t in e["rank_tests_per_swarm_best"].items():
            L.append(f"| {k} | {t['n'][0]} vs {t['n'][1]} | "
                     f"{f(t['median'][0])} vs {f(t['median'][1])} | "
                     f"{t['p_two_sided']:.3g} | {t['better']} |")
    L.append("")

L.append("## Across-anchor, swarm count equalised by subsampling\n")
L.append(f"| run | native N | boot N | boot best RMSD median (Å) | 5–95% | P(best <3Å) | exact? |")
L.append("|---|---|---|---|---|---|---|")
for r in rows:
    if "boot_best_median" not in r:
        continue
    L.append(f"| {r['run']} | {r['swarms_scored']} | {r['bootN']} | "
             f"{f(r['boot_best_median'])} | {f(r['boot_best_p5'])}–{f(r['boot_best_p95'])} | "
             f"{r['boot_frac_under_3A']:.3f} | "
             f"{'yes (N==bootN)' if r['boot_exact_no_subsampling'] else 'no'} |")

L.append("\n## Is the 3 Å gate reachable at all? (intrinsic rigid-docking floor)\n")
L.append("The docked ligand is an **apo NMR conformer treated as rigid**. It cannot be "
         "superimposed on the bound conformer, so the optimally-superposed (Kabsch) "
         "deviation between the two is a hard lower bound on any pose's RMSD. A gate "
         "threshold at or below that floor is unreachable no matter how good the docking "
         "is.\n")
L.append("| run | floor all-atom (Å) | best all-atom (Å) | excess | floor backbone (Å) | best backbone (Å) | excess | swarms <3Å backbone | floor CA (Å) | best CA (Å) |")
L.append("|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    if r["intrinsic_floor_allatom_A"] is None:
        continue
    L.append(f"| {r['run']} | {f(r['intrinsic_floor_allatom_A'])} | {f(r['best_rmsd'])} | "
             f"{f(r['best_rmsd_excess_over_floor_A'])} | "
             f"{f(r['intrinsic_floor_backbone_A'])} | {f(r['best_rmsd_backbone'])} | "
             f"{f(r['best_backbone_excess_A'])} | {r['swarms_under_3A_backbone']} | "
             f"{f(r['intrinsic_floor_ca_A'])} | {f(r['best_rmsd_ca'])} |")
L.append("\n`excess` = best achieved − floor, i.e. how much of the error the docking is "
         "responsible for. A run whose floor already exceeds the gate threshold cannot "
         "pass it, and its failure says nothing about the docking.\n")

L.append("\n## Gate\n")
q = verdict["questions"]
L.append(f"**(a) Any run reaching <3 Å best ligand-RMSD:** "
         f"{', '.join(q['a_sub3A_reached']['runs_under_3A']) or '**NONE**'}")
L.append(f"\n**(a') Full (+)/(−) interface anchor A2 (7EGR), best RMSD by rung:** "
         + ", ".join(f"{k} {f(v)} Å" for k, v in
                     q["a_sub3A_full_interface_A2"]["best_rmsd_by_rung"].items()))
L.append("\n**(b) Monotonic in disclosed information?**\n")
if mono:
    L.append("By best RMSD (a single order statistic — noisy):\n")
    L.append("| anchor | BLIND (Å) | XDON (Å) | NATIVE (Å) | NATIVE<XDON | XDON<BLIND | monotonic |")
    L.append("|---|---|---|---|---|---|---|")
    for anchor, v in mono.items():
        b = v["best_rmsd"]
        L.append(f"| {anchor} | {f(b['BLIND'])} | {f(b['XDON'])} | {f(b['NATIVE'])} | "
                 f"{v['NATIVE_better_than_XDON']} | {v['XDON_better_than_BLIND']} | "
                 f"**{v['monotonic_in_disclosure']}** |")
    L.append("\nBy per-swarm-best distribution (hundreds of values per rung — the reliable read):\n")
    L.append("| anchor | statistic | BLIND (Å) | XDON (Å) | NATIVE (Å) | monotonic |")
    L.append("|---|---|---|---|---|---|")
    for anchor, v in mono.items():
        for lab, key, flag in (("median", "per_swarm_best_median", "monotonic_by_median"),
                               ("p5", "per_swarm_best_p5", "monotonic_by_p5")):
            d = v[key]
            L.append(f"| {anchor} | {lab} | {f(d.get('BLIND'))} | {f(d.get('XDON'))} | "
                     f"{f(d.get('NATIVE'))} | **{v[flag]}** |")
    L.append("\n**(c) Gap attributable to disclosure (BLIND − NATIVE):**\n")
    for anchor, v in mono.items():
        L.append(f"- {anchor}: **{v['gap_BLIND_minus_NATIVE_A']:+.2f} Å** "
                 f"(BLIND−XDON {v['gap_BLIND_minus_XDON_A']:+.2f}, "
                 f"XDON−NATIVE {v['gap_XDON_minus_NATIVE_A']:+.2f})")
else:
    L.append("_not all three rungs available yet_")
L.append("\n## Can the scoring function find the pose the sampler found?\n")
L.append("Sign convention: LightDock scoring is **maximised** "
         "(`write_ranking_to_file` sorts `reverse=True`), so a scoring function "
         "that discriminates correctly gives a **negative** Spearman ρ against "
         "RMSD. Percentile 100 = best-scoring pose in the run.\n")
L.append("| run | R | poses | best RMSD (Å) | rank of best pose by score | **score percentile of best pose** | top-1 by score (Å) | best in top-10 | best in top-100 | ρ pose-level | ρ swarm-level | median ρ within swarm | frac. swarms ρ<0 |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for sr in score_rows:
    if "error" in sr:
        L.append(f"| {sr['run']} | — | — | — | — | — | — | — | — | — | — | — | — |")
        continue
    L.append(f"| {sr['run']} | {sr['restraints_R']} | {sr['poses']} | "
             f"{f(sr['best_rmsd'])} | {sr['best_pose_score_rank']} | "
             f"**{f(sr['best_pose_score_percentile'])}** | "
             f"{f(sr['rmsd_of_top1_by_score'])} | {f(sr['best_rmsd_in_top10_by_score'])} | "
             f"{f(sr['best_rmsd_in_top100_by_score'])} | "
             f"{f(sr['spearman_pose_level'], 4)} | {f(sr['spearman_swarm_level'], 4)} | "
             f"{f(sr['median_within_swarm_spearman'], 4)} | "
             f"{f(sr['frac_swarms_spearman_negative'], 3)} |")

# per anchor x rung aggregation of the same block
agg = {}
for r in runs:
    if not r["score_row"] or "error" in r["score_row"]:
        continue
    agg.setdefault((r["anchor"], r["rung"]), []).append(r["score_row"])
if agg:
    L.append("\nAggregated per anchor × rung:\n")
    L.append("| anchor | rung | R | reps | best RMSD (Å) | score percentile of best pose | top-1 by score (Å) | ρ pose-level | ρ swarm-level |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for (anchor, rung) in sorted(agg, key=lambda k: (k[0], RUNG_ORDER[k[1]])):
        v = agg[(anchor, rung)]
        j = lambda key, nd=2: " / ".join(f(x[key], nd) for x in v)
        L.append(f"| {anchor} | {rung} | {v[0]['restraints_R']} | {len(v)} | "
                 f"{j('best_rmsd')} | {j('best_pose_score_percentile')} | "
                 f"{j('rmsd_of_top1_by_score')} | {j('spearman_pose_level', 4)} | "
                 f"{j('spearman_swarm_level', 4)} |")

if verdict["incomplete_runs"]:
    L.append("\n⚠ **Incomplete runs — numbers above are lower bounds on coverage:**\n")
    for r in verdict["incomplete_runs"]:
        L.append(f"- {r['run']}: {r['scored']}/{r['total']} swarms")
open(out_prefix + ".md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
