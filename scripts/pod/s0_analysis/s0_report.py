#!/usr/bin/env python3
"""
Final S0 summary: the tables that go to HQ.

Reads everything s0_ladder.py and the diagnostic scripts wrote and emits one
markdown file. Run from the s0_analysis directory.

Usage: s0_report.py <s0_analysis_dir> <out_md>
"""
import json
import os
import sys

AN = os.path.abspath(sys.argv[1])
out_md = sys.argv[2]

RUNG_ORDER = {"BLIND": 0, "XDON": 1, "NATIVE": 2}


def load(name, default=None):
    p = os.path.join(AN, name)
    if not os.path.exists(p):
        return default
    try:
        return json.load(open(p))
    except Exception:
        return default


lad = load("ladder.json")
if lad is None:
    print("ladder.json not found — nothing to report")
    sys.exit(1)
floors = load("ensemble_floor.json", {})
cov = {r["run"]: r for r in (load("swarm_coverage.json", []) or [])}
probe = {r["run"]: r for r in (load("scoring_probe.json", []) or [])}
obo = load("offbyone_check.json", {})

rows = lad["rows"]
by_run = {r["run"]: r for r in rows}
score = {s["run"]: s for s in lad.get("score_vs_rmsd", []) if "error" not in s}
GATE = 3.0


def f(x, nd=2):
    return "—" if x is None else f"{x:.{nd}f}"


L = []
L.append("# S0 3-rung ladder — results\n")
n_done = sum(1 for r in rows if r["complete"])
L.append(f"Runs scored: **{n_done}/12 complete**"
         + ("" if n_done == 12 else "  ⚠ *partial — incomplete runs flagged below*") + "\n")
L.append("Anchors, from live RCSB queries on 2026-07-30:\n")
L.append("| | A1 | A2 |")
L.append("|---|---|---|")
L.append("| PDB | **9NX0** | **7EGR** |")
L.append("| complex | human α7 nAChR (CHRNA7) + α-conotoxin ImI | *Aplysia californica* AChBP + RgIA |")
L.append("| method / resolution | cryo-EM 3.06 Å | X-ray 2.5 Å |")
L.append("| released | 2026-02-04 | 2022-04-06 |")
L.append("| docked receptor | chains D,E — 392 res each (1–478, **ECD + TMD**), 6326 heavy atoms | chains H,G — 206 res each (19–224, **ECD only**), 3276 heavy atoms |")
L.append("| docked ligand (apo, rigid) | ImI, 1IMI model 1, 90 atoms | RgIA, 2JUT model 1, 105 atoms |")
L.append("| native reference | 9NX0 chain I | 7EGR chain P |")
L.append("| matched atoms | 90 / 90 | 105 / 105 |")
L.append("\nBoth anchors are subunit **pairs cut from a homopentamer**, so both present a "
         "complete (+)/(−) interface. A2 is *not* a human receptor.\n")

# ---------------- gate
L.append("\n## (a) Kill gate — does any run reach < 3 Å best ligand-RMSD?\n")
sub3 = [r["run"] for r in rows if r["best_rmsd"] is not None and r["best_rmsd"] < GATE]
L.append(f"**All-atom: {'YES — ' + ', '.join(sub3) if sub3 else 'NO. Nothing under 3 Å.'}**\n")
bb3 = [r["run"] for r in rows if r.get("best_rmsd_backbone") is not None
       and r["best_rmsd_backbone"] < GATE]
L.append(f"**Backbone: {'YES — ' + ', '.join(bb3) if bb3 else 'NO. Nothing under 3 Å.'}**\n")
L.append("But the threshold has to be read against the **intrinsic floor** — the docked "
         "ligand is an apo conformer treated as rigid, and the optimally superposed "
         "(Kabsch) deviation from the bound conformer is a hard lower bound on any pose:\n")
L.append("| run | floor all-atom | best all-atom | **excess** | floor backbone | best backbone | **excess** | floor CA | best CA |")
L.append("|---|---|---|---|---|---|---|---|---|")
for r in rows:
    if r.get("intrinsic_floor_allatom_A") is None:
        continue
    L.append(f"| {r['run']}{'' if r['complete'] else ' ⚠'} | "
             f"{f(r['intrinsic_floor_allatom_A'])} | {f(r['best_rmsd'])} | "
             f"**{f(r['best_rmsd_excess_over_floor_A'])}** | "
             f"{f(r['intrinsic_floor_backbone_A'])} | {f(r['best_rmsd_backbone'])} | "
             f"**{f(r['best_backbone_excess_A'])}** | "
             f"{f(r['intrinsic_floor_ca_A'])} | {f(r['best_rmsd_ca'])} |")
L.append("\n`excess = best − floor` is the part the docking is responsible for.\n")
L.append("**A2's all-atom floor is 4.453 Å — above the gate. The all-atom gate could not "
         "be passed at A2 no matter how good the docking was**; the minimum over all 20 "
         "deposited 2JUT models is still 3.248 Å. The informative number at A2 is the "
         "excess, and on backbone (floor 2.397 Å) the gate *was* reachable and still failed.\n")

# ---------------- ladder
L.append("\n## (b) Is accuracy monotonic in disclosed information?\n")
for anchor in sorted(lad["ladder"]):
    e = lad["ladder"][anchor]
    L.append(f"### {anchor} — {e['pdb']}\n")
    L.append("| rung | R restraints | swarms | best all-atom (Å) rep1 / rep2 | pooled best | best backbone rep1 / rep2 | per-swarm-best median (Å) | top-1 by score (Å) |")
    L.append("|---|---|---|---|---|---|---|---|")
    for rung in ("BLIND", "XDON", "NATIVE"):
        v = e["rungs"].get(rung)
        if not v:
            continue
        reps = [by_run[f"{anchor}_{e['pdb']}_{rung}_rep{i}"]
                for i in v["reps"] if f"{anchor}_{e['pdb']}_{rung}_rep{i}" in by_run]
        per = " / ".join(f(x) for x in v["best_rmsd_per_rep"])
        bbs = " / ".join(f(r.get("best_rmsd_backbone")) for r in reps)
        sw = " / ".join(f"{a}/{b}" for a, b in zip(v["swarms_scored"], v["swarms_total"]))
        t1 = " / ".join(f(x) for x in v["top1_by_score_rmsd_per_rep"])
        L.append(f"| {rung} | {v['restraints_R']} | {sw} | {per} | "
                 f"**{f(v['best_rmsd_pooled'])}** | {bbs} | "
                 f"{f(v['pooled_per_swarm_best_median'])} | {t1} |")
    mono = lad["verdict"]["questions"]["b_monotonic_in_information"].get(anchor)
    if mono:
        b = mono["best_rmsd"]
        inverted = b["BLIND"] < b["XDON"] < b["NATIVE"]
        if inverted:
            verdict_txt = ("**strictly monotonic in the WRONG direction** — more disclosed "
                           "information, worse best pose")
        elif b["BLIND"] < b["NATIVE"]:
            verdict_txt = "**inverted but not strictly ordered** — BLIND best, rungs not clean"
        else:
            verdict_txt = "**improves with disclosure**"
        L.append(f"\nBest-RMSD ordering: BLIND {f(b['BLIND'])} / XDON {f(b['XDON'])} / "
                 f"NATIVE {f(b['NATIVE'])} Å → {verdict_txt} "
                 f"(gap BLIND−NATIVE = {mono['gap_BLIND_minus_NATIVE_A']:+.3f} Å; the "
                 f"hypothesised direction would make this positive)")
    L.append("")

L.append("Two things move in **opposite** directions, which is the substance of the answer:\n")
L.append("- **Best achievable pose gets worse with more information.** Restrained swarm "
         "placement concentrates sampling on the restrained patch; where that patch only "
         "partly covers the true footprint it pulls the search to a partly-wrong location. "
         "Measured: nearest swarm centre to the native ligand centroid is 6.08 Å for A2 "
         "BLIND but 10.67 Å for A2 XDON.")
L.append("- **The *typical* pose gets markedly better with more information**, consistently "
         "at both anchors: per-swarm-best median A1 43.77 → 20.55 Å (BLIND → XDON), "
         "A2 27.64 → 14.57 / 17.86 Å. Restraints do concentrate sampling near the site; "
         "they just do not improve its best outcome.")
L.append("- **The score-selected pose shows no consistent pattern.** A2 improves then flattens "
         "(top-1 by score 25.25 → 12.12 / 11.82 → 12.76 / 12.96 Å) while A1 is non-monotone "
         "(11.47 → 15.84 → 11.10 Å), with XDON the worst rung at A1 and the best at A2. "
         "⚠ Independently of that, the column is not comparable across rungs at all: "
         "receptor restraints multiply the score by (1 + f_R), giving BLIND a ×1 ceiling and "
         "restrained rungs up to ×2 (`FINDINGS.md` F4). Re-scored on one common scale the A2 "
         "ordering reverses outright — see the rung-comparable table below. **This row must "
         "not be quoted as evidence about pose selection.**\n")

# ---------------- (c)
L.append("\n## (c) How much of the 2.14 Å figure is attributable to disclosed information?\n")
qc = lad["verdict"]["questions"].get("c_gap_magnitude", {})
for anchor, v in qc.items():
    L.append(f"- **{anchor}: {v['gap_BLIND_minus_NATIVE_A']:+.3f} Å.** "
             "A negative number means full disclosure made the best pose *worse*.")
if len(qc) == 2:
    _bm = lad["verdict"]["questions"]["b_monotonic_in_information"]
    (_na, _va), (_nb, _vb) = list(qc.items())
    _ga, _gb = _va["gap_BLIND_minus_NATIVE_A"], _vb["gap_BLIND_minus_NATIVE_A"]
    _ra = round(_bm[_na]["best_rmsd"]["BLIND"], 2) - round(_bm[_na]["best_rmsd"]["NATIVE"], 2)
    _rb = round(_bm[_nb]["best_rmsd"]["BLIND"], 2) - round(_bm[_nb]["best_rmsd"]["NATIVE"], 2)
    L.append(
        f"\n⚠ **The two anchors landing on the same 2-decimal value ({_ga:+.2f} / {_gb:+.2f} Å) "
        f"is a coincidence, not a copied cell.** Each gap is computed independently from its own "
        f"anchor's pooled best-RMSD: {_na} = {_bm[_na]['best_rmsd']['BLIND']:.4f} − "
        f"{_bm[_na]['best_rmsd']['NATIVE']:.4f} = {_ga:+.4f} Å; {_nb} = "
        f"{_bm[_nb]['best_rmsd']['BLIND']:.4f} − {_bm[_nb]['best_rmsd']['NATIVE']:.4f} = "
        f"{_gb:+.4f} Å. They differ by {abs(_ga - _gb):.4f} Å, which is why they are quoted to "
        f"three decimals here. Differencing the 2-decimal values shown in the tables above "
        f"instead gives {_na} {_ra:+.2f} and {_nb} {_rb:+.2f} Å — the difference of rounded "
        f"numbers, not the rounded difference; the values quoted here are the latter.")
L.append("\nSo on this ladder, **none of it.** Disclosure does not buy best-pose accuracy; "
         "it costs it. The 2.14 Å therefore cannot be attributed to the restraints.\n")
L.append("What the ladder does **not** isolate, stated plainly: the WP1 2.14 Å run differs "
         "from S0's NATIVE rung in three ways at once, and S0 changed all three. Per "
         "AGENT.md/STATUS.md (the original run directory was lost with the previous pod and "
         "could not be re-inspected):\n")
L.append("| | WP1 2.14 Å run | S0 A2 NATIVE |")
L.append("|---|---|---|")
L.append("| receptor | 6HY7 chain A — **one α9 ECD subunit**, 1701 atoms, 210 res (2–213): half an interface | 7EGR chains H,G — complete (+)/(−) interface |")
L.append("| ligand | bound coordinates (bound–bound self-redock) | **apo NMR**, rigid |")
L.append("| disclosure | 12 restraints, 12/12 identical to 6HY7 native contacts | 33 restraints = 7EGR's own native contacts |")
L.append("\nDisclosure is maximal in both, so it is the one factor the comparison rules "
         "*out*. The remaining candidates are the apo-vs-bound ligand and the half-vs-full "
         "interface. Isolating them needs one run per factor, not more ladder rungs.\n")

# ---------------- score vs rmsd
L.append("\n## Can the scoring function find the pose the sampler found?\n")
L.append("LightDock scoring is **maximised**, so a correctly discriminating function gives "
         "a **negative** Spearman ρ against RMSD. Percentile 100 = best-scoring pose.\n")
L.append("| run | R | poses | best RMSD (Å) | rank of best pose | **score percentile of best pose** | top-1 by score (Å) | best in top-10 | best in top-100 | ρ pose | ρ swarm | frac. swarms ρ<0 |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    s = score.get(r["run"])
    if not s:
        continue
    L.append(f"| {s['run']} | {s['restraints_R']} | {s['poses']} | {f(s['best_rmsd'])} | "
             f"{s['best_pose_score_rank']} | **{f(s['best_pose_score_percentile'])}** | "
             f"{f(s['rmsd_of_top1_by_score'])} | {f(s['best_rmsd_in_top10_by_score'])} | "
             f"{f(s['best_rmsd_in_top100_by_score'])} | {f(s['spearman_pose_level'], 3)} | "
             f"{f(s['spearman_swarm_level'], 3)} | "
             f"{f(s['frac_swarms_spearman_negative'], 2)} |")
agg = {}
for r in rows:
    s = score.get(r["run"])
    if s:
        agg.setdefault((r["anchor"], r["rung"]), []).append(s)
if agg:
    L.append("\nPer anchor × rung:\n")
    L.append("| anchor | rung | R | best RMSD (Å) | score percentile of best pose | top-1 by score (Å) | ρ pose | ρ swarm |")
    L.append("|---|---|---|---|---|---|---|---|")
    for k in sorted(agg, key=lambda k: (k[0], RUNG_ORDER[k[1]])):
        v = agg[k]
        j = lambda key, nd=2: " / ".join(f(x[key], nd) for x in v)
        L.append(f"| {k[0]} | {k[1]} | {v[0]['restraints_R']} | {j('best_rmsd')} | "
                 f"{j('best_pose_score_percentile')} | {j('rmsd_of_top1_by_score')} | "
                 f"{j('spearman_pose_level', 3)} | {j('spearman_swarm_level', 3)} |")
rs = {}
rdir = os.path.join(AN, "rescore")
if os.path.isdir(rdir):
    for fn in sorted(os.listdir(rdir)):
        if fn.endswith(".meta.json"):
            m = json.load(open(os.path.join(rdir, fn)))
            rs[m["run"]] = m
if rs:
    L.append("\n### Same question with a rung-comparable score\n")
    L.append("The published `Scoring` column carries the `(1 + f_R)` restraint factor, which "
             "is identically 1 for BLIND and up to 2 for restrained rungs, so it cannot be "
             "compared across rungs (`FINDINGS.md` F4). Re-scoring the same poses with the "
             "restraints removed puts every rung on one scale. Percentile and ρ are estimated "
             "from an unbiased random sub-sample of poses; the estimator was validated on "
             "A2 BLIND rep1, where it gives 7.05 / −0.073 against 7.64 / −0.080 from the full "
             "ranking.\n")
    L.append("| run | R | median(restraint-free / published) | top-1 by **published** score → RMSD | top-1 by **restraint-free** score → RMSD | shift | published rank of restraint-free top-1 | best-pose percentile (restraint-free) |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in rows:
        m = rs.get(r["run"])
        if not m:
            continue
        c, v = m["rung_comparable"], m["vs_file_column"]
        pub = score.get(r["run"], {}).get("rmsd_of_top1_by_score")
        free = c["rmsd_of_top1_by_restraintfree_score"]
        shift = (free - pub) if pub is not None else None
        L.append(f"| {r['run']} | {m['restraints_R']} | "
                 f"{f(v['median_ratio_recomputed_over_file'], 4)} | {f(pub)} | {f(free)} | "
                 f"{('+' if shift and shift > 0 else '')}{f(shift)} | "
                 f"{c['stale_rank_of_corrected_top1']} | "
                 f"{f(c['best_rmsd_pose_score_percentile_est'])} |")
    L.append("\nThe `median(restraint-free / published)` column is the measurement of F4: "
             "≈1.0 wherever R = 0 and <1 wherever restraints are declared.\n")

L.append("\n## ⚠ Which of these metrics actually reproduce\n")
L.append("Each rung was run twice, identical except for the GSO seed, so rep1-vs-rep2 measures "
         "how much of each metric is noise (`FINDINGS.md` F7).\n")

# computed, not hardcoded, so it stays correct as runs land
pairs = {}
for r in rows:
    if r["complete"]:
        pairs.setdefault(r["run"].rsplit("_rep", 1)[0], []).append(r)
pairs = {k: sorted(v, key=lambda x: x["rep"]) for k, v in pairs.items() if len(v) >= 2}
cd_extent, cd_top1 = {}, {}
for k, v in pairs.items():
    for r in v:
        p = os.path.join(AN, "rmsd", r["run"] + ".npz")
        if os.path.exists(p):
            import numpy as _np
            d = _np.load(p)
            cd = d["centroid_disp"]
            cd_extent[r["run"]] = float(cd.max())
            cd_top1[r["run"]] = float(cd[int(_np.argmax(d["scoring"]))])

METRICS = [("per_swarm_best_median", "per-swarm-best **median**", 1),
           ("per_swarm_best_p5", "per-swarm-best **p5**", 1),
           ("best_rmsd", "`best_rmsd`", 1),
           ("rmsd_of_top1_by_score", "**top-1-by-score → RMSD**", 1),
           ("best_pose_score_pct", "percentile of best pose", 1)]
if pairs:
    order = sorted(pairs)
    L.append("| metric | " + " | ".join(k.replace("_9NX0", "").replace("_7EGR", "")
                                        for k in order) + " | max Δ | verdict |")
    L.append("|---" * (len(order) + 3) + "|")
    for key, label, _ in METRICS:
        cells, deltas = [], []
        for k in order:
            a, b = pairs[k]
            if key == "best_pose_score_pct":
                va = score.get(a["run"], {}).get("best_pose_score_percentile")
                vb = score.get(b["run"], {}).get("best_pose_score_percentile")
            else:
                va, vb = a.get(key), b.get(key)
            if va is None or vb is None:
                cells.append("—")
                continue
            cells.append(f"{va:.2f} / {vb:.2f}")
            deltas.append(abs(va - vb))
        mx = max(deltas) if deltas else None
        if mx is None:
            verdict = "—"
        elif mx <= 1.2:
            verdict = "**reproduces**"
        elif mx <= 2.0:
            verdict = "usable, state a band"
        else:
            verdict = "**unusable**"
        L.append(f"| {label} | " + " | ".join(cells) + f" | {f(mx)} | {verdict} |")
    L.append("\nThe two failing metrics are each defined by a **single pose** out of "
             "26 400–78 000, so nothing averages the noise away.")
    a1b = [r["run"] for k, v in pairs.items() if k.startswith("A1") and "BLIND" in k for r in v]
    if len(a1b) == 2 and all(x in cd_top1 for x in a1b):
        L.append(f"\nThe A1 BLIND blow-up is mechanical: the top-scoring pose's centroid sits "
                 f"**{cd_top1[a1b[0]]:.1f} Å** from native in rep1 and **{cd_top1[a1b[1]]:.1f} Å** "
                 f"in rep2. That is possible because A1's receptor is the full-length pair "
                 f"including the TMD — poses reach "
                 f"≈{max(cd_extent[x] for x in a1b):.0f} Å of centroid displacement at A1 "
                 + (f"versus ≈{max(v for k, v in cd_extent.items() if k.startswith('A2')):.0f} Å "
                    f"across A2" if any(k.startswith("A2") for k in cd_extent) else "")
                 + ". On a receptor that large the best-vdW pose can sit almost anywhere on it, "
                   "and the seed decides where.\n")
L.append("**So: quote the per-swarm-best distribution and `best_rmsd` (with a ±1.5 Å band); do "
         "not interpret any score-selected single-pose number.** AGENT.md §2 already prescribed "
         "the per-swarm-best distribution with rank tests over raw maxima — this is an "
         "independent confirmation of it.\n")
L.append("Spearman ρ has thousands of pairs behind it and is stable within a run, but its "
         "**sign is not uniform across runs**: negative (correct direction, |ρ| 0.06–0.46) in "
         "most, while A1 XDON rep1 is −0.003 pose-level and **+0.115 swarm-level**, and both A1 "
         "restrained runs go positive on the restraint-free score (+0.074, +0.041). "
         "Discrimination collapses specifically in A1's restrained rungs, not at A2. Quote ρ "
         "per run, never as \"always negative\".\n")
L.append("**None of this touches the (a)–(c) verdicts above**, which rest on `best_rmsd`, the "
         "intrinsic floors and the per-swarm-best distribution.\n")

# ---------------- bootstrap
L.append("\n## Across-anchor comparison, swarm count equalised\n")
L.append(f"`best_rmsd` is a minimum, so it grows with swarm count and A1 (390–392) cannot "
         f"be compared raw to A2 (132–134). Every run's per-swarm-best array is subsampled "
         f"without replacement to the common minimum "
         f"(**N = {lad['bootstrap_common_N']}**), B = {lad['bootstrap_B']}, "
         f"seed {lad['bootstrap_seed']}, and the minimum retaken inside each subsample.\n")
L.append("| run | native N | boot best RMSD median (Å) | 5–95 % | P(best < 3 Å) | subsampled? |")
L.append("|---|---|---|---|---|---|")
for r in rows:
    if "boot_best_median" not in r:
        continue
    L.append(f"| {r['run']} | {r['swarms_scored']} | {f(r['boot_best_median'])} | "
             f"{f(r['boot_best_p5'])}–{f(r['boot_best_p95'])} | "
             f"{r['boot_frac_under_3A']:.3f} | "
             f"{'no (already at N)' if r['boot_exact_no_subsampling'] else 'yes'} |")

L.append("\n## Instrument findings\n")
L.append("Separate from the docking result, and the part that survives it — full detail with "
         "how each was established in `FINDINGS.md`:\n")
L.append("1. **The 3 Å all-atom gate was unreachable at A2 by construction** (floor 4.453 Å; "
         "3.248 Å for the best of 20 apo models). Two RgIA copies inside 7EGR itself differ "
         "by 2.526 Å all-atom, so the threshold sits barely above the ground truth's own noise.")
L.append("2. **Apo-source choice alone can decide a 3 Å benchmark.** For the same toxin ImI, "
         "the achievable floor spans **1.02–4.57 Å** across three deposited apo entries "
         "(1G2G 1.02–2.13, 1IMI 1.76–2.74, 1CNL 4.02–4.57) — a spread wider than the "
         "threshold.")
L.append("3. **`gso_<step>.out` pairs post-move coordinates with pre-move scores.** Proven by "
         "a consecutive-step experiment: E(P₁₀) matches gso_11's `Scoring` at median "
         "|Δ| 0.0081 (186/200 within 0.05) while the same-file pairing fails at median 5.68. "
         "`lgd_rank_swarm.py` sorts on that column, so LightDock's own ranking tool inherits it.")
L.append("4. **Receptor restraints enter the score**, giving BLIND a ×1 ceiling and restrained "
         "rungs up to ×2. Verified from each run's `lightdock.info` and by recomputation "
         "(median recomputed/file = 1.0000 with restraints wired, 0.9257 / 0.9694 without). "
         "Raw scores are not comparable across rungs.")
_esc = [r["escape_fraction"] for r in rows if r.get("escape_fraction") is not None]
L.append(f"5. **Swarm placement biases but does not bound the search** — "
         f"{100 * min(_esc):.0f}–{100 * max(_esc):.0f} % of glowworms end outside their "
         f"initial swarm sphere, lowest in the restrained rungs.")

inc = lad["verdict"].get("incomplete_runs") or []
if inc:
    L.append("\n⚠ **Incomplete runs — their numbers are lower bounds on coverage:**\n")
    for r in inc:
        L.append(f"- {r['run']}: {r['scored']}/{r['total']} swarms")

open(out_md, "w").write("\n".join(L) + "\n")
print(f"wrote {out_md} ({len(L)} lines, {n_done}/12 runs complete)")
