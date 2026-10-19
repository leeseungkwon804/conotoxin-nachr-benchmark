#!/usr/bin/env python3
"""
Supplementary tables for the round-2 additions, built from local primary data only.

Two of the four round-2 tables can be built here; the other two (the A2 half of the interface
inventory, and the CAPRI/DockQ per-run panel) need per-pose arrays that live on the pod and are
not mirrored locally, so they are not faked here — see the TODO note printed at the end.

  Table S18  Native interface inventory and per-method reproduction at anchor A1
             <- manuscript/figures_v3/panels/{native,recovery}_contacts_A1_9NX0.csv
  Table S21  Co-folding confidence by condition, Boltz-2 and Chai-1, with exact permutation tests
             <- manuscript/figures_v3/panels/panel_data.json

The permutation tests are exact rather than sampled: 5 vs 5 has only C(10,5)/2 = 126 distinct
labellings, so the whole null is enumerable and the smallest attainable two-sided p is 1/126.

Usage: python3 scripts/make_supp_tables_round2.py
"""
import csv
import itertools
import json
import os
import statistics

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
PAN = os.path.join(ROOT, "manuscript", "figures_v3", "panels")
SUP = os.path.join(ROOT, "manuscript", "supplementary")

LOOP_ORDER = ["loopA", "loopB", "loopC", "loopD", "loopE", "loopF", "outside_A-F"]
LOOP_LABEL = {"loopA": "A", "loopB": "B", "loopC": "C", "loopD": "D", "loopE": "E",
              "loopF": "F", "outside_A-F": "outside A–F"}
AA3 = {"ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys", "GLN": "Gln",
       "GLU": "Glu", "GLY": "Gly", "HIS": "His", "ILE": "Ile", "LEU": "Leu", "LYS": "Lys",
       "MET": "Met", "PHE": "Phe", "PRO": "Pro", "SER": "Ser", "THR": "Thr", "TRP": "Trp",
       "TYR": "Tyr", "VAL": "Val"}


def res(name, num):
    return "%s%s" % (AA3.get(name, name.title()), num)


# ------------------------------------------------------------------ Table S18

def table_s18():
    nat = {(r["rec_chain"], r["rec_resnum"], r["lig_resnum"]): r
           for r in csv.DictReader(open(os.path.join(PAN, "native_contacts_A1_9NX0.csv")))}
    rows = list(csv.DictReader(open(os.path.join(PAN, "recovery_contacts_A1_9NX0.csv"))))

    rungs = [("BLIND", ["A1_9NX0_BLIND_rep1", "A1_9NX0_BLIND_rep2"]),
             ("XDON", ["A1_9NX0_XDON_rep1", "A1_9NX0_XDON_rep2"]),
             ("NATIVE", ["A1_9NX0_NATIVE_rep1", "A1_9NX0_NATIVE_rep2"])]
    # The two poses a user actually receives. Figure 5a uses exactly this pair, and the table is
    # the figure's source, so it must not quietly pool a third. The oracle best-fnat pose is
    # reported separately in footnote [d] instead of being folded in here.
    poses = ["best_rmsd", "top1_score"]

    rows.sort(key=lambda r: (LOOP_ORDER.index(r["loop"]), int(r["rec_resnum"]),
                             int(r["lig_resnum"])))

    out = []
    for r in rows:
        n = nat.get((r["rec_chain"], r["rec_resnum"], r["lig_resnum"]), {})
        kinds = []
        for key, lab in (("n_hbond", "H-bond"), ("n_saltbridge", "ionic"),
                         ("n_hydrophobic", "hydrophobic"), ("n_aromatic", "aromatic"),
                         ("n_cation_pi", "cation–π")):
            c = int(n.get(key, 0) or 0)
            if c:
                kinds.append("%s ×%d" % (lab, c) if c > 1 else lab)
        b2 = sum(int(r["boltz2_m%d" % i]) for i in range(5))
        cells = [LOOP_LABEL[r["loop"]],
                 "%s(%s)" % (res(r["rec_resname"], r["rec_resnum"]), r["rec_chain"]),
                 res(r["lig_resname"], r["lig_resnum"]),
                 r["n_atom_pairs"], r["min_dist_A"],
                 ", ".join(kinds) or "—",
                 "%d" % b2]
        for _, runs in rungs:
            hit = any(int(r["%s|%s" % (run, p)]) for run in runs for p in poses)
            cells.append("yes" if hit else "no")
        out.append(cells)

    head = ["Loop", "Receptor residue", "ImI residue", "Atom pairs < 5 Å", "Min. distance (Å)",
            "Interaction type(s)", "Boltz-2 samples (of 5)", "BLIND", "XDON", "NATIVE"]

    tot = {k: sum(int(v[k]) for v in nat.values())
           for k in ("n_hbond", "n_saltbridge", "n_hydrophobic", "n_aromatic", "n_cation_pi")}

    # Hydrogen bonds with a cysteine partner are separated out: at this interface every cysteine
    # involved is disulfide-bonded (the receptor's loop C vicinal pair and ImI's framework), so
    # the geometry is real but the chemistry is doubtful. Counted, not asserted.
    cys_rows = [v for v in nat.values()
                if int(v["n_hbond"]) and "CYS" in (v["rec_resname"], v["lig_resname"])]
    cys_hb = sum(int(v["n_hbond"]) for v in cys_rows)
    cys_desc = ["%s(%s)–%s" % (res(v["rec_resname"], v["rec_resnum"]), v["rec_chain"],
                               res(v["lig_resname"], v["lig_resnum"])) for v in cys_rows]
    per_sample = [sum(int(r["boltz2_m%d" % i]) for r in rows) for i in range(5)]
    per_rung, oracle = {}, {}
    for lab, runs in rungs:
        per_rung[lab] = sum(1 for r in rows
                            if any(int(r["%s|%s" % (run, p)]) for run in runs for p in poses))
        oracle[lab] = sum(1 for r in rows
                          if any(int(r["%s|%s" % (run, p)]) for run in runs
                                 for p in poses + ["best_fnat"]))

    cap = (
        "**Table S18.** Native receptor–ligand residue contacts at anchor **A1 · 9NX0** "
        "(human α7 nAChR + α-conotoxin ImI) and their reproduction by each method. "
        "One row per native contact, ordered by binding-site loop and then by receptor residue "
        "number; %d contacts in total. A residue pair counts as in contact when any heavy-atom "
        "pair falls within 5.0 Å, the same definition used for fnat throughout (Methods §4.4, "
        "§4.7). Residue numbering is that of the deposited coordinates. The Boltz-2 column gives "
        "how many of the five diffusion samples reproduced the contact. The three rung columns "
        "record whether **any** reported pose of that rung reproduced it, pooling both GSO "
        "replicates and both user-facing poses per run (best-RMSD and top-1-by-score) — "
        "a deliberately generous criterion, so a \"no\" means no pose the search reported "
        "contained that contact. Interaction types follow the geometric criteria in footnote [a]. "
        "This table is the source of Figure 5a,b. Machine-readable companion: "
        "`Table_S18_interface_inventory_A1.csv`."
        % len(rows))

    foot = [
        "",
        "[a] Geometric criteria, applied to the deposited coordinates. **H-bond**: donor–acceptor "
        "heavy-atom distance ≤ 3.5 Å with donor-H···acceptor angle ≥ 120°, hydrogens placed by the "
        "structure-preparation step of Methods §4.1. **Ionic**: carboxylate/guanidinium or "
        "ammonium nitrogen–oxygen pair ≤ 4.0 Å. **Hydrophobic**: carbon–carbon pair ≤ 4.5 Å "
        "between apolar side-chain atoms. **Aromatic** (π-stacking) and **cation–π** were "
        "evaluated on ring-centroid geometry; neither occurs at this interface, and both columns "
        "are therefore zero throughout.",
        "",
        "[b] Native interaction totals at this interface: %d hydrogen bonds, %d ionic, "
        "%d hydrophobic, %d aromatic, %d cation–π. %d of the hydrogen bonds, spread over %d "
        "contacts (%s), have a disulfide-bonded cysteine as one partner; these are geometrically "
        "identified but chemically questionable, so hydrogen-bond counts in the main text are "
        "quoted both with and without them (%d with, %d without)."
        % (tot["n_hbond"], tot["n_saltbridge"], tot["n_hydrophobic"],
           tot["n_aromatic"], tot["n_cation_pi"], cys_hb, len(cys_rows),
           "; ".join(cys_desc), tot["n_hbond"], tot["n_hbond"] - cys_hb),
        "",
        "[c] Contacts reproduced, totalled down each column: Boltz-2 samples 0–4 reproduced "
        "%s of %d respectively; BLIND %d, XDON %d, NATIVE %d."
        % ("/".join(str(v) for v in per_sample), len(rows),
           per_rung["BLIND"], per_rung["XDON"], per_rung["NATIVE"]),
        "",
        "[d] The rung columns pool across poses and replicates and so overstate what any single "
        "docked pose achieved. No individual docked pose at this anchor exceeded fnat 0.245 "
        "(13 of %d contacts). Adding the oracle best-fnat pose of each run — a pose no user could "
        "select, since it is chosen by knowing the answer — raises the pooled totals only to "
        "BLIND %d, XDON %d, NATIVE %d." % (len(rows), oracle["BLIND"], oracle["XDON"],
                                           oracle["NATIVE"]),
    ]

    write_md(os.path.join(SUP, "Table_S18_interface_inventory_A1.md"), cap, head, out, foot)
    with open(os.path.join(SUP, "Table_S18_interface_inventory_A1.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerows(out)
    print("  Table S18: %d contact rows" % len(out))


# ------------------------------------------------------------------ Table S21

def perm_p(a, b):
    """Exact two-sided permutation p on the difference of means, 5 vs 5 -> 126 labellings."""
    pool = list(a) + list(b)
    obs = abs(statistics.mean(a) - statistics.mean(b))
    n = len(a)
    seen, hits, total = set(), 0, 0
    for idx in itertools.combinations(range(len(pool)), n):
        comp = tuple(sorted(set(range(len(pool))) - set(idx)))
        key = min((idx, comp))
        if key in seen:
            continue
        seen.add(key)
        total += 1
        g1 = [pool[i] for i in idx]
        g2 = [pool[i] for i in comp]
        if abs(statistics.mean(g1) - statistics.mean(g2)) >= obs - 1e-12:
            hits += 1
    return hits / total, total


def table_s21():
    d = json.load(open(os.path.join(PAN, "panel_data.json")))

    # (block, condition label, boltz key or None, chai key or None)
    species = [
        ("Vc1.1 · full mature receptor (ECD+TMD)", "Vc1.1 / full mature ECD+TMD", None),
        ("Vc1.1 · ECD only", "Vc1.1 / ECD only", "VC11"),
        ("Vc1.1 · ECD only, explicit disulfides", "Vc1.1 + explicit SS / ECD only", None),
        ("RgIA · ECD only", "RgIA / ECD only", "RGIA"),
    ]
    e2 = [("ImI (true binder)", "E2_TRUE_IMI"),
          ("ImI, non-cysteine positions shuffled", "E2_SHUF"),
          ("RgIA (α9α10 toxin, off-target here)", "E2_NB_RGIA"),
          ("MVIIA (Ca_V2.2 toxin, off-target here)", "E2_NB_MVIIA")]

    rows = []
    for label, bkey, ckey in species:
        for model, vals in (("Boltz-2", d["boltz_species"].get(bkey)),
                            ("Chai-1", None if not ckey else
                             {"human": d["chai"]["E1_HUM_ECD_%s" % ckey] if ckey == "VC11"
                              else d["chai"]["E1b_HUM_ECD_%s" % ckey],
                              "rat": d["chai"]["E1_RAT_ECD_%s" % ckey] if ckey == "VC11"
                              else d["chai"]["E1b_RAT_ECD_%s" % ckey]})):
            if not vals:
                continue
            hu, ra = vals["human"], vals["rat"]
            p, ntot = perm_p(hu, ra)
            rows.append(["α9α10 orthologue", label, model,
                         "%.4f–%.4f" % (min(hu), max(hu)), "%.4f" % statistics.mean(hu),
                         "%.4f–%.4f" % (min(ra), max(ra)), "%.4f" % statistics.mean(ra),
                         "%+.4f" % (statistics.mean(hu) - statistics.mean(ra)),
                         "%.4f" % p])

    true_b2 = d["boltz_e2"]["E2_TRUE_IMI"]
    true_ch = d["chai"]["E2_TRUE_IMI"]
    for label, key in e2:
        for model, src, true in (("Boltz-2", d["boltz_e2"], true_b2),
                                 ("Chai-1", d["chai"], true_ch)):
            v = src[key]
            cell = ["binder vs non-binder", label, model,
                    "%.4f–%.4f" % (min(v), max(v)), "%.4f" % statistics.mean(v), "—", "—",
                    "%+.4f" % (statistics.mean(v) - statistics.mean(true)),
                    "—" if key == "E2_TRUE_IMI" else "%.4f" % perm_p(v, true)[0]]
            rows.append(cell)

    head = ["Experiment", "Condition", "Model", "Human range", "Human mean",
            "Rat range", "Rat mean", "Δ mean (vs comparator)", "Exact p"]

    cap = (
        "**Table S21.** Co-folding confidence by condition in two architecturally independent "
        "models, all five diffusion samples per condition. Boltz-2 values are its composite "
        "confidence; Chai-1 values are its aggregate score. The two scales are **not** "
        "comparable to one another and are never pooled — each model is read only against its own "
        "conditions. The upper block is the α9α10 orthologue experiment (comparator: the rat "
        "construct of the same row); the lower block is the binder/non-binder experiment at the "
        "A1 α7 site (comparator: the true binder ImI in the same model). *p* is a two-sided "
        "exact permutation test on the difference of means over all %d distinct 5-vs-5 "
        "labellings, so the smallest attainable value is 1/%d = %.4f; no multiplicity correction "
        "is applied because no comparison here is used to support a positive claim. This table is "
        "the source of Figure 7. Machine-readable companion: "
        "`Table_S21_cofolding_confidence.csv`." % (126, 126, 1 / 126))

    foot = [
        "",
        "[a] Both models were run with matched settings (5 diffusion samples, 200 timesteps, "
        "3 recycles, server-generated MSAs), so a difference between them is not a difference in "
        "sampling budget. Inference provenance — software versions, weight checksums and run "
        "dates — is recorded in Table S12; the raw per-sample outputs are in File S1.",
        "",
        "[b] The measured pharmacology this experiment is set against runs the other way: Vc1.1 "
        "(ACV1) is roughly two orders of magnitude weaker at human α9α10 than at rat. A "
        "confidence gap that is small, model-dependent in sign, or of the wrong sign therefore "
        "does not recover the species difference, whatever its *p*.",
        "",
        "[c] The receptor construct matters for the Boltz-2 rows: the same toxin and the same "
        "receptor sequences give different gaps depending on whether the transmembrane domain is "
        "present, which is why the full-receptor and ECD-only rows are reported separately "
        "rather than pooled.",
    ]

    write_md(os.path.join(SUP, "Table_S21_cofolding_confidence.md"), cap, head, rows, foot)
    with open(os.path.join(SUP, "Table_S21_cofolding_confidence.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerows(rows)
    print("  Table S21: %d condition rows" % len(rows))


# ------------------------------------------------------------------ output

def write_md(path, caption, head, rows, footnotes):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(caption + "\n\n")
        fh.write("| " + " | ".join(head) + " |\n")
        fh.write("|" + "---|" * len(head) + "\n")
        for r in rows:
            fh.write("| " + " | ".join(str(c) for c in r) + " |\n")
        for line in footnotes:
            fh.write(line + "\n")


if __name__ == "__main__":
    table_s18()
    table_s21()
    print("\nStill to build, and blocked on pod access (per-pose arrays are not mirrored locally):")
    print("  Table S19  interface inventory at A2 · 7EGR   <- round2 recovery arrays")
    print("  Table S20  CAPRI class and DockQ components, per run and per sample")
    print("  Table S22  extracellular-domain control ladder at A1, per run")
