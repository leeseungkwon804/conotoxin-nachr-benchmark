#!/usr/bin/env python3
"""
The three round-2 supplementary tables that needed pod data, now mirrored under pod_mirror_r2r3/.

  Table S19  Native interface inventory and per-method reproduction at anchor A2 (7EGR)
  Table S20  CAPRI class and DockQ components for every method at both anchors
  Table S22  Extracellular-domain control ladder at anchor A1, per run

S18 and S21 are built by make_supp_tables_round2.py from data that was already local; this script
is separate only because its inputs arrived later, not because the two sets differ in kind.

Usage: python3 scripts/make_supp_tables_round3.py
"""
import csv
import glob
import json
import os

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
MIR = os.path.join(ROOT, "pod_mirror_r2r3")
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


def write_md(path, caption, head, rows, footnotes):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(caption + "\n\n")
        fh.write("| " + " | ".join(head) + " |\n")
        fh.write("|" + "---|" * len(head) + "\n")
        for r in rows:
            fh.write("| " + " | ".join(str(c) for c in r) + " |\n")
        for line in footnotes:
            fh.write(line + "\n")
    with open(path[:-3] + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        w.writerows(rows)


# ------------------------------------------------------------------ Table S19

def table_s19():
    ifc = os.path.join(MIR, "round2", "iface")
    nat = {(r["rec_chain"], r["rec_resnum"], r["lig_resnum"]): r
           for r in csv.DictReader(open(os.path.join(ifc, "native_contacts_A2_7EGR.csv")))}
    rows = list(csv.DictReader(open(os.path.join(ifc, "recovery_contacts_A2_7EGR.csv"))))
    inter = list(csv.DictReader(open(os.path.join(ifc, "native_interactions_A2_7EGR.csv"))))

    rungs = [("BLIND", ["A2_7EGR_BLIND_rep1", "A2_7EGR_BLIND_rep2"]),
             ("XDON", ["A2_7EGR_XDON_rep1", "A2_7EGR_XDON_rep2"]),
             ("NATIVE", ["A2_7EGR_NATIVE_rep1", "A2_7EGR_NATIVE_rep2"])]
    poses = ["best_rmsd", "top1_score"]          # the two a user receives; matches Figure 5a's rule

    rows.sort(key=lambda r: (LOOP_ORDER.index(r["loop"]) if r["loop"] in LOOP_ORDER
                             else len(LOOP_ORDER), int(r["rec_resnum"]), int(r["lig_resnum"])))

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
        cells = [LOOP_LABEL.get(r["loop"], r["loop"]),
                 "%s(%s)" % (res(r["rec_resname"], r["rec_resnum"]), r["rec_chain"]),
                 res(r["lig_resname"], r["lig_resnum"]),
                 r["n_atom_pairs"], r["min_dist_A"],
                 ", ".join(kinds) or "—",
                 "%d" % sum(int(r["boltz2_m%d" % i]) for i in range(5)),
                 "%d" % sum(int(r["boltz1_m%d" % i]) for i in range(5))]
        for _, runs in rungs:
            cells.append("yes" if any(int(r["%s|%s" % (run, p)]) for run in runs for p in poses)
                         else "no")
        out.append(cells)

    head = ["Loop", "Receptor residue", "RgIA residue", "Atom pairs < 5 Å", "Min. distance (Å)",
            "Interaction type(s)", "Boltz-2 (of 5)", "Boltz-1 (of 5)", "BLIND", "XDON", "NATIVE"]

    tot = {k: sum(int(v[k]) for v in nat.values())
           for k in ("n_hbond", "n_saltbridge", "n_hydrophobic", "n_aromatic", "n_cation_pi")}
    sg = [r for r in inter if r["type"] == "hbond"
          and (r["rec_atom"] == "SG" or r["lig_atom"] == "SG")]
    sg_acc = [r for r in sg if (r["detail"] == "rec_donor" and r["lig_atom"] == "SG")
              or (r["detail"] == "lig_donor" and r["rec_atom"] == "SG")]
    b2 = [sum(int(r["boltz2_m%d" % i]) for r in rows) for i in range(5)]
    b1 = [sum(int(r["boltz1_m%d" % i]) for r in rows) for i in range(5)]
    per_rung = {lab: sum(1 for r in rows
                         if any(int(r["%s|%s" % (run, p)]) for run in runs for p in poses))
                for lab, runs in rungs}

    cap = (
        "**Table S19.** Native receptor–ligand residue contacts at anchor **A2 · 7EGR** "
        "(*Aplysia* acetylcholine-binding protein + α-conotoxin RgIA) and their reproduction by "
        "each method — the companion to Table S18, built identically so the two anchors can be "
        "read side by side. %d contacts in total, one per row, ordered by binding-site loop and "
        "then by receptor residue number. Contact definition, residue numbering, rung-column "
        "criterion and interaction-type criteria are exactly as in Table S18. A2 additionally "
        "carries the Boltz-1 temporal control, whose five samples are given their own column. "
        "Machine-readable companion: `Table_S19_interface_inventory_A2.csv`." % len(rows))

    foot = [
        "",
        "[a] Native interaction totals at this interface: %d hydrogen bonds, %d ionic, "
        "%d hydrophobic, %d aromatic, %d cation–π. %d hydrogen bonds involve a cysteine sulfur, "
        "of which %d have that sulfur as the acceptor; as at A1, main-text hydrogen-bond counts "
        "are quoted both with and without those."
        % (tot["n_hbond"], tot["n_saltbridge"], tot["n_hydrophobic"], tot["n_aromatic"],
           tot["n_cation_pi"], len(sg), len(sg_acc)),
        "",
        "[b] Contacts reproduced, totalled down each column: Boltz-2 samples 0–4, %s of %d; "
        "Boltz-1 samples 0–4, %s of %d; BLIND %d, XDON %d, NATIVE %d."
        % ("/".join(str(v) for v in b2), len(rows), "/".join(str(v) for v in b1), len(rows),
           per_rung["BLIND"], per_rung["XDON"], per_rung["NATIVE"]),
        "",
        "[c] This anchor is where the comparison is not one-directional. Four contacts are missed "
        "by all five Boltz-2 samples; the closest of them, the Thr53(G)–Arg409 hydrogen bond at "
        "2.58 \u00c5 — the third-closest approach in the whole interface — is reproduced by all five "
        "Boltz-1 samples and by docked poses at two rungs. One contact, Asp176(G)–Cys412, is "
        "reached by no model of any kind.",
        "",
        "[d] A2 is the leakage-exposed anchor (Section 2.4, Table S14). Boltz-1, trained before "
        "the relevant depositions, is reported alongside precisely so that this table is not read "
        "as evidence about Boltz-2 alone.",
    ]

    write_md(os.path.join(SUP, "Table_S19_interface_inventory_A2.md"), cap, head, out, foot)
    print("  Table S19: %d contact rows" % len(out))


# ------------------------------------------------------------------ Table S20

CAPRI_ORDER = {"High": 0, "Medium": 1, "Acceptable": 2, "Incorrect": 3}


def table_s20():
    src = list(csv.DictReader(open(os.path.join(MIR, "round2", "dockq",
                                                "CAPRI_DockQ_all_methods.csv"))))
    std = {(r["method"], r["anchor"], r["run"], r["pose"]): r
           for r in src if r["peptide_mode"] == "False"}
    pep = {(r["method"], r["anchor"], r["run"], r["pose"]): r
           for r in src if r["peptide_mode"] == "True"}

    def sort_key(k):
        method, anchor, run, pose = k
        mrank = {"LightDock": 0, "Boltz-1": 1, "Boltz-2": 2}.get(method, 3)
        return (anchor, mrank, run, pose)

    out = []
    for k in sorted(std, key=sort_key):
        r = std[k]
        p = pep.get(k, {})
        label = r["run"] if r["method"] == "LightDock" else "%s %s" % (r["method"], r["pose"])
        out.append([
            r["anchor"], r["method"], label,
            {"best_rmsd": "best RMSD", "top1_score": "top 1 by score"}.get(r["pose"], r["pose"]),
            "%s/%s" % (r["n_shared_contacts"], r["n_native_contacts"]),
            "%.3f" % float(r["fnat"]), "%.3f" % float(r["fnonnat"]),
            "%.3f" % float(r["LRMS_A"]), "%.3f" % float(r["iRMS_A"]),
            "%.4f" % float(r["DockQ"]), r["CAPRI_class"],
            p.get("CAPRI_class_peptide_bands", "—"),
        ])

    head = ["Anchor", "Method", "Run / model", "Pose selected by", "Contacts shared / native",
            "fnat", "fnon-nat", "L-RMSD (Å)", "I-RMSD (Å)", "DockQ",
            "CAPRI class (standard bands)", "CAPRI class (peptide bands)"]

    cap = (
        "**Table S20.** CAPRI classification and DockQ components for every model evaluated in "
        "this work, at both anchors. %d rows: the twelve LightDock runs contribute two poses each "
        "(the best-RMSD pose in the run and the pose its scoring function ranked first), Boltz-2 "
        "contributes five samples per anchor, Boltz-1 five at A2, and the rigid-placement ceiling "
        "is listed for reference. fnat, L-RMSD and I-RMSD follow the CAPRI definitions and are "
        "combined into DockQ by its published definition [18]; the implementation reproduces "
        "DockQ v2.1.3 across every row (Methods §4.8). Both band systems are reported: the "
        "standard protein–protein CAPRI bands, and the stricter peptide-target bands appropriate "
        "to a twelve-to-thirteen-residue ligand. The choice of band changes the reading, which is "
        "why neither is given alone. This table is the source of Figure 6a. Machine-readable "
        "companion: `Table_S20_capri_dockq_all_methods.csv`." % len(out))

    def rng(method, anchor, pose=None):
        v = [float(r["DockQ"]) for kk, r in std.items()
             if kk[0] == method and kk[1] == anchor and (pose is None or kk[3] == pose)]
        return (min(v), max(v), len(v)) if v else (float("nan"),) * 2 + (0,)

    foot = [
        "",
        "[a] DockQ ranges by method, standard bands. A1: LightDock best-RMSD %.3f–%.3f, "
        "LightDock top-1 %.3f–%.3f, Boltz-2 %.3f–%.3f. A2: LightDock best-RMSD %.3f–%.3f, "
        "LightDock top-1 %.3f–%.3f, Boltz-2 %.3f–%.3f, Boltz-1 %.3f–%.3f."
        % (rng("LightDock", "A1", "best_rmsd")[:2] + rng("LightDock", "A1", "top1_score")[:2]
           + rng("Boltz-2", "A1")[:2] + rng("LightDock", "A2", "best_rmsd")[:2]
           + rng("LightDock", "A2", "top1_score")[:2] + rng("Boltz-2", "A2")[:2]
           + rng("Boltz-1", "A2")[:2]),
        "",
        "[b] The receptor is held at its deposited coordinates throughout the docking runs, and "
        "the receptor file used for docking is byte-identical to the one used to build the native "
        "complex. The CAPRI receptor superposition is therefore the identity (Kabsch fit RMSD "
        "0.0000 Å in every row), and CAPRI L-RMSD coincides with the un-superposed ligand "
        "backbone RMSD. This is a property of this benchmark's construction, not a general one.",
        "",
        "[c] The best-RMSD pose of a run is an oracle selection — it requires knowing the answer. "
        "The top-1-by-score row is what a user actually receives. Reading the two together is the "
        "point of the table: at both anchors the oracle row is one to three CAPRI classes better "
        "than the delivered row.",
    ]

    write_md(os.path.join(SUP, "Table_S20_capri_dockq_all_methods.md"), cap, head, out, foot)
    print("  Table S20: %d model rows" % len(out))


# ------------------------------------------------------------------ Table S22

def table_s22():
    runs = []
    for f in sorted(glob.glob(os.path.join(MIR, "ecd_rmsd", "*.json"))):
        runs.append(json.load(open(f)))
    ecd = json.load(open(os.path.join(MIR, "ecd_definition.json")))

    order = {"BLIND": 0, "XDON": 1, "NATIVE": 2}
    runs.sort(key=lambda d: (order[d["run"].split("_")[1]], d["run"]))

    out = []
    for d in runs:
        rung = d["run"].split("_")[1]
        out.append([rung, d["run"].split("_")[2], d["swarms_total"], "{:,}".format(d["poses_scored"]),
                    d["restraints_R"], d["restraints_L"],
                    "%.3f" % d["best_rmsd"], "%.3f" % d["rmsd_of_top1_by_score"],
                    d["best_rmsd_swarm"], d["best_rmsd_score_rank"]])

    head = ["Rung", "GSO replicate", "Swarms", "Poses scored", "Receptor restraints",
            "Ligand restraints", "Best ligand RMSD (Å)", "RMSD of top-1-by-score pose (Å)",
            "Swarm holding best pose", "Score rank of best pose"]

    means = {}
    for rung in ("BLIND", "XDON", "NATIVE"):
        v = [d["best_rmsd"] for d in runs if d["run"].split("_")[1] == rung]
        means[rung] = sum(v) / len(v)

    cap = (
        "**Table S22.** The three-rung disclosure ladder repeated at anchor A1 on a truncated "
        "receptor carrying the extracellular domain only, one row per run. This is the control "
        "for the difference in receptor extent between the two anchors: A1 presents both the "
        "extracellular and transmembrane domains while A2 presents an extracellular domain only, "
        "so their ladders are not otherwise obtained on comparable surfaces (Methods §4.10). "
        "Toxin, protocol, scoring function, swarm construction and step count are unchanged from "
        "the full-receptor A1 runs; only the receptor file differs. The domain boundary was not "
        "asserted from sequence annotation: A1's two docked chains were structurally aligned to "
        "A2's extracellular-domain chains by the same routine that builds the cross-donated "
        "restraints, and the retained region is the aligned one (%d and %d residues aligned on "
        "the (+) and (\u2212) faces, refined RMSD %.2f and %.2f \u00c5, sequence identity %.0f%% and "
        "%.0f%% over the alignment). The truncated receptor retains %s of the %s heavy atoms of "
        "the full construct, against %s in A2 \u2014 so the control brings A1 to A2's scale rather "
        "than to an arbitrary one. Machine-readable companion: "
        "`Table_S22_ecd_control_ladder.csv`."
        % (ecd["pairs"][0]["aligned"], ecd["pairs"][1]["aligned"],
           ecd["pairs"][0]["refined_rmsd"], ecd["pairs"][1]["refined_rmsd"],
           100 * ecd["pairs"][0]["seq_ident_over_aligned"],
           100 * ecd["pairs"][1]["seq_ident_over_aligned"],
           "{:,}".format(ecd["atoms_a1_ecd"]), "{:,}".format(ecd["atoms_a1_full"]),
           "{:,}".format(ecd["atoms_a2_for_scale"])))

    foot = [
        "",
        "[a] Rung means of the best ligand RMSD: BLIND %.2f Å, XDON %.2f Å, NATIVE %.2f Å. "
        "The ladder improves monotonically with disclosure here, the opposite of the "
        "full-receptor A1 ladder, and both orderings are replicate-consistent. The sign of the "
        "disclosure effect therefore belongs to the receptor construct rather than to the ligand "
        "class — the finding this table exists to establish."
        % (means["BLIND"], means["XDON"], means["NATIVE"]),
        "",
        "[b] Swarm counts differ between rungs because restraints restrict where swarms may be "
        "placed; this is LightDock's own behaviour and is not a setting that was varied. The "
        "poses-scored column shows the search budget actually spent in each run, so the rungs can "
        "be compared on effort as well as on outcome.",
        "",
        "[c] The top-1-by-score column is included to show that the scoring failure documented at "
        "the full receptor (Section 2.3) is unchanged by truncation: no rung delivers a "
        "near-native pose as its first-ranked result, whichever way the best-pose ladder runs.",
    ]

    write_md(os.path.join(SUP, "Table_S22_ecd_control_ladder.md"), cap, head, out, foot)
    print("  Table S22: %d run rows (BLIND %.2f / XDON %.2f / NATIVE %.2f Å)"
          % (len(out), means["BLIND"], means["XDON"], means["NATIVE"]))


if __name__ == "__main__":
    table_s19()
    table_s20()
    table_s22()
