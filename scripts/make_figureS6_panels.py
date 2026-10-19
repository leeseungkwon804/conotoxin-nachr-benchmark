#!/usr/bin/env python3
"""
Figure S6 raw panels -- residue-level interface view, deposited vs Boltz-2, both
anchors. Rendered headless with PyMOL 3.x.

Run:  /opt/homebrew/bin/pymol -cq scripts/make_figureS6_panels.py

This follows the receptor-CA-superposition method established in
scripts/make_structure_panels.py::panel_4a() for the prediction-overlay panels --
homomeric receptor chains are re-paired by exhaustive permutation, keeping the
mapping that minimises CA RMSD, because a naive cmd.align() mis-paired chains
here before and inflated RMSD from 0.8 A to 17.7 A (see that script's docstring
and pod_mirror_n2/eval_cofold.py's docstring, which records the same failure mode
independently: 0.825 A -> 18.271 A on A2).

Model selection (median by interface all-atom ligand RMSD, not by confidence):
  A1 9NX0  -> pod_mirror_boltz/eval_ms/summary.csv           -> model 1
  A2 7EGR  -> pod_mirror_n2/eval_7EGR_b1_vs_b2.json, Boltz-2  -> model 4
  (verified against those two files directly; A2's Boltz-2 rows in that JSON are
  the SAME 5 poses as pod_mirror_boltz/eval_ms's A2_7EGR set, index-permuted --
  eval_7EGR_b1_vs_b2.json is the canonical source per project rule, so the A2
  prediction PDB is read from pod_mirror_n2/outputs_b2_7EGR/, not pod_mirror_boltz/)

Native/deposited coordinates: pod_mirror/structures/lightdock_inputs/ (deposited
frame; pod_mirror_n2/native/ holds byte-identical copies of the same two files).

Contact definition used to pick which residues get sticks + labels: receptor atom
within 5.0 A of a ligand atom (the same cutoff eval_cofold.py and eval_iface.py
use for fnat).
"""
import itertools
import os

from pymol import cmd

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
LD = os.path.join(ROOT, "pod_mirror/structures/lightdock_inputs")
OUT = os.path.join(ROOT, "manuscript", "figures_v2", "panels")
os.makedirs(OUT, exist_ok=True)

W, H, DPI = 1700, 1700, 400

C_PLUS = "0x6FA8DC"    # principal (+) subunit, light blue
C_MINUS = "0xB4B4B4"   # complementary (-) subunit, grey
C_TOX = "0xE69F00"     # deposited toxin, orange
C_TOX2 = "0xCC79A7"    # Boltz-2 predicted toxin, purple
C_IFACE = "0x3D3D3D"   # interface receptor side-chain carbons, dark grey

ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def label_expr(resn, resi):
    return "%s%s" % (ONE.get(resn, resn), resi)


def base_setup():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("stick_radius", 0.15)
    cmd.set("valence", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.15)
    cmd.set("label_size", 16)
    cmd.set("label_font_id", 7)  # Sans Bold
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_shadow_mode", 0)
    cmd.set("dash_gap", 0.35)
    cmd.set("dash_length", 0.18)
    cmd.set("dash_width", 2.0)
    cmd.set("dash_color", "red")
    cmd.set("dash_radius", 0.035)


def render(name):
    cmd.png(os.path.join(OUT, name + ".png"), width=W, height=H, dpi=DPI, ray=1)
    print("  rendered", name)


def label_residues(sel):
    """Label the CA of every residue in sel with one-letter code + resi."""
    model = cmd.get_model(sel + " and name CA")
    for at in model.atom:
        cmd.label("%s and chain %s and resi %s and name CA" % (sel, at.chain, at.resi),
                   '"%s"' % label_expr(at.resn, at.resi))


def show_interface(rec_obj, lig_obj, rec_sel_extra="", lig_sel_extra=""):
    """rec_obj / lig_obj: object (or object+chain) selections already loaded.
    Selects receptor residues within 5 A of the ligand, shows sticks + labels on
    both sides, and draws dashed polar contacts between them."""
    rec_sel = rec_obj + (" and (%s)" % rec_sel_extra if rec_sel_extra else "")
    lig_sel = lig_obj + (" and (%s)" % lig_sel_extra if lig_sel_extra else "")
    cmd.select("iface_rec", "byres (%s within 5.0 of (%s))" % (rec_sel, lig_sel))
    n = cmd.count_atoms("iface_rec and name CA")
    cmd.show("sticks", "iface_rec and not (name C+N+O)")
    cmd.color(C_IFACE, "iface_rec and elem C")
    cmd.show("sticks", lig_sel)
    label_residues("iface_rec")
    label_residues(lig_sel)
    cmd.distance("hb", "iface_rec and (elem N+O)", "(%s) and (elem N+O)" % lig_sel,
                 3.6, mode=2)
    cmd.hide("labels", "hb")
    cmd.delete("iface_rec")
    return n


# --------------------------------------------------------------- deposited column
def deposited_panel(tag, rec_pdb, lig_pdb, lig_color):
    base_setup()
    cmd.load(os.path.join(LD, rec_pdb), "rec")
    cmd.load(os.path.join(LD, lig_pdb), "lig")
    cmd.hide("everything")
    cmd.show("cartoon", "rec")
    chains = cmd.get_chains("rec")
    cmd.color(C_PLUS, "rec and chain %s" % chains[0])
    cmd.color(C_MINUS, "rec and chain %s" % chains[1])
    cmd.set("cartoon_transparency", 0.45, "rec")
    cmd.show("cartoon", "lig")
    cmd.color(lig_color, "lig")
    n = show_interface("rec", "lig")
    print("  %s (deposited): %d receptor residues within 5 A of the toxin" % (tag, n))
    cmd.orient("lig")
    cmd.zoom("lig", 3.2)
    render(tag)


# --------------------------------------------------------------- prediction column
def prediction_panel(tag, rec_pdb, lig_pdb, pred_pdb, ligch, lig_color):
    if not os.path.exists(pred_pdb):
        print("  SKIP %s (prediction not found: %s)" % (tag, pred_pdb))
        return
    base_setup()
    # load the deposited complex too, ONLY to fix an identical camera to the
    # deposited-column panel for this anchor -- it is hidden, never rendered.
    cmd.load(os.path.join(LD, rec_pdb), "rec")
    cmd.load(os.path.join(LD, lig_pdb), "lig")
    cmd.load(pred_pdb, "pred")
    cmd.hide("everything")
    cmd.select("pred_lig", "pred and chain %s" % ligch)

    pred_ch = [c for c in cmd.get_chains("pred") if c != ligch]
    nat_ch = cmd.get_chains("rec")
    best = None
    for perm in itertools.permutations(pred_ch):
        cmd.create("probe", "pred")
        pairs = list(zip(perm, nat_ch))
        m = cmd.pair_fit(*[x for pc, nc in pairs
                           for x in ("probe and chain %s and name CA" % pc,
                                     "rec and chain %s and name CA" % nc)])
        cmd.delete("probe")
        if best is None or m < best[0]:
            best = (m, pairs)
    rms, pairs = best
    cmd.pair_fit(*[x for pc, nc in pairs
                   for x in ("pred and chain %s and name CA" % pc,
                             "rec and chain %s and name CA" % nc)])
    print("  %s (Boltz-2): receptor CA superposition RMSD %.3f A  (chain map %s)"
          % (tag, rms, {p_: n_ for p_, n_ in pairs}))

    cmd.hide("everything")
    cmd.show("cartoon", "pred and not chain %s" % ligch)
    for pc, nc in pairs:
        cmd.color(C_PLUS if nc == nat_ch[0] else C_MINUS, "pred and chain %s" % pc)
    cmd.set("cartoon_transparency", 0.45, "pred and not chain %s" % ligch)
    cmd.show("cartoon", "pred_lig")
    cmd.color(lig_color, "pred_lig")
    n = show_interface("pred and not chain %s" % ligch, "pred_lig")
    print("  %s (Boltz-2): %d receptor residues within 5 A of the predicted toxin"
          % (tag, n))
    # camera fixed from the deposited ligand coordinates -> identical framing to
    # the deposited-column panel of this anchor
    cmd.orient("lig")
    cmd.zoom("lig", 3.2)
    cmd.disable("rec"); cmd.disable("lig")
    cmd.hide("everything", "rec or lig")
    render(tag)


PANELS = [
    dict(tag="S6_A1_9NX0_deposited", rec="9NX0_receptorPAIR_DE.pdb",
         lig="9NX0_ligand_I.pdb",
         pred=os.path.join(ROOT, "pod_mirror_boltz/outputs_ms/boltz_results_A1_9NX0/"
                                  "predictions/A1_9NX0/A1_9NX0_model_1.pdb"),
         ligch="I"),
    dict(tag="S6_A2_7EGR_deposited", rec="7EGR_receptorPAIR_HG.pdb",
         lig="7EGR_ligand_P.pdb",
         pred=os.path.join(ROOT, "pod_mirror_n2/outputs_b2_7EGR/boltz_results_A2_7EGR/"
                                  "predictions/A2_7EGR/A2_7EGR_model_4.pdb"),
         ligch="P"),
]

if __name__ in ("__main__", "pymol"):
    for p in PANELS:
        deposited_panel(p["tag"], p["rec"], p["lig"], C_TOX)
        prediction_panel(p["tag"].replace("_deposited", "_boltz2"),
                          p["rec"], p["lig"], p["pred"], p["ligch"], C_TOX2)
    print("\nrenders in", OUT)
