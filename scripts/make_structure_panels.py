#!/usr/bin/env python3
"""
Structural panels, rendered headless with PyMOL 3.x.

Run:  /opt/homebrew/bin/pymol -cq scripts/make_structure_panels.py

Produces the raw renders that Figure 1 and Figure 4a are assembled from. Every
selection here is taken from the deposited coordinates in
pod_mirror/structures/lightdock_inputs/, which are in the DEPOSITED frame -- the
LightDock working files are origin-centred and must never be mixed with these.

Conventions followed (see figs_research/00_피겨관행_조사.md):
  * receptor as cartoon, toxin as cartoon + sticks so the disulfide scaffold reads;
  * the two receptor subunits get distinct neutral colours so the (+)/(-) interface
    is visible without relying on red/green;
  * white background, ray-traced, ambient occlusion off for panel clarity;
  * a consistent view matrix per anchor so panels can be compared.
"""
import os

from pymol import cmd

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
LD = os.path.join(ROOT, "pod_mirror/structures/lightdock_inputs")
OUT = os.path.join(ROOT, "manuscript", "figures_v2", "panels")
os.makedirs(OUT, exist_ok=True)

W, H, DPI = 1700, 1500, 400

# colour-vision-safe, no red/green pairing
C_PLUS = "0x6FA8DC"    # principal (+) subunit, light blue
C_MINUS = "0xB4B4B4"   # complementary (-) subunit, grey
C_TOX = "0xE69F00"     # toxin, orange
C_TOX2 = "0xCC79A7"    # second toxin / predicted, purple
C_APO = "0x009E73"     # apo conformer, teal


def base_setup():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("stick_radius", 0.16)
    cmd.set("valence", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.15)


def render(name):
    cmd.png(os.path.join(OUT, name + ".png"), width=W, height=H, dpi=DPI, ray=1)
    print("  rendered", name)


# --------------------------------------------------------------- 1a  whole assemblies
def panel_1a():
    for tag, rec, lig, ligname in (
            ("1a_A1_9NX0", "9NX0_receptor_ABCDE.pdb", "9NX0_ligand_I.pdb", "ImI"),
            ("1a_A2_7EGR", "7EGR_receptor_FGHIJ.pdb", "7EGR_ligand_P.pdb", "RgIA")):
        base_setup()
        cmd.load(os.path.join(LD, rec), "rec")
        cmd.load(os.path.join(LD, lig), "tox")
        cmd.hide("everything")
        cmd.show("cartoon", "rec")
        cmd.color(C_MINUS, "rec")
        # highlight the two subunits that form the modelled site
        pair = "chain D or chain E" if tag.endswith("9NX0") else "chain H or chain G"
        cmd.color(C_PLUS, "rec and (%s)" % pair)
        cmd.show("cartoon", "tox")
        cmd.show("sticks", "tox")
        cmd.color(C_TOX, "tox")
        cmd.set("cartoon_transparency", 0.55, "rec")
        cmd.set("cartoon_transparency", 0.0, "rec and (%s)" % pair)
        # orient on the ECD ring only, so the TMD (present in 9NX0) does not dominate
        # the principal axes and tip the pentamer on its side.
        # The pseudo-5-fold axis is the principal axis of the whole pentamer. orient()
        # puts the LONGEST axis horizontal, so after orient() the 5-fold runs left-right
        # and a -90 turn about y brings us to look straight down it.
        cmd.orient("rec")
        cmd.zoom("rec", 1.0)
        render(tag + "_side")
        cmd.turn("y", -90)
        cmd.zoom("rec", 1.0)
        render(tag + "_top")


# --------------------------------------------------------------- 1b  interface close-up
def panel_1b():
    for tag, rec, lig in (("1b_A1_9NX0", "9NX0_receptorPAIR_DE.pdb", "9NX0_ligand_I.pdb"),
                          ("1b_A2_7EGR", "7EGR_receptorPAIR_HG.pdb", "7EGR_ligand_P.pdb")):
        base_setup()
        cmd.load(os.path.join(LD, rec), "rec")
        cmd.load(os.path.join(LD, lig), "tox")
        chains = cmd.get_chains("rec")
        cmd.hide("everything")
        cmd.show("cartoon", "rec")
        cmd.color(C_PLUS, "rec and chain %s" % chains[0])
        cmd.color(C_MINUS, "rec and chain %s" % chains[1])
        cmd.show("cartoon", "tox"); cmd.show("sticks", "tox"); cmd.color(C_TOX, "tox")
        # native contact residues within 5 A of the toxin
        cmd.select("iface", "byres (rec within 5.0 of tox)")
        cmd.show("sticks", "iface and not (name C+N+O)")
        cmd.color("0x444444", "iface and elem C")
        n = cmd.count_atoms("iface and name CA")
        print("  %s: %d receptor residues within 5 A of the toxin" % (tag, n))
        cmd.set("cartoon_transparency", 0.35, "rec")
        cmd.orient("tox")
        cmd.zoom("tox", 3.0)
        render(tag)


# --------------------------------------------------------------- 1d  apo vs bound
def panel_1d():
    """The rigid-body limit: the apo conformer that was docked, superposed on the
    bound conformer. The residual after optimal superposition is the floor."""
    import glob
    for tag, native, apo_glob, floor in (
            ("1d_A1_ImI", "9NX0_ligand_I.pdb",
             os.path.join(ROOT, "pod_mirror/apo_ligands/ImI/ImI_1IMI_m01.pdb"), 2.439),
            ("1d_A2_RgIA", "7EGR_ligand_P.pdb",
             os.path.join(ROOT, "pod_mirror/apo_ligands/RgIA/RgIA_2JUT_m01.pdb"), 4.453)):
        if not os.path.exists(apo_glob):
            cand = sorted(glob.glob(os.path.join(os.path.dirname(apo_glob), "*_m01.pdb")))
            if not cand:
                print("  SKIP %s (no apo conformer found at %s)" % (tag, apo_glob))
                continue
            apo_glob = cand[0]
        base_setup()
        cmd.load(os.path.join(LD, native), "bound")
        cmd.load(apo_glob, "apo")
        cmd.hide("everything")
        for o, c in (("bound", C_TOX), ("apo", C_APO)):
            cmd.show("cartoon", o); cmd.show("sticks", o); cmd.color(c, o)
        rms = cmd.align("apo", "bound")[0]
        print("  %s: PyMOL align RMSD %.3f A (reported rigid-body floor %.3f A)" % (tag, rms, floor))
        cmd.orient("bound")
        cmd.zoom("bound", 1.2)
        render(tag)


# --------------------------------------------------------------- 4a  prediction overlay
def panel_4a():
    preds = [("4a_A1_9NX0", "9NX0_receptorPAIR_DE.pdb", "9NX0_ligand_I.pdb",
              os.path.join(ROOT, "pod_mirror_boltz/outputs_ms/boltz_results_A1_9NX0/"
                                 "predictions/A1_9NX0/A1_9NX0_model_4.pdb"), "I"),
             ("4a_A2_7EGR", "7EGR_receptorPAIR_HG.pdb", "7EGR_ligand_P.pdb",
              os.path.join(ROOT, "pod_mirror_boltz/outputs_ms/boltz_results_A2_7EGR/"
                                 "predictions/A2_7EGR/A2_7EGR_model_1.pdb"), "P")]
    for tag, rec, lig, pred, ligch in preds:
        if not os.path.exists(pred):
            print("  SKIP %s (prediction not found)" % tag)
            continue
        base_setup()
        cmd.load(os.path.join(LD, rec), "rec")
        cmd.load(os.path.join(LD, lig), "native")
        cmd.load(pred, "pred")
        # superpose the predicted receptor onto the native receptor, then carry the
        # predicted toxin along -- exactly what the RMSD evaluation does.
        cmd.select("pred_rec", "pred and not chain %s" % ligch)
        cmd.select("pred_lig", "pred and chain %s" % ligch)
        import itertools
        pred_ch = [c for c in cmd.get_chains("pred") if c != ligch]
        nat_ch = cmd.get_chains("rec")
        best = None
        for perm in itertools.permutations(pred_ch):
            cmd.create("probe", "pred")
            pairs = list(zip(perm, nat_ch))
            # pair-wise fit on the first mapped chain, then score over all mapped chains
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
        print("  %s: receptor CA superposition RMSD %.3f A  (chain map %s)"
              % (tag, rms, {p_: n_ for p_, n_ in pairs}))
        cmd.hide("everything")
        cmd.show("cartoon", "rec")
        chains = cmd.get_chains("rec")
        cmd.color(C_PLUS, "rec and chain %s" % chains[0])
        cmd.color(C_MINUS, "rec and chain %s" % chains[1])
        cmd.set("cartoon_transparency", 0.55, "rec")
        for o, c in (("native", C_TOX), ("pred_lig", C_TOX2)):
            cmd.show("cartoon", o); cmd.show("sticks", o); cmd.color(c, o)
        cmd.orient("native")
        cmd.zoom("native", 3.0)
        render(tag)


if __name__ in ("__main__", "pymol"):
    print("panel 1a -- whole assemblies")
    panel_1a()
    print("panel 1b -- interface close-up")
    panel_1b()
    print("panel 1d -- apo vs bound conformer")
    panel_1d()
    print("panel 4a -- prediction vs deposited")
    panel_4a()
    print("\nrenders in", OUT)
