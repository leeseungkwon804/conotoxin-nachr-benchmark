#!/usr/bin/env python3
"""
Figure 5c structural panel, rendered headless with PyMOL 3.x.
Run:  /opt/homebrew/bin/pymol -cq scripts/make_structure_panels_v3.py

Shows, on one receptor, where each method's pose actually landed at anchor A1 (9NX0, human alpha7):
the deposited ImI, the best Boltz-2 sample, the best pose the docking search reached, and the pose
the scoring function returned. The last of these sits on the transmembrane domain, 88 A away, and
that is the panel's whole point — so the view has to include the whole receptor rather than a
close-up of the orthosteric site.

Conventions follow scripts/make_structure_panels.py (which produced the Figure 1 and 4a renders):
white background, ray-traced, shadows and depth cue off, the two receptor subunits in distinct
neutral greys so the (+)/(-) interface reads without colour coding, cartoon + sticks on the toxins.

Colours match the plotted panels: deposited = grey, Boltz-2 = green (#009E73),
docking best pose = blue (#0072B2), docking top-1 by score = vermillion (#D55E00).

All coordinates are in the DEPOSITED frame. The docked poses were written there by
working/poses_fig5c/ (LightDock working files are origin-centred and must never be mixed with
deposited coordinates); the Boltz model is superposed onto the deposited receptor here.
"""
import os

from pymol import cmd

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
LD = os.path.join(ROOT, "pod_mirror/structures/lightdock_inputs")
POSES = os.path.join(ROOT, "working/poses_fig5c")
BOLTZ = os.path.join(ROOT, "pod_mirror_boltz/outputs_ms/boltz_results_A1_9NX0/predictions/A1_9NX0")
OUT = os.path.join(ROOT, "manuscript", "figures_v3", "panels")
os.makedirs(OUT, exist_ok=True)

W, H, DPI = 2000, 1700, 400

def setup():
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("stick_radius", 0.16)
    cmd.set("valence", 0)
    cmd.set("depth_cue", 0)
    cmd.set("specular", 0.15)
    cmd.set("cartoon_transparency", 0.0)
    # reinitialize() clears user colours, so they are defined after it, not before
    cmd.set_color("c_dep", [0.30, 0.30, 0.30])   # deposited ligand: near-black, not the receptor grey
    cmd.set_color("c_boltz", [0.0, 0.62, 0.45])       # #009E73
    cmd.set_color("c_dockbest", [0.0, 0.447, 0.698])  # #0072B2
    cmd.set_color("c_docktop1", [0.835, 0.369, 0.0])  # #D55E00


setup()

# receptor: the docked subunit pair, in the deposited frame
cmd.load(os.path.join(LD, "9NX0_receptorPAIR_DE.pdb"), "rec")
cmd.hide("everything", "rec")
cmd.show("cartoon", "rec")
cmd.color("grey90", "rec and chain D")
cmd.color("grey80", "rec and chain E")
cmd.set("cartoon_transparency", 0.70, "rec")

# deposited ImI
cmd.load(os.path.join(LD, "9NX0_ligand_I.pdb"), "dep")
cmd.hide("everything", "dep")
cmd.show("sticks", "dep")
cmd.set("stick_radius", 0.26, "dep")
cmd.color("c_dep", "dep")

# best Boltz-2 sample, superposed on the receptor by its own receptor chains
cmd.load(os.path.join(BOLTZ, "A1_9NX0_model_0.pdb"), "b2")
# the toxin is the short chain; align on everything else
cmd.select("b2_rec", "b2 and polymer and not (byres (b2 and polymer and chain I))")
try:
    cmd.align("b2_rec", "rec")
except Exception:
    cmd.align("b2 and name CA", "rec and name CA")
cmd.hide("everything", "b2")
cmd.select("b2_lig", "b2 and chain I")
if cmd.count_atoms("b2_lig") == 0:
    # Boltz keeps input chain ids, but fall back to the shortest chain if it did not
    chains = cmd.get_chains("b2")
    shortest = min(chains, key=lambda c: cmd.count_atoms("b2 and chain %s" % c))
    cmd.select("b2_lig", "b2 and chain %s" % shortest)
cmd.show("sticks", "b2_lig")
cmd.set("stick_radius", 0.17, "b2_lig")
cmd.color("c_boltz", "b2_lig")

# docked poses, already in the deposited frame
cmd.load(os.path.join(POSES, "A1_9NX0_BLIND_rep2_best_rmsd.pdb"), "dock_best")
cmd.load(os.path.join(POSES, "A1_9NX0_BLIND_rep2_top1_score.pdb"), "dock_top1")
# docking poses drawn thicker than the co-folded one so they are not hidden behind it
for obj, col, rad in (("dock_best", "c_dockbest", 0.22), ("dock_top1", "c_docktop1", 0.30)):
    cmd.hide("everything", obj)
    cmd.show("sticks", obj)
    cmd.set("stick_radius", rad, obj)
    cmd.color(col, obj)

cmd.deselect()
cmd.orient("rec")
cmd.zoom("rec or dock_top1", 4)
cmd.turn("x", -12)

cmd.png(os.path.join(OUT, "Figure5_c_structure.png"), width=W, height=H, dpi=DPI, ray=1)
print("wrote Figure5_c_structure.png")

# a second view, rotated, so the 88 A separation is not an artefact of one projection
cmd.turn("y", 90)
cmd.png(os.path.join(OUT, "Figure5_c_structure_rot90.png"), width=W, height=H, dpi=DPI, ray=1)
print("wrote Figure5_c_structure_rot90.png")

# close-up of the orthosteric site: deposited vs Boltz-2 only, to show they coincide
cmd.zoom("dep", 6)
cmd.disable("dock_top1")
cmd.png(os.path.join(OUT, "Figure5_c_site_closeup.png"), width=W, height=H, dpi=DPI, ray=1)
print("wrote Figure5_c_site_closeup.png")
