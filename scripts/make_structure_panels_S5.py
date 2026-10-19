#!/usr/bin/env python3
"""
Figure S5 structural panels -- swarm-centre geometry around the native binding site,
rendered headless with PyMOL 3.x.

Run:  /opt/homebrew/bin/pymol -cq scripts/make_structure_panels_S5.py

Conventions are inherited from scripts/make_structure_panels.py (same colour palette,
same base_setup, same white/ray-traced/no-shadow style) so this panel reads as part of
the same figure family as Figure 1 / Figure 4a. Nothing here invents a new colour code:
the rung colours (BLIND/XDON/NATIVE) are the exact RUNG_C triplet already used for the
rung comparisons in Figure 1c/2/3 (scripts/make_figures_v2.py), just carried into 3D.

*** COORDINATE-FRAME NOTE (read before touching this file) ***
pod_mirror/structures/lightdock_inputs/*.pdb are in the DEPOSITED frame. LightDock's own
working files (receptor.pdb copies inside each run directory, and every glowworm
coordinate written to swarm_N/gso_*.out) are in a RECEPTOR-CENTRED frame -- the receptor
CA centroid sits at the origin there, not at its deposited-frame position. This script
never mixes the two without an explicit translation: swarm-centre coordinates were
reconstructed by scripts/reconstruct_swarm_centres_S5.py as
approx = mean(200 generation-0 glowworm translations per swarm) + deposited-frame
receptor CA centroid, and written to
results/s0_ladder/swarm_centres_reconstructed_S5.json (already in the deposited frame --
this script does no further translation). That reconstruction was validated against
results/s0_ladder/swarm_coverage.json's own min_dist_centre_to_native_centroid_A and
n_centres_within_2radius_of_centroid for all six runs (agreement to 0.1-0.7 A on the
distance and within 1-2 counts on the shell population -- expected from averaging ~200
noisy samples per swarm, since raw per-swarm centre coordinates were never persisted:
setup.json has write_starting_positions:false for the production runs). The NUMBERS
drawn on the figure are the exact canonical values from swarm_coverage.json / Table S10;
only the marker POSITIONS are the (disclosed, validated) reconstruction.
"""
import json
import math
import os

from pymol import cmd
from pymol.cgo import ALPHA, BEGIN, COLOR, CYLINDER, END, LINE_LOOP, LINEWIDTH, SPHERE, VERTEX

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
LD = os.path.join(ROOT, "pod_mirror/structures/lightdock_inputs")
OUT = os.path.join(ROOT, "manuscript", "figures_v2", "panels")
os.makedirs(OUT, exist_ok=True)

_RECON = json.load(open(os.path.join(
    ROOT, "results/s0_ladder/swarm_centres_reconstructed_S5.json")))
# reshape scripts/reconstruct_swarm_centres_S5.py's flat {"A1_BLIND": {...}, ...} output
# into the {"A1": {"BLIND": {...}, "BLIND_cloud_xyz": [...]}, ...} shape this script uses
DATA = {}
for _tag, _d in _RECON.items():
    _anc, _rung = _tag.split("_", 1)
    DATA.setdefault(_anc, {})[_rung] = {
        "n_swarms": _d["n_swarms"],
        "nearest_point_xyz": _d["nearest_point_xyz_deposited"],
        "reconstructed_min_dist": _d["reconstructed_min_dist_A"],
        "canonical_min_dist": _d["canonical_min_dist_A_swarm_coverage_json"],
    }
    if _rung == "BLIND":
        DATA[_anc]["BLIND_cloud_xyz"] = _d["centers_xyz_deposited"]

W, H, DPI = 1700, 1700, 400

# -- exact palette reused from scripts/make_structure_panels.py --
C_PLUS = "0x6FA8DC"
C_MINUS = "0xB4B4B4"
C_TOX = "0xE69F00"

# -- exact rung palette reused from scripts/make_figures_v2.py RUNG_C --
RUNG_RGB = {
    "BLIND": (0x00 / 255, 0x72 / 255, 0xB2 / 255),
    "XDON": (0xE6 / 255, 0x9F / 255, 0x00 / 255),
    "NATIVE": (0xCC / 255, 0x79 / 255, 0xA7 / 255),
}
CLOUD_RGB = (0.55, 0.55, 0.55)


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


def sphere_cgo(xyz, r, rgb, alpha=1.0):
    obj = [COLOR, rgb[0], rgb[1], rgb[2]]
    if alpha < 1.0:
        obj = [ALPHA, alpha] + obj
    obj += [SPHERE, xyz[0], xyz[1], xyz[2], r]
    return obj


def cylinder_cgo(a, b, r, rgb):
    return [CYLINDER, a[0], a[1], a[2], b[0], b[1], b[2], r,
            rgb[0], rgb[1], rgb[2], rgb[0], rgb[1], rgb[2]]


def wire_circle_cgo(center, r, plane, rgb, n=72, width=2.0):
    """A closed line loop of radius r around `center`, lying in the given plane
    ('xy', 'xz', or 'yz'). Used as a scale ring instead of a solid CGO sphere --
    solid alpha spheres do not ray-trace as translucent reliably at this radius."""
    obj = [LINEWIDTH, width, COLOR, rgb[0], rgb[1], rgb[2], BEGIN, LINE_LOOP]
    for i in range(n):
        t = 2 * math.pi * i / n
        c, s = r * math.cos(t), r * math.sin(t)
        if plane == "xy":
            v = (center[0] + c, center[1] + s, center[2])
        elif plane == "xz":
            v = (center[0] + c, center[1], center[2] + s)
        else:
            v = (center[0], center[1] + c, center[2] + s)
        obj += [VERTEX, v[0], v[1], v[2]]
    obj += [END]
    return obj


def read_all_atom_centroid(path):
    xs, ys, zs, n = 0.0, 0.0, 0.0, 0
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                xs += float(line[30:38]); ys += float(line[38:46]); zs += float(line[46:54])
                n += 1
    return (xs / n, ys / n, zs / n)


def panel(tag, recfile, ligfile, anchor_key, closest_rung):
    base_setup()
    cmd.load(os.path.join(LD, recfile), "rec")
    cmd.load(os.path.join(LD, ligfile), "tox")
    cmd.hide("everything")
    cmd.show("cartoon", "rec")
    # receptor rendered as a single grey cartoon per the task spec (not the +/- two-tone
    # used in Fig 1/4 -- this panel's point is swarm geometry, not interface chemistry)
    cmd.color(C_MINUS, "rec")
    cmd.set("cartoon_transparency", 0.25, "rec")
    cmd.show("cartoon", "tox")
    cmd.show("sticks", "tox")
    cmd.color(C_TOX, "tox")

    lig_centroid = read_all_atom_centroid(os.path.join(LD, ligfile))

    cgo = []
    # background swarm cloud (BLIND rung only, real reconstructed positions, ~390/134 pts)
    for p in DATA[anchor_key]["BLIND_cloud_xyz"]:
        cgo += sphere_cgo(p, 0.9, CLOUD_RGB, alpha=0.55)

    # per-rung nearest swarm-centre marker + arrow from ligand centroid. (A same-
    # centre "white halo" sphere was tried here to separate the XDON marker
    # (rung colour = same orange as the native-ligand sticks, C_TOX) from the
    # ligand, but a same-centre larger sphere fully occludes the smaller one
    # from every viewing angle -- it does not read as an outline. Disambiguation
    # instead relies on the arrow (a cylinder back to the ligand centroid, absent
    # from the real ligand sticks) plus the colour-matched distance label placed
    # under each panel in the matplotlib composite.)
    for rung in ("BLIND", "XDON", "NATIVE"):
        d = DATA[anchor_key][rung]
        p = d["nearest_point_xyz"]
        rgb = RUNG_RGB[rung]
        cgo += cylinder_cgo(lig_centroid, p, 0.5, rgb)
        r = 2.7 if rung == closest_rung else 2.1
        cgo += sphere_cgo(p, r, rgb, alpha=1.0)

    # one 10 A radius wire-frame ring (3 orthogonal great circles) for scale,
    # centred on the overall-nearest swarm centre
    scale_p = DATA[anchor_key][closest_rung]["nearest_point_xyz"]
    for plane in ("xy", "xz", "yz"):
        cgo += wire_circle_cgo(scale_p, 10.0, plane, (0.35, 0.35, 0.35), width=1.6)

    cmd.load_cgo(cgo, "swarms")

    # "swarms" is a CGO object, not an atom selection -- it cannot be combined with
    # "rec" in a boolean selection expression, but "all" and a bare orient() do
    # include CGO extents, so zoom on those instead of "rec or swarms".
    cmd.orient()
    cmd.zoom("all", 8.0)
    render(tag)

    # detail inset: same view direction, zoomed on the ligand with enough buffer to
    # keep all three rung markers and the scale ring in frame -- the whole-receptor
    # scale (13-15 A at A1) makes the rung markers too small to tell apart above.
    cmd.zoom("tox", 26.0)
    render(tag + "_detail")


if __name__ in ("__main__", "pymol"):
    print("S5 panel -- A1 9NX0")
    panel("S5_A1_9NX0", "9NX0_receptorPAIR_DE.pdb", "9NX0_ligand_I.pdb",
          "A1", closest_rung="NATIVE")
    print("S5 panel -- A2 7EGR")
    panel("S5_A2_7EGR", "7EGR_receptorPAIR_HG.pdb", "7EGR_ligand_P.pdb",
          "A2", closest_rung="BLIND")
    print("\nrenders in", OUT)
