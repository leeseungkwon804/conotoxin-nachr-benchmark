#!/usr/bin/env /usr/bin/python3
"""
Figure S5 -- swarm-centre geometry around the native binding site.

Run with /usr/bin/python3 -- see scripts/make_figures_v2.py header for why
(Homebrew python3's matplotlib is broken here).

Composites the two PyMOL renders from scripts/make_structure_panels_S5.py
(S5_A1_9NX0, S5_A2_7EGR) into a 2-panel figure (one panel per anchor), with a
grouped bar strip below showing escape fraction by rung (both anchors kept in
separate bars per project rule -- 7EGR/9NX0 leakage status differs and is never
pooled, though escape fraction itself is a docking-geometry quantity, not a
leakage-relevant one; the labelling is kept anyway for consistency with every
other cross-anchor panel in this project).

Source data (verify-on-open, not re-derived numbers):
  - results/s0_ladder/swarm_coverage.json           (min_dist_centre_to_native_centroid_A,
                                                       n_swarm_centres, n_centres_within_2radius,
                                                       n_centres_touching_native_ligand)
  - results/s0_ladder/ladder.json  rows[]            (escape_fraction, per run/rep)
  - results/s0_ladder/METHOD.md                      (corrected escape-fraction narrative,
                                                       cross-checked against ladder.json directly)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
    "legend.fontsize": 6.3, "xtick.labelsize": 7, "ytick.labelsize": 7,
})

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
# figures_v2 is frozen; set FIGOUT to render a revised version elsewhere without touching it
OUT = os.environ.get("FIGOUT") or os.path.join(ROOT, "manuscript", "figures_v2")
PAN = os.path.join(OUT, "panels")
MM = 1 / 25.4

C_BLIND, C_XDON, C_NATIVE = "#0072B2", "#E69F00", "#CC79A7"
C_A1, C_A2 = "#0072B2", "#D55E00"
C_GREY = "#5A5A5A"
RUNG_C = {"BLIND": C_BLIND, "XDON": C_XDON, "NATIVE": C_NATIVE}
RUNGS = ["BLIND", "XDON", "NATIVE"]

# canonical Table S10 numbers (results/s0_ladder/swarm_coverage.json,
# min_dist_centre_to_native_centroid_A -- re-verified against the raw JSON on open,
# and cross-checked in results/s0_ladder/METHOD.md's own table, 2026-08 pass)
NEAREST_A = {
    "A1": {"BLIND": 13.879, "XDON": 14.615, "NATIVE": 13.625},
    "A2": {"BLIND": 6.076, "XDON": 10.673, "NATIVE": 8.813},
}
# escape_fraction, both reps, from results/s0_ladder/ladder.json rows[]
ESCAPE = {
    "A1": {"BLIND": [0.98177, 0.98185], "XDON": [0.97047, 0.96988], "NATIVE": [0.95495, 0.95311]},
    "A2": {"BLIND": [0.97183, 0.96918], "XDON": [0.82153, 0.81877], "NATIVE": [0.89549, 0.89583]},
}
ANCHOR_LABEL = {"A1": "A1  9NX0 (human $\\alpha$7 nAChR + ImI)",
                "A2": "A2  7EGR ($\\it{Aplysia}$ AChBP + RgIA)"}


def panel_label(ax, s, dx=-0.06, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")


def struct_panel(ax, anchor, letter):
    name = "S5_%s_%s" % (anchor, "9NX0" if anchor == "A1" else "7EGR")
    p = os.path.join(PAN, name + ".png")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if os.path.exists(p):
        ax.imshow(mpimg.imread(p))
    else:
        ax.text(0.5, 0.5, "[%s]\nrender missing" % name, transform=ax.transAxes,
                ha="center", va="center", fontsize=7, color="#999999")
    ax.set_title(ANCHOR_LABEL[anchor], fontsize=8, pad=4)
    panel_label(ax, letter)

    # detail inset: same camera direction, zoomed on the ligand -- at whole-receptor
    # scale the three rung markers (13-15 A apart at A1) are too close together to
    # tell apart, so a cropped close-up is composited into the corner.
    dp = os.path.join(PAN, name + "_detail.png")
    if os.path.exists(dp):
        iax = ax.inset_axes([0.66, 0.02, 0.34, 0.34])
        iax.imshow(mpimg.imread(dp))
        iax.set_xticks([]); iax.set_yticks([])
        for s in iax.spines.values():
            s.set_visible(True); s.set_color("#666666"); s.set_linewidth(0.7)
        iax.set_title("detail", fontsize=5.6, pad=1.5, color="#444444")

    # per-rung nearest-swarm-centre distance, colour-matched to the arrows in the render
    y0 = -0.035
    xs = [0.16, 0.5, 0.84]
    for x, r in zip(xs, RUNGS):
        ax.text(x, y0, "%s\n%.2f Å" % (r, NEAREST_A[anchor][r]), transform=ax.transAxes,
                ha="center", va="top", fontsize=6.6, fontweight="bold", color=RUNG_C[r])
    ax.text(0.5, y0 - 0.115, "nearest swarm centre → native ligand centroid, by rung",
            transform=ax.transAxes, ha="center", va="top", fontsize=6, color="#333333")


def escape_panel(ax):
    x = np.arange(3); w = 0.36
    for off, anc, c in ((-w / 2, "A1", C_A1), (w / 2, "A2", C_A2)):
        means = [np.mean(ESCAPE[anc][r]) for r in RUNGS]
        los = [means[i] - min(ESCAPE[anc][r]) for i, r in enumerate(RUNGS)]
        his = [max(ESCAPE[anc][r]) - means[i] for i, r in enumerate(RUNGS)]
        ax.bar(x + off, means, w, color=c, alpha=0.9,
               label="A1 (9NX0)" if anc == "A1" else "A2 (7EGR)", zorder=2)
        ax.errorbar(x + off, means, yerr=[los, his], fmt="none", ecolor="black",
                     elinewidth=0.8, capsize=2, zorder=3)
        for xi, m in zip(x + off, means):
            ax.text(xi, m + 0.012, "%.3f" % m, ha="center", fontsize=5.6)
    ax.set_xticks(x)
    ax.set_xticklabels(["BLIND\nno restraints", "XDON\ncross-donated\nrestraints",
                        "NATIVE\ntarget's own\nrestraints"])
    ax.set_ylabel("Escape fraction\n(final centroid farther from its own\nswarm centre than 10 Å)")
    ax.set_ylim(0.70, 1.06)
    ax.axhspan(0.819, 0.982, color=C_GREY, alpha=0.08, zorder=0)
    ax.legend(frameon=False, loc="lower left", handletextpad=0.5, ncol=2)
    ax.set_title("Most glowworms leave their initial 10 Å swarm sphere in every run "
                 "(range 0.819–0.982 across all 12 runs, shaded)", fontsize=7.5, pad=5)


def main():
    fig = plt.figure(figsize=(170 * MM, 205 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 0.75],
                          hspace=0.55, wspace=0.10,
                          left=0.075, right=0.98, top=0.95, bottom=0.315)

    ax = fig.add_subplot(gs[0, 0]); struct_panel(ax, "A1", "A")
    ax = fig.add_subplot(gs[0, 1]); struct_panel(ax, "A2", "B")
    ax = fig.add_subplot(gs[1, :]); escape_panel(ax); panel_label(ax, "C", dx=-0.045)

    legend_y = 0.255
    swatches = [(C_GREY, "swarm centres, BLIND rung (n=390 / 134, illustrative positions -- see caption)"),
                (C_BLIND, "nearest swarm centre, BLIND"),
                (C_XDON, "nearest swarm centre, XDON"),
                (C_NATIVE, "nearest swarm centre, NATIVE"),
                ("#E69F00", "native ligand (deposited)")]
    x0 = 0.045
    for c, lab in swatches:
        ax0 = fig.add_axes([x0, legend_y - 0.008, 0.014, 0.014])
        ax0.set_xticks([]); ax0.set_yticks([])
        for s in ax0.spines.values():
            s.set_visible(False)
        ax0.scatter([0.5], [0.5], s=42, color=c)
        fig.text(x0 + 0.018, legend_y, lab, fontsize=6.2, va="center", ha="left")
        x0 += 0.018 + 0.0092 * len(lab) + 0.014

    caption = (
        "Figure S5. Swarm-centre geometry around the native binding site. (a,b) Receptor and native "
        "ligand (deposited frame), A1 (9NX0), A2 (7EGR); grey spheres = BLIND rung's swarm-centre "
        "cloud (n = 390 A1, 134 A2; positions reconstructed, illustrative -- exact distances in "
        "Table S10); rung-coloured spheres = each rung's nearest swarm centre to native centroid "
        "(Table S10); ring = 10 Å swarm radius (Methods §4.2). Inset: detail. (c) Escape fraction "
        "by rung, both anchors (Table S10); bars = mean, error bars = min/max, 2 replicate runs; "
        "shading = 0.819-0.982 range, 12 runs."
    )
    # The caption is drawn into the figure for standalone review copies. Set FIGNOCAP=1 for
    # the submission render, where the caption lives in the document and a second copy
    # baked into the image would be a duplicate.
    if not os.environ.get("FIGNOCAP"):
        fig.text(0.045, 0.225, caption, fontsize=5.85, va="top", ha="left", wrap=True,
                 linespacing=1.42, transform=fig.transFigure)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS5." + e)
        fig.savefig(p, dpi=int(os.environ.get("FIGDPI", 400)), bbox_inches="tight")
    plt.close(fig)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS5." + e)
        print("wrote %-40s %8.1f KB" % (p, os.path.getsize(p) / 1024))


if __name__ == "__main__":
    main()
