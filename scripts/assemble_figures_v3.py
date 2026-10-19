#!/usr/bin/env python3
"""
Rough composites of Figures 5-7 (and Figure 2 with its new panel e) from the individual panels.

These are a starting layout to be rearranged in PowerPoint, not the final figures — the point is
that DK has something with the right panels in roughly the right place rather than a bag of PNGs.
Panel letters are drawn here so the composite is readable on its own; the PPTX version will carry
its own.

Usage: python3 assemble_figures_v3.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
V3 = os.path.join(ROOT, "manuscript", "figures_v3")
PAN = os.path.join(V3, "panels")
MM = 1 / 25.4

matplotlib.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42,
                            "ps.fonttype": 42, "font.size": 9})


def place(ax, name, letter, dx=-0.02, dy=1.0):
    f = os.path.join(PAN, name + ".png")
    if os.path.exists(f):
        ax.imshow(mpimg.imread(f))
    else:
        ax.text(0.5, 0.5, "%s\n(not yet made)" % name, ha="center", va="center",
                fontsize=8, color="#B00020", transform=ax.transAxes)
        ax.set_facecolor("#FAFAFA")
    ax.axis("off")
    if letter:
        ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="top", ha="left")


def build(name, layout, figsize):
    """layout: list of (gridspec-slice, panel-name, letter)."""
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(*layout["grid"], hspace=layout.get("hspace", 0.10),
                          wspace=layout.get("wspace", 0.06))
    for spec, pname, letter in layout["items"]:
        place(fig.add_subplot(gs[spec]), pname, letter)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(V3, "%s.%s" % (name, ext)), dpi=400,
                    bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("  wrote", name)


print("assembling into", V3)

# Figure 5 — contact map is tall, so it takes the left column full height
build("Figure5", {
    "grid": (2, 2), "wspace": 0.04, "hspace": 0.14,
    "items": [((slice(0, 2), 0), "Figure5_a", "a"),
              ((0, 1), "Figure5_b", "b"),
              ((1, 1), "Figure5_c", "c")],
}, (180 * MM, 150 * MM))

# Figure 6 — two wide panels stacked
build("Figure6", {
    "grid": (2, 1), "hspace": 0.18,
    "items": [((0, 0), "Figure6_a", "a"), ((1, 0), "Figure6_b", "b")],
}, (170 * MM, 130 * MM))

# Figure 7 — two wide panels stacked
build("Figure7", {
    "grid": (2, 1), "hspace": 0.20,
    "items": [((0, 0), "Figure7_a", "a"), ((1, 0), "Figure7_b", "b")],
}, (170 * MM, 130 * MM))

# Figure 2 — existing a-d plus the new e
build("Figure2_with_e", {
    "grid": (3, 2), "hspace": 0.16, "wspace": 0.06,
    "items": [((0, 0), "Figure2_a", "a"), ((0, 1), "Figure2_b", "b"),
              ((1, 0), "Figure2_c", "c"), ((1, 1), "Figure2_d", "d"),
              ((2, slice(0, 2)), "Figure2_e", "e")],
}, (180 * MM, 200 * MM))
print("done")
