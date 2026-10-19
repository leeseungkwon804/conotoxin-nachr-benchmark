#!/usr/bin/env /usr/bin/python3
"""
Figure S6 -- residue-level view of the co-folded interface at both anchors.

Run with /usr/bin/python3 -- see scripts/make_figures_v2.py header for why
(Homebrew python3's matplotlib is broken here).

Composites the four PyMOL renders from scripts/make_figureS6_panels.py
(S6_A1_9NX0_deposited / _boltz2, S6_A2_7EGR_deposited / _boltz2) into a 2x2
grid: rows = anchor (A1 9NX0 top, A2 7EGR bottom), columns = (deposited
structure, Boltz-2 prediction). Every panel carries an explicit leakage-risk
stamp (project rule: 7EGR and 9NX0 must never be shown side by side without
one) -- HIGH / positive control for 7EGR, MEDIUM for 9NX0 -- because 7EGR
(released 2022-04-06) sits inside the Boltz-2 training window and 9NX0
(released 2026-02-04) does not, but the co-folding pipeline treats both the
same way, so the panel itself has to carry the warning.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.8,
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8,
})

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
# figures_v2 is frozen; set FIGOUT to render a revised version elsewhere without touching it
OUT = os.environ.get("FIGOUT") or os.path.join(ROOT, "manuscript", "figures_v2")
PAN = os.path.join(OUT, "panels")
MM = 1 / 25.4

C_HIGH = "#B22222"   # HIGH leakage stamp -- firebrick
C_MED = "#B8860B"    # MEDIUM leakage stamp -- dark goldenrod
C_TOX = "#E69F00"
C_TOX2 = "#CC79A7"
C_PLUS = "#6FA8DC"
C_MINUS = "#B4B4B4"

# receptor CA superposition RMSD after exhaustive homomeric chain-pair
# enumeration, printed by make_figureS6_panels.py -- read here only to keep
# the caption's own numbers traceable to that run's stdout, not re-derived.
SUPERPOSITION_RMSD = {"A1": 7.235, "A2": 0.832}
N_IFACE_RES = {("A1", "deposited"): 28, ("A1", "boltz2"): 25,
               ("A2", "deposited"): 33, ("A2", "boltz2"): 32}

ROWS = [
    dict(key="A1", pdb="9NX0", tox="ImI", iface_label="$\\alpha$7(+)$\\alpha$7(−) interface",
         leak="MEDIUM", leak_full="leakage: MEDIUM", stamp_c=C_MED,
         released="released 2026-02-04"),
    dict(key="A2", pdb="7EGR", tox="RgIA", iface_label="AChBP(+)(−) interface",
         leak="HIGH", leak_full="leakage: HIGH — positive control", stamp_c=C_HIGH,
         released="released 2022-04-06"),
]
COLS = [
    dict(key="deposited", title="Deposited structure"),
    dict(key="boltz2", title="Boltz-2 prediction (median model)"),
]


def stamp(ax, text, color):
    """Colored leakage-risk badge, top-left corner, inside the panel."""
    ax.text(0.02, 0.985, text, transform=ax.transAxes, ha="left", va="top",
             fontsize=6.6, fontweight="bold", color="white",
             bbox=dict(boxstyle="round,pad=0.28", facecolor=color,
                       edgecolor="none", alpha=0.92), zorder=10)


def panel_tag(ax, s):
    ax.text(0.985, 0.985, s, transform=ax.transAxes, ha="right", va="top",
             fontsize=8.5, fontweight="bold", color="black",
             bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                       edgecolor="black", linewidth=0.6, alpha=0.85), zorder=10)


def main():
    fig = plt.figure(figsize=(170 * MM, 190 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.10, wspace=0.04,
                          left=0.045, right=0.99, top=0.93, bottom=0.20)

    letters = iter("ABCD")
    for i, row in enumerate(ROWS):
        for j, col in enumerate(COLS):
            ax = fig.add_subplot(gs[i, j])
            name = "S6_%s_%s_%s" % (row["key"], row["pdb"], col["key"])
            p = os.path.join(PAN, name + ".png")
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
            if os.path.exists(p):
                ax.imshow(mpimg.imread(p))
            else:
                ax.text(0.5, 0.5, "[%s]\nrender missing" % name, transform=ax.transAxes,
                        ha="center", va="center", fontsize=7, color="#999999")
            stamp(ax, row["leak_full"], row["stamp_c"])
            panel_tag(ax, next(letters))
            n = N_IFACE_RES.get((row["key"], col["key"]))
            sub = "%d receptor residues within 5 Å of the toxin" % n if n else ""
            if i == 0:
                ax.set_title(col["title"], fontsize=8.5, pad=6)
            if j == 0:
                ax.text(-0.045, 0.5, "%s   %s\n%s" % (row["key"], row["pdb"], row["iface_label"]),
                        transform=ax.transAxes, rotation=90, ha="center", va="center",
                        fontsize=8.2, fontweight="bold")
            if sub:
                ax.text(0.5, -0.025, sub, transform=ax.transAxes, ha="center", va="top",
                        fontsize=6.3, color="#333333")

    legend_y = 0.155
    swatches = [(C_PLUS, "receptor (+) principal subunit"),
                (C_MINUS, "receptor (−) complementary subunit"),
                (C_TOX, "deposited toxin"),
                (C_TOX2, "Boltz-2 predicted toxin"),
                ("#CC0000", "polar contact ≤ 3.6 Å (dashed)")]
    x0 = 0.06
    for c, lab in swatches:
        if lab.startswith("polar"):
            fig.add_artist(plt.Line2D([x0, x0 + 0.018], [legend_y, legend_y],
                                      transform=fig.transFigure, color=c, lw=1.3,
                                      linestyle=(0, (2, 1.3))))
        else:
            fig.patches.append(plt.Rectangle((x0, legend_y - 0.006), 0.018, 0.012,
                                             transform=fig.transFigure, facecolor=c,
                                             edgecolor="none"))
        fig.text(x0 + 0.024, legend_y, lab, fontsize=6.4, va="center", ha="left")
        x0 += 0.024 + 0.012 * len(lab) + 0.012

    caption = (
        "Figure S6. Co-folded interface at both anchors (A1 = 9NX0/ImI/$\\alpha$7 nAChR; "
        "A2 = 7EGR/RgIA/AChBP; see Supplementary Materials Guide §0). Rows: anchor (top A1, "
        "bottom A2); columns: deposited structure vs. Boltz-2 prediction, shown as the "
        "median-RMSD (not top-confidence) model: A1 1/5, A2 4/5. A2's five poses are "
        "index-permuted duplicates, not an independent n = 10 (pod_mirror_boltz/eval_ms, "
        "pod_mirror_n2/eval_7EGR_b1_vs_b2.json). Interface residues and receptor "
        "C$\\alpha$ superposition follow Methods §4.4 (RMSD 7.235 Å at A1, 0.832 Å at A2); "
        "dashed lines mark polar contacts $\\leq$3.6 Å. Leakage-risk stamps per Guide §0; "
        "timeline in Figure 4d."
    )
    # The caption is drawn into the figure for standalone review copies. Set FIGNOCAP=1 for
    # the submission render, where the caption lives in the document and a second copy
    # baked into the image would be a duplicate.
    if not os.environ.get("FIGNOCAP"):
        fig.text(0.045, 0.125, caption, fontsize=6.1, va="top", ha="left", wrap=True,
                 linespacing=1.45, transform=fig.transFigure)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS6." + e)
        fig.savefig(p, dpi=int(os.environ.get("FIGDPI", 400)), bbox_inches="tight")
    plt.close(fig)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS6." + e)
        print("wrote %-40s %8.1f KB" % (p, os.path.getsize(p) / 1024))


if __name__ == "__main__":
    main()
