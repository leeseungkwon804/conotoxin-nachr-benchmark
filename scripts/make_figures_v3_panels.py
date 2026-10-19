#!/usr/bin/env python3
"""
Panels for the new Figures 5-7, exported individually so they can be laid out in PowerPoint.

Conventions are inherited from scripts/make_figures_v2.py deliberately — same rcParams, same
Okabe-Ito assignments, same 400 dpi png+pdf pair, same panels/ naming — so a v3 panel dropped
beside a v2 panel in the same deck does not read as a different figure set.

Colour assignments (validated: the three-rung set and the two-anchor set each pass the
lightness-band, chroma, CVD-separation, normal-vision and contrast checks):
    BLIND #0072B2   XDON #E69F00   NATIVE #CC79A7      A1 #0072B2   A2 #D55E00
Orange and pink sit below 3:1 against white, so every mark carrying them is directly labelled
rather than relying on the legend alone.

Numbers are hard-coded from the measured outputs rather than recomputed here; each block cites
the file it came from so a value can be traced back. Nothing in this script derives a result.

Usage: python3 make_figures_v3_panels.py
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 7.5,
    "legend.fontsize": 6.3, "xtick.labelsize": 7, "ytick.labelsize": 7,
})

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
PAN = os.path.join(ROOT, "manuscript", "figures_v3", "panels")
os.makedirs(PAN, exist_ok=True)

C_BLIND, C_XDON, C_NATIVE = "#0072B2", "#E69F00", "#CC79A7"
C_GREY, C_FLOOR = "#5A5A5A", "#009E73"
C_A1, C_A2 = "#0072B2", "#D55E00"
RUNG_C = {"BLIND": C_BLIND, "XDON": C_XDON, "NATIVE": C_NATIVE}
MM = 1 / 25.4

# Measured values, exported from the analysis JSONs on the pods (see panel_data.json's
# "provenance" block). Hard-coding per-sample numbers into a plotting script is how a figure
# ends up showing something the data never said, so nothing below invents a value.
import json
DATA = json.load(open(os.path.join(PAN, "panel_data.json")))


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(PAN, "%s.%s" % (name, ext)), dpi=400,
                    bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  wrote", name)


# ----------------------------------------------------------------- Figure 5b
# Interaction-type recovery. Source: round2 contact track.
# A1 LightDock BLIND best-RMSD pose recovers none of any type; that zero is the panel.
def fig5b():
    # Counts are 0-42 with one category at 1, so drawing the native total as a background bar
    # piles three labels on top of each other in that column. The total goes in the tick label
    # instead, and the bars carry only what each method reproduced.
    cats = ["Hydrogen bonds\n(13 in native)", "Ionic\n(1 in native)", "Hydrophobic\n(42 in native)"]
    boltz_lo, boltz_hi = [4, 0, 25], [7, 1, 28]

    fig, ax = plt.subplots(figsize=(70 * MM, 46 * MM))
    x = np.arange(len(cats))
    w = 0.34

    ax.bar(x - w / 2, boltz_lo, width=w, color=C_FLOOR, edgecolor="white", linewidth=0.8, zorder=3)
    for xi, lo, hi in zip(x - w / 2, boltz_lo, boltz_hi):
        if hi > lo:
            ax.plot([xi, xi], [lo, hi], color=C_FLOOR, lw=1.6, solid_capstyle="butt", zorder=4)
        ax.text(xi, hi + 0.9, "%d–%d" % (lo, hi) if hi > lo else "%d" % lo,
                ha="center", va="bottom", fontsize=6.2, color=C_FLOOR, zorder=5)

    ax.bar(x + w / 2, [0, 0, 0], width=w, color=C_BLIND, edgecolor="white", linewidth=0.8, zorder=3)
    for xi in x + w / 2:
        ax.text(xi, 0.9, "0", ha="center", va="bottom", fontsize=8,
                color=C_BLIND, fontweight="bold", zorder=5)

    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=6.3, linespacing=1.2)
    ax.set_ylabel("Native interactions reproduced")
    ax.set_ylim(0, 34)
    ax.set_title("A1 (9NX0)  human $\\alpha$7 – ImI", pad=6)
    ax.legend(handles=[Line2D([], [], color=C_FLOOR, lw=5, label="Boltz-2 (n = 5, range)"),
                       Line2D([], [], color=C_BLIND, lw=5, label="LightDock BLIND, best pose")],
              loc="upper center", bbox_to_anchor=(0.5, 1.0), frameon=False, handlelength=1.1)
    save(fig, "Figure5_b")


# ----------------------------------------------------------------- Figure 6a
# DockQ by method, anchors separate. Source: round2/dockq/, round3/ceiling/.
def fig6a():
    fig, axes = plt.subplots(1, 2, figsize=(120 * MM, 48 * MM), sharey=True)
    bands = [(0.80, 1.00, "High"), (0.49, 0.80, "Medium"), (0.23, 0.49, "Acceptable"),
             (0.00, 0.23, "Incorrect")]
    data = {
        "A1 (9NX0) human $\\alpha$7": [
            ("docking\nbest", 0.249, 0.387, C_BLIND),
            ("docking\ntop-1", 0.004, 0.237, C_GREY),
            ("Boltz-2", 0.826, 0.872, C_FLOOR),
        ],
        "A2 (7EGR) AChBP": [
            ("docking\nbest", 0.300, 0.604, C_BLIND),
            ("docking\ntop-1", 0.004, 0.250, C_GREY),
            ("Boltz-2", 0.909, 0.939, C_FLOOR),
            ("Boltz-1\ntemporal ctrl", 0.808, 0.821, "#56B4E9"),
        ],
    }
    for ax, (title, rows) in zip(axes, data.items()):
        for lo, hi, lab in bands:
            ax.axhspan(lo, hi, color="#F4F4F4" if lab in ("Medium", "Incorrect") else "#FBFBFB",
                       zorder=0, lw=0)
            # Band labels are drawn once, on the left panel only: the y-axis is shared, so the
            # bands are identical in both, and right-aligning them in each panel put "High" on
            # top of the A2 Boltz-1 value label (overlap measured at 3.5 x 5.0 px). Left edge is
            # clear of every bar because xlim starts at -0.6 while the first bar sits at 0.
            if ax is axes[0]:
                ax.text(0.015, (lo + hi) / 2, lab, transform=ax.get_yaxis_transform(),
                        fontsize=5.8, color="#9A9A9A", va="center", ha="left", zorder=1)
        for i, (lab, lo, hi, c) in enumerate(rows):
            ax.plot([i, i], [lo, hi], color=c, lw=5, solid_capstyle="round", zorder=3)
            ax.text(i, hi + 0.03, "%.2f" % hi, ha="center", va="bottom", fontsize=6.2, color=c)
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([r[0] for r in rows], fontsize=6.3, linespacing=1.15)
        ax.set_xlim(-0.6, len(rows) - 0.4)
        ax.set_title(title, pad=6)
    axes[0].set_ylabel("DockQ")
    axes[0].set_ylim(0, 1.02)
    save(fig, "Figure6_a")


# ----------------------------------------------------------------- Figure 6b
# fnat: achieved vs rigid-placement ceiling vs weakest co-folded sample.
# Source: round3/ceiling/search_*.json, round2/dockq/.
def fig6b():
    fig, ax = plt.subplots(figsize=(78 * MM, 48 * MM))
    groups = ["A1 (9NX0)", "A2 (7EGR)"]
    achieved = [13 / 53, 31 / 56]
    ceiling = [47 / 53, 43 / 56]
    boltz_weak = [np.nan, 47 / 56]

    x = np.arange(2)
    w = 0.26
    ax.bar(x - w, achieved, width=w, color=C_BLIND, edgecolor="white", lw=0.8, label="docking, best of all poses")
    ax.bar(x, ceiling, width=w, color=C_GREY, edgecolor="white", lw=0.8, label="rigid-placement ceiling")
    ax.bar(x + w, boltz_weak, width=w, color=C_FLOOR, edgecolor="white", lw=0.8, label="Boltz-2, weakest sample")

    for xi, v, n, tot in [(x[0] - w, achieved[0], 13, 53), (x[1] - w, achieved[1], 31, 56),
                          (x[0], ceiling[0], 47, 53), (x[1], ceiling[1], 43, 56),
                          (x[1] + w, boltz_weak[1], 47, 56)]:
        ax.text(xi, v + 0.02, "%d/%d" % (n, tot), ha="center", va="bottom", fontsize=6.2, color="#333333")

    # the A2 split: what a rigid search could have found vs what it could never reach.
    # Drawn as a bracket to the right of the bars so it never crosses a value label.
    xb = x[1] + w + 0.20
    ax.plot([xb, xb], [43 / 56, 47 / 56], color=C_A2, lw=1.0)
    for yv in (43 / 56, 47 / 56):
        ax.plot([xb - 0.05, xb + 0.05], [yv, yv], color=C_A2, lw=1.0)
    ax.text(xb + 0.09, (47 / 56 + 43 / 56) / 2, "4 contacts no rigid\nplacement can reach",
            fontsize=6.0, color=C_A2, va="center", linespacing=1.2)

    ax.set_xticks(x); ax.set_xticklabels(groups)
    ax.set_ylabel("fnat (fraction of native contacts)")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(-0.55, 1.95)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False,
              handlelength=1.1, ncol=3, columnspacing=1.2)
    save(fig, "Figure6_b")


# ----------------------------------------------------------------- Figure 7a
# Species: per-sample confidence, both models, both ligands. Source: vc11 eval/, cofold2 analysis/.
def fig7a():
    bs = DATA["boltz_species"]
    boltz = {"Vc1.1": (bs["Vc1.1 / ECD only"]["human"], bs["Vc1.1 / ECD only"]["rat"]),
             "RgIA": (bs["RgIA / ECD only"]["human"], bs["RgIA / ECD only"]["rat"])}
    ch = DATA["chai"]
    chai = {"Vc1.1": (ch["E1_HUM_ECD_VC11"], ch["E1_RAT_ECD_VC11"]),
            "RgIA": (ch["E1b_HUM_ECD_RGIA"], ch["E1b_RAT_ECD_RGIA"])}

    fig, axes = plt.subplots(1, 2, figsize=(118 * MM, 52 * MM))
    notes = ["no confidence metric separates\nthe species (1 of 28 tests\nnominal, = chance)",
             "separates both ligands, but the\neffect size runs opposite to the\nmeasured potency gap"]
    for ax, (model, dat, ylab), note in zip(
            axes, [("Boltz-2", boltz, "composite confidence"),
                   ("Chai-1", chai, "aggregate score")], notes):
        pos, labels = [], []
        for k, (lig, (hum, rat)) in enumerate(dat.items()):
            for j, (vals, c) in enumerate([(hum, C_A2), (rat, C_A1)]):
                pp = k * 2.4 + j
                pos.append(pp); labels.append("%s\n%s" % (lig, "human" if j == 0 else "rat"))
                jit = np.linspace(-0.10, 0.10, len(vals))
                ax.scatter(pp + jit, vals, s=13, color=c, zorder=3, linewidths=0)
                ax.plot([pp - 0.22, pp + 0.22], [np.mean(vals)] * 2, color=c, lw=1.5, zorder=4)
        ax.set_xticks(pos); ax.set_xticklabels(labels, fontsize=6.2, linespacing=1.2)
        ax.set_title(model, pad=6)
        ax.set_ylabel(ylab)
        ax.set_xlim(-0.7, max(pos) + 0.7)
        # headroom for the note so it never lands on a point
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + (hi - lo) * 0.34)
        ax.text(0.5, 0.99, note, transform=ax.transAxes, ha="center", va="top",
                fontsize=6.0, color=C_GREY, linespacing=1.25)
    axes[0].legend(handles=[Line2D([], [], marker="o", ls="", color=C_A2, ms=4, label="human"),
                            Line2D([], [], marker="o", ls="", color=C_A1, ms=4, label="rat")],
                   loc="lower left", frameon=False, handlelength=0.8, ncol=2, columnspacing=1.0)
    save(fig, "Figure7_a")


# ----------------------------------------------------------------- Figure 7b
# Binder discrimination at A1: the control that outscores the true binder.
def fig7b():
    order = [("E2_TRUE_IMI", "ImI\ntrue"), ("E2_SHUF", "ImI\nshuffled"),
             ("E2_NB_RGIA", "RgIA\noff-target"), ("E2_NB_MVIIA", "MVIIA\noff-target")]
    src = [("Boltz-2", DATA["boltz_e2"], "composite confidence"),
           ("Chai-1", DATA["chai"], "aggregate score")]
    fig, axes = plt.subplots(1, 2, figsize=(118 * MM, 52 * MM))
    for ax, (model, d, ylab) in zip(axes, src):
        for i, (key, lab) in enumerate(order):
            vals = d.get(key)
            if not vals:
                continue
            c = C_GREY if i == 0 else (C_A2 if i == 2 else "#9C9C9C")
            jit = np.linspace(-0.10, 0.10, len(vals))
            ax.scatter([i] * len(vals) + jit, vals, s=13, color=c, zorder=3, linewidths=0)
            ax.plot([i - 0.22, i + 0.22], [np.mean(vals)] * 2, color=c, lw=1.5, zorder=4)
        tv = d.get("E2_TRUE_IMI")
        if tv:
            ax.axhspan(min(tv), max(tv), color=C_GREY, alpha=0.12, zorder=0, lw=0)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([l for _, l in order], fontsize=6.2, linespacing=1.2)
        ax.set_ylabel(ylab); ax.set_title(model, pad=6)
        ax.set_xlim(-0.6, len(order) - 0.4)
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * 0.22)
    # the load-bearing observation: an off-target ligand outscores the true binder
    bt = DATA["boltz_e2"]
    axes[0].annotate("RgIA outscores\nthe true binder", xy=(2, max(bt["E2_NB_RGIA"])),
                     xytext=(0.30, 0.93), textcoords="axes fraction",
                     fontsize=6.0, color=C_A2, linespacing=1.2, ha="left", va="top",
                     arrowprops=dict(arrowstyle="->", color=C_A2, lw=0.8))
    save(fig, "Figure7_b")


# ----------------------------------------------------------------- Figure 2e
# Disclosure ladder by receptor construct — the sign inversion. Source: round3/ecd_rmsd/, s0_analysis/rmsd/.
def fig2e():
    ladders = {
        "A1 full-length\n(ECD + TMD)": {"BLIND": [8.745, 7.647], "XDON": [9.150, 9.007], "NATIVE": [10.640, 9.834]},
        "A1 ECD only": {"BLIND": [8.445, 10.900], "XDON": [8.624, 8.507], "NATIVE": [7.278, 6.922]},
        "A2 (AChBP, ECD)": {"BLIND": [5.990, 7.168], "XDON": [9.349, 7.807], "NATIVE": [8.181, 8.858]},
    }
    fig, ax = plt.subplots(figsize=(84 * MM, 50 * MM))
    order = ["BLIND", "XDON", "NATIVE"]
    styles = {"A1 full-length\n(ECD + TMD)": ("-", "o"), "A1 ECD only": ("-", "s"),
              "A2 (AChBP, ECD)": ("--", "^")}
    cols = {"A1 full-length\n(ECD + TMD)": C_A1, "A1 ECD only": C_NATIVE, "A2 (AChBP, ECD)": C_A2}
    for name, d in ladders.items():
        means = [np.mean(d[r]) for r in order]
        ls, mk = styles[name]
        ax.plot(range(3), means, ls, marker=mk, ms=5, lw=1.6, color=cols[name], zorder=3,
                markeredgecolor="white", markeredgewidth=0.7)
        for i, r in enumerate(order):
            ax.scatter([i] * len(d[r]), d[r], s=8, color=cols[name], alpha=0.45, zorder=2, linewidths=0)
        ax.text(2.08, means[-1], name, fontsize=6.2, color=cols[name], va="center", linespacing=1.15)
    ax.set_xticks(range(3)); ax.set_xticklabels(order)
    ax.set_xlim(-0.25, 3.55)
    ax.set_xlabel("interface information disclosed to the search")
    # axis inverted so better poses sit higher; say so in the label rather than with a
    # floating annotation, which collides with the y-label at this figure width
    ax.set_ylabel("best ligand RMSD (Å), mean of two GSO\nreplicates — higher is better")
    ax.invert_yaxis()
    save(fig, "Figure2_e")


# ----------------------------------------------------------------- Figure 5a
# Contact map: every native contact x what each method reproduced.
# Source: round2/iface/recovery_contacts_A1_9NX0.csv (one row per native contact).
def fig5a():
    import csv
    rows = list(csv.DictReader(open(os.path.join(PAN, "recovery_contacts_A1_9NX0.csv"))))
    # loop order down the page: principal face first, then complementary, then the rest
    LOOP_ORDER = ["loopA", "loopB", "loopC", "loopD", "loopE", "loopF"]
    def loop_key(r):
        l = r["loop"]
        return (LOOP_ORDER.index(l) if l in LOOP_ORDER else len(LOOP_ORDER), int(r["rec_resnum"]))
    rows.sort(key=loop_key)

    bcols = ["boltz2_m%d" % i for i in range(5)]
    dcols = [c for c in rows[0] if c.startswith("A1_9NX0_") and c.endswith("|best_rmsd")]
    ycols = [c for c in rows[0] if c.startswith("A1_9NX0_") and c.endswith("|top1_score")]

    n = len(rows)
    fig, ax = plt.subplots(figsize=(96 * MM, max(60, n * 2.0) * MM))
    for i, r in enumerate(rows):
        y = n - 1 - i
        nb = sum(int(r[c]) for c in bcols)
        # Boltz-2: shaded by how many of the five samples reproduced it
        ax.add_patch(plt.Rectangle((0, y - 0.42), 1, 0.84,
                                   facecolor=C_FLOOR, alpha=0.18 + 0.82 * nb / 5,
                                   edgecolor="white", lw=0.6))
        if nb:
            ax.text(0.5, y, str(nb), ha="center", va="center", fontsize=5.6,
                    color="white" if nb >= 3 else "#2A2A2A")
        # docking: any of the reported poses, per rung
        for j, rung in enumerate(["BLIND", "XDON", "NATIVE"]):
            got = any(int(r[c]) for c in dcols + ycols if "_%s_" % rung in c)
            ax.add_patch(plt.Rectangle((1.15 + j, y - 0.42), 1, 0.84,
                                       facecolor=RUNG_C[rung] if got else "#F0F0F0",
                                       edgecolor="white", lw=0.6))
        lab = "%s%s–%s%s" % (r["rec_resname"].capitalize(), r["rec_resnum"],
                             r["lig_resname"].capitalize(), r["lig_resnum"])
        ax.text(-0.15, y, lab, ha="right", va="center", fontsize=5.4)
        lp = r["loop"].replace("loop", "")
        if lp.startswith("outside"):
            lp = "–"
        ax.text(4.35, y, lp, ha="left", va="center", fontsize=5.4, color=C_GREY)

    ax.set_xlim(-2.6, 5.4); ax.set_ylim(-0.8, n + 1.4)
    ax.set_xticks([0.5, 1.65, 2.65, 3.65])
    ax.set_xticklabels(["Boltz-2\n(of 5)", "BLIND", "XDON", "NATIVE"], fontsize=6.0, linespacing=1.2)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.text(2.65, n + 1.05, "LightDock", ha="center", fontsize=6.3, color=C_GREY)
    ax.plot([1.15, 4.15], [n + 0.75] * 2, color=C_GREY, lw=0.6)
    ax.text(4.35, n + 0.35, "loop", ha="left", va="center", fontsize=5.6, color=C_GREY)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    save(fig, "Figure5_a")


# ----------------------------------------------------------------- Figure 5c
# What sampling reached vs what scoring handed back, per A1 run.
# A structural superposition would show this too, but the gap is the point and it is a number:
# the pose a user receives is up to 88 A away, on the transmembrane domain.
# Source: s0_analysis/rmsd/A1_9NX0_*.json (best_rmsd, rmsd_of_top1_by_score).
def fig5c():
    runs = [("BLIND rep1", 8.745, 11.472), ("BLIND rep2", 7.647, 87.990),
            ("XDON rep1", 9.150, 15.842), ("XDON rep2", 9.007, 15.691),
            ("NATIVE rep1", 10.640, 11.103), ("NATIVE rep2", 9.834, 21.116)]
    fig, ax = plt.subplots(figsize=(84 * MM, 48 * MM))
    y = np.arange(len(runs))[::-1]
    for yi, (lab, best, top1) in zip(y, runs):
        ax.plot([best, top1], [yi, yi], color="#D5D5D5", lw=1.4, zorder=1,
                solid_capstyle="butt")
        ax.scatter([best], [yi], s=26, color=C_FLOOR, zorder=3, linewidths=0)
        ax.scatter([top1], [yi], s=26, color=C_A2, zorder=3, linewidths=0)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in runs], fontsize=6.3)
    ax.set_xscale("log")
    ax.set_xlabel("ligand RMSD from the native pose (Å, log scale)")
    ax.axvline(3, color=C_GREY, ls=":", lw=0.8, zorder=0)
    ax.text(3, len(runs) - 0.35, " 3 Å", fontsize=5.8, color=C_GREY, va="top")
    ax.annotate("on the transmembrane\ndomain, no shared\ncontact residue",
                xy=(87.99, y[1]), xytext=(0.42, 0.30), textcoords="axes fraction",
                fontsize=6.0, color=C_A2, linespacing=1.2, ha="left", va="center",
                arrowprops=dict(arrowstyle="->", color=C_A2, lw=0.8))
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C_FLOOR, ms=4,
                              label="best pose the search found"),
                       Line2D([], [], marker="o", ls="", color=C_A2, ms=4,
                              label="pose the scoring function returned")],
              loc="upper center", bbox_to_anchor=(0.5, -0.30), frameon=False,
              handlelength=0.8, ncol=1)
    ax.set_title("A1 (9NX0)  all six runs", pad=6, fontsize=7)
    save(fig, "Figure5_c")


if __name__ == "__main__":
    print("writing panels to", PAN)
    fig5a(); fig5b(); fig5c(); fig6a(); fig6b(); fig7a(); fig7b(); fig2e()
    print("done")
