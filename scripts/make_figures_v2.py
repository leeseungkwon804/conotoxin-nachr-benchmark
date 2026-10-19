#!/usr/bin/env /usr/bin/python3
"""
Main figures for the conopeptide modelling paper, v2.

Run with /usr/bin/python3 -- the Homebrew python 3.14 on this machine has a broken
pyexpat and cannot import matplotlib.

Structural panels come from scripts/make_structure_panels.py (PyMOL, headless) and
are composited in here; if a render is missing the panel degrades to a labelled
placeholder rather than failing the whole figure.

Design rules applied throughout:
  * Arial, pdf.fonttype 42 so the vector output stays editable.
  * No red/green pairing anywhere; the categorical palette is Okabe-Ito derived and
    survives greyscale.
  * n = 12 runs and n = 5 samples are small, so every distribution shows the actual
    points on top of the summary.
  * Docking accuracy is always drawn against its rigid-body floor. Raw RMSD is not
    interpretable alone here: A2's all-atom floor (4.453 A) already exceeds the 3 A
    acceptance threshold, so at that anchor the threshold was unreachable by
    construction, and any panel that hid the floor would misrepresent it.
  * LightDock scoring is MAXIMISED, so a scoring function that works gives a NEGATIVE
    Spearman rho against RMSD. Panels involving rho state the sign convention.
  * 7EGR and 9NX0 are never pooled: their data-leakage status differs and every panel
    that shows both keeps them in separate axes or separate groups, labelled.
"""
import csv
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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
OUT = os.path.join(ROOT, "manuscript", "figures_v2")
PAN = os.path.join(OUT, "panels")
os.makedirs(OUT, exist_ok=True)

C_BLIND, C_XDON, C_NATIVE = "#0072B2", "#E69F00", "#CC79A7"
C_GREY, C_FLOOR = "#5A5A5A", "#009E73"
C_A1, C_A2 = "#0072B2", "#D55E00"
RUNG_C = {"BLIND": C_BLIND, "XDON": C_XDON, "NATIVE": C_NATIVE}
RUNGS = ["BLIND", "XDON", "NATIVE"]
MM = 1 / 25.4


def panel_label(ax, s, dx=-0.16, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")


def export_all_panels(fig, gs, panels, prefix, pad=0.06):
    """Save every named panel in `panels` (dict: name -> (axes_list, (row, col)))
    as its own standalone, tightly cropped PNG+PDF.

    Each panel's own natural tight bbox (title, tick labels, panel_label() letter,
    any companion artist such as a colorbar's Axes) is used as-is on any edge that
    has no neighbouring panel -- so a two-line title or a long barh category label
    is never truncated just because it spills into the inter-panel gap; that gap
    exists precisely so such content has somewhere to go.

    Two panels are only clipped against each other where their RAW bboxes
    actually overlap across the gap -- e.g. a colorbar's rotated axis label that
    genuinely reaches into the neighbour's territory. In that case both are
    trimmed to the exact midpoint of the gap, so the same pixel is never drawn in
    both exports. Panels that stay within their own half are left untouched."""
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bottoms, tops, lefts, rights = gs.get_grid_positions(fig)
    nrows, ncols = len(bottoms), len(lefts)
    wpx = fig.get_figwidth() * renderer.dpi
    hpx = fig.get_figheight() * renderer.dpi

    raw = {}
    for name, (axes_list, rc) in panels.items():
        boxes = [a.get_tightbbox(renderer) for a in axes_list]
        raw[name] = (Bbox.union(boxes), rc)
    bnds = {name: list(b.extents) for name, (b, rc) in raw.items()}  # [x0,y0,x1,y1] px
    clipped = {name: [False, False, False, False] for name in panels}  # x0,y0,x1,y1
    by_cell = {}
    for name, (b, rc) in raw.items():
        by_cell.setdefault(rc, []).append(name)

    for name, (b, (row, col)) in raw.items():
        if col + 1 < ncols:  # right neighbour
            mid_x = (rights[col] + lefts[col + 1]) / 2 * wpx
            for other in by_cell.get((row, col + 1), []):
                ob = raw[other][0]
                # clip only where the two RAW bboxes genuinely intersect as
                # rectangles -- crossing the midpoint alone is not enough: a
                # panel's bbox can cross into the next cell's column range
                # while sitting at a row the neighbour's own content never
                # reaches (or vice versa), and bbox_inches crops the shared
                # canvas pixel-for-pixel, so an unclipped panel whose bbox
                # merely overlaps the neighbour's cell -- without the
                # neighbour's own content being there -- is still safe.
                if Bbox.intersection(b, ob) is not None:
                    bnds[name][2] = min(bnds[name][2], mid_x); clipped[name][2] = True
                    bnds[other][0] = max(bnds[other][0], mid_x); clipped[other][0] = True
        if row + 1 < nrows:  # row index increases downward; display y increases upward
            mid_y = (bottoms[row] + tops[row + 1]) / 2 * hpx
            for other in by_cell.get((row + 1, col), []):
                ob = raw[other][0]
                if Bbox.intersection(b, ob) is not None:
                    bnds[name][1] = max(bnds[name][1], mid_y); clipped[name][1] = True
                    bnds[other][3] = min(bnds[other][3], mid_y); clipped[other][3] = True

    # pad_inches is a no-op when bbox_inches is a Bbox instance rather than the
    # literal string "tight" (matplotlib only applies it in "tight" mode), so the
    # padding has to be baked into the bbox extents by hand here -- but only on
    # edges that kept their natural, unclipped extent. Padding a collision-clipped
    # edge would push the crop back past the exact midpoint that was computed to
    # keep the two panels' pixels from ever overlapping.
    pad_px = pad * renderer.dpi
    signs = (-1, -1, 1, 1)
    for name in panels:
        ext = [v + s * pad_px if not c else v
               for v, c, s in zip(bnds[name], clipped[name], signs)]
        full_in = Bbox.from_extents(*ext).transformed(fig.dpi_scale_trans.inverted())
        for ext_ in ("png", "pdf"):
            fig.savefig("%s_%s.%s" % (prefix, name, ext_), dpi=400, bbox_inches=full_in)


def show_render(ax, name, caption=None):
    """Composite a PyMOL render; degrade to a placeholder if it is not there yet."""
    p = os.path.join(PAN, name + ".png")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    if os.path.exists(p):
        ax.imshow(mpimg.imread(p))
    else:
        ax.text(0.5, 0.5, "[%s]\nrender pending" % name, transform=ax.transAxes,
                ha="center", va="center", fontsize=6.5, color="#999999")
    if caption:
        ax.set_title(caption, fontsize=7, pad=3)
    return os.path.exists(p)


def load_main():
    rows = list(csv.DictReader(open(os.path.join(ROOT, "results/s0_ladder/s0_12run_main.csv"))))
    for r in rows:
        for k, v in list(r.items()):
            if k not in ("run", "anchor", "rung", "complete"):
                try:
                    r[k] = float(v)
                except (TypeError, ValueError):
                    pass
    return rows


def per_swarm_best(run):
    z = np.load(os.path.join(ROOT, "pod_mirror/s0_analysis/rmsd", run + ".npz"))
    sw, rm = z["swarm"], z["rmsd"]
    o = np.argsort(sw, kind="mergesort")
    sw, rm = sw[o], rm[o]
    return np.array([s.min() for s in np.split(rm, np.flatnonzero(np.diff(sw)) + 1)], float)


# =============================================================================== Fig 1
def figure1():
    fig = plt.figure(figsize=(170 * MM, 210 * MM))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 1.0, 0.95],
                          hspace=0.34, wspace=0.16,
                          left=0.07, right=0.98, top=0.96, bottom=0.06)

    PANELS = {}

    # (a) each anchor in Top and Side orientation (2x2 sub-block: rows 0-1 x anchor columns)
    ax = fig.add_subplot(gs[0, 0])
    show_render(ax, "1a_A1_9NX0_top", "A1  9NX0 — human $\\alpha$7 nAChR + $\\alpha$-conotoxin ImI\ncryo-EM 3.06 A, released 2026-02-04\nTop")
    panel_label(ax, "a", dx=-0.05)
    PANELS["a_A1_top"] = ([ax], (0, 0))
    ax = fig.add_subplot(gs[0, 1])
    show_render(ax, "1a_A2_7EGR_top", "A2  7EGR — $\\it{Aplysia}$ AChBP + $\\alpha$-conotoxin RgIA\nX-ray 2.50 A, released 2022-04-06\nTop")
    PANELS["a_A2_top"] = ([ax], (0, 1))
    ax = fig.add_subplot(gs[1, 0])
    show_render(ax, "1a_A1_9NX0_side", "Side")
    PANELS["a_A1_side"] = ([ax], (1, 0))
    ax = fig.add_subplot(gs[1, 1])
    show_render(ax, "1a_A2_7EGR_side", "Side")
    PANELS["a_A2_side"] = ([ax], (1, 1))

    ax = fig.add_subplot(gs[2, 0])
    show_render(ax, "1b_A1_9NX0", "Orthosteric site, (+)/(−) subunit pair\n28 receptor residues within 5 A")
    panel_label(ax, "b", dx=-0.05)
    PANELS["b_A1_9NX0"] = ([ax], (2, 0))
    ax = fig.add_subplot(gs[2, 1])
    show_render(ax, "1b_A2_7EGR", "Orthosteric site, (+)/(−) subunit pair\n33 receptor residues within 5 A")
    PANELS["b_A2_7EGR"] = ([ax], (2, 1))

    # (c) the three information rungs and how much of the true interface each supplies
    ax = fig.add_subplot(gs[3, 0])
    recall = {"BLIND": (0.0, 0.0), "XDON": (0.893, 0.364), "NATIVE": (1.0, 1.0)}
    x = np.arange(3); w = 0.36
    ax.bar(x - w / 2, [recall[r][0] for r in RUNGS], w, color=C_A1, label="A1 (9NX0)")
    ax.bar(x + w / 2, [recall[r][1] for r in RUNGS], w, color=C_A2, label="A2 (7EGR)")
    for i, r in enumerate(RUNGS):
        for j, off in enumerate((-w / 2, w / 2)):
            v = recall[r][j]
            ax.text(i + off, v + 0.03, "%.2f" % v if v else "0", ha="center", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(["BLIND\nnone", "XDON\nfrom another\ncomplex", "NATIVE\nthe target's\nown contacts"])
    ax.set_ylabel("Fraction of the native contact\nset supplied as restraints")
    ax.set_ylim(0, 1.18)
    ax.legend(frameon=False, loc="upper left", handletextpad=0.5)
    ax.set_title("Information supplied to the search", fontsize=7.5, pad=4)
    panel_label(ax, "c")
    PANELS["c"] = ([ax], (3, 0))

    # (d) apo vs bound conformer -- the rigid-body limit
    ax = fig.add_subplot(gs[3, 1])
    ok = show_render(ax, "1d_A2_RgIA",
                     "Docked apo conformer (teal) vs bound (orange)\nthe residual limits any rigid-body pose")
    if ok:
        ax.text(0.5, -0.06, "A1 floor 2.44 A   ·   A2 floor 4.45 A  (all-atom)",
                transform=ax.transAxes, ha="center", va="top", fontsize=6.5)
    panel_label(ax, "d", dx=-0.05, dy=1.24)
    PANELS["d"] = ([ax], (3, 1))

    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, "Figure1." + e), dpi=400, bbox_inches="tight")
    export_all_panels(fig, gs, PANELS, os.path.join(PAN, "Figure1"))
    plt.close(fig)
    print("Figure1 written (%d individual panels)" % len(PANELS))


# =============================================================================== Fig 2
def figure2():
    rows = load_main()
    fig = plt.figure(figsize=(170 * MM, 148 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.30,
                          left=0.10, right=0.97, top=0.92, bottom=0.09)
    PANELS = {}

    # (a) excess over floor, drawn as a stem from the floor to the achieved value
    ax = fig.add_subplot(gs[0, 0])
    xt, xl, x = [], [], 0.0
    for anc, pdb in (("A1", "9NX0"), ("A2", "7EGR")):
        fl = [r for r in rows if r["anchor"] == anc][0]["floor_allatom_A"]
        xs0 = x
        for rung in RUNGS:
            rs = [r for r in rows if r["anchor"] == anc and r["rung"] == rung]
            for r in rs:
                ax.vlines(x, fl, r["best_allatom_A"], color=RUNG_C[rung], lw=2.4,
                          alpha=0.85, zorder=2)
                ax.plot([x], [r["best_allatom_A"]], "o", ms=4.2, color=RUNG_C[rung],
                        mec="white", mew=0.6, zorder=3)
                x += 0.44
            xt.append(x - 0.44 * len(rs) / 2 - 0.22); xl.append(rung)
            x += 0.50
        ax.hlines(fl, xs0 - 0.30, x - 0.72, color=C_FLOOR, lw=1.6, zorder=4)
        ax.text(xs0 - 0.28, fl - 0.42, "floor %.2f" % fl, color=C_FLOOR, fontsize=6,
                va="top", ha="left")
        x += 1.25
    ax.axhline(3.0, color=C_GREY, ls=(0, (4, 3)), lw=1.0, zorder=1)
    ax.text(x - 1.3, 3.0, "3 A ", color=C_GREY, fontsize=6, va="bottom", ha="right")
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=6, rotation=30, ha="right")
    ax.set_ylabel("Best ligand RMSD, all-atom (A)")
    ax.set_ylim(0, 12.8); ax.set_xlim(-0.8, x - 1.5)
    ax.text(np.mean(xt[:3]), 12.4, "A1  9NX0", ha="center", fontsize=6.8, fontweight="bold")
    ax.text(np.mean(xt[3:]), 12.4, "A2  7EGR", ha="center", fontsize=6.8, fontweight="bold")
    ax.text(0.015, 0.90, "bar height = excess over\nthe rigid-body floor",
            transform=ax.transAxes, fontsize=6, va="top", color="#333333")
    panel_label(ax, "a")
    PANELS["a"] = ([ax], (0, 0))

    # (b) per-swarm-best distributions
    ax = fig.add_subplot(gs[0, 1])
    data, cols, xt2, xl2, pos = [], [], [], [], 0
    for anc, pdb in (("A1", "9NX0"), ("A2", "7EGR")):
        for rung in RUNGS:
            data.append(np.concatenate([per_swarm_best("%s_%s_%s_rep%d" % (anc, pdb, rung, k))
                                        for k in (1, 2)]))
            cols.append(RUNG_C[rung]); xt2.append(pos); xl2.append(rung); pos += 1
        pos += 0.8
    parts = ax.violinplot(data, positions=xt2, widths=0.80, showextrema=False)
    for b, c in zip(parts["bodies"], cols):
        b.set_facecolor(c); b.set_alpha(0.35); b.set_edgecolor(c); b.set_linewidth(0.7)
    ax.boxplot(data, positions=xt2, widths=0.19, showfliers=False, patch_artist=True,
               medianprops=dict(color="black", lw=1.0),
               boxprops=dict(facecolor="white", lw=0.7),
               whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    ax.set_xticks(xt2); ax.set_xticklabels(xl2, fontsize=6, rotation=30, ha="right")
    ax.set_yscale("log")
    ax.set_ylabel("Per-swarm best ligand RMSD (A)")
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi * 2.1)
    ax.text(np.mean(xt2[:3]), hi * 1.45, "A1  9NX0", ha="center", fontsize=6.8, fontweight="bold")
    ax.text(np.mean(xt2[3:]), hi * 1.45, "A2  7EGR", ha="center", fontsize=6.8, fontweight="bold")
    ax.text(0.5, 0.015, "not one of the 3,144 scored swarms reached 3 A",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=6, color="#333333")
    panel_label(ax, "b")
    PANELS["b"] = ([ax], (0, 1))

    # (c) BLIND vs NATIVE, paired
    ax = fig.add_subplot(gs[1, 0])
    for anc, pdb in (("A1", "9NX0"), ("A2", "7EGR")):
        for rep in (1, 2):
            for basis, key in (("all-atom", "best_allatom_A"), ("backbone", "best_backbone_A")):
                b = [r for r in rows if r["run"] == "%s_%s_BLIND_rep%d" % (anc, pdb, rep)][0][key]
                n = [r for r in rows if r["run"] == "%s_%s_NATIVE_rep%d" % (anc, pdb, rep)][0][key]
                c = C_A1 if anc == "A1" else C_A2
                ls = "-" if basis == "all-atom" else (0, (2, 1.5))
                ax.plot([0, 1], [b, n], color=c, lw=1.0, ls=ls, alpha=0.9, zorder=2)
                ax.plot([0, 1], [b, n], "o", ms=3.3, color=c, mec="white", mew=0.4, zorder=3)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["BLIND\nno interface\ninformation", "NATIVE\nthe target's own\ncontacts supplied"])
    ax.set_xlim(-0.40, 1.40)
    ax.set_ylabel("Best ligand RMSD (A)")
    ax.set_title("Supplying the interface never improved the best pose\n8 of 8 comparisons, margin 1.36–2.99 A",
                 fontsize=7.5, pad=5)
    ax.legend(handles=[Line2D([], [], color=C_A1, lw=1.2, label="A1 (9NX0)"),
                       Line2D([], [], color=C_A2, lw=1.2, label="A2 (7EGR)"),
                       Line2D([], [], color="k", lw=1.0, label="all-atom"),
                       Line2D([], [], color="k", lw=1.0, ls=(0, (2, 1.5)), label="backbone")],
              loc="upper left", frameon=False, ncol=2, columnspacing=0.7, handletextpad=0.4)
    panel_label(ax, "c")
    PANELS["c"] = ([ax], (1, 0))

    # (d) apo entry floors
    ax = fig.add_subplot(gs[1, 1])
    rowsE = json.load(open(os.path.join(ROOT, "pod_mirror_t0/ensemble_floor_full.json")))["9NX0_ligand_I.pdb"]
    groups = {k: [r["allatom"] for r in rowsE if "_%s_m" % k in r["file"]]
              for k in ("1G2G", "1IMI", "1CNL")}
    order2 = ["1G2G", "1IMI", "1CNL"]; cols2 = ["#0072B2", "#5A5A5A", "#D55E00"]
    rng = np.random.default_rng(0)
    for i, k in enumerate(order2):
        v = np.array(groups[k], float)
        ax.scatter(np.full(len(v), i) + rng.uniform(-0.17, 0.17, len(v)), v,
                   s=10, color=cols2[i], alpha=0.8, linewidths=0, zorder=3)
        ax.hlines(np.median(v), i - 0.28, i + 0.28, color="black", lw=1.3, zorder=4)
        ax.text(i, 5.15, "n=%d" % len(v), ha="center", fontsize=6)
    ax.axhline(3.0, color=C_GREY, ls=(0, (4, 3)), lw=1.0)
    ax.text(2.48, 3.05, "3 A", color=C_GREY, fontsize=6, va="bottom", ha="right")
    ax.annotate("", xy=(2.30, 1.022), xytext=(2.30, 4.574),
                arrowprops=dict(arrowstyle="<->", lw=0.9, color="#333333"))
    ax.text(2.24, 2.8, "1.02–4.57 A", rotation=90, fontsize=6, ha="right", va="center")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["1G2G\n$\\alpha$7-recognition", "1IMI\nused here", "1CNL"])
    ax.set_ylabel("Rigid-body floor, all-atom (A)")
    ax.set_ylim(0, 5.6); ax.set_xlim(-0.5, 2.6)
    ax.set_title("The apo NMR entry alone moves the achievable floor\nby more than the threshold itself",
                 fontsize=7.5, pad=5)
    panel_label(ax, "d")
    PANELS["d"] = ([ax], (1, 1))

    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, "Figure2." + e), dpi=400, bbox_inches="tight")
    export_all_panels(fig, gs, PANELS, os.path.join(PAN, "Figure2"))
    plt.close(fig)
    print("Figure2 written (%d individual panels)" % len(PANELS))


# =============================================================================== Fig 3
def figure3():
    panel = list(csv.DictReader(open(os.path.join(ROOT, "pod_mirror_t0/multiscore_panel.csv"))))
    summ = list(csv.DictReader(open(os.path.join(ROOT, "pod_mirror_t0/multiscore_summary.csv"))))
    for r in panel:
        for k in ("spearman_random_subsample", "rmsd_of_top1_by_score_A", "best_rmsd_in_set_A",
                  "best_rmsd_in_top10_A", "best_rmsd_in_top100_A"):
            r[k] = float(r[k])
    funcs = [s["scoring_function"] for s in sorted(summ, key=lambda s: float(s["rho_median"]))]
    runs = ["%s_%s_%s_rep%d" % (a, p, g, k)
            for a, p in (("A1", "9NX0"), ("A2", "7EGR")) for g in RUNGS for k in (1, 2)]

    fig = plt.figure(figsize=(170 * MM, 150 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.65,
                          left=0.14, right=0.96, top=0.90, bottom=0.10)
    PANELS = {}

    # (a) rho heatmap
    ax = fig.add_subplot(gs[0, 0])
    M = np.full((len(funcs), len(runs)), np.nan)
    for r in panel:
        M[funcs.index(r["scoring_function"]), runs.index(r["run"])] = r["spearman_random_subsample"]
    v = np.nanmax(np.abs(M))
    im = ax.imshow(M, cmap="PuOr_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels([r.split("_")[2][0] + r[-1] for r in runs], fontsize=5.5)
    ax.set_yticks(range(len(funcs))); ax.set_yticklabels(funcs, fontsize=6.5)
    ax.axvline(5.5, color="white", lw=1.6)
    ax.text(2.5, -1.35, "A1 (9NX0)", ha="center", fontsize=6.5, fontweight="bold")
    ax.text(8.5, -1.35, "A2 (7EGR)", ha="center", fontsize=6.5, fontweight="bold")
    ax.set_xlabel("run   (B/X/N = BLIND/XDON/NATIVE, 1–2 = seed)", fontsize=6.3)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Spearman $\\rho$  (score vs RMSD)", fontsize=6.3)
    cb.ax.tick_params(labelsize=6)
    cb.ax.annotate("works", xy=(1.0, 0.13), xytext=(3, 0), textcoords="offset points",
                   xycoords="axes fraction", fontsize=6, color=C_GREY, rotation=90,
                   ha="left", va="center")
    ax.set_title("LightDock scores are maximised, so a working function\ngives negative $\\rho$ (purple)",
                 fontsize=7.5, pad=20)
    panel_label(ax, "a", dx=-0.34)
    PANELS["a"] = ([ax, cb.ax], (0, 0))

    # (b) top-1 RMSD by function
    ax = fig.add_subplot(gs[0, 1])
    data = [[r["rmsd_of_top1_by_score_A"] for r in panel if r["scoring_function"] == f] for f in funcs]
    ax.boxplot(data, vert=False, widths=0.6, showfliers=False, patch_artist=True,
               medianprops=dict(color="black", lw=1.0),
               boxprops=dict(facecolor="#E4E4E4", lw=0.7),
               whiskerprops=dict(lw=0.7), capprops=dict(lw=0.7))
    rng = np.random.default_rng(1)
    for i, f in enumerate(funcs, start=1):
        sub = [r for r in panel if r["scoring_function"] == f]
        ax.scatter([r["rmsd_of_top1_by_score_A"] for r in sub],
                   np.full(len(sub), i) + rng.uniform(-0.2, 0.2, len(sub)),
                   s=8, c=[C_A1 if r["anchor"] == "A1" else C_A2 for r in sub],
                   alpha=0.85, linewidths=0, zorder=3)
    ax.axvline(7.647, color=C_A1, ls=(0, (3, 2)), lw=1.0)
    ax.axvline(5.990, color=C_A2, ls=(0, (3, 2)), lw=1.0)
    ax.set_yticks(range(1, len(funcs) + 1)); ax.set_yticklabels(funcs, fontsize=6.5)
    ax.set_xscale("log"); ax.set_xlabel("RMSD of the top-scoring pose (A)")
    ax.set_title("Dashed lines mark the best pose that was actually\npresent in each anchor's scored set",
                 fontsize=7.5, pad=5)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=C_A1, label="A1 run", ms=4),
                       Line2D([], [], marker="o", ls="", color=C_A2, label="A2 run", ms=4)],
              loc="lower right", frameon=False, handletextpad=0.4)
    panel_label(ax, "b", dx=-0.30)
    PANELS["b"] = ([ax], (0, 1))

    # (c) enrichment: does looking deeper into the ranking find the good poses?
    ax = fig.add_subplot(gs[1, 0])
    depths = [1, 10, 100]
    RUN_C = "A2_7EGR_BLIND_rep1"   # the run that produced the campaign's best pose
    for f in funcs:
        sub = [r for r in panel if r["scoring_function"] == f and r["run"] == RUN_C]
        if not sub:
            continue
        r0 = sub[0]
        y = [r0["rmsd_of_top1_by_score_A"], r0["best_rmsd_in_top10_A"], r0["best_rmsd_in_top100_A"]]
        hi = f in ("fastdfire", "dfire", "cpydock")
        ax.plot(depths, y, "-o", ms=3.4, lw=1.5 if hi else 0.8,
                color={"fastdfire": C_A2, "dfire": C_A2, "cpydock": C_NATIVE}.get(f, "#BBBBBB"),
                alpha=1.0 if hi else 0.8, zorder=3 if hi else 2,
                label=("dfire / fastdfire (coincide)" if f == "fastdfire"
                       else (None if f == "dfire" else f)) if hi else None)
    smin = float([r["best_rmsd_in_set_A"] for r in panel if r["run"] == RUN_C][0])
    ax.axhline(smin, color=C_FLOOR, lw=1.4)
    ax.text(1.05, smin * 1.04, "best pose this run produced (%.2f A)" % smin,
            fontsize=6, color=C_FLOOR, va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(depths); ax.set_xticklabels(["top 1", "top 10", "top 100"])
    ax.set_xlabel("How deep into the score ranking we look")
    ax.set_ylabel("Best ligand RMSD found (A)")
    ax.set_title("A2 BLIND rep1: going 100 poses deep still does not\nreach the good pose this run had already produced",
                 fontsize=7.5, pad=5)
    ax.legend(frameon=False, loc="upper right", handletextpad=0.5)
    panel_label(ax, "c")
    PANELS["c"] = ([ax], (1, 0))

    # (d) rank-1 instability
    ax = fig.add_subplot(gs[1, 1])
    rf = list(csv.DictReader(open(os.path.join(ROOT, "results/s0_rescore12/s0_restraintfree_12run.csv"))))
    for r in rf:
        r["restraints_R"] = int(r["restraints_R"])
        r["published_rank_of_restraint_free_winner"] = int(r["published_rank_of_restraint_free_winner"])
    rf.sort(key=lambda r: r["published_rank_of_restraint_free_winner"])
    ys = np.arange(len(rf))
    ax.barh(ys, [r["published_rank_of_restraint_free_winner"] for r in rf],
            color=[C_BLIND if r["restraints_R"] == 0 else C_NATIVE for r in rf], height=0.66)
    ax.set_yticks(ys)
    ax.set_yticklabels([r["run"].replace("_9NX0", "").replace("_7EGR", "") for r in rf], fontsize=5.8)
    ax.set_xscale("log"); ax.axvline(1, color="black", lw=1.0)
    ax.set_xlabel("Published rank of the restraint-free winner")
    ax.set_title("Removing one non-physical term changes which pose\nranks first in every restrained run",
                 fontsize=7.5, pad=5)
    ax.legend(handles=[Patch(facecolor=C_BLIND, label="BLIND (no restraints)"),
                       Patch(facecolor=C_NATIVE, label="restrained (XDON / NATIVE)")],
              loc="lower right", frameon=False, handletextpad=0.5)
    panel_label(ax, "d", dx=-0.30)
    PANELS["d"] = ([ax], (1, 1))

    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, "Figure3." + e), dpi=400, bbox_inches="tight")
    export_all_panels(fig, gs, PANELS, os.path.join(PAN, "Figure3"))
    plt.close(fig)
    print("Figure3 written (%d individual panels)" % len(PANELS))


# =============================================================================== Fig 4
def figure4():
    ms = list(csv.DictReader(open(os.path.join(ROOT, "pod_mirror_boltz/eval_ms/summary.csv"))))
    for r in ms:
        for k in ("lig_allatom_iface_A", "fnat", "confidence"):
            r[k] = float(r[k])
    b12 = json.load(open(os.path.join(ROOT, "pod_mirror_n2/eval_7EGR_b1_vs_b2.json")))
    rows = load_main()
    pan = list(csv.DictReader(open(os.path.join(ROOT, "pod_mirror_t0/multiscore_panel.csv"))))

    fig = plt.figure(figsize=(170 * MM, 152 * MM))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.60,
                          left=0.10, right=0.97, top=0.93, bottom=0.09)
    PANELS = {}

    # (a) prediction vs deposited
    ax = fig.add_subplot(gs[0, 0])
    show_render(ax, "4a_A2_7EGR",
                "Boltz-2 prediction (purple) on the deposited\ntoxin (orange), after receptor superposition")
    panel_label(ax, "a", dx=-0.05)
    PANELS["a"] = ([ax], (0, 0))

    # (b) method comparison, anchors kept apart
    ax = fig.add_subplot(gs[0, 1])
    y, yt, ylab = 0.0, [], []
    for anc, lab in (("A2", "A2  7EGR"), ("A1", "A1  9NX0")):
        dock = min(r["best_allatom_A"] for r in rows if r["anchor"] == anc)
        best_sf = min(float(r["rmsd_of_top1_by_score_A"]) for r in pan if r["anchor"] == anc)
        items = [("docking, best of 12 runs", dock, C_GREY, None),
                 ("best of 9 scoring functions", best_sf, "#9E9E9E", None)]
        if anc == "A2":
            v = [d["lig_allatom_iface_A"] for d in b12 if d["model_set"] == "Boltz-1"]
            items.append(("Boltz-1 co-folding", float(np.median(v)), C_FLOOR, v))
        if anc == "A2":
            # same five poses as pod_mirror_boltz/eval_ms (verified: 0.006 A after
            # rigid superposition, indices permuted). Taking the copy that was scored
            # by the same script as Boltz-1 keeps panels b and c on one footing, and
            # avoids any suggestion that this is an independent repeat -- n = 5, not 10.
            v2 = [d["lig_allatom_iface_A"] for d in b12 if d["model_set"] == "Boltz-2"]
        else:
            v2 = [r["lig_allatom_iface_A"] for r in ms if r["anchor"].startswith(anc)]
        items.append(("Boltz-2 co-folding", float(np.median(v2)), C_A2 if anc == "A2" else C_A1, v2))
        for name, val, c, pts in items:
            ax.barh(y, val, color=c, height=0.62, zorder=2)
            if pts:
                ax.scatter(pts, np.full(len(pts), y), s=8, color="black", zorder=4,
                           alpha=0.85, linewidths=0)
            yt.append(y); ylab.append(name); y += 1
        ax.text(10.8, y - len(items) / 2 - 0.5, lab, fontsize=6.8, fontweight="bold",
                rotation=90, va="center", ha="center")
        y += 0.9
    ax.axvline(3.0, color=C_GREY, ls=(0, (4, 3)), lw=1.0, zorder=1)
    ax.set_yticks(yt); ax.set_yticklabels(ylab, fontsize=6.2)
    ax.set_xlabel("Ligand RMSD (A)"); ax.set_xlim(0, 11.6)
    ax.set_title("Anchors are never pooled: their leakage status differs\nBoltz-2 n = 5 samples per anchor",
                 fontsize=7.5, pad=5)
    panel_label(ax, "b", dx=-0.42)
    PANELS["b"] = ([ax], (0, 1))

    # (c) temporal control
    ax = fig.add_subplot(gs[1, 0])
    for i, label in enumerate(["Boltz-1", "Boltz-2"]):
        v = [d["lig_allatom_iface_A"] for d in b12 if d["model_set"] == label]
        ax.scatter(np.full(len(v), i) + np.linspace(-0.10, 0.10, len(v)), v, s=20,
                   color=C_FLOOR if label == "Boltz-1" else C_A2, zorder=3, linewidths=0)
        ax.hlines(np.median(v), i - 0.22, i + 0.22, color="black", lw=1.3, zorder=4)
    dock2 = min(r["best_allatom_A"] for r in rows if r["anchor"] == "A2")
    ax.axhline(dock2, color=C_GREY, lw=1.1)
    ax.text(-0.42, dock2 + 0.10, "best of 3,144 docking swarms (%.2f A)" % dock2,
            fontsize=6, color=C_GREY, va="bottom")
    ax.axhline(3.0, color=C_GREY, ls=(0, (4, 3)), lw=1.0)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Boltz-1\ncutoff 2021-09-30\ncould not have seen 7EGR",
                        "Boltz-2\ncutoff 2023-06-01\n7EGR inside the window"], fontsize=6.2)
    ax.set_xlim(-0.5, 1.5); ax.set_ylim(0, 7.0)
    ax.set_ylabel("Ligand RMSD (A)")
    ax.set_title("Temporal control at 7EGR, 5 samples each", fontsize=7.5, pad=5)
    panel_label(ax, "c")
    PANELS["c"] = ([ax], (1, 0))

    # (d) leakage timeline
    ax = fig.add_subplot(gs[1, 1])
    ents = [("2BR8", 2005.43, "PnIA"), ("2BYP", 2005.76, "ImI"), ("2C9T", 2006.12, "ImI"),
            ("2UZ6", 2007.60, "TxIA"), ("4EZ1", 2013.35, "BuIA"), ("5JME", 2016.36, "PeIA"),
            ("5CO5", 2016.55, "GIC"), ("5XGL", 2018.45, "LvIA"), ("6HY7", 2019.39, "RgIA"),
            ("7N43", 2021.92, "OmIA"), ("7EGR", 2022.26, "RgIA"), ("9NX0", 2026.09, "ImI")]
    for i, (pdb, yr, tox) in enumerate(ents):
        tgt = pdb in ("7EGR", "9NX0")
        c = C_A2 if pdb == "7EGR" else (C_A1 if pdb == "9NX0" else C_GREY)
        ax.scatter([yr], [i], s=36 if tgt else 15, color=c, zorder=4,
                   linewidths=0.7 if tgt else 0, edgecolors="black")
        left = yr > 2021.0
        gap = 0.55 if tgt else 0.35
        ax.text(yr + (-gap if left else gap), i, "%s (%s)" % (pdb, tox), fontsize=5.6,
                va="center", ha="right" if left else "left",
                fontweight="bold" if tgt else "normal",
                color="black" if tgt else "#555555", zorder=5)
    ax.axvline(2021.75, color=C_FLOOR, lw=1.3, zorder=1)
    ax.axvline(2023.42, color="#8B0000", lw=1.3, zorder=1)
    # the two cutoffs sit only 1.67 yr apart, too narrow for two centred multi-line
    # labels -- point each label outward (away from the other line) instead.
    ax.text(2021.75 - 0.25, len(ents) + 0.5, "Boltz-1\ncutoff", fontsize=5.8, color=C_FLOOR,
            ha="right", va="bottom")
    ax.text(2023.42 + 0.25, len(ents) + 0.5, "Boltz-2\ncutoff", fontsize=5.8, color="#8B0000",
            ha="left", va="bottom")
    ax.set_yticks([]); ax.set_xlabel("PDB initial release date")
    ax.set_xlim(2002.6, 2028.4); ax.set_ylim(-1.1, len(ents) + 2.6)
    ax.set_title("Every conotoxin–receptor complex except our two targets\npredates the Boltz-2 cutoff",
                 fontsize=7.5, pad=5)
    panel_label(ax, "d", dx=-0.08)
    PANELS["d"] = ([ax], (1, 1))

    for e in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, "Figure4." + e), dpi=400, bbox_inches="tight")
    export_all_panels(fig, gs, PANELS, os.path.join(PAN, "Figure4"))
    plt.close(fig)
    print("Figure4 written (%d individual panels)" % len(PANELS))


if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4()
    print("\noutput:", OUT)
    for f in sorted(glob.glob(os.path.join(OUT, "Figure*.png"))):
        print("  %-16s %6.0f KB" % (os.path.basename(f), os.path.getsize(f) / 1024))
