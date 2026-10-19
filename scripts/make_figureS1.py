#!/usr/bin/env /usr/bin/python3
"""
Figure S1 -- Score-RMSD rank correlation computed on the full re-scored subset.

This is the supplementary companion to main-text Figure 3A: same 9 (scoring
function) x 12 (run) grid, same run ordering, same diverging colormap and sign
convention, but the cell value is `spearman_all` (computed on the full re-scored
subset per run, n = 2,742-2,961 poses) instead of `spearman_random_subsample`
(the unbiased-estimator column Figure 3A uses, n = 2,000 seeded random poses/run).

Run with /usr/bin/python3 -- see scripts/make_figures_v2.py header for why.

Design rules inherited from make_figures_v2.py (read there for the full
rationale): Arial + pdf.fonttype 42, PuOr_r diverging colormap, LightDock scores
are MAXIMISED so a working function gives a NEGATIVE Spearman rho (purple),
7EGR (A2, HIGH leakage risk) and 9NX0 (A1, MEDIUM leakage risk) are kept in
separate blocks of the same axis, never pooled into one statistic.

Source: pod_mirror_t0/multiscore_panel.csv, columns spearman_all and
spearman_random_subsample; pod_mirror_t0/multiscore_summary.csv, column
rho_median (used only to fix the y-axis order, identical to Figure 3A).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

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
os.makedirs(OUT, exist_ok=True)
C_GREY = "#5A5A5A"
MM = 1 / 25.4
RUNGS = ["BLIND", "XDON", "NATIVE"]


def main():
    panel = list(csv.DictReader(open(os.path.join(ROOT, "pod_mirror_t0/multiscore_panel.csv"))))
    summ = list(csv.DictReader(open(os.path.join(ROOT, "pod_mirror_t0/multiscore_summary.csv"))))
    for r in panel:
        for k in ("spearman_all", "spearman_random_subsample", "n_scored"):
            r[k] = float(r[k])

    # identical y-axis order to Figure 3A: ascending median of spearman_random_subsample
    funcs = [s["scoring_function"] for s in sorted(summ, key=lambda s: float(s["rho_median"]))]
    runs = ["%s_%s_%s_rep%d" % (a, p, g, k)
            for a, p in (("A1", "9NX0"), ("A2", "7EGR")) for g in RUNGS for k in (1, 2)]

    n_lo = min(r["n_scored"] for r in panel)
    n_hi = max(r["n_scored"] for r in panel)
    print("n_scored range across the 12 runs: %d-%d" % (n_lo, n_hi))

    M_all = np.full((len(funcs), len(runs)), np.nan)
    M_sub = np.full((len(funcs), len(runs)), np.nan)
    for r in panel:
        i, j = funcs.index(r["scoring_function"]), runs.index(r["run"])
        M_all[i, j] = r["spearman_all"]
        M_sub[i, j] = r["spearman_random_subsample"]

    assert np.isfinite(M_all).all() and np.isfinite(M_sub).all(), "missing cell in the 9x12 grid"

    # verify the two counts this figure's caption depends on, directly from the data
    n_cells = M_all.size
    n_sign_flip = int(np.sum(np.sign(M_all) != np.sign(M_sub)))
    n_larger = int(np.sum(np.abs(M_all) > np.abs(M_sub)))
    print("cells total: %d" % n_cells)
    print("sign flips vs Figure 3A (spearman_random_subsample): %d/%d" % (n_sign_flip, n_cells))
    print("|rho| larger on the full subset than on the random subsample: %d/%d" % (n_larger, n_cells))

    v = float(np.nanmax(np.abs(M_all)))
    print("max |spearman_all| in this panel (colour-scale limit): %.4f" % v)

    fig = plt.figure(figsize=(150 * MM, 108 * MM))
    ax = fig.add_axes([0.28, 0.20, 0.62, 0.62])

    im = ax.imshow(M_all, cmap="PuOr_r", vmin=-v, vmax=v, aspect="auto")
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels([r.split("_")[2][0] + r[-1] for r in runs], fontsize=6.5)
    ax.set_yticks(range(len(funcs)))
    ax.set_yticklabels(funcs, fontsize=7)
    ax.axvline(5.5, color="white", lw=1.8)
    ax.text(2.5, -1.05, "A1 (9NX0) -- MEDIUM leakage risk", ha="center", fontsize=6.5,
            fontweight="bold")
    ax.text(8.5, -1.05, "A2 (7EGR) -- HIGH leakage risk", ha="center", fontsize=6.5,
            fontweight="bold")
    ax.set_xlabel("run   (B/X/N = BLIND/XDON/NATIVE, 1-2 = GSO seed rep)", fontsize=6.8)

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Spearman $\\rho$ (score vs. RMSD), full re-scored subset\n"
                  "(n = %d-%d poses per run)" % (n_lo, n_hi), fontsize=6.3)
    cb.ax.tick_params(labelsize=6)
    cb.ax.annotate("works", xy=(1.9, 0.13), xycoords="axes fraction", fontsize=6,
                   color=C_GREY, rotation=90, ha="center", va="center")

    ax.set_title("LightDock scores are maximised, so a working function gives\n"
                  "negative $\\rho$ (purple) -- full-subset estimate, cf. Figure 3A",
                  fontsize=7.5, pad=22)

    caption = (
        "Figure S1. Score–RMSD rank correlation on the full re-scored subset. "
        "Same grid, run order, and sign convention as Figure 3A; here the cell value is "
        "spearman_all, the full re-scored subset per run (subset construction: Methods §4.3), "
        "not the random-subsample column Figure 3A uses. n = %d–%d poses/run (Table S8). "
        "%d of %d cells differ in sign from Figure 3A; $|\\rho|$ is larger here in %d of %d "
        "cells. A1/A2, BLIND/XDON/NATIVE: Supplementary Materials Guide §0. Not a leakage "
        "assessment (Figure 4d)."
    ) % (n_lo, n_hi, n_sign_flip, n_cells, n_larger, n_cells)
    fig.text(0.03, -0.06, caption, fontsize=6.1, va="top", ha="left", wrap=True,
             transform=fig.transFigure)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS1." + e)
        fig.savefig(p, dpi=400, bbox_inches="tight")
    plt.close(fig)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS1." + e)
        print("wrote %-40s %8.1f KB" % (p, os.path.getsize(p) / 1024))

    return n_sign_flip, n_cells, n_larger, n_lo, n_hi


if __name__ == "__main__":
    main()
