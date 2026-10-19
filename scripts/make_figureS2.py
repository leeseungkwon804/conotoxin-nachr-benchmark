#!/usr/bin/env /usr/bin/python3
"""
Figure S2 -- Score-RMSD funnels for the remaining seven scoring functions on a
single run.

Companion to main-text Figure 3's DFIRE/cpydock design (spec_fig3.md section 5,
"panel C"): the same 2,919 poses of run A2_7EGR_BLIND_rep1, scored by each of
the seven scoring functions that are NOT the two singled out in that design
(DFIRE, the most negative-rho / best-behaved function, and cpydock, the most
positive-rho / worst-behaved function -- see spec_fig3.md lines 171-175, 313-320
for the rho_median ordering and the exact DFIRE-left/cpydock-right pairing).

Function set for this figure, in ascending rho_median order (identical
ordering convention to Figure S1 / main-text Figure 3A):
    fastDFIRE, vdW, PISA, MJ3h, SD, TOBI, DFIRE2

NOTE ON A SOURCE-LIST DISCREPANCY (documented, not silently fixed):
figs_research/spec_supplementary.md line 729 lists the seven source npz files
as {cpydock, dfire2, mj3h, pisa, sd, tobi, vdw} -- i.e. it includes cpydock and
omits fastdfire. That cannot be right for a "remaining seven" figure: cpydock is
one of the two functions the companion main-text design already shows (spec_fig3.md
Sec.5), and fastdfire/dfire are near-duplicate potentials that "coincide" in the
currently-rendered main Figure 3 panel c (see scripts/make_figures_v2.py, the
`hi = f in ("fastdfire","dfire","cpydock")` line) -- of that pair, only DFIRE is
the one spec_fig3.md's panel-C scatter design actually plots. Swapping cpydock
for fastdfire is also the only reading that yields exactly seven panels (9 total
scoring functions minus the 2 shown = 7; excluding all three of fastdfire/dfire/
cpydock would leave only 6). This script therefore plots
{fastdfire, vdw, pisa, mj3h, sd, tobi, dfire2}, not the literal spec_supplementary.md
line-729 list.

Design rules inherited from make_figures_v2.py / make_figureS1.py (Arial,
pdf.fonttype 42, LightDock scores are MAXIMISED so a working function gives a
NEGATIVE Spearman rho, 7EGR is HIGH leakage risk and never pooled with 9NX0).
All numbers below (rho, best-pose score, top-1-by-score pose) were read
directly from the *.npz files and cross-checked against
pod_mirror_t0/multiscore_panel.csv; see the printed verification block.

Source: pod_mirror_t0/multiscore/A2_7EGR_BLIND_rep1__{function}.npz
        (keys: rmsd, score, is_random_sample)
Run with /usr/bin/python3.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator

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
# figures_v2 is frozen; set FIGOUT to render a revised version elsewhere without touching it
OUT = os.environ.get("FIGOUT") or os.path.join(ROOT, "manuscript", "figures_v2")
MS_DIR = os.path.join(ROOT, "pod_mirror_t0", "multiscore")
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4

RUN = "A2_7EGR_BLIND_rep1"
BEST_RMSD = 5.990  # best pose actually present in this run's scored set

C_RANDOM = "#9a9a9a"
C_ENRICH = "#c9c9c9"
C_STAR = "#2a78d6"
C_DIAMOND = "#eb6834"

# ascending rho_median order, excluding dfire and cpydock (the two functions
# shown in the companion main-text design -- see module docstring)
FUNCS = ["fastdfire", "vdw", "pisa", "mj3h", "sd", "tobi", "dfire2"]
DISPLAY = {
    "fastdfire": "fastDFIRE", "vdw": "vdW", "pisa": "PISA", "mj3h": "MJ3h",
    "sd": "SD", "tobi": "TOBI", "dfire2": "DFIRE2",
}
# 'sd' has an extreme heavy tail (steric-clash penalty on badly overlapping
# poses, score_min ~ -2.28e9 vs score_max ~ +344) -- a linear axis would
# collapse >99% of the points into a hairline; give it a symlog y-axis instead
# and say so.
SYMLOG = {"sd"}


def manual_spearman(x, y):
    """Spearman rho with average-rank tie handling, no scipy dependency."""
    def rankdata_avg(a):
        order = np.argsort(a, kind="mergesort")
        ranks = np.empty(len(a), dtype=float)
        sorted_a = a[order]
        i, n = 0, len(a)
        while i < n:
            j = i
            while j + 1 < n and sorted_a[j + 1] == sorted_a[i]:
                j += 1
            ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return ranks
    rx, ry = rankdata_avg(x), rankdata_avg(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    return float(np.sum(rx * ry) / np.sqrt(np.sum(rx ** 2) * np.sum(ry ** 2)))


def main():
    data = {}
    for f in FUNCS:
        z = np.load(os.path.join(MS_DIR, f"{RUN}__{f}.npz"))
        rmsd, score, irs = z["rmsd"], z["score"], z["is_random_sample"].astype(bool)
        assert len(rmsd) == 2919, f"{f}: expected n=2919, got {len(rmsd)}"
        i_best = int(np.argmin(rmsd))
        i_top1 = int(np.argmax(score))
        assert abs(rmsd[i_best] - BEST_RMSD) < 0.001, \
            f"{f}: best-rmsd pose is {rmsd[i_best]:.3f}, not {BEST_RMSD}"
        rho_sub = manual_spearman(score[irs], rmsd[irs])
        rho_all = manual_spearman(score, rmsd)
        data[f] = dict(rmsd=rmsd, score=score, irs=irs,
                        best_rmsd=float(rmsd[i_best]), best_score=float(score[i_best]),
                        best_is_random=bool(irs[i_best]),
                        top1_rmsd=float(rmsd[i_top1]), top1_score=float(score[i_top1]),
                        n_random=int(irs.sum()), n_enrich=int((~irs).sum()),
                        rho_sub=rho_sub, rho_all=rho_all)
        print(f"{f:10s}  n_random={data[f]['n_random']}  n_enrich={data[f]['n_enrich']}  "
              f"best-pose score={data[f]['best_score']:.3f} (is_random_sample={data[f]['best_is_random']})  "
              f"top1-by-score: rmsd={data[f]['top1_rmsd']:.3f} score={data[f]['top1_score']:.3f}  "
              f"rho_sub={rho_sub:+.4f}  rho_all={rho_all:+.4f}")

    fig = plt.figure(figsize=(195 * MM, 108 * MM))
    gs = fig.add_gridspec(2, 4, hspace=0.55, wspace=0.42,
                           left=0.055, right=0.985, top=0.90, bottom=0.11)

    letters = "ABCDEFG"
    for k, f in enumerate(FUNCS):
        row, col = divmod(k, 4)
        ax = fig.add_subplot(gs[row, col])
        d = data[f]
        irs = d["irs"]

        ax.scatter(d["rmsd"][irs], d["score"][irs], s=3, c=C_RANDOM, alpha=0.35,
                   linewidths=0, zorder=2)
        ax.scatter(d["rmsd"][~irs], d["score"][~irs], s=3, facecolors="none",
                   edgecolors=C_ENRICH, alpha=0.25, linewidths=0.4, zorder=2)

        ax.axvline(BEST_RMSD, color=C_STAR, ls=(0, (3, 2)), lw=0.8, alpha=0.55, zorder=1)

        ax.scatter([d["best_rmsd"]], [d["best_score"]], marker="*", s=60,
                   c=C_STAR, edgecolors="black", linewidths=0.4, zorder=5)
        ax.scatter([d["top1_rmsd"]], [d["top1_score"]], marker="D", s=32,
                   c=C_DIAMOND, edgecolors="black", linewidths=0.4, zorder=5)

        if f in SYMLOG:
            ax.set_yscale("symlog", linthresh=1000)
            ax.set_ylim(d["score"].min() * 1.08, 3000)
            ax.set_yticks([-1e9, -1e7, -1e5, -1e3, 0, 1e3])
            ax.yaxis.set_minor_locator(NullLocator())
        else:
            ymin = min(d["score"].min(), d["best_score"], d["top1_score"])
            ymax = max(d["score"].max(), d["best_score"], d["top1_score"])
            pad = 0.12 * (ymax - ymin)
            ax.set_ylim(ymin - pad, ymax + pad)

        ax.set_xlim(5, 62)
        # Text box anchored bottom-right: the star is always at the left edge
        # (x = 5.990 A) and the diamond is always at the top of the panel (it
        # is the argmax score, and "up" is "better" throughout this figure),
        # so bottom-right is the one corner guaranteed clear of both special
        # markers regardless of which function's panel this is.
        ax.text(0.97, 0.04, r"$\rho$ = %+.3f" % d["rho_sub"] + "\n(random subsample,\nn = 2,000)",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=5.6,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#bbbbbb", lw=0.5))
        ax.set_title(DISPLAY[f] + (" (y: symlog)" if f in SYMLOG else ""),
                     fontsize=7.3, pad=3)
        ax.set_ylabel(DISPLAY[f] + " score\n(maximised)", fontsize=6.3, labelpad=2)
        ax.text(-0.19, 1.0, u"↑ better", transform=ax.transAxes, fontsize=5.3,
                color="#5A5A5A", ha="left", va="bottom")
        ax.tick_params(labelsize=6)
        if row == 1 or (row == 0 and col == 3 and False):
            ax.set_xlabel("ligand all-atom RMSD (Å)", fontsize=6.3)
        ax.text(-0.30, 1.14, letters[k], transform=ax.transAxes, fontsize=9,
                fontweight="bold", va="top", ha="left")

    # bottom axis label also needed under the 3 top-row panels that have no
    # panel directly below them at same column position -- here all 7 panels
    # sit in a 4-wide x 2-row grid with row0 fully populated (cols 0-3) and
    # row1 partially populated (cols 0-2, col 3 = legend), so row-0 panels in
    # col 3 has nothing below it and needs its own x label.
    fig.axes[3].set_xlabel("ligand all-atom RMSD (Å)", fontsize=6.3)

    # legend panel in the unused (row1, col3) slot -- hand-placed, not
    # plt.legend(), so title and entries cannot collide regardless of how
    # many lines each label wraps to.
    axl = fig.add_subplot(gs[1, 3])
    axl.axis("off")
    axl.set_xlim(0, 1)
    axl.set_ylim(0, 1)
    axl.text(0.5, 0.99, "run %s\nn = 2,919 poses" % RUN, ha="center", va="top",
              fontsize=6.5, fontweight="bold")

    legend_items = [
        ("o", C_RANDOM, "none", 4.5, "random subsample\n(n = 2,000)", 1),
        ("o", "none", C_ENRICH, 4.5,
         "enrichment-selected\n(top-500 by score +\ntop-500 by RMSD; n = 919)", 2),
        ("*", C_STAR, "black", 9, "best pose in set,\n%.3f Å" % BEST_RMSD, 1),
        ("D", C_DIAMOND, "black", 6,
         "rank-1 by score\n(RMSD varies\nby function)", 2),
    ]
    y = 0.80
    for marker, fc, ec, ms, label, extra_lines in legend_items:
        axl.plot([0.09], [y], marker=marker, ms=ms,
                  mfc=fc, mec=ec, mew=0.5, ls="", clip_on=False)
        axl.text(0.20, y + 0.025, label, ha="left", va="top", fontsize=5.8)
        y -= 0.155 + 0.058 * extra_lines

    axl.plot([0.03, 0.15], [y, y], color=C_STAR, ls=(0, (3, 2)), lw=0.9, alpha=0.6,
              clip_on=False)
    axl.text(0.20, y + 0.025,
              "common reference line\nat %.3f Å in every panel" % BEST_RMSD,
              ha="left", va="top", fontsize=5.8)

    caption = (
        "Figure S2. Score–RMSD funnels for the seven scoring functions not shown in "
        "Figure 3C (DFIRE, cpydock): fastDFIRE, vdW, PISA, MJ3h, SD, TOBI, DFIRE2, on the "
        "same n = 2,919 poses (see Methods §4.3) from run A2_7EGR_BLIND_rep1. Dashed line "
        "marks the best pose in the set, %.3f Å. Pure re-scoring of stored docking output; "
        "data-leakage status (Figure 4d) not applicable."
    ) % BEST_RMSD
    # The caption is drawn into the figure for standalone review copies. Set FIGNOCAP=1 for
    # the submission render, where the caption lives in the document and a second copy
    # baked into the image would be a duplicate.
    if not os.environ.get("FIGNOCAP"):
        fig.text(0.03, -0.07, caption, fontsize=6.1, va="top", ha="left", wrap=True,
                 transform=fig.transFigure)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS2." + e)
        fig.savefig(p, dpi=int(os.environ.get("FIGDPI", 400)), bbox_inches="tight")
    plt.close(fig)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS2." + e)
        print("wrote %-70s %8.1f KB" % (p, os.path.getsize(p) / 1024))

    return data


if __name__ == "__main__":
    main()
