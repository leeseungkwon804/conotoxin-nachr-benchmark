#!/usr/bin/env /usr/bin/python3
"""
Figure S3 -- Per-swarm-best distributions for each of the twelve runs separately.

Supplementary companion to main-text Figure 2B, which pools all 3,144 scored
swarms into six rung-level (BLIND/XDON/NATIVE) violins, one per anchor. This
figure keeps every one of the 12 runs (anchor x rung x seed rep) as its own
row so seed-to-seed reproducibility of the *distribution shape* can be read
directly -- the paper's central claim about the docking swarms is that the
distribution reproduces across the two GSO seeds while single-pose statistics
(e.g. which swarm is "best") do not (cf. Figure 3D, rank-1 instability).

Run with /usr/bin/python3 -- see scripts/make_figures_v2.py header for why
(Homebrew python 3.14's broken pyexpat cannot import matplotlib).

A1 (9NX0, MEDIUM leakage risk) and A2 (7EGR, HIGH leakage risk) are drawn as
two separate blocks (top six rows / bottom six rows) on two separate Axes
with two separate x-axis scales, and are never pooled into one statistic or
one shared scale: the per-swarm-best rigid-body floor differs by nearly a
factor of two between anchors (2.439 A vs 4.453 A) and the pooled maxima
differ by a factor of ~2.6 (131.127 A vs 50.284 A, both verified below from
the raw per-pose arrays) -- a shared axis would visually compress A2's true
spread relative to A1's, or vice versa.

Source: pod_mirror/s0_analysis/rmsd/<run>.npz (12 files), keys 'swarm' and
'rmsd' (all-atom ligand RMSD per pose, matching the 'best_allatom_A' /
'floor_allatom_A' columns of results/s0_ladder/s0_12run_main.csv). Per-swarm
minima are computed here by sorting by swarm id, splitting into groups, and
taking the min per group -- identical logic to per_swarm_best() in
scripts/make_figures_v2.py (this file does not import that module so it has
no hidden dependency on it, but the function body is copied verbatim and
verified below to reproduce the file's own precomputed 'per_swarm_best' key
and the summary columns in s0_12run_main.csv).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

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

C_BLIND, C_XDON, C_NATIVE = "#0072B2", "#E69F00", "#CC79A7"
C_GREY, C_FLOOR = "#5A5A5A", "#009E73"
RUNG_C = {"BLIND": C_BLIND, "XDON": C_XDON, "NATIVE": C_NATIVE}
RUNGS = ["BLIND", "XDON", "NATIVE"]
MM = 1 / 25.4

ANCHORS = [
    dict(anc="A1", pdb="9NX0", leak="MEDIUM leakage risk", floor=2.439, xlim=(1.9, 165)),
    dict(anc="A2", pdb="7EGR", leak="HIGH leakage risk", floor=4.453, xlim=(2.1, 68)),
]


def per_swarm_best(run):
    """Per-swarm minimum RMSD -- identical grouping logic to per_swarm_best()
    in scripts/make_figures_v2.py: sort by swarm id, split on swarm-id
    boundaries, take the min per group."""
    z = np.load(os.path.join(ROOT, "pod_mirror/s0_analysis/rmsd", run + ".npz"))
    sw, rm = z["swarm"], z["rmsd"]
    o = np.argsort(sw, kind="mergesort")
    sw, rm = sw[o], rm[o]
    return np.array([s.min() for s in np.split(rm, np.flatnonzero(np.diff(sw)) + 1)], float)


def raincloud_row(ax, data, y, color):
    """One row of a horizontal raincloud: half-violin ("cloud", upper half
    only) + thin boxplot (centre line) + jittered strip ("rain", lower half)."""
    parts = ax.violinplot([data], positions=[y], vert=False, widths=0.86,
                          showextrema=False)
    body = parts["bodies"][0]
    body.set_facecolor(color); body.set_edgecolor(color)
    body.set_alpha(0.40); body.set_linewidth(0.7); body.set_zorder(2)
    # Clip to the upper half of the row to turn the (symmetric) violin into a
    # half-violin "cloud". The clip rectangle's x-extent is expressed in axes
    # fraction (0-1), not data coordinates: on a log-scale x-axis, a data-coordinate
    # rectangle wide enough to guarantee full coverage (e.g. 1e-3 to 1e6) makes
    # matplotlib's tight-bbox machinery blow up the saved canvas to ~10x the
    # intended size, because it (incorrectly) folds the clip path's own extent
    # into the axes' tight bounding box. Axes-fraction x-coordinates are bounded
    # by construction and sidestep that entirely.
    clip_trans = blended_transform_factory(ax.transAxes, ax.transData)
    clip = Rectangle((0.0, y), 1.0, 0.44, transform=clip_trans)
    body.set_clip_path(clip)

    bp = ax.boxplot([data], positions=[y], vert=False, widths=0.13,
                    showfliers=False, patch_artist=True, zorder=5)
    for box in bp["boxes"]:
        box.set_facecolor("white"); box.set_edgecolor("black"); box.set_linewidth(0.7)
    for key in ("whiskers", "caps"):
        for ln in bp[key]:
            ln.set_color("black"); ln.set_linewidth(0.7)
    for ln in bp["medians"]:
        ln.set_color("black"); ln.set_linewidth(1.0)

    rng = np.random.default_rng(int(abs(y) * 1000) + len(data))
    jitter = y - 0.08 - rng.uniform(0, 0.30, size=len(data))
    ax.scatter(data, jitter, s=2.8, color=color, alpha=0.32, linewidths=0, zorder=1)


def block(ax, anc, pdb, leak, floor, xlim):
    rows = [(rung, rep) for rung in RUNGS for rep in (1, 2)]
    ys = list(range(len(rows) - 1, -1, -1))  # 5..0, BLIND rep1 at top
    stats = []
    for (rung, rep), y in zip(rows, ys):
        run = "%s_%s_%s_rep%d" % (anc, pdb, rung, rep)
        d = per_swarm_best(run)
        n, med, p5, mn = len(d), float(np.median(d)), float(np.percentile(d, 5)), float(d.min())
        stats.append((run, rung, rep, n, med, p5, mn, d))
        raincloud_row(ax, d, y, RUNG_C[rung])

    ax.axvline(3.0, color=C_GREY, ls=(0, (4, 3)), lw=1.0, zorder=6)
    ax.axvline(floor, color=C_FLOOR, lw=1.3, zorder=6)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.75, 5.85)
    ax.set_yticks(ys)
    ax.set_yticklabels(["%s rep %d" % (rung, rep) for (rung, rep), y in zip(rows, ys)],
                        fontsize=6.5)
    ax.set_xlabel("Per-swarm-best ligand RMSD, all-atom (A), log scale")
    ax.set_title("%s  %s -- %s" % (anc, pdb, leak), fontsize=7.8, pad=6,
                 fontweight="bold", loc="left")

    # landmark labels
    ytop = 5.85
    ax.text(3.0, ytop, " 3 A gate", color=C_GREY, fontsize=6, rotation=90,
            ha="left", va="top")
    ax.text(floor, ytop, " floor %.3f A" % floor, color=C_FLOOR, fontsize=6,
            rotation=90, ha="left", va="top")

    # per-row n / median / p5 / min, in the axes' blended (data-y, axes-x) frame
    trans = ax.get_yaxis_transform()
    for (run, rung, rep, n, med, p5, mn, d), y in zip(stats, ys):
        ax.text(1.015, y, "n=%d\nmed %.1f\np5 %.1f\nmin %.2f" % (n, med, p5, mn),
                transform=trans, fontsize=5.3, va="center", ha="left",
                linespacing=1.15, color="#222222")

    return stats


def main():
    fig = plt.figure(figsize=(180 * MM, 225 * MM))
    gs = fig.add_gridspec(2, 1, hspace=0.20, left=0.16, right=0.82, top=0.91, bottom=0.045)

    print("run                        rung    rep    n   median     p5      min")
    print("-" * 72)
    all_stats = {}
    for spec in ANCHORS:
        ax = fig.add_subplot(gs[ANCHORS.index(spec), 0])
        s = block(ax, spec["anc"], spec["pdb"], spec["leak"], spec["floor"], spec["xlim"])
        all_stats[spec["anc"]] = s
        for run, rung, rep, n, med, p5, mn, d in s:
            print("%-26s %-6s %2d  %4d  %7.3f  %6.3f  %6.3f" % (run, rung, rep, n, med, p5, mn))

    fig.legend(handles=[Line2D([], [], color=C_BLIND, lw=5, alpha=0.6, label="BLIND (no restraints)"),
                        Line2D([], [], color=C_XDON, lw=5, alpha=0.6, label="XDON (cross-donated)"),
                        Line2D([], [], color=C_NATIVE, lw=5, alpha=0.6, label="NATIVE (own contacts)"),
                        Line2D([], [], color=C_GREY, lw=1.2, ls=(0, (4, 3)), label="3 A gate"),
                        Line2D([], [], color=C_FLOOR, lw=1.4, label="anchor's rigid-body floor")],
              loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.49, 1.012),
              handletextpad=0.5, columnspacing=1.0)

    # ---- sanity checks against the file's own precomputed key and the CSV summary ----
    ok_key = True
    for anc_stats in all_stats.values():
        for run, rung, rep, n, med, p5, mn, d in anc_stats:
            z = np.load(os.path.join(ROOT, "pod_mirror/s0_analysis/rmsd", run + ".npz"))
            ok_key &= bool(np.allclose(np.sort(d), np.sort(z["per_swarm_best"])))
    print("\ncomputed per_swarm_best matches file's stored 'per_swarm_best' key for all 12 runs:", ok_key)

    main_csv = {r["run"]: r for r in csv.DictReader(
        open(os.path.join(ROOT, "results/s0_ladder/s0_12run_main.csv")))}
    ok_csv = True
    for anc_stats in all_stats.values():
        for run, rung, rep, n, med, p5, mn, d in anc_stats:
            row = main_csv[run]
            ok_csv &= abs(med - float(row["per_swarm_best_median_A"])) < 0.005
            ok_csv &= abs(p5 - float(row["per_swarm_best_p5_A"])) < 0.005
            ok_csv &= abs(mn - float(row["best_allatom_A"])) < 0.005
    print("computed median/p5/min matches results/s0_ladder/s0_12run_main.csv for all 12 runs:", ok_csv)
    assert ok_key and ok_csv, "Figure S3 data does not reproduce the checked-in summary -- stop."

    pooled_A1 = np.concatenate([d for _, _, _, _, _, _, _, d in all_stats["A1"]])
    pooled_A2 = np.concatenate([d for _, _, _, _, _, _, _, d in all_stats["A2"]])
    print("A1 pooled per-swarm-best: n=%d min=%.3f max=%.3f" % (len(pooled_A1), pooled_A1.min(), pooled_A1.max()))
    print("A2 pooled per-swarm-best: n=%d min=%.3f max=%.3f" % (len(pooled_A2), pooled_A2.min(), pooled_A2.max()))
    ns = [n for s in all_stats.values() for _, _, _, n, _, _, _, _ in s]
    print("n range across the 12 runs: %d-%d" % (min(ns), max(ns)))

    caption = (
        "Figure S3. Per-swarm-best ligand RMSD for each of the 12 runs individually "
        "(n = %d–%d swarms/run), rather than pooled by rung as in Figure 2b; rows are "
        "annotated with n, median, 5th percentile, and minimum (see Methods §4.6). A1 and "
        "A2 are plotted on separate, non-shared axes; pooled per-swarm-best RMSD reaches "
        "%.1f Å at A1 and %.1f Å at A2. The 3 Å gate and each anchor's rigid-body floor are "
        "marked on each panel."
    ) % (min(ns), max(ns), pooled_A1.max(), pooled_A2.max())
    fig.text(0.03, -0.025, caption, fontsize=6.1, va="top", ha="left", wrap=True,
             transform=fig.transFigure)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS3." + e)
        fig.savefig(p, dpi=400, bbox_inches="tight")
    plt.close(fig)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS3." + e)
        print("wrote %-46s %8.1f KB" % (p, os.path.getsize(p) / 1024))


if __name__ == "__main__":
    main()
