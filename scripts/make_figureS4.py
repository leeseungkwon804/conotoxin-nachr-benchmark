#!/usr/bin/env /usr/bin/python3
"""
Figure S4 -- Conformational floor of every deposited apo NMR model of the two
toxins (ImI at anchor A1/9NX0, RgIA at anchor A2/7EGR).

Companion figure to Table S6 (manuscript/supplementary/Table_S06_...): that
table gives the full 82-row per-model listing; this figure shows the same
numbers as individual points so the *shape* of each entry's distribution is
visible, not just its range -- in particular that all 10 deposited 1CNL
models sit above the 3 A gate and (separately, panel B) all 20 deposited
2JUT models do too, which the headline range alone (1.02-4.57 A) does not
convey.

Run with /usr/bin/python3 -- the Homebrew python 3.14 on this machine has a
broken pyexpat and cannot import matplotlib (see make_figures_v2.py header).

Design rules inherited from make_figures_v2.py / make_figureS3.py (read there
for full rationale, restated here only where this figure adds to them):
  * Arial, pdf.fonttype 42.
  * 7EGR (A2, HIGH leakage risk) and 9NX0 (A1, MEDIUM leakage risk) are two
    separate panels with two separate y-ranges of the SAME metric (all-atom
    Kabsch RMSD, A) -- never pooled into one distribution, both panels
    labelled with anchor + leakage risk in the title.
  * C_A1 (blue, #0072B2) and C_A2 (vermillion, #D55E00) are this project's
    established anchor colours (make_figures_v2.py C_A1/C_A2) -- panel A's
    three ImI entries are drawn as a light-to-dark ORDINAL ramp of the A1
    blue hue (entries have no inherent categorical identity beyond "which
    PDB deposition"; ramping preserves "blue = A1" while still separating
    1G2G/1IMI/1CNL by lightness), panel B's single RgIA entry (2JUT) is
    solid C_A2 vermillion.
  * C_FLOOR (green, #009E73) marks the one conformer S0 actually docked from
    in each ensemble (1IMI m01, 2JUT m01) -- identical colour role to the
    "anchor's rigid-body floor" line in Figure S3, since these two points
    ARE that floor value.
  * C_GREY (#5A5A5A) is reused for the 3 A gate line (as in Figure S3) and,
    as a diamond marker (shape carries the identity, not colour alone), for
    the bound-copy reference control in each panel -- this is not a
    deposited apo model, so it is deliberately drawn as a different marker
    shape in a desaturated colour rather than a fourth competing hue.

Source data (read directly, not from memory):
  - pod_mirror_t0/ensemble_floor_full.json, key "9NX0_ligand_I.pdb"
    (61 records: ImI_1G2G_m01-30, ImI_1IMI_m01-20, ImI_1CNL_m01-10,
    ImI_bound_9NX0_F)
  - results/s0_ladder/ensemble_floor.json, key "7EGR_ligand_P.pdb"
    (21 records: RgIA_2JUT_m01-20, RgIA_bound_7EGR_Q)
  - results/s0_ladder/intrinsic_floor.json (independent recomputation of the
    two S0-docked floors, used here only as a cross-check)
  - manuscript/supplementary/Table_S06_conformational_floors_all_apo_models.csv
    (the checked-in companion table for this figure; used only as a second
    independent cross-check, not as the primary source)
"""
import csv
import json
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
# figures_v2 is frozen; set FIGOUT to render a revised version elsewhere without touching it
OUT = os.environ.get("FIGOUT") or os.path.join(ROOT, "manuscript", "figures_v2")
os.makedirs(OUT, exist_ok=True)
MM = 1 / 25.4

C_A1, C_A2 = "#0072B2", "#D55E00"
C_GREY, C_FLOOR = "#5A5A5A", "#009E73"
# ordinal blue ramp for the three ImI entries, light->dark (validated with
# scripts/../dataviz skill's validate_palette.js --ordinal: monotone L, all
# adjacent steps >=0.06 OKLCH L apart, light end clears 3:1 vs white surface)
C_1G2G, C_1IMI, C_1CNL = "#5A9BC7", C_A1, "#00374F"


def swarm_offsets(y, bin_width, step):
    """Deterministic histogram-based beeswarm: bin y-values, spread points
    within a bin symmetrically about the category centre. Good enough for
    n<=30 per category; avoids a random-jitter dependency on RNG seed."""
    y = np.asarray(y, float)
    bin_id = np.floor(y / bin_width).astype(int)
    counts = {}
    for b in bin_id:
        counts[b] = counts.get(b, 0) + 1
    seen = {}
    x = np.zeros(len(y))
    for i, b in enumerate(bin_id):
        k = seen.get(b, 0)
        seen[b] = k + 1
        n = counts[b]
        x[i] = (k - (n - 1) / 2.0) * step
    return x


def load_imi():
    d = json.load(open(os.path.join(ROOT, "pod_mirror_t0/ensemble_floor_full.json")))
    recs = d["9NX0_ligand_I.pdb"]
    assert len(recs) == 61, "expected 61 ImI records (60 apo + 1 bound), got %d" % len(recs)
    entries = {"1G2G": [], "1IMI": [], "1CNL": []}
    bound = None
    for r in recs:
        f = r["file"]
        if f.startswith("ImI_bound"):
            bound = r
        else:
            entry = f.split("_")[1]
            entries[entry].append(r)
    assert [len(entries[e]) for e in ("1G2G", "1IMI", "1CNL")] == [30, 20, 10], \
        "ImI entry counts do not match spec (1G2G 30 / 1IMI 20 / 1CNL 10)"
    assert bound is not None and abs(bound["allatom"] - 0.001) < 1e-9, \
        "ImI bound control should be 0.001 A, got %r" % (bound,)
    m01 = next(r for r in entries["1IMI"] if r["file"] == "ImI_1IMI_m01.pdb")
    assert abs(m01["allatom"] - 2.439) < 1e-9, "1IMI m01 floor should be 2.439 A, got %.3f" % m01["allatom"]
    return entries, bound, m01


def load_rgia():
    d = json.load(open(os.path.join(ROOT, "results/s0_ladder/ensemble_floor.json")))
    recs = d["7EGR_ligand_P.pdb"]
    assert len(recs) == 21, "expected 21 RgIA records (20 apo + 1 bound), got %d" % len(recs)
    apo = [r for r in recs if r["file"].startswith("RgIA_2JUT")]
    bound = next(r for r in recs if r["file"].startswith("RgIA_bound"))
    assert len(apo) == 20, "expected 20 2JUT models, got %d" % len(apo)
    assert abs(bound["allatom"] - 2.526) < 1e-9 and abs(bound["backbone"] - 0.697) < 1e-9, \
        "RgIA bound control should be 2.526/0.697 A, got %r" % (bound,)
    m01 = next(r for r in apo if r["file"] == "RgIA_2JUT_m01.pdb")
    assert abs(m01["allatom"] - 4.453) < 1e-9, "2JUT m01 floor should be 4.453 A, got %.3f" % m01["allatom"]
    return apo, bound, m01


def cross_check(entries_imi, bound_imi, apo_rgia, bound_rgia):
    """Independent cross-checks against intrinsic_floor.json (recomputed by a
    different script) and Table_S06 (built by a different agent from the same
    two JSON sources) -- both must agree before this figure is trusted."""
    intrinsic = json.load(open(os.path.join(ROOT, "results/s0_ladder/intrinsic_floor.json")))
    by_run = {r["run"]: r for r in intrinsic}
    a1 = by_run["A1_9NX0_BLIND_rep1"]
    a2 = by_run["A2_7EGR_BLIND_rep1"]
    assert abs(a1["kabsch_rmsd_allatom_A"] - 2.439) < 1e-9
    assert abs(a2["kabsch_rmsd_allatom_A"] - 4.453) < 1e-9
    print("cross-check vs intrinsic_floor.json: A1 2.439, A2 4.453 -- OK")

    csv_path = os.path.join(ROOT, "manuscript/supplementary/Table_S06_conformational_floors_all_apo_models.csv")
    if os.path.exists(csv_path):
        rows = list(csv.DictReader(open(csv_path)))
        assert len(rows) == 82, "Table S6 CSV should have 82 rows, has %d" % len(rows)
        n_imi_apo = sum(1 for r in rows if r["toxin"] == "ImI" and r["is_bound_control"] == "FALSE")
        n_imi_bound = sum(1 for r in rows if r["toxin"] == "ImI" and r["is_bound_control"] == "TRUE")
        n_rgia_apo = sum(1 for r in rows if r["toxin"] == "RgIA" and r["is_bound_control"] == "FALSE")
        n_rgia_bound = sum(1 for r in rows if r["toxin"] == "RgIA" and r["is_bound_control"] == "TRUE")
        assert (n_imi_apo, n_imi_bound, n_rgia_apo, n_rgia_bound) == (60, 1, 20, 1), \
            "Table S6 row composition mismatch: %r" % ((n_imi_apo, n_imi_bound, n_rgia_apo, n_rgia_bound),)
        docked = [r for r in rows if r["s0_docked_model"] == "TRUE"]
        docked_vals = sorted(round(float(r["allatom_rmsd_A"]), 3) for r in docked)
        assert docked_vals == [2.439, 4.453], "Table S6 s0_docked_model rows should be 2.439/4.453, got %r" % docked_vals
        print("cross-check vs Table_S06 CSV (82 rows, entry counts, docked-model floors): OK")
    else:
        print("NOTE: Table_S06 CSV not found at build time -- skipping that cross-check "
              "(figure still trusted on ensemble_floor_full.json / ensemble_floor.json / "
              "intrinsic_floor.json agreement above).")


def rank_of(apo_rgia, m01):
    """1IMI m01 and 2JUT m01 are docked-from conformers, not necessarily the
    ensemble's best floor -- report 2JUT m01's rank among its own 20 models,
    used in the panel-B annotation (matches Table S6 footnote c)."""
    vals = sorted(r["allatom"] for r in apo_rgia)
    return vals.index(m01["allatom"]) + 1, len(vals)


def draw_panel(ax, groups, gate=3.0, ylim=(-0.32, 5.35)):
    """groups: list of dicts with keys x, label, color, vals (list of allatom
    RMSD), n, highlight (dict with 'file','allatom' or None), bound (record
    or None), bin_width, step."""
    xt, xtl = [], []
    for g in groups:
        xt.append(g["x"])
        xtl.append(g["label"])
        if g["vals"] is not None and len(g["vals"]) > 0:
            vals = np.array(g["vals"])
            xs = g["x"] + swarm_offsets(vals, g["bin_width"], g["step"])
            is_hl = g.get("highlight_mask")
            base = ~is_hl if is_hl is not None else np.ones(len(vals), bool)
            ax.scatter(xs[base], vals[base], s=10, color=g["color"], alpha=0.78,
                       edgecolors="white", linewidths=0.3, zorder=3)
            med = float(np.median(vals))
            ax.plot([g["x"] - 0.30, g["x"] + 0.30], [med, med], color="black",
                     lw=1.1, zorder=4, solid_capstyle="butt")
            if is_hl is not None and is_hl.any():
                hx, hy = xs[is_hl][0], vals[is_hl][0]
                ax.scatter([hx], [hy], marker="*", s=170, facecolor=C_FLOOR,
                           edgecolor="black", linewidths=0.7, zorder=6)
                ax.annotate(g["highlight_label"], xy=(hx, hy),
                            xytext=g["highlight_xytext"], textcoords="data",
                            fontsize=5.6, ha=g.get("highlight_ha", "left"),
                            va="center",
                            arrowprops=dict(arrowstyle="-", color="black", lw=0.6,
                                             shrinkA=2, shrinkB=3))
        if g.get("bound") is not None:
            bx, by = g["x"], g["bound"]["allatom"]
            ax.scatter([bx], [by], marker="D", s=26, facecolor=C_GREY,
                       edgecolor="black", linewidths=0.5, zorder=5)
            va = g.get("bound_va", "bottom")
            ax.annotate(g["bound_label"], xy=(bx, by), xytext=(bx, g["bound_ytext"]),
                        textcoords="data", fontsize=5.6, ha="center", va=va,
                        linespacing=1.35,
                        arrowprops=dict(arrowstyle="-", color=C_GREY, lw=0.6,
                                         shrinkA=2, shrinkB=3))

    # separator between the deposited apo entries and the bound-copy control --
    # the control is not a deposited apo model and must read as visually distinct
    if len(groups) >= 2 and groups[-1].get("bound") is not None:
        sep_x = (groups[-2]["x"] + groups[-1]["x"]) / 2.0
        ax.axvline(sep_x, color="#c9c9c9", lw=0.7, ls=":", zorder=1)

    ax.axhline(gate, color=C_GREY, ls=(0, (4, 3)), lw=2.3, zorder=2)
    ax.text(xt[0] - 0.52, gate + 0.10, "3 Å gate", color=C_GREY, fontsize=6.5, ha="left")
    ax.set_xticks(xt)
    ax.set_xticklabels(xtl, fontsize=6.6)
    ax.set_xlim(xt[0] - 0.62, xt[-1] + 0.62)
    ax.set_ylim(*ylim)
    return xt


def main():
    entries_imi, bound_imi, m01_imi = load_imi()
    apo_rgia, bound_rgia, m01_rgia = load_rgia()
    cross_check(entries_imi, bound_imi, apo_rgia, bound_rgia)
    rank, n_rgia = rank_of(apo_rgia, m01_rgia)
    print("2JUT m01 (S0-docked) all-atom floor rank: %d of %d (worse = higher rank)" % (rank, n_rgia))

    all_1cnl_above_gate = all(r["allatom"] > 3.0 for r in entries_imi["1CNL"])
    all_2jut_above_gate = all(r["allatom"] > 3.0 for r in apo_rgia)
    all_1g2g_below_gate = all(r["allatom"] < 3.0 for r in entries_imi["1G2G"])
    all_1imi_below_gate = all(r["allatom"] < 3.0 for r in entries_imi["1IMI"])
    print("every 1CNL model above 3 A gate:", all_1cnl_above_gate)
    print("every 2JUT model above 3 A gate:", all_2jut_above_gate)
    print("every 1G2G model below 3 A gate:", all_1g2g_below_gate)
    print("every 1IMI model below 3 A gate:", all_1imi_below_gate)
    assert all_1cnl_above_gate and all_2jut_above_gate and all_1g2g_below_gate and all_1imi_below_gate

    fig = plt.figure(figsize=(180 * MM, 100 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.05, 1.55], wspace=0.10,
                          left=0.075, right=0.985, top=0.86, bottom=0.225)

    # ---------------- panel a: ImI / 9NX0 (A1, MEDIUM leakage risk) ----------------
    axA = fig.add_subplot(gs[0, 0])
    vals_1imi = np.array([r["allatom"] for r in entries_imi["1IMI"]])
    hl_mask = np.array([r["file"] == "ImI_1IMI_m01.pdb" for r in entries_imi["1IMI"]])
    groupsA = [
        dict(x=0, label="1G2G\n(n=30)", color=C_1G2G,
             vals=[r["allatom"] for r in entries_imi["1G2G"]],
             bin_width=0.05, step=0.075),
        dict(x=1, label="1IMI\n(n=20)", color=C_1IMI, vals=vals_1imi,
             bin_width=0.05, step=0.075, highlight_mask=hl_mask,
             highlight_label="S0 docked (m01)\nfloor 2.439 Å",
             highlight_xytext=(1.42, 3.55), highlight_ha="left"),
        dict(x=2, label="1CNL\n(n=10)", color=C_1CNL,
             vals=[r["allatom"] for r in entries_imi["1CNL"]],
             bin_width=0.05, step=0.075),
        dict(x=3.15, label="bound copy\n(n=1)", color=C_GREY, vals=None,
             bound=bound_imi, bound_label="chain F vs I (NCS)\n0.001 Å",
             bound_ytext=1.05),
    ]
    draw_panel(axA, groupsA)
    axA.set_ylabel("All-atom Kabsch RMSD to native bound toxin (Å)\n[rigid-body floor -- no docking]")
    axA.set_title("A   ImI — 9NX0 anchor (A1, MEDIUM leakage risk)",
                  fontsize=7.8, loc="left", fontweight="bold", pad=8)

    # ---------------- panel b: RgIA / 7EGR (A2, HIGH leakage risk) ----------------
    axB = fig.add_subplot(gs[0, 1])
    vals_2jut = np.array([r["allatom"] for r in apo_rgia])
    hl_mask_b = np.array([r["file"] == "RgIA_2JUT_m01.pdb" for r in apo_rgia])
    groupsB = [
        dict(x=0, label="2JUT\n(n=20)", color=C_A2, vals=vals_2jut,
             bin_width=0.07, step=0.075, highlight_mask=hl_mask_b,
             highlight_label="S0 docked (m01)\nfloor 4.453 Å\n(rank %d/%d)" % (rank, n_rgia),
             highlight_xytext=(0.55, 4.0), highlight_ha="left"),
        dict(x=1.15, label="bound copy\n(n=1)", color=C_GREY, vals=None,
             bound=bound_rgia,
             bound_label="chain Q vs P (2nd copy)\n2.526 Å (bb 0.697 Å)",
             bound_ytext=1.55, bound_va="top"),
    ]
    draw_panel(axB, groupsB)
    axB.set_yticklabels([])
    axB.set_title("B   RgIA — 7EGR anchor (A2, HIGH leakage risk)",
                  fontsize=7.8, loc="left", fontweight="bold", pad=8)

    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", color=C_GREY, label="deposited apo NMR model", ms=5),
        Line2D([], [], marker="*", ls="", markerfacecolor=C_FLOOR, markeredgecolor="black",
               label="conformer S0 actually docked from (= anchor floor, cf. Figure S3)", ms=11),
        Line2D([], [], marker="D", ls="", markerfacecolor=C_GREY, markeredgecolor="black",
               label="bound-copy reference (2nd crystal/cryo-EM copy, not a deposited apo model)", ms=6),
        Line2D([], [], color="black", lw=1.1, label="per-entry median"),
        Line2D([], [], color=C_GREY, lw=2.3, ls=(0, (4, 3)), label="3 Å docking gate"),
    ], loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.0),
        handletextpad=0.6, columnspacing=1.6, fontsize=6.2)

    caption = (
        "Figure S4. Apo-ensemble conformational floor of the two docked toxins, all-atom "
        "Kabsch RMSD to native (see Methods §4.1). (a) ImI-9NX0: PDB 1G2G (n=30), 1IMI "
        "(n=20), 1CNL (n=10); S0-docked 1IMI m01 floor %.3f Å; bound-copy control (chain F "
        "vs I) %.3f Å. (b) RgIA-7EGR: PDB 2JUT (n=20); S0-docked m01 floor %.3f Å (rank %d "
        "of %d); bound-copy control (chain Q vs P) %.3f Å all-atom / %.3f Å backbone."
    ) % (m01_imi["allatom"], bound_imi["allatom"], m01_rgia["allatom"], rank, n_rgia,
         bound_rgia["allatom"], bound_rgia["backbone"])
    # The caption is drawn into the figure for standalone review copies. Set FIGNOCAP=1 for
    # the submission render, where the caption lives in the document and a second copy
    # baked into the image would be a duplicate.
    if not os.environ.get("FIGNOCAP"):
        fig.text(0.02, -0.16, caption, fontsize=6.1, va="top", ha="left", wrap=True,
                 transform=fig.transFigure)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS4." + e)
        fig.savefig(p, dpi=int(os.environ.get("FIGDPI", 400)), bbox_inches="tight")
    plt.close(fig)

    for e in ("png", "pdf"):
        p = os.path.join(OUT, "FigureS4." + e)
        print("wrote %-40s %8.1f KB" % (p, os.path.getsize(p) / 1024))


if __name__ == "__main__":
    main()
