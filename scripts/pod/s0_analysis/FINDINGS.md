# S0 — findings that are about the *instrument*, not the result

Everything here was measured on this pod on 2026-07-30. Each item states how it was
established. These are separable from the ladder result and are the part of the work
that survives regardless of how the docking turns out.

---

## F1. The 3 Å gate was unreachable at anchor A2 by construction

S0 docks an **apo NMR conformer as a rigid body**. A rigid body cannot be superimposed
on a different conformer, so the optimally superposed (Kabsch) deviation between the
docked conformer and the native bound conformer is a **hard lower bound on any pose's
RMSD**. Measured, on the matched heavy atoms actually used in the analysis:

| anchor | toxin | docked apo conformer | floor all-atom | floor backbone | floor CA |
|---|---|---|---|---|---|
| A1 | ImI | 1IMI model 1 | **2.439 Å** | 1.056 Å | 1.073 Å |
| A2 | RgIA | 2JUT model 1 | **4.453 Å** | 2.397 Å | 2.117 Å |

**A2's floor exceeds the gate threshold.** And it is not a bad-model artefact — the
minimum over all 20 deposited 2JUT models is still **3.248 Å** all-atom (model 19):

| 2JUT ensemble (20 models) | min | median | max |
|---|---|---|---|
| all-atom | 3.248 (m19) | 4.337 | 4.852 |
| backbone | 1.298 (m18) | 2.276 | 2.518 |
| CA | 1.277 (m18) | 2.035 | 2.236 |

So **no available apo conformer of RgIA can reach 3 Å all-atom against 7EGR chain P**,
whatever the docking does. A negative gate result at A2 on the all-atom basis carries no
information about docking quality.

Worse for the metric: the two RgIA copies **inside 7EGR itself** (chains P and Q) differ
by **2.526 Å all-atom** (backbone 0.697 Å, CA 0.677 Å). The gate threshold is barely above
the copy-to-copy variation of the ground truth. All-atom ligand-RMSD is not a well-posed
metric for a peptide with four arginines — the deviation lives almost entirely in the long
side chains (worst atoms after superposition: NH2 13.78 Å, NH1 12.62 Å, OH 12.46 Å, CZ
11.84 Å).

At A1 the picture differs: the five ImI copies in 9NX0 are NCS-identical
(chain F vs chain I = **0.001 Å** on every basis), so there is no ground-truth noise, and
the floor of 2.439 Å leaves 0.56 Å of headroom under 3 Å — reachable, but only just.

**Consequence for reporting.** The gate must be quoted on a basis where it is reachable
(backbone or CA), with the floor stated, and the all-atom number reported alongside as
uninterpretable at A2. Every run is therefore scored on all three bases and the
`excess = best − floor` is reported, which is the part the docking is responsible for.

## F2. Apo-source choice alone can decide a 3 Å benchmark

Three apo NMR entries for the *same* toxin (ImI) are in the workspace. Their achievable
floors against 9NX0 chain I:

| entry | models | min all-atom | median | max | min backbone | min CA |
|---|---|---|---|---|---|---|
| **1G2G** (α7-recognition conformation, 2000-11-08) | 30 | **1.022** (m24) | 1.403 | 2.131 | 0.702 | 0.630 |
| 1IMI (used in S0) | 20 | 1.760 (m04) | 2.192 | 2.735 | 0.986 | 0.985 |
| 1CNL | 10 | 4.018 (m03) | 4.478 | 4.574 | 2.157 | 2.191 |

The spread across entries is **1.02–4.57 Å, wider than the 3 Å threshold itself.** Picking
1CNL would make the gate unpassable; picking 1G2G would make it comfortable. 1G2G model 1
is already 1.096 Å, so the effect is not model cherry-picking — it is entry choice.

A benchmark that docks "the apo structure" without specifying and justifying the entry and
model is reporting its own selection, not the tool's accuracy. **This applies to the
published conopeptide docking literature as much as to us.**

*Caveat on using this to rescue a result:* selecting the conformer by RMSD-to-native is
answer-selection and reintroduces the leakage the S0 ladder exists to remove. Any
conformer swap must be justified by a rule fixed in advance (e.g. "model 1 of the entry
determined as the α7-recognition conformation") and never by proximity to the answer.

## F3. `gso_<step>.out` pairs post-move coordinates with pre-move scores

LightDock 0.9.4's per-step output file writes each glowworm's position **after** that
step's move next to the objective value evaluated **before** it.

Mechanism, from source: `GSO.run()` executes `update_luciferin()` → `movement_phase()` →
`save()` per step, and `Glowworm.compute_luciferin()` re-evaluates the objective only
`if self.moved or self.step == 0`. So `scoring` is the objective at the position held at
the *start* of the step, while the saved coordinates are the position after the move.

Established four independent ways:

1. **Signature.** `gso_0.out` carries `Scoring = 0.00000000` for every glowworm and
   `Luciferin = 5.0` (the initial value) — `save(0)` runs before any evaluation.
2. **No-movement control.** A purpose-built 1-step run: at step 1 every glowworm has the
   same luciferin, so `is_neighbor` (which requires `self.luciferin < other.luciferin`)
   finds no neighbours and **0 of 200 glowworms move**. There the correctly-specified
   scoring function reproduces the file to **median |Δ| = 0.0083** — i.e. the machinery is
   exact when the position is unchanged.
3. **Consecutive-step test (decisive).** An 11-step run saves both step 10 and step 11,
   with 180 of 200 glowworms moving between them:

   | pairing tested | median &#124;Δ&#124; | max | within 0.05 |
   |---|---|---|---|
   | **E(P₁₀) vs `Scoring` in gso_11.out** (pre-move) | **0.0081** | 0.188 | **186/200** |
   | E(P₁₁) vs `Scoring` in gso_11.out (post-move) | 5.31 | 216.1 | 24/200 |
   | E(P₁₀) vs `Scoring` in gso_10.out (same-file) | 5.68 | 203.3 | 24/200 |

   The 24/200 residual matches the 20 non-movers plus coincidences.
4. **At production step count.** On the real 100-step runs, recomputing the correctly
   specified score at the written coordinates disagrees with the file for **0 of 1000**
   sampled poses, median |Δ| ≈ 0.42 (the mismatch shrinks as the swarm converges: ≈5.3 at
   step 11, ≈0.42 at step 100).

**Impact.** `lgd_rank_swarm.py` calls `write_ranking_to_file(..., order_by="scoring")`,
which sorts on this column. LightDock's own ranking tool therefore ranks each pose by a
neighbouring pose's energy. Anyone who selects "the top pose" from LightDock output
inherits this. It does **not** affect the search itself — the GSO uses correctly paired
values internally — only the post-hoc reading of the output.

## F4. Receptor restraints enter the score, so BLIND and restrained rungs have different
score ceilings

The vdW driver returns `score = energy × (1 + f_R + f_L)` with `energy = −vdw_energy`, and
`ScoringFunction.restraints_satisfied()` returns `0.0` for an empty restraint set.

Measured per run from the run's own `lightdock.info`, which dumps `vars(args)` at
simulation start:

| run | `restraints` | `receptor_restraints` | `ligand_restraints` |
|---|---|---|---|
| A2_7EGR_BLIND_rep1 | `None` | `None` | `None` |
| A2_7EGR_XDON_rep1 | `restraints.list` | `{'active': 12, 'passive': 0, 'blocked': 0}` | `{'active': 0, …}` |
| A2_7EGR_NATIVE_rep1 | `restraints.list` | `{'active': 33, …}` | `{'active': 0, …}` |

Confirmed numerically by recomputing the score at the written coordinates two ways —
median(recomputed / file):

| run | restraint-free | restraints wired in |
|---|---|---|
| A2_7EGR_BLIND_rep1 (R=0) | **1.0000** | 1.0000 (identical by construction) |
| A2_7EGR_XDON_rep1 (R=12) | 0.9257 | **1.0000** |
| A2_7EGR_NATIVE_rep1 (R=33) | 0.9694 | **1.0000** |

The BLIND row is the control: with no restraints the restraint-free recomputation matches,
so the systematic 7.4 % / 3.1 % shortfall in the restrained rows is the missing `(1 + f_R)`
factor (median f_R ≈ 0.080 and 0.032).

**So `f_L ≡ 0` everywhere — the brief's claim holds — but `f_R` is live in XDON/NATIVE and
identically 0 in BLIND.** The achievable multiplier is ×1 for BLIND versus up to ×2 for the
restrained rungs.

**Consequences.**
- Raw `Scoring` values must **not** be compared across rungs.
- This is the same defect class as AGENT.md §2 item 1 (the ×3-vs-×2 *ligand*-restraint
  asymmetry). The S0 design closed the L half — L = 0 verified in all 12 runs — but not
  the R half.

### Measured size of the effect

The same poses were re-scored with the restraints removed (an adapter built exactly as
`lgd_calculate_scoring.py` does), putting every rung on one scale. Per run, ≈2 800 poses:
the 500 best by published score, the 500 best by RMSD, and a seeded random 2 000. The
random block alone is used for percentiles and ρ, so those are unbiased; the estimator was
validated on A2 BLIND rep1, where it returns 7.05 / −0.073 against 7.64 / −0.080 from the
full ranking.

| run | R | published top-1 → RMSD | restraint-free top-1 → RMSD | shift | published rank of restraint-free top-1 |
|---|---|---|---|---|---|
| A2 BLIND rep1 | 0 | 25.251 Å | 25.251 Å | +0.000 | 1 |
| A2 BLIND rep2 | 0 | 25.249 Å | 25.544 Å | +0.294 | 2 |
| A2 XDON rep1 | 12 | 12.123 Å | 15.236 Å | +3.113 | 33 |
| A2 XDON rep2 | 12 | 11.817 Å | 11.811 Å | −0.006 | 11 |
| A2 NATIVE rep1 | 33 | 12.763 Å | **29.101 Å** | **+16.338** | 156 |
| A2 NATIVE rep2 | 33 | 12.960 Å | **29.267 Å** | **+16.307** | 68 |

Three things follow, and they are **anchor-specific** — see the next subsection:

1. **The BLIND control lands where it must** (+0.000, +0.294), which validates the procedure.
2. **At A2, NATIVE's published top-1 is inflated by ≈16.3 Å in both reps.** On the common
   scale A2 NATIVE's score-selected pose is 29.1 / 29.3 Å — the **worst of all six A2 runs,
   worse than BLIND's 25.3 / 25.5**. There the restraint term does manufacture apparent
   scoring competence, which makes an "the A2 NATIVE rung ranks poses well" claim circular:
   the bonus uses the answer. **This does not generalise — A1 NATIVE shows +0.219 Å.**
3. **A2 XDON's shift does not reproduce** (+3.113 vs −0.006) and must not be cited on its own.

*Limit of this test, stated plainly:* the re-scored objective is not the one the search
optimised, so these numbers mean "among the poses this run produced, which does pure vdW
prefer" — not "what would restraint-free docking have produced" (that is the BLIND run).
For a fair rung-to-rung comparison of pose *selection*, this is the right quantity.

### The reshuffle is real; its size is erratic and does NOT track restraint coverage

Restraint-set coverage was recomputed here rather than taken from the project notes, and
agrees with them exactly (native = residues within 5.0 Å of the native ligand):

| restraint set | R | native footprint | overlap | **recall of native** | precision |
|---|---|---|---|---|---|
| A1 XDON (9NX0 ← 2BYP) | 29 | 28 | 25 | **0.893** | 0.862 |
| A1 NATIVE | 28 | 28 | 28 | **1.000** | 1.000 |
| A2 XDON (7EGR ← 6HY7) | 12 | 33 | 12 | **0.364** | 1.000 |
| A2 NATIVE | 33 | 33 | 33 | **1.000** | 1.000 |

Full result, `swap` = where the restraint-free winner sits in the published ranking:

| run | R | recall | published top-1 | restraint-free top-1 | change | swap |
|---|---|---|---|---|---|---|
| A2 BLIND rep1 | 0 | 0.000 | 25.251 | 25.251 | +0.000 | 1 |
| A2 BLIND rep2 | 0 | 0.000 | 25.249 | 25.544 | +0.294 | 2 |
| A1 BLIND rep1 | 0 | 0.000 | 11.472 | 11.537 | +0.066 | 13 |
| A2 XDON rep2 | 12 | 0.364 | 11.817 | 11.811 | −0.006 | 11 |
| A2 XDON rep1 | 12 | 0.364 | 12.123 | 15.236 | +3.113 | 33 |
| **A1 NATIVE rep1** | 28 | **1.000** | 11.103 | 11.322 | **+0.219** | 10 |
| A1 XDON rep1 | 29 | 0.893 | 15.842 | 26.239 | **+10.397** | 287 |
| A2 NATIVE rep2 | 33 | 1.000 | 12.960 | 29.267 | **+16.307** | 68 |
| A2 NATIVE rep1 | 33 | 1.000 | 12.763 | 29.101 | **+16.338** | 156 |

An earlier draft of this file claimed a dose–response in recall, based on the first 8 runs.
**A1 NATIVE rep1 refutes it.** It has recall 1.000 — identical to A2 NATIVE, which shows
+16.3 Å — and changes by only **+0.219 Å**. A prediction of +10 to +16 Å was registered in
advance for this run and was wrong. The apparent monotonicity was an accident of which runs
had been scored.

What does hold:

1. **In every restrained run the restraint term changes which pose ranks first.** The
   restraint-free winner sits at published rank 10–287, never 1. The BLIND controls sit at
   1, 2 and 13, i.e. essentially unchanged, which is the required behaviour when f_R ≡ 0.
2. **The RMSD consequence of that reshuffle is erratic — +0.22 Å to +16.34 Å — and is not
   predictable from the restraint set.** It tracks how deep the swap goes (swaps inside the
   top ~13 cost ≤0.3 Å; swaps at 68–287 cost +10 to +16 Å) far better than it tracks
   coverage, but with nine points and A1 XDON out of order (swap 287 → +10.4, below
   swap 156 → +16.3) that is a description, not a law.
3. So the operative statement is a **robustness** one: *the published pose ranking is not
   stable under removing a non-physical term from the objective, and the error this induces
   cannot be bounded from the restraint set.* Any "top pose" taken from a restrained
   LightDock run inherits that.

### Sign of ρ is not uniform across runs

*(Updated to all 12 runs. An earlier draft of this paragraph read "7 of 8 negative, with
A1_9NX0_XDON_rep1 the sole exception"; with the remaining four runs scored that understates
it — the exception reproduces across seeds and is a property of A1's restrained rungs, not
of one run.)*

Across all 12 runs, Spearman ρ against RMSD is negative (the correct direction) in
**11/12** pose-level published, **10/12** swarm-level published and **9/12** pose-level
restraint-free. Every departure sits at **A1's restrained rungs**:

| run | pose-level published | swarm-level published | pose-level restraint-free |
|---|---|---|---|
| A1_9NX0_XDON_rep1 | −0.0027 (≈ 0) | **+0.1146** | **+0.0736** |
| A1_9NX0_XDON_rep2 | **+0.0062** | **+0.0547** | **+0.0832** |
| A1_9NX0_NATIVE_rep1 | −0.0585 | −0.1112 | **+0.0412** |
| A1_9NX0_NATIVE_rep2 | −0.0803 | −0.1577 | −0.0009 (≈ 0) |

All eight A2 runs and both A1 BLIND runs are negative on all three bases (|ρ| 0.065–0.460).
So the correct statement is not "one odd run": **score–RMSD discrimination collapses
specifically in A1's restrained rungs, and does so reproducibly across both GSO seeds.**
ρ is stable *within* a run, having thousands of pairs behind it, but is not uniformly signed
*across* runs, and must be reported per run rather than summarised as "always negative".

## F7. Every single-pose metric is unstable across GSO seeds; only distribution statistics reproduce

Each rung was run twice, identical in every respect except the GSO seed (`-sg` 324324 vs
20260730). That makes rep1-vs-rep2 a direct measurement of how much each candidate metric is
just noise. Measured over the four rungs that have both reps:

| metric | A1 BLIND | A2 BLIND | A2 XDON | A2 NATIVE | verdict |
|---|---|---|---|---|---|
| per-swarm-best **median** | 43.77 / 44.05 (Δ0.28) | 28.54 / 27.50 (Δ1.04) | 14.67 / 14.52 (Δ0.15) | 17.81 / 17.92 (Δ0.11) | **reproduces** |
| per-swarm-best **p5** | 19.20 / 18.67 (Δ0.54) | 9.84 / 10.54 (Δ0.70) | 10.28 / 10.56 (Δ0.28) | 10.18 / 11.27 (Δ1.08) | **reproduces** |
| `best_rmsd` | 8.74 / 7.65 (Δ1.10) | 5.99 / 7.17 (Δ1.18) | 9.35 / 7.81 (Δ1.54) | 8.18 / 8.86 (Δ0.68) | usable with ±1.5 Å |
| **top-1-by-score → RMSD** | 11.47 / **87.99** (Δ**76.5**) | 25.25 / 25.25 (Δ0.00) | 12.12 / 11.82 (Δ0.31) | 12.76 / 12.96 (Δ0.20) | **unusable** |
| percentile of best pose | 91.04 / 96.69 | 7.64 / 8.03 | 6.78 / 3.03 | **88.10 / 0.85** | **unusable** |

Two metrics fail, and both fail for the same structural reason: they are defined by a
**single pose** out of 26 400–78 000, so nothing averages the noise out.

The A1 BLIND blow-up is worth understanding rather than dismissing. The centroid displacement
of the top-scoring pose is **4.2 Å in rep1 and 87.6 Å in rep2** — in rep2 the best-vdW pose is
on the far side of the receptor. That is possible because A1's docked receptor is the
full-length subunit pair including the transmembrane domain: the maximum pose-centroid
displacement reached is **≈140 Å for A1** versus **≈59 Å for A2 BLIND** and ≈31–37 Å for the
A2 restrained rungs. With a receptor of that extent the highest-vdW-score pose can sit almost
anywhere on it, and which "anywhere" wins is decided by the search seed. This is the same
receptor-scope asymmetry noted earlier, now showing up as metric instability.

**Recommendation.** Report the **per-swarm-best distribution** (median, p5, counts under
thresholds) plus rank tests, and `best_rmsd` with an explicit ±1.5 Å reproducibility band.
Treat every score-selected single-pose statistic as uninterpretable. This is not a new
prescription — AGENT.md §2 already states that the correct statistic is the per-swarm-best
distribution with rank tests rather than raw maxima; these measurements independently confirm
it, and the analysis drifted from it before this check was run.

Importantly, **the ladder verdict does not depend on the unstable metrics**: question (b) rests
on `best_rmsd` and the per-swarm-best median, both of which reproduce.

### The percentile statistic, in detail

A natural-looking statistic — where the best-RMSD pose ranks by score — turns out to be
noise. Within the **same rung**, changing only the GSO seed:

| run | published percentile of best-RMSD pose |
|---|---|
| A2 NATIVE rep1 | **88.10** |
| A2 NATIVE rep2 | **0.85** |

Both runs have identical receptor, ligand, restraints, swarm placement and swarm count; they
differ only in `-sg`. The statistic swings from the 88th to the 0.85th percentile. The reason
is structural: there is essentially one near-native pose among ~26 400, so its score rank is
a **single order statistic** with no averaging behind it.

This retrospectively voids two readings that looked like results earlier in this analysis —
"the near-native pose sits in the bottom 8 % by score" (A2 BLIND rep1, 7.64) and "in some
runs it sits in the top 10 %" (88.10 / 91.04). Neither is evidence about the scoring
function.

**Use instead:** the Spearman ρ values reported **per run** — they have thousands of pairs
behind them and are stable within a run, but their sign is not uniform across runs (see F4,
last subsection: 11/12, 10/12 and 9/12 negative on the three bases, every departure at A1's
restrained rungs).

⚠ *An earlier draft of this paragraph recommended `top-1-by-score → RMSD` as the stable
replacement, on the strength of the A2 runs alone (A2 BLIND 25.25 / 25.25, XDON 12.12 / 11.82,
NATIVE 12.76 / 12.96). **That recommendation is withdrawn, and is contradicted by the table
above:** A1 BLIND gives 11.47 / 87.99 across nothing but a seed change (Δ 76.52 Å).
`top-1-by-score → RMSD` is a single-pose statistic and must not be quoted.*

## F5. Swarm placement is a bias, not a bound

Nearest swarm centre to the native ligand centroid, with `swarm_radius = 10.0 Å`:

| runs | centres | nearest centre → native centroid | centres within 10 Å of a native atom | centroid inside a swarm sphere |
|---|---|---|---|---|
| A1 BLIND | 390 | 13.88 Å | 4 | no |
| A1 XDON | 390 | 14.61 Å | 8 | no |
| A1 NATIVE | 392 | 13.62 Å | 12 | no |
| A2 BLIND | 134 | **6.08 Å** | 15 | **yes** (5) |
| A2 XDON | 134 | 10.67 Å | 57 | no |
| A2 NATIVE | 132 | 8.81 Å | 43 | **yes** (1) |

Read alone this would say A1 cannot reach the site. It does not hold: **82–98 % of
glowworms end outside their own swarm sphere**, so initial placement biases the search
without bounding it. What matters is the achieved displacement, and that is reported per run.

*(An earlier draft of this paragraph gave the range as "97–98 % / 0.970–0.982", measured on
the BLIND runs only. Over all 12 completed runs the escape fraction is 0.819–0.982, and the
spread is systematic rather than noise — restraints hold glowworms nearer their spheres:)*

| rung | escape fraction | |
|---|---|---|
| BLIND | 0.969–0.982 | least confined |
| NATIVE | 0.895–0.955 | |
| XDON | 0.819–0.970 | most confined (A2 XDON 0.819 / 0.822) |

The conclusion is unaffected — a 82 % escape rate is still a large majority leaving the
sphere — but the correct statement is "82–98 %", and the restrained rungs are the low end.

Note the counter-intuitive entry: **A2 BLIND's nearest swarm centre (6.08 Å) is closer to
the native site than A2 XDON's (10.67 Å)**. Restrained placement concentrates swarms on
the *restrained patch*, and the cross-donated patch covers only part of the true footprint
(recall 12/33), so it pulls sampling to a partly-wrong location while uniform blind
placement happens to land a centre on the site. This is the mechanism behind the ladder
inversion.

## F6. Integrity checks that passed

- **Native reference is the correct ligand copy.** Both anchors are homopentamers with
  several bound toxin copies; the docked receptor is one subunit pair. Verified by contact
  counting rather than filename: 9NX0 chain I contacts chain D (328 atom pairs, 16
  residues) and chain E (121, 12), and has **zero** contacts to A/B/C (20.6–30.4 Å away).
  7EGR chain P contacts chain G (203, 16) and chain H (305, 17), zero to F/I/J
  (17.1–32.2 Å). The 5 Å contact-residue totals — **28** for 9NX0 and **33** for 7EGR —
  independently reproduce the NATIVE restraint-set sizes.
- **Rungs differ only in restraints.** Per anchor, `receptor.pdb` and
  `lightdock_ligand.pdb` are md5-identical across all 6 runs; `restraints.list` md5 is
  identical between rep1 and rep2 of a rung; `use_anm` false, `glowworms` 200,
  `surface_density` 50.0, `swarm_radius` 10.0, `starting_points_seed` 324324 in all 12.
  Reps differ only in the GSO seed (`-sg` 324324 vs 20260730).
- **Pose reconstruction is exact.** The quaternion convention was determined, not assumed:
  of four candidates only `wxyz` reproduces `lgd_generate_conformations.py`
  (max per-atom deviation 0.000808 Å at A1, 0.000814 Å at A2; the others are off by
  18–25 Å). An independent path through real PDB files agrees to 6.8e-5 Å (A1) and
  7.7e-5 Å (A2) over 200 poses each.
- **Atom correspondence is complete.** A1 90/90 matched (the only unmatched native atom is
  the NH2 cap nitrogen of residue 13, stripped from the docked ligand); A2 105/105, all
  atoms on both sides, numbering offset +400. Paired by Needleman–Wunsch on residue names,
  not residue-number arithmetic.

---

## What this means for the write-up

The benchmark-audit contribution is stronger than the docking result. F1 and F2 say the
*threshold* used across this literature is not well posed for flexible peptides and can be
decided by an undocumented choice of apo entry. F3 and F4 are defects in the tool's output
and scoring that affect anyone using LightDock with restraints. None of these required the
docking to succeed.
