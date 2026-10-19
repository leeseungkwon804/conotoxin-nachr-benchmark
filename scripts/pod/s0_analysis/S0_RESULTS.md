# S0 3-rung ladder — results

Runs scored: **12/12 complete**

Anchors, from live RCSB queries on 2026-07-30:

| | A1 | A2 |
|---|---|---|
| PDB | **9NX0** | **7EGR** |
| complex | human α7 nAChR (CHRNA7) + α-conotoxin ImI | *Aplysia californica* AChBP + RgIA |
| method / resolution | cryo-EM 3.06 Å | X-ray 2.5 Å |
| released | 2026-02-04 | 2022-04-06 |
| docked receptor | chains D,E — 392 res each (1–478, **ECD + TMD**), 6326 heavy atoms | chains H,G — 206 res each (19–224, **ECD only**), 3276 heavy atoms |
| docked ligand (apo, rigid) | ImI, 1IMI model 1, 90 atoms | RgIA, 2JUT model 1, 105 atoms |
| native reference | 9NX0 chain I | 7EGR chain P |
| matched atoms | 90 / 90 | 105 / 105 |

Both anchors are subunit **pairs cut from a homopentamer**, so both present a complete (+)/(−) interface. A2 is *not* a human receptor.


## (a) Kill gate — does any run reach < 3 Å best ligand-RMSD?

**All-atom: NO. Nothing under 3 Å.**

**Backbone: NO. Nothing under 3 Å.**

But the threshold has to be read against the **intrinsic floor** — the docked ligand is an apo conformer treated as rigid, and the optimally superposed (Kabsch) deviation from the bound conformer is a hard lower bound on any pose:

| run | floor all-atom | best all-atom | **excess** | floor backbone | best backbone | **excess** | floor CA | best CA |
|---|---|---|---|---|---|---|---|---|
| A1_9NX0_BLIND_rep1 | 2.44 | 8.74 | **6.31** | 1.06 | 8.20 | **7.14** | 1.07 | 8.19 |
| A1_9NX0_BLIND_rep2 | 2.44 | 7.65 | **5.21** | 1.06 | 6.32 | **5.27** | 1.07 | 6.32 |
| A1_9NX0_XDON_rep1 | 2.44 | 9.15 | **6.71** | 1.06 | 8.00 | **6.95** | 1.07 | 7.91 |
| A1_9NX0_XDON_rep2 | 2.44 | 9.01 | **6.57** | 1.06 | 7.48 | **6.43** | 1.07 | 7.66 |
| A1_9NX0_NATIVE_rep1 | 2.44 | 10.64 | **8.20** | 1.06 | 9.56 | **8.50** | 1.07 | 9.63 |
| A1_9NX0_NATIVE_rep2 | 2.44 | 9.83 | **7.39** | 1.06 | 9.31 | **8.26** | 1.07 | 9.52 |
| A2_7EGR_BLIND_rep1 | 4.45 | 5.99 | **1.54** | 2.40 | 4.58 | **2.19** | 2.12 | 4.39 |
| A2_7EGR_BLIND_rep2 | 4.45 | 7.17 | **2.71** | 2.40 | 5.75 | **3.35** | 2.12 | 5.63 |
| A2_7EGR_XDON_rep1 | 4.45 | 9.35 | **4.90** | 2.40 | 7.95 | **5.56** | 2.12 | 7.89 |
| A2_7EGR_XDON_rep2 | 4.45 | 7.81 | **3.35** | 2.40 | 6.90 | **4.50** | 2.12 | 6.79 |
| A2_7EGR_NATIVE_rep1 | 4.45 | 8.18 | **3.73** | 2.40 | 6.93 | **4.53** | 2.12 | 6.66 |
| A2_7EGR_NATIVE_rep2 | 4.45 | 8.86 | **4.41** | 2.40 | 7.44 | **5.04** | 2.12 | 7.63 |

`excess = best − floor` is the part the docking is responsible for.

**A2's all-atom floor is 4.453 Å — above the gate. The all-atom gate could not be passed at A2 no matter how good the docking was**; the minimum over all 20 deposited 2JUT models is still 3.248 Å. The informative number at A2 is the excess, and on backbone (floor 2.397 Å) the gate *was* reachable and still failed.


## (b) Is accuracy monotonic in disclosed information?

### A1 — 9NX0

| rung | R restraints | swarms | best all-atom (Å) rep1 / rep2 | pooled best | best backbone rep1 / rep2 | per-swarm-best median (Å) | top-1 by score (Å) |
|---|---|---|---|---|---|---|---|
| BLIND | 0 | 390/390 / 390/390 | 8.74 / 7.65 | **7.65** | 8.20 / 6.32 | 43.89 | 11.47 / 87.99 |
| XDON | 29 | 390/390 / 390/390 | 9.15 / 9.01 | **9.01** | 8.00 / 7.48 | 20.53 | 15.84 / 15.69 |
| NATIVE | 28 | 392/392 / 392/392 | 10.64 / 9.83 | **9.83** | 9.56 / 9.31 | 20.70 | 11.10 / 21.12 |

Best-RMSD ordering: BLIND 7.65 / XDON 9.01 / NATIVE 9.83 Å → **strictly monotonic in the WRONG direction** — more disclosed information, worse best pose (gap BLIND−NATIVE = -2.187 Å; the hypothesised direction would make this positive)

### A2 — 7EGR

| rung | R restraints | swarms | best all-atom (Å) rep1 / rep2 | pooled best | best backbone rep1 / rep2 | per-swarm-best median (Å) | top-1 by score (Å) |
|---|---|---|---|---|---|---|---|
| BLIND | 0 | 134/134 / 134/134 | 5.99 / 7.17 | **5.99** | 4.58 / 5.75 | 27.64 | 25.25 / 25.25 |
| XDON | 12 | 134/134 / 134/134 | 9.35 / 7.81 | **7.81** | 7.95 / 6.90 | 14.57 | 12.12 / 11.82 |
| NATIVE | 33 | 132/132 / 132/132 | 8.18 / 8.86 | **8.18** | 6.93 / 7.44 | 17.86 | 12.76 / 12.96 |

Best-RMSD ordering: BLIND 5.99 / XDON 7.81 / NATIVE 8.18 Å → **strictly monotonic in the WRONG direction** — more disclosed information, worse best pose (gap BLIND−NATIVE = -2.191 Å; the hypothesised direction would make this positive)

Two things move in **opposite** directions, which is the substance of the answer:

- **Best achievable pose gets worse with more information.** Restrained swarm placement concentrates sampling on the restrained patch; where that patch only partly covers the true footprint it pulls the search to a partly-wrong location. Measured: nearest swarm centre to the native ligand centroid is 6.08 Å for A2 BLIND but 10.67 Å for A2 XDON.
- **The *typical* pose gets markedly better with more information**, consistently at both anchors: per-swarm-best median A1 43.77 → 20.55 Å (BLIND → XDON), A2 27.64 → 14.57 / 17.86 Å. Restraints do concentrate sampling near the site; they just do not improve its best outcome.
- **The score-selected pose shows no consistent pattern.** A2 improves then flattens (top-1 by score 25.25 → 12.12 / 11.82 → 12.76 / 12.96 Å) while A1 is non-monotone (11.47 → 15.84 → 11.10 Å), with XDON the worst rung at A1 and the best at A2. ⚠ Independently of that, the column is not comparable across rungs at all: receptor restraints multiply the score by (1 + f_R), giving BLIND a ×1 ceiling and restrained rungs up to ×2 (`FINDINGS.md` F4). Re-scored on one common scale the A2 ordering reverses outright — see the rung-comparable table below. **This row must not be quoted as evidence about pose selection.**


## (c) How much of the 2.14 Å figure is attributable to disclosed information?

- **A1: -2.187 Å.** A negative number means full disclosure made the best pose *worse*.
- **A2: -2.191 Å.** A negative number means full disclosure made the best pose *worse*.

⚠ **The two anchors landing on the same 2-decimal value (-2.19 / -2.19 Å) is a coincidence, not a copied cell.** Each gap is computed independently from its own anchor's pooled best-RMSD: A1 = 7.6466 − 9.8340 = -2.1874 Å; A2 = 5.9901 − 8.1808 = -2.1907 Å. They differ by 0.0034 Å, which is why they are quoted to three decimals here. Differencing the 2-decimal values shown in the tables above instead gives A1 -2.18 and A2 -2.19 Å — the difference of rounded numbers, not the rounded difference; the values quoted here are the latter.

So on this ladder, **none of it.** Disclosure does not buy best-pose accuracy; it costs it. The 2.14 Å therefore cannot be attributed to the restraints.

What the ladder does **not** isolate, stated plainly: the WP1 2.14 Å run differs from S0's NATIVE rung in three ways at once, and S0 changed all three. Per AGENT.md/STATUS.md (the original run directory was lost with the previous pod and could not be re-inspected):

| | WP1 2.14 Å run | S0 A2 NATIVE |
|---|---|---|
| receptor | 6HY7 chain A — **one α9 ECD subunit**, 1701 atoms, 210 res (2–213): half an interface | 7EGR chains H,G — complete (+)/(−) interface |
| ligand | bound coordinates (bound–bound self-redock) | **apo NMR**, rigid |
| disclosure | 12 restraints, 12/12 identical to 6HY7 native contacts | 33 restraints = 7EGR's own native contacts |

Disclosure is maximal in both, so it is the one factor the comparison rules *out*. The remaining candidates are the apo-vs-bound ligand and the half-vs-full interface. Isolating them needs one run per factor, not more ladder rungs.


## Can the scoring function find the pose the sampler found?

LightDock scoring is **maximised**, so a correctly discriminating function gives a **negative** Spearman ρ against RMSD. Percentile 100 = best-scoring pose.

| run | R | poses | best RMSD (Å) | rank of best pose | **score percentile of best pose** | top-1 by score (Å) | best in top-10 | best in top-100 | ρ pose | ρ swarm | frac. swarms ρ<0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1_9NX0_BLIND_rep1 | 0 | 78000 | 8.74 | 6991 | **91.04** | 11.47 | 11.38 | 11.38 | -0.205 | -0.275 | 0.79 |
| A1_9NX0_BLIND_rep2 | 0 | 78000 | 7.65 | 2579 | **96.69** | 87.99 | 87.76 | 15.64 | -0.180 | -0.252 | 0.78 |
| A1_9NX0_XDON_rep1 | 29 | 78000 | 9.15 | 297 | **99.62** | 15.84 | 10.50 | 10.32 | -0.003 | 0.115 | 0.85 |
| A1_9NX0_XDON_rep2 | 29 | 78000 | 9.01 | 153 | **99.81** | 15.69 | 10.15 | 9.18 | 0.006 | 0.055 | 0.83 |
| A1_9NX0_NATIVE_rep1 | 28 | 78400 | 10.64 | 364 | **99.54** | 11.10 | 11.10 | 10.96 | -0.058 | -0.111 | 0.82 |
| A1_9NX0_NATIVE_rep2 | 28 | 78400 | 9.83 | 3200 | **95.92** | 21.12 | 20.96 | 11.43 | -0.080 | -0.158 | 0.81 |
| A2_7EGR_BLIND_rep1 | 0 | 26800 | 5.99 | 24753 | **7.64** | 25.25 | 25.07 | 24.79 | -0.080 | -0.274 | 0.60 |
| A2_7EGR_BLIND_rep2 | 0 | 26800 | 7.17 | 24648 | **8.03** | 25.25 | 25.21 | 24.97 | -0.077 | -0.387 | 0.60 |
| A2_7EGR_XDON_rep1 | 12 | 26800 | 9.35 | 24985 | **6.78** | 12.12 | 11.89 | 11.86 | -0.143 | -0.248 | 0.76 |
| A2_7EGR_XDON_rep2 | 12 | 26800 | 7.81 | 25990 | **3.03** | 11.82 | 11.77 | 11.76 | -0.125 | -0.289 | 0.74 |
| A2_7EGR_NATIVE_rep1 | 33 | 26400 | 8.18 | 3142 | **88.10** | 12.76 | 12.63 | 9.40 | -0.237 | -0.460 | 0.72 |
| A2_7EGR_NATIVE_rep2 | 33 | 26400 | 8.86 | 26177 | **0.85** | 12.96 | 12.78 | 9.13 | -0.217 | -0.363 | 0.67 |

Per anchor × rung:

| anchor | rung | R | best RMSD (Å) | score percentile of best pose | top-1 by score (Å) | ρ pose | ρ swarm |
|---|---|---|---|---|---|---|---|
| A1 | BLIND | 0 | 8.74 / 7.65 | 91.04 / 96.69 | 11.47 / 87.99 | -0.205 / -0.180 | -0.275 / -0.252 |
| A1 | XDON | 29 | 9.15 / 9.01 | 99.62 / 99.81 | 15.84 / 15.69 | -0.003 / 0.006 | 0.115 / 0.055 |
| A1 | NATIVE | 28 | 10.64 / 9.83 | 99.54 / 95.92 | 11.10 / 21.12 | -0.058 / -0.080 | -0.111 / -0.158 |
| A2 | BLIND | 0 | 5.99 / 7.17 | 7.64 / 8.03 | 25.25 / 25.25 | -0.080 / -0.077 | -0.274 / -0.387 |
| A2 | XDON | 12 | 9.35 / 7.81 | 6.78 / 3.03 | 12.12 / 11.82 | -0.143 / -0.125 | -0.248 / -0.289 |
| A2 | NATIVE | 33 | 8.18 / 8.86 | 88.10 / 0.85 | 12.76 / 12.96 | -0.237 / -0.217 | -0.460 / -0.363 |

### Same question with a rung-comparable score

The published `Scoring` column carries the `(1 + f_R)` restraint factor, which is identically 1 for BLIND and up to 2 for restrained rungs, so it cannot be compared across rungs (`FINDINGS.md` F4). Re-scoring the same poses with the restraints removed puts every rung on one scale. Percentile and ρ are estimated from an unbiased random sub-sample of poses; the estimator was validated on A2 BLIND rep1, where it gives 7.05 / −0.073 against 7.64 / −0.080 from the full ranking.

| run | R | median(restraint-free / published) | top-1 by **published** score → RMSD | top-1 by **restraint-free** score → RMSD | shift | published rank of restraint-free top-1 | best-pose percentile (restraint-free) |
|---|---|---|---|---|---|---|---|
| A1_9NX0_BLIND_rep1 | 0 | 1.0000 | 11.47 | 11.54 | +0.07 | 13 | 89.00 |
| A1_9NX0_BLIND_rep2 | 0 | 1.0000 | 87.99 | 87.91 | -0.08 | 18 | 96.30 |
| A1_9NX0_XDON_rep1 | 29 | 0.9686 | 15.84 | 26.24 | +10.40 | 287 | 95.55 |
| A1_9NX0_XDON_rep2 | 29 | 0.9695 | 15.69 | 24.87 | +9.18 | 152 | 99.25 |
| A1_9NX0_NATIVE_rep1 | 28 | 0.9632 | 11.10 | 11.32 | +0.22 | 10 | 97.85 |
| A1_9NX0_NATIVE_rep2 | 28 | 0.9580 | 21.12 | 21.02 | -0.09 | 11 | 89.20 |
| A2_7EGR_BLIND_rep1 | 0 | 0.9999 | 25.25 | 25.25 | 0.00 | 1 | 7.05 |
| A2_7EGR_BLIND_rep2 | 0 | 1.0000 | 25.25 | 25.54 | +0.29 | 2 | 6.75 |
| A2_7EGR_XDON_rep1 | 12 | 0.9231 | 12.12 | 15.24 | +3.11 | 33 | 3.85 |
| A2_7EGR_XDON_rep2 | 12 | 0.9231 | 11.82 | 11.81 | -0.01 | 11 | 2.35 |
| A2_7EGR_NATIVE_rep1 | 33 | 0.9334 | 12.76 | 29.10 | +16.34 | 156 | 51.80 |
| A2_7EGR_NATIVE_rep2 | 33 | 0.9361 | 12.96 | 29.27 | +16.31 | 68 | 0.80 |

The `median(restraint-free / published)` column is the measurement of F4: ≈1.0 wherever R = 0 and <1 wherever restraints are declared.


## ⚠ Which of these metrics actually reproduce

Each rung was run twice, identical except for the GSO seed, so rep1-vs-rep2 measures how much of each metric is noise (`FINDINGS.md` F7).

| metric | A1_BLIND | A1_NATIVE | A1_XDON | A2_BLIND | A2_NATIVE | A2_XDON | max Δ | verdict |
|---|---|---|---|---|---|---|---|---|
| per-swarm-best **median** | 43.77 / 44.05 | 20.90 / 20.52 | 20.55 / 20.45 | 28.54 / 27.50 | 17.81 / 17.92 | 14.67 / 14.52 | 1.04 | **reproduces** |
| per-swarm-best **p5** | 19.20 / 18.67 | 12.77 / 12.77 | 12.92 / 12.51 | 9.84 / 10.54 | 10.18 / 11.27 | 10.28 / 10.56 | 1.08 | **reproduces** |
| `best_rmsd` | 8.74 / 7.65 | 10.64 / 9.83 | 9.15 / 9.01 | 5.99 / 7.17 | 8.18 / 8.86 | 9.35 / 7.81 | 1.54 | usable, state a band |
| **top-1-by-score → RMSD** | 11.47 / 87.99 | 11.10 / 21.12 | 15.84 / 15.69 | 25.25 / 25.25 | 12.76 / 12.96 | 12.12 / 11.82 | 76.52 | **unusable** |
| percentile of best pose | 91.04 / 96.69 | 99.54 / 95.92 | 99.62 / 99.81 | 7.64 / 8.03 | 88.10 / 0.85 | 6.78 / 3.03 | 87.25 | **unusable** |

The two failing metrics are each defined by a **single pose** out of 26 400–78 000, so nothing averages the noise away.

The A1 BLIND blow-up is mechanical: the top-scoring pose's centroid sits **4.2 Å** from native in rep1 and **87.6 Å** in rep2. That is possible because A1's receptor is the full-length pair including the TMD — poses reach ≈140 Å of centroid displacement at A1 versus ≈59 Å across A2. On a receptor that large the best-vdW pose can sit almost anywhere on it, and the seed decides where.

**So: quote the per-swarm-best distribution and `best_rmsd` (with a ±1.5 Å band); do not interpret any score-selected single-pose number.** AGENT.md §2 already prescribed the per-swarm-best distribution with rank tests over raw maxima — this is an independent confirmation of it.

Spearman ρ has thousands of pairs behind it and is stable within a run, but its **sign is not uniform across runs**: negative (correct direction, |ρ| 0.06–0.46) in most, while A1 XDON rep1 is −0.003 pose-level and **+0.115 swarm-level**, and both A1 restrained runs go positive on the restraint-free score (+0.074, +0.041). Discrimination collapses specifically in A1's restrained rungs, not at A2. Quote ρ per run, never as "always negative".

**None of this touches the (a)–(c) verdicts above**, which rest on `best_rmsd`, the intrinsic floors and the per-swarm-best distribution.


## Across-anchor comparison, swarm count equalised

`best_rmsd` is a minimum, so it grows with swarm count and A1 (390–392) cannot be compared raw to A2 (132–134). Every run's per-swarm-best array is subsampled without replacement to the common minimum (**N = 132**), B = 20000, seed 20260730, and the minimum retaken inside each subsample.

| run | native N | boot best RMSD median (Å) | 5–95 % | P(best < 3 Å) | subsampled? |
|---|---|---|---|---|---|
| A1_9NX0_BLIND_rep1 | 390 | 10.15 | 8.74–14.45 | 0.000 | yes |
| A1_9NX0_BLIND_rep2 | 390 | 8.88 | 7.65–15.25 | 0.000 | yes |
| A1_9NX0_XDON_rep1 | 390 | 10.29 | 9.15–10.97 | 0.000 | yes |
| A1_9NX0_XDON_rep2 | 390 | 9.57 | 9.01–10.94 | 0.000 | yes |
| A1_9NX0_NATIVE_rep1 | 392 | 10.78 | 10.64–11.79 | 0.000 | yes |
| A1_9NX0_NATIVE_rep2 | 392 | 10.41 | 9.83–11.83 | 0.000 | yes |
| A2_7EGR_BLIND_rep1 | 134 | 5.99 | 5.99–5.99 | 0.000 | yes |
| A2_7EGR_BLIND_rep2 | 134 | 7.17 | 7.17–7.17 | 0.000 | yes |
| A2_7EGR_XDON_rep1 | 134 | 9.35 | 9.35–9.35 | 0.000 | yes |
| A2_7EGR_XDON_rep2 | 134 | 7.81 | 7.81–7.81 | 0.000 | yes |
| A2_7EGR_NATIVE_rep1 | 132 | 8.18 | 8.18–8.18 | 0.000 | no (already at N) |
| A2_7EGR_NATIVE_rep2 | 132 | 8.86 | 8.86–8.86 | 0.000 | no (already at N) |

## Instrument findings

Separate from the docking result, and the part that survives it — full detail with how each was established in `FINDINGS.md`:

1. **The 3 Å all-atom gate was unreachable at A2 by construction** (floor 4.453 Å; 3.248 Å for the best of 20 apo models). Two RgIA copies inside 7EGR itself differ by 2.526 Å all-atom, so the threshold sits barely above the ground truth's own noise.
2. **Apo-source choice alone can decide a 3 Å benchmark.** For the same toxin ImI, the achievable floor spans **1.02–4.57 Å** across three deposited apo entries (1G2G 1.02–2.13, 1IMI 1.76–2.74, 1CNL 4.02–4.57) — a spread wider than the threshold.
3. **`gso_<step>.out` pairs post-move coordinates with pre-move scores.** Proven by a consecutive-step experiment: E(P₁₀) matches gso_11's `Scoring` at median |Δ| 0.0081 (186/200 within 0.05) while the same-file pairing fails at median 5.68. `lgd_rank_swarm.py` sorts on that column, so LightDock's own ranking tool inherits it.
4. **Receptor restraints enter the score**, giving BLIND a ×1 ceiling and restrained rungs up to ×2. Verified from each run's `lightdock.info` and by recomputation (median recomputed/file = 1.0000 with restraints wired, 0.9257 / 0.9694 without). Raw scores are not comparable across rungs.
5. **Swarm placement biases but does not bound the search** — 82–98 % of glowworms end outside their initial swarm sphere, lowest in the restrained rungs.
