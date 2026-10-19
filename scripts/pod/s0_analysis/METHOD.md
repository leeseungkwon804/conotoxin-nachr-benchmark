# S0 3-rung ladder — ligand-RMSD analysis method

Written 2026-07-30 by the P1 analysis worker. Everything here was measured on this pod;
nothing is quoted from memory.

## What is being measured

Each of the 12 S0 runs docks an **apo NMR** conotoxin onto a **rigid receptor pair** taken
from a deposited complex. The question is how close the docked ligand gets to the ligand's
**native bound pose in that same deposited complex**.

The reported number is an **interface-frame ligand-RMSD with no superposition**. The
receptor never moves (`use_anm: false` in every `setup.json`, verified), so the docked
ligand and the native ligand share a coordinate frame once LightDock's centring is undone.
Superposing would hide a translated pose, which is exactly the failure mode the gate is
looking for.

## Pose algebra

`lightdock3_setup.py` writes `lightdock_receptor.pdb` / `lightdock_ligand.pdb` with each
body translated so its own centroid is at the origin. Measured, both anchors:

| | receptor atoms | shift (Å) | component spread |
|---|---|---|---|
| A1 `A1_9NX0_*` | 6326 | `[151.029, 148.350, 140.669]` | 2.84e-14 |
| A2 `A2_7EGR_*` | 3276 | `[-44.865, -51.886, 16.485]` | 1.42e-14 |

The spread confirms it is a pure translation — the receptor was not rotated. So:

```
pose_deposited_frame = R(q) @ lightdock_ligand_xyz + t + receptor_shift
```

where `(t, q)` are read per glowworm from `swarm_*/gso_100.out`.

`t` is therefore the pose's ligand **centroid** in the receptor-centred frame, which is used
below as a sampling diagnostic.

## Quaternion convention — determined, not assumed

Four candidate conventions were tested by reproducing the output of LightDock's own
`lgd_generate_conformations.py` on a real swarm. Max per-atom deviation:

| convention | A1 (swarm_0, 12 poses) | A2 (swarm_0, 8 poses) |
|---|---|---|
| **`wxyz`** (file order = w,x,y,z; v' = q v q*) | **0.000808 Å** | **0.000814 Å** |
| `wxyz` transposed | 18.42 Å | 23.16 Å |
| `xyzw` | 21.01 Å | 24.73 Å |
| `xyzw` transposed | 20.08 Å | 22.78 Å |

`wxyz` wins by four orders of magnitude; the residual 8e-4 Å is PDB's 3-decimal coordinate
field. Stored in `_calib/*/calibration.json`.

## Independent cross-check of the RMSD itself

`s0_verify_rmsd.py` regenerates a full swarm (200 poses) as **real PDB files** via
`lgd_generate_conformations.py`, then recomputes RMSD by a code path that shares nothing
with the analytic one except the receptor shift (atoms matched by reading the generated PDB
text).

| anchor | poses compared | atom pairs (PDB path / analytic) | analytic best | PDB-path best | max &#124;diff&#124; |
|---|---|---|---|---|---|
| A1 swarm_0 | 200 | 90 / 90 | 47.9908 Å | 47.9908 Å | 6.8e-5 Å |
| A2 swarm_0 | 200 | 105 / 105 | 33.7573 Å | 33.7573 Å | 7.7e-5 Å |

Both **AGREE**.

## Atom correspondence — docked apo ligand vs native bound ligand

Residues are paired by Needleman–Wunsch on 3-letter residue names (not by residue-number
arithmetic, so a numbering offset or an amidation cap cannot silently shift the mapping),
then atoms are paired by name within each residue pair.

| anchor | docked ligand | native ligand | matched atoms | dropped |
|---|---|---|---|---|
| A1 | `apo_ligands/ImI/ImI_1IMI_m01.pdb`, NH2 stripped → chain A res 1–12, **90 atoms** | `9NX0_ligand_I.pdb`, chain I res 1–13, 91 atoms | **90 / 90 (100 % of docked)** | 1 native-only atom: the `NH2` cap N of res 13 |
| A2 | `apo_ligands/RgIA/RgIA_2JUT_m01.pdb` → chain A res 1–13, **105 atoms** | `7EGR_ligand_P.pdb`, chain P res 401–413, 105 atoms | **105 / 105 (all atoms both sides)** | none |

Residue names agree position-for-position in both anchors
(A1 `GLY CYS CYS SER ASP PRO ARG CYS ALA TRP ARG CYS`;
A2 `GLY CYS CYS SER ASP PRO ARG CYS ARG TYR ARG CYS ARG`), and A2's numbering offset is a
clean +400. So the "apo vs bound atom count may differ" worry does **not** bite here: no
side chain is missing on either side, and the only unmatched atom in the whole analysis is
one amide nitrogen.

## Anchor identities — fetched live from RCSB on 2026-07-30

| | A1 | A2 |
|---|---|---|
| PDB | **9NX0** | **7EGR** |
| title | Alpha7-nicotinic acetylcholine receptor bound to conotoxin ImI | Co-crystal structure of Ac-AChBPP in complex with RgIA |
| method / resolution | cryo-EM, **3.06 Å** (FSC 0.143) | X-ray, **2.5 Å** |
| released | **2026-02-04** | **2022-04-06** |
| receptor entity | Neuronal acetylcholine receptor subunit alpha-7, ***Homo sapiens***, CHRNA7, chains A–E | Soluble acetylcholine receptor (AChBP), ***Aplysia californica***, chains A,D,G,H,I,J |
| ligand entity | alpha-conotoxin ImI, *Conus imperialis*, chains F–J (5 copies) | RgIA |
| docked receptor pair | `9NX0_receptorPAIR_DE.pdb` — chains D,E, **392 residues each, range 1–478 (ECD + TMD, full length)**, 6326 heavy atoms | `7EGR_receptorPAIR_HG.pdb` — chains H,G, **206 residues each, range 19–224 (ECD only)**, 3276 heavy atoms |

**Two corrections to the dispatch brief.**

1. The brief frames A2 as the anchor with the complete (+)/(−) interface. Both anchors are
   subunit **pairs cut from a homopentamer**, so both present a complete (+)/(−)
   interface. A1 additionally carries the transmembrane domain.
2. A2 is **not a human receptor** — 7EGR is *Aplysia californica* AChBP, the soluble
   surrogate that this project's own literature audit criticised Gulsevin 2020 for relying
   on. A1 (9NX0) is the human α7 cryo-EM anchor. A2 remains a legitimate docking-difficulty
   test (real experimental coordinates, complete interface, orthosteric ligand) but must not
   be written up as a human-receptor result.

The receptor-scope difference is also the real cause of the swarm-count asymmetry (390 vs
134): A1's blind search must cover the TMD surface as well. Within an anchor the rungs are
matched, so the ladder comparison is unaffected; **across** anchors only the subsampled
numbers are quoted.

## Metrics

| metric | why |
|---|---|
| `best_rmsd` | min over all poses in the run. Sampling ceiling — "could the search reach the native pose at all". **This is the kill-gate number.** |
| `rmsd_of_top1_by_score` | the pose you would actually hand to a biologist. Separates sampling from scoring. |
| `best_rmsd_in_top_by_score[N]` | how deep into the score-ranked list you must dig. |
| `per_swarm_best_rmsd` | one value per swarm; the unit for the unequal-N bootstrap and for the rank tests. Per AGENT.md §2 this, not `z=(max−mean)/sd`, is the correct statistic. |
| `sampling_floor` | see below. |

## Sampling floor — separating a placement artefact from a physics result

RMSD can never be smaller than the pose's centroid displacement from the native centroid,
so that displacement is a hard floor. Both are recorded.

Swarm-centre coverage of the native site, all 12 runs (`swarm_coverage.json`,
`swarm_radius = 10.0 Å`):

| run group | swarm centres | nearest centre → native centroid | centres within 10 Å of a native **atom** | native centroid inside a swarm sphere |
|---|---|---|---|---|
| A1 BLIND | 390 | 13.88 Å | 4 | no |
| A1 XDON | 390 | 14.61 Å | 8 | no |
| A1 NATIVE | 392 | 13.62 Å | 12 | no |
| A2 BLIND | 134 | 6.08 Å | 15 | **yes** (5 centres) |
| A2 XDON | 134 | 10.67 Å | 57 | no |
| A2 NATIVE | 132 | 8.81 Å | 43 | **yes** (1 centre) |

On its own this would look like A1 cannot reach the native site. **It does not hold**, because
GSO glowworms leave their initial sphere: measured escape fraction (final centroid farther
from its own swarm centre than `swarm_radius`) is **0.982 in A1_BLIND_rep1 and 0.972 in
A2_BLIND_rep1**, and **0.819–0.982 across all 12 runs** (the low end is the restrained rungs,
which confine the search more). Initial placement is therefore a bias, not a bound, and the
achieved floor is what matters — now over complete swarm coverage:

- `A1_9NX0_BLIND_rep1` (390/390 swarms): min pose-centroid displacement **3.934 Å**
- `A2_7EGR_BLIND_rep1` (134/134 swarms): min pose-centroid displacement **2.145 Å**

A2 blind docking does place the peptide centroid within 2.15 Å of where it belongs while
still scoring 5.99 Å RMSD — i.e. right neighbourhood, wrong orientation. A1 blind reaches
3.93 Å of centroid displacement while scoring 8.74 Å. That distinction is worth carrying
into the write-up.

*(Correction: an earlier draft of this section reported escape fractions of 0.970 / 0.971 and
an A1 sampling floor of 6.74 Å. Those came from partial scans — 99 of 390 swarms for A1 —
taken while the runs were still being scored. The values above are from the completed
`ladder.json` at full swarm coverage.)*

## Unequal-N correction

`best_rmsd` is a minimum, so it grows with swarm count and A1/A2 raw values are not
comparable (390 vs 134). Every run's `per_swarm_best` array is subsampled **without
replacement** to the common minimum swarm count, `B = 20000`, seed `20260730`, and the min is
retaken inside each subsample. This yields a distribution of "best RMSD you would have seen
with N swarms". Runs already at the common N are flagged `boot_exact_no_subsampling` — for
them the distribution is a point mass and no correction is possible or needed.

Within an anchor the rungs were built at matched swarm counts by design
(A1 390/390/392, A2 134/134/132), so raw best-RMSD is compared directly and backed by a
Mann–Whitney test (tie-corrected, normal approximation — scipy is not installed here) on the
per-swarm-best distributions.

## Files

```
s0_analysis/
  bin/s0_lib.py             parsing, pose algebra, NW residue alignment, quaternion candidates
  bin/s0_calibrate.py       pins the quaternion convention against lgd_generate_conformations
  bin/s0_verify_rmsd.py     independent RMSD cross-check through real PDB files
  bin/s0_rmsd_run.py        scores one run  -> rmsd/<run>.json + .npz
  bin/s0_swarm_coverage.py  swarm-centre coverage of the native site
  bin/s0_ladder.py          assembles the ladder, bootstrap, rank tests, verdict -> ladder.{json,md}
  bin/s0_watch.sh           polls queue_v3_A*.log, scores each run the moment it exits
  _calib/                   calibration + verification scratch
  rmsd/                     per-run results (written as runs finish)
  swarm_coverage.json       all 12 runs
  ladder.json ladder.md     the verdict
  watch.log                 watcher trail
```

Completion is taken from `--- <run> rc=<n> END` in `s0_gate/queue_v3_A*.log`, never from
counting `gso_100.out` files, so a half-written file can never be scored.
