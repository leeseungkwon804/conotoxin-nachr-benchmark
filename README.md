# Benchmarking structural models of α-conotoxin–nAChR complexes

Analysis code for:

> **Benchmarking Structural Models of α-Conotoxin–Nicotinic Acetylcholine Receptor
> Complexes: Apo-Conformer Selection, Interface Information, and Deep-Learning
> Co-Folding.** *Toxins* (submitted 2026).

## What this repository contains

The full benchmark pipeline: restraint-driven docking at three levels of interface
information, nine rescoring functions, deep-learning co-folding, and the temporal
control that separates generalisation from memorisation of a training-set structure.

```
scripts/00_–03_    sequence retrieval, conotoxin curation, MSA
scripts/10_–16_    anchor contacts, PDB/UniProt numbering checks, nAChR ECD cropping,
                   Boltz input generation, scoring against the benchmark
scripts/make_*     figures and supplementary figures
scripts/make_supp_tables_round2.py, round3.py    supplementary tables
scripts/build_manuscript.py, assemble_figures_v3*.py   document and figure assembly
scripts/Run_Boltz2_nAChR_axis.ipynb   co-folding runs
results/           machine-readable outputs underlying Tables S1–S23 (JSON/CSV)
```

## Structures and tools

No structural coordinates are redistributed here. Everything is public:

| Anchor | PDB | Method |
|---|---|---|
| A1 — human α7 nAChR + α-conotoxin ImI | [9NX0](https://www.rcsb.org/structure/9NX0) | cryo-EM, 3.06 Å, released 2026-02-04 |
| A2 — *Aplysia* AChBP + α-conotoxin RgIA | [7EGR](https://www.rcsb.org/structure/7EGR) | X-ray, released 2022-04-06 |

Cross-donated restraint sets use 2BYP (A1) and 6HY7 (A2).

Tool versions, fixed because results depend on them: LightDock 0.9.4 (pip package;
its `__version__` string "3.0" is a branding artifact), Boltz-2 v2.2.1, Boltz-1
(training cutoff 2021-09-30, used as the temporal control), Chai-1 v0.6.1.

## Reproducing

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy pandas biopython matplotlib lightdock==0.9.4
python3 scripts/11_verify_pdb_uniprot_numbering.py
python3 scripts/12_crop_nachr_ecd.py
python3 scripts/10_anchor_contacts.py
python3 scripts/13_make_boltz_yaml.py
python3 scripts/15_score_vs_benchmark.py
```

Co-folding needs a GPU; inference settings (three recycling steps, 200 sampling steps,
five diffusion samples, structural kernels disabled, no MSA subsampling) are set in
`13_make_boltz_yaml.py` and the notebook, and are held identical across models.

## Licence

MIT (see `LICENSE`). Please cite the paper if you use this code.
