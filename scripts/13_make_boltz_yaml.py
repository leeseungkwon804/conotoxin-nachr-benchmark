#!/usr/bin/env python3
"""
Build Boltz-2 input YAMLs for the alpha-conotoxin / nAChR benchmark axis.

IMPORTANT (2026-07-21, see Decisions D11): Boltz-2's `properties: affinity` block
is documented (docs/prediction.md) as valid ONLY for small-molecule ligand chains
("must be a ligand chain, not a protein, DNA or RNA"). Conotoxins are peptides
(protein chains), so NO affinity property is emitted here -- these YAMLs request
structure co-folding only. Confidence metrics (ipTM / interface pLDDT / PAE) from
the Boltz output are the ranking proxy for axis (b), not a predicted IC50.

Toxin sequences: mature peptide only (propeptide/signal cleaved), verified via
UniProt Peptide/Chain features or, for anchor toxins, the actual modeled ATOM
records of the real complex (not the entity_poly string, which for ImI included
a trailing unmodeled cap shown as "X" -- confirmed removed by checking coordinates).
PTMs (C-terminal amidation on all 5; 4-hydroxyproline + gamma-carboxyglutamate on
Vc1.1) are NOT encoded in the fold -- documented as a stated limitation rather than
guessed CCD codes.

Channel side: cropped nAChR ECDs from data/human_channels/ecd_cropped/ (task 2).
Homomeric/heteromeric targets use a minimal 2-subunit interface (matching the real
9NX0 pattern: each alpha-conotoxin contacts exactly 2 flanking subunits at one
subunit-subunit interface), NOT the full homopentamer -- documented per-pair where
the specific bound interface for a heteromeric subtype (e.g. alpha6/alpha3beta2beta3)
is a modeling simplification (the native receptor has multiple non-equivalent
interfaces and which one Vc1.1 prefers is not independently confirmed here).

Usage:
    python3 scripts/13_make_boltz_yaml.py

Output:
    results/boltz_inputs/nachr_axis/<pair_id>.yaml
    results/boltz_inputs/nachr_axis/manifest.csv
"""

import csv
from pathlib import Path

ECD_DIR = Path("data/human_channels/ecd_cropped")
OUT_DIR = Path("results/boltz_inputs/nachr_axis")

# mature peptide sequences (verified 2026-07-21; see script docstring for provenance)
TOXINS = {
    "ImI":   "GCCSDPRCAWRC",       # Conus imperialis, anchor 9NX0 (real complex, alpha7)
    "ImII":  "ACCSDRRCRWRC",       # Conus imperialis, anchor 9NX1 (real complex, alpha7)
    "RgIA":  "GCCSDPRCRYRCR",      # Conus regius, anchor 6HY7 (real complex, alpha9 ECD)
    "Vc1.1": "GCCSDPRCNYDHPEIC",   # Conus victoriae (UniProt P69747, Peptide 50-65)
    "PeIA":  "GCCSHPACSVNHPELC",   # Conus pergrandis (UniProt Q1L777, Peptide 22-37)
}

# target_name -> list of ECD-cropped subunit files forming the minimal docking interface
# (2 copies for homomeric alpha7/alpha9; 1+1 heteromeric interface otherwise)
TARGETS = {
    "alpha7":              (["a7nAChR", "a7nAChR"], "homomeric interface (matches 9NX0/9NX1 anchor pattern)"),
    "alpha9":              (["a9nAChR", "a9nAChR"], "homomeric interface"),
    "alpha9alpha10":       (["a9nAChR", "a10nAChR"], "heteromeric interface, minimal pair"),
    "alpha3beta2":         (["a3nAChR", "b2nAChR"], "heteromeric interface, minimal pair"),
    "alpha3beta4":         (["a3nAChR", "b4nAChR"], "heteromeric interface, minimal pair"),
    "alpha6_alpha3beta2beta3": (["a6nAChR", "b3nAChR"], "SIMPLIFICATION: native receptor is a heteropentamer with multiple non-equivalent interfaces (alpha6/alpha3/beta2/beta3); only one representative alpha6-beta3 interface is modeled here, not independently confirmed as the Vc1.1-preferred site"),
    "alpha6_alpha3beta4":  (["a6nAChR", "b4nAChR"], "SIMPLIFICATION: see alpha6_alpha3beta2beta3 note"),
}

# (toxin, target, source) -- source = benchmark row(s) this validates, or "anchor" for
# real-structure pose-validation cases (see Lab notebook/Structural anchors,
# Lab notebook/Benchmark set)
PAIRS = [
    ("ImI",   "alpha7",              "anchor 9NX0 (pose validation, no measured IC50 in our benchmark set)"),
    ("ImII",  "alpha7",              "anchor 9NX1 (pose validation)"),
    ("RgIA",  "alpha9",              "anchor 6HY7 (pose validation) + benchmark RgIA->a9a10 rat 5.2 nM (imperfect subtype match)"),
    ("RgIA",  "alpha9alpha10",       "benchmark: RgIA -> a9a10 (rat) IC50 5.2 nM"),
    ("Vc1.1", "alpha9",              "benchmark: Vc1.1 -> a9 (human) IC50 160.4 nM"),
    ("Vc1.1", "alpha9alpha10",       "benchmark: Vc1.1 -> a9a10 (human) IC50 906.9 nM; also rat 19 nM"),
    ("Vc1.1", "alpha3beta2",         "benchmark: Vc1.1 -> a3b2 (human) IC50 232.5 nM"),
    ("Vc1.1", "alpha3beta4",         "benchmark: Vc1.1 -> a3b4 (human) IC50 10.95 uM"),
    ("Vc1.1", "alpha6_alpha3beta2beta3", "benchmark: Vc1.1 -> a6/a3b2b3 (human) IC50 451.2 nM"),
    ("Vc1.1", "alpha6_alpha3beta4",  "benchmark: Vc1.1 -> a6/a3b4 (human) IC50 3.72 uM"),
    ("PeIA",  "alpha9alpha10",       "benchmark: PeIA -> a9a10 (rat) IC50 54.9 nM"),
    ("PeIA",  "alpha3beta2",         "benchmark: PeIA -> a3b2 (rat) IC50 97.5 nM"),
]


def read_ecd_seq(name):
    path = ECD_DIR / f"{name}_ECD.fasta"
    lines = path.read_text().splitlines()
    return "".join(l.strip() for l in lines[1:])


def make_yaml(toxin_id, toxin_seq, target_subunits, target_seqs):
    """Two-entity co-folding YAML: toxin (chain T) + channel interface (chains C1, C2...).
    No `properties: affinity` -- see module docstring / Decisions D11."""
    lines = ["version: 1", "sequences:"]
    lines.append("  - protein:")
    lines.append("      id: T")
    lines.append(f"      sequence: {toxin_seq}")
    for i, (subunit, seq) in enumerate(zip(target_subunits, target_seqs), 1):
        chain = chr(ord("A") + i - 1)
        lines.append(f"  - protein:")
        lines.append(f"      id: {chain}")
        lines.append(f"      sequence: {seq}")
    return "\n".join(lines) + "\n"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for toxin_name, target_name, source in PAIRS:
        toxin_seq = TOXINS[toxin_name]
        subunits, target_note = TARGETS[target_name]
        target_seqs = [read_ecd_seq(s) for s in subunits]

        pair_id = f"{toxin_name}__{target_name}".replace(".", "")
        yaml_text = make_yaml(pair_id, toxin_seq, subunits, target_seqs)
        out_path = OUT_DIR / f"{pair_id}.yaml"
        out_path.write_text(yaml_text)

        rows.append({
            "pair_id": pair_id, "toxin": toxin_name, "toxin_length": len(toxin_seq),
            "target": target_name, "target_subunits": "+".join(subunits),
            "target_note": target_note, "benchmark_source": source,
            "yaml_file": str(out_path),
        })
        print(f"{pair_id:35s}  toxin={len(toxin_seq)}aa  target={'+'.join(subunits)}")

    manifest_path = OUT_DIR / "manifest.csv"
    with open(manifest_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\n[done] {len(rows)} YAML files -> {OUT_DIR}")
    print(f"manifest -> {manifest_path}")
    print("\nNOTE: no `properties: affinity` block emitted (invalid for peptide binders "
          "per Boltz docs -- see Decisions D11). These YAMLs request structure co-folding "
          "only; run with `boltz predict --use_msa_server` (see script 14).")


if __name__ == "__main__":
    main()
