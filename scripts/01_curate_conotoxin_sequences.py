#!/usr/bin/env python3
"""
Curate conotoxin sequences manually from ConoServer.

ConoServer (conopeptide.org) does not have a clean API, so this script provides
helper functions for manual curation workflow:

1. Visit conoserver.org
2. For each Conus species (geographus, striatus, tulipa, catus, purpurascens)
   - Filter by toxin type (alpha, mu, delta, omega)
   - Download FASTA sequences
   - Save to data/conotoxins/raw/conus_[species]_[type].fasta

3. Run this script to consolidate and QC

Usage:
    # First: manually download from ConoServer
    # Then: python3 scripts/01_curate_conotoxin_sequences.py consolidate

Output:
    - data/conotoxins/conotoxins_master.fasta (consolidated, deduplicated)
    - data/conotoxins/conotoxin_metadata.csv (sequence info table)
"""

import os
import re
from pathlib import Path
from collections import defaultdict

CONUS_SPECIES = ['geographus', 'striatus', 'tulipa', 'catus', 'purpurascens']
TOXIN_TYPES = ['mu', 'delta', 'omega', 'alpha']

RAW_DIR = Path('data/conotoxins/raw')
OUTPUT_DIR = Path('data/conotoxins')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def parse_fasta(fasta_file):
    """Parse FASTA file, return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_header:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences

def validate_conotoxin_cysteine_pattern(seq):
    """
    Conotoxins have characteristic disulfide patterns.
    Most contain 4–6 cysteines in conserved frameworks.
    Return True if sequence looks like a conotoxin (has ≥4 cysteines, 12–30 aa typical).
    """
    cys_count = seq.count('C')
    return 12 <= len(seq) <= 35 and cys_count >= 4

def consolidate_sequences():
    """Read all raw conotoxin files, deduplicate, output master FASTA."""

    if not RAW_DIR.exists():
        print(f"❌ No {RAW_DIR} found. Please download sequences from ConoServer first.")
        print(f"   Expected files like: data/conotoxins/raw/conus_geographus_mu.fasta")
        return

    master_seqs = {}  # {seq_id: (header, sequence, species, toxin_type)}
    all_files = list(RAW_DIR.glob('conus_*.fasta'))

    if not all_files:
        print(f"⚠ No FASTA files found in {RAW_DIR}")
        return

    for fasta_file in sorted(all_files):
        print(f"Reading {fasta_file.name}...", end=" ")
        try:
            seqs = parse_fasta(fasta_file)
            print(f"({len(seqs)} sequences)")

            for header, seq in seqs:
                if not seq or seq.count('X') > len(seq) * 0.1:
                    continue  # Skip empty or >10% unknown residues

                if not validate_conotoxin_cysteine_pattern(seq):
                    print(f"  ⚠ Questionable: {header[:60]} (only {seq.count('C')} cysteines, len={len(seq)})")
                    continue

                seq_upper = seq.upper()
                if seq_upper not in master_seqs:
                    master_seqs[seq_upper] = (header, seq_upper)

        except Exception as e:
            print(f"Error: {e}")
            continue

    # Write consolidated master file
    master_file = OUTPUT_DIR / 'conotoxins_master.fasta'
    with open(master_file, 'w') as f:
        for seq_id, (header, seq) in enumerate(master_seqs.items(), 1):
            f.write(f">{header}\n{seq}\n")

    print(f"\n✓ Consolidated {len(master_seqs)} unique sequences → {master_file}")

def download_instructions():
    """Print instructions for manual ConoServer download."""
    print("""
╔════════════════════════════════════════════════════════════════╗
║          ConoServer Manual Download Instructions               ║
╚════════════════════════════════════════════════════════════════╝

1. Go to https://conopeptide.org

2. For each Conus species, download toxin sequences:

   Species list:
   - Conus geographus
   - Conus striatus
   - Conus tulipa
   - Conus catus
   - Conus purpurascens

3. For each species, filter by toxin type and download:
   - μ-conotoxin (site 1 blockers)
   - δ-conotoxin (site 6 blockers)
   - ω-conotoxin (site 3 blockers) — esp. geographus/magus for MVIIA benchmark

4. Save files as:
   data/conotoxins/raw/conus_[species]_[type].fasta

   Example:
   data/conotoxins/raw/conus_geographus_mu.fasta
   data/conotoxins/raw/conus_geographus_delta.fasta
   data/conotoxins/raw/conus_geographus_omega.fasta
   data/conotoxins/raw/conus_striatus_mu.fasta
   ... etc

5. Return here and run:
   python3 scripts/01_curate_conotoxin_sequences.py consolidate

Questions? ConoServer has detailed species pages and sequence information.
If a species has no toxins of a certain type, that's OK — just skip it.

Target: 20–50 μ-sequences, 20–50 δ-sequences, 5–10 ω-sequences total.
    """)

if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == 'consolidate':
        consolidate_sequences()
    else:
        download_instructions()
