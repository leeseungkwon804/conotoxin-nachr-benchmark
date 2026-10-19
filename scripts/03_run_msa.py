#!/usr/bin/env python3
"""
Run multiple sequence alignments (MSA) on conotoxins and fish NaV proteins.

Prerequisites:
  - MAFFT installed (conda install -c bioconda mafft)

Usage:
    python3 scripts/03_run_msa.py conotoxins mu     # Align all μ-conotoxins
    python3 scripts/03_run_msa.py conotoxins delta  # Align all δ-conotoxins
    python3 scripts/03_run_msa.py nav              # Align fish + human NaV

Output:
    - results/alignments/conotoxin_mu_aligned.fasta
    - results/alignments/conotoxin_delta_aligned.fasta
    - results/alignments/nav_fish_human_aligned.fasta
    - results/alignments/*.html (interactive visualization)
"""

import subprocess
import os
import sys
from pathlib import Path
from Bio import AlignIO, SeqIO
import tempfile

ALIGNMENTS_DIR = Path('results/alignments')
ALIGNMENTS_DIR.mkdir(parents=True, exist_ok=True)

def check_mafft():
    """Verify MAFFT is installed."""
    try:
        subprocess.run(['mafft', '--version'], capture_output=True, check=True)
        print("✓ MAFFT found")
        return True
    except FileNotFoundError:
        print("❌ MAFFT not found. Install with: conda install -c bioconda mafft")
        return False

def run_mafft(input_fasta, output_afa, algorithm='auto'):
    """
    Align sequences with MAFFT.
    """
    print(f"Aligning {input_fasta}...", end=" ")
    try:
        cmd = [
            'mafft',
            f'--{algorithm}',  # 'auto' or 'linsi' (slow, accurate)
            '--thread', '4',
            input_fasta,
        ]
        with open(output_afa, 'w') as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
        print(f"✓")
        return output_afa
    except subprocess.CalledProcessError as e:
        print(f"❌ Error: {e}")
        return None

def align_conotoxins(toxin_type):
    """
    Align conotoxins of given type (mu, delta, omega).
    """
    input_file = Path(f'data/conotoxins/conotoxins_master.fasta')
    if not input_file.exists():
        print(f"❌ {input_file} not found. Run: python3 scripts/01_curate_conotoxin_sequences.py consolidate")
        return

    # Filter by toxin type
    toxin_type_lower = toxin_type.lower()
    filtered_file = ALIGNMENTS_DIR / f"conotoxin_{toxin_type_lower}_input.fasta"

    count = 0
    with open(filtered_file, 'w') as out:
        for record in SeqIO.parse(input_file, 'fasta'):
            # Heuristic: check if header contains the toxin type
            if toxin_type_lower in record.description.lower():
                out.write(f">{record.id} {record.description}\n{str(record.seq)}\n")
                count += 1

    if count == 0:
        print(f"⚠ No {toxin_type} conotoxins found in master FASTA")
        return

    output_afa = ALIGNMENTS_DIR / f"conotoxin_{toxin_type_lower}_aligned.afa"
    run_mafft(str(filtered_file), str(output_afa), algorithm='linsi')

    print(f"  {count} {toxin_type}-conotoxins aligned")
    generate_alignment_stats(str(output_afa), toxin_type)

def align_nav_sequences():
    """
    Align fish NaV (curated) + human NaV sequences.
    """
    # Collect all NaV sequences
    combined_file = ALIGNMENTS_DIR / "nav_combined_input.fasta"

    with open(combined_file, 'w') as out:
        # Add human NaV
        for record in SeqIO.parse('data/human_channels/all_human_channels.fasta', 'fasta'):
            out.write(f">Human_{record.id}\n{str(record.seq)}\n")

        # Add fish NaV (from curated files)
        for curated_file in Path('results/nav_extraction').glob('*_curated.fasta'):
            if curated_file.exists():
                for record in SeqIO.parse(curated_file, 'fasta'):
                    out.write(f">{record.id}\n{str(record.seq)}\n")

    # Run alignment with slower, more accurate algorithm for proteins
    output_afa = ALIGNMENTS_DIR / "nav_fish_human_aligned.afa"
    run_mafft(str(combined_file), str(output_afa), algorithm='linsi')

    generate_alignment_stats(str(output_afa), 'NaV')

def generate_alignment_stats(afa_file, name):
    """Generate summary statistics for alignment."""
    alignment = AlignIO.read(afa_file, 'fasta')
    print(f"\n  Alignment stats ({name}):")
    print(f"    Sequences: {len(alignment)}")
    print(f"    Length: {alignment.get_alignment_length()} bp")

    # Conservation score (% columns with same residue in ≥80% of sequences)
    conserved_cols = 0
    for col_idx in range(alignment.get_alignment_length()):
        col = alignment[:, col_idx]
        res_counts = {}
        for res in col:
            res_counts[res] = res_counts.get(res, 0) + 1

        max_count = max(res_counts.values()) if res_counts else 0
        if max_count / len(alignment) >= 0.8:
            conserved_cols += 1

    conservation_pct = 100 * conserved_cols / alignment.get_alignment_length()
    print(f"    Conservation (≥80%): {conservation_pct:.1f}%")

def main():
    if not check_mafft():
        sys.exit(1)

    if len(sys.argv) < 2:
        print("""Usage:
    python3 scripts/03_run_msa.py conotoxins mu
    python3 scripts/03_run_msa.py conotoxins delta
    python3 scripts/03_run_msa.py conotoxins omega
    python3 scripts/03_run_msa.py nav
        """)
        return

    target = sys.argv[1]

    if target == 'conotoxins' and len(sys.argv) >= 3:
        toxin_type = sys.argv[2]
        align_conotoxins(toxin_type)
    elif target == 'nav':
        align_nav_sequences()
    else:
        print(f"Unknown target: {target}")

if __name__ == '__main__':
    main()
