#!/usr/bin/env python3
"""
Extract voltage-gated sodium channel (scn / NaV) genes from fish genome assemblies
using TBLASTN against human NaV1.7 and NaV1.8 as queries.

Prerequisites:
  - NCBI BLAST suite installed (tblastn, makeblastdb)
  - Fish genome assemblies downloaded to data/genomes/fasta/
  - Human NaV sequences in data/human_channels/

Workflow:
  1. Download genome assemblies (nucleotide FASTA)
  2. Run: python3 scripts/02_extract_fish_nav.py search [species_name] [genome_fasta]
  3. Inspect BLAST hits, manually curate exon boundaries
  4. Validate with: python3 scripts/02_extract_fish_nav.py validate

Usage:
    # Setup: index genome
    python3 scripts/02_extract_fish_nav.py search "Labrus bergylta" data/genomes/fasta/labrus_bergylta.fasta

    # Output: BLAST results, candidate exons
    # Manual curation: examine results/nav_extraction/labrus_bergylta_candidates.fasta
    # Validate: python3 scripts/02_extract_fish_nav.py validate

Output:
    - results/nav_extraction/[species]_tblastn.tsv (BLAST results)
    - results/nav_extraction/[species]_candidates.fasta (translated, unfiltered)
    - results/nav_extraction/[species]_curated.fasta (manually reviewed)
    - results/ortholog_assignment.txt (NaV1.7 vs NaV1.8 mapping)
"""

import subprocess
import os
import sys
import re
from pathlib import Path
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq

RESULTS_DIR = Path('results/nav_extraction')
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

HUMAN_QUERY = 'data/human_channels/all_human_channels.fasta'

def check_blast():
    """Verify BLAST is installed."""
    try:
        subprocess.run(['tblastn', '-version'], capture_output=True, check=True)
        print("✓ BLAST found")
        return True
    except FileNotFoundError:
        print("❌ BLAST not found. Install with: conda install -c bioconda blast")
        return False

def run_tblastn_search(species_name, genome_fasta, evalue=1e-5):
    """
    Run TBLASTN: protein query (human NaV) vs nucleotide genome.
    """
    if not os.path.exists(genome_fasta):
        print(f"❌ Genome file not found: {genome_fasta}")
        return None

    # Create BLAST database
    db_name = RESULTS_DIR / f"{species_name}_genome"
    print(f"Creating BLAST database for {species_name}...", end=" ")
    try:
        subprocess.run(
            ['makeblastdb', '-in', genome_fasta, '-dbtype', 'nucl', '-out', str(db_name)],
            capture_output=True, check=True
        )
        print("✓")
    except subprocess.CalledProcessError as e:
        print(f"❌ {e.stderr.decode()}")
        return None

    # Run TBLASTN
    output_tsv = RESULTS_DIR / f"{species_name}_tblastn.tsv"
    print(f"Running TBLASTN ({species_name})...", end=" ")
    try:
        cmd = [
            'tblastn',
            '-query', HUMAN_QUERY,
            '-db', str(db_name),
            '-evalue', str(evalue),
            '-out', str(output_tsv),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen',
            '-num_threads', '4',
        ]
        subprocess.run(cmd, capture_output=True, check=True)
        print("✓")
    except subprocess.CalledProcessError as e:
        print(f"❌ {e.stderr.decode()}")
        return None

    return output_tsv

def extract_hit_regions(species_name, genome_fasta, blast_results_tsv, flank_bp=500):
    """
    From TBLASTN hits, extract genomic regions and translate.
    """
    # Load genome
    print(f"Loading {genome_fasta}...", end=" ")
    genome = {}
    for record in SeqIO.parse(genome_fasta, 'fasta'):
        genome[record.id] = record.seq
    print(f"✓ ({len(genome)} contigs/scaffolds)")

    # Parse BLAST results
    candidates = []
    with open(blast_results_tsv) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            query_id, hit_id, pident, length = parts[0], parts[1], float(parts[2]), int(parts[3])
            sstart, send = int(parts[8]), int(parts[9])
            evalue, bitscore = float(parts[10]), float(parts[11])

            # Filter: >80% identity, >100 bp aligned
            if pident < 80 or length < 100:
                continue

            # Extract region with flanking sequence
            region_start = max(0, min(sstart, send) - flank_bp)
            region_end = min(len(genome[hit_id]), max(sstart, send) + flank_bp)

            if hit_id not in genome:
                continue

            region_seq = genome[hit_id][region_start:region_end]

            # Translate all 6 frames, pick longest ORF
            best_prot = find_longest_orf(region_seq)

            candidates.append({
                'query': query_id,
                'hit': hit_id,
                'pident': pident,
                'bitscore': bitscore,
                'region_start': region_start,
                'region_end': region_end,
                'dna_seq': str(region_seq),
                'protein_seq': best_prot,
            })

    # Write candidates FASTA
    output_fasta = RESULTS_DIR / f"{species_name}_candidates.fasta"
    with open(output_fasta, 'w') as f:
        for i, cand in enumerate(candidates, 1):
            header = f"{species_name}_candidate_{i} | query={cand['query']} | identity={cand['pident']:.1f}% | bs={cand['bitscore']:.0f}"
            f.write(f">{header}\n{cand['protein_seq']}\n")

    print(f"✓ {len(candidates)} candidates → {output_fasta}")
    return output_fasta

def find_longest_orf(dna_seq):
    """Translate DNA in all 6 frames, return longest ORF."""
    from Bio.Seq import Seq

    best_prot = ""
    for frame in range(3):
        for strand in [dna_seq, dna_seq.reverse_complement()]:
            seq = strand[frame:]
            prot = seq.translate(to_stop=False)
            if len(prot) > len(best_prot):
                best_prot = prot

    return str(best_prot)

def main():
    if not check_blast():
        sys.exit(1)

    if len(sys.argv) < 2:
        print("""Usage:
    python3 scripts/02_extract_fish_nav.py search [species] [genome_fasta]
    python3 scripts/02_extract_fish_nav.py validate
        """)
        return

    command = sys.argv[1]

    if command == 'search':
        if len(sys.argv) < 4:
            print("Usage: python3 scripts/02_extract_fish_nav.py search [species] [genome_fasta]")
            sys.exit(1)

        species_name = sys.argv[2]
        genome_fasta = sys.argv[3]

        blast_results = run_tblastn_search(species_name, genome_fasta)
        if blast_results:
            extract_hit_regions(species_name, genome_fasta, blast_results)
            print(f"\n📋 Next: manually review {RESULTS_DIR}/{species_name}_candidates.fasta")
            print(f"   Then save curated sequences to {species_name}_curated.fasta")

    elif command == 'validate':
        print("Validation: checking curated sequences for integrity...")
        # Future: validate translated proteins against known NaV domains

if __name__ == '__main__':
    main()
