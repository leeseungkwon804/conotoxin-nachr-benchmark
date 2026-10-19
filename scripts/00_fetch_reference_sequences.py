#!/usr/bin/env python3
"""
Fetch reference human ion-channel sequences from UniProt.

All accessions below were verified live against the UniProt REST API
(rest.uniprot.org) on 2026-07-21. Do NOT edit them from memory --
an earlier version of this file had NaV1.7 = P08735, which is actually
maize glyceraldehyde-3-phosphate dehydrogenase. Verify with:
    curl -s https://rest.uniprot.org/uniprotkb/<ACC>.fasta | head -1

Usage:
    python3 scripts/00_fetch_reference_sequences.py

Output:
    data/human_channels/<name>.fasta  (one per channel)
    data/human_channels/all_human_channels.fasta  (combined)
"""

import os
import sys
import urllib.request
from pathlib import Path

# name -> (UniProt accession, gene, role in project)
CHANNELS = {
    # Primary drug targets (pain / epilepsy)
    "NaV1.7": ("Q15858", "SCN9A",   "primary target (pain)"),
    "NaV1.8": ("Q9Y5Y9", "SCN10A",  "primary target (pain)"),
    # Off-target / selectivity counter-screen
    "NaV1.5": ("Q14524", "SCN5A",   "off-target (cardiac) - safety"),
    "NaV1.2": ("Q99250", "SCN2A",   "off-target (CNS) + has KIIIA complex 6J8E"),
    "NaV1.6": ("Q9UQD0", "SCN8A",   "off-target (CNS)"),
    # Benchmark anchor (ziconotide / MVIIA -> CaV2.2)
    "CaV2.2": ("Q00975", "CACNA1B", "benchmark anchor (ziconotide) 7VFU/7MIX"),
    "CaV2.1": ("O00555", "CACNA1A", "off-target (P/Q) + has MVIIC complex 8X91"),
    # nAChR targets (alpha-conotoxin arm, if included)
    "a7nAChR":   ("P36544", "CHRNA7",  "nAChR target"),
    "a9nAChR":   ("Q9UGM1", "CHRNA9",  "nAChR target (analgesic)"),
    "a10nAChR":  ("Q9GZZ6", "CHRNA10", "nAChR target (analgesic)"),
    # nAChR subtypes needed for the Vc1.1 6-subtype selectivity benchmark panel
    # (accessions confirmed live via UniProt search API 2026-07-21, not from memory)
    "a3nAChR":   ("P32297", "CHRNA3",  "Vc1.1 panel: a3b2/a3b4 subtypes"),
    "b2nAChR":   ("P17787", "CHRNB2",  "Vc1.1 panel: a3b2 subtype"),
    "b4nAChR":   ("P30926", "CHRNB4",  "Vc1.1 panel: a3b4 subtype"),
    "a6nAChR":   ("Q15825", "CHRNA6",  "Vc1.1 panel: a6/a3b2b3, a6/a3b4 subtypes"),
    "b3nAChR":   ("Q05901", "CHRNB3",  "Vc1.1 panel: a6/a3b2b3 subtype"),
    # anchor 9K6X target (not a planned drug target, but a real anchor complex)
    "NET":       ("P23975", "SLC6A2",  "anchor 9K6X (AoIA-NET complex)"),
}

OUTPUT_DIR = Path("data/human_channels")


def fetch_fasta(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    combined = []
    print(f"Fetching {len(CHANNELS)} human channel sequences from UniProt REST...\n")
    for name, (acc, gene, role) in CHANNELS.items():
        try:
            fasta = fetch_fasta(acc)
        except Exception as e:
            print(f"  {name:10s} {acc}  FAILED: {e}")
            continue
        header = fasta.splitlines()[0]
        # sanity check: header should be the human gene, not something random
        ok = gene in header and "HUMAN" in header.upper()
        flag = "OK " if ok else "!! CHECK"
        seqlen = len("".join(fasta.splitlines()[1:]))
        print(f"  {name:10s} {acc}  {flag}  ({seqlen} aa)  {role}")
        if not ok:
            print(f"      -> unexpected header: {header}")
        out = OUTPUT_DIR / f"{name}.fasta"
        out.write_text(fasta)
        combined.append(fasta if fasta.endswith("\n") else fasta + "\n")

    (OUTPUT_DIR / "all_human_channels.fasta").write_text("".join(combined))
    print(f"\nWrote {len(combined)} sequences -> {OUTPUT_DIR/'all_human_channels.fasta'}")


if __name__ == "__main__":
    main()
