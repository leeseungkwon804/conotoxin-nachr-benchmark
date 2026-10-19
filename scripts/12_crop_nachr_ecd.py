#!/usr/bin/env python3
"""
Crop each nAChR subunit's canonical UniProt sequence to its extracellular
ligand-binding domain (ECD), using boundaries fetched LIVE per-subunit from the
UniProt features API on 2026-07-21 (NOT assumed uniform -- signal peptides ranged
21-31 aa and ECD end boundaries ranged 233-250 aa across the 8 subunits checked).

This crop is valid ONLY for the alpha-conotoxin / nAChR axis. The NaV/CaV
pore-blocker axis (mu/omega/iota-conotoxins) CANNOT be cropped this way -- see
Decisions D5: the toxin binding site there is formed by 4 homologous
transmembrane repeats folding together, and even the cryo-EM-ordered region is
1200-1800 aa.

Usage:
    python3 scripts/12_crop_nachr_ecd.py

Output:
    data/human_channels/ecd_cropped/<name>_ECD.fasta
    results/nav_extraction/nachr_ecd_boundaries.csv
"""

import csv
from pathlib import Path

HUMAN_DIR = Path("data/human_channels")
OUT_DIR = HUMAN_DIR / "ecd_cropped"
RESULTS_DIR = Path("results/nav_extraction")

# name -> (fasta filename, ECD start, ECD end)  -- boundaries = UniProt "Topological
# domain: Extracellular" feature immediately following the "Signal" feature,
# fetched live via rest.uniprot.org/uniprotkb/<acc>.json on 2026-07-21.
ECD_BOUNDARIES = {
    "a7nAChR":  ("a7nAChR.fasta",  23, 233),   # P36544, signal 1-22
    "a9nAChR":  ("a9nAChR.fasta",  23, 238),   # Q9UGM1, signal 1-22
    "a10nAChR": ("a10nAChR.fasta", 25, 237),   # Q9GZZ6, signal 1-24
    "a3nAChR":  ("a3nAChR.fasta",  32, 250),   # P32297, signal 1-31
    "a6nAChR":  ("a6nAChR.fasta",  26, 239),   # Q15825, signal 1-25
    "b2nAChR":  ("b2nAChR.fasta",  26, 238),   # P17787, signal 1-25
    "b3nAChR":  ("b3nAChR.fasta",  22, 237),   # Q05901, signal 1-21
    "b4nAChR":  ("b4nAChR.fasta",  22, 239),   # P30926, signal 1-21
}


def read_fasta(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0]
    seq = "".join(l.strip() for l in lines[1:])
    return header, seq


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for name, (fasta_name, start, end) in ECD_BOUNDARIES.items():
        path = HUMAN_DIR / fasta_name
        if not path.exists():
            print(f"{name}: !! missing {path}")
            continue
        header, seq = read_fasta(path)
        full_len = len(seq)
        ecd_seq = seq[start - 1:end]  # UniProt 1-indexed -> python 0-indexed slice
        out_path = OUT_DIR / f"{name}_ECD.fasta"
        out_path.write_text(f">{name}_ECD | {header[1:]} | residues {start}-{end}\n{ecd_seq}\n")

        flag = "OK" if len(ecd_seq) == (end - start + 1) else "!! LENGTH MISMATCH"
        print(f"{name:10s} full={full_len:4d} aa -> ECD {start}-{end} ({len(ecd_seq)} aa)  {flag}")

        rows.append({
            "name": name, "fasta_source": fasta_name, "full_length": full_len,
            "ecd_start": start, "ecd_end": end, "ecd_length": len(ecd_seq),
        })

    out_csv = RESULTS_DIR / "nachr_ecd_boundaries.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n-> {out_csv}\n-> {OUT_DIR}/*.fasta ({len(rows)} files)")


if __name__ == "__main__":
    main()
