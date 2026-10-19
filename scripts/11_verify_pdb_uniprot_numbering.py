#!/usr/bin/env python3
"""
Verify that PDB auth_seq_id numbering for each anchor's channel chain matches the
downloaded UniProt canonical sequence numbering -- do NOT assume they match. Cryo-EM
constructs commonly carry N-terminal purification tags (Strep/FLAG) which shift
auth_seq_id relative to native UniProt numbering; 6J8E's NaV1.2 entity sequence
visibly starts with such a tag ("MASWSHPQFEKGGG...FDYKDDDDK...").

Method: take the OBSERVED (auth_seq_id, resname) pairs straight from ATOM records
(not the entity_poly sequence, which may not line up 1:1 with auth_seq_id if loops
are unmodeled) and scan candidate integer offsets against the UniProt sequence,
picking the offset that maximizes exact residue-identity agreement. Reports the
match rate so a wrong/ambiguous mapping is visible rather than silently accepted.

Usage:
    python3 scripts/11_verify_pdb_uniprot_numbering.py

Output:
    results/nav_extraction/numbering_offsets.csv
"""

import csv
import glob
from pathlib import Path
import numpy as np

PDB_CACHE = Path("results/nav_extraction/pdb_cache")
HUMAN_DIR = Path("data/human_channels")

AA3to1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U",
}

# anchor -> (structure_format-agnostic file stem, channel_chain, uniprot_fasta_file)
CHANNEL_ANCHORS = [
    ("6J8E", "A", "NaV1.2.fasta"),
    ("7VFU", "A", "CaV2.2.fasta"),
    ("7MIX", "A", "CaV2.2.fasta"),
    ("8X91", "A", "CaV2.1.fasta"),
    ("9NX0", "A", "a7nAChR.fasta"),
    ("9NX1", "B", "a7nAChR.fasta"),
    ("6HY7", "A", "a9nAChR.fasta"),
    ("25II", "A", "NaV1.6.fasta"),
    ("9K6X", "A", "NET.fasta"),  # NET not in current fetch list -- flagged below if missing
]


def read_fasta_seq(path):
    lines = Path(path).read_text().splitlines()
    return "".join(l.strip() for l in lines if l and not l.startswith(">"))


def observed_residues_pdb(path, chain):
    """(auth_seq_id, one-letter resname) per residue, from legacy PDB CA atoms only."""
    out = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] != chain:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            resname = line[17:20].strip()
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            out[resnum] = AA3to1.get(resname, "X")
    return out


def observed_residues_mmcif(path, chain):
    out = {}
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            parts = line.split()
            if len(parts) < 21:
                continue
            if parts[3] != "CA":
                continue
            if parts[18] != chain:
                continue
            try:
                resnum = int(parts[16])
            except ValueError:
                continue
            resname = parts[17]
            out[resnum] = AA3to1.get(resname, "X")
    return out


def find_offset(observed, uniprot_seq, offset_range=range(-400, 401)):
    """observed: dict auth_seq_id -> one-letter aa. Try offsets so that
    uniprot_seq[auth_seq_id - offset - 1] == observed aa (1-indexed UniProt).
    Returns (best_offset, match_rate, n_compared)."""
    items = list(observed.items())
    best = (None, -1.0, 0)
    L = len(uniprot_seq)
    for off in offset_range:
        n_ok = 0
        n_cmp = 0
        for resnum, aa in items:
            idx = resnum - off - 1  # 0-indexed into uniprot_seq
            if 0 <= idx < L:
                n_cmp += 1
                if uniprot_seq[idx] == aa:
                    n_ok += 1
        if n_cmp >= 20:  # need enough overlap to trust the rate
            rate = n_ok / n_cmp
            if rate > best[1]:
                best = (off, rate, n_cmp)
    return best


def main():
    rows = []
    for pdb_id, chain, fasta_name in CHANNEL_ANCHORS:
        fasta_path = HUMAN_DIR / fasta_name
        if not fasta_path.exists():
            print(f"{pdb_id} chain {chain}: !! no local FASTA {fasta_path} -- SKIPPED, need to fetch")
            rows.append({"pdb_id": pdb_id, "chain": chain, "uniprot_file": fasta_name,
                         "offset": None, "match_rate": None, "n_compared": None, "note": "FASTA missing"})
            continue
        uniprot_seq = read_fasta_seq(fasta_path)

        pdb_path = PDB_CACHE / f"{pdb_id}.pdb"
        cif_path = PDB_CACHE / f"{pdb_id}.cif"
        if pdb_path.exists():
            observed = observed_residues_pdb(pdb_path, chain)
        elif cif_path.exists():
            observed = observed_residues_mmcif(cif_path, chain)
        else:
            print(f"{pdb_id}: !! no cached structure file found")
            continue

        offset, rate, n_cmp = find_offset(observed, uniprot_seq)
        flag = "OK" if (offset is not None and rate >= 0.9) else "!! CHECK"
        print(f"{pdb_id} chain {chain}: offset={offset}  match_rate={rate:.1%} (n={n_cmp})  {flag}")
        rows.append({"pdb_id": pdb_id, "chain": chain, "uniprot_file": fasta_name,
                     "offset": offset, "match_rate": round(rate, 4) if offset is not None else None,
                     "n_compared": n_cmp, "note": flag})

    out_csv = Path("results/nav_extraction/numbering_offsets.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n-> {out_csv}")


if __name__ == "__main__":
    main()
