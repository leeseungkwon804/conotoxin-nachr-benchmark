#!/usr/bin/env python3
"""
Extract toxin<->channel contact residues from the ~9 real conopeptide-human-channel
anchor structures, using ONLY legacy-PDB coordinate parsing + numpy (no Biopython --
it fails to install in this environment due to a broken system expat/pyexpat, see
Lab notebook/Tooling.md). Runs locally (CPU, small files, seconds) -- NOT a GPU job,
see Decisions D-compute-boundary discussion 2026-07-21.

Chain/entity mapping verified live via RCSB Data API on 2026-07-21 (see ANCHORS below;
each entry: pdb_id, toxin chains, channel (pharmacological target) chains -- auxiliary
structural subunits like beta-subunits are intentionally excluded from "channel").

Usage:
    python3 scripts/10_anchor_contacts.py

Output:
    results/nav_extraction/anchor_contacts.csv   (all anchors, all contacts <= 5 Angstrom)
    results/nav_extraction/pdb_cache/*.pdb        (cached downloads)
"""

import os
import csv
import urllib.request
import urllib.error
from pathlib import Path
import numpy as np

CUTOFF_A = 5.0

# pdb_id -> (toxin_name, toxin_chains, channel_name, channel_chains)
# verified via data.rcsb.org/rest/v1/core/polymer_entity/<pdb>/<entity> on 2026-07-21
ANCHORS = {
    "6J8E": ("Mu-conotoxin KIIIA",       ["D"],           "NaV1.2 (SCN2A)",      ["A"]),
    "7VFU": ("Omega-conotoxin MVIIA",    ["E"],           "CaV2.2 (CACNA1B)",    ["A"]),
    "7MIX": ("Omega-conotoxin MVIIA",    ["B"],           "CaV2.2 (CACNA1B)",    ["A"]),
    "8X91": ("Omega-conotoxin MVIIC",    ["X"],           "CaV2.1 (CACNA1A)",    ["A"]),
    "9NX0": ("Alpha-conotoxin ImI",      ["F","G","H","I","J"], "alpha7 nAChR (CHRNA7)", ["A","B","C","D","E"]),
    "9NX1": ("Alpha-conotoxin ImII",     ["A"],           "alpha7 nAChR (CHRNA7)", ["B","C","D","E","F"]),
    "6HY7": ("Alpha-conotoxin RgIA",     ["B"],           "alpha9 nAChR ECD (CHRNA9)", ["A"]),
    "25II": ("Iota-conotoxin RXIA",      ["D"],           "NaV1.6 (SCN8A)",      ["A"]),
    "9K6X": ("conotoxin AoIA",           ["C","D"],       "NET (SLC6A2)",        ["A","B"]),
}

RESULTS_DIR = Path("results/nav_extraction")
PDB_CACHE = RESULTS_DIR / "pdb_cache"


def download_structure(pdb_id):
    """Try legacy PDB first; some large/many-chain cryo-EM entries (e.g. 9K6X) are
    404 in legacy format on RCSB, so fall back to mmCIF."""
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    pdb_out = PDB_CACHE / f"{pdb_id}.pdb"
    cif_out = PDB_CACHE / f"{pdb_id}.cif"
    if pdb_out.exists() and pdb_out.stat().st_size > 0:
        return pdb_out, "pdb"
    if cif_out.exists() and cif_out.stat().st_size > 0:
        return cif_out, "cif"

    try:
        print(f"  downloading {pdb_id}.pdb ...", end=" ")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", pdb_out)
        print(f"{pdb_out.stat().st_size} bytes")
        return pdb_out, "pdb"
    except urllib.error.HTTPError as e:
        print(f"legacy PDB unavailable ({e.code}), falling back to mmCIF ...", end=" ")
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.cif", cif_out)
        print(f"{cif_out.stat().st_size} bytes")
        return cif_out, "cif"


def parse_atoms_pdb(pdb_path):
    """Legacy PDB fixed-width columns. Returns chain_id -> [(resnum, resname, x, y, z)]."""
    chains = {}
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain_id = line[21]
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            resname = line[17:20].strip()
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            chains.setdefault(chain_id, []).append((resnum, resname, x, y, z))
    return chains


def parse_atoms_mmcif(cif_path):
    """mmCIF _atom_site loop, whitespace-delimited. Column order confirmed live for
    RCSB-served files on 2026-07-21: group_PDB id type_symbol label_atom_id label_alt_id
    label_comp_id label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x
    Cartn_y Cartn_z occupancy B_iso_or_equiv pdbx_formal_charge auth_seq_id auth_comp_id
    auth_asym_id auth_atom_id pdbx_PDB_model_num  (21 columns; auth_* used, matching the
    auth_asym_ids reported by the RCSB polymer_entity API used to build ANCHORS above)."""
    chains = {}
    with open(cif_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            parts = line.split()
            if len(parts) < 21:
                continue
            try:
                x = float(parts[10]); y = float(parts[11]); z = float(parts[12])
                resnum = int(parts[16])
            except ValueError:
                continue
            resname = parts[17]
            chain_id = parts[18]
            chains.setdefault(chain_id, []).append((resnum, resname, x, y, z))
    return chains


def parse_atoms(struct_path, fmt):
    return parse_atoms_pdb(struct_path) if fmt == "pdb" else parse_atoms_mmcif(struct_path)


def residue_atom_arrays(chain_atoms):
    """Group a chain's atom list into per-residue coordinate arrays.
    Returns list of (resnum, resname, coords[N,3])."""
    by_res = {}
    order = []
    for resnum, resname, x, y, z in chain_atoms:
        key = resnum
        if key not in by_res:
            by_res[key] = (resname, [])
            order.append(key)
        by_res[key][1].append((x, y, z))
    return [(k, by_res[k][0], np.array(by_res[k][1])) for k in order]


def contacts_between(toxin_residues, channel_residues, cutoff=CUTOFF_A):
    """For every toxin residue, find channel residues within cutoff (min atom-atom dist)."""
    rows = []
    for t_resnum, t_resname, t_coords in toxin_residues:
        for c_resnum, c_resname, c_coords in channel_residues:
            # pairwise distance matrix between the two small coord sets
            diff = t_coords[:, None, :] - c_coords[None, :, :]
            dmat = np.sqrt((diff ** 2).sum(axis=-1))
            mind = dmat.min()
            if mind <= cutoff:
                rows.append((t_resnum, t_resname, c_resnum, c_resname, round(float(mind), 2)))
    return rows


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = RESULTS_DIR / "anchor_contacts.csv"
    all_rows = []

    for pdb_id, (toxin_name, toxin_chains, channel_name, channel_chains) in ANCHORS.items():
        print(f"\n=== {pdb_id}  ({toxin_name} -> {channel_name}) ===")
        struct_path, fmt = download_structure(pdb_id)
        chains = parse_atoms(struct_path, fmt)

        missing = [c for c in toxin_chains + channel_chains if c not in chains]
        if missing:
            print(f"  !! WARNING: chains not found in parsed structure: {missing} "
                  f"(available: {sorted(chains.keys())})")

        for t_chain in toxin_chains:
            if t_chain not in chains:
                continue
            t_res = residue_atom_arrays(chains[t_chain])
            for c_chain in channel_chains:
                if c_chain not in chains:
                    continue
                c_res = residue_atom_arrays(chains[c_chain])
                rows = contacts_between(t_res, c_res)
                for t_resnum, t_resname, c_resnum, c_resname, mind in rows:
                    all_rows.append({
                        "pdb_id": pdb_id,
                        "toxin_name": toxin_name,
                        "toxin_chain": t_chain,
                        "toxin_resnum": t_resnum,
                        "toxin_resname": t_resname,
                        "channel_name": channel_name,
                        "channel_chain": c_chain,
                        "channel_resnum": c_resnum,
                        "channel_resname": c_resname,
                        "min_dist_A": mind,
                    })
                n_contact_res = len(set((r["channel_chain"], r["channel_resnum"]) for r in all_rows
                                         if r["pdb_id"] == pdb_id and r["toxin_chain"] == t_chain and r["channel_chain"] == c_chain))
                print(f"  toxin chain {t_chain} vs channel chain {c_chain}: "
                      f"{n_contact_res} contact residues (<= {CUTOFF_A} A)")

    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n[done] {len(all_rows)} contact rows across {len(ANCHORS)} anchors -> {out_csv}")


if __name__ == "__main__":
    main()
