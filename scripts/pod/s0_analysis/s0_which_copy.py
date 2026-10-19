#!/usr/bin/env python3
"""
Which ligand copy actually sits at the receptor pair we docked onto?

Both anchors are homopentamers with several bound toxin copies (9NX0: ImI in
chains F-J; 7EGR: RgIA in several chains).  The docked receptor is ONE subunit
pair.  If the native reference PDB is a different copy -- one bound at a
different interface of the same pentamer -- then every RMSD in the analysis is
measured against a toxin sitting tens of angstroms away and the whole result is
meaningless.  This checks it by contact counting rather than trusting a filename.

The correct copy is the one with 5 A contacts to BOTH chains of the pair, since
alpha-conotoxins bind the (+)/(-) subunit interface.

Usage: s0_which_copy.py <full_receptor_pdb> <pair_chain_ids> <ligand_multichain_pdb> [<ref_pdb>]
   e.g. s0_which_copy.py 9NX0_receptor_ABCDE.pdb DE 9NX0_all_ligands.pdb 9NX0_ligand_I.pdb
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import parse_pdb_atoms

rec_pdb, pair, lig_pdb = sys.argv[1], sys.argv[2].upper(), sys.argv[3]
ref_pdb = sys.argv[4] if len(sys.argv) > 4 else None
CUT = 5.0

rec = parse_pdb_atoms(rec_pdb)
lig = parse_pdb_atoms(lig_pdb)
rec_by_chain = {}
for a in rec:
    rec_by_chain.setdefault(a["chain"], []).append(a["xyz"])
rec_by_chain = {k: np.asarray(v) for k, v in rec_by_chain.items()}
lig_by_chain = {}
for a in lig:
    lig_by_chain.setdefault(a["chain"], []).append(a["xyz"])
lig_by_chain = {k: np.asarray(v) for k, v in lig_by_chain.items()}

print(f"receptor {os.path.basename(rec_pdb)}: chains "
      f"{ {k: len(v) for k, v in sorted(rec_by_chain.items())} }")
print(f"ligand   {os.path.basename(lig_pdb)}: chains "
      f"{ {k: len(v) for k, v in sorted(lig_by_chain.items())} }")
print(f"docked receptor pair: {' + '.join(pair)}   contact cutoff {CUT} A\n")

rows = []
for lc, L in sorted(lig_by_chain.items()):
    counts = {}
    for rc, R in sorted(rec_by_chain.items()):
        d = np.sqrt(((L[:, None, :] - R[None, :, :]) ** 2).sum(-1))
        counts[rc] = int((d < CUT).sum())
    pair_total = sum(counts.get(c, 0) for c in pair)
    both = all(counts.get(c, 0) > 0 for c in pair)
    other = sum(v for k, v in counts.items() if k not in pair)
    rows.append((lc, counts, pair_total, both, other))
    flag = "  <== at the docked pair" if both and pair_total > 0 and other < pair_total else ""
    print(f"  ligand chain {lc}: contacts {counts}  "
          f"| to pair {pair}: {pair_total}  both-chains: {both}  elsewhere: {other}{flag}")

best = max(rows, key=lambda r: (r[3], r[2]))
print(f"\nCopy bound at pair {pair}: chain {best[0]} "
      f"({best[2]} atom contacts to the pair, {best[4]} elsewhere)")

if ref_pdb:
    ref = parse_pdb_atoms(ref_pdb)
    ref_chains = sorted({a["chain"] for a in ref})
    print(f"reference file {os.path.basename(ref_pdb)} holds chain(s) {ref_chains}, "
          f"{len(ref)} atoms")
    ok = len(ref_chains) == 1 and ref_chains[0] == best[0]
    if not ok:
        # a renamed extraction is still fine if the coordinates match
        R = np.asarray([a["xyz"] for a in ref])
        for lc, L in lig_by_chain.items():
            if len(L) == len(R) and np.allclose(np.sort(L, axis=0), np.sort(R, axis=0),
                                                atol=1e-3):
                print(f"  -> coordinates identical to ligand chain {lc}")
                ok = (lc == best[0])
                break
    print("VERDICT:", "reference is the correct copy" if ok
          else f"*** MISMATCH: reference is not the copy at pair {pair} ***")
