#!/usr/bin/env python3
"""
Does any swarm even sit on the native site?

A BLIND run can only reach the native pose if some swarm centre lies within
reach of it.  Glowworms are seeded inside a sphere of setup.json:swarm_radius
around their centre, so if the nearest swarm centre is far outside that radius
the run's best-RMSD reports a sampling gap, not a scoring or physics result.
This distinguishes the two, which the kill gate needs.

Usage: s0_swarm_coverage.py <native_ligand_pdb> <run_dir> [<run_dir> ...]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s0_lib import parse_pdb_atoms, receptor_shift, xyz_array

native_pdb = sys.argv[1]
nat_xyz = xyz_array(parse_pdb_atoms(native_pdb))
nat_centroid = nat_xyz.mean(0)

out = []
for run_dir in sys.argv[2:]:
    run_dir = os.path.abspath(run_dir)
    name = os.path.basename(run_dir)
    setup = json.load(open(os.path.join(run_dir, "setup.json")))
    shift, _, _ = receptor_shift(run_dir)
    centres_pdb = os.path.join(run_dir, "init", "swarm_centers.pdb")
    if not os.path.exists(centres_pdb):
        print(f"{name}: no init/swarm_centers.pdb")
        continue
    c = xyz_array(parse_pdb_atoms(centres_pdb)) + shift    # to deposited frame
    d_centroid = np.linalg.norm(c - nat_centroid, axis=1)
    # distance from each swarm centre to the closest native ligand atom
    d_atoms = np.sqrt(((c[:, None, :] - nat_xyz[None, :, :]) ** 2).sum(-1)).min(1)
    R = float(setup["swarm_radius"])
    rec = {
        "run": name,
        "n_swarm_centres": int(len(c)),
        "swarm_radius_A": R,
        "min_dist_centre_to_native_centroid_A": float(d_centroid.min()),
        "min_dist_centre_to_nearest_native_atom_A": float(d_atoms.min()),
        "n_centres_within_radius_of_centroid": int((d_centroid <= R).sum()),
        "n_centres_within_2radius_of_centroid": int((d_centroid <= 2 * R).sum()),
        "n_centres_touching_native_ligand": int((d_atoms <= R).sum()),
        "native_site_sampled": bool(d_centroid.min() <= R),
    }
    out.append(rec)
    print(f"{name}: centres {rec['n_swarm_centres']}, radius {R} A | "
          f"nearest centre to native centroid {rec['min_dist_centre_to_native_centroid_A']:.2f} A "
          f"({rec['n_centres_within_radius_of_centroid']} within radius, "
          f"{rec['n_centres_within_2radius_of_centroid']} within 2R) | "
          f"{rec['n_centres_touching_native_ligand']} centres within {R} A of a native atom | "
          f"native site sampled: {rec['native_site_sampled']}")

if out:
    dest = os.environ.get("COVERAGE_OUT")
    if dest:
        prev = json.load(open(dest)) if os.path.exists(dest) else []
        names = {r["run"] for r in out}
        json.dump([r for r in prev if r["run"] not in names] + out,
                  open(dest, "w"), indent=2)
        print("wrote", dest)
