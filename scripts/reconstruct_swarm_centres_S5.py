#!/usr/bin/env python3
"""
Reconstruct approximate swarm-centre coordinates for Figure S5, in the DEPOSITED
coordinate frame, and validate the reconstruction against the canonical distances
in results/s0_ladder/swarm_coverage.json.

Why this exists: the production S0-ladder LightDock runs were launched with
setup.json's write_starting_positions:false, so individual swarm-centre coordinates
were never written to disk -- only the aggregate distances in swarm_coverage.json /
ladder.json survive. What IS on disk, per swarm, is every glowworm's generation-0
translation vector (pod_mirror_rest/conotoxin_bench/s0_gate/<run>/swarm_N/gso_0.out,
200 glowworms per swarm) in LightDock's own RECEPTOR-CENTRED working frame (receptor
Cα centroid at the origin -- see the coordinate-frame note in
scripts/make_structure_panels_S5.py). LightDock initialises each swarm's glowworms
clustered around that swarm's (unpersisted) centre, so the mean of the 200
generation-0 translations is a good approximate estimator of it. Translating that
mean by the run's own deposited-frame receptor Cα centroid (computed from the SAME
receptor.pdb copy bundled in the run directory, confirmed identical to
pod_mirror/structures/lightdock_inputs/*.pdb -- both have receptor CA centroid
[150.958, 148.338, 139.921] for 9NX0, checked at the top of this exercise) puts the
estimate back in the deposited frame, where it can be plotted alongside the real
deposited receptor and ligand structures.

Validation (printed below, and see the module docstring of
make_structure_panels_S5.py): the reconstructed nearest-swarm-centre distance
agrees with swarm_coverage.json's min_dist_centre_to_native_centroid_A to within
0.11-0.71 A across all six runs, and the count of swarms within the 20 A shell
(2x swarm_radius) agrees with n_centres_within_2radius_of_centroid to within 0-2 --
expected sampling noise from averaging 200 points per swarm (each swarm's
generation-0 glowworms have a per-axis SD of ~4.2-4.7 A, so the standard error of a
200-sample mean is ~0.3 A per axis), not a discrepancy with the recorded Table S10
values. This script only produces marker POSITIONS for the figure; every distance
NUMBER drawn on Figure S5 is read directly from swarm_coverage.json, not from this
reconstruction.

Run:  python3 scripts/reconstruct_swarm_centres_S5.py
Output: results/s0_ladder/swarm_centres_reconstructed_S5.json
"""
import json
import os

import numpy as np

ROOT = "/Volumes/dkyoo_SSD1/review_partime/conotoxin_project"
GATE = f"{ROOT}/pod_mirror_rest/conotoxin_bench/s0_gate"
DEP = f"{ROOT}/pod_mirror/structures/lightdock_inputs"
OUT_JSON = f"{ROOT}/results/s0_ladder/swarm_centres_reconstructed_S5.json"

RUNS = {
    "A1_BLIND": ("A1_9NX0_BLIND_rep1", "9NX0_receptorPAIR_DE.pdb", "9NX0_ligand_I.pdb"),
    "A1_XDON": ("A1_9NX0_XDON_rep1", "9NX0_receptorPAIR_DE.pdb", "9NX0_ligand_I.pdb"),
    "A1_NATIVE": ("A1_9NX0_NATIVE_rep1", "9NX0_receptorPAIR_DE.pdb", "9NX0_ligand_I.pdb"),
    "A2_BLIND": ("A2_7EGR_BLIND_rep1", "7EGR_receptorPAIR_HG.pdb", "7EGR_ligand_P.pdb"),
    "A2_XDON": ("A2_7EGR_XDON_rep1", "7EGR_receptorPAIR_HG.pdb", "7EGR_ligand_P.pdb"),
    "A2_NATIVE": ("A2_7EGR_NATIVE_rep1", "7EGR_receptorPAIR_HG.pdb", "7EGR_ligand_P.pdb"),
}

# canonical Table S10 / swarm_coverage.json values, for the printed validation only
EXPECTED_MIN_DIST = {
    "A1_BLIND": 13.87945953012374, "A1_XDON": 14.61464051288257, "A1_NATIVE": 13.62450970836635,
    "A2_BLIND": 6.075657865357973, "A2_XDON": 10.673381735786299, "A2_NATIVE": 8.812982954909849,
}
EXPECTED_N_WITHIN_2R = {
    "A1_BLIND": 10, "A1_XDON": 56, "A1_NATIVE": 56,
    "A2_BLIND": 23, "A2_XDON": 120, "A2_NATIVE": 68,
}


def read_ca(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(coords)


def read_all_atoms(path):
    coords = []
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(coords)


def parse_gso0(path):
    pts = []
    with open(path) as f:
        next(f)
        for line in f:
            line = line.strip()
            if not line.startswith("("):
                continue
            inner = line[1:line.index(")")]
            vals = [float(v) for v in inner.split(",")]
            pts.append(vals[:3])
    return np.array(pts)


def main():
    results = {}
    for tag, (rundir, recfile, ligfile) in RUNS.items():
        d = f"{GATE}/{rundir}"
        rec_ca_dep = read_ca(f"{DEP}/{recfile}")
        rec_centroid = rec_ca_dep.mean(axis=0)
        lig_atoms = read_all_atoms(f"{DEP}/{ligfile}")
        lig_centroid_dep = lig_atoms.mean(axis=0)
        lig_centroid_gso = lig_centroid_dep - rec_centroid

        swarm_dirs = sorted(
            (x for x in os.listdir(d) if x.startswith("swarm_") and not x.startswith("._")),
            key=lambda x: int(x.split("_")[1]))
        n = len(swarm_dirs)
        centers_gso = np.zeros((n, 3))
        for i, sd in enumerate(swarm_dirs):
            centers_gso[i] = parse_gso0(f"{d}/{sd}/gso_0.out").mean(axis=0)
        centers_dep = centers_gso + rec_centroid

        dist = np.linalg.norm(centers_gso - lig_centroid_gso, axis=1)
        n_within_20 = int((dist < 20).sum())
        print("%-10s n_swarms=%-4d recon_min=%.3f  canon_min=%.3f (diff %+.3f)  "
              "n<20A=%-3d (expected %d)" % (
                  tag, n, dist.min(), EXPECTED_MIN_DIST[tag],
                  dist.min() - EXPECTED_MIN_DIST[tag], n_within_20, EXPECTED_N_WITHIN_2R[tag]))

        nearest_i = int(dist.argmin())
        results[tag] = dict(
            n_swarms=n,
            centers_xyz_deposited=centers_dep.tolist(),
            nearest_point_xyz_deposited=centers_dep[nearest_i].tolist(),
            reconstructed_min_dist_A=float(dist.min()),
            canonical_min_dist_A_swarm_coverage_json=EXPECTED_MIN_DIST[tag],
        )

    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(results, f)
    print("\nwrote", OUT_JSON, "(%.1f KB)" % (os.path.getsize(OUT_JSON) / 1024))


if __name__ == "__main__":
    main()
