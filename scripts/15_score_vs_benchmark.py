#!/usr/bin/env python3
"""
Score Boltz-2 nAChR-axis predictions against the measured affinity benchmark.

IMPORTANT SCOPE NOTE (Decisions D11): Boltz-2 does not output a predicted IC50
for peptide-protein pairs (that module is small-molecule-ligand-only). What we
actually have is Boltz's own STRUCTURE-CONFIDENCE metrics (confidence_score,
iptm, complex_plddt -- see docs/prediction.md output schema) for each toxin-
channel co-fold. This script tests them as an INDIRECT ranking proxy against
measured affinity. Confidence == "how sure is the model of its own structure",
NOT == "how strong is the binding" -- these are conceptually different, and a
weak/absent correlation here is a legitimate, reportable benchmark finding, not
a bug in this script.

Usage:
    python3 scripts/15_score_vs_benchmark.py [--mock]

    --mock   run on synthetic confidence values (since no Colab run exists yet)
             to validate the parsing/join/correlation logic end-to-end.
    (no flag) read real confidence_*.json files from
             results/boltz_outputs/nachr_axis/<pair_id>/predictions/<pair_id>/

Output:
    results/nachr_axis_scoring.csv
"""

import csv
import json
import re
import sys
from pathlib import Path

MANIFEST = Path("results/boltz_inputs/nachr_axis/manifest.csv")
BENCHMARK = Path("data/benchmark_affinity.csv")
BOLTZ_OUT = Path("results/boltz_outputs/nachr_axis")
SCORE_OUT = Path("results/nachr_axis_scoring.csv")

# maps our YAML "target" names -> the (target, target_species) spelling used in
# data/benchmark_affinity.csv, since the two files were built independently
TARGET_ALIASES = {
    "alpha9": [("alpha9 nAChR (homomeric)", "human")],
    "alpha9alpha10": [("alpha9alpha10 nAChR", "human"), ("alpha9alpha10 nAChR", "rat")],
    "alpha3beta2": [("alpha3beta2 nAChR", "human"), ("alpha3beta2 nAChR", "rat")],
    "alpha3beta4": [("alpha3beta4 nAChR", "human")],
    "alpha6_alpha3beta2beta3": [("alpha6/alpha3beta2beta3 nAChR", "human")],
    "alpha6_alpha3beta4": [("alpha6/alpha3beta4 nAChR", "human")],
}


def parse_affinity_to_pIC50(row):
    """Convert a benchmark_affinity.csv row's value/unit to pIC50 (-log10 M),
    higher = more potent. Returns None if not an IC50-type row we can convert
    (e.g. pIC50 already, or Kd -- handled separately)."""
    value_str = row["value"]
    unit = row["unit"].lower()
    m = re.match(r"[~]?([\d.]+)", value_str.replace("+/-", " ").strip())
    if not m:
        return None
    value = float(m.group(1))

    if "pic50" in unit:
        # value is already -log10(M); if unit string also has an approx molar
        # conversion e.g. "(~5 pM)" we ignore that and trust the pIC50 number.
        return value

    scale = {"nm": 1e-9, "um": 1e-6, "µm": 1e-6, "m": 1.0, "pm": 1e-12}
    for key, factor in scale.items():
        if unit.startswith(key):
            molar = value * factor
            if molar <= 0:
                return None
            import math
            return -math.log10(molar)
    return None


def normalize_toxin_name(name):
    """'Vc1.1 (ACV1)' -> 'vc1.1'; 'MVIIA (ziconotide)' -> 'mviia'. Strip any
    parenthetical synonym and whitespace/punctuation, keep only the base name."""
    base = re.sub(r"\(.*?\)", "", name).strip().lower()
    return base


def load_benchmark_lookup():
    """target_name (our alias) + species -> best-matching pIC50 (first match wins;
    if multiple species rows exist for one target we keep the human one when present)."""
    rows = list(csv.DictReader(open(BENCHMARK)))
    lookup = {}
    for our_name, aliases in TARGET_ALIASES.items():
        for bench_target, species in aliases:
            for r in rows:
                if r["target"].strip().lower() == bench_target.lower() and \
                   r["species"].strip().lower() == species.lower():
                    pic50 = parse_affinity_to_pIC50(r)
                    if pic50 is not None:
                        lookup.setdefault(our_name, []).append(
                            {"toxin": normalize_toxin_name(r["toxin"]), "species": species,
                             "pIC50": pic50, "raw_value": r["value"], "raw_unit": r["unit"]}
                        )
    return lookup


def load_confidence(pair_id):
    """Read the real Boltz confidence JSON for a pair, if it exists.

    NOTE (2026-07-22, verified against real Colab output): the actual nesting is
    one level deeper than the generic example in Boltz's docs/prediction.md --
    <out_dir>/<pair_id>/boltz_results_<pair_id>/predictions/<pair_id>/confidence_*.json
    -- so this searches recursively instead of assuming an exact fixed depth."""
    pair_dir = BOLTZ_OUT / pair_id
    if not pair_dir.exists():
        return None
    jsons = sorted(pair_dir.rglob("confidence_*_model_0.json"))
    if not jsons:
        return None
    return json.loads(jsons[0].read_text())


def mock_confidence(pair_id, seed_val):
    """Synthetic confidence values for pipeline validation only -- NOT real
    predictions. Deterministic function of pair_id so --mock is reproducible."""
    h = sum(ord(c) for c in pair_id) + seed_val
    return {
        "confidence_score": 0.5 + 0.4 * ((h % 97) / 97),
        "iptm": 0.4 + 0.5 * ((h % 89) / 89),
        "complex_plddt": 0.6 + 0.3 * ((h % 83) / 83),
    }


def spearman(xs, ys):
    """Plain-Python Spearman rank correlation (no scipy dependency needed)."""
    n = len(xs)
    if n < 3:
        return None
    def rank(vals):
        idx = sorted(range(len(vals)), key=lambda i: vals[i])
        ranks = [0] * len(vals)
        i = 0
        while i < len(idx):
            j = i
            while j + 1 < len(idx) and vals[idx[j + 1]] == vals[idx[i]]:
                j += 1
            avg_rank = (i + j) / 2 + 1
            for k in range(i, j + 1):
                ranks[idx[k]] = avg_rank
            i = j + 1
        return ranks
    rx, ry = rank(xs), rank(ys)
    mean_rx, mean_ry = sum(rx) / n, sum(ry) / n
    cov = sum((a - mean_rx) * (b - mean_ry) for a, b in zip(rx, ry))
    var_x = sum((a - mean_rx) ** 2 for a in rx)
    var_y = sum((b - mean_ry) ** 2 for b in ry)
    if var_x == 0 or var_y == 0:
        return None
    return cov / (var_x ** 0.5 * var_y ** 0.5)


def main():
    use_mock = "--mock" in sys.argv
    manifest = list(csv.DictReader(open(MANIFEST)))
    bench_lookup = load_benchmark_lookup()

    rows = []
    for m in manifest:
        pair_id, toxin, target = m["pair_id"], m["toxin"], m["target"]

        bench_matches = bench_lookup.get(target, [])
        toxin_norm = normalize_toxin_name(toxin)
        toxin_matches = [b for b in bench_matches if b["toxin"] == toxin_norm]
        measured = toxin_matches[0]["pIC50"] if toxin_matches else None

        conf = load_confidence(pair_id) if not use_mock else mock_confidence(pair_id, seed_val=0)
        if conf is None:
            rows.append({"pair_id": pair_id, "toxin": toxin, "target": target,
                         "status": "NO BOLTZ OUTPUT YET", "confidence_score": None,
                         "iptm": None, "complex_plddt": None,
                         "measured_pIC50": round(measured, 3) if measured is not None else None})
            continue

        rows.append({
            "pair_id": pair_id, "toxin": toxin, "target": target,
            "status": "mock" if use_mock else "real",
            "confidence_score": round(conf["confidence_score"], 4),
            "iptm": round(conf["iptm"], 4),
            "complex_plddt": round(conf["complex_plddt"], 4),
            "measured_pIC50": round(measured, 3) if measured is not None else None,
        })

    SCORE_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(SCORE_OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"{'pair_id':35s} {'status':20s} {'iptm':>8s} {'pLDDT':>8s} {'measured_pIC50':>15s}")
    for r in rows:
        print(f"{r['pair_id']:35s} {r['status']:20s} "
              f"{str(r['iptm']):>8s} {str(r['complex_plddt']):>8s} {str(r['measured_pIC50']):>15s}")

    scored = [r for r in rows if r["measured_pIC50"] is not None and r["iptm"] is not None]
    print(f"\n{len(scored)}/{len(rows)} pairs have both a confidence score and a measured affinity.")
    if len(scored) >= 3:
        rho_iptm = spearman([r["iptm"] for r in scored], [r["measured_pIC50"] for r in scored])
        rho_conf = spearman([r["confidence_score"] for r in scored], [r["measured_pIC50"] for r in scored])
        print(f"Spearman rho (ipTM vs measured pIC50):              {rho_iptm:.3f}" if rho_iptm is not None else "rho (ipTM): undefined (no variance)")
        print(f"Spearman rho (confidence_score vs measured pIC50):  {rho_conf:.3f}" if rho_conf is not None else "rho (confidence_score): undefined")
    else:
        print("Not enough scored pairs yet for a correlation (need >=3).")

    if use_mock:
        print("\n*** THESE ARE MOCK/SYNTHETIC CONFIDENCE VALUES -- for validating this "
              "script's parsing/join/correlation logic only. Run without --mock once "
              "script 14 has produced real Boltz output on Colab. ***")
    else:
        n_missing = sum(1 for r in rows if r["status"] == "NO BOLTZ OUTPUT YET")
        if n_missing:
            print(f"\n{n_missing}/{len(rows)} pairs have no Boltz output yet -- run script 14 on Colab first.")


if __name__ == "__main__":
    main()
