#!/usr/bin/env python3
"""
Colab-ready runner for the nAChR-axis Boltz-2 benchmark batch.

*** THIS SCRIPT IS FOR DK TO RUN HIMSELF IN COLAB, ON GPU. ***
It was authored and reviewed locally (no GPU here, and per the established
local-vs-Colab compute boundary -- see Lab notebook/Tooling.md -- Boltz-2
inference is a "heavy" job that must run in Colab). No predictions have
actually been executed as part of preparing this script.

Setup in a Colab notebook (one cell each):

    !pip install boltz[cuda] -U

    from google.colab import drive
    drive.mount('/content/drive')

Then upload/sync the YAML batch to Drive:
    /content/drive/MyDrive/Colab Notebooks/review_part_time/conotoxin/boltz_inputs/nachr_axis/*.yaml
  (= the contents of results/boltz_inputs/nachr_axis/ produced by script 13, copied
   or git-pulled onto Drive)

Then run this script as a Colab cell (or `!python 14_colab_boltz_runner.py`):

    python3 14_colab_boltz_runner.py

Behavior:
  - Iterates every .yaml in INPUT_DIR.
  - Skips a pair if its output directory already contains a completed structure
    (resumable across Colab session timeouts/disconnects).
  - Runs `boltz predict <yaml> --use_msa_server --out_dir <OUTPUT_DIR>/<pair_id>`.
  - Uses the MSA server (--use_msa_server) since no precomputed .a3m is provided
    in the YAMLs (see script 13 docstring).
  - Logs a manifest row (pair_id, status, wall time, output path) to
    OUTPUT_DIR/run_log.csv incrementally, so a crash mid-batch loses at most the
    row for the pair that was in flight.

No `properties: affinity` is used (see Decisions D11 -- invalid for peptide
binders); this only produces STRUCTURE predictions + Boltz's own confidence
metrics (confidence_*.json: ptm, iptm, complex_plddt, etc.), which script 15
uses as the axis-(b) ranking proxy against measured affinity.
"""

import csv
import subprocess
import sys
import time
from pathlib import Path

# Adjust this if running outside the expected Colab Drive layout.
DRIVE_ROOT = Path("/content/drive/MyDrive/Colab Notebooks/review_part_time/conotoxin")
INPUT_DIR = DRIVE_ROOT / "boltz_inputs" / "nachr_axis"
OUTPUT_DIR = DRIVE_ROOT / "boltz_outputs" / "nachr_axis"


def already_done(pair_out_dir: Path) -> bool:
    """A Boltz run is done if it produced at least one predicted structure file."""
    if not pair_out_dir.exists():
        return False
    return any(pair_out_dir.rglob("*.cif")) or any(pair_out_dir.rglob("*.pdb"))


def run_one(yaml_path: Path, out_dir: Path) -> tuple[bool, str]:
    cmd = [
        "boltz", "predict", str(yaml_path),
        "--use_msa_server",
        "--out_dir", str(out_dir),
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            return False, f"boltz exit {result.returncode}: {result.stderr[-500:]}"
        return True, "ok"
    except subprocess.TimeoutExpired:
        return False, "timeout (>1h)"
    except FileNotFoundError:
        return False, "boltz CLI not found -- did you run: pip install boltz[cuda] -U ?"


def main():
    if not INPUT_DIR.exists():
        print(f"!! {INPUT_DIR} not found.")
        print("   Upload/sync results/boltz_inputs/nachr_axis/*.yaml from the local")
        print("   project to this Drive path before running.")
        sys.exit(1)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log_path = OUTPUT_DIR / "run_log.csv"
    log_exists = log_path.exists()
    log_file = open(log_path, "a", newline="")
    log_writer = csv.writer(log_file)
    if not log_exists:
        log_writer.writerow(["pair_id", "status", "message", "seconds", "out_dir"])
        log_file.flush()

    yamls = sorted(INPUT_DIR.glob("*.yaml"))
    print(f"Found {len(yamls)} YAML inputs in {INPUT_DIR}")

    for yaml_path in yamls:
        pair_id = yaml_path.stem
        pair_out = OUTPUT_DIR / pair_id

        if already_done(pair_out):
            print(f"  [skip] {pair_id} -- already has output")
            continue

        print(f"  [run]  {pair_id} ...", end=" ", flush=True)
        t0 = time.time()
        ok, msg = run_one(yaml_path, pair_out)
        dt = time.time() - t0
        status = "ok" if ok else "FAILED"
        print(f"{status} ({dt:.0f}s)")
        if not ok:
            print(f"         {msg}")

        log_writer.writerow([pair_id, status, msg, round(dt, 1), str(pair_out)])
        log_file.flush()  # incremental save -- survives a Colab disconnect

    log_file.close()
    print(f"\n[done] log -> {log_path}")


if __name__ == "__main__":
    main()
