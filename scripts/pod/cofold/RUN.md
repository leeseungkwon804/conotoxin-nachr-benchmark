# Boltz-2 co-folding run — exact commands used 2026-08-03 (pod, RTX A4500 20 GB)

Environment as used: `boltz-2.2.1`, torch 2.4.1+cu124, python 3.11, CUDA on A4500.

⚠ **`--no_kernels` is REQUIRED on this GPU.** Without it boltz-2.2.1 aborts with
`ModuleNotFoundError: No module named 'cuequivariance_torch'` inside
`boltz/model/layers/triangular_mult.py:22`. This is the documented non-H100 kernel problem.

⚠ Install in a **venv**, not system python. Installing into system python downgraded
shared packages on the shared pod (scipy 1.17.1 -> 1.13.1, biopython 1.87 -> 1.84,
gemmi 0.7.5 -> 0.6.5). See STATUS.md 2026-08-03 incident entry.

```bash
python3 -m venv ~/venvs/boltz && ~/venvs/boltz/bin/pip install boltz
```

Single prediction (A2 took ~40 s on the A4500; first run also downloads ~7.5 GB of
weights + CCD into --cache):

```bash
boltz predict inputs/A2_7EGR.yaml --out_dir outputs --use_msa_server \
    --output_format pdb --cache <CACHE> --no_kernels
```

Robustness run (5 diffusion samples, fixed seed) — matches the 2-seed design on the docking side:

```bash
boltz predict inputs/A2_7EGR.yaml --out_dir outputs_ms --use_msa_server \
    --output_format pdb --cache <CACHE> --no_kernels --diffusion_samples 5 --seed 20260803
```

MSA comes from the public ColabFold server (api.colabfold.com). Inputs are already-public
deposited PDB sequences.

## Evaluation

```bash
python3 eval_cofold.py <pred.pdb> \
   <bench>/structures/lightdock_inputs/7EGR_receptorPAIR_HG.pdb \
   <bench>/structures/lightdock_inputs/7EGR_ligand_P.pdb  P  A2_eval.json
python3 eval_iface.py  <same four args> A2_eval_iface.json [CUT=10.0]
```

⚠ **Use `structures/lightdock_inputs/*_receptorPAIR_*.pdb`, NOT
`s0_gate/<run>/lightdock_receptor.pdb`.** The latter is in the origin-centred LightDock
frame (centroid 0,0,0) while the native ligand is in the deposited frame
(7EGR ligand centroid -46.74, -53.87, 3.01). Mixing them gives 0 native contacts.

⚠ The receptor is a **homodimer**, so chain assignment is ambiguous and both scripts
enumerate permutations and keep the minimum. Wrong assignment moved receptor CA-RMSD
from 0.825 A to 18.271 A at A2.
