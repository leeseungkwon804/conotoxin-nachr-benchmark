#!/usr/bin/env python3
"""
Evaluate a co-folded complex against the deposited complex on a basis that can
sit in the same table as our rigid-docking numbers.

Docking numbers are ligand-RMSD in the receptor frame with NO superposition --
valid there because the receptor is held fixed and IS the native receptor.
A co-folding model predicts the receptor too, so the equivalent quantity is
ligand-RMSD AFTER superposing the predicted receptor onto the native receptor
(CAPRI L-RMSD). We also report receptor CA-RMSD and fnat so the two rows are not
silently compared on different footings.

Receptor chains are homomeric here, so the chain assignment is ambiguous; we
enumerate permutations and keep the one minimising receptor CA-RMSD (the same
optimal-mapping rule CAPRI/DockQ use). Getting this wrong inflated A2's receptor
CA-RMSD from 0.825 A to 18.271 A.

IMPORTANT: the rigid-body conformational floor that bounds the docking runs does
NOT bound a co-folding model -- it predicts the ligand conformer too. A1 floor is
2.439 A all-atom and Boltz-2 scored 1.130 A; A2 floor is 4.453 A and Boltz-2
scored 2.408 A. That is not a contradiction, it is the point, and any table
carrying both methods must footnote it.

Usage: eval_cofold.py <pred.pdb> <native_receptor.pdb> <native_ligand.pdb> <lig_chain> <out.json>
"""
import collections, itertools, json, sys
import numpy as np


def read(path):
    at = []
    for L in open(path):
        if not L.startswith(("ATOM", "HETATM")):
            continue
        if L[76:78].strip() == "H":
            continue
        at.append(dict(ch=L[21], resseq=L[22:27].strip(), resn=L[17:20].strip(),
                       name=L[12:16].strip(),
                       xyz=(float(L[30:38]), float(L[38:46]), float(L[46:54]))))
    return at


def residues(at):
    out = collections.OrderedDict()
    for a in at:
        out.setdefault(a["ch"], collections.OrderedDict()).setdefault(a["resseq"], []).append(a)
    return out


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))]) @ Wt
    return R, Pc, Qc, float(np.sqrt((((P - Pc) @ R - (Q - Qc)) ** 2).sum(1).mean()))


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))


pred_f, natrec_f, natlig_f, ligch, out_f = sys.argv[1:6]
pred, natrec, natlig = read(pred_f), read(natrec_f), read(natlig_f)
pr, nr = residues(pred), residues(natrec)
pred_rc = [c for c in pr if c != ligch]
nat_rc = list(nr.keys())


def ca_pairs(mapping):
    P, Q = [], []
    for pc, nc in mapping:
        prl, nrl = list(pr[pc].values()), list(nr[nc].values())
        for i in range(min(len(prl), len(nrl))):
            pa = {a["name"]: a["xyz"] for a in prl[i]}
            na = {a["name"]: a["xyz"] for a in nrl[i]}
            if "CA" in pa and "CA" in na:
                P.append(pa["CA"]); Q.append(na["CA"])
    return np.array(P), np.array(Q)


best = None
for perm in itertools.permutations(pred_rc):
    m = list(zip(perm, nat_rc))
    P, Q = ca_pairs(m)
    if len(P) == 0:
        continue
    R, Pc, Qc, r = kabsch(P, Q)
    if best is None or r < best[0]:
        best = (r, m, R, Pc, Qc, len(P))
rec_rmsd, mapping, R, Pc, Qc, n_ca = best
xform = lambda v: (np.asarray(v) - Pc) @ R + Qc

plig = list(pr[ligch].values())
nlig_map = list(residues(natlig).values())[0]
nlig = list(nlig_map.values())
nat_resseq = list(nlig_map.keys())
Pl, Ql, names = [], [], []
for i in range(min(len(plig), len(nlig))):
    pa = {a["name"]: a["xyz"] for a in plig[i]}
    na = {a["name"]: a["xyz"] for a in nlig[i]}
    for nm in pa:
        if nm in na:
            Pl.append(pa[nm]); Ql.append(na[nm]); names.append(nm)
Pl, Ql = np.array(Pl), np.array(Ql)
Plt = xform(Pl)
bb = np.array([nm in ("N", "CA", "C", "O") for nm in names])
ca = np.array([nm == "CA" for nm in names])

# fnat on native chain labels / native ligand residue numbering
nat_lig_at = [a for r in nlig for a in r]
Rn = np.array([a["xyz"] for a in natrec]); Ln = np.array([a["xyz"] for a in nat_lig_at])
rk = [(a["ch"], a["resseq"]) for a in natrec]
lkeys = []
for i, r in enumerate(nlig):
    lkeys += [nat_resseq[i]] * len(r)
d = np.linalg.norm(Rn[:, None, :] - Ln[None, :, :], axis=-1)
ii, jj = np.where(d < 5.0)
nat_c = set((rk[i], lkeys[j]) for i, j in zip(ii, jj))

pmap = {}
for pc, nc in mapping:
    prl, nrl = list(pr[pc].keys()), list(nr[nc].keys())
    for i in range(min(len(prl), len(nrl))):
        pmap[(pc, prl[i])] = (nc, nrl[i])
pred_rec_at = [a for c in pred_rc for r in pr[c].values() for a in r]
pred_lig_at, plkeys = [], []
for i, r in enumerate(plig):
    if i < len(nat_resseq):
        pred_lig_at += r; plkeys += [nat_resseq[i]] * len(r)
Rp = np.array([a["xyz"] for a in pred_rec_at]); Lp = np.array([a["xyz"] for a in pred_lig_at])
dp = np.linalg.norm(Rp[:, None, :] - Lp[None, :, :], axis=-1)
ii, jj = np.where(dp < 5.0)
pred_c = set()
for i, j in zip(ii, jj):
    a = pred_rec_at[i]
    k = pmap.get((a["ch"], a["resseq"]))
    if k:
        pred_c.add((k, plkeys[j]))

res = dict(pred=pred_f,
           chain_mapping={p: n for p, n in mapping},
           n_receptor_CA_matched=int(n_ca),
           receptor_CA_rmsd_after_superposition_A=round(rec_rmsd, 3),
           n_ligand_atoms_matched=int(len(Pl)),
           ligand_rmsd_allatom_A=round(rmsd(Plt, Ql), 3),
           ligand_rmsd_backbone_A=round(rmsd(Plt[bb], Ql[bb]), 3),
           ligand_rmsd_CA_A=round(rmsd(Plt[ca], Ql[ca]), 3),
           n_native_contacts=int(len(nat_c)), n_pred_contacts=int(len(pred_c)),
           n_shared_contacts=int(len(nat_c & pred_c)),
           fnat=round(len(nat_c & pred_c) / max(1, len(nat_c)), 3),
           note=("ligand RMSD is after superposing predicted receptor CA onto native receptor CA "
                 "(CAPRI L-RMSD), with optimal homomeric chain assignment. The rigid-body "
                 "conformational floor of the docking runs does NOT apply to a co-folding model."))
json.dump(res, open(out_f, "w"), indent=1)
for k, v in res.items():
    if k != "note":
        print("%-45s %s" % (k, v))
