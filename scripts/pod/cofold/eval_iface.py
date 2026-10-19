#!/usr/bin/env python3
"""Interface-local variant of the co-folding evaluation.

Superposing on the WHOLE receptor is unfair when the receptor carries a domain
that is irrelevant to the binding site and mobile. 9NX0 (A1) has a TMD and the
toxin binds the ECD interface, so a whole-receptor fit charges the ligand for TMD
placement error: A1 whole-receptor CA-RMSD is 7.896 A and the ligand then scores
3.099 A, while on interface-local superposition the same prediction gives
interface CA-RMSD 0.826 A and ligand 1.130 A.

So we also superpose on receptor residues that are actually at the site -- any
native receptor residue with a heavy atom within CUT A of the NATIVE ligand --
and re-measure. Both numbers get reported; neither replaces the other.

Usage: eval_iface.py <pred.pdb> <native_receptor.pdb> <native_ligand.pdb> <lig_chain> <out.json> [CUT]
"""
import collections, itertools, json, sys
import numpy as np

CUT = float(sys.argv[6]) if len(sys.argv) > 6 else 10.0


def read(p):
    at = []
    for L in open(p):
        if L.startswith(("ATOM", "HETATM")) and L[76:78].strip() != "H":
            at.append(dict(ch=L[21], rs=L[22:27].strip(), nm=L[12:16].strip(),
                           xyz=(float(L[30:38]), float(L[38:46]), float(L[46:54]))))
    return at


def res(at):
    o = collections.OrderedDict()
    for a in at:
        o.setdefault(a["ch"], collections.OrderedDict()).setdefault(a["rs"], []).append(a)
    return o


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ Wt))]) @ Wt
    return R, Pc, Qc, float(np.sqrt((((P - Pc) @ R - (Q - Qc)) ** 2).sum(1).mean()))


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))


pred_f, natrec_f, natlig_f, ligch, out_f = sys.argv[1:6]
pred, natrec, natlig = read(pred_f), read(natrec_f), read(natlig_f)
pr, nr = res(pred), res(natrec)
prc = [c for c in pr if c != ligch]
nrc = list(nr.keys())

# native interface residues: within CUT of any native ligand heavy atom
Ln = np.array([a["xyz"] for a in natlig])
iface = set()
for ch, rd in nr.items():
    for rs, atoms in rd.items():
        X = np.array([a["xyz"] for a in atoms])
        if np.linalg.norm(X[:, None, :] - Ln[None, :, :], axis=-1).min() < CUT:
            iface.add((ch, rs))


def pairs(mapping, only_iface):
    P, Q = [], []
    for pc, nc in mapping:
        pl, nl = list(pr[pc].values()), list(nr[nc].items())
        for i in range(min(len(pl), len(nl))):
            nrs, natoms = nl[i]
            if only_iface and (nc, nrs) not in iface:
                continue
            pa = {a["nm"]: a["xyz"] for a in pl[i]}
            na = {a["nm"]: a["xyz"] for a in natoms}
            if "CA" in pa and "CA" in na:
                P.append(pa["CA"]); Q.append(na["CA"])
    return np.array(P), np.array(Q)


best = None
for perm in itertools.permutations(prc):
    m = list(zip(perm, nrc))
    P, Q = pairs(m, True)
    if len(P) < 3:
        continue
    R, Pc, Qc, r = kabsch(P, Q)
    if best is None or r < best[0]:
        best = (r, m, R, Pc, Qc, len(P))
ir, mapping, R, Pc, Qc, n_if = best
xf = lambda v: (np.asarray(v) - Pc) @ R + Qc

plig = list(pr[ligch].values())
nlm = list(res(natlig).values())[0]
nlig = list(nlm.values())
Pl, Ql, nms = [], [], []
for i in range(min(len(plig), len(nlig))):
    pa = {a["nm"]: a["xyz"] for a in plig[i]}
    na = {a["nm"]: a["xyz"] for a in nlig[i]}
    for nm in pa:
        if nm in na:
            Pl.append(pa[nm]); Ql.append(na[nm]); nms.append(nm)
Pl, Ql = np.array(Pl), np.array(Ql)
Plt = xf(Pl)
bb = np.array([n in ("N", "CA", "C", "O") for n in nms])
ca = np.array([n == "CA" for n in nms])

out = dict(pred=pred_f,
           superposition="native interface receptor CA only (<%.1f A of native ligand)" % CUT,
           chain_mapping={p: n for p, n in mapping}, n_interface_CA_used=int(n_if),
           interface_CA_rmsd_A=round(ir, 3),
           ligand_rmsd_allatom_A=round(rmsd(Plt, Ql), 3),
           ligand_rmsd_backbone_A=round(rmsd(Plt[bb], Ql[bb]), 3),
           ligand_rmsd_CA_A=round(rmsd(Plt[ca], Ql[ca]), 3))
json.dump(out, open(out_f, "w"), indent=1)
for k, v in out.items():
    print("%-42s %s" % (k, v))
