#!/usr/bin/env python3
"""
Contact-weighted BindPred: mean-pool the RECEPTOR side's ESM2 embedding only over
residues confirmed as toxin contacts (<=5A) in a REAL solved anchor structure,
instead of naive whole-ECD mean-pooling.

Contact positions (ECD-relative, 1-indexed) verified from anchor_contacts.csv
(filtered to real amino acids only -- an earlier version of that file included
water/HOH and ethylene-glycol/EDO crystallization additives as spurious "contacts",
caught and excluded here), converted via the PDB->UniProt offsets already verified
in scripts/11_verify_pdb_uniprot_numbering.py (+23 for a7 anchors 9NX0/9NX1, +25 for
a9 anchor 6HY7), then to ECD-crop-relative coordinates.

SCOPE, stated honestly: only alpha7 and alpha9 have a real solved anchor complex in
our set. alpha10, alpha3, alpha6, beta2, beta3, beta4 have NO direct structural
contact evidence here -- for heteromeric pairs, the non-anchored partner subunit
still uses naive whole-sequence pooling (unavoidable without a homology-transfer
step, which is a separate, not-yet-done extension). This script therefore only
covers the 9 pairs where at least one side has real contact data, for a fair,
honestly-scoped comparison against the naive whole-sequence baseline on the SAME
pairs (not the full 17).

For alpha7, 6 of 39 raw PDB contacts fall in the TM region (outside our ECD crop,
since 9NX0/9NX1 are full-length structures but our crop is ECD-only) and are
necessarily dropped -- our contact-weighted embedding is therefore based on the
ECD-visible portion of the real interface only, not the complete one.
"""

import csv
import torch
import esm

TOXINS = {
    "ImI":   "GCCSDPRCAWRC",
    "ImII":  "ACCSDRRCRWRC",
    "RgIA":  "GCCSDPRCRYRCR",
    "Vc1.1": "GCCSDPRCNYDHPEIC",
    "PeIA":  "GCCSHPACSVNHPELC",
    "Mr1.1": "GCCSHPACSVNNPDIC",
    "AuIB":  "GCCSYPPCFATNPDC",
}

SUBUNITS = {
    "a7nAChR":  "GEFQRKLYKELVKNYNPLERPVANDSQPLTVYFSLSLLQIMDVDEKNQVLTTNIWLQMSWTDHYLQWNVSEYPGVKTVRFPDGQIWKPDILLYNSADERFDATFHTNVLVNSSGHCQYLPPGIFKSSCYIDVRWFPFDVQHCKLKFGSWSYGGWSLDLQMQEADISGYIPNGEWDLVGIPGKRSERFYECCKEPYPDVTFTVTMRRRTLYY",
    "a9nAChR":  "AETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSS",
    "a10nAChR": "AEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAA",
}

# ECD-relative (1-indexed) contact positions, derived from real anchor structures
CONTACT_POS = {
    "a7nAChR": [36, 46, 47, 55, 57, 76, 77, 93, 109, 111, 112, 115, 117, 119, 145,
                147, 148, 149, 150, 151, 152, 164, 166, 186, 188, 189, 190, 191,
                192, 193, 195, 196, 197],  # from 9NX0 (all 5 copies) + 9NX1 (B-F), union
    "a9nAChR": [98, 152, 153, 154, 155, 156, 157, 195, 200, 202, 203, 204],  # from 6HY7
}

# only pairs where at least one side has real contact data (honest scope)
PAIRS = [
    ("ImI",   ["a7nAChR"]),
    ("ImII",  ["a7nAChR"]),
    ("RgIA",  ["a9nAChR"]),
    ("RgIA",  ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a9nAChR"]),
    ("Vc1.1", ["a9nAChR", "a10nAChR"]),
    ("PeIA",  ["a9nAChR", "a10nAChR"]),
    ("Mr1.1", ["a9nAChR", "a10nAChR"]),
    ("AuIB",  ["a7nAChR"]),
]


def per_residue_embedding(model, alphabet, batch_converter, seq, device):
    _, _, tokens = batch_converter([("x", seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[33], return_contacts=False)
    reps = out["representations"][33][0]
    return reps[1:len(seq) + 1]  # (L, 1280), excl BOS/EOS


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    cache = {}
    def get_reps(seq):
        if seq not in cache:
            cache[seq] = per_residue_embedding(model, alphabet, batch_converter, seq, device)
        return cache[seq]

    from catboost import CatBoostRegressor
    cb = CatBoostRegressor()
    cb.load_model("/root/.cache/huggingface/hub/models--hbp5181--BindPred/snapshots/"
                  "89a595d1c161a95783d904f5c5b3cce565e52848/ESM2_BindPred.cbm", format="cbm")

    import numpy as np
    rows = []
    for toxin_name, subunit_names in PAIRS:
        toxin_seq = TOXINS[toxin_name]
        lig_reps = get_reps(toxin_seq)
        lig_emb = lig_reps.mean(0).cpu().numpy()  # toxin: whole-sequence pool (unchanged, it's already all "contact")

        # receptor: per-subunit, contact-weighted where available, else whole-seq
        sub_embs = []
        contact_coverage = []
        for s in dict.fromkeys(subunit_names):
            reps = get_reps(SUBUNITS[s])
            if s in CONTACT_POS:
                idx = [p - 1 for p in CONTACT_POS[s]]  # 0-indexed
                sub_emb = reps[idx].mean(0).cpu().numpy()
                contact_coverage.append(f"{s}:contact({len(idx)}res)")
            else:
                sub_emb = reps.mean(0).cpu().numpy()
                contact_coverage.append(f"{s}:naive_full({reps.shape[0]}res)")
            sub_embs.append(sub_emb)
        rec_emb = np.mean(sub_embs, axis=0)  # average across subunits if heteromeric

        X = np.concatenate([lig_emb, rec_emb]).reshape(1, -1)
        pred_log10_kd = float(cb.predict(X)[0])
        pred_pKd = -pred_log10_kd

        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        cov_str = "+".join(contact_coverage)
        print(f"{pair_id:30s} pKd={pred_pKd:.3f}  [{cov_str}]")
        rows.append({
            "pair_id": pair_id, "toxin": toxin_name,
            "receptor_subunits": "+".join(subunit_names),
            "predicted_pKd_contact_weighted": round(pred_pKd, 4),
            "contact_coverage": cov_str,
        })

    with open("bindpred_contact_weighted_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> bindpred_contact_weighted_predictions.csv")


if __name__ == "__main__":
    main()
