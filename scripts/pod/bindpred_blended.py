#!/usr/bin/env python3
"""
Blended BindPred: receptor embedding = mean(whole-sequence pool, contact-only pool),
instead of full replacement. Tests whether a compromise retains naive's advantage
while adding some structural signal, given pure contact-replacement made things
WORSE (not just neutral) -- likely because it pushes the input far from the
statistical distribution BindPred's pretrained CatBoost regressor was calibrated on
(whole-protein-averaged ESM2 embeddings), see bindpred_contact_weighted.py findings.
"""

import csv
import numpy as np
import torch
import esm
from catboost import CatBoostRegressor

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

CONTACT_POS = {
    "a7nAChR": [36, 46, 47, 55, 57, 76, 77, 93, 109, 111, 112, 115, 117, 119, 145,
                147, 148, 149, 150, 151, 152, 164, 166, 186, 188, 189, 190, 191,
                192, 193, 195, 196, 197],
    "a9nAChR": [98, 152, 153, 154, 155, 156, 157, 195, 200, 202, 203, 204],
}

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
    return reps[1:len(seq) + 1]


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    cache = {}
    def get_reps(seq):
        if seq not in cache:
            cache[seq] = per_residue_embedding(model, alphabet, batch_converter, seq, device)
        return cache[seq]

    cb = CatBoostRegressor()
    cb.load_model("/root/.cache/huggingface/hub/models--hbp5181--BindPred/snapshots/"
                  "89a595d1c161a95783d904f5c5b3cce565e52848/ESM2_BindPred.cbm", format="cbm")

    rows = []
    for toxin_name, subunit_names in PAIRS:
        toxin_seq = TOXINS[toxin_name]
        lig_emb = get_reps(toxin_seq).mean(0).cpu().numpy()

        sub_embs = []
        for s in dict.fromkeys(subunit_names):
            reps = get_reps(SUBUNITS[s])
            whole_emb = reps.mean(0).cpu().numpy()
            if s in CONTACT_POS:
                idx = [p - 1 for p in CONTACT_POS[s]]
                contact_emb = reps[idx].mean(0).cpu().numpy()
                blended = (whole_emb + contact_emb) / 2.0
                sub_embs.append(blended)
            else:
                sub_embs.append(whole_emb)
        rec_emb = np.mean(sub_embs, axis=0)

        X = np.concatenate([lig_emb, rec_emb]).reshape(1, -1)
        pred_pKd = -float(cb.predict(X)[0])

        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        print(f"{pair_id:30s} blended_pKd={pred_pKd:.3f}")
        rows.append({"pair_id": pair_id, "toxin": toxin_name,
                     "receptor_subunits": "+".join(subunit_names),
                     "predicted_pKd_blended": round(pred_pKd, 4)})

    with open("bindpred_blended_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> bindpred_blended_predictions.csv")


if __name__ == "__main__":
    main()
