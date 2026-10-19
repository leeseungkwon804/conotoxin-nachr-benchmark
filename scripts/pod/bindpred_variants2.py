#!/usr/bin/env python3
"""
Two more engineering variants, both trying to sidestep BindPred's CatBoost
distributional-mismatch problem found in the contact-only/blended tests:

A) Pure embedding cosine similarity (toxin whole-seq vs receptor, contact-weighted
   by inverse-distance) -- bypasses the pretrained CatBoost regressor entirely, so
   whatever distributional calibration it has can't distort the signal.
B) Distance-weighted contact pooling (soft weights = 1/min_dist_A, not binary
   in/out) fed into CatBoost -- refines the contact-only idea with a softer,
   more informative signal than a hard mask.
"""

import csv
import json
import numpy as np
import torch
import esm
from catboost import CatBoostRegressor

TOXINS = {
    "ImI":   "GCCSDPRCAWRC", "ImII":  "ACCSDRRCRWRC", "RgIA":  "GCCSDPRCRYRCR",
    "Vc1.1": "GCCSDPRCNYDHPEIC", "PeIA":  "GCCSHPACSVNHPELC",
    "Mr1.1": "GCCSHPACSVNNPDIC", "AuIB":  "GCCSYPPCFATNPDC",
}
SUBUNITS = {
    "a7nAChR":  "GEFQRKLYKELVKNYNPLERPVANDSQPLTVYFSLSLLQIMDVDEKNQVLTTNIWLQMSWTDHYLQWNVSEYPGVKTVRFPDGQIWKPDILLYNSADERFDATFHTNVLVNSSGHCQYLPPGIFKSSCYIDVRWFPFDVQHCKLKFGSWSYGGWSLDLQMQEADISGYIPNGEWDLVGIPGKRSERFYECCKEPYPDVTFTVTMRRRTLYY",
    "a9nAChR":  "AETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSS",
    "a10nAChR": "AEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAA",
}

with open("contact_dist.json") as f:
    CONTACT_DIST = json.load(f)  # {"a7nAChR": {"36": 3.84, ...}, "a9nAChR": {...}}

PAIRS = [
    ("ImI",   ["a7nAChR"]), ("ImII",  ["a7nAChR"]),
    ("RgIA",  ["a9nAChR"]), ("RgIA",  ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a9nAChR"]), ("Vc1.1", ["a9nAChR", "a10nAChR"]),
    ("PeIA",  ["a9nAChR", "a10nAChR"]), ("Mr1.1", ["a9nAChR", "a10nAChR"]),
    ("AuIB",  ["a7nAChR"]),
]


def per_residue_embedding(model, alphabet, batch_converter, seq, device):
    _, _, tokens = batch_converter([("x", seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[33], return_contacts=False)
    return out["representations"][33][0][1:len(seq) + 1]


def soft_weighted_pool(reps, contact_dist_map):
    """weight_i = 1/dist_i for contact positions, 0 elsewhere; normalized."""
    L = reps.shape[0]
    weights = torch.zeros(L)
    for pos_str, dist in contact_dist_map.items():
        pos = int(pos_str) - 1
        if 0 <= pos < L:
            weights[pos] = 1.0 / max(dist, 0.5)
    if weights.sum() == 0:
        return reps.mean(0)
    weights = weights / weights.sum()
    return (reps.cpu() * weights.unsqueeze(-1)).sum(0)


def cosine(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-9))


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
        tox_reps = get_reps(toxin_seq)
        tox_whole = tox_reps.mean(0).cpu().numpy()

        soft_embs, whole_embs = [], []
        for s in dict.fromkeys(subunit_names):
            reps = get_reps(SUBUNITS[s])
            whole_embs.append(reps.mean(0).cpu().numpy())
            if s in CONTACT_DIST:
                soft_embs.append(soft_weighted_pool(reps, CONTACT_DIST[s]).numpy())
            else:
                soft_embs.append(reps.mean(0).cpu().numpy())
        rec_soft = np.mean(soft_embs, axis=0)
        rec_whole = np.mean(whole_embs, axis=0)

        # Variant A: cosine similarity, toxin vs distance-weighted receptor contact embedding
        cos_soft = cosine(tox_whole, rec_soft)
        cos_whole = cosine(tox_whole, rec_whole)  # baseline: cosine with naive whole-seq too

        # Variant B: distance-weighted contact pooling -> CatBoost
        X = np.concatenate([tox_whole, rec_soft]).reshape(1, -1)
        pred_pKd_distw = -float(cb.predict(X)[0])

        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        print(f"{pair_id:30s} cos_soft={cos_soft:+.4f}  cos_whole={cos_whole:+.4f}  "
              f"distw_catboost_pKd={pred_pKd_distw:.3f}")
        rows.append({
            "pair_id": pair_id, "toxin": toxin_name,
            "receptor_subunits": "+".join(subunit_names),
            "cosine_contact_weighted": round(cos_soft, 4),
            "cosine_whole_seq": round(cos_whole, 4),
            "catboost_dist_weighted_pKd": round(pred_pKd_distw, 4),
        })

    with open("bindpred_variants2_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> bindpred_variants2_predictions.csv")


if __name__ == "__main__":
    main()
