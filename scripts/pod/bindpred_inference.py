#!/usr/bin/env python3
"""
BindPred zero-shot inference for the conotoxin nAChR-axis benchmark.

Reproduces the exact feature construction from hbp5181/BindPred's train.py
(fetched live from HuggingFace 2026-07-22): mean-pooled ESM2 (esm2_t33_650M_UR50D,
final layer) embedding per protein, concatenated as [ligand_features, receptor_features]
(1280+1280=2560-dim), fed to the pretrained CatBoostRegressor (ESM2_BindPred.cbm),
which outputs predicted log10(Kd in M).

STATED APPROXIMATION: ESM2 is a single-sequence model with no multimer/cross-chain
attention. For heteromeric receptor pairs (e.g. alpha9+alpha10), the two subunit
sequences are concatenated into one string before embedding (common practical
workaround for applying single-sequence PLMs to multi-chain complexes). For
homomeric receptors (e.g. alpha7+alpha7), only one copy is embedded (redundant to
embed the identical sequence twice). Toxin = "ligand", channel = "receptor",
matching the paper's smaller/larger molecule convention.
"""

import csv
import json
import torch
import esm
from catboost import CatBoostRegressor

MODEL_PATH = "/root/.cache/huggingface/hub/models--hbp5181--BindPred/snapshots/89a595d1c161a95783d904f5c5b3cce565e52848/ESM2_BindPred.cbm"

TOXINS = {
    "ImI":   "GCCSDPRCAWRC",
    "ImII":  "ACCSDRRCRWRC",
    "RgIA":  "GCCSDPRCRYRCR",
    "Vc1.1": "GCCSDPRCNYDHPEIC",
    "PeIA":  "GCCSHPACSVNHPELC",
    "Mr1.1": "GCCSHPACSVNNPDIC",
    "AuIB":  "GCCSYPPCFATNPDC",
}

# ECD-cropped subunit sequences (from data/human_channels/ecd_cropped/, script 12)
SUBUNITS = {
    "a7nAChR":  "GEFQRKLYKELVKNYNPLERPVANDSQPLTVYFSLSLLQIMDVDEKNQVLTTNIWLQMSWTDHYLQWNVSEYPGVKTVRFPDGQIWKPDILLYNSADERFDATFHTNVLVNSSGHCQYLPPGIFKSSCYIDVRWFPFDVQHCKLKFGSWSYGGWSLDLQMQEADISGYIPNGEWDLVGIPGKRSERFYECCKEPYPDVTFTVTMRRRTLYY",
    "a9nAChR":  "AETADGKYAQKLFNDLFEDYSNALRPVEDTDKVLNVTLQITLSQIKDMDERNQILTAYLWIRQIWHDAYLTWDRDQYDGLDSIRIPSDLVWRPDIVLYNKADDESSEPVNTNVVLRYDGLITWDAPAITKSSCVVDVTYFPFDNQQCNLTFGSWTYNGNQVDIFNALDSGDLSDFIEDVEWEVHGMPAVKNVISYGCCSEPYPDVTFTLLLKRRSS",
    "a10nAChR": "AEGRLALKLFRDLFANYTSALRPVADTDQTLNVTLEVTLSQIIDMDERNQVLTLYLWIRQEWTDAYLRWDPNAYGGLDAIRIPSSLVWRPDIVLYNKADAQPPGSASTNVVLRHDGAVRWDAPAITRSSCRVDVAAFPFDAQHCGLTFGSWTHGGHQLDVRPRGAAASLADFVENVEWRVLGMPARRRVLTYGCCSEPYPDVTFTLLLRRRAA",
    "a3nAChR":  "SEAEHRLFERLFEDYNEIIRPVANVSDPVIIHFEVSMSQLVKVDEVNQIMETNLWLKQIWNDYKLKWNPSDYGGAEFMRVPAQKIWKPDIVLYNNAVGDFQVDDKTKALLKYTGEVTWIPPAIFKSSCKIDVTYFPFDYQNCTMKFGSWSYDKAKIDLVLIGSSMNLKDYWESGEWAIIKAPGYKHDIKYNCCEEIYPDITYSLYIRRLPLFYTINLII",
    "b2nAChR":  "TDTEERLVEHLLDPSRYNKLIRPATNGSELVTVQLMVSLAQLISVHEREQIMTTNVWLTQEWEDYRLTWKPEEFDNMKKVRLPSKHIWLPDVVLYNNADGMYEVSFYSNAVVSYDGSIFWLPPAIYKSACKIEVKHFPFDQQNCTMKFRSWTYDRTEIDLVLKSEVASLDDFTPSGEWDIVALPGRRNENPDDSTYVDITYDFIIRRKPLFYT",
    "b4nAChR":  "RVANAEEKLMDDLLNKTRYNNLIRPATSSSQLISIKLQLSLAQLISVNEREQIMTTNVWLKQEWTDYRLTWNSSRYEGVNILRIPAKRIWLPDIVLYNNADGTYEVSVYTNLIVRSNGSVLWLPPAIYKSACKIEVKYFPFDQQNCTLKFRSWTYDHTEIDMVLMTPTASMDDFTPSGEWDIVALPGRRTVNPQDPSYVDVTYDFIIKRKPLFYTINL",
    "a6nAChR":  "KGCVGCATEERLFHKLFSHYNQFIRPVENVSDPVTVHFEVAITQLANVDEVNQIMETNLWLRHIWNDYKLRWDPMEYDGIETLRVPADKIWKPDIVLYNNAVGDFQVEGKTKALLKYNGMITWTPPAIFKSSCPMDITFFPFDHQNCSLKFGSWTYDKAEIDLLIIGSKVDMNDFWENSEWEIIDASGYKHDIKYNCCEEIYTDITYSFYIRRL",
    "b3nAChR":  "FNSIAENEDALLRHLFQGYQKWVRPVLHSNDTIKVYFGLKISQLVDVDEKNQLMTTNVWLKQEWTDHKLRWNPDDYGGIHSIKVPSESLWLPDIVLFENADGRFEGSLMTKVIVKSNGTVVWTPPASYKSSCTMDVTFFPFDRQNCSMKFGSWTYDGTMVDLILINENVDRKDFFDNGEWEILNAKGMKGNRRDGVYSYPFITYSFVLRRLPLFYT",
}

PAIRS = [
    ("ImI",   ["a7nAChR"]),
    ("ImII",  ["a7nAChR"]),
    ("RgIA",  ["a9nAChR"]),
    ("RgIA",  ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a9nAChR"]),
    ("Vc1.1", ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a3nAChR", "b2nAChR"]),
    ("Vc1.1", ["a3nAChR", "b4nAChR"]),
    ("Vc1.1", ["a6nAChR", "b3nAChR"]),
    ("Vc1.1", ["a6nAChR", "b4nAChR"]),
    ("PeIA",  ["a9nAChR", "a10nAChR"]),
    ("PeIA",  ["a3nAChR", "b2nAChR"]),
    ("Mr1.1", ["a9nAChR", "a10nAChR"]),
    ("Mr1.1", ["a3nAChR", "b2nAChR"]),
    ("AuIB",  ["a3nAChR", "b4nAChR"]),
    ("AuIB",  ["a3nAChR", "b2nAChR"]),
    ("AuIB",  ["a7nAChR"]),
]


def mean_pool_embedding(model, alphabet, batch_converter, seq, device):
    _, _, tokens = batch_converter([("x", seq)])
    tokens = tokens.to(device)
    with torch.no_grad():
        out = model(tokens, repr_layers=[33], return_contacts=False)
    reps = out["representations"][33][0]  # (L+2, 1280) incl BOS/EOS
    # mean-pool over real residues only (exclude BOS at 0 and EOS at len+1)
    return reps[1:len(seq) + 1].mean(0).cpu().numpy()


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}")

    print("loading ESM2 esm2_t33_650M_UR50D ...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device).eval()
    batch_converter = alphabet.get_batch_converter()

    print("loading BindPred CatBoost model ...")
    cb = CatBoostRegressor()
    cb.load_model(MODEL_PATH, format="cbm")

    # cache embeddings so each unique sequence is only embedded once
    cache = {}
    def get_emb(seq):
        if seq not in cache:
            cache[seq] = mean_pool_embedding(model, alphabet, batch_converter, seq, device)
        return cache[seq]

    rows = []
    for toxin_name, subunit_names in PAIRS:
        toxin_seq = TOXINS[toxin_name]
        # heteromeric -> concatenate subunit sequences; homomeric -> single copy
        receptor_seq = "".join(SUBUNITS[s] for s in dict.fromkeys(subunit_names))

        lig_emb = get_emb(toxin_seq)
        rec_emb = get_emb(receptor_seq)
        import numpy as np
        X = np.concatenate([lig_emb, rec_emb]).reshape(1, -1)
        pred_log10_kd = float(cb.predict(X)[0])
        pred_pKd = -pred_log10_kd

        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        print(f"{pair_id:30s} predicted log10(Kd)={pred_log10_kd:+.3f}  pKd={pred_pKd:.3f}")
        rows.append({
            "pair_id": pair_id, "toxin": toxin_name,
            "receptor_subunits": "+".join(subunit_names),
            "predicted_log10_Kd": round(pred_log10_kd, 4),
            "predicted_pKd": round(pred_pKd, 4),
        })

    with open("bindpred_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> bindpred_predictions.csv")


if __name__ == "__main__":
    main()
