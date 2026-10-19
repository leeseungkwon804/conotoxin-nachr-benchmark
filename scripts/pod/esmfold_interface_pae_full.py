#!/usr/bin/env python3
"""
Extend interface-PAE scoring to all 17 pairs: alpha7/alpha9 pairs use real anchor
contacts (as before); alpha10/alpha3/alpha6/beta2/beta3/beta4 pairs use the
homology-transferred contact positions (from homology_transfer_contacts.py,
alpha7 -> each subunit via BLOSUM62 pairwise alignment, ~40-44% identity).
For heteromeric pairs, interface PAE is computed against BOTH subunits'
contact/transferred positions combined.
"""

import csv
import json
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

TOXINS = {
    "ImI":   "GCCSDPRCAWRC", "ImII":  "ACCSDRRCRWRC", "RgIA":  "GCCSDPRCRYRCR",
    "Vc1.1": "GCCSDPRCNYDHPEIC", "PeIA":  "GCCSHPACSVNHPELC",
    "Mr1.1": "GCCSHPACSVNNPDIC", "AuIB":  "GCCSYPPCFATNPDC",
}
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

with open("contact_dist.json") as f:
    CONTACT_DIST = json.load(f)
ANCHOR_CONTACT_POS = {k: [int(p) for p in v.keys()] for k, v in CONTACT_DIST.items()}

with open("homology_contact_positions.json") as f:
    HOMOLOGY = json.load(f)
HOMOLOGY_CONTACT_POS = {k: v["positions"] for k, v in HOMOLOGY.items()}

ALL_CONTACT_POS = {**ANCHOR_CONTACT_POS, **HOMOLOGY_CONTACT_POS}

PAIRS = [
    ("ImI",   ["a7nAChR"]), ("ImII",  ["a7nAChR"]),
    ("RgIA",  ["a9nAChR"]), ("RgIA",  ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a9nAChR"]), ("Vc1.1", ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a3nAChR", "b2nAChR"]), ("Vc1.1", ["a3nAChR", "b4nAChR"]),
    ("Vc1.1", ["a6nAChR", "b3nAChR"]), ("Vc1.1", ["a6nAChR", "b4nAChR"]),
    ("PeIA",  ["a9nAChR", "a10nAChR"]), ("PeIA",  ["a3nAChR", "b2nAChR"]),
    ("Mr1.1", ["a9nAChR", "a10nAChR"]), ("Mr1.1", ["a3nAChR", "b2nAChR"]),
    ("AuIB",  ["a3nAChR", "b4nAChR"]), ("AuIB",  ["a3nAChR", "b2nAChR"]),
    ("AuIB",  ["a7nAChR"]),
]
LINKER = "G" * 25


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
    model = model.to(device)
    model.esm = model.esm.half()
    model.eval()

    rows = []
    for toxin_name, subunit_names in PAIRS:
        toxin_seq = TOXINS[toxin_name]
        receptor_seq = "".join(SUBUNITS[s] for s in dict.fromkeys(subunit_names))
        complex_seq = toxin_seq + LINKER + receptor_seq

        tokenized = tokenizer([complex_seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.to(device) for k, v in tokenized.items()}
        with torch.no_grad():
            out = model(**tokenized)
        pae = out["predicted_aligned_error"][0]
        n_tox = len(toxin_seq)
        n_link = len(LINKER)
        tox_idx = list(range(n_tox))

        # gather contact indices across all subunits in this pair (real anchor or homology-transferred)
        offset = n_tox + n_link
        contact_idx = []
        source_tags = []
        for s in dict.fromkeys(subunit_names):
            src = "anchor" if s in ANCHOR_CONTACT_POS else "homology"
            positions = ALL_CONTACT_POS.get(s, [])
            contact_idx += [offset + p - 1 for p in positions]
            source_tags.append(f"{s}:{src}({len(positions)})")
            offset += len(SUBUNITS[s])

        if contact_idx:
            sub_pae = pae[tox_idx][:, contact_idx]
            sub_pae2 = pae[contact_idx][:, tox_idx]
            interface_pae = float(((sub_pae.mean() + sub_pae2.mean()) / 2).item())
        else:
            interface_pae = None

        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        tag_str = "+".join(source_tags)
        print(f"{pair_id:32s} interface_PAE={interface_pae:.3f}  [{tag_str}]")
        rows.append({"pair_id": pair_id, "toxin": toxin_name,
                     "receptor_subunits": "+".join(subunit_names),
                     "interface_PAE": round(interface_pae, 4) if interface_pae is not None else "",
                     "contact_source": tag_str})

    with open("esmfold_interface_pae_full_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> esmfold_interface_pae_full_predictions.csv")


if __name__ == "__main__":
    main()
