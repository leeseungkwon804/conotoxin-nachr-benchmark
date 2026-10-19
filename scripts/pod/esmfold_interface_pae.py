#!/usr/bin/env python3
"""
Extract ESMFold's interface PAE (predicted aligned error) specifically between the
toxin and the real anchor-derived contact-region receptor residues, as a genuinely
different signal from the ESM2-embedding family (BindPred variants) already
exhausted. Uses the same poly-glycine-linker pseudo-multimer approach validated
earlier (esmfold_inference.py).

interface_pae = mean PAE between toxin positions and contact-region receptor
positions specifically (not the whole receptor, not the whole complex) -- lower
PAE = higher confidence in their relative position, the closest available proxy to
a true interface confidence without native multimer attention.
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
}

with open("contact_dist.json") as f:
    CONTACT_DIST = json.load(f)
CONTACT_POS = {k: [int(p) for p in v.keys()] for k, v in CONTACT_DIST.items()}

PAIRS = [
    ("ImI",   ["a7nAChR"]), ("ImII",  ["a7nAChR"]),
    ("RgIA",  ["a9nAChR"]), ("RgIA",  ["a9nAChR", "a10nAChR"]),
    ("Vc1.1", ["a9nAChR"]), ("Vc1.1", ["a9nAChR", "a10nAChR"]),
    ("PeIA",  ["a9nAChR", "a10nAChR"]), ("Mr1.1", ["a9nAChR", "a10nAChR"]),
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
        # only score against the FIRST subunit that has contact data (a7 or a9 side)
        anchor_subunit = next((s for s in subunit_names if s in CONTACT_POS), None)
        receptor_seq = "".join(SUBUNITS[s] for s in dict.fromkeys(subunit_names))
        complex_seq = toxin_seq + LINKER + receptor_seq

        tokenized = tokenizer([complex_seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.to(device) for k, v in tokenized.items()}
        with torch.no_grad():
            out = model(**tokenized)

        pae = out["predicted_aligned_error"][0]  # (L, L)
        n_tox = len(toxin_seq)
        n_link = len(LINKER)

        # offset of the anchor subunit within the concatenated receptor block
        offset = n_tox + n_link
        for s in dict.fromkeys(subunit_names):
            if s == anchor_subunit:
                break
            offset += len(SUBUNITS[s])

        if anchor_subunit:
            contact_idx = [offset + p - 1 for p in CONTACT_POS[anchor_subunit]]
            tox_idx = list(range(n_tox))
            # interface PAE: toxin residues <-> contact-region receptor residues, both directions
            sub_pae = pae[tox_idx][:, contact_idx]
            sub_pae2 = pae[contact_idx][:, tox_idx]
            interface_pae = float(((sub_pae.mean() + sub_pae2.mean()) / 2).item())
        else:
            interface_pae = None

        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        print(f"{pair_id:30s} interface_PAE={interface_pae:.3f}" if interface_pae is not None
              else f"{pair_id:30s} interface_PAE=N/A")
        rows.append({"pair_id": pair_id, "toxin": toxin_name,
                     "receptor_subunits": "+".join(subunit_names),
                     "interface_PAE": round(interface_pae, 4) if interface_pae is not None else ""})

    with open("esmfold_interface_pae_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> esmfold_interface_pae_predictions.csv")


if __name__ == "__main__":
    main()
