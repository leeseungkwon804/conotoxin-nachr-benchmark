#!/usr/bin/env python3
"""
ESMFold zero-training structure-confidence scoring for the conotoxin nAChR-axis
benchmark, as a third method (after Boltz-2 structure-confidence, BindPred
sequence-embedding regression).

Uses HuggingFace's facebook/esmfold_v1 (verified working on this pod 2026-07-24;
the fair-esm package's own esmfold_v1 loads the SAME 3B-parameter model via a
different checkpoint format -- reusing the already-downloaded HF checkpoint avoids
a redundant ~15GB download).

KNOWN LIMITATION, STATED not hidden: ESMFold has no native cross-chain attention
(single-sequence model). Multimer input via a literal ":" chain-break character
FAILS on this HF port's tokenizer (confirmed bug, reproduces a prior failure log
from a separate project on this same pod). Standard workaround (as used in Meta's
own ESMFold multimer examples): join chains with a poly-glycine linker (25 x 'G'),
fold as one sequence, then mask out the linker region before computing per-chain
confidence. This is an approximation of a true multimer prediction, not a substitute
for real cross-chain modeling (same caveat class as Boltz-2's homomeric/heteromeric
minimal-interface simplification, see Structural anchors notes).
"""

import csv
import torch
from transformers import AutoTokenizer, EsmForProteinFolding

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

LINKER = "G" * 25


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"device={device}")
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

        plddt = out["plddt"][0].mean(-1)  # per-residue, mean over atoms
        n_tox = len(toxin_seq)
        n_link = len(LINKER)

        toxin_plddt = plddt[:n_tox].mean().item()
        receptor_plddt = plddt[n_tox + n_link:].mean().item()
        overall_ptm = out["ptm"].item()
        # interface proxy: pLDDT of toxin residues only (how confidently is the
        # SHORT chain folded in the presence of the receptor context) -- closest
        # analogue available to Boltz's ipTM without true multimer attention
        pair_id = f"{toxin_name}__{'_'.join(subunit_names)}"
        print(f"{pair_id:30s} ptm={overall_ptm:.3f}  toxin_pLDDT={toxin_plddt:.3f}  "
              f"receptor_pLDDT={receptor_plddt:.3f}")

        rows.append({
            "pair_id": pair_id, "toxin": toxin_name,
            "receptor_subunits": "+".join(subunit_names),
            "ptm": round(overall_ptm, 4),
            "toxin_pLDDT": round(toxin_plddt, 4),
            "receptor_pLDDT": round(receptor_plddt, 4),
        })

    with open("esmfold_predictions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("\n-> esmfold_predictions.csv")


if __name__ == "__main__":
    main()
