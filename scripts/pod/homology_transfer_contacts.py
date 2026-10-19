#!/usr/bin/env python3
"""
Homology-transfer the alpha7 contact positions (33 residues, from two independent
real complexes 9NX0+9NX1 -- the most robust contact set we have) onto the six
nAChR subunits with NO direct anchor structure (alpha10, alpha3, alpha6, beta2,
beta3, beta4), via pairwise sequence alignment (BLOSUM62, Biopython PairwiseAligner
-- available on this Linux pod; broken on the local Mac's expat, hence deferred
here).

STATED APPROXIMATION, disclosed not hidden: this assumes the toxin-facing surface
is structurally homologous across all Cys-loop nAChR subunits (a reasonable but
unverified assumption for this project -- alpha and beta subunits share the same
fold but the specific loop-C region that dominates principal-face pharmacology is
the most sequence-variable part, exactly where this approximation is weakest).
Positions that fall in a gap of the pairwise alignment (no homologous residue)
are dropped and logged, not silently guessed.
"""

import json
from Bio.Align import PairwiseAligner, substitution_matrices

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

A7_CONTACT_POS = [36, 46, 47, 55, 57, 76, 77, 93, 109, 111, 112, 115, 117, 119, 145,
                  147, 148, 149, 150, 151, 152, 164, 166, 186, 188, 189, 190, 191,
                  192, 193, 195, 196, 197]

TARGET_SUBUNITS = ["a10nAChR", "a3nAChR", "a6nAChR", "b2nAChR", "b3nAChR", "b4nAChR"]


def align_and_transfer(ref_seq, ref_positions, target_seq):
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "global"

    alignment = aligner.align(ref_seq, target_seq)[0]
    ref_aln, tgt_aln = alignment[0], alignment[1]

    # build ref_pos (1-indexed, ungapped) -> alignment_column index
    ref_pos_to_col = {}
    ref_counter = 0
    for col, ch in enumerate(ref_aln):
        if ch != "-":
            ref_counter += 1
            ref_pos_to_col[ref_counter] = col

    # build alignment_column -> target_pos (1-indexed, ungapped)
    col_to_tgt_pos = {}
    tgt_counter = 0
    for col, ch in enumerate(tgt_aln):
        if ch != "-":
            tgt_counter += 1
            col_to_tgt_pos[col] = tgt_counter

    transferred, dropped = [], []
    for p in ref_positions:
        col = ref_pos_to_col.get(p)
        if col is None:
            dropped.append(p)
            continue
        tgt_pos = col_to_tgt_pos.get(col)
        if tgt_pos is None:
            dropped.append(p)  # aligned to a gap in target
        else:
            transferred.append(tgt_pos)

    identity = sum(1 for a, b in zip(ref_aln, tgt_aln) if a == b and a != "-") / \
               sum(1 for a in ref_aln if a != "-")
    return sorted(set(transferred)), dropped, identity


def main():
    a7_seq = SUBUNITS["a7nAChR"]
    results = {}
    for target in TARGET_SUBUNITS:
        transferred, dropped, identity = align_and_transfer(a7_seq, A7_CONTACT_POS, SUBUNITS[target])
        results[target] = {"positions": transferred, "dropped_from_a7": dropped, "seq_identity_to_a7": round(identity, 3)}
        print(f"{target:12s} identity_to_a7={identity:.1%}  transferred={len(transferred)}/{len(A7_CONTACT_POS)}  dropped={dropped}")
        # sanity: print residues at transferred positions
        residues = [SUBUNITS[target][p-1] for p in transferred]
        print(f"             residues: {residues}")

    with open("homology_contact_positions.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\n-> homology_contact_positions.json")


if __name__ == "__main__":
    main()
