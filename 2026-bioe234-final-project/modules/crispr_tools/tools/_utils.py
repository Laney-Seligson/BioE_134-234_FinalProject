"""Shared constants and utilities for crispr_tools tools.

This file contains constants (like the codon table) that multiple tools use.
Students can add additional shared utilities here.
"""

from __future__ import annotations

# Standard genetic code (DNA codons -> amino acids)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Valid sequence characters (DNA/RNA + IUPAC ambiguity codes)
# A, T, C, G = standard DNA bases
# U = RNA uracil
# R = A or G (purine)
# Y = C or T (pyrimidine)
# S = G or C (strong)
# W = A or T (weak)
# K = G or T (keto)
# M = A or C (amino)
# N = any base
VALID_SEQUENCE_CHARS = set("ATUCGRSYKWMN")


def _score_guide(protospacer: str) -> float:
    gc = (protospacer.count("G") + protospacer.count("C")) / len(protospacer)
    score = 100 - abs(gc - 0.5) * 200   # 100 at 50% GC, 0 at 0% or 100%
    if "TTTT" in protospacer:            score -= 30   # Pol III termination risk
    if any(b * 5 in protospacer for b in "ACGT"): score -= 20   # homopolymer
    return round(score, 1)


def rank_guides(guides: list) -> list:
    """Score, print, and return guides sorted best-first."""
    scored = sorted(((g, _score_guide(g["protospacer"])) for g in guides),
                    key=lambda x: x[1], reverse=True)
    print(f"\n{'Rank':<5} {'Score':<7} {'GC%':<6} {'PAM':<6} Protospacer")
    for rank, (g, score) in enumerate(scored, 1):
        p = g["protospacer"]
        gc_pct = f"{(p.count('G')+p.count('C'))/len(p)*100:.0f}%"
        print(f"  #{rank:<3} {score:<7} {gc_pct:<6} {g.get('pam_site',''):<6} {p}")
    return [g | {"score": s} for g, s in scored]
