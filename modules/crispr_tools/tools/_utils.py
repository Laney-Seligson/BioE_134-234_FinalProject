"""Shared constants and utilities for crispr_tools tools.

This file contains constants (like the codon table) that multiple tools use.
Students can add additional shared utilities here.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Organism normalization
# ---------------------------------------------------------------------------

ORGANISM_ALIASES: dict[str, str] = {
    # E. coli variants
    "ecoli":                    "Escherichia coli",
    "e.coli":                   "Escherichia coli",
    "e. coli":                  "Escherichia coli",
    "escherichia coli":         "Escherichia coli",
    "bacteria":                 "Escherichia coli",
    # Human
    "human":                    "Homo sapiens",
    "homo sapiens":             "Homo sapiens",
    # Mouse
    "mouse":                    "Mus musculus",
    "mus musculus":             "Mus musculus",
    # Yeast
    "yeast":                    "Saccharomyces cerevisiae",
    "s. cerevisiae":            "Saccharomyces cerevisiae",
    "saccharomyces cerevisiae": "Saccharomyces cerevisiae",
    # C. elegans
    "worm":                     "Caenorhabditis elegans",
    "c. elegans":               "Caenorhabditis elegans",
    "caenorhabditis elegans":   "Caenorhabditis elegans",
    # Fly
    "fly":                      "Drosophila melanogaster",
    "drosophila":               "Drosophila melanogaster",
    "drosophila melanogaster":  "Drosophila melanogaster",
    # Zebrafish
    "zebrafish":                "Danio rerio",
    "danio rerio":              "Danio rerio",
    # Arabidopsis
    "arabidopsis":              "Arabidopsis thaliana",
    "arabidopsis thaliana":     "Arabidopsis thaliana",
    # Rat
    "rat":                      "Rattus norvegicus",
    "rattus norvegicus":        "Rattus norvegicus",
    # Plants
    "tobacco":                  "Nicotiana tabacum",
    "nicotiana tabacum":        "Nicotiana tabacum",
    "rice":                     "Oryza sativa",
    "oryza sativa":             "Oryza sativa",
    "corn":                     "Zea mays",
    "maize":                    "Zea mays",
    "zea mays":                 "Zea mays",
}


def normalize_organism(name: str) -> str:
    """Return the canonical scientific name for a common organism alias.

    Unrecognised names are returned unchanged so callers never silently drop
    a valid organism the user provided.
    """
    return ORGANISM_ALIASES.get(name.strip().lower(), name.strip())


# ---------------------------------------------------------------------------
# DNA sequence utilities
# ---------------------------------------------------------------------------

VALID_DNA: set[str] = set("ATGC")

_COMPLEMENT: dict[str, str] = {"A": "T", "T": "A", "G": "C", "C": "G"}


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(_COMPLEMENT[b] for b in reversed(seq.upper()))


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


def _score_guide(protospacer: str, pam_site: str = "") -> float:
    """
    Score a guide protospacer for on-target quality.

    Criteria and sources:
    - TTTT disqualifier: Doench et al. 2016 Nat Biotechnol doi:10.1038/nbt.3437
      explicitly excluded guides with 4+ consecutive T's (Pol III termination);
      CRISPRdirect (Naito et al. 2015 Bioinformatics doi:10.1093/bioinformatics/btu743)
      flags the same motif.
    - GC content 40–60%: practical guideline from Doench 2016 (GC count is a
      model feature; <10 or >10 of 20 nt are binary flags) and widely adopted
      in CRISPOR / CRISPRdirect documentation.
    - Homopolymer penalty: implied by the TTTT reasoning; runs ≥5 of any base
      may impair transcription or folding.
    - Cas12a TTTV PAM bonus: TTTA > TTTC > TTTG ordering from
      Zetsche et al. 2015 Cell doi:10.1016/j.cell.2015.09.038 and
      Kim et al. 2017 Nat Biotechnol doi:10.1038/nbt.3803.
      No equivalent bonus for SpCas9 NGG: within-canonical-PAM N-base effects
      are encoded as complex learned dinucleotide weights in Rule Set 2 and
      are too small to justify a simple heuristic ordering.
    """
    p = protospacer.upper()
    n = len(p)

    # Priority 1: TTTT run terminates Pol III transcription — hard disqualifier
    if "TTTT" in p:
        return 0.0

    score = 100.0

    # Priority 2: GC content; 40–60% optimal, penalise outside that window
    gc = (p.count("G") + p.count("C")) / n
    if gc < 0.40:
        score -= (0.40 - gc) * 250   # −25 per 10% below 40%
    elif gc > 0.60:
        score -= (gc - 0.60) * 250   # −25 per 10% above 60%

    # Homopolymer run of 5+ identical bases
    if any(b * 5 in p for b in "ACGT"):
        score -= 20

    # Cas12a TTTV PAM tie-breaker (Zetsche 2015, Kim 2017)
    pam = pam_site.upper()
    if len(pam) == 4 and pam.startswith("TTT"):
        score += {"A": 5, "C": 3, "G": 1}.get(pam[3], 0)

    return round(max(0.0, score), 1)


def _guide_flags(protospacer: str) -> str:
    p = protospacer.upper()
    flags = []
    if "TTTT" in p:
        flags.append("TTTT!")
    gc = (p.count("G") + p.count("C")) / len(p)
    if gc < 0.40:
        flags.append("low-GC")
    elif gc > 0.60:
        flags.append("high-GC")
    if any(b * 5 in p for b in "ACGT"):
        flags.append("homopolymer")
    return ",".join(flags)


def rank_guides(guides: list) -> list:
    """Score and return guides sorted best-first."""
    ranked = []

    for g in guides:
        p = g["protospacer"]
        score = _score_guide(p, g.get("pam_site", ""))
        gc_pct = round((p.count("G") + p.count("C")) / len(p) * 100, 1)

        ranked.append(g | {
            "score": score,
            "gc_percent": gc_pct,
            "flags": _guide_flags(p),
        })

    return sorted(ranked, key=lambda g: g["score"], reverse=True)
