"""Shared CRISPR defaults used by the local design tools."""


E_COLI_PCRISPR_DEFAULTS = {
    "vector_key": "pcrispr",
    "assembly_method": "golden_gate",
    "enzyme": "BsaI",
    "top_overhang": "AAAC",
    "bottom_overhang": "AAAAC",
    "u6_requires_5prime_g": False,
    "cell_strain": "HME63 or MG1655 carrying pCas9",
    "selection": "Kan",
    "notes": (
        "E. coli pCRISPR guide-array plasmid used with a separate pCas9 "
        "plasmid carrying tracrRNA and Cas9. Jiang et al. used this "
        "two-plasmid system with Lambda Red recombineering in HME63 for "
        "efficient genome editing."
    ),
    "source": (
        "Jiang et al. Nat Biotechnol 2013, doi:10.1038/nbt.2508; "
        "pCas9 carries tracrRNA/Cas9 and pCRISPR carries the spacer array. "
        "Supplementary Fig. 9 and Supplementary Table 2 describe BsaI "
        "insertion of spacers into pCRISPR using AAAC/AAAAC-style oligos."
    ),
}


DEFAULT_VECTOR_BY_REFERENCE = {
    "ecoli_rpsl": E_COLI_PCRISPR_DEFAULTS["vector_key"],
}
