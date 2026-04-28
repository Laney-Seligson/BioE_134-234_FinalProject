"""
Shared citation registry for the verification / scoring tools authored by
Jillian (rank_guides, predict_offtargets, predict_editing_efficiency,
colony_calculator, interpret_ice_tide, verify_edit).

Matches the Citation/format_citations pattern established in
design_cloning_oligos.py (Laney's tool) so the citation shape is
consistent across the whole project:

    {"label": "...", "reference": "...", "claim": "..."}

Each tool imports the citations it needs by key and assembles a
per-call list reflecting which rules actually fired.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Citation:
    """One structured literature citation. Same shape as the one in
    design_cloning_oligos.py to keep output schema consistent."""
    label: str            # e.g. "Doench et al. 2016, Nat Biotechnol"
    url_or_reference: str # DOI or stable URL
    claim: str            # which rule/recommendation this citation supports


def format_citations(citations) -> list[dict]:
    """Render a sequence of Citation objects to the project's standard
    list-of-dicts format. Mirrors design_cloning_oligos._format_citations."""
    return [
        {"label": c.label, "reference": c.url_or_reference, "claim": c.claim}
        for c in citations
    ]


# ---------------------------------------------------------------------------
# Citation registry — keyed by short symbolic name. Tools pick by key.
# ---------------------------------------------------------------------------

CITATIONS: dict[str, Citation] = {
    # On-target efficiency (Doench Rule Set 2)
    "doench_2016": Citation(
        label="Doench et al. 2016, Nat Biotechnol 34:184-191",
        url_or_reference="https://doi.org/10.1038/nbt.3437",
        claim="Rule Set 2: position-specific nucleotide weights, GC content "
              "optimum 40-65%, NGGT > NGGA PAM context, and CFD off-target "
              "scoring (Suppl. Table 19).",
    ),

    # Cas9 mechanism / NGG PAM
    "jinek_2012": Citation(
        label="Jinek et al. 2012, Science 337:816-821",
        url_or_reference="https://doi.org/10.1126/science.1225829",
        claim="SpCas9 cuts 3bp upstream of the NGG PAM producing a blunt "
              "double-strand break.",
    ),

    # Cas9 seed region + off-target
    "hsu_2013": Citation(
        label="Hsu et al. 2013, Nat Biotechnol 31:827-832",
        url_or_reference="https://doi.org/10.1038/nbt.2647",
        claim="SpCas9 seed region is positions 1-12 from the PAM-proximal end; "
              "mismatches in the seed reduce cutting much more than PAM-distal "
              "mismatches. NAG is a weak alternate PAM (~10x lower cutting).",
    ),

    # Cas12a mechanism / TTTV PAM
    "zetsche_2015": Citation(
        label="Zetsche et al. 2015, Cell 163(3):759-771",
        url_or_reference="https://doi.org/10.1016/j.cell.2015.09.038",
        claim="LbCas12a recognizes a TTTV (V=A/C/G) PAM 5' of a 23bp "
              "protospacer and produces a staggered cut with a 5-nt 5' overhang.",
    ),

    # Cas12a seed
    "kim_2016_cas12a": Citation(
        label="Kim et al. 2016, Nat Biotechnol 34:863-868",
        url_or_reference="https://doi.org/10.1038/nbt.3609",
        claim="Cas12a seed region is positions 1-10 from the PAM-proximal "
              "end; Cas12a is more mismatch-intolerant than Cas9 across the "
              "full guide.",
    ),

    # Cas12a position-specific scoring
    "kim_2018_cas12a": Citation(
        label="Kim et al. 2018, Nat Biotechnol 36:239-241",
        url_or_reference="https://doi.org/10.1038/nbt.4061",
        claim="Position-specific nucleotide preferences for Cas12a on-target "
              "activity (DeepCpf1 model).",
    ),

    # HDR efficiency
    "paquet_2016": Citation(
        label="Paquet et al. 2016, Nature 533:125-129",
        url_or_reference="https://doi.org/10.1038/nature17664",
        claim="HDR knockin efficiency is typically ~5-10% in mammalian cells, "
              "~1/10 of NHEJ cutting efficiency.",
    ),

    # Bacterial CRISPR efficiency benchmark
    "jiang_2013": Citation(
        label="Jiang et al. 2013, Nat Biotechnol 31:233-239",
        url_or_reference="https://doi.org/10.1038/nbt.2508",
        claim="CRISPR-Cas9 editing in E. coli typically reaches 30-70% "
              "efficiency with pCRISPR/lambda Red recombineering systems.",
    ),

    # RNP delivery efficiency
    "kim_2014_rnp": Citation(
        label="Kim et al. 2014, Genome Res 24:1012-1019",
        url_or_reference="https://doi.org/10.1101/gr.171322.113",
        claim="Cas9 RNP delivery to mammalian cells achieves ~50-80% editing "
              "efficiency, higher than plasmid transfection.",
    ),

    # Lin RNP
    "lin_2014": Citation(
        label="Lin et al. 2014, eLife 3:e04766",
        url_or_reference="https://doi.org/10.7554/eLife.04766",
        claim="RNP delivery improves editing efficiency and reduces "
              "off-target effects vs. plasmid expression.",
    ),

    # Base editors
    "komor_2016": Citation(
        label="Komor et al. 2016, Nature 533:420-424",
        url_or_reference="https://doi.org/10.1038/nature17946",
        claim="Cytidine base editors achieve ~30-50% editing efficiency at "
              "target Cs without double-strand breaks.",
    ),

    # Prime editors
    "anzalone_2019": Citation(
        label="Anzalone et al. 2019, Nature 576:149-157",
        url_or_reference="https://doi.org/10.1038/s41586-019-1711-4",
        claim="Prime editing achieves ~10-30% editing efficiency for diverse "
              "edit types without requiring HDR.",
    ),

    # ICE
    "hsiau_2019_ice": Citation(
        label="Hsiau et al. 2019, bioRxiv 251082",
        url_or_reference="https://doi.org/10.1101/251082",
        claim="ICE (Synthego) deconvolves Sanger trace mixed peaks into indel "
              "alleles and reports a KO score; R^2 below 0.80 indicates the "
              "fit is unreliable.",
    ),

    # TIDE
    "brinkman_2014_tide": Citation(
        label="Brinkman et al. 2014, Nucleic Acids Res 42(22):e168",
        url_or_reference="https://doi.org/10.1093/nar/gku936",
        claim="TIDE quantifies CRISPR indel frequencies by linear "
              "decomposition of mixed Sanger trace signal downstream of "
              "the cut site.",
    ),

    # Wallace primer Tm
    "wallace_1979": Citation(
        label="Wallace et al. 1979, Nucleic Acids Res 6(11):3543-3557",
        url_or_reference="https://doi.org/10.1093/nar/6.11.3543",
        claim="Primer melting temperature for short oligos approximated by "
              "Tm = 2*(A+T) + 4*(G+C) degC.",
    ),

    # PCR primer design rules
    "dieffenbach_1993": Citation(
        label="Dieffenbach et al. 1993, PCR Methods Appl 3(3):S30-37",
        url_or_reference="https://doi.org/10.1101/gr.3.3.S30",
        claim="PCR primer design guidelines: Tm 50-65 degC, GC 40-60%, "
              "no poly-N runs, primer pair |delta-Tm| < 5 degC.",
    ),

    # Pol III termination
    "polIII_termination": Citation(
        label="Bogenhagen & Brown 1981, Cell 24:261-270; Nielsen et al. 2013",
        url_or_reference="https://doi.org/10.1016/0092-8674(81)90522-5",
        claim="Pol III terminates transcription at runs of 4+ T's on the "
              "non-template strand; TTTT in a U6-driven sgRNA truncates "
              "the guide RNA before it is fully transcribed.",
    ),

    # Fonfara Cas12a cut staggered
    "fonfara_2016": Citation(
        label="Fonfara et al. 2016, Nature 532:517-521",
        url_or_reference="https://doi.org/10.1038/nature17945",
        claim="LbCas12a makes a staggered double-strand break: non-template "
              "strand between positions 18-19, template strand between 23-24.",
    ),

    # ---- Wet-lab protocol sources (used by _protocols.py / lab_sheet) ----

    # Golden Gate one-pot assembly
    "engler_2008_goldengate": Citation(
        label="Engler et al. 2008, PLoS ONE 3(11):e3647",
        url_or_reference="https://doi.org/10.1371/journal.pone.0003647",
        claim="Golden Gate cloning: type IIs restriction enzyme (e.g. BsaI) "
              "and T4 DNA ligase in a single reaction, with thermocycling "
              "between 37C (cut) and 16C (ligate).",
    ),

    # NEB Golden Gate / BsaI-HFv2 product manual
    "neb_e1601": Citation(
        label="NEB Golden Gate Assembly Kit (BsaI-HFv2), product E1601",
        url_or_reference="https://www.neb.com/en-us/products/e1601-neb-golden-gate-assembly-kit-bsai-hfv2",
        claim="NEB recommended Golden Gate reaction: 1 uL T4 ligase buffer, "
              "0.5 uL T4 ligase, 0.5 uL BsaI-HFv2, ~1 uL DNA mix, ddH2O to "
              "10 uL; cycle 37C/16C.",
    ),

    # Gibson Assembly original paper
    "gibson_2009": Citation(
        label="Gibson et al. 2009, Nat Methods 6:343-345",
        url_or_reference="https://doi.org/10.1038/nmeth.1318",
        claim="Gibson Assembly: exonuclease, polymerase, and ligase in a "
              "single isothermal 50C reaction join overlapping DNA fragments.",
    ),

    # NEB Gibson Master Mix manual
    "neb_e2611": Citation(
        label="NEB Gibson Assembly Master Mix, product E2611",
        url_or_reference="https://www.neb.com/en-us/products/e2611-gibson-assembly-master-mix",
        claim="NEB Gibson 2X Master Mix: 5 uL master mix + up to 5 uL DNA "
              "(0.02-0.5 pmol per fragment); incubate 50C for 15-60 min.",
    ),

    # NEB Q5 polymerase manual
    "neb_m0491_q5": Citation(
        label="NEB Q5 High-Fidelity DNA Polymerase, product M0491",
        url_or_reference="https://www.neb.com/en-us/products/m0491-q5-high-fidelity-dna-polymerase",
        claim="Q5 PCR (50 uL): 10 uL 5X Q5 buffer, 1 uL 10 mM dNTPs, "
              "2.5 uL each 10 uM primer, template, 0.5 uL Q5 polymerase, "
              "ddH2O to 50 uL. Extension 30 sec/kb at 72C.",
    ),

    # Chemical transformation - Inoue method
    "inoue_1990": Citation(
        label="Inoue et al. 1990, Gene 96(1):23-28",
        url_or_reference="https://doi.org/10.1016/0378-1119(90)90336-p",
        claim="Chemically competent E. coli prepared with the Inoue method "
              "achieve >1e8 cfu/ug; standard heat-shock protocol: 30 min on "
              "ice, 42C 30-45 sec, 2 min on ice, recover in SOC at 37C.",
    ),

    # Sambrook canonical reference
    "sambrook_2001": Citation(
        label="Sambrook & Russell 2001, Molecular Cloning: A Laboratory Manual, 3rd ed.",
        url_or_reference="https://www.cshlpress.com/default.tpl?action=full&--eqskudatarq=525",
        claim="Standard transformation, plating, and antibiotic selection "
              "protocols; rescue/recovery in SOC or 2YT for 45-60 min before "
              "plating on selection.",
    ),

    # Zymo miniprep
    "zymo_d4015": Citation(
        label="Zymo Research Plasmid Miniprep-Classic, product D4015",
        url_or_reference="https://www.zymoresearch.com/products/zyppy-plasmid-miniprep-kit",
        claim="Silica-column plasmid miniprep: lyse 600 uL overnight culture, "
              "neutralize, bind to column, wash, elute in 30-50 uL elution "
              "buffer or ddH2O.",
    ),

    # Sanger sequencing premix (UC Berkeley DNA Sequencing Facility, standard)
    "sanger_premix_standard": Citation(
        label="UC Berkeley DNA Sequencing Facility submission guidelines",
        url_or_reference="https://mcb.berkeley.edu/barker/dnaseq/",
        claim="Standard Sanger premix submission: ~10 uL final volume "
              "containing 200-500 ng plasmid DNA and 3.2 pmol primer "
              "(typical: 4 uL miniprep DNA + 3 uL of 2.66 uM primer + ddH2O).",
    ),

    # IDT Alt-R Cas9 RNP assembly (vendor protocol; widely used for
    # mammalian electroporation/lipofection experiments).
    "idt_altr_rnp": Citation(
        label="IDT Alt-R CRISPR-Cas9 System: RNP assembly protocol",
        url_or_reference="https://www.idtdna.com/pages/products/crispr-genome-editing/alt-r-crispr-cas9-system",
        claim="Anneal Alt-R crRNA + tracrRNA at equimolar (1 uL each, 100 uM "
              "stock) in 1 uL Duplex Buffer; complex with Cas9 protein at 1:1.2 "
              "molar ratio (RNA:protein) at room temp 20 min before delivery.",
    ),

    # Lonza/Neon electroporation (foundational Lonza Nucleofector protocol).
    "lonza_nucleofection": Citation(
        label="Lonza Nucleofector (4D) CRISPR/Cas9 RNP delivery protocol",
        url_or_reference="https://bioscience.lonza.com/lonza_bs/CH/en/Cell-and-Transfection",
        claim="20 uL Nucleocuvette: 100,000-500,000 cells in P3/SF/SE solution, "
              "add 5 uL pre-formed RNP complex, run program (e.g. CA-137 for "
              "K562, CM-138 for HEK293, EH-100 for primary T cells); transfer "
              "to pre-warmed media within 10 min.",
    ),

    # Lipofection of CRISPRMAX (ThermoFisher) for plasmid/RNP transfection.
    "thermo_crisprmax": Citation(
        label="ThermoFisher Lipofectamine CRISPRMAX user guide",
        url_or_reference="https://www.thermofisher.com/order/catalog/product/CMAX00001",
        claim="24-well lipofection: dilute 7.5 pmol Cas9 RNP + 1 uL Cas9 Plus "
              "Reagent in 25 uL Opti-MEM; in parallel dilute 1.5 uL CRISPRMAX "
              "in 25 uL Opti-MEM; combine 1:1, incubate 5 min RT, add 50 uL "
              "complex to cells in 500 uL antibiotic-free media.",
    ),

    # NEB EnGen sgRNA Synthesis Kit / IVT (E3322).
    "neb_e3322_ivt": Citation(
        label="NEB EnGen sgRNA Synthesis Kit (E3322)",
        url_or_reference="https://www.neb.com/en-us/products/e3322-engen-sgrna-synthesis-kit-s-pyogenes",
        claim="One-pot in vitro transcription of sgRNA from a target-specific "
              "DNA oligo: 1 uL target oligo (1 uM) + 2 uL sgRNA scaffold mix + "
              "10 uL NTP buffer mix + 2 uL T7 RNA polymerase mix + 5 uL ddH2O; "
              "incubate 37C 30 min; DNase I treat 15 min; column or LiCl purify.",
    ),
}


def cites(*keys: str) -> list[Citation]:
    """Convenience: pull a list of Citation objects from CITATIONS by key.

    Raises KeyError if a key is unknown — fail fast rather than silently
    drop a citation and ship a paper-less recommendation.
    """
    return [CITATIONS[k] for k in keys]
