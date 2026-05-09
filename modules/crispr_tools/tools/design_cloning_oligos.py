"""
CRISPR Cloning Workflow Designer
=================================
A guided wet-lab simulation tool for four CRISPR guide-cloning workflows:

  1. TypeIISOligoCloning  — BbsI/BsmBI/BsaI annealed-oligo sticky-end ligation
  2. RestrictionLigation  — conventional RE cloning (e.g. SpeI for pTargetF)
  3. GibsonAssembly       — overlap-based seamless assembly
  4. GoldenGateAssembly   — Type IIS PCR-based multi-fragment assembly with
                            user-designed 4-nt overhangs; emits all four primer
                            fields required by create_construction_file's
                            assembly_strategy="GoldenGate" mode.

Behavior
--------
- If all required inputs are present → design the oligos/primers and return results.
- If inputs are missing → return a "needs_user_input" dict with human-readable
  questions instead of raising an error.
- Hard errors are reserved for invalid DNA strings and genuinely broken biology.

construction_file_inputs compatibility
---------------------------------------
TypeIISOligoCloning, RestrictionLigation, and GibsonAssembly emit
assembly_strategy="DirectSynthesis" because create_construction_file's Gibson
and GoldenGate strategies require all four PCR-primer fields including vector
linearisation primers those branches do not design.

GoldenGateAssembly is the exception: it designs all four primers (insert
forward/reverse + vector forward/reverse) and therefore emits
assembly_strategy="GoldenGate", which is fully compatible with
create_construction_file's GoldenGate step.

DISCLAIMER: For educational planning only. Validate designs against the actual
plasmid map before ordering reagents.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from modules.crispr_tools.tools.fetch_addgene_vector import FetchAddgeneVector as _AddgeneFetcher
from modules.crispr_tools.tools._utils import normalize_organism, reverse_complement as _reverse_complement


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

TOOL_DIR        = Path(__file__).resolve().parent
CRISPR_TOOLS_DIR = TOOL_DIR.parent
BUNDLED_DATA_DIR = CRISPR_TOOLS_DIR / "data"


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Citation:
    """One structured literature or Addgene citation."""
    label: str            # e.g. "Cong et al. Science 2013"
    url_or_reference: str # DOI, Addgene URL, or ISBN
    claim: str            # what this citation supports


@dataclass(frozen=True)
class VectorSpec:
    """
    Describes one cloning vector preset.

    cloning_method controls which branch .run() dispatches to:
      "TypeIISOligoCloning"  — annealed-oligo guide insertion via Type IIS RE
      "RestrictionLigation"  — conventional RE cloning with compatible sticky ends
      "GibsonAssembly"       — overlap-based seamless assembly

    dna_source is the physical form of the insert in a real experiment:
      "annealed_oligos"      — two synthetic oligos annealed together
      "pcr_product"          — PCR-amplified with RE-site-tailed primers
      "overlap_fragment"     — synthesised fragment carrying designed overlaps
    """
    name: str
    cloning_method: str        # see docstring
    dna_source: str
    enzyme: str
    promoter: str
    nuclease_system: str
    guide_type: str
    scaffold_in_vector: bool

    # TypeIIS-specific (empty string when not applicable)
    recognition_site: str = ""
    top_overhang: str = ""
    bottom_overhang: str = ""
    # "prefers" not "requires" — U6 efficiency is higher with 5′G but the
    # guide still works without it.
    u6_prefers_5prime_g: bool = False

    # RestrictionLigation-specific
    restriction_site_sequence: str = ""  # RE recognition seq added to PCR primer tails

    # GibsonAssembly-specific
    recommended_overlap_bp: int = 20
    gibson_left_context: str = ""   # known vector sequence 5′ of guide insertion site
    gibson_right_context: str = ""  # known vector sequence 3′ of guide insertion site

    # Resource / strain metadata
    backbone_resource: Optional[str] = None  # key in BACKBONE_RESOURCES
    cell_strain: str = ""
    selection: str = ""
    notes: str = ""
    citations: tuple[Citation, ...] = field(default_factory=tuple)


# ---------------------------------------------------------------------------
# File-system registries
# ---------------------------------------------------------------------------

BACKBONE_RESOURCES: dict[str, Path] = {
    # pX330 sequence source: NovoPro pX330-U6-Chimeric_BB-CBh-hSpCas9
    # GenBank Map (.gb), 8484 bp, Cat. No. V005940.
    # Plasmid identity and wet-lab metadata are cross-checked against
    # Addgene plasmid #42230 (pX330-U6-Chimeric_BB-CBh-hSpCas9).
    "px330":         BUNDLED_DATA_DIR / "pX330-U6-Chimeric_BB-CBh-hSpCas9.gb",
    "pET28a":        BUNDLED_DATA_DIR / "pET28a.gb",
    "pBR322":        BUNDLED_DATA_DIR / "pBR322.gb",
    "pCRISPR_rpsL":  BUNDLED_DATA_DIR / "pCRISPR_rpsL.gbk",
}


# ---------------------------------------------------------------------------
# Vector preset registry
# ---------------------------------------------------------------------------

VECTOR_SPECS: dict[str, VectorSpec] = {

    # ── TypeIIS / Annealed-Oligo Cloning ─────────────────────────────────────
    # The vector backbone encodes a Type IIS enzyme site adjacent to the guide
    # scaffold.  Digest leaves 4-nt 5′ overhangs whose sequence is vector-specific.
    # Two synthetic oligos with matching overhangs are annealed and ligated in.
    # No PCR of the insert is required.

    "px330": VectorSpec(
        name="pX330",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BbsI",
        promoter="U6",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="GAAGAC",
        top_overhang="CACC",
        bottom_overhang="AAAC",
        u6_prefers_5prime_g=True,
        backbone_resource="px330",
        cell_strain="Any mammalian",
        selection="Amp",
        notes=(
            "All-in-one SpCas9 + sgRNA vector.  BbsI digest leaves CACC/AAAC "
            "overhangs.  5′G is prepended when absent because U6 transcription "
            "efficiency is higher when the +1 position is a guanosine. Addgene "
            "#42230 lists pX330 as a mammalian CRISPR plasmid with pUC ori "
            "backbone, Ampicillin resistance, Stbl3 growth at 37°C, and high "
            "copy number. The local backbone sequence is the public NovoPro "
            "GenBank map for pX330-U6-Chimeric_BB-CBh-hSpCas9."
        ),
        citations=(
            Citation("Cong et al. Science 2013",
                     "https://doi.org/10.1126/science.1231143",
                     "Introduced pX330 for mammalian CRISPR-Cas9 editing"),
            Citation("Ran et al. Nat Protoc 2013",
                     "https://doi.org/10.1038/nprot.2013.143",
                     "BbsI CACC/AAAC annealed-oligo cloning protocol"),
            Citation("Addgene #42230",
                     "https://www.addgene.org/42230/",
                     "pX330 plasmid repository record: pUC ori backbone; Ampicillin 100 μg/mL; Stbl3; 37°C; high copy"),
            Citation("NovoPro V005940 GenBank map",
                     "https://www.novoprolabs.com/vector/Vgy3dima",
                     "Public 8484 bp pX330-U6-Chimeric_BB-CBh-hSpCas9 GenBank sequence used as local backbone resource"),
        ),
    ),

    "lenticrispr_v2": VectorSpec(
        name="lentiCRISPR v2",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BsmBI",
        promoter="U6",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="CGTCTC",
        top_overhang="CACC",
        bottom_overhang="AAAC",
        u6_prefers_5prime_g=True,
        cell_strain="Any mammalian (lentiviral delivery)",
        selection="Puro",
        notes=(
            "Lentiviral all-in-one vector for stable Cas9 + guide integration.  "
            "BsmBI is a Type IIS enzyme; lentiCRISPR v2 uses CACC/AAAC overhangs "
            "that match the pX330-BbsI convention, so the oligo design is the same "
            "for these two specific vectors. Other BsmBI vectors may use different overhangs."
        ),
        citations=(
            Citation("Sanjana et al. Nat Methods 2014",
                     "https://doi.org/10.1038/nmeth.3047",
                     "lentiCRISPR v2 BsmBI CACC/AAAC cloning protocol"),
            Citation("Addgene #52961",
                     "https://www.addgene.org/52961/",
                     "lentiCRISPR v2 plasmid repository record"),
        ),
    ),

    "pdr274": VectorSpec(
        name="pDR274",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BsaI",
        promoter="T7",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="GGTCTC",
        top_overhang="TAGG",
        bottom_overhang="AAAC",
        # T7 requires G at +1, but the TAGG overhang already provides it as its
        # last base — no additional prepend is applied (u6_prefers_5prime_g=False).
        u6_prefers_5prime_g=False,
        cell_strain="Zebrafish embryo (microinjection)",
        selection="Kan",
        notes=(
            "Zebrafish sgRNA-only vector (guide RNA cassette; no Cas9).  T7 RNA "
            "polymerase drives in vitro transcription; SpCas9 mRNA or protein is "
            "co-injected separately.  TAGG overhang encodes the required +1 G; "
            "no extra 5′G is prepended."
        ),
        citations=(
            Citation("Hwang et al. PLoS One 2013",
                     "https://doi.org/10.1371/journal.pone.0068708",
                     "pDR274 BsaI TAGG/AAAC zebrafish sgRNA cloning"),
            Citation("Varshney et al. Nat Protoc 2016",
                     "https://doi.org/10.1038/nprot.2016.099",
                     "Zebrafish guide RNA cloning and microinjection protocol"),
            Citation("Addgene #42250",
                     "https://www.addgene.org/42250/",
                     "pDR274 plasmid repository record"),
        ),
    ),

    "pcrispr": VectorSpec(
        name="pCRISPR",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BsaI",
        promoter="J23119",
        nuclease_system="SpCas9",      # Cas9 lives on companion pCas9 plasmid
        guide_type="spacer_array",
        scaffold_in_vector=True,
        recognition_site="GGTCTC",
        top_overhang="AAAC",
        bottom_overhang="AAAAC",       # 5-nt bottom overhang is specific to this vector
        u6_prefers_5prime_g=False,
        backbone_resource="pCRISPR_rpsL",
        cell_strain="HME63 or MG1655 carrying pCas9",
        selection="Kan",
        notes=(
            "Two-plasmid E. coli system: pCRISPR carries the spacer array under "
            "J23119; pCas9 provides Cas9 and tracrRNA.  BsaI TypeIIS cloning.  "
            "The AAAAC bottom overhang (5 nt) is wider than the mammalian 4-nt AAAC."
        ),
        citations=(
            Citation("Jiang et al. Nat Biotechnol 2013",
                     "https://doi.org/10.1038/nbt.2508",
                     "pCas9/pCRISPR two-plasmid E. coli system; BsaI guide insertion"),
            Citation("Addgene #44505",
                     "https://www.addgene.org/44505/",
                     "pCRISPR plasmid repository record"),
        ),
    ),
"prc11_lbcpf1": VectorSpec(
    name="pRC11-U6-DR-crRNA-BsmBI(x2); EFS-Puro-WPRE",
    cloning_method="TypeIISOligoCloning",
    dna_source="annealed_oligos",
    enzyme="BsmBI",
    promoter="U6",
    nuclease_system="Cas12a",
    guide_type="crRNA",
    scaffold_in_vector=True,
    recognition_site="CGTCTC",
    top_overhang="",
    bottom_overhang="",
    u6_prefers_5prime_g=False,
    cell_strain="Mammalian lentiviral delivery",
    selection="Puro",
    notes=(
        "Lentiviral crRNA-only vector for use with LbCas12a (crRNA cassette; no "
        "Cas12a/Cpf1 nuclease). LbCas12a is provided by a separate construct or "
        "stable cell line. Uses a U6 direct-repeat crRNA cassette and BsmBI cloning. "
        "Overhangs must be verified from the plasmid map/protocol."
    ),
    citations=(
        Citation("Chow et al Nat Methods. 2019 May;16(5):405–408",
                 "https://doi.org/10.1038/s41592-019-0371-5",
                 "U6-driven crRNA cassette; LbCpf1 direct repeat; BsmBI cloning cassette; lentiviral delivery; puromycin selection"),
        Citation("Addgene #123360",
                "https://www.addgene.org/123360/",
                "Plasmid map verification"),
    ),
),

"bpk4446_fncas12a": VectorSpec(
    name="pUC19-U6-FnCas12a crRNA-BsmBI cassette (BPK4446)",
    cloning_method="TypeIISOligoCloning",
    dna_source="annealed_oligos",
    enzyme="BsmBI",
    promoter="U6",
    nuclease_system="Cas12a",
    guide_type="crRNA",
    scaffold_in_vector=True,
    recognition_site="CGTCTC",
    top_overhang="",
    bottom_overhang="",
    u6_prefers_5prime_g=False,
    cell_strain="Mammalian / cloning entry vector",
    selection="Amp",
    notes=(
        "crRNA-only entry vector for use with FnCas12a (crRNA cassette; no Cas12a "
        "nuclease). FnCas12a is provided by a separate construct. U6 promoter; "
        "spacer oligos are cloned into a BsmBI cassette. Overhangs must be verified "
        "from the plasmid map/protocol."
    ),
    citations=(
        Citation(
            "Kleinstiver et al. Nat Biotechnol 2019",
            "https://doi.org/10.1038/s41587-018-0011-0",
            "FnCas12a crRNA cloning vector using BsmBI spacer insertion"),
        Citation(
            "Addgene #114087",
            "https://www.addgene.org/114087/",
            "FnCas12a crRNA entry vector; clone spacer oligos into BsmBI cassette"),
    ),
),

    # ── Standard Restriction Ligation ────────────────────────────────────────
    # The insert is PCR-amplified with primers whose 5′ tails carry the RE site.
    # Digest of both PCR product and vector yields compatible sticky ends for
    # T4 ligation.  The enzyme cuts the same site in insert and vector — the
    # recognition site is disrupted at the ligation junction.
    # This is biologically distinct from Type IIS cloning.

    "ptargetf": VectorSpec(
        name="pTargetF",
        cloning_method="RestrictionLigation",
        dna_source="pcr_product",
        enzyme="SpeI",
        promoter="J23119",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=False,      # full cassette (promoter+guide+scaffold) is the insert
        recognition_site="ACTAGT",     # SpeI: A↓CTAGT → 5′-CTAG-3′ overhang
        top_overhang="CTAGT",
        bottom_overhang="A",
        u6_prefers_5prime_g=False,
        restriction_site_sequence="ACTAGT",
        cell_strain="E. coli MG1655 or equivalent",
        selection="Amp",
        notes=(
            "Used with COMPANION pCas9-CR4 plasmid for E. coli editing.  The "
            "sgRNA cassette (J23119 + guide + scaffold) is PCR-amplified with "
            "SpeI-site-tailed primers, then restriction-ligated into pTargetF.  "
            "SpeI and XbaI leave compatible CTAG overhangs — the ligation product "
            "is NOT re-cuttable by either enzyme alone."
        ),
        citations=(
            Citation("Jiang et al. Appl Environ Microbiol. 2015 Apr;81(7):2506-14. doi: 10.1128/AEM.04023-14. Epub 2015 Jan 30.",
                     "https://journals.asm.org/doi/10.1128/aem.04023-14",
                     "pTargetF/pCas9-CR4 two-plasmid system; SpeI restriction ligation"),
            Citation("Addgene #62226",
                     "https://www.addgene.org/62226/",
                     "pTargetF plasmid repository record"),
            Citation("NEB SpeI product page",
                     "https://www.neb.com/en-us/products/r0133-spei-hf",
                     "SpeI recognition ACTAGT; leaves 5′-CTAG-3′ overhang"),
        ),
    ),

    # ── Gibson Assembly ──────────────────────────────────────────────────────
    # Isothermal single-tube reaction (T5 exonuclease + Phusion polymerase +
    # Taq ligase).  Adjacent fragments share 15–30 bp of homologous sequence.
    # No RE sites or sticky ends in the final product.

    "px458_gibson": VectorSpec(
        name="pX458 (Gibson)",
        cloning_method="GibsonAssembly",
        dna_source="overlap_fragment",
        enzyme="Gibson",
        promoter="U6",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        u6_prefers_5prime_g=True,
        recommended_overlap_bp=20,
        # U6 promoter flanking sequences at the BbsI guide insertion site
        # (positions 251-267 in pX330 GenBank; shared with pX458).
        gibson_left_context="ATATATCTTGTGGAAAGGACGAAACACCGG",
        gibson_right_context="GTTTTAGAGCTAGAAATAGCAAGTTAAAAT",
        cell_strain="Any mammalian",
        selection="Amp",
        notes=(
            "pX458 can accept a guide insert by Gibson assembly if the insert "
            "carries overlaps matching the sequences flanking the linearisation "
            "site. Vector linearisation primers are NOT designed by this tool."
        ),
        citations=(
            Citation("Ran et al. Nat Protoc 2013",
                     "https://doi.org/10.1038/nprot.2013.143",
                     "pX458 plasmid for Cas9 + GFP reporter"),
            Citation("Addgene #48138",
                     "https://www.addgene.org/48138/",
                     "pX458 plasmid repository record"),
            Citation("Gibson et al. Nat Methods 2009",
                     "https://doi.org/10.1038/nmeth.1318",
                     "Gibson assembly isothermal multi-fragment cloning"),
        ),
    ),

    # ── Plant binary vectors (Agrobacterium T-DNA delivery) ─────────────────
    # AtU6-26 driven sgRNA cloned via BsaI TypeIIS sites into a binary vector
    # backbone.  BsaI leaves ACCG (top) and AAAC (bottom) 4-nt overhangs.
    # Verify overhangs against the plasmid map before ordering oligos.

    "pkse401": VectorSpec(
        name="pKSE401",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BsaI",
        promoter="AtU6-26",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="GGTCTC",
        top_overhang="ACCG",
        bottom_overhang="AAAC",
        u6_prefers_5prime_g=True,
        cell_strain="Arabidopsis thaliana (Agrobacterium floral dip)",
        selection="Spec",
        notes=(
            "Binary vector for Arabidopsis CRISPR. SpCas9 under 35S; sgRNA under "
            "AtU6-26. BsaI digest yields ACCG/AAAC overhangs. Verified against "
            "Addgene #62203 map before ordering oligos."
        ),
        citations=(
            Citation("Xing et al. BMC Plant Biol 2014",
                     "https://doi.org/10.1186/s12870-014-0327-y",
                     "pKSE401 BsaI ACCG/AAAC Arabidopsis sgRNA cloning protocol"),
            Citation("Addgene #62203",
                     "https://www.addgene.org/62203/",
                     "pKSE401 plasmid repository record"),
        ),
    ),

    "phee401e": VectorSpec(
        name="pHEE401E",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BsaI",
        promoter="AtU6-26",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="GGTCTC",
        top_overhang="ACCG",
        bottom_overhang="AAAC",
        u6_prefers_5prime_g=True,
        cell_strain="Arabidopsis thaliana (Agrobacterium floral dip)",
        selection="Hyg",
        notes=(
            "High-efficiency Arabidopsis CRISPR via egg cell-specific EC1.2en-EC1.1p "
            "promoter driving Cas9; improves heritable T1 editing rates over 35S-Cas9. "
            "sgRNA cloning (BsaI ACCG/AAAC) is identical to pKSE401."
        ),
        citations=(
            Citation("Wang et al. Plant Cell 2015",
                     "https://doi.org/10.1105/tpc.15.00454",
                     "pHEE401E egg cell-specific Cas9 for heritable Arabidopsis editing"),
            Citation("Addgene #71288",
                     "https://www.addgene.org/71288/",
                     "pHEE401E plasmid repository record"),
        ),
    ),

    # ── Yeast (Saccharomyces cerevisiae) ─────────────────────────────────────
    # pML104/pML107: all-in-one S. cerevisiae vectors. SpCas9 under pTDH3;
    # sgRNA under SNR52 RNA Pol III promoter. The guide is cloned into a
    # BclI-SwaI cassette: BclI leaves a 5' GATC overhang and SwaI is blunt.

    "pml104": VectorSpec(
        name="pML104",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BclI-SwaI",
        promoter="SNR52",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="TGATCA",   # BclI recognition site (Dam-methylation-sensitive)
        top_overhang="GATC",         # BclI 5′ overhang after digest
        bottom_overhang="",          # SwaI produces a blunt end
        u6_prefers_5prime_g=False,   # SNR52 does not require a 5′G
        cell_strain="Saccharomyces cerevisiae (lithium acetate transformation)",
        selection="URA3 (yeast) / Ampicillin (bacteria)",
        notes=(
            "All-in-one S. cerevisiae CRISPR vector (pRSII426 backbone, 2μ, high copy). "
            "Cas9 under pTDH3 (GAP) promoter; sgRNA under pSNR52 RNA Pol III promoter. "
            "Guide spacer cloned by annealing two oligos into BclI/SwaI-digested vector; "
            "the top oligo has a 5′ GATC BclI-compatible overhang and the opposite end is blunt. "
            "IMPORTANT: BclI is Dam-methylation-sensitive — purify plasmid from a "
            "Dam- E. coli strain (e.g. SCS110) before digestion or BclI will not cut. "
            "Include GTTTTAGAGCTAG after the 20 nt guide to rebuild the 5′ sgRNA segment."
        ),
        citations=(
            Citation("Laughery et al. Yeast 2015",
                     "https://doi.org/10.1002/yea.3098",
                     "pML104/pML107 all-in-one S. cerevisiae CRISPR vectors; BclI-SwaI guide spacer cloning"),
            Citation("Addgene #67638",
                     "https://www.addgene.org/67638/",
                     "pML104 plasmid repository record; pRSII426 backbone; URA3 selection; BclI-SwaI guide cloning"),
        ),
    ),

    "pml107": VectorSpec(
        name="pML107",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BclI-SwaI",
        promoter="SNR52",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="TGATCA",
        top_overhang="GATC",
        bottom_overhang="",
        u6_prefers_5prime_g=False,
        cell_strain="Saccharomyces cerevisiae (lithium acetate transformation)",
        selection="LEU2 (yeast) / Ampicillin (bacteria)",
        notes=(
            "All-in-one S. cerevisiae CRISPR vector (pRSII425 backbone, 2μ, high copy). "
            "Cas9 under pTDH3 (GAP) promoter; sgRNA under pSNR52 RNA Pol III promoter. "
            "Guide spacer cloned by annealing two oligos into BclI/SwaI-digested vector; "
            "the top oligo has a 5′ GATC BclI-compatible overhang and the opposite end is blunt. "
            "IMPORTANT: BclI is Dam-methylation-sensitive — purify plasmid from a "
            "Dam- E. coli strain (e.g. SCS110) before digestion or BclI will not cut. "
            "Include GTTTTAGAGCTAG after the 20 nt guide to rebuild the 5′ sgRNA segment."
        ),
        citations=(
            Citation("Laughery et al. Yeast 2015",
                     "https://doi.org/10.1002/yea.3098",
                     "pML104/pML107 all-in-one S. cerevisiae CRISPR vectors; BclI-SwaI guide spacer cloning"),
            Citation("Addgene #67639",
                     "https://www.addgene.org/67639/",
                     "pML107 plasmid repository record; pRSII425 backbone; LEU2 selection; BclI-SwaI guide cloning"),
        ),
    ),

    # ── Caenorhabditis elegans ────────────────────────────────────────────────
    # pDD162: Peft-3::Cas9 drives Cas9 in soma and germline; U6 promoter drives
    # sgRNA.  Guide inserted by Gibson assembly (Ligation Independent Cloning).
    # Vector linearisation primers must be supplied separately — see Addgene
    # #47549 protocol and the Goldstein Lab CRISPR collection notes.

    "pdd162": VectorSpec(
        name="pDD162",
        cloning_method="GibsonAssembly",
        dna_source="overlap_fragment",
        enzyme="Gibson",
        promoter="U6",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recommended_overlap_bp=20,
        u6_prefers_5prime_g=False,
        cell_strain="Caenorhabditis elegans (germline microinjection)",
        selection="Amp",
        notes=(
            "C. elegans all-in-one vector.  Peft-3 drives Cas9 in soma and "
            "germline; U6 promoter drives the sgRNA.  Guide is inserted by Gibson "
            "assembly (Ligation Independent Cloning) into a pre-built sgRNA "
            "scaffold that is missing its 19-nt targeting sequence.  Provide "
            "left_overlap_context and right_overlap_context from the plasmid map "
            "flanking the empty targeting slot.  Vector linearisation primers are "
            "NOT designed by this tool — add them manually before calling "
            "create_construction_file with assembly_strategy='Gibson'."
        ),
        citations=(
            Citation("Dickinson et al. Nat Methods 2013",
                     "https://doi.org/10.1038/nmeth.2641",
                     "pDD162 C. elegans Cas9+sgRNA vector; Gibson LIC guide insertion"),
            Citation("Addgene #47549",
                     "https://www.addgene.org/47549/",
                     "pDD162 plasmid repository record; cloning method: Gibson LIC"),
        ),
    ),

    # ── Drosophila melanogaster ───────────────────────────────────────────────
    # pCFD3: dU6:3-driven single-guide vector.  SpCas9 is supplied in trans
    # from a stable nos-cas9 or act5c-cas9 transgenic line.  BbsI TypeIIS
    # cloning.  Overhangs must be verified from the plasmid map.

    "pcfd3": VectorSpec(
        name="pCFD3",
        cloning_method="TypeIISOligoCloning",
        dna_source="annealed_oligos",
        enzyme="BbsI",
        promoter="U6:3",
        nuclease_system="SpCas9",
        guide_type="sgRNA",
        scaffold_in_vector=True,
        recognition_site="GAAGAC",
        top_overhang="",    # verify from Addgene #49410 plasmid map / protocol
        bottom_overhang="", # verify from Addgene #49410 plasmid map / protocol
        u6_prefers_5prime_g=False,  # dU6:3 does not enforce a 5′G requirement
        cell_strain="Drosophila melanogaster (embryo injection)",
        selection="Amp",
        notes=(
            "Single-guide Drosophila CRISPR vector.  dU6:3 RNA Pol III promoter "
            "drives the sgRNA; SpCas9 is provided in trans by a nos-cas9 or "
            "act5c-cas9 transgenic line.  BbsI TypeIIS cloning — verify exact "
            "4-nt overhangs from the plasmid map/protocol before ordering oligos."
        ),
        citations=(
            Citation("Port et al. PNAS 2014",
                     "https://doi.org/10.1073/pnas.1405500111",
                     "pCFD3 single-guide Drosophila vector; BbsI cloning; Cas9 in trans"),
            Citation("Addgene #49410",
                     "https://www.addgene.org/49410/",
                     "pCFD3 plasmid repository record"),
        ),
    ),

}

# Workflow-level citations (method references independent of any vector)
WORKFLOW_CITATIONS: dict[str, tuple[Citation, ...]] = {
    "TypeIISOligoCloning": (
        Citation("Ran et al. Nat Protoc 2013",
                 "https://doi.org/10.1038/nprot.2013.143",
                 "General annealed-oligo Type IIS guide cloning protocol"),
        Citation("Addgene, Cloning of Oligos for sgRNA/shRNA (Zhang Lab, Sept 2015)",
                 "https://www.addgene.org/protocols/cloning-of-oligos/",
                 "Standard protocol for BbsI/BsmBI annealed-oligo ligation into Cas9 vectors"),
        Citation("NEB FAQ: Which restriction enzymes are used in Golden Gate assembly?",
                 "https://www.neb.com/en-us/faqs/which-restriction-enzymes-are-used-in-golden-gate-assembly",
                 "BbsI, BsmBI, and BsaI cut outside their recognition site leaving 4-nt overhangs used in annealed-oligo guide cloning"),
        Citation("Hu et al. 2020, The Crop Journal 8(3):403-407",
                 "https://doi.org/10.1016/j.cj.2019.06.007",
                 "Cas12a crRNA oligo cloning via BbsI/BsmBI into array vectors"),
    ),
    "RestrictionLigation": (
        Citation("Sambrook & Russell, Molecular Cloning 4th ed.",
                 "ISBN 978-1-936113-42-2",
                 "Standard restriction enzyme ligation cloning methodology"),
    ),
    "GibsonAssembly": (
        Citation("Gibson et al. Nat Methods 2009",
                 "https://doi.org/10.1038/nmeth.1318",
                 "Isothermal Gibson assembly of overlapping DNA fragments"),
        Citation("NEB HiFi Assembly usage guidelines",
                 "https://www.neb.com/en-us/tools-and-resources/usage-guidelines/nebuilder-hifi-dna-assembly",
                 "Practical overlap design guidelines: 15–30 bp, Tm ≥ 55 °C"),
    ),
    "GoldenGateAssembly": (
        Citation("Engler et al. PLoS One 2008",
                 "https://doi.org/10.1371/journal.pone.0003647",
                 "Original Golden Gate assembly using Type IIS enzymes with user-designed overhangs"),
        Citation("NEB Golden Gate Assembly",
                 "https://www.neb.com/applications/cloning/golden-gate-assembly",
                 "NEB Golden Gate reagents and overhang uniqueness design guidelines"),
    ),
}

# ---------------------------------------------------------------------------
# Golden Gate enzyme specifications
# recognition: Type IIS recognition sequence (written 5'→3')
# spacer_len:  number of bases between recognition site and the cut point
#              (= number of N's in the primer tail between site and overhang)
# ---------------------------------------------------------------------------

GOLDEN_GATE_ENZYME_SPECS: dict[str, dict] = {
    "BsaI":  {"recognition": "GGTCTC", "spacer_len": 1},
    "BsmBI": {"recognition": "CGTCTC", "spacer_len": 1},
    "Esp3I": {"recognition": "CGTCTC", "spacer_len": 1},  # BsmBI isoschizomer
    "BbsI":  {"recognition": "GAAGAC", "spacer_len": 2},
}

# ---------------------------------------------------------------------------
# Addgene API fallback: common TypeIIS overhangs for well-known pX330/lentiCRISPR-
# family vectors. CACC/AAAC is the predominant convention for BbsI and BsmBI in
# those lineages, but is NOT universal — other vectors using the same enzyme may
# differ. Always verify against the actual plasmid protocol before ordering oligos.
# BsaI overhangs vary by vector; always ask the user.
# ---------------------------------------------------------------------------

_TYPEII_ENZYME_OVERHANGS: dict[str, dict] = {
    "BbsI":  {"top": "CACC", "bottom": "AAAC", "recognition": "GAAGAC"},
    "BsmBI": {"top": "CACC", "bottom": "AAAC", "recognition": "CGTCTC"},
}


def _parse_addgene_id(vector: str) -> Optional[int]:
    """Return numeric Addgene ID if vector looks like one, else None."""
    v = vector.strip()
    if v.lower().startswith("addgene:"):
        v = v[8:].strip()
    return int(v) if v.isdigit() else None


def _addgene_data_to_vector_spec(
    data: dict,
) -> "tuple[Optional[VectorSpec], Optional[dict]]":
    """
    Try to build a VectorSpec from a FetchAddgeneVector result dict.

    Returns (VectorSpec, None) on success, or (None, needs_user_input_dict)
    when critical fields can't be inferred (e.g. BsaI overhangs are vector-
    specific and are not in the standard lookup table).
    """
    addgene_id = data.get("addgene_id", "")
    name = data.get("name") or f"Addgene #{addgene_id}"
    url = data.get("url") or f"https://www.addgene.org/{addgene_id}/"
    description = (data.get("description") or "").lower()

    # Map raw clone_method string to our cloning method keys
    raw = (data.get("clone_method_raw") or "").lower()
    if "gibson" in raw:
        cloning_method = "GibsonAssembly"
        dna_source = "overlap_fragment"
    elif "golden gate" in raw:
        cloning_method = "GoldenGateAssembly"
        dna_source = "overlap_fragment"
    elif "restriction" in raw and "type ii" not in raw:
        cloning_method = "RestrictionLigation"
        dna_source = "pcr_product"
    else:
        # Covers "Type IIS Restriction Enzyme Cloning" and unknowns
        cloning_method = "TypeIISOligoCloning"
        dna_source = "annealed_oligos"

    enzyme = data.get("enzyme") or ""
    promoter = data.get("promoter") or ""

    # Infer nuclease system from description / name
    if any(kw in description for kw in ("cas12a", "cpf1", "cas12a")):
        nuclease_system, guide_type = "Cas12a", "crRNA"
    elif "sacas9" in description:
        nuclease_system, guide_type = "SaCas9", "sgRNA"
    else:
        nuclease_system, guide_type = "SpCas9", "sgRNA"

    # Selection markers
    all_markers = [
        m for m in (data.get("resistance_markers") or []) + [data.get("bacterial_resistance") or ""]
        if m
    ]
    selection = ", ".join(sorted(set(all_markers)))
    cell_strain = data.get("growth_strain") or ""

    doi = data.get("article_doi") or ""
    citations: tuple[Citation, ...] = tuple(filter(None, [
        Citation(f"Article (DOI: {doi})", f"https://doi.org/{doi}", f"Original publication for {name}") if doi else None,
        Citation(f"Addgene #{addgene_id}", url, f"{name} plasmid repository record"),
    ]))
    notes = (
        f"Dynamically fetched from Addgene #{addgene_id}. "
        "Verify overhangs and cloning parameters against the plasmid protocol before ordering oligos."
    )

    if cloning_method == "TypeIISOligoCloning":
        enz_info = _TYPEII_ENZYME_OVERHANGS.get(enzyme, {})
        top_overhang = enz_info.get("top", "")
        bottom_overhang = enz_info.get("bottom", "")
        recognition_site = enz_info.get("recognition", "")

        if not top_overhang or not bottom_overhang:
            return None, {
                "status": "needs_user_input",
                "addgene_id": addgene_id,
                "plasmid_name": name,
                "cloning_method": cloning_method,
                "enzyme": enzyme,
                "promoter": promoter,
                "missing_fields": ["top_overhang", "bottom_overhang"],
                "questions": [
                    f"Vector '{name}' (Addgene #{addgene_id}) uses {enzyme} — "
                    f"standard overhangs for {enzyme} are not in the lookup table "
                    f"(they vary by vector). Please provide top and bottom 5′ "
                    f"overhangs from the plasmid protocol at {url}."
                ],
                "partial_spec": {
                    "name": name,
                    "cloning_method": cloning_method,
                    "enzyme": enzyme,
                    "promoter": promoter,
                    "nuclease_system": nuclease_system,
                },
            }

        spec = VectorSpec(
            name=name,
            cloning_method=cloning_method,
            dna_source=dna_source,
            enzyme=enzyme,
            promoter=promoter,
            nuclease_system=nuclease_system,
            guide_type=guide_type,
            scaffold_in_vector=True,
            recognition_site=recognition_site,
            top_overhang=top_overhang,
            bottom_overhang=bottom_overhang,
            u6_prefers_5prime_g="U6" in promoter.upper(),
            cell_strain=cell_strain,
            selection=selection,
            notes=notes,
            citations=citations,
        )
        return spec, None

    if cloning_method == "GibsonAssembly":
        spec = VectorSpec(
            name=name,
            cloning_method=cloning_method,
            dna_source=dna_source,
            enzyme="Gibson",
            promoter=promoter,
            nuclease_system=nuclease_system,
            guide_type=guide_type,
            scaffold_in_vector=True,
            recommended_overlap_bp=20,
            cell_strain=cell_strain,
            selection=selection,
            notes=notes,
            citations=citations,
        )
        return spec, None

    # RestrictionLigation / GoldenGate / unknown — partial info only
    return None, {
        "status": "needs_user_input",
        "addgene_id": addgene_id,
        "plasmid_name": name,
        "cloning_method": cloning_method,
        "enzyme": enzyme,
        "missing_fields": ["cloning_parameters"],
        "questions": [
            f"Vector '{name}' (Addgene #{addgene_id}) uses {cloning_method}. "
            f"Please supply the full cloning parameters (overhangs, restriction sites) "
            f"from the protocol at {url}."
        ],
    }


# ---------------------------------------------------------------------------
# Known promoter and scaffold sequences for auto-cassette assembly
# Used by RestrictionLigation when scaffold_in_vector=False and the promoter
# and nuclease are both known — avoids asking the user for the full cassette.
# Sequences sourced from published protocols; verify before ordering.
# ---------------------------------------------------------------------------

PROMOTER_SEQUENCES: dict[str, str] = {
    # Anderson J23119 constitutive promoter (iGEM standard part BBa_J23119)
    "J23119": "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC",
}

SCAFFOLD_SEQUENCES: dict[str, str] = {
    # Standard SpCas9 sgRNA scaffold (Jinek et al. 2012 / Cong et al. 2013)
    "SpCas9": "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGT"
              "GGCACCGAGTCGGTGC",
}


# ---------------------------------------------------------------------------
# Organism compatibility map
# Maps lowercase target-organism keywords to substrings expected in
# VectorSpec.cell_strain when the vector is compatible with that organism.
# ---------------------------------------------------------------------------

_ORGANISM_COMPAT: dict[str, list[str]] = {
    "escherichia coli":          ["e. coli", "mg1655", "hme63", "mach1", "dh5"],
    "e. coli":                   ["e. coli", "mg1655", "hme63", "mach1", "dh5"],
    "bacteria":                  ["e. coli", "mg1655", "hme63", "mach1", "dh5"],
    "mammalian":                 ["mammalian"],
    "human":                     ["mammalian"],
    "homo sapiens":              ["mammalian"],
    "mouse":                     ["mammalian"],
    "mus musculus":              ["mammalian"],
    "zebrafish":                 ["zebrafish"],
    "danio rerio":               ["zebrafish"],
    "arabidopsis":               ["arabidopsis"],
    "arabidopsis thaliana":      ["arabidopsis"],
    "plant":                     ["arabidopsis", "agrobacterium"],
    "saccharomyces cerevisiae":  ["saccharomyces"],
    "s. cerevisiae":             ["saccharomyces"],
    "yeast":                     ["saccharomyces"],
    "caenorhabditis elegans":    ["caenorhabditis"],
    "c. elegans":                ["caenorhabditis"],
    "drosophila":                ["drosophila"],
    "drosophila melanogaster":   ["drosophila"],
}


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

_YEAST_BCLI_SWAI_SGRNA_PREFIX = "GTTTTAGAGCTAG"
_IUPAC_AMBIGUITY = set("RYSWKMBDHVN")


def _validate_dna(seq: str, label: str) -> str:
    seq = seq.upper().strip()
    if not seq:
        raise ValueError(f"{label} must not be empty.")
    invalid = sorted(set(seq) - set("ATGC"))
    if invalid:
        raise ValueError(
            f"Invalid base(s) in {label}: {invalid}. Only A, T, G, C are accepted."
        )
    return seq


def _normalize_resource_dna(seq: str, label: str) -> str:
    """Normalize backbone resource sequences; ambiguous IUPAC bases become N."""
    seq = seq.upper().strip()
    if not seq:
        raise ValueError(f"{label} must not be empty.")
    normalized = "".join("N" if base in _IUPAC_AMBIGUITY else base for base in seq)
    invalid = sorted(set(normalized) - set("ATGCN"))
    if invalid:
        raise ValueError(
            f"Invalid base(s) in {label}: {invalid}. Only A, T, G, C, N and IUPAC ambiguity codes are accepted in backbone resources."
        )
    return normalized



def _uses_yeast_bcli_swai_cloning(spec: Optional[VectorSpec]) -> bool:
    return (
        spec is not None
        and spec.promoter.upper() == "SNR52"
        and spec.enzyme == "BclI-SwaI"
        and spec.top_overhang == "GATC"
    )


def _read_sequence_file(path: Path) -> str:
    """Return uppercase DNA string from a FASTA or GenBank file."""
    if not path.exists():
        raise ValueError(f"Required local sequence file is missing: {path}")
    text = path.read_text()

    if path.suffix.lower() in {".fa", ".fasta", ".fna"}:
        seq = "".join(
            line.strip()
            for line in text.splitlines()
            if line.strip() and not line.startswith(">")
        )
        return _normalize_resource_dna(seq, path.name)

    if path.suffix.lower() in {".gb", ".gbk", ".genbank"}:
        in_origin, parts = False, []
        for line in text.splitlines():
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin and line.startswith("//"):
                break
            if in_origin:
                parts.extend(c for c in line if c.isalpha())
        return _normalize_resource_dna("".join(parts), path.name)

    raise ValueError(f"Unsupported sequence file type: {path.suffix}")


def _format_citations(citations: tuple[Citation, ...]) -> list[dict]:
    return [{"label": c.label, "reference": c.url_or_reference, "claim": c.claim}
            for c in citations]


def _needs_user_input(
    cloning_method: str,
    vector: str,
    missing_fields: list[str],
    questions: list[str],
    note: str = "",
) -> dict:
    """Construct a standard 'needs_user_input' response."""
    return {
        "status": "needs_user_input",
        "cloning_method": cloning_method,
        "vector": vector,
        "missing_fields": missing_fields,
        "questions": questions,
        "notes": note,
    }


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class CRISPRCloningDesigner:
    """
    Description:
        Designs oligos or primers for CRISPR guide insertion via TypeIIS ligation,
        restriction ligation, Gibson, or Golden Gate assembly.

    Input:
        vector (str): Known vector key (e.g. "px330") or "custom".
        protospacer (str): 20-nt guide for TypeIIS cloning.
        insert_sequence (str): Insert DNA for Gibson/Golden Gate.
        left_overlap_context (str): Vector sequence 5′ of insertion site.
        right_overlap_context (str): Vector sequence 3′ of insertion site.
        cloning_method (str): Required when vector="custom".

    Output:
        dict: status="ready" with oligo/primer sequences, or
              status="needs_user_input" with missing fields and questions.

    Tests:
        - Case:
            Input: vector="px330", protospacer="ACAGAAACCTGCCAGTTTGC"
            Expected Output: status == "ready", top_oligo starts with "CACC"
            Description: TypeIIS cloning with a known vector.
        - Case:
            Input: vector="px330" (no protospacer)
            Expected Output: status == "needs_user_input"
            Description: Missing protospacer returns structured prompt.
        - Case:
            Input: vector="px330", protospacer="INVALID!!"
            Expected Exception: ValueError
            Description: Non-DNA protospacer raises ValueError.
    """

    def initiate(self) -> None:
        """No mutable state; exists for JSON tool harness compatibility."""
        pass

    # ── Vector resolution ────────────────────────────────────────────────────
    cas12a_aliases = {"cpf1", "cas12a", "cas12a_custom", "cpf1_custom"}

    def resolve_vector(self, vector: Optional[str]) -> Optional[VectorSpec]:
        """
        Return the VectorSpec for a known vector key, or None for "custom".
        Raises ValueError only for a non-empty, non-custom, unrecognised name.
        """
        if not vector or vector.lower().strip() == "custom":
            return None
        key = vector.lower().strip()
        if key in self.cas12a_aliases:
            return None 
        if key not in VECTOR_SPECS:
            raise ValueError(
                f"Unknown vector '{vector}'. "
                f"Known vectors: {', '.join(sorted(VECTOR_SPECS))}. "
                "Pass vector='custom' to design with a vector not in this list, "
                "or pass an Addgene plasmid ID (e.g. '42230') to fetch from the Addgene API."
            )
        return VECTOR_SPECS[key]

    # ── Organism compatibility check ─────────────────────────────────────────

    def _check_organism_compatibility(
        self, spec: VectorSpec, target_organism: str
    ) -> Optional[dict]:
        """
        Return a needs_user_input dict if the vector's cell_strain is incompatible
        with target_organism, or None if compatible or organism is unrecognised.
        """
        org_lower = normalize_organism(target_organism).lower()
        strain_lower = spec.cell_strain.lower()

        compatible_substrings: Optional[list[str]] = None
        for key, substrings in _ORGANISM_COMPAT.items():
            if key in org_lower or org_lower in key:
                compatible_substrings = substrings
                break

        if compatible_substrings is None:
            return None  # Unrecognised organism — don't block

        if any(s in strain_lower for s in compatible_substrings):
            return None  # Compatible

        suggestions = [
            k for k, v in VECTOR_SPECS.items()
            if v.nuclease_system not in ("N/A",)
            and any(s in v.cell_strain.lower() for s in compatible_substrings)
        ]

        suggestion_str = (
            f"Compatible preset vectors for {target_organism}: "
            + ", ".join(suggestions) + "."
            if suggestions else
            "No preset vectors are configured for this organism — "
            "pass vector='custom' and supply cloning parameters manually."
        )

        return {
            "status": "needs_user_input",
            "cloning_method": spec.cloning_method,
            "vector": spec.name,
            "missing_fields": ["vector"],
            "questions": [
                f"Vector '{spec.name}' is designed for '{spec.cell_strain}', "
                f"but the target organism is '{target_organism}'. "
                + suggestion_str
            ],
            "notes": (
                f"Organism mismatch: vector cell_strain='{spec.cell_strain}' "
                f"vs target_organism='{target_organism}'. "
                "To skip this check, pass vector='custom'."
            ),
        }

    # ── Requirements assessment ──────────────────────────────────────────────

    def assess_requirements(
        self,
        cloning_method: Optional[str],
        spec: Optional[VectorSpec],
        # TypeIIS inputs
        protospacer: Optional[str],
        top_overhang: Optional[str],
        bottom_overhang: Optional[str],
        enzyme: Optional[str],
        u6_prefers_5prime_g: Optional[bool],
        scaffold_in_vector: Optional[bool],
        promoter: Optional[str],
        # RE ligation inputs
        guide_cassette_sequence: Optional[str],
        restriction_site_sequence: Optional[str],
        # Gibson / Golden Gate inputs
        insert_sequence: Optional[str],
        left_overlap_context: Optional[str],
        right_overlap_context: Optional[str],
        # Golden Gate-specific inputs
        left_overhang: Optional[str],
        right_overhang: Optional[str],
    ) -> dict:
        """
        Inspect available inputs and decide whether we can proceed.

        Returns {"status": "ready"} if all required fields are present.
        Returns {"status": "needs_user_input", ...} with missing fields and
        human-readable questions if any required inputs are absent.

        This method never raises — it is the place where missing information
        is surfaced as structured guidance rather than hard errors.
        """
        # Infer cloning_method from spec if not supplied directly
        method = cloning_method or (spec.cloning_method if spec else None)

        if not method:
            return _needs_user_input(
                cloning_method="unknown",
                vector="unknown",
                missing_fields=["cloning_method"],
                questions=[
                    "Which cloning method do you want to use? "
                    "Options: TypeIISOligoCloning, RestrictionLigation, GibsonAssembly, GoldenGateAssembly."
                ],
                note=(
                    "No vector was recognised and no cloning_method was specified. "
                    "Provide a known vector name or set vector='custom' and specify "
                    "the cloning_method."
                ),
            )

        vector_name = spec.name if spec else "custom"

        if method == "TypeIISOligoCloning":
            return self._assess_typeIIS(
                spec, vector_name, protospacer, top_overhang, bottom_overhang,
                enzyme, u6_prefers_5prime_g, scaffold_in_vector, promoter,
            )

        if method == "RestrictionLigation":
            return self._assess_restriction(
                spec, vector_name, guide_cassette_sequence,
                restriction_site_sequence, enzyme,
            )

        if method == "GibsonAssembly":
            return self._assess_gibson(
                vector_name, insert_sequence,
                left_overlap_context, right_overlap_context,
            )

        if method == "GoldenGateAssembly":
            return self._assess_golden_gate(
                vector_name, insert_sequence,
                left_overlap_context, right_overlap_context,
                left_overhang, right_overhang, enzyme,
            )

        return _needs_user_input(
            cloning_method=method,
            vector=vector_name,
            missing_fields=["cloning_method"],
            questions=[
                f"'{method}' is not a recognised cloning method. "
                "Use: TypeIISOligoCloning, RestrictionLigation, GibsonAssembly, or GoldenGateAssembly."
            ],
        )

    def _assess_typeIIS(
        self,
        spec: Optional[VectorSpec],
        vector_name: str,
        protospacer: Optional[str],
        top_overhang: Optional[str],
        bottom_overhang: Optional[str],
        enzyme: Optional[str],
        u6_prefers_5prime_g: Optional[bool],
        scaffold_in_vector: Optional[bool],
        promoter: Optional[str],
    ) -> dict:
        missing, questions = [], []

        if not protospacer:
            missing.append("protospacer")
            questions.append(
                "What is the protospacer sequence? "
                "(20 nt for SpCas9 / SaCas9, 23 nt for Cas12a; A/T/G/C only)"
            )

        # Spec values are authoritative for known vectors — user-supplied overrides
        # are only used when the spec field is absent or empty (e.g. pcfd3).
        resolved_top_overhang = (spec.top_overhang if spec and spec.top_overhang else top_overhang) or None
        if not resolved_top_overhang:
            missing.append("top_overhang")
            questions.append(
                "Type IIS cloning overhangs are vector-specific and are not determined by the guide sequence. "
                "Please provide the top-strand 5′ overhang from the plasmid map or protocol."
            )

        resolved_bottom_overhang = (spec.bottom_overhang if spec and spec.bottom_overhang else bottom_overhang) or None
        if not resolved_bottom_overhang and not _uses_yeast_bcli_swai_cloning(spec):
            missing.append("bottom_overhang")
            questions.append(
                "Type IIS cloning overhangs are vector-specific and are not determined by the guide sequence. "
                "Please provide the bottom-strand 5′ overhang from the plasmid map or protocol."
                )
        resolved_enzyme = (spec.enzyme if spec and spec.enzyme else enzyme) or None
        if not resolved_enzyme:
            missing.append("enzyme")
            questions.append(
                "Which Type IIS restriction enzyme does your vector use? "
                "(e.g. BbsI, BsmBI, BsaI)"
            )

        resolved_scaffold = scaffold_in_vector if scaffold_in_vector is not None else (
            spec.scaffold_in_vector if spec else None
        )
        if resolved_scaffold is None:
            missing.append("scaffold_in_vector")
            questions.append(
                "Is the guide scaffold or CRISPR RNA direct repeat already encoded in the vector? "
                "(True for most guide-cloning vectors; False if you are cloning "
                "the entire guide RNA/crRNA cassette as the insert)"
            )

        resolved_promoter = promoter or (spec.promoter if spec else None)
        if resolved_promoter is None:
            missing.append("promoter")
            questions.append(
                "What promoter does your vector use to drive the guide RNA or CRISPR RNA? "
                "(e.g. U6, H1, T7 — affects whether a 5′G is needed)"
            )

        resolved_prepend_g = u6_prefers_5prime_g if u6_prefers_5prime_g is not None else (
            spec.u6_prefers_5prime_g if spec else None
        )
        if resolved_prepend_g is None:
            missing.append("u6_prefers_5prime_g")
            questions.append(
                "Should a 5′G be prepended when the protospacer does not already "
                "start with G? (often used for U6/H1-driven guide RNA or CRISPR RNA expression; not needed for T7)"
            )
        if missing: 
            return _needs_user_input(
                cloning_method="TypeIISOligoCloning",
                vector=vector_name,
                missing_fields=missing,
                questions=questions,
            )  
        return {"status": "ready"}
    def _assess_restriction(
        self,
        spec: Optional[VectorSpec],
        vector_name: str,
        guide_cassette_sequence: Optional[str],
        restriction_site_sequence: Optional[str],
        enzyme: Optional[str],
    ) -> dict:
        missing, questions = [], []

        if not guide_cassette_sequence:
            missing.append("guide_cassette_sequence")
            questions.append(
                "What is the full guide cassette sequence to clone? "
                "This should include the promoter, spacer, and guide scaffold or "
                "CRISPR RNA direct repeat if those elements are not already in the vector. "
                "(A/T/G/C only)"
            )

        if spec is None:
            if not enzyme:
                missing.append("enzyme")
                questions.append(
                    "Which restriction enzyme will you use to linearise the vector "
                    "and digest the insert? (e.g. SpeI, EcoRI, HindIII)"
                )
            if not restriction_site_sequence:
                missing.append("restriction_site_sequence")
                questions.append(
                    "What is the restriction enzyme recognition sequence to add to "
                    "your PCR primer 5′ tails? "
                    "(e.g. ACTAGT for SpeI, GAATTC for EcoRI)"
                )

        if missing:
            return _needs_user_input(
                cloning_method="RestrictionLigation",
                vector=vector_name,
                missing_fields=missing,
                questions=questions,
                note=(
                    "For restriction ligation, provide the full guide cassette "
                    "sequence. This tool will design PCR primers with RE-site tails."
                ),
            )
        return {"status": "ready"}

    def _assess_gibson(
        self,
        vector_name: str,
        insert_sequence: Optional[str],
        left_overlap_context: Optional[str],
        right_overlap_context: Optional[str],
    ) -> dict:
        missing, questions = [], []

        if not insert_sequence:
            missing.append("insert_sequence")
            questions.append(
                "What is the insert sequence? "
                "(The sequence you want to introduce into the vector, "
                "without overlap regions — those are computed from the context.)"
            )
        if not left_overlap_context:
            missing.append("left_overlap_context")
            questions.append(
                "What is the vector sequence immediately 5′ of the insertion site? "
                "Provide at least 20–30 bp so an overlap region can be defined. "
                "(Copy from your plasmid map upstream of where the insert goes.)"
            )
        if not right_overlap_context:
            missing.append("right_overlap_context")
            questions.append(
                "What is the vector sequence immediately 3′ of the insertion site? "
                "Provide at least 20–30 bp for the right overlap region."
            )

        if missing:
            return _needs_user_input(
                cloning_method="GibsonAssembly",
                vector=vector_name,
                missing_fields=missing,
                questions=questions,
                note=(
                    "Gibson assembly joins fragments by designed sequence overlaps. "
                    "The left and right overlap contexts come from your plasmid map, "
                    "not from the guide sequence itself."
                ),
            )
        return {"status": "ready"}

    def _assess_golden_gate(
        self,
        vector_name: str,
        insert_sequence: Optional[str],
        left_vector_context: Optional[str],
        right_vector_context: Optional[str],
        left_overhang: Optional[str],
        right_overhang: Optional[str],
        enzyme: Optional[str],
    ) -> dict:
        missing, questions = [], []

        if not insert_sequence:
            missing.append("insert_sequence")
            questions.append(
                "What is the insert sequence to clone in? "
                "(A/T/G/C only; do not include the overhang bases — those are added by the primers.)"
            )
        if not left_vector_context:
            missing.append("left_vector_context")
            questions.append(
                "What is the vector sequence immediately 5′ of the insertion site? "
                "Provide at least 20 bp so a vector reverse primer can be designed."
            )
        if not right_vector_context:
            missing.append("right_vector_context")
            questions.append(
                "What is the vector sequence immediately 3′ of the insertion site? "
                "Provide at least 20 bp so a vector forward primer can be designed."
            )
        if not left_overhang:
            missing.append("left_overhang")
            questions.append(
                "What 4-nt sequence should be used for the left junction overhang? "
                "(e.g. ATCG — this is the sequence in the assembled product reading 5′→3′ "
                "on the top strand at the left boundary of the insert.)"
            )
        if not right_overhang:
            missing.append("right_overhang")
            questions.append(
                "What 4-nt sequence should be used for the right junction overhang? "
                "(e.g. GCTA — must be different from the left overhang and, ideally, "
                "not the reverse complement of any other overhang in the assembly.)"
            )
        if not enzyme:
            missing.append("enzyme")
            questions.append(
                "Which Type IIS enzyme will you use? "
                f"Supported: {', '.join(sorted(GOLDEN_GATE_ENZYME_SPECS))}."
            )

        if missing:
            return _needs_user_input(
                cloning_method="GoldenGateAssembly",
                vector=vector_name,
                missing_fields=missing,
                questions=questions,
                note=(
                    "Golden Gate assembly uses a Type IIS enzyme to expose freely designable "
                    "4-nt overhangs. All four primers (insert fwd/rev + vector fwd/rev) are "
                    "designed so create_construction_file can use assembly_strategy='GoldenGate'."
                ),
            )
        return {"status": "ready"}

    def _resolve_backbone(self, spec: Optional[VectorSpec]) -> tuple[str, str]:
        """Return (backbone_name, backbone_sequence_or_placeholder)."""
        if spec is None or spec.backbone_resource is None:
            name = spec.name if spec else "custom_vector"
            return name, "N"
        key = spec.backbone_resource
        if key not in BACKBONE_RESOURCES:
            raise ValueError(f"No local backbone file configured for '{key}'.")
        return key, _read_sequence_file(BACKBONE_RESOURCES[key])

    # ── Workflow branch A: TypeIIS / Annealed-Oligo Cloning ──────────────────

    def design_typeIIS_oligos(
        self,
        protospacer: str,
        spec: Optional[VectorSpec],
        # Custom-vector overrides (ignored when spec is not None, unless provided)
        top_overhang_override: Optional[str] = None,
        bottom_overhang_override: Optional[str] = None,
        prepend_g_override: Optional[bool] = None,
        enzyme_override: Optional[str] = None,
        construct_name: Optional[str] = None,
    ) -> dict:
        """
        Design top and bottom annealed oligos for Type IIS guide cloning.

        Wet-lab steps simulated
        -----------------------
        1. Digest vector with Type IIS enzyme (BbsI, BsmBI, or BsaI).
           The enzyme recognition site is in the backbone; cutting downstream
           leaves 4-nt 5′ overhangs whose sequence is fixed by the vector.
        2. Design two synthetic oligos:
               top    = 5′-[top_overhang][protospacer]-3′
               bottom = 5′-[bottom_overhang][RC(protospacer)]-3′
        3. Anneal oligos (95 °C → slow-cool) → dsDNA insert with 4-nt overhangs.
        4. Ligate annealed insert into digested vector (T4 ligase, 16 °C o/n).

        The overhangs are properties of the vector, not the guide. Only the
        central portion of each oligo changes when you change the protospacer.

        construction_file_inputs uses assembly_strategy="DirectSynthesis" because
        annealed oligos are synthesized directly — no backbone PCR occurs.
        """
        # Spec overhangs and enzyme are authoritative for known vectors.
        top_oh = (spec.top_overhang if spec and spec.top_overhang else top_overhang_override) or ""
        bot_oh = (spec.bottom_overhang if spec and spec.bottom_overhang else bottom_overhang_override) or ""
        use_g  = prepend_g_override if prepend_g_override is not None else (
            spec.u6_prefers_5prime_g if spec else False
        )
        enz = (spec.enzyme if spec and spec.enzyme else enzyme_override) or "N/A"

        top_oh = _validate_dna(top_oh, "top_overhang")
        bot_oh = _validate_dna(bot_oh, "bottom_overhang") if bot_oh else ""

        g_prepended = False
        original    = protospacer
        if use_g and protospacer[0] != "G":
            protospacer = "G" + protospacer
            g_prepended = True

        if _uses_yeast_bcli_swai_cloning(spec):
            sgRNA_prefix = _YEAST_BCLI_SWAI_SGRNA_PREFIX
            top_oligo = top_oh + protospacer + sgRNA_prefix
            bot_oligo = _reverse_complement(sgRNA_prefix) + _reverse_complement(protospacer)
            insert_sequence = protospacer + sgRNA_prefix
            end_structure = "BclI-compatible 5' GATC overhang and SwaI blunt end"
        else:
            top_oligo = top_oh + protospacer
            bot_oligo = bot_oh + _reverse_complement(protospacer)
            insert_sequence = protospacer
            end_structure = "two sticky-end overhangs"

        slug       = (spec.name if spec else "custom").lower().replace(" ", "_")
        top_name   = f"{slug}_guide_top"
        bot_name   = f"{slug}_guide_bottom"
        ins_name   = f"{slug}_annealed_guide_insert"
        final_name = construct_name or f"{slug}_guide_construct"

        backbone_name, backbone_seq = self._resolve_backbone(spec)
        notes = self._typeIIS_notes(spec, enz, g_prepended, original, protospacer)

        construction_file_inputs = {
            "construct_name": final_name,
            "assembly_strategy": "TypeIISOligoCloning",
            "backbone_name": backbone_name,
            "backbone_sequence": backbone_seq,
            "insert_name": ins_name,
            "insert_sequence": insert_sequence,
            # Type IIS annealed-oligo cloning uses dedicated top/bottom oligo
            # fields; reusing the insert primer fields creates duplicate parts.
            "insert_forward_primer_name": "",
            "insert_forward_primer_sequence": "",
            "insert_reverse_primer_name": "",
            "insert_reverse_primer_sequence": "",
            "top_oligo_name": top_name,
            "top_oligo_sequence": top_oligo,
            "bottom_oligo_name": bot_name,
            "bottom_oligo_sequence": bot_oligo,
            "top_overhang": top_oh,
            "bottom_overhang": bot_oh,
            "vector_forward_primer_name": "",
            "vector_forward_primer_sequence": "",
            "vector_reverse_primer_name": "",
            "vector_reverse_primer_sequence": "",
            "enzyme": "" if enz == "N/A" else enz,
            "cell_strain": spec.cell_strain if spec else "",
            "selection": spec.selection if spec else "",
            "temperature_c": 37,
            "notes": " ".join(notes) + " This insert is a guide RNA/crRNA protospacer oligo (annealed oligo cloning), not a gene or multi-fragment insert.",
        }

        all_citations = (spec.citations if spec else ()) + WORKFLOW_CITATIONS["TypeIISOligoCloning"]

        return {
            "status": "ready",
            "cloning_method": "TypeIISOligoCloning",
            "vector": spec.name if spec else "custom",
            "enzyme": enz,
            "citations": _format_citations(all_citations),
            "top_overhang": top_oh,
            "bottom_overhang": bot_oh,
            "top_oligo_name": top_name,
            "bottom_oligo_name": bot_name,
            "top_oligo": top_oligo,
            "bottom_oligo": bot_oligo,
            "end_structure": end_structure,
            "g_prepended": g_prepended,
            "final_protospacer": protospacer,
            "construction_file_inputs": construction_file_inputs,
        }

    # ── Workflow branch B: Standard Restriction Ligation ─────────────────────

    def design_restriction_insert(
        self,
        guide_cassette_sequence: str,
        spec: Optional[VectorSpec],
        enzyme_override: Optional[str] = None,
        restriction_site_override: Optional[str] = None,
        construct_name: Optional[str] = None,
    ) -> dict:
        """
        Design PCR primers for standard restriction-enzyme-based guide cloning.

        This is biologically distinct from Type IIS cloning:
          • The RE site is added to PCR primer 5′ tails, not encoded in the backbone.
          • After digest, sticky ends are fixed by the enzyme's cut pattern —
            independent of the guide sequence.
          • The enzyme recognition site is disrupted at the ligation junction.

        Wet-lab steps simulated
        -----------------------
        1. Design PCR primers:
               fwd = [4-nt clamp][RE site][first 20 nt of cassette]
               rev = [4-nt clamp][RE site][RC of last 20 nt of cassette]
        2. PCR-amplify the guide cassette.
        3. Digest PCR product + vector with RE.
        4. Gel-purify; dephosphorylate vector (CIP).
        5. Ligate (T4, 16 °C overnight); transform; screen.

        SpeI / XbaI note: both leave a 5′-CTAG-3′ overhang. Their ligation
        product is not re-cuttable by either enzyme alone.

        construction_file_inputs uses assembly_strategy="DirectSynthesis" because
        create_construction_file has no RestrictionLigation strategy and the
        Gibson/GoldenGate strategies require vector PCR primers not designed here.
        """
        re_site = (
            restriction_site_override
            or (spec.restriction_site_sequence if spec else "")
        )
        enz = enzyme_override or (spec.enzyme if spec else "")

        if not re_site:
            raise ValueError(
                "restriction_site_sequence is required for RestrictionLigation design "
                "but was not provided and is not set on the vector spec."
            )

        clamp = "AAAA"  # 4-nt clamp improves RE digestion efficiency near PCR termini
        fwd_primer = clamp + re_site + guide_cassette_sequence[:20]
        rv_primer  = clamp + re_site + _reverse_complement(guide_cassette_sequence[-20:])

        slug       = (spec.name if spec else "custom").lower().replace(" ", "_")
        fwd_name   = f"{slug}_fwd_re_primer"
        rv_name    = f"{slug}_rv_re_primer"
        final_name = construct_name or f"{slug}_re_construct"

        backbone_name, backbone_seq = self._resolve_backbone(spec)

        note_lines = [
            f"Restriction ligation into {spec.name if spec else 'custom vector'} using {enz}.",
            f"RE site ({enz}: {re_site}) added to primer 5′ tails with {len(clamp)}-nt clamp.",
            f"Digest both PCR product and vector with {enz} before ligation.",
        ]

        construction_file_inputs = {
            "construct_name": final_name,
            "assembly_strategy": "RestrictionLigation",
            "backbone_name": backbone_name,
            "backbone_sequence": backbone_seq,
            "insert_name": f"{slug}_guide_cassette",
            "insert_sequence": guide_cassette_sequence,
            "insert_forward_primer_name": fwd_name,
            "insert_forward_primer_sequence": fwd_primer,
            "insert_reverse_primer_name": rv_name,
            "insert_reverse_primer_sequence": rv_primer,
            "vector_forward_primer_name": "",
            "vector_forward_primer_sequence": "",
            "vector_reverse_primer_name": "",
            "vector_reverse_primer_sequence": "",
            "enzyme": enz,
            "cell_strain": spec.cell_strain if spec else "",
            "selection": spec.selection if spec else "",
            "temperature_c": 37,
            "notes": " ".join(note_lines),
        }

        all_citations = (spec.citations if spec else ()) + WORKFLOW_CITATIONS["RestrictionLigation"]

        return {
            "status": "ready",
            "cloning_method": "RestrictionLigation",
            "vector": spec.name if spec else "custom",
            "enzyme": enz,
            "restriction_site": re_site,
            "citations": _format_citations(all_citations),
            "forward_primer_name": fwd_name,
            "reverse_primer_name": rv_name,
            "forward_primer": fwd_primer,
            "reverse_primer": rv_primer,
            "guide_cassette_sequence": guide_cassette_sequence,
            "construction_file_inputs": construction_file_inputs,
        }

    # ── Workflow branch C: Gibson Assembly ───────────────────────────────────

    def design_gibson_fragment(
        self,
        insert_sequence: str,
        left_overlap_context: str,
        right_overlap_context: str,
        spec: Optional[VectorSpec],
        overlap_bp: Optional[int] = None,
        construct_name: Optional[str] = None,
    ) -> dict:
        """
        Design a Gibson assembly insert with overlap regions.

        Gibson assembly is not restriction-enzyme cloning:
          • No sticky ends or RE scars in the final product.
          • Fragment order is set entirely by overlap sequences.
          • Multiple fragments can be joined in one tube.

        Inputs
        ------
        insert_sequence       — sequence to introduce (without overlap regions)
        left_overlap_context  — vector sequence immediately 5′ of insertion site
        right_overlap_context — vector sequence immediately 3′ of insertion site

        Outputs include:
          assembled_fragment   — left_overlap + insert + right_overlap
          overlap_forward_primer — left_overlap tail + first 20 nt of insert
          overlap_reverse_primer — right_overlap tail + RC of last 20 nt of insert

        construction_file_inputs uses assembly_strategy="DirectSynthesis" because
        the Gibson strategy in create_construction_file requires all four primer
        fields including vector linearisation primers this tool does not design.
        To use assembly_strategy="Gibson", supply vector linearisation primers
        separately and change the field before calling create_construction_file.

        Verify overlap Tm ≥ 55 °C and avoid repetitive sequences in the overlap.
        """
        obs = overlap_bp or (spec.recommended_overlap_bp if spec else 20)

        if len(left_overlap_context) < obs:
            raise ValueError(
                f"left_overlap_context ({len(left_overlap_context)} bp) is shorter "
                f"than the requested overlap ({obs} bp)."
            )
        if len(right_overlap_context) < obs:
            raise ValueError(
                f"right_overlap_context ({len(right_overlap_context)} bp) is shorter "
                f"than the requested overlap ({obs} bp)."
            )

        left_overlap  = left_overlap_context[-obs:]
        right_overlap = right_overlap_context[:obs]
        assembled     = left_overlap + insert_sequence + right_overlap

        fwd_primer = left_overlap  + insert_sequence[:20]
        rv_primer  = right_overlap + _reverse_complement(insert_sequence[-20:])

        slug       = (spec.name if spec else "custom").lower().replace(" ", "_")
        fwd_name   = f"{slug}_gibson_fwd"
        rv_name    = f"{slug}_gibson_rv"
        final_name = construct_name or f"{slug}_gibson_construct"

        backbone_name, backbone_seq = self._resolve_backbone(spec)

        note_lines = [
            f"Gibson assembly into {spec.name if spec else 'custom vector'}.  "
            f"Overlap: {obs} bp.",
            f"Left overlap: {left_overlap}  Right overlap: {right_overlap}.",
            "Incubate with NEBuilder HiFi or Gibson Assembly kit at 50 °C, 60 min.",
            "Vector linearisation primers are NOT included — add them manually to "
            "use assembly_strategy='Gibson' with create_construction_file.",
        ]

        construction_file_inputs = {
            "construct_name": final_name,
            "assembly_strategy": "Gibson",
            "backbone_name": backbone_name,
            "backbone_sequence": backbone_seq,
            "insert_name": f"{slug}_gibson_fragment",
            "insert_sequence": assembled,
            "insert_forward_primer_name": fwd_name,
            "insert_forward_primer_sequence": fwd_primer,
            "insert_reverse_primer_name": rv_name,
            "insert_reverse_primer_sequence": rv_primer,
            "vector_forward_primer_name": "",
            "vector_forward_primer_sequence": "",
            "vector_reverse_primer_name": "",
            "vector_reverse_primer_sequence": "",
            "enzyme": "",
            "cell_strain": spec.cell_strain if spec else "",
            "selection": spec.selection if spec else "",
            "temperature_c": 50,
            "notes": " ".join(note_lines),
        }

        all_citations = (spec.citations if spec else ()) + WORKFLOW_CITATIONS["GibsonAssembly"]

        return {
            "status": "ready",
            "cloning_method": "GibsonAssembly",
            "vector": spec.name if spec else "custom",
            "enzyme": "Gibson (isothermal)",
            "overlap_bp": obs,
            "citations": _format_citations(all_citations),
            "left_overlap": left_overlap,
            "right_overlap": right_overlap,
            "assembled_fragment": assembled,
            "overlap_forward_primer_name": fwd_name,
            "overlap_reverse_primer_name": rv_name,
            "overlap_forward_primer": fwd_primer,
            "overlap_reverse_primer": rv_primer,
            "insert_sequence": insert_sequence,
            "construction_file_inputs": construction_file_inputs,
        }

    # ── Workflow branch D: Golden Gate Assembly ───────────────────────────────

    def design_golden_gate_fragment(
        self,
        insert_sequence: str,
        left_vector_context: str,
        right_vector_context: str,
        left_overhang: str,
        right_overhang: str,
        enzyme: str,
        spec: Optional[VectorSpec],
        construct_name: Optional[str] = None,
    ) -> dict:
        """
        Design all four primers for Golden Gate assembly of an insert into a vector.

        Primer design
        -------------
        For a Type IIS enzyme (e.g. BsaI: GGTCTCN↓) the recognition site is placed
        in each primer's 5′ tail, oriented so the enzyme cuts toward the junction and
        exposes the user-defined 4-nt overhang.  Primer structures:

          insert_fwd : 5'-[clamp][recognition][spacer][left_overhang][insert_first_20nt]-3'
                       → exposes 5′-left_overhang on insert top strand

          insert_rev : 5'-[clamp][recognition][spacer][RC(right_overhang)][RC(insert_last_20nt)]-3'
                       → exposes 5′-RC(right_overhang) on insert bottom strand

          vector_fwd : 5'-[clamp][recognition][spacer][right_overhang][right_context_first_20nt]-3'
                       → exposes 5′-right_overhang on vector top strand (right side)

          vector_rev : 5'-[clamp][recognition][spacer][RC(left_overhang)][RC(left_context_last_20nt)]-3'
                       → exposes 5′-RC(left_overhang) on vector bottom strand (left side)

        After digest, insert and vector PCR products have complementary 5′ overhangs at
        both junctions and ligate directionally in a single T4-ligase or GG-mix reaction.

        construction_file_inputs uses assembly_strategy="GoldenGate" because all four
        primer fields are provided — unlike the other branches that use DirectSynthesis.

        Wet-lab steps
        -------------
        1. PCR insert with insert_fwd / insert_rev.
        2. PCR vector backbone with vector_fwd / vector_rev (amplifies entire backbone).
        3. Digest both PCR products with the Type IIS enzyme.
        4. Ligate (T4 ligase 16 °C overnight, or single-tube NEB Golden Gate mix).
        5. Transform; select; confirm by sequencing across both junctions.

        Overhang uniqueness: in a two-fragment assembly (insert + linearised vector)
        only two overhangs are required.  For multi-fragment assemblies, ALL overhangs
        (and their reverse complements) must be unique — verify with the NEB Ligation
        Fidelity Viewer or a similar tool before ordering primers.
        """
        enz_spec = GOLDEN_GATE_ENZYME_SPECS.get(enzyme)
        if enz_spec is None:
            known = ", ".join(sorted(GOLDEN_GATE_ENZYME_SPECS))
            raise ValueError(
                f"Unknown Golden Gate enzyme '{enzyme}'. "
                f"Supported enzymes: {known}."
            )

        recognition = enz_spec["recognition"]
        spacer      = "A" * enz_spec["spacer_len"]
        clamp       = "AAAA"

        left_oh  = _validate_dna(left_overhang,  "left_overhang")
        right_oh = _validate_dna(right_overhang, "right_overhang")

        if len(left_oh) != 4:
            raise ValueError(f"left_overhang must be exactly 4 nt, got {len(left_oh)}.")
        if len(right_oh) != 4:
            raise ValueError(f"right_overhang must be exactly 4 nt, got {len(right_oh)}.")
        if left_oh == right_oh:
            raise ValueError("left_overhang and right_overhang must be different.")
        if len(left_vector_context) < 20:
            raise ValueError(
                f"left_vector_context must be at least 20 bp, got {len(left_vector_context)}."
            )
        if len(right_vector_context) < 20:
            raise ValueError(
                f"right_vector_context must be at least 20 bp, got {len(right_vector_context)}."
            )
        if len(insert_sequence) < 20:
            raise ValueError(
                f"insert_sequence must be at least 20 bp, got {len(insert_sequence)}."
            )

        tail        = clamp + recognition + spacer
        insert_fwd  = tail + left_oh  + insert_sequence[:20]
        insert_rev  = tail + _reverse_complement(right_oh) + _reverse_complement(insert_sequence[-20:])
        vector_fwd  = tail + right_oh + right_vector_context[:20]
        vector_rev  = tail + _reverse_complement(left_oh)  + _reverse_complement(left_vector_context[-20:])

        slug         = (spec.name if spec else "custom").lower().replace(" ", "_")
        ins_fwd_name = f"{slug}_gg_insert_fwd"
        ins_rev_name = f"{slug}_gg_insert_rev"
        vec_fwd_name = f"{slug}_gg_vector_fwd"
        vec_rev_name = f"{slug}_gg_vector_rev"
        final_name   = construct_name or f"{slug}_gg_construct"

        backbone_name, backbone_seq = self._resolve_backbone(spec)

        note_lines = [
            f"Golden Gate assembly using {enzyme} (recognition site: {recognition}; spacer: {len(spacer)} nt).",
            f"Left junction overhang: {left_oh}.  Right junction overhang: {right_oh}.",
            "PCR insert with insert_fwd/rev; PCR vector backbone with vector_fwd/rev.",
            f"Digest both PCR products with {enzyme}, then ligate (T4 or NEB Golden Gate Mix).",
            "For multi-fragment assemblies verify all overhang sequences are unique "
            "(use NEB Ligation Fidelity Viewer).",
        ]

        construction_file_inputs = {
            "construct_name":                 final_name,
            "assembly_strategy":              "GoldenGate",
            "backbone_name":                  backbone_name,
            "backbone_sequence":              backbone_seq,
            "insert_name":                    f"{slug}_insert",
            "insert_sequence":                insert_sequence,
            "insert_forward_primer_name":     ins_fwd_name,
            "insert_forward_primer_sequence": insert_fwd,
            "insert_reverse_primer_name":     ins_rev_name,
            "insert_reverse_primer_sequence": insert_rev,
            "vector_forward_primer_name":     vec_fwd_name,
            "vector_forward_primer_sequence": vector_fwd,
            "vector_reverse_primer_name":     vec_rev_name,
            "vector_reverse_primer_sequence": vector_rev,
            "enzyme":                         enzyme,
            "cell_strain":                    spec.cell_strain if spec else "",
            "selection":                      spec.selection   if spec else "",
            "temperature_c":                  37,
            "notes":                          " ".join(note_lines),
        }

        all_citations = (spec.citations if spec else ()) + WORKFLOW_CITATIONS["GoldenGateAssembly"]

        return {
            "status":                      "ready",
            "cloning_method":              "GoldenGateAssembly",
            "vector":                      spec.name if spec else "custom",
            "enzyme":                      enzyme,
            "recognition_site":            recognition,
            "left_overhang":               left_oh,
            "right_overhang":              right_oh,
            "citations":                   _format_citations(all_citations),
            "insert_forward_primer_name":  ins_fwd_name,
            "insert_reverse_primer_name":  ins_rev_name,
            "vector_forward_primer_name":  vec_fwd_name,
            "vector_reverse_primer_name":  vec_rev_name,
            "insert_forward_primer":       insert_fwd,
            "insert_reverse_primer":       insert_rev,
            "vector_forward_primer":       vector_fwd,
            "vector_reverse_primer":       vector_rev,
            "construction_file_inputs":    construction_file_inputs,
        }

    # ── Public dispatcher ─────────────────────────────────────────────────────

    def run(
        self,
        # ── universal ─────────────────────────────────────────────────────
        vector: Optional[str] = None,
        cloning_method: Optional[str] = None,
        construct_name: Optional[str] = None,
        target_organism: Optional[str] = None,
        # ── TypeIIS branch ────────────────────────────────────────────────
        protospacer: Optional[str] = None,
        top_overhang: Optional[str] = None,
        bottom_overhang: Optional[str] = None,
        prepend_g: Optional[bool] = None,
        enzyme: Optional[str] = None,
        scaffold_in_vector: Optional[bool] = None,
        promoter: Optional[str] = None,
        u6_prefers_5prime_g: Optional[bool] = None,
        # ── RestrictionLigation branch ────────────────────────────────────
        guide_cassette_sequence: Optional[str] = None,
        restriction_site_sequence: Optional[str] = None,
        # ── GibsonAssembly branch ─────────────────────────────────────────
        insert_sequence: Optional[str] = None,
        left_overlap_context: Optional[str] = None,
        right_overlap_context: Optional[str] = None,
        overlap_bp: Optional[int] = None,
        # ── GoldenGateAssembly branch ─────────────────────────────────────
        left_overhang: Optional[str] = None,
        right_overhang: Optional[str] = None,
    ) -> dict:
        """
        Public entry point.

        Returns one of:
          - A "needs_user_input" dict if required information is missing.
          - A design result dict (with status="ready") if all inputs are present.

        Dispatch is by the cloning_method encoded in the resolved VectorSpec,
        or by the cloning_method argument when vector="custom".

        For known vectors: pass a vector name and the workflow-specific inputs.
        For custom vectors: pass vector="custom", cloning_method, and all
          workflow-specific parameters for your vector.

        """
        # Addgene fallback: if the vector string looks like an Addgene ID and is
        # not in the local preset registry, fetch metadata from the Addgene API
        # and try to build a VectorSpec from it.
        addgene_id = _parse_addgene_id(vector) if vector else None
        if addgene_id is not None and (not vector or vector.lower().strip() not in VECTOR_SPECS):
            fetcher = _AddgeneFetcher()
            fetcher.initiate()
            fetch_result = fetcher.run(addgene_id)
            if fetch_result.get("status") == "needs_user_input":
                return fetch_result  # API key missing — surface to user
            addgene_spec, needs_input = _addgene_data_to_vector_spec(fetch_result)
            if needs_input:
                return needs_input  # Overhangs unknown — ask user
            spec = addgene_spec
        else:
            # Resolve vector spec (None for custom)
            spec = self.resolve_vector(vector)

        # Organism compatibility guard
        if spec is not None and target_organism:
            compat_check = self._check_organism_compatibility(spec, target_organism)
            if compat_check:
                return compat_check

        # Determine effective cloning method
        method = cloning_method or (spec.cloning_method if spec else None)

        # Auto-assemble guide cassette for RestrictionLigation vectors where
        # scaffold_in_vector=False and promoter/nuclease sequences are known.
        if (
            method == "RestrictionLigation"
            and not guide_cassette_sequence
            and protospacer
            and spec is not None
            and not spec.scaffold_in_vector
        ):
            promoter_seq  = PROMOTER_SEQUENCES.get(spec.promoter, "")
            scaffold_seq  = SCAFFOLD_SEQUENCES.get(spec.nuclease_system, "")
            if promoter_seq and scaffold_seq:
                guide_cassette_sequence = promoter_seq + protospacer + scaffold_seq

        # For Gibson presets with known flanking sequences, auto-fill overlap
        # contexts and insert sequence so the user doesn't have to supply them.
        if method == "GibsonAssembly" and spec is not None:
            if not insert_sequence and protospacer and spec.scaffold_in_vector:
                insert_sequence = protospacer
            if not left_overlap_context and spec.gibson_left_context:
                left_overlap_context = spec.gibson_left_context
            if not right_overlap_context and spec.gibson_right_context:
                right_overlap_context = spec.gibson_right_context

        # Requirements check — may return early with needs_user_input
        check = self.assess_requirements(
            cloning_method=method,
            spec=spec,
            protospacer=protospacer,
            top_overhang=top_overhang,
            bottom_overhang=bottom_overhang,
            enzyme=enzyme,
            u6_prefers_5prime_g=u6_prefers_5prime_g,
            scaffold_in_vector=scaffold_in_vector,
            promoter=promoter,
            guide_cassette_sequence=guide_cassette_sequence,
            restriction_site_sequence=restriction_site_sequence,
            insert_sequence=insert_sequence,
            left_overlap_context=left_overlap_context,
            right_overlap_context=right_overlap_context,
            left_overhang=left_overhang,
            right_overhang=right_overhang,
        )
        if check["status"] == "needs_user_input":
            return check

        # All required inputs present — dispatch to design branch
        if method == "TypeIISOligoCloning":
            clean = _validate_dna(protospacer, "protospacer")
            return self.design_typeIIS_oligos(
                protospacer=clean,
                spec=spec,
                top_overhang_override=top_overhang,
                bottom_overhang_override=bottom_overhang,
                prepend_g_override=prepend_g,
                enzyme_override=enzyme,
                construct_name=construct_name,
            )

        if method == "RestrictionLigation":
            clean = _validate_dna(guide_cassette_sequence, "guide_cassette_sequence")
            return self.design_restriction_insert(
                guide_cassette_sequence=clean,
                spec=spec,
                enzyme_override=enzyme,
                restriction_site_override=restriction_site_sequence,
                construct_name=construct_name,
            )

        if method == "GibsonAssembly":
            return self.design_gibson_fragment(
                insert_sequence=_validate_dna(insert_sequence, "insert_sequence"),
                left_overlap_context=_validate_dna(left_overlap_context, "left_overlap_context"),
                right_overlap_context=_validate_dna(right_overlap_context, "right_overlap_context"),
                spec=spec,
                overlap_bp=overlap_bp,
                construct_name=construct_name,
            )

        if method == "GoldenGateAssembly":
            return self.design_golden_gate_fragment(
                insert_sequence=_validate_dna(insert_sequence, "insert_sequence"),
                left_vector_context=_validate_dna(left_overlap_context, "left_overlap_context"),
                right_vector_context=_validate_dna(right_overlap_context, "right_overlap_context"),
                left_overhang=_validate_dna(left_overhang, "left_overhang"),
                right_overhang=_validate_dna(right_overhang, "right_overhang"),
                enzyme=enzyme,
                spec=spec,
                construct_name=construct_name,
            )

        # Should not reach here after assess_requirements, but be explicit
        raise ValueError(
            f"Unrecognised cloning_method '{method}'. "
            "Expected: TypeIISOligoCloning, RestrictionLigation, GibsonAssembly, "
            "or GoldenGateAssembly."
        )

    # ── Private helpers ───────────────────────────────────────────────────────

    @staticmethod
    def _typeIIS_notes(
        spec: Optional[VectorSpec],
        enzyme: str,
        g_prepended: bool,
        original: str,
        final: str,
    ) -> list[str]:
        name  = spec.name if spec else "custom vector"
        notes = [f"Vector '{name}' — {enzyme} (TypeIIS) cloning."]
        if spec:
            notes.append(spec.notes)
        if g_prepended:
            notes.append(
                "5′G prepended for promoter compatibility — "
                "cloned guide is one base longer than the original protospacer."
            )
            notes.append(f"Original: {original}  Modified: {final}")
        else:
            notes.append("No 5′G added.")
        return notes


# ---------------------------------------------------------------------------
# Module-level singleton (harness compatibility)
# ---------------------------------------------------------------------------

_instance = CRISPRCloningDesigner()
_instance.initiate()
design_cloning_oligos = _instance.run


# ---------------------------------------------------------------------------
# Example usage
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("── 1. TypeIIS pX330 — complete inputs ──────────────────────────────")
    r = design_cloning_oligos(vector="px330", protospacer="ACAGAAACCTGCCAGTTTGC")
    print(f"  status       : {r['status']}")
    print(f"  top_oligo    : {r['top_oligo']}")
    print(f"  bottom_oligo : {r['bottom_oligo']}")
    print(f"  g_prepended  : {r['g_prepended']}")
    print()

    print("── 2. TypeIIS — missing protospacer → needs_user_input ─────────────")
    r2 = design_cloning_oligos(vector="px330")
    print(f"  status          : {r2['status']}")
    print(f"  missing_fields  : {r2['missing_fields']}")
    print(f"  question        : {r2['questions'][0]}")
    print()

    print("── 3. Custom TypeIIS — missing overhangs → needs_user_input ────────")
    r3 = design_cloning_oligos(
        vector="custom",
        cloning_method="TypeIISOligoCloning",
        protospacer="ACAGAAACCTGCCAGTTTGC",
    )
    print(f"  status          : {r3['status']}")
    print(f"  missing_fields  : {r3['missing_fields']}")
    print()

    print("── 4. Custom TypeIIS — all inputs provided ─────────────────────────")
    r4 = design_cloning_oligos(
        vector="custom",
        cloning_method="TypeIISOligoCloning",
        protospacer="ACAGAAACCTGCCAGTTTGC",
        top_overhang="CACC",
        bottom_overhang="AAAC",
        enzyme="BbsI",
        scaffold_in_vector=True,
        promoter="U6",
        u6_prefers_5prime_g=True,
    )
    print(f"  status       : {r4['status']}")
    print(f"  top_oligo    : {r4['top_oligo']}")
    print()

    print("── 5. RestrictionLigation pTargetF ─────────────────────────────────")
    cassette = (
        "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCGTTTTACAGAAACCTGCCAGTTTGC"
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAG"
    )
    r5 = design_cloning_oligos(
        vector="ptargetf",
        guide_cassette_sequence=cassette,
    )
    print(f"  status          : {r5['status']}")
    print(f"  forward_primer  : {r5['forward_primer']}")
    print(f"  restriction_site: {r5['restriction_site']}")
    print()

    print("── 6. GibsonAssembly — missing contexts → needs_user_input ─────────")
    r6 = design_cloning_oligos(
        vector="px458_gibson",
        insert_sequence="ACAGAAACCTGCCAGTTTGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG",
    )
    print(f"  status         : {r6['status']}")
    print(f"  missing_fields : {r6['missing_fields']}")
    print()

    print("── 7. GibsonAssembly — all inputs ──────────────────────────────────")
    r7 = design_cloning_oligos(
        vector="px458_gibson",
        insert_sequence="ACAGAAACCTGCCAGTTTGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG",
        left_overlap_context="GGCCGGCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTT",
        right_overlap_context="GCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTA",
    )
    print(f"  status         : {r7['status']}")
    print(f"  left_overlap   : {r7['left_overlap']}")
    print(f"  overlap_bp     : {r7['overlap_bp']}")
    print()
