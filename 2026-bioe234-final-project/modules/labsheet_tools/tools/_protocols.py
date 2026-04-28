"""
Canonical wet-lab protocol registry.

Every recipe rendered into a lab sheet should resolve to one of these
entries. Volumes and conditions are transcribed from manufacturer manuals
or the original primary literature, never invented. Each entry carries a
`source` key that resolves into citations.CITATIONS, so the lab_sheet
output can ship a `protocol_sources` field that is auditable end-to-end.

If you need to update a recipe (e.g. NEB revises a kit insert), edit the
entry here — every downstream lab sheet picks up the change automatically.

Schema:
    Protocol.name        — human-readable label rendered into the sheet
    Protocol.reagents    — ordered list of (reagent, volume) tuples; the
                           sum of volumes is the total reaction volume
    Protocol.program     — thermocycler program / incubation conditions
    Protocol.source      — citation key (must exist in citations.CITATIONS)
    Protocol.notes       — optional protocol-specific notes
"""

from __future__ import annotations

from dataclasses import dataclass, field

from modules.crispr_tools.tools.citations import CITATIONS


@dataclass(frozen=True)
class Protocol:
    name: str
    reagents: tuple[tuple[str, str], ...]
    program: str
    source: str            # key into CITATIONS
    notes: str = ""
    extra_sources: tuple[str, ...] = field(default_factory=tuple)


# ---------------------------------------------------------------------------
# Registry. Six canonical protocols covering the steps lab_sheet emits.
# ---------------------------------------------------------------------------

PROTOCOLS: dict[str, Protocol] = {

    # Q5 PCR (NEB M0491) — standard 50 uL reaction
    "pcr_q5": Protocol(
        name="Q5 High-Fidelity PCR (50 uL)",
        reagents=(
            ("5X Q5 Reaction Buffer", "10 uL"),
            ("10 mM dNTPs", "1 uL"),
            ("10 uM forward primer", "2.5 uL"),
            ("10 uM reverse primer", "2.5 uL"),
            ("Template DNA (1-10 ng plasmid / 50-250 ng gDNA)", "1 uL"),
            ("Q5 High-Fidelity DNA Polymerase", "0.5 uL"),
            ("ddH2O", "to 50 uL"),
        ),
        program="98C 30s; [98C 10s, 60C* 30s, 72C 30s/kb] x 30; 72C 2min; 4C hold",
        source="neb_m0491_q5",
        notes="*Annealing temperature: use NEB Tm Calculator for primer pair. "
              "Extension 30 sec/kb (longer for amplicons >6kb).",
    ),

    # Golden Gate assembly (Engler 2008 + NEB E1601) — 10 uL reaction
    "goldengate_bsai": Protocol(
        name="Golden Gate Assembly (BsaI-HFv2, 10 uL)",
        reagents=(
            ("ddH2O", "7 uL"),
            ("T4 DNA Ligase Buffer (10X)", "1 uL"),
            ("DNA Mix (5 uL per fragment, ~75 ng each)", "1 uL"),
            ("T4 DNA Ligase", "0.5 uL"),
            ("BsaI-HFv2", "0.5 uL"),
        ),
        program="[37C 5min, 16C 5min] x 30; 60C 5min; 80C 5min; 4C hold",
        source="engler_2008_goldengate",
        extra_sources=("neb_e1601",),
        notes="Use BsmBI-v2 instead of BsaI-HFv2 for MoClo Level 2 reactions.",
    ),

    # Gibson Assembly (Gibson 2009 + NEB E2611) — 10 uL reaction
    "gibson_neb_2x": Protocol(
        name="Gibson Assembly (NEB 2X Master Mix, 10 uL)",
        reagents=(
            ("ddH2O", "4 uL"),
            ("DNA Mix (5 uL per fragment, 0.02-0.5 pmol each)", "1 uL"),
            ("Gibson Assembly Master Mix (2X)", "5 uL"),
        ),
        program="50C 15-60 min (15 min for 2-3 fragments, 60 min for 4-6); 4C hold",
        source="gibson_2009",
        extra_sources=("neb_e2611",),
        notes="For >6 fragments or fragments <200 bp, use NEBuilder HiFi instead.",
    ),

    # Chemical transformation (Inoue method + Sambrook protocols)
    "transformation_chemical": Protocol(
        name="Chemical Transformation (heat shock)",
        reagents=(
            ("Chemically competent E. coli", "50 uL"),
            ("Assembly product or plasmid (1-100 ng)", "1-5 uL"),
            ("SOC or 2YT recovery medium", "950 uL"),
        ),
        program="Mix on ice 30 min; 42C heat shock 30-45 sec; ice 2 min; "
                "add recovery medium; 37C shaking 45-60 min; plate on antibiotic.",
        source="inoue_1990",
        extra_sources=("sambrook_2001",),
        notes="Recovery time depends on antibiotic: Amp = 30 min sufficient; "
              "Kan/Spec/Cm = 60 min recommended for full marker expression.",
    ),

    # Plasmid miniprep (Zymo silica-column)
    "miniprep_zymo": Protocol(
        name="Plasmid Miniprep (Zymo silica column)",
        reagents=(
            ("Overnight bacterial culture", "600 uL"),
            ("7X Lysis Buffer (blue)", "100 uL"),
            ("Cold Neutralization Buffer (yellow)", "350 uL"),
            ("DNA Wash Buffer", "200 uL (x2)"),
            ("DNA Elution Buffer or ddH2O", "30-50 uL"),
        ),
        program="Lyse 1-2 min (color blue->yellow on neutralize); centrifuge "
                "16,000 x g 4 min; bind to column; wash 2x; elute warm.",
        source="zymo_d4015",
        notes="For low-copy plasmids (e.g. BAC, pSC101), use 5-10 mL culture "
              "and a midiprep kit instead.",
    ),

    # Type IIS oligo cloning (single-pot digestion-ligation of annealed
    # oligo duplex into a Type IIS-cut vector, e.g. pX330/pSpCas9 with BbsI
    # or pCRISPR/Lenti-Guide with BsmBI). Engler 2008 + NEB E1601 recipe.
    "typeiis_oligo_cloning": Protocol(
        name="Type IIS Oligo Cloning (anneal + digest-ligate, 10 uL)",
        reagents=(
            ("100 uM top oligo", "1 uL"),
            ("100 uM bottom oligo", "1 uL"),
            ("T4 DNA Ligase Buffer (10X) — also serves as anneal buffer", "1 uL"),
            ("Vector (50 ng, undigested)", "1 uL"),
            ("Type IIS enzyme (BbsI-HF / BsmBI-v2)", "0.5 uL"),
            ("T4 DNA Ligase", "0.5 uL"),
            ("ddH2O", "5 uL"),
        ),
        program=(
            "Phase 1 - Anneal oligos: combine 1 uL top + 1 uL bottom + 1 uL "
            "T4 ligase buffer + 7 uL ddH2O in a PCR tube. Heat 95C 5 min, "
            "ramp -0.1C/sec to 25C (or boil and slow-cool in a beaker). "
            "Dilute 1:200 in ddH2O before use.\n"
            "Phase 2 - Digest-ligate: in a fresh tube combine the recipe "
            "above using 1 uL of the diluted oligo duplex in place of the "
            "anneal mix. Run [37C 5 min, 16C 5 min] x 30; 60C 5 min; "
            "80C 5 min; 4C hold."
        ),
        source="engler_2008_goldengate",
        extra_sources=("neb_e1601",),
        notes=(
            "Use BbsI-HF (NEB R3539) for pX330/pSpCas9-style Cas9 vectors "
            "and BsmBI-v2 (NEB R0739) for Lenti-Guide-Puro / pCRISPR. "
            "Vector overhangs (e.g. CACC/AAAC) come from "
            "design_cloning_oligos. Resulting plasmid carries the protospacer "
            "in-frame with the sgRNA scaffold."
        ),
    ),

    # CRISPR RNP assembly (IDT Alt-R)
    "crispr_rnp_assembly": Protocol(
        name="Cas9 RNP Assembly (IDT Alt-R, 5 uL)",
        reagents=(
            ("Alt-R crRNA (100 uM)", "1 uL"),
            ("Alt-R tracrRNA (100 uM)", "1 uL"),
            ("Nuclease-Free Duplex Buffer", "1 uL"),
            ("Alt-R S.p. Cas9 Nuclease V3 (10 ug/uL = ~62 uM)", "1 uL"),
            ("PBS or Opti-MEM (delivery buffer)", "1 uL"),
        ),
        program="Anneal crRNA+tracrRNA: 95C 5 min, slow cool to RT (-1C/sec). "
                "Complex with Cas9 at 1:1.2 (RNA:protein) molar ratio: "
                "incubate at room temperature 20 min. Use within 1 hr.",
        source="idt_altr_rnp",
        notes="Final RNP concentration ~10 uM. Scale up linearly for "
              "multiple electroporations or larger lipofection volumes.",
    ),

    # CRISPR delivery — electroporation (Lonza Nucleofector)
    "crispr_electroporation": Protocol(
        name="RNP Electroporation (Lonza 4D Nucleofector, 20 uL cuvette)",
        reagents=(
            ("Pre-formed Cas9 RNP (10 uM, from crispr_rnp_assembly)", "5 uL"),
            ("Cells in P3/SF/SE Nucleofection Solution", "15 uL"),
            ("(target: 1e5 - 5e5 cells per cuvette)", ""),
        ),
        program="Resuspend cells in Lonza solution; add RNP; transfer to "
                "Nucleocuvette; run cell-line-specific program (e.g. CA-137 "
                "for K562, CM-138 for HEK293, EH-100 for primary T cells); "
                "transfer to pre-warmed media within 10 min.",
        source="lonza_nucleofection",
        notes="Always check Lonza's published optimization tables for the "
              "exact program for your cell line. RNP > plasmid for editing "
              "efficiency and lower off-target rate.",
    ),

    # CRISPR delivery — lipofection (Lipofectamine CRISPRMAX)
    "crispr_lipofection": Protocol(
        name="RNP Lipofection (Lipofectamine CRISPRMAX, 24-well)",
        reagents=(
            ("Cas9 RNP (7.5 pmol)", "1.5 uL of 5 uM RNP"),
            ("Cas9 Plus Reagent", "1 uL"),
            ("Opti-MEM (mix 1)", "25 uL"),
            ("Lipofectamine CRISPRMAX", "1.5 uL"),
            ("Opti-MEM (mix 2)", "25 uL"),
        ),
        program="Combine RNP + Plus Reagent + Opti-MEM (mix 1). In parallel, "
                "combine CRISPRMAX + Opti-MEM (mix 2). Incubate each 5 min RT. "
                "Combine 1:1; incubate 5 min RT. Add 50 uL complex to cells "
                "in 500 uL antibiotic-free media. Refresh media at 24 hr.",
        source="thermo_crisprmax",
        notes="Use for adherent cell lines that tolerate lipofection. For "
              "primary cells or hard-to-transfect lines, prefer "
              "crispr_electroporation instead.",
    ),

    # In-vitro transcription of sgRNA (NEB EnGen E3322)
    "crispr_ivt_sgrna": Protocol(
        name="sgRNA In Vitro Transcription (NEB EnGen E3322, 20 uL)",
        reagents=(
            ("Target-specific DNA oligo (1 uM, with T7 promoter + spacer)", "1 uL"),
            ("sgRNA Scaffold DNA Mix", "2 uL"),
            ("NTP Buffer Mix", "10 uL"),
            ("EnGen sgRNA Enzyme Mix (T7 RNAP + others)", "2 uL"),
            ("Nuclease-free water", "5 uL"),
        ),
        program="37C 30 min; add 1 uL DNase I, 37C 15 min; LiCl precipitate "
                "or column purify (Monarch RNA Cleanup). Quantify; store -80C.",
        source="neb_e3322_ivt",
        notes="Use IVT sgRNA when you need many guides cheaply (e.g. screening). "
              "Synthetic Alt-R guides (used in crispr_rnp_assembly) give higher "
              "and more reproducible editing in mammalian cells.",
    ),

    # Sanger sequencing premix (UC Berkeley standard)
    "sanger_sequencing": Protocol(
        name="Sanger Sequencing Premix Submission",
        reagents=(
            ("Miniprep plasmid DNA (200-500 ng)", "4 uL"),
            ("Sequencing primer (2.66 uM working stock)", "3 uL"),
            ("ddH2O", "3 uL"),
        ),
        program="No thermocycling — submit premix tubes to sequencing facility.",
        source="sanger_premix_standard",
        notes="Working primer stock prep: dilute 100 uM stock 1:37.6 in ddH2O "
              "(13.3 uL stock + 487 uL ddH2O => 2.66 uM, 500 uL). "
              "Total premix volume ~10 uL per reaction.",
    ),
}


def get(key: str) -> Protocol:
    """Look up a protocol by key. Raises KeyError if unknown — fail fast
    rather than silently emit a recipe with no source."""
    return PROTOCOLS[key]


def reagent_block(key: str) -> str:
    """Render a protocol's reagents as a multi-line string suitable for
    pasting into the lab sheet's plain-text view."""
    proto = PROTOCOLS[key]
    return "\n".join(f"{vol} {reagent}" for reagent, vol in proto.reagents)


def protocol_source_record(key: str) -> dict:
    """Render a protocol's source(s) as a dict suitable for inclusion in a
    lab sheet's `protocol_sources` field. Includes both the primary source
    and any extra_sources (e.g. NEB kit insert in addition to the original
    method paper)."""
    proto = PROTOCOLS[key]
    keys = (proto.source,) + tuple(proto.extra_sources)
    sources = []
    for k in keys:
        cit = CITATIONS[k]
        sources.append({
            "label": cit.label,
            "reference": cit.url_or_reference,
            "claim": cit.claim,
        })
    return {
        "step": proto.name,
        "sources": sources,
    }
