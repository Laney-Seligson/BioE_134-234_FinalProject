"""
Tests for design_cloning_oligos covering:
  - Organism compatibility guard
  - Auto-cassette assembly for pTargetF
  - New vector presets: pml104, pdd162, pcfd3
  - pCRISPR rpsL cleanup
  - Misc cloning method routing
"""

import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.design_cloning_oligos import (
    design_cloning_oligos,
    VECTOR_SPECS,
    PROMOTER_SEQUENCES,
    SCAFFOLD_SEQUENCES,
)

PROTOSPACER = "CAGCTGGCGTAATAGCGAAG"


# ---------------------------------------------------------------------------
# Organism compatibility guard
# ---------------------------------------------------------------------------

class TestOrganismCompatibility:

    def test_mammalian_vector_blocked_for_ecoli(self):
        r = design_cloning_oligos(
            vector="px330", protospacer=PROTOSPACER, target_organism="Escherichia coli"
        )
        assert r["status"] == "needs_user_input"
        assert "vector" in r["missing_fields"]
        q = r["questions"][0]
        assert "pX330" in q
        assert "Escherichia coli" in q or "ecoli" in q.lower()

    def test_mammalian_vector_blocked_for_ecoli_suggests_compatible(self):
        r = design_cloning_oligos(
            vector="px330", protospacer=PROTOSPACER, target_organism="Escherichia coli"
        )
        q = r["questions"][0]
        assert "pcrispr" in q or "ptargetf" in q

    def test_ecoli_vector_passes_for_ecoli(self):
        r = design_cloning_oligos(
            vector="pcrispr", protospacer=PROTOSPACER, target_organism="Escherichia coli"
        )
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_mammalian_vector_passes_for_mammalian(self):
        r = design_cloning_oligos(
            vector="px330", protospacer=PROTOSPACER, target_organism="Homo sapiens"
        )
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_unknown_organism_does_not_block(self):
        r = design_cloning_oligos(
            vector="px330", protospacer=PROTOSPACER, target_organism="Rattus norvegicus"
        )
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_yeast_vector_blocked_for_mammalian(self):
        r = design_cloning_oligos(
            vector="pml104", protospacer=PROTOSPACER, target_organism="Homo sapiens"
        )
        assert r["status"] == "needs_user_input"
        assert "vector" in r["missing_fields"]

    def test_yeast_vector_passes_for_yeast(self):
        r = design_cloning_oligos(
            vector="pml104", protospacer=PROTOSPACER, target_organism="Saccharomyces cerevisiae"
        )
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_drosophila_vector_blocked_for_ecoli(self):
        r = design_cloning_oligos(
            vector="pcfd3", protospacer=PROTOSPACER, target_organism="Escherichia coli"
        )
        assert r["status"] == "needs_user_input"
        assert "vector" in r["missing_fields"]

    def test_no_target_organism_skips_check(self):
        r = design_cloning_oligos(vector="px330", protospacer=PROTOSPACER)
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_px330_uses_real_backbone_sequence(self):
        r = design_cloning_oligos(vector="px330", protospacer=PROTOSPACER)
        inputs = r["construction_file_inputs"]
        assert inputs["backbone_name"] == "px330"
        assert inputs["backbone_sequence"] != "N"
        assert len(inputs["backbone_sequence"]) == 8484


# ---------------------------------------------------------------------------
# Auto-cassette assembly (pTargetF RestrictionLigation)
# ---------------------------------------------------------------------------

class TestAutoCassetteAssembly:

    def test_ptargetf_protospacer_only_returns_ready(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert r["status"] == "ready"

    def test_cassette_starts_with_j23119(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert r["guide_cassette_sequence"].startswith(PROMOTER_SEQUENCES["J23119"])

    def test_cassette_contains_protospacer(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert PROTOSPACER in r["guide_cassette_sequence"]

    def test_cassette_ends_with_scaffold(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert r["guide_cassette_sequence"].endswith(SCAFFOLD_SEQUENCES["SpCas9"])

    def test_cassette_length_is_promoter_plus_spacer_plus_scaffold(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        expected = (
            len(PROMOTER_SEQUENCES["J23119"])
            + len(PROTOSPACER)
            + len(SCAFFOLD_SEQUENCES["SpCas9"])
        )
        assert len(r["guide_cassette_sequence"]) == expected

    def test_forward_primer_contains_spei_site(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert "ACTAGT" in r["forward_primer"]

    def test_reverse_primer_contains_spei_site(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert "ACTAGT" in r["reverse_primer"]

    def test_scaffold_in_vector_true_still_routes_correctly(self):
        r = design_cloning_oligos(
            vector="pcrispr", cloning_method="RestrictionLigation", protospacer=PROTOSPACER
        )
        assert r["status"] in ("needs_user_input", "ready")


# ---------------------------------------------------------------------------
# New vector presets
# ---------------------------------------------------------------------------

class TestNewVectorPresets:

    def test_pml104_in_registry(self):
        assert "pml104" in VECTOR_SPECS

    def test_pml104_cloning_method(self):
        assert VECTOR_SPECS["pml104"].cloning_method == "TypeIISOligoCloning"

    def test_pml104_enzyme(self):
        assert VECTOR_SPECS["pml104"].enzyme == "BclI-SwaI"

    def test_pml104_cell_strain_is_yeast(self):
        assert "Saccharomyces" in VECTOR_SPECS["pml104"].cell_strain

    def test_pml104_has_citations(self):
        assert len(VECTOR_SPECS["pml104"].citations) > 0

    def test_pml104_design_uses_bcli_swai_blunt_insert(self):
        r = design_cloning_oligos(
            vector="pml104", protospacer=PROTOSPACER, target_organism="yeast"
        )
        assert r["status"] == "ready"
        assert r["top_overhang"] == "GATC"
        assert r["bottom_overhang"] == ""
        assert r["top_oligo"] == "GATCCAGCTGGCGTAATAGCGAAGGTTTTAGAGCTAG"
        assert r["bottom_oligo"] == "CTAGCTCTAAAACCTTCGCTATTACGCCAGCTG"
        assert r["end_structure"] == "BclI-compatible 5' GATC overhang and SwaI blunt end"

    def test_pml107_in_registry(self):
        assert "pml107" in VECTOR_SPECS

    def test_pml107_marker_is_leu2(self):
        assert "LEU2" in VECTOR_SPECS["pml107"].selection

    def test_pdd162_in_registry(self):
        assert "pdd162" in VECTOR_SPECS

    def test_pdd162_cloning_method_is_gibson(self):
        assert VECTOR_SPECS["pdd162"].cloning_method == "GibsonAssembly"

    def test_pdd162_enzyme_is_gibson(self):
        assert VECTOR_SPECS["pdd162"].enzyme == "Gibson"

    def test_pdd162_cell_strain_is_celegans(self):
        assert "elegans" in VECTOR_SPECS["pdd162"].cell_strain.lower()

    def test_pdd162_has_citations(self):
        assert len(VECTOR_SPECS["pdd162"].citations) > 0

    def test_pdd162_missing_contexts_asks_for_them(self):
        r = design_cloning_oligos(
            vector="pdd162", insert_sequence=PROTOSPACER * 3, target_organism="C. elegans"
        )
        assert r["status"] == "needs_user_input"
        assert "left_overlap_context" in r["missing_fields"]
        assert "right_overlap_context" in r["missing_fields"]

    def test_pdd162_routes_to_gibson_when_complete(self):
        context = "ATGCATGCATGCATGCATGCATGCATGC"  # >= 20 bp
        r = design_cloning_oligos(
            vector="pdd162",
            insert_sequence=PROTOSPACER * 3,
            left_overlap_context=context,
            right_overlap_context=context,
            target_organism="C. elegans",
        )
        assert r["status"] == "ready"
        assert r["cloning_method"] == "GibsonAssembly"

    def test_pcfd3_in_registry(self):
        assert "pcfd3" in VECTOR_SPECS

    def test_pcfd3_cloning_method(self):
        assert VECTOR_SPECS["pcfd3"].cloning_method == "TypeIISOligoCloning"

    def test_pcfd3_enzyme(self):
        assert VECTOR_SPECS["pcfd3"].enzyme == "BbsI"

    def test_pcfd3_cell_strain_is_drosophila(self):
        assert "Drosophila" in VECTOR_SPECS["pcfd3"].cell_strain

    def test_pcfd3_has_citations(self):
        assert len(VECTOR_SPECS["pcfd3"].citations) > 0

    def test_pcfd3_missing_overhangs_ask_for_them(self):
        r = design_cloning_oligos(
            vector="pcfd3", protospacer=PROTOSPACER, target_organism="Drosophila melanogaster"
        )
        assert r["status"] == "needs_user_input"
        assert "top_overhang" in r["missing_fields"]


# ---------------------------------------------------------------------------
# pCRISPR rpsL cleanup
# ---------------------------------------------------------------------------

class TestPCRISPRCleanup:

    def test_pcrispr_name_has_no_rpsl(self):
        assert "rpsL" not in VECTOR_SPECS["pcrispr"].name

    def test_pcrispr_name_is_pcrispr(self):
        assert VECTOR_SPECS["pcrispr"].name == "pCRISPR"

    def test_pcrispr_notes_have_no_rpsl(self):
        assert "rpsL" not in VECTOR_SPECS["pcrispr"].notes

    def test_pcrispr_citation_claims_have_no_rpsl(self):
        for citation in VECTOR_SPECS["pcrispr"].citations:
            assert "rpsL" not in citation.claim


# ---------------------------------------------------------------------------
# pET28a / pBR322 removed from VECTOR_SPECS
# ---------------------------------------------------------------------------

def test_pet28a_not_in_vector_specs():
    assert "pet28a" not in VECTOR_SPECS


def test_pbr322_not_in_vector_specs():
    assert "pbr322" not in VECTOR_SPECS


# ---------------------------------------------------------------------------
# pcrispr oligo output
# ---------------------------------------------------------------------------

def test_pcrispr_ecoli_oligo_output():
    result = design_cloning_oligos(
        protospacer="TACTTTACGCAGCGCGGAGT",
        vector="pcrispr",
        target_organism="E. coli",
    )
    assert result["vector"].lower() == "pcrispr"
    assert result["enzyme"] == "BsaI"
    assert any(c["label"].startswith("Jiang et al.") for c in result["citations"])
    assert result["top_oligo"] == "AAACTACTTTACGCAGCGCGGAGT"
    assert result["bottom_oligo"] == "AAAACACTCCGCGCTGCGTAAAGTA"
    assert result["construction_file_inputs"]["backbone_name"] == "pCRISPR_rpsL"
    assert result["construction_file_inputs"]["cell_strain"] == "HME63 or MG1655 carrying pCas9"
    assert result["construction_file_inputs"]["selection"] == "Kan"


def test_pcrispr_zebrafish_returns_needs_user_input():
    result = design_cloning_oligos(
        protospacer="CCGGATGCTCCTCAGCTCTG",
        vector="pcrispr",
        target_organism="zebrafish",
    )
    assert result["status"] == "needs_user_input"
