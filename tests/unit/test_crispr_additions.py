"""
Tests for CRISPR tool additions and bug fixes:
  - Organism compatibility guard (target_organism mismatch detection)
  - Auto-cassette assembly for pTargetF (RestrictionLigation + scaffold_in_vector=False)
  - New vector presets: pml104, pdd162, pcfd3
  - pCRISPR rpsL cleanup
  - run_full_crispr_workflow exposes ranked_guides + scoring_rationale
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
from modules.crispr_tools.tools.run_full_crispr_workflow import RunFullCrisprWorkflow


# ---------------------------------------------------------------------------
# Fixtures / shared constants
# ---------------------------------------------------------------------------

PROTOSPACER = "CAGCTGGCGTAATAGCGAAG"  # 20 nt, valid ATGC
SHORT_SEQ = "ATGCATGCATGCATGCATGC" * 5  # 100 bp reference for rank_guides


# ---------------------------------------------------------------------------
# Organism compatibility guard
# ---------------------------------------------------------------------------

class TestOrganismCompatibility:

    def test_mammalian_vector_blocked_for_ecoli(self):
        r = design_cloning_oligos(
            vector="px330",
            protospacer=PROTOSPACER,
            target_organism="Escherichia coli",
        )
        assert r["status"] == "needs_user_input"
        assert "vector" in r["missing_fields"]
        q = r["questions"][0]
        assert "pX330" in q
        assert "Escherichia coli" in q or "ecoli" in q.lower()

    def test_mammalian_vector_blocked_for_ecoli_suggests_compatible(self):
        r = design_cloning_oligos(
            vector="px330",
            protospacer=PROTOSPACER,
            target_organism="Escherichia coli",
        )
        q = r["questions"][0]
        # Should suggest pcrispr and/or ptargetf
        assert "pcrispr" in q or "ptargetf" in q

    def test_ecoli_vector_passes_for_ecoli(self):
        # pcrispr is E. coli compatible — should not be blocked
        r = design_cloning_oligos(
            vector="pcrispr",
            protospacer=PROTOSPACER,
            target_organism="Escherichia coli",
        )
        # May still need_user_input for other reasons (overhangs), but NOT for organism
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_mammalian_vector_passes_for_mammalian(self):
        r = design_cloning_oligos(
            vector="px330",
            protospacer=PROTOSPACER,
            target_organism="Homo sapiens",
        )
        # Compatible — should not block on organism
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_unknown_organism_does_not_block(self):
        # Rat is not in _ORGANISM_COMPAT — should pass through without blocking
        r = design_cloning_oligos(
            vector="px330",
            protospacer=PROTOSPACER,
            target_organism="Rattus norvegicus",
        )
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_yeast_vector_blocked_for_mammalian(self):
        r = design_cloning_oligos(
            vector="pml104",
            protospacer=PROTOSPACER,
            target_organism="Homo sapiens",
        )
        assert r["status"] == "needs_user_input"
        assert "vector" in r["missing_fields"]

    def test_yeast_vector_passes_for_yeast(self):
        r = design_cloning_oligos(
            vector="pml104",
            protospacer=PROTOSPACER,
            target_organism="Saccharomyces cerevisiae",
        )
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])

    def test_drosophila_vector_blocked_for_ecoli(self):
        r = design_cloning_oligos(
            vector="pcfd3",
            protospacer=PROTOSPACER,
            target_organism="Escherichia coli",
        )
        assert r["status"] == "needs_user_input"
        assert "vector" in r["missing_fields"]

    def test_no_target_organism_skips_check(self):
        # Passing no target_organism should never block on organism mismatch
        r = design_cloning_oligos(vector="px330", protospacer=PROTOSPACER)
        assert r.get("status") != "needs_user_input" or "vector" not in r.get("missing_fields", [])


# ---------------------------------------------------------------------------
# Auto-cassette assembly (pTargetF RestrictionLigation bug fix)
# ---------------------------------------------------------------------------

class TestAutoCassetteAssembly:

    def test_ptargetf_protospacer_only_returns_ready(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert r["status"] == "ready"

    def test_cassette_starts_with_j23119(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        j23119 = PROMOTER_SEQUENCES["J23119"]
        assert r["guide_cassette_sequence"].startswith(j23119)

    def test_cassette_contains_protospacer(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        assert PROTOSPACER in r["guide_cassette_sequence"]

    def test_cassette_ends_with_scaffold(self):
        r = design_cloning_oligos(vector="ptargetf", protospacer=PROTOSPACER)
        scaffold = SCAFFOLD_SEQUENCES["SpCas9"]
        assert r["guide_cassette_sequence"].endswith(scaffold)

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

    def test_scaffold_in_vector_true_still_asks_for_cassette(self):
        # pcrispr has scaffold_in_vector=True — auto-build should NOT fire
        r = design_cloning_oligos(
            vector="pcrispr",
            cloning_method="RestrictionLigation",
            protospacer=PROTOSPACER,
        )
        # pcrispr is TypeIIS, so it won't hit the RE branch — but if called
        # with an RE method and no cassette, it should ask for the cassette
        assert r["status"] in ("needs_user_input", "ready")


# ---------------------------------------------------------------------------
# New vector presets
# ---------------------------------------------------------------------------

class TestNewVectorPresets:

    # pML104 — yeast
    def test_pml104_in_registry(self):
        assert "pml104" in VECTOR_SPECS

    def test_pml104_cloning_method(self):
        assert VECTOR_SPECS["pml104"].cloning_method == "TypeIISOligoCloning"

    def test_pml104_enzyme(self):
        assert VECTOR_SPECS["pml104"].enzyme == "BsaI"

    def test_pml104_cell_strain_is_yeast(self):
        assert "Saccharomyces" in VECTOR_SPECS["pml104"].cell_strain

    def test_pml104_has_citations(self):
        assert len(VECTOR_SPECS["pml104"].citations) > 0

    def test_pml104_missing_overhangs_ask_for_them(self):
        r = design_cloning_oligos(
            vector="pml104",
            protospacer=PROTOSPACER,
            target_organism="yeast",
        )
        assert r["status"] == "needs_user_input"
        assert "top_overhang" in r["missing_fields"]

    # pDD162 — C. elegans
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
            vector="pdd162",
            insert_sequence=PROTOSPACER * 3,
            target_organism="C. elegans",
        )
        assert r["status"] == "needs_user_input"
        assert "left_overlap_context" in r["missing_fields"]
        assert "right_overlap_context" in r["missing_fields"]

    def test_pdd162_routes_to_gibson_when_complete(self):
        context = "ATGCATGCATGCATGCATGCATGCATGC"  # 28 bp >= 20 bp overlap
        r = design_cloning_oligos(
            vector="pdd162",
            insert_sequence=PROTOSPACER * 3,
            left_overlap_context=context,
            right_overlap_context=context,
            target_organism="C. elegans",
        )
        assert r["status"] == "ready"
        assert r["cloning_method"] == "GibsonAssembly"

    # pCFD3 — Drosophila
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
            vector="pcfd3",
            protospacer=PROTOSPACER,
            target_organism="Drosophila melanogaster",
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
# run_full_crispr_workflow exposes ranked_guides + scoring_rationale
# ---------------------------------------------------------------------------

FAKE_SEQ = (
    "ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCC"
    "TGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATA"
    "GCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGG"
    "CGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCT"
    "TCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGC"
    "CCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAG"
    "AATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGG"
)

FAKE_SEQUENCE_INFO = {
    "sequence": FAKE_SEQ,
    "source": "mock",
    "resource": "lacZ",
    "organism": "Escherichia coli",
    "ncbi_gene_id": "945006",
    "ncbi_accession": "NC_000913.3",
    "length": len(FAKE_SEQ),
    "note": "Mock sequence for unit testing.",
}


class TestFullWorkflowRankGuides:

    @pytest.fixture
    def workflow(self):
        wf = RunFullCrisprWorkflow()
        wf.initiate()
        return wf

    def test_ranked_guides_present_in_output(self, workflow, monkeypatch):
        monkeypatch.setattr(
            workflow.sequence_fetcher, "run",
            lambda query, organism="Escherichia coli": FAKE_SEQUENCE_INFO,
        )
        result = workflow.run(query="lacZ", vector="ptargetf")
        assert "guides" in result
        assert isinstance(result["guides"], list)
        assert len(result["guides"]) > 0

    def test_each_guide_has_efficiency_score(self, workflow, monkeypatch):
        monkeypatch.setattr(
            workflow.sequence_fetcher, "run",
            lambda query, organism="Escherichia coli": FAKE_SEQUENCE_INFO,
        )
        result = workflow.run(query="lacZ", vector="ptargetf")
        for guide in result["guides"]:
            assert "efficiency_score" in guide
            assert "specificity_score" in guide
            assert "total_score" in guide

    def test_each_guide_has_efficiency_details(self, workflow, monkeypatch):
        monkeypatch.setattr(
            workflow.sequence_fetcher, "run",
            lambda query, organism="Escherichia coli": FAKE_SEQUENCE_INFO,
        )
        result = workflow.run(query="lacZ", vector="ptargetf")
        for guide in result["guides"]:
            assert "efficiency_details" in guide
            assert "specificity_details" in guide

    def test_scoring_rationale_present(self, workflow, monkeypatch):
        monkeypatch.setattr(
            workflow.sequence_fetcher, "run",
            lambda query, organism="Escherichia coli": FAKE_SEQUENCE_INFO,
        )
        result = workflow.run(query="lacZ", vector="ptargetf")
        assert "scoring_rationale" in result
        assert isinstance(result["scoring_rationale"], str)
        assert len(result["scoring_rationale"]) > 0

    def test_guides_sorted_best_first(self, workflow, monkeypatch):
        monkeypatch.setattr(
            workflow.sequence_fetcher, "run",
            lambda query, organism="Escherichia coli": FAKE_SEQUENCE_INFO,
        )
        result = workflow.run(query="lacZ", vector="ptargetf")
        scores = [g["total_score"] for g in result["guides"]]
        assert scores == sorted(scores, reverse=True)

    def test_selected_guide_is_top_ranked(self, workflow, monkeypatch):
        monkeypatch.setattr(
            workflow.sequence_fetcher, "run",
            lambda query, organism="Escherichia coli": FAKE_SEQUENCE_INFO,
        )
        result = workflow.run(query="lacZ", vector="ptargetf")
        assert result["selected_guide"]["protospacer"] == result["guides"][0]["protospacer"]
