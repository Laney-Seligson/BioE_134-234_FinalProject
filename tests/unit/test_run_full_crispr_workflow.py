import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.run_full_crispr_workflow import RunFullCrisprWorkflow

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


@pytest.fixture
def workflow():
    wf = RunFullCrisprWorkflow()
    wf.initiate()
    return wf


def _mock_fetcher(info: dict):
    """Return a sequence_fetcher.run replacement that accepts target_type."""
    def _run(query, organism="Escherichia coli", target_type="genomic_locus"):
        return info
    return _run


def test_ranked_guides_present_in_output(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        confirmed=True, target_type="genomic_locus",
    )
    assert "guides" in result
    assert isinstance(result["guides"], list)
    assert len(result["guides"]) > 0


def test_each_guide_has_efficiency_score(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        confirmed=True, target_type="genomic_locus",
    )
    for guide in result["guides"]:
        assert "efficiency_score" in guide
        assert "specificity_score" in guide
        assert "total_score" in guide


def test_each_guide_has_efficiency_details(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        confirmed=True, target_type="genomic_locus",
    )
    for guide in result["guides"]:
        assert "efficiency_details" in guide
        assert "specificity_details" in guide


def test_scoring_rationale_present(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        confirmed=True, target_type="genomic_locus",
    )
    assert "scoring_rationale" in result
    assert isinstance(result["scoring_rationale"], str)
    assert len(result["scoring_rationale"]) > 0


def test_guides_sorted_best_first(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        confirmed=True, target_type="genomic_locus",
    )
    scores = [g["total_score"] for g in result["guides"]]
    assert scores == sorted(scores, reverse=True)


def test_selected_guide_is_top_ranked(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        confirmed=True, target_type="genomic_locus",
    )
    assert result["selected_guide"]["protospacer"] == result["guides"][0]["protospacer"]


def test_empty_query_raises(workflow):
    with pytest.raises(ValueError):
        workflow.run(query="", vector="ptargetf")


def test_gene_symbol_without_organism_asks_for_organism(workflow):
    result = workflow.run(query="CFTR", vector="px330", confirmed=True)
    assert result["status"] == "needs_user_input"
    assert "organism" in result["missing_fields"]
    assert result["continue_with"]["organism"] == "<selected_organism>"
    assert result["workflow_trace"] == []


def test_gene_query_without_target_type_asks_for_target_type(workflow):
    # Added after the target_type gate was introduced in a0157b2.
    # Without target_type the workflow must stop and ask before fetching.
    result = workflow.run(query="lacZ", organism="Escherichia coli", vector="ptargetf")
    assert result["status"] == "needs_user_input"
    assert "target_type" in result["missing_fields"]
    assert "continue_with" in result


def test_raw_dna_sequence_skips_target_type_gate(workflow, monkeypatch):
    # A raw DNA query is not a gene name, so target_type is not required.
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query=FAKE_SEQ, organism="Escherichia coli", vector="ptargetf", confirmed=True
    )
    assert result.get("missing_fields") != ["target_type"]


def test_upstream_selected_gene_requires_explicit_gene_confirmation(workflow):
    result = workflow.run(
        query="EGFR",
        organism="Homo sapiens",
        vector="px330",
        source_query="genes in Homo sapiens associated with lung cancer",
        confirmed=True,
    )
    assert result["status"] == "needs_user_input"
    assert "gene_confirmed" in result["missing_fields"]
    assert result["continue_with"]["query"] == "EGFR"
    assert result["continue_with"]["source_query"] == (
        "genes in Homo sapiens associated with lung cancer"
    )
    assert result["continue_with"]["gene_confirmed"] is True
    assert result["workflow_trace"] == []


def test_upstream_selected_gene_confirmation_then_continues(workflow, monkeypatch):
    monkeypatch.setattr(
        workflow.sequence_fetcher, "run",
        _mock_fetcher({**FAKE_SEQUENCE_INFO, "resource": "EGFR", "organism": "Homo sapiens"}),
    )
    result = workflow.run(
        query="EGFR",
        organism="Homo sapiens",
        vector="",
        source_query="genes in Homo sapiens associated with lung cancer",
        gene_confirmed=True,
        target_type="genomic_locus",
    )
    assert result["status"] == "needs_user_input"
    assert "vector" in result["missing_fields"]
    assert "cas_recommendation" in result
    assert result["continue_with"]["query"] == "EGFR"
    assert result["continue_with"]["gene_confirmed"] is True
    assert any(step["tool"] == "crispr_cas_selector" for step in result["workflow_trace"])


def test_missing_vector_returns_needs_user_input(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="",
        target_type="genomic_locus",
    )
    assert result["status"] == "needs_user_input"
    assert "vector" in result["missing_fields"]
    assert "cas_recommendation" in result
    assert result["continue_with"]["organism"] == "Escherichia coli"
    assert any(step["tool"] == "crispr_cas_selector" for step in result["workflow_trace"])


def test_missing_vector_for_celegans_recommends_pdd162_after_cas_selection(workflow, monkeypatch):
    monkeypatch.setattr(
        workflow.sequence_fetcher, "run",
        _mock_fetcher({
            **FAKE_SEQUENCE_INFO,
            "resource": "dna-2",
            "organism": "Caenorhabditis elegans",
        }),
    )
    result = workflow.run(
        query="dna-2",
        organism="Caenorhabditis elegans",
        vector="",
        target_type="genomic_locus",
    )
    assert result["status"] == "needs_user_input"
    assert "cas_recommendation" in result
    assert result["vector_recommendations"][0]["vector_key"] == "pdd162"


def test_confirmation_gate_happens_after_cas_selection(workflow, monkeypatch):
    monkeypatch.setattr(workflow.sequence_fetcher, "run", _mock_fetcher(FAKE_SEQUENCE_INFO))
    result = workflow.run(
        query="lacZ", organism="Escherichia coli", vector="ptargetf",
        target_type="genomic_locus",
    )
    assert result["status"] == "needs_user_input"
    assert "confirmed" in result["missing_fields"]
    assert "cas_recommendation" in result
    assert "recommended_system" in result
    assert result["continue_with"]["organism"] == "Escherichia coli"
    assert any(step["tool"] == "crispr_cas_selector" for step in result["workflow_trace"])
    assert "guides" not in result
