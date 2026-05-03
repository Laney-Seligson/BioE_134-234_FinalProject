import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.rank_guides import rank_guides


def test_rank_guides_scores_single_guide_cas9():
    # GC=50%, no TTTT → efficiency max (2); on-target → specificity 1
    protospacer = "ATGCATGCATGCATGCATGG"
    reference = protospacer + "AGG" + "C" * 60
    result = rank_guides(
        guides=[{"protospacer": protospacer}],
        reference=reference,
        nuclease="cas9",
    )
    top = result["ranked_guides"][0]
    assert top["efficiency_score"] == 2
    assert top["specificity_score"] == 1
    assert top["total_score"] == 3
    assert top["efficiency_details"]["gc_content_ok"] is True
    assert top["efficiency_details"]["no_polyt_run"] is True
    assert "g_at_pam_proximal" not in top["efficiency_details"]
    assert result["best_guide"] is top
    assert "scoring_rationale" in result


def test_rank_guides_penalizes_polyt_and_low_gc():
    bad  = "AAAAAAAAAAAAAAAAAAAA"   # GC=0%, fails
    good = "ATGCATGCATGCATGCATGG"  # GC=50%, passes
    reference = good + "AGG" + "C" * 40 + bad + "AGG" + "C" * 40
    result = rank_guides(
        guides=[{"protospacer": bad}, {"protospacer": good}],
        reference=reference,
        nuclease="cas9",
    )
    assert result["best_guide"]["protospacer"] == good
    assert result["ranked_guides"][0]["efficiency_score"] >= result["ranked_guides"][1]["efficiency_score"]


def test_rank_guides_empty_list_raises():
    with pytest.raises(ValueError):
        rank_guides(guides=[], reference="ATGATGATG", nuclease="cas9")


def test_rank_guides_invalid_nuclease_raises():
    with pytest.raises(ValueError):
        rank_guides(
            guides=[{"protospacer": "ATGATGATGATGATGATGAG"}],
            reference="ATGATGATGATGATGATGAGAGG",
            nuclease="cas13",
        )


def test_rank_guides_missing_protospacer_raises():
    with pytest.raises(ValueError):
        rank_guides(
            guides=[{"not_a_protospacer": "ATGATGATGATGATGATGAG"}],
            reference="ATGATGATGATGATGATGAGAGG",
            nuclease="cas9",
        )


def test_rank_guides_cas12a_max_efficiency_is_two():
    protospacer = "ATGCATGCATGCATGCATGCATG"  # 23 bp, GC=52%, no TTTT
    reference = "TTTA" + protospacer + "A" * 40
    result = rank_guides(
        guides=[{"protospacer": protospacer}],
        reference=reference,
        nuclease="cas12a",
    )
    top = result["ranked_guides"][0]
    assert top["efficiency_score"] <= 2
    assert top["specificity_score"] == 1
    assert top["total_score"] == 3
    assert "g_at_pam_proximal" not in top["efficiency_details"]


def test_rank_guides_sorted_best_first():
    bad  = "AAAAAAAAAAAAAAAAAAAA"
    good = "ATGCATGCATGCATGCATGG"
    reference = good + "AGG" + "C" * 40 + bad + "AGG" + "C" * 40
    result = rank_guides(
        guides=[{"protospacer": bad}, {"protospacer": good}],
        reference=reference,
        nuclease="cas9",
    )
    scores = [g["total_score"] for g in result["ranked_guides"]]
    assert scores == sorted(scores, reverse=True)


def test_rank_guides_emits_citations():
    result = rank_guides(
        guides=[{"protospacer": "ATGCATGCATGCATGCATGG"}],
        reference="ATGCATGCATGCATGCATGGAGG" + "C" * 50,
        nuclease="cas9",
    )
    assert "citations" in result
    assert isinstance(result["citations"], list)
    assert len(result["citations"]) > 0
    for c in result["citations"]:
        assert set(c.keys()) == {"label", "reference", "claim"}
