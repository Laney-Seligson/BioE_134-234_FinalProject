import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.design_cas9_grna import DesignCas9Grna

_SEQ_WITH_NGG = (
    "CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTT"
    "TAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
)
_SEQ_NO_NGG = "ATGATGATGATGATGATGATGATG"


@pytest.fixture
def designer():
    d = DesignCas9Grna()
    d.initiate()
    return d


def test_returns_list_of_guides(designer):
    guides = designer.run(_SEQ_WITH_NGG)
    assert isinstance(guides, list)
    assert len(guides) > 0


def test_protospacer_is_20_bp(designer):
    guides = designer.run(_SEQ_WITH_NGG)
    for g in guides:
        assert len(g["protospacer"]) == 20


def test_pam_is_ngg(designer):
    guides = designer.run(_SEQ_WITH_NGG)
    for g in guides:
        assert g["pam_site"][1] == "G" and g["pam_site"][2] == "G"


def test_guide_contains_protospacer(designer):
    guides = designer.run(_SEQ_WITH_NGG)
    for g in guides:
        assert g["protospacer"].replace("U", "T") in g["grna_sequence"].replace("U", "T")


def test_empty_sequence_raises(designer):
    with pytest.raises(ValueError):
        designer.run("")


def test_no_ngg_pam_raises(designer):
    with pytest.raises(ValueError):
        designer.run(_SEQ_NO_NGG)


def test_returns_at_most_10_guides(designer):
    long_seq = _SEQ_WITH_NGG * 10
    guides = designer.run(long_seq)
    assert len(guides) <= 10
