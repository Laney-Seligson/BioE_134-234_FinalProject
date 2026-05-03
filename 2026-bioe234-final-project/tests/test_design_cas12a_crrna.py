import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.design_cas12a_crrna import DesignCas12aCrrna

_SEQ_WITH_TTTV = (
    "CCCTAGATGCCTTTTAGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTT"
    "TAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
)
_SEQ_NO_TTTV = "ATGATGATGATGATGATGATGATG"


@pytest.fixture
def designer():
    d = DesignCas12aCrrna()
    d.initiate()
    return d


def test_returns_list_of_guides(designer):
    guides = designer.run(_SEQ_WITH_TTTV)
    assert isinstance(guides, list)
    assert len(guides) > 0


def test_protospacer_is_23_bp(designer):
    guides = designer.run(_SEQ_WITH_TTTV)
    for g in guides:
        assert len(g["protospacer"]) == 23


def test_pam_is_tttv(designer):
    guides = designer.run(_SEQ_WITH_TTTV)
    for g in guides:
        pam = g["pam_site"]
        assert pam[:3] == "TTT"
        assert pam[3] in "ACG"  # V = not T


def test_empty_sequence_raises(designer):
    with pytest.raises(ValueError):
        designer.run("")


def test_no_tttv_pam_raises(designer):
    with pytest.raises(ValueError):
        designer.run(_SEQ_NO_TTTV)


def test_returns_at_most_10_guides(designer):
    long_seq = _SEQ_WITH_TTTV * 10
    guides = designer.run(long_seq)
    assert len(guides) <= 10
