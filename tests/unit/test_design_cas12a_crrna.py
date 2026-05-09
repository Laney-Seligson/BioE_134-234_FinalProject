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
    # Zetsche et al. 2015, Cell doi:10.1016/j.cell.2015.09.038
    # Cas12a scans the target strand for TTTV PAMs and returns each downstream 23 nt as a crRNA candidate.
    guides = designer.run(_SEQ_WITH_TTTV)
    assert isinstance(guides, list)
    assert len(guides) > 0


def test_protospacer_is_23_bp(designer):
    # Zetsche et al. 2015, Cell doi:10.1016/j.cell.2015.09.038
    # Cas12a cleaves with a 23 nt protospacer, 3 nt longer than SpCas9, contributing to its specificity profile.
    guides = designer.run(_SEQ_WITH_TTTV)
    for g in guides:
        assert len(g["protospacer"]) == 23


def test_pam_is_tttv(designer):
    # Zetsche et al. 2015, Cell doi:10.1016/j.cell.2015.09.038
    # Cas12a requires a 5'-TTTV PAM (V = A, C, or G) on the non-template strand, upstream of the protospacer.
    guides = designer.run(_SEQ_WITH_TTTV)
    for g in guides:
        pam = g["pam_site"]
        assert pam[:3] == "TTT"
        assert pam[3] in "ACG"  # V = not T


def test_empty_sequence_raises(designer):
    # No PAM sites can exist in an empty sequence; raising ValueError enforces valid input at the boundary.
    with pytest.raises(ValueError):
        designer.run("")


def test_no_tttv_pam_raises(designer):
    # Zetsche et al. 2015, Cell doi:10.1016/j.cell.2015.09.038
    # Without a TTTV PAM, Cas12a cannot bind; no valid crRNA can be returned.
    with pytest.raises(ValueError):
        designer.run(_SEQ_NO_TTTV)


def test_returns_at_most_10_guides(designer):
    # Zetsche et al. 2015, Cell doi:10.1016/j.cell.2015.09.038
    # Returning the top-ranked subset (≤10) mirrors standard crRNA selection practice for Cas12a experiments.
    long_seq = _SEQ_WITH_TTTV * 10
    guides = designer.run(long_seq)
    assert len(guides) <= 10
