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
    # Jinek et al. 2012, Science doi:10.1126/science.1225829
    # SpCas9 scans the target strand for NGG PAMs and returns each upstream 20 nt as a candidate guide.
    guides = designer.run(_SEQ_WITH_NGG)
    assert isinstance(guides, list)
    assert len(guides) > 0


def test_protospacer_is_20_bp(designer):
    # Jinek et al. 2012, Science doi:10.1126/science.1225829
    # SpCas9 base-pairs with a 20 nt protospacer immediately upstream of the NGG PAM.
    guides = designer.run(_SEQ_WITH_NGG)
    for g in guides:
        assert len(g["protospacer"]) == 20


def test_pam_is_ngg(designer):
    # Jinek et al. 2012, Science doi:10.1126/science.1225829
    # SpCas9 activity requires an NGG PAM on the non-template strand immediately 3' of the protospacer.
    guides = designer.run(_SEQ_WITH_NGG)
    for g in guides:
        assert g["pam_site"][1] == "G" and g["pam_site"][2] == "G"


def test_guide_contains_protospacer(designer):
    # Cong et al. 2013, Science doi:10.1126/science.1231143
    # The sgRNA spacer sequence is identical to the protospacer; cleavage depends on Watson-Crick complementarity.
    guides = designer.run(_SEQ_WITH_NGG)
    for g in guides:
        assert g["protospacer"].replace("U", "T") in g["grna_sequence"].replace("U", "T")


def test_empty_sequence_raises(designer):
    # No PAM sites can exist in an empty sequence; raising ValueError enforces valid input at the boundary.
    with pytest.raises(ValueError):
        designer.run("")


def test_no_ngg_pam_raises(designer):
    # Jinek et al. 2012, Science doi:10.1126/science.1225829
    # Without an NGG PAM, SpCas9 cannot bind; no valid guide can be returned.
    with pytest.raises(ValueError):
        designer.run(_SEQ_NO_NGG)


def test_returns_at_most_10_guides(designer):
    # Mali et al. 2013, Science doi:10.1126/science.1232033
    # Returning the top-ranked subset (≤10) mirrors standard sgRNA selection practice.
    long_seq = _SEQ_WITH_NGG * 10
    guides = designer.run(long_seq)
    assert len(guides) <= 10
