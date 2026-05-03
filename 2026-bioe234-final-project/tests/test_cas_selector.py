import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.cas_selector import CasSelector

_GC_RICH  = "GCGCGCGCGCGG"    # many NGG PAMs → Cas9 wins on guide count
# TTTV PAM + valid 23-bp spacer; no GG anywhere → zero valid Cas9 guides
_AT_RICH  = "TTTACTGATGAACAGTGACGACAGCGT" * 3


@pytest.fixture
def selector():
    s = CasSelector()
    s.initiate()
    return s


def test_gc_rich_recommends_cas9(selector):
    result = selector.run(seq=_GC_RICH, repair_template=False)
    assert result["recommendation"] == "Cas9"


def test_at_rich_recommends_cas12a(selector):
    result = selector.run(seq=_AT_RICH, repair_template=False)
    assert result["recommendation"] == "Cas12a"


def test_hdr_does_not_change_recommendation(selector):
    result_no_hdr = selector.run(seq=_GC_RICH, repair_template=False)
    result_hdr    = selector.run(seq=_GC_RICH, repair_template=True)
    assert result_no_hdr["recommendation"] == result_hdr["recommendation"]


def test_multiplexing_overrides_guide_count(selector):
    result = selector.run(seq=_GC_RICH, repair_template=False, num_targets=3)
    assert result["recommendation"] == "Cas12a"


def test_result_has_rationale(selector):
    result = selector.run(seq=_GC_RICH, repair_template=False)
    assert "rationale" in result
    assert isinstance(result["rationale"], str)
    assert len(result["rationale"]) > 0


def test_empty_sequence_raises(selector):
    with pytest.raises(ValueError):
        selector.run(seq="", repair_template=False)


def test_zero_num_targets_raises(selector):
    with pytest.raises(ValueError):
        selector.run(seq=_GC_RICH, repair_template=False, num_targets=0)


def test_margin_threshold_below_one_raises(selector):
    with pytest.raises(ValueError):
        selector.run(seq=_GC_RICH, repair_template=False, margin_threshold=0.5)
