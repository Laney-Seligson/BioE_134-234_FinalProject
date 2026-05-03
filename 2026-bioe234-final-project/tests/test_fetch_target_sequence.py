import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.crispr_tools.tools.fetch_target_sequence import FetchTargetSequence

_RAW_DNA = "ATGCATGCATGCATGCATGCATGC"


@pytest.fixture
def fetcher():
    f = FetchTargetSequence()
    f.initiate()
    return f


def test_raw_dna_returned_as_is(fetcher):
    result = fetcher.run(query=_RAW_DNA)
    assert result["sequence"] == _RAW_DNA
    assert result["source"] == "raw_input"


def test_raw_dna_is_uppercase(fetcher):
    result = fetcher.run(query=_RAW_DNA.lower())
    assert result["sequence"] == _RAW_DNA


def test_invalid_raw_dna_raises(fetcher):
    with pytest.raises(ValueError):
        fetcher.run(query="ATGCXYZ")


def test_empty_query_raises(fetcher):
    with pytest.raises(ValueError):
        fetcher.run(query="")


def test_result_has_required_keys(fetcher):
    result = fetcher.run(query=_RAW_DNA)
    for key in ("sequence", "source"):
        assert key in result
