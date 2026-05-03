import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.construction_file_tools.tools.create_construction_file import create_construction_file
from modules.construction_file_tools.tools.validate_construction_file import validate_construction_record

VALID_CONSTRUCTION_INPUT = {
    "construct_name": "pET28a_REP24",
    "assembly_strategy": "GoldenGate",
    "backbone_name": "pET28a",
    "backbone_sequence": (
        "AGATCTCGATCCCGCGAAATTAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAAC"
        "TTTAAGAAGGAGATATACCATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCATATGGCTAGC"
        "ATGACTGGTGGACAGCAAATGGGTCGCGGATCCGAATTCGAGCTCCGTCGACAAGCTTGCGGCCGCACTCGAGCACCACCACCACCA"
        "CCACTGAGATCCGGCTGCTAACAAAGCCCGAAAGGAAGCTGAGTTGGCTGCTGCCACCGCTGAGCAATAACTAGCATAACCCCTTGGG"
        "GCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGGAGGAACTATATCCGG"
    ),
    "insert_name": "REP24",
    "insert_sequence": (
        "ATGAAAAATGTTTTAATGGTTACTACTTCTCATGATGTTATGGGTAATTCTAATGAAAAAACTGGTTTATGGTTATCTGAATTAACTC"
        "ATCCTTATTATTCTATTATTGATAAAAATATTAATATTGATATTGTTTCTATTATGGGTGGTGAAATTCCTATTGATCCTAATTCTGT"
        "TGCTCAAGAAGATTATTATAATGATAAATTTTTAGCTGATGATAATTTAAAAAATATTATGAAAAATTCTACTTCTTTACGTGATGTT"
        "AATATTAAAGAATATGATGCTATTATTTTTGCTGGTGGTCATGGTACTATGTGGGATTTTCCTAATAATGCTAATATTCATTCTAAAG"
    ),
    "insert_forward_primer_name": "my_insert_fwd",
    "insert_forward_primer_sequence": "CCATAGGTCTCAATGAAAAATGTTTTAATGGTTACTA",
    "insert_reverse_primer_name": "my_insert_rev",
    "insert_reverse_primer_sequence": "CAGATGGTCTCACGAGATTAAATAAAGATTGTTTAAT",
    "vector_forward_primer_name": "vector_left_primer",
    "vector_forward_primer_sequence": "CCATAGGTCTCACTCGAGCACCACCACCACCACCACT",
    "vector_reverse_primer_name": "vector_right_primer",
    "vector_reverse_primer_sequence": "CAGATGGTCTCATCATGGTATATCTCCTTCTTAAAGT",
    "enzyme": "BsaI",
    "cell_strain": "Mach1",
    "selection": "Kan",
    "temperature_c": 37,
    "notes": "valid real-world example",
}


def _build_record():
    return create_construction_file(**VALID_CONSTRUCTION_INPUT)["structured_construction_file"]


def test_validate_passes_valid_record():
    report = validate_construction_record(_build_record(), strict=False)
    assert report["construct_name"] == "pET28a_REP24"
    assert len(report["step_results"]) >= 1


def test_validate_step_types_match_record():
    record = _build_record()
    report = validate_construction_record(record, strict=False)
    step_types = {s["step_type"] for s in report["step_results"]}
    assert "PCR" in step_types
    assert "GoldenGate" in step_types
    assert "Transform" in step_types


def test_validate_fails_with_bad_primer():
    bad_input = dict(VALID_CONSTRUCTION_INPUT)
    bad_input["insert_reverse_primer_sequence"] = "AAAAAAAAAAAAAAAAAAAA"
    record = create_construction_file(**bad_input)["structured_construction_file"]
    report = validate_construction_record(record, strict=False)
    assert report["is_valid"] is False
    assert len(report["errors"]) >= 1


def test_validate_pcr_step_is_first(  ):
    report = validate_construction_record(_build_record(), strict=False)
    assert report["step_results"][0]["step_type"] == "PCR"
