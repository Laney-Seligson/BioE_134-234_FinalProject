import sys
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.construction_file_tools.tools.create_construction_file import create_construction_file

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
        "TTTTAGATATTTATGCTAAAAATGGTGTTATTGGTGCTATTTGTCATGGTGTTGCTGCTTTAATTAATGTTAAAGATAATAATGGTC"
        "AAAATATTATTCGTGATAAAGAAGTTACTGGTTTTTCTAATAATGAAGAAAAAATTGTTGGTTTAACTGATGTTGTTCCTTTTTCTTA"
        "TAGAAGATTCTTTAGTTGAAGCTGGTGCTAAATATTCTTCTGCTTCTGAATGGCAATCTTATGTTAAATCTGATTCTAAAATTATTAC"
        "TGCTCAAAATCCTCAATCTGCTACTGATTTTGCTAAAGCTATTAAACAATCTTTATTTAAT"
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


def test_create_construction_file_returns_required_keys():
    result = create_construction_file(**VALID_CONSTRUCTION_INPUT)
    assert result["construct_name"] == "pET28a_REP24"
    assert result["assembly_strategy"] == "GoldenGate"
    assert "structured_construction_file" in result
    assert "construction_file_txt" in result


def test_create_construction_file_structured_record_has_operations():
    result = create_construction_file(**VALID_CONSTRUCTION_INPUT)
    record = result["structured_construction_file"]
    assert "operations" in record
    assert len(record["operations"]) > 0


def test_create_construction_file_has_pcr_and_transform_steps():
    result = create_construction_file(**VALID_CONSTRUCTION_INPUT)
    step_types = {op["step_type"] for op in result["structured_construction_file"]["operations"]}
    assert "PCR" in step_types
    assert "Transform" in step_types


def test_create_construction_file_goldengate_step_present():
    result = create_construction_file(**VALID_CONSTRUCTION_INPUT)
    step_types = [op["step_type"] for op in result["structured_construction_file"]["operations"]]
    assert "GoldenGate" in step_types
