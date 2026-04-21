"""
Unit tests for the seq_basics example tools.

Each tool is now a class following the Python Function Object Pattern
(initiate / run). Tests cover both the canonical class interface AND the
module-level alias (for example: `reverse_complement = _instance.run`) so that
direct imports continue to work for students who prefer that style.
"""

import pytest
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from modules.seq_basics.tools.translate import translate
from modules.seq_basics.tools.reverse_complement import reverse_complement
from modules.crispr_tools.tools.create_construction_file import create_construction_file
from modules.crispr_tools.tools.design_cloning_oligos import design_cloning_oligos
from modules.crispr_tools.tools.construction_file_validation import (
    validate_construction_record,
)
from modules.crispr_tools.tools.predict_offtargets import predict_offtargets
from modules.crispr_tools.tools.verify_edit import verify_edit


def test_reverse_complement_basic():
    assert reverse_complement("ATGC") == "GCAT"


def test_reverse_complement_ambiguity_codes():
    # Should not error for supported IUPAC subset
    assert reverse_complement("ATRYSWKMN")


def test_translate_basic():
    assert translate("ATGGCT") == "MA"


def test_translate_frame_validation():
    with pytest.raises(ValueError):
        translate("ATGGCT", frame=0)
    with pytest.raises(ValueError):
        translate("ATGGCT", frame=4)


def test_translate_with_coordinates_and_frame():
    # sequence: A ATG GCT AAA
    # start=1 -> ATGGCTAAA
    # frame=1 -> ATG GCT AAA -> M A K
    assert translate("AATGGCTAAA", start=1, end=None, frame=1) == "MAK"


VALID_CONSTRUCTION_INPUT = {
    "construct_name": "pET28a_REP24",
    "assembly_strategy": "GoldenGate",
    "backbone_name": "pET28a",
    "backbone_sequence": (
        "AGATCTCGATCCCGCGAAATTAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACCATGGGCAGCAGCCATCATCATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCATATGGCTAGCATGACTGGTGGACAGCAAATGGGTCGCGGATCCGAATTCGAGCTCCGTCGACAAGCTTGCGGCCGCACTCGAGCACCACCACCACCACCACTGAGATCCGGCTGCTAACAAAGCCCGAAAGGAAGCTGAGTTGGCTGCTGCCACCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGGAGGAACTATATCCGGATTGGCGAATGGGACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTACGCGCAGCGTGACCGCTACACTTGCCAGCGCCCTAGCGCCCGCTCCTTTCGCTTTCTTCCCTTCCTTTCTCGCCACGTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAAAAACTTGATTAGGGTGATGGTTCACGTAGTGGGCCATCGCCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCTTTAATAGTGGACTCTTGTTCCAAACTGGAACAACACTCAACCCTATCTCGGTCTATTCTTTTGATTTATAAGGGATTTTGCCGATTTCGGCCTATTGGTTAAAAAATGAGCTGATTTAACAAAAATTTAACGCGAATTTTAACAAAATATTAACGTTTACAATTTCAGGTGGCACTTTTCGGGGAAATGTGCGCGGAACCCCTATTTGTTTATTTTTCTAAATACATTCAAATATGTATCCGCTCATGAATTAATTCTTAGAAAAACTCATCGAGCATCAAATGAAACTGCAATTTATTCATATCAGGATTATCAATACCATATTTTTGAAAAAGCCGTTTCTGTAATGAAGGAGAAAACTCACCGAGGCAGTTCCATAGGATGGCAAGATCCTGGTATCGGTCTGCGATTCCGACTCGTCCAACATCAATACAACCTATTAATTTCCCCTCGTCAAAAATAAGGTTATCAAGTGAGAAATCACCATGAGTGACGACTGAATCCGGTGAGAATGGCAAAAGTTTATGCATTTCTTTCCAGACTTGTTCAACAGGCCAGCCATTACGCTCGTCATCAAAATCACTCGCATCAACCAAACCGTTATTCATTCGTGATTGCGCCTGAGCGAGACGAAATACGCGATCGCTGTTAAAAGGACAATTACAAACAGGAATCGAATGCAACCGGCGCAGGAACACTGCCAGCGCATCAACAATATTTTCACCTGAATCAGGATATTCTTCTAATACCTGGAATGCTGTTTTCCCGGGGATCGCAGTGGTGAGTAACCATGCATCATCAGGAGTACGGATAAAATGCTTGATGGTCGGAAGAGGCATAAATTCCGTCAGCCAGTTTAGTCTGACCATCTCATCTGTAACATCATTGGCAACGCTACCTTTGCCATGTTTCAGAAACAACTCTGGCGCATCGGGCTTCCCATACAATCGATAGATTGTCGCACCTGATTGCCCGACATTATCGCGAGCCCATTTATACCCATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTAGAGCAAGACGTTTCCCGTTGAATATGGCTCATAACACCCCTTGTATTACTGTTTATGTAAGCAGACAGTTTTATTGTTCATGACCAAAATCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCGTAGAAAAGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATCAAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGTGTAGCCGTAGTTAGGCCACCACTTCAAGAACTCTGTAGCACCGCCTACATACCTCGCTCTGCTAATCCTGTTACCAGTGGCTGCTGCCAGTGGCGATAAGTCGTGTCTTACCGGGTTGGACTCAAGACGATAGTTACCGGATAAGGCGCAGCGGTCGGGCTGAACGGGGGGTTCGTGCACACAGCCCAGCTTGGAGCGAACGACCTACACCGAACTGAGATACCTACAGCGTGAGCTATGAGAAAGCGCCACGCTTCCCGAAGGGAGAAAGGCGGACAGGTATCCGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCCAGGGGGAAACGCCTGGTATCTTTATAGTCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGCCTATGGAAAAACGCCAGCAACGCGGCCTTTTTACGGTTCCTGGCCTTTTGCTGGCCTTTTGCTCACATGTTCTTTCCTGCGTTATCCCCTGATTCTGTGGATAACCGTATTACCGCCTTTGAGTGAGCTGATACCGCTCGCCGCAGCCGAACGACCGAGCGCAGCGAGTCAGTGAGCGAGGAAGCGGAAGAGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAGTTAAGCCAGTATACACTCCGCTATCGCTACGTGACTGGGTCATGGCTGCGCCCCGACACCCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGGCATCCGCTTACAGACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTCACCGTCATCACCGAAACGCGCGAGGCAGCTGCGGTAAAGCTCATCAGCGTGGTCGTGAAGCGATTCACAGATGTCTGCCTGTTCATCCGCGTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCACTGATGCCTCCGTGTAAGGGGGATTTCTGTTCATGGGGGTAATGATACCGATGAAACGAGAGAGGATGCTCACGATACGGGTTACTGATGATGAACATGCCCGGTTACTGGAACGTTGTGAGGGTAAACAACTGGCGGTATGGATGCGGCGGGACCAGAGAAAAATCACTCAGGGTCAATGCCAGCGCTTCGTTAATACAGATGTAGGTGTTCCACAGGGTAGCCAGCAGCATCCTGCGATGCAGATCCGGAACATAATGGTGCAGGGCGCTGACTTCCGCGTTTCCAGACTTTACGAAACACGGAAACCGAAGACCATTCATGTTGTTGCTCAGGTCGCAGACGTTTTGCAGCAGCAGTCGCTTCACGTTCGCTCGCGTATCGGTGATTCATTCTGCTAACCAGTAAGGCAACCCCGCCAGCCTAGCCGGGTCCTCAACGACAGGAGCACGATCATGCGCACCCGTGGGGCCGCCATGCCGGCGATAATGGCCTGCTTCTCGCCGAAACGTTTGGTGGCGGGACCAGTGACGAAGGCTTGAGCGAGGGCGTGCAAGATTCCGAATACCGCAAGCGACAGGCCGATCATCGTCGCGCTCCAGCGAAAGCGGTCCTCGCCGAAAATGACCCAGAGCGCTGCCGGCACCTGTCCTACGAGTTGCATGATAAAGAAGACAGTCATAAGTGCGGCGACGATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTCTCAAGGGCATCGGTCGAGATCCCGGTGCCTAATGAGTGAGCTAACTTACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCCAGGGTGGTTTTTCTTTTCACCAGTGAGACGGGCAACAGCTGATTGCCCTTCACCGCCTGGCCCTGAGAGAGTTGCAGCAAGCGGTCCACGCTGGTTTGCCCCAGCAGGCGAAAATCCTGTTTGATGGTGGTTAACGGCGGGATATAACATGAGCTGTCTTCGGTATCGTCGTATCCCACTACCGAGATATCCGCACCAACGCGCAGCCCGGACTCGGTAATGGCGCGCATTGCGCCCAGCGCCATCTGATCGTTGGCAACCAGCATCGCAGTGGGAACGATGCCCTCATTCAGCATTTGCATGGTTTGTTGAAAACCGGACATGGCACTCCAGTCGCCTTCCCGTTCCGCTATCGGCTGAATTTGATTGCGAGTGAGATATTTATGCCAGCCAGCCAGACGCAGACGCGCCGAGACAGAACTTAATGGGCCCGCTAACAGCGCGATTTGCTGGTGACCCAATGCGACCAGATGCTCCACGCCCAGTCGCGTACCGTCTTCATGGGAGAAAATAATACTGTTGATGGGTGTCTGGTCAGAGACATCAAGAAATAACGCCGGAACATTAGTGCAGGCAGCTTCCACAGCAATGGCATCCTGGTCATCCAGCGGATAGTTAATGATCAGCCCACTGACGCGTTGCGCGAGAAGATTGTGCACCGCCGCTTTACAGGCTTCGACGCCGCTTCGTTCTACCATCGACACCACCACGCTGGCACCCAGTTGATCGGCGCGAGATTTAATCGCCGCGACAATTTGCGACGGCGCGTGCAGGGCCAGACTGGAGGTGGCAACGCCAATCAGCAACGACTGTTTGCCCGCCAGTTGTTGTGCCACGCGGTTGGGAATGTAATTCAGCTCCGCCATCGCCGCTTCCACTTTTTCCCGCGTTTTCGCAGAAACGTGGCTGGCCTGGTTCACCACGCGGGAAACGGTCTGATAAGAGACACCGGCATACTCTGCGACATCGTATAACGTTACTGGTTTCACATTCACCACCCTGAATTGACTCTCTTCCGGGCGCTATCATGCCATACCGCGAAAGGTTTTGCGCCATTCGATGGTGTCCGGGATCTCGACGCTCTCCCTTATGCGACTCCTGCATTAGGAAGCAGCCCAGTAGTAGGTTGAGGCCGTTGAGCACCGCCGCCGCAAGGAATGGTGCATGCAAGGAGATGGCGCCCAACAGTCCCCCGGCCACGGGGCCTGCCACCATACCCACGCCGAAACAAGCGCTCATGAGCCCGAAGTGGCGAGCCCGATCTTCCCCATCGGTGATGTCGGCGATATAGGCGCCAGCAACCGCACCTGTGGCGCCGGTGATGCCGGCCACGATGCGTCCGGCGTAGAGGATCG"
    ),
    "insert_name": "REP24",
    "insert_sequence": (
        "ATGAAAAATGTTTTAATGGTTACTACTTCTCATGATGTTATGGGTAATTCTAATGAAAAAACTGGTTTATGGTTATCTGAATTAACTCATCCTTATTATTCTATTATTGATAAAAATATTAATATTGATATTGTTTCTATTATGGGTGGTGAAATTCCTATTGATCCTAATTCTGTTGCTCAAGAAGATTATTATAATGATAAATTTTTAGCTGATGATAATTTAAAAAATATTATGAAAAATTCTACTTCTTTACGTGATGTTAATATTAAAGAATATGATGCTATTATTTTTGCTGGTGGTCATGGTACTATGTGGGATTTTCCTAATAATGCTAATATTCATTCTAAAGTTTTAGATATTTATGCTAAAAATGGTGTTATTGGTGCTATTTGTCATGGTGTTGCTGCTTTAATTAATGTTAAAGATAATAATGGTCAAAATATTATTCGTGATAAAGAAGTTACTGGTTTTTCTAATAATGAAGAAAAAATTGTTGGTTTAACTGATGTTGTTCCTTTTTCTTTAGAAGATTCTTTAGTTGAAGCTGGTGCTAAATATTCTTCTGCTTCTGAATGGCAATCTTATGTTAAATCTGATTCTAAAATTATTACTGCTCAAAATCCTCAATCTGCTACTGATTTTGCTAAAGCTATTAAACAATCTTTATTTAAT"
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
    "notes": "valid real-world example with non-fixed names",
}


def test_create_construction_file_basic():
    result = create_construction_file(**VALID_CONSTRUCTION_INPUT)

    assert result["construct_name"] == "pET28a_REP24"
    assert result["assembly_strategy"] == "GoldenGate"
    assert "structured_construction_file" in result
    assert "construction_file_txt" in result


def test_design_cloning_oligos_defaults_to_ecoli_pcrispr_reference():
    result = design_cloning_oligos(
        protospacer="TACTTTACGCAGCGCGGAGT",
        organism="E. coli",
        target_reference="ecoli_rpsl",
    )

    assert result["vector"] == "pcrispr"
    assert result["enzyme"] == "BsaI"
    assert result["source"].startswith("Jiang et al. Nat Biotechnol 2013")
    assert result["target_verification"]["verified"] is True
    assert result["target_verification"]["reference"] == "ecoli_rpsl"
    assert result["top_oligo"] == "AAACTACTTTACGCAGCGCGGAGT"
    assert result["bottom_oligo"] == "AAAACACTCCGCGCTGCGTAAAGTA"
    assert result["construction_file_inputs"]["cell_strain"] == "HME63 or MG1655 carrying pCas9"
    assert result["construction_file_inputs"]["selection"] == "Kan"


def test_design_cloning_oligos_rejects_unsupported_organism_verification():
    with pytest.raises(ValueError, match="downloaded local references"):
        design_cloning_oligos(
            protospacer="CCGGATGCTCCTCAGCTCTG",
            vector="pET28a",
            organism="zebrafish",
        )


def test_validate_construction_record_name_agnostic():
    result = create_construction_file(**VALID_CONSTRUCTION_INPUT)
    record = result["structured_construction_file"]

    report = validate_construction_record(record, strict=False)

    assert report["construct_name"] == "pET28a_REP24"
    assert len(report["step_results"]) == 4

    first_step = report["step_results"][0]
    assert first_step["step_type"] == "PCR"
    assert first_step["output_name"] == record["operations"][0]["output"]

    second_step = report["step_results"][1]
    assert second_step["step_type"] == "PCR"
    assert second_step["output_name"] == record["operations"][1]["output"]

    third_step = report["step_results"][2]
    assert third_step["step_type"] == "GoldenGate"

    fourth_step = report["step_results"][3]
    assert fourth_step["step_type"] == "Transform"


def test_validate_construction_record_fails_with_bad_primer():
    bad_input = dict(VALID_CONSTRUCTION_INPUT)
    bad_input["insert_reverse_primer_sequence"] = "AAAAAAAAAAAAAAAAAAAA"

    result = create_construction_file(**bad_input)
    record = result["structured_construction_file"]

    report = validate_construction_record(record, strict=False)

    assert report["is_valid"] is False
    assert report["step_results"][0]["step_type"] == "PCR"
    assert report["step_results"][0]["is_valid"] is False
    assert len(report["errors"]) >= 1


# ---------------------------------------------------------------------------
# predict_offtargets tests
# ---------------------------------------------------------------------------

def test_predict_offtargets_exact_match_high_risk():
    # exact match + NGG PAM → should be flagged HIGH risk
    result = predict_offtargets(
        protospacer="ATGATGATGATGATGATGAT",
        reference="ATGATGATGATGATGATGATAGG",
    )
    assert len(result["offtarget_sites"]) >= 1
    top = result["offtarget_sites"][0]
    assert top["mismatches"] == 0
    assert top["has_pam"] is True
    assert top["risk"] == "HIGH"


def test_predict_offtargets_no_match_returns_empty():
    # no similarity → no sites within default 3-mismatch threshold
    result = predict_offtargets(
        protospacer="ATGATGATGATGATGATGAT",
        reference="CCCCCCCCCCCCCCCCCCCCAGG",
    )
    # all sites should have mismatches > 0 (or list is empty)
    for site in result["offtarget_sites"]:
        assert site["mismatches"] > 0


def test_predict_offtargets_empty_protospacer_raises():
    with pytest.raises(ValueError):
        predict_offtargets(protospacer="", reference="ATGATGATG")


def test_predict_offtargets_empty_reference_raises():
    with pytest.raises(ValueError):
        predict_offtargets(protospacer="ATGATGATGATGATGATGAT", reference="")


def test_predict_offtargets_max_mismatches_zero():
    # with max_mismatches=0, only the perfect-match site should come back
    result = predict_offtargets(
        protospacer="ATGATGATGATGATGATGAT",
        reference="ATGATGATGATGATGATGATAGG",
        max_mismatches=0,
    )
    for site in result["offtarget_sites"]:
        assert site["mismatches"] == 0


# ---------------------------------------------------------------------------
# verify_edit tests
# ---------------------------------------------------------------------------

_VERIFY_REFERENCE = (
    "CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTT"
    "AGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
)
_VERIFY_PROTOSPACER = "TCAGAAACCTGCCAGTTTGC"


def test_verify_edit_standard_case():
    # protospacer at position 15, PAM = TGG, cut at position 32
    result = verify_edit(protospacer=_VERIFY_PROTOSPACER, reference=_VERIFY_REFERENCE)
    assert isinstance(result["cut_position"], int)
    assert result["cut_position"] >= 0
    assert result["pam_sequence"] == "TGG"
    assert result["strand"] == "+"
    assert result["protospacer_position"] == 15
    assert result["cut_position"] == 32


def test_verify_edit_protospacer_not_in_reference_raises():
    with pytest.raises(ValueError):
        verify_edit(
            protospacer="AAAAAAAAAAAAAAAAAAAA",
            reference="ATGCATGCATGC",
        )


def test_verify_edit_empty_protospacer_raises():
    with pytest.raises(ValueError):
        verify_edit(protospacer="", reference="ATGCATGCATGC")


def test_verify_edit_empty_reference_raises():
    with pytest.raises(ValueError):
        verify_edit(protospacer=_VERIFY_PROTOSPACER, reference="")


def test_verify_edit_reverse_strand():
    # RC of _VERIFY_PROTOSPACER is GCAAACTGGCAGGTTTCTGA
    # Build a reference that contains the RC (so the tool must search the minus strand)
    rc_protospacer = "GCAAACTGGCAGGTTTCTGA"
    # CCN before RC protospacer = PAM on minus strand (NGG when RC'd)
    # prepend CCT (RC = AGG, an NGG) so the PAM is valid
    reference = "AAAAAAAAAAAAAAAAAAAAACCT" + rc_protospacer + "AAAAAAAAAAAAAAAAAAAAAA"
    result = verify_edit(protospacer=_VERIFY_PROTOSPACER, reference=reference)
    assert result["strand"] == "-"
    assert result["pam_sequence"][1] == "G"
    assert result["pam_sequence"][2] == "G"
    assert isinstance(result["cut_position"], int)


# ---------------------------------------------------------------------------
# Additional edge case tests
# ---------------------------------------------------------------------------

def test_predict_offtargets_no_pam_at_end_of_reference():
    # protospacer fits at position 0 but there is no room for a PAM after it
    protospacer = "ATGATGATGATGATGATGAT"
    reference = protospacer  # exactly 20 bp, no trailing bases
    result = predict_offtargets(protospacer=protospacer, reference=reference, max_mismatches=0)
    # the site should be found but has_pam must be False (no bases after the window)
    assert len(result["offtarget_sites"]) >= 1
    assert result["offtarget_sites"][0]["has_pam"] is False


def test_predict_offtargets_circular_wraps_origin():
    # Arrange so the protospacer straddles the circular origin:
    #   reference = [protospacer[10:]] [AGG PAM] [filler] [protospacer[:10]]
    # When the reference wraps, position 33 onward reads protospacer[:10] + protospacer[10:]
    # = the full protospacer, with AGG at position 10 as the PAM.
    protospacer = "ATGATGATGATGATGATGAT"
    reference = protospacer[10:] + "AGG" + "CCCCCCCCCCCCCCCCCCCC" + protospacer[:10]
    result_circular = predict_offtargets(
        protospacer=protospacer, reference=reference, max_mismatches=0, is_circular=True
    )
    result_linear = predict_offtargets(
        protospacer=protospacer, reference=reference, max_mismatches=0, is_circular=False
    )
    # circular scan should find the wrapped site; linear scan should miss it
    assert len(result_circular["offtarget_sites"]) > len(result_linear["offtarget_sites"])
