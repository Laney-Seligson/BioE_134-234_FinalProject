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
from modules.crispr_tools.tools.rank_guides import rank_guides
from modules.crispr_tools.tools.colony_calculator import colony_calculator
from modules.crispr_tools.tools.interpret_ice_tide import interpret_ice_tide
from modules.crispr_tools.tools.predict_editing_efficiency import predict_editing_efficiency
from modules.crispr_tools.tools.verify_edit import verify_edit
from modules.crispr_tools.tools.lab_sheet import lab_sheet
from modules.crispr_tools.tools.lookup_gene_sequence import LookupGeneSequence


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
    assert result["construction_file_inputs"]["backbone_name"] == "pCRISPR_rpsL"
    assert result["construction_file_inputs"]["backbone_sequence"] != "N"
    assert len(result["construction_file_inputs"]["backbone_sequence"]) == 2433
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

# ---------------------------------------------------------------------------
# Cas12a tests for predict_offtargets
# ---------------------------------------------------------------------------

def test_predict_offtargets_cas12a_exact_match_high_risk():
    # TTTA PAM (TTTV where V=A) before 23bp protospacer → HIGH risk
    protospacer = "ATGATGATGATGATGATGATGAT"
    reference = "TTTA" + protospacer + "AAAAAAAAAA"
    result = predict_offtargets(protospacer=protospacer, reference=reference, nuclease="cas12a")
    assert len(result["offtarget_sites"]) >= 1
    top = result["offtarget_sites"][0]
    assert top["mismatches"] == 0
    assert top["has_pam"] is True
    assert top["risk"] == "HIGH"


def test_predict_offtargets_cas12a_no_pam_returns_low_risk():
    # 23bp protospacer present but no TTTV before it
    protospacer = "ATGATGATGATGATGATGATGAT"
    reference = "AAAA" + protospacer + "AAAAAAAAAA"
    result = predict_offtargets(
        protospacer=protospacer, reference=reference, nuclease="cas12a", max_mismatches=0
    )
    for site in result["offtarget_sites"]:
        assert site["has_pam"] is False
        assert site["risk"] == "LOW"


def test_predict_offtargets_cas12a_invalid_nuclease_raises():
    with pytest.raises(ValueError):
        predict_offtargets(
            protospacer="ATGATGATGATGATGATGAT",
            reference="ATGATGATG",
            nuclease="talenuclease",
        )


def test_predict_offtargets_cas12a_result_includes_nuclease_key():
    protospacer = "ATGATGATGATGATGATGATGAT"
    reference = "TTTA" + protospacer + "AAAAAAAAAA"
    result = predict_offtargets(protospacer=protospacer, reference=reference, nuclease="cas12a")
    assert result["nuclease"] == "cas12a"


# ---------------------------------------------------------------------------
# Cas12a tests for verify_edit
# ---------------------------------------------------------------------------

def test_verify_edit_cas12a_forward_strand():
    # TTTA (TTTV) PAM before 23bp protospacer; cut at ps_position + 18
    protospacer = "ATGATGATGATGATGATGATGAT"
    # ps_position = 4, expected cut = 4 + 18 = 22
    reference = "TTTA" + protospacer + "A" * 200
    result = verify_edit(protospacer=protospacer, reference=reference, nuclease="cas12a")
    assert result["strand"] == "+"
    assert result["pam_sequence"] == "TTTA"
    assert result["cut_position"] == 22
    assert result["nuclease"] == "cas12a"


def test_verify_edit_cas12a_reverse_strand():
    # RC protospacer on + strand, TAAA after it (= RC of TTTA = TTTV PAM on minus strand)
    protospacer = "ATGATGATGATGATGATGATGAT"
    from modules.crispr_tools.tools.verify_edit import _reverse_complement
    rc = _reverse_complement(protospacer)
    # forward strand: [padding][RC protospacer][TAAA][padding]
    reference = "A" * 20 + rc + "TAAA" + "A" * 200
    result = verify_edit(protospacer=protospacer, reference=reference, nuclease="cas12a")
    assert result["strand"] == "-"
    assert result["pam_sequence"] == "TTTA"


def test_verify_edit_cas12a_wrong_pam_raises():
    # AAAA before protospacer — not a valid TTTV PAM
    protospacer = "ATGATGATGATGATGATGATGAT"
    reference = "AAAA" + protospacer + "A" * 200
    with pytest.raises(ValueError, match="TTTV"):
        verify_edit(protospacer=protospacer, reference=reference, nuclease="cas12a")


def test_verify_edit_invalid_nuclease_raises():
    with pytest.raises(ValueError):
        verify_edit(
            protospacer="TCAGAAACCTGCCAGTTTGC",
            reference=_VERIFY_REFERENCE,
            nuclease="talenuclease",
        )


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


# ---------------------------------------------------------------------------
# lab_sheet tests
# ---------------------------------------------------------------------------

def _make_record(strategy="GoldenGate", selection="Kan"):
    ops = [
        {
            "step_number": 1,
            "step_type": "PCR",
            "inputs": ["ins_fwd", "ins_rev", "insert"],
            "parameters": {"forward_primer": "ins_fwd", "reverse_primer": "ins_rev", "template": "insert"},
            "output": "ins_pcr",
            "description": "",
        },
        {
            "step_number": 2,
            "step_type": "PCR",
            "inputs": ["vec_fwd", "vec_rev", "backbone"],
            "parameters": {"forward_primer": "vec_fwd", "reverse_primer": "vec_rev", "template": "backbone"},
            "output": "vec_pcr",
            "description": "",
        },
    ]
    if strategy == "GoldenGate":
        ops.append({
            "step_number": 3,
            "step_type": "GoldenGate",
            "inputs": ["vec_pcr", "ins_pcr"],
            "parameters": {"enzyme": "BsaI"},
            "output": "my_construct",
            "description": "",
        })
    else:
        ops.append({
            "step_number": 3,
            "step_type": "Gibson",
            "inputs": ["vec_pcr", "ins_pcr"],
            "parameters": {"reagent": "GibsonMix", "overlap_bp": 20},
            "output": "my_construct",
            "description": "",
        })
    ops.append({
        "step_number": 4,
        "step_type": "Transform",
        "inputs": ["my_construct"],
        "parameters": {"cells": "Mach1", "selection": selection, "temperature_c": 37},
        "output": "my_construct_e",
        "description": "",
    })
    return {
        "construct_name": "my_construct",
        "assembly_strategy": strategy,
        "parts": [],
        "operations": ops,
        "notes": "",
    }


def test_lab_sheet_empty_record_raises():
    with pytest.raises(ValueError):
        lab_sheet({})


def test_lab_sheet_returns_required_keys():
    result = lab_sheet(_make_record())
    assert "construct_name" in result
    assert "assembly_strategy" in result
    assert "lab_sheet_text" in result
    assert "step_count" in result


def test_lab_sheet_pcr_section_format():
    result = lab_sheet(_make_record())
    text = result["lab_sheet_text"]
    assert "A: PCR" in text
    assert "samples:" in text
    assert "source:" in text
    assert "primer1" in text
    assert "destination: thermocycler1A" in text
    assert "program: Q5/Q5-4K" in text


def test_lab_sheet_gel_dpni_section():
    result = lab_sheet(_make_record())
    text = result["lab_sheet_text"]
    assert "A: Gel and DpnI" in text
    assert "DpnI" in text
    assert "BstEII ladder" in text


def test_lab_sheet_zymo_section():
    result = lab_sheet(_make_record())
    text = result["lab_sheet_text"]
    assert "A: Zymo" in text
    assert "elution_volume" in text
    assert "50 uL" in text


def test_lab_sheet_goldengate_assemble_section():
    result = lab_sheet(_make_record(strategy="GoldenGate"))
    text = result["lab_sheet_text"]
    assert "A: Assemble" in text
    assert "DNA Mix:" in text
    assert "BsaI" in text
    assert "T4 DNA ligase" in text
    assert "program: main/GG1" in text


def test_lab_sheet_gibson_assemble_section():
    result = lab_sheet(_make_record(strategy="Gibson"))
    text = result["lab_sheet_text"]
    assert "A: Assemble" in text
    assert "2X Gibson Mix" in text
    assert "program: main/GIB2" in text


def test_lab_sheet_transform_no_rescue_for_amp():
    result = lab_sheet(_make_record(selection="Amp"))
    text = result["lab_sheet_text"]
    assert "rescue_required: no" in text
    assert "rescue step" not in text


def test_lab_sheet_transform_rescue_for_spec():
    result = lab_sheet(_make_record(selection="Spec"))
    text = result["lab_sheet_text"]
    assert "rescue_required: yes" in text
    assert "rescue step" in text


def test_lab_sheet_thread_prefix():
    result = lab_sheet(_make_record(), thread="K")
    text = result["lab_sheet_text"]
    assert "K: PCR" in text
    assert "K: Assemble" in text
    assert "K: Transform" in text
    assert "K1" in text


def test_lab_sheet_step_count_includes_all_sections():
    result = lab_sheet(_make_record())
    # PCR + Gel+DpnI + Zymo + Assemble + Transform + Pick + Miniprep + Sequencing = 8
    assert result["step_count"] == 8


def test_lab_sheet_pick_section():
    result = lab_sheet(_make_record())
    text = result["lab_sheet_text"]
    assert "A: Pick" in text
    assert "number" in text
    assert "labels" in text
    assert "A1A" in text


def test_lab_sheet_miniprep_section():
    result = lab_sheet(_make_record())
    text = result["lab_sheet_text"]
    assert "A: Miniprep" in text
    assert "culture" in text
    assert "location" in text


def test_lab_sheet_sequencing_section():
    result = lab_sheet(_make_record())
    text = result["lab_sheet_text"]
    assert "A: Sequencing" in text
    assert "L4440" in text
    assert "2.66 uM" in text
    assert "Stanley Hall" in text


def test_lab_sheet_direct_synthesis():
    record = {
        "construct_name": "synth_insert",
        "assembly_strategy": "DirectSynthesis",
        "parts": [],
        "operations": [
            {
                "step_number": 1,
                "step_type": "DirectSynthesis",
                "inputs": ["synth_insert"],
                "parameters": {},
                "output": "synth_insert",
                "description": "",
            }
        ],
        "notes": "",
    }
    result = lab_sheet(record)
    text = result["lab_sheet_text"]
    assert "DirectSynthesis" in text
    assert "vendor" in text


def test_lab_sheet_notes_appended_when_present():
    record = _make_record()
    record["notes"] = "Handle with care."
    result = lab_sheet(record, include_notes=True)
    assert "Handle with care." in result["lab_sheet_text"]


def test_lab_sheet_notes_suppressed_when_flag_false():
    record = _make_record()
    record["notes"] = "Handle with care."
    result = lab_sheet(record, include_notes=False)
    assert "Handle with care." not in result["lab_sheet_text"]


# ---------------------------------------------------------------------------
# lookup_gene_sequence tests (NCBI calls are mocked)
# ---------------------------------------------------------------------------

def _make_lookup_instance():
    inst = LookupGeneSequence()
    inst.initiate(email="test@example.com")
    return inst


def test_lookup_gene_sequence_empty_gene_name_raises():
    inst = _make_lookup_instance()
    with pytest.raises(ValueError, match="gene_name"):
        inst.run(gene_name="", organism="E. coli")


def test_lookup_gene_sequence_empty_organism_raises():
    inst = _make_lookup_instance()
    with pytest.raises(ValueError, match="organism"):
        inst.run(gene_name="lacZ", organism="")


def test_lookup_gene_sequence_alias_resolved(monkeypatch):
    """Common alias 'E. coli' maps to 'Escherichia coli' before hitting NCBI."""
    captured = {}

    def fake_esearch(db, term, retmax):
        captured["term"] = term
        # simulate NCBI returning one hit
        class FakeHandle:
            def close(self): pass
        return FakeHandle()

    def fake_read(handle):
        if "IdList" not in captured:
            return {"IdList": ["945006"]}
        return {"IdList": []}

    monkeypatch.setattr("Bio.Entrez.esearch", fake_esearch)
    monkeypatch.setattr("Bio.Entrez.read", fake_read)

    inst = _make_lookup_instance()
    # will fail at _fetch_cds_for_gene since we only mocked esearch, but
    # we can check the query term was resolved correctly before the error
    try:
        inst.run(gene_name="lacZ", organism="E. coli")
    except (ValueError, Exception):
        pass

    assert "Escherichia coli" in captured.get("term", "")


def test_lookup_gene_sequence_no_results_raises(monkeypatch):
    """When NCBI returns no gene IDs, ValueError is raised."""
    class FakeHandle:
        def close(self): pass

    monkeypatch.setattr("Bio.Entrez.esearch", lambda **kw: FakeHandle())
    monkeypatch.setattr("Bio.Entrez.read", lambda h: {"IdList": []})

    inst = _make_lookup_instance()
    with pytest.raises(ValueError, match="No gene found"):
        inst.run(gene_name="xyzzy_fake_9999", organism="E. coli")


def test_lookup_gene_sequence_cds_extraction(monkeypatch):
    """Full happy-path mock: returns a CDS sequence from a GenBank record."""
    import io

    fake_gb = (
        "LOCUS       NM_000001               10 bp    mRNA    linear   BCT 01-JAN-2020\n"
        "DEFINITION  Escherichia coli lacZ mRNA.\n"
        "ACCESSION   NM_000001\n"
        "VERSION     NM_000001.1\n"
        "FEATURES             Location/Qualifiers\n"
        "     CDS             1..10\n"
        '                     /product="beta-galactosidase"\n'
        "ORIGIN\n"
        "        1 atgaaatttg\n"
        "//\n"
    )

    class FakeHandle:
        def close(self): pass

    read_responses = iter([
        {"IdList": ["945006"]},
        [{"LinkSetDb": [{"Link": [{"Id": "12345"}]}]}],
    ])

    monkeypatch.setattr("Bio.Entrez.esearch", lambda **kw: FakeHandle())
    monkeypatch.setattr("Bio.Entrez.read", lambda h: next(read_responses))
    monkeypatch.setattr("Bio.Entrez.elink", lambda **kw: FakeHandle())
    monkeypatch.setattr(
        "Bio.Entrez.efetch",
        lambda **kw: io.StringIO(fake_gb) if kw.get("db") == "nucleotide" else FakeHandle(),
    )

    inst = _make_lookup_instance()
    result = inst.run(gene_name="lacZ", organism="E. coli")

    assert result["gene_name"] == "lacZ"
    assert result["organism"] == "Escherichia coli"
    assert result["product"] == "beta-galactosidase"
    assert result["source"] == "CDS annotation"
    assert len(result["sequence"]) > 0


# ---------------------------------------------------------------------------
# verify_edit primer Tm validation tests
# ---------------------------------------------------------------------------

def test_verify_edit_returns_primer_tm_fields():
    protospacer = "TCAGAAACCTGCCAGTTTGC"
    reference = "CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC" + "A" * 200
    result = verify_edit(protospacer=protospacer, reference=reference, nuclease="cas9")
    assert "forward_primer_tm" in result
    assert "reverse_primer_tm" in result
    assert "tm_difference" in result
    assert "primer_warnings" in result
    assert isinstance(result["primer_warnings"], list)
    # Tm should be a reasonable Wallace-rule value for 20bp primer (somewhere 40-80)
    assert 30 <= result["forward_primer_tm"] <= 90
    assert 30 <= result["reverse_primer_tm"] <= 90
    # tm_difference should be non-negative
    assert result["tm_difference"] >= 0


def test_verify_edit_flags_high_gc_primer_warning():
    # Build a reference where primers fall in a GC-rich zone -> should warn
    # Put the cut site in the middle of a GC-rich reference
    protospacer = "ATGATGATGATGATGATGAT"
    reference = "G" * 160 + protospacer + "AGG" + "G" * 160
    result = verify_edit(protospacer=protospacer, reference=reference, nuclease="cas9")
    # Warnings should contain something about GC or Tm being too high
    assert len(result["primer_warnings"]) > 0


def test_verify_edit_balanced_primers_have_no_warnings():
    # Build a balanced reference with ~50% GC around the cut site
    flank = "ATGCATGCATGCATGCATGC" * 10  # 200bp, 50% GC, no poly-N
    protospacer = "ATGATGATGATGATGATGAT"
    reference = flank + protospacer + "AGG" + flank
    result = verify_edit(protospacer=protospacer, reference=reference, nuclease="cas9")
    # Tms should be ~60 (Wallace: 10 AT * 2 + 10 GC * 4 = 60)
    assert 55 <= result["forward_primer_tm"] <= 65
    assert 55 <= result["reverse_primer_tm"] <= 65
    assert result["tm_difference"] <= 5


# ---------------------------------------------------------------------------
# rank_guides tests
# ---------------------------------------------------------------------------

def test_rank_guides_scores_single_guide_cas9():
    # GC = 50%, no TTTT, ends in G -> efficiency max (3)
    protospacer = "ATGCATGCATGCATGCATGG"
    reference = protospacer + "AGG" + "C" * 60
    result = rank_guides(
        guides=[{"protospacer": protospacer}],
        reference=reference,
        nuclease="cas9",
    )
    top = result["ranked_guides"][0]
    assert top["efficiency_score"] == 3
    assert top["efficiency_details"]["gc_content_ok"] is True
    assert top["efficiency_details"]["no_polyt_run"] is True
    assert top["efficiency_details"]["g_at_pam_proximal"] is True
    assert result["best_guide"] is top
    assert "scoring_rationale" in result


def test_rank_guides_penalizes_polyt_and_low_gc():
    # all A's: GC = 0% (fail), no TTTT but no G at PAM-proximal
    bad = "AAAAAAAAAAAAAAAAAAAA"
    # balanced good guide
    good = "ATGCATGCATGCATGCATGG"
    reference = good + "AGG" + "C" * 40 + bad + "AGG" + "C" * 40
    result = rank_guides(
        guides=[{"protospacer": bad}, {"protospacer": good}],
        reference=reference,
        nuclease="cas9",
    )
    # good should rank first
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
    # Cas12a has no PAM-proximal G bonus -> max efficiency = 2
    protospacer = "ATGATGATGATGATGATGATGAT"  # 23 bp, GC ~ 33% actually -> fails GC
    # Use a balanced 23-mer: GC 40-70%, no TTTT
    protospacer = "ATGCATGCATGCATGCATGCATG"  # GC = 52%
    reference = "TTTA" + protospacer + "A" * 40
    result = rank_guides(
        guides=[{"protospacer": protospacer}],
        reference=reference,
        nuclease="cas12a",
    )
    top = result["ranked_guides"][0]
    assert top["efficiency_score"] <= 2
    assert "g_at_pam_proximal" not in top["efficiency_details"]


# ---------------------------------------------------------------------------
# colony_calculator tests
# ---------------------------------------------------------------------------

def test_colony_calculator_50pct_efficiency_one_clone():
    # P(X>=1) = 1 - 0.5^n >= 0.95 -> n=5 (1-0.5^5=0.96875)
    result = colony_calculator(editing_efficiency=0.5, desired_clones=1, confidence=0.95)
    assert result["colonies_to_pick"] == 5
    assert result["probability_at_chosen_n"] >= 0.95


def test_colony_calculator_low_efficiency_high_burden():
    # 5% efficiency, 1 clone, 95% conf -> 59 colonies
    # log(0.05)/log(0.95) = -2.9957/-0.0513 = 58.4 -> 59
    result = colony_calculator(editing_efficiency=0.05, desired_clones=1, confidence=0.95)
    assert result["colonies_to_pick"] == 59


def test_colony_calculator_zero_efficiency_raises():
    with pytest.raises(ValueError):
        colony_calculator(editing_efficiency=0, desired_clones=1)


def test_colony_calculator_efficiency_above_one_raises():
    with pytest.raises(ValueError):
        colony_calculator(editing_efficiency=1.5)


def test_colony_calculator_preset_hdr_mammalian():
    result = colony_calculator(preset="hdr_mammalian", desired_clones=3, confidence=0.95)
    # 3 HDR clones at 5% efficiency requires many colonies
    assert result["colonies_to_pick"] > 60
    assert result["editing_efficiency"] == 0.05


def test_colony_calculator_unknown_preset_raises():
    with pytest.raises(ValueError):
        colony_calculator(preset="not_a_real_preset")


def test_colony_calculator_no_inputs_raises():
    with pytest.raises(ValueError):
        colony_calculator()


def test_colony_calculator_safety_margin_is_1_5x():
    result = colony_calculator(editing_efficiency=0.5, desired_clones=1, confidence=0.95)
    # 5 * 1.5 = 7.5 -> ceil -> 8
    assert result["safety_margin_recommendation"] == 8


# ---------------------------------------------------------------------------
# interpret_ice_tide tests
# ---------------------------------------------------------------------------

def test_interpret_ice_tide_excellent_result():
    result = interpret_ice_tide(editing_pct=85, r_squared=0.95, tool="ice")
    assert result["efficiency_classification"] == "EXCELLENT"
    assert result["fit_quality"] == "HIGH"
    assert result["is_reliable"] is True
    assert result["warnings"] == []
    assert len(result["next_steps"]) >= 1


def test_interpret_ice_tide_failed_editing():
    result = interpret_ice_tide(editing_pct=5, r_squared=0.95)
    assert result["efficiency_classification"] == "FAILED"
    assert result["is_reliable"] is True


def test_interpret_ice_tide_unreliable_fit():
    result = interpret_ice_tide(editing_pct=85, r_squared=0.5, tool="ice")
    assert result["is_reliable"] is False
    assert any("R^2" in w for w in result["warnings"])
    # next_steps should suggest re-sequencing
    assert any("re-sequence" in step.lower() or "re-pcr" in step.lower() for step in result["next_steps"])


def test_interpret_ice_tide_marginal_efficiency():
    result = interpret_ice_tide(editing_pct=20, r_squared=0.92)
    assert result["efficiency_classification"] == "MARGINAL"
    assert result["is_reliable"] is True


def test_interpret_ice_tide_invalid_editing_pct_raises():
    with pytest.raises(ValueError):
        interpret_ice_tide(editing_pct=-5, r_squared=0.9)
    with pytest.raises(ValueError):
        interpret_ice_tide(editing_pct=150, r_squared=0.9)


def test_interpret_ice_tide_invalid_r_squared_raises():
    with pytest.raises(ValueError):
        interpret_ice_tide(editing_pct=50, r_squared=1.5)


def test_interpret_ice_tide_invalid_tool_raises():
    with pytest.raises(ValueError):
        interpret_ice_tide(editing_pct=50, r_squared=0.9, tool="bogus")


def test_interpret_ice_tide_indel_distribution_dominant():
    result = interpret_ice_tide(
        editing_pct=70,
        r_squared=0.95,
        indel_distribution={"+1": 45.0, "-3": 15.0, "0": 30.0, "+2": 10.0},
    )
    assert result["dominant_indel"] == "+1"
    assert result["dominant_indel_pct"] == 45.0


def test_interpret_ice_tide_unedited_dominant_warns():
    result = interpret_ice_tide(
        editing_pct=30,
        r_squared=0.95,
        indel_distribution={"+1": 15.0, "0": 60.0, "-3": 10.0, "+2": 15.0},
    )
    assert result["dominant_indel"] == "0"
    assert any("unedited" in w.lower() for w in result["warnings"])


# ---------------------------------------------------------------------------
# predict_editing_efficiency tests
# ---------------------------------------------------------------------------

def test_predict_editing_efficiency_good_guide_high_score():
    result = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG",  # GC=50%, ends in G
        pam="AGG",
        nuclease="cas9",
        delivery="rnp",
        outcome="nhej",
    )
    # Good guide + RNP + NHEJ should predict moderate-to-high efficiency
    assert result["on_target_efficiency_pct"] > 50
    assert result["warnings"] == []
    assert "feature_contributions" in result
    assert len(result["confidence_range"]) == 2


def test_predict_editing_efficiency_polyt_guide_dies():
    result = predict_editing_efficiency(
        protospacer="ATGCATTTTTGCATGCATGG",  # contains TTTT
        pam="AGG",
        nuclease="cas9",
    )
    # Poly-T should kill the score and produce a warning
    assert result["on_target_efficiency_pct"] < 25
    assert any("TTTT" in w or "Pol III" in w for w in result["warnings"])


def test_predict_editing_efficiency_weak_pam_warns():
    result = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG",
        pam="AAG",  # NAG is weak
        nuclease="cas9",
    )
    assert any("PAM" in w for w in result["warnings"])


def test_predict_editing_efficiency_hdr_lower_than_nhej():
    nhej = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG", pam="AGG", outcome="nhej"
    )
    hdr = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG", pam="AGG", outcome="hdr"
    )
    # HDR is ~10% of NHEJ
    assert hdr["on_target_efficiency_pct"] < nhej["on_target_efficiency_pct"] / 5


def test_predict_editing_efficiency_rnp_better_than_plasmid():
    rnp = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG", pam="AGG", delivery="rnp"
    )
    plasmid = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG", pam="AGG", delivery="plasmid"
    )
    assert rnp["on_target_efficiency_pct"] > plasmid["on_target_efficiency_pct"]


def test_predict_editing_efficiency_empty_protospacer_raises():
    with pytest.raises(ValueError):
        predict_editing_efficiency(protospacer="", pam="AGG")


def test_predict_editing_efficiency_wrong_length_raises():
    with pytest.raises(ValueError):
        # 19bp instead of 20
        predict_editing_efficiency(protospacer="ATGCATGCATGCATGCATG", pam="AGG", nuclease="cas9")


def test_predict_editing_efficiency_invalid_delivery_raises():
    with pytest.raises(ValueError):
        predict_editing_efficiency(
            protospacer="ATGCATGCATGCATGCATGG", pam="AGG", delivery="magic_beans"
        )


# ---------------------------------------------------------------------------
# CFD score tests for predict_offtargets
# ---------------------------------------------------------------------------

def test_predict_offtargets_cfd_perfect_match_is_one():
    protospacer = "ATGATGATGATGATGATGAT"
    reference = protospacer + "AGG" + "C" * 50
    result = predict_offtargets(protospacer=protospacer, reference=reference, nuclease="cas9")
    # First site should be the on-target with CFD = 1.0
    on = result["offtarget_sites"][0]
    assert on["mismatches"] == 0
    assert on["cfd_score"] == 1.0


def test_predict_offtargets_aggregate_cfd_excludes_on_target():
    protospacer = "ATGATGATGATGATGATGAT"
    # only on-target in reference
    reference = protospacer + "AGG" + "C" * 50
    result = predict_offtargets(protospacer=protospacer, reference=reference, nuclease="cas9")
    # aggregate CFD excludes the on-target -> 0 if no other matches
    assert result["aggregate_offtarget_cfd"] == 0.0


def test_predict_offtargets_pam_distal_mismatch_higher_cfd_than_seed_mismatch():
    # seed mismatch (PAM-proximal) penalizes harder than PAM-distal mismatch
    # protospacer ATGATGATGATGATGATGAT
    # PAM-distal mismatch (position 20 from PAM-proximal end = position 1 in seq):
    #   change first base
    seed_mismatch_off = "ATGATGATGATGATGATGAA"  # last base differs (PAM-proximal pos 1)
    distal_mismatch_off = "TTGATGATGATGATGATGAT"  # first base differs (PAM-distal pos 20)
    # Both with PAM
    reference = seed_mismatch_off + "AGG" + "C" * 30 + distal_mismatch_off + "AGG" + "C" * 30
    result = predict_offtargets(
        protospacer="ATGATGATGATGATGATGAT",
        reference=reference,
        nuclease="cas9",
        max_mismatches=1,
    )
    # find the two off-target entries
    seed_site = next(s for s in result["offtarget_sites"] if s["sequence"] == seed_mismatch_off)
    distal_site = next(s for s in result["offtarget_sites"] if s["sequence"] == distal_mismatch_off)
    # Distal mismatch should have HIGHER CFD (less harmful) than seed mismatch
    assert distal_site["cfd_score"] > seed_site["cfd_score"]


def test_predict_offtargets_no_pam_cfd_is_zero():
    protospacer = "ATGATGATGATGATGATGAT"
    # no NGG after the protospacer
    reference = protospacer + "ACC" + "C" * 50
    result = predict_offtargets(protospacer=protospacer, reference=reference, nuclease="cas9")
    for s in result["offtarget_sites"]:
        if not s["has_pam"]:
            assert s["cfd_score"] == 0.0


# ---------------------------------------------------------------------------
# Citation-field smoke tests — every Jillian tool must emit a "citations"
# field shaped like Laney's: list of {"label", "reference", "claim"} dicts.
# ---------------------------------------------------------------------------

def _assert_well_formed_citations(citations):
    assert isinstance(citations, list)
    assert len(citations) > 0
    for c in citations:
        assert isinstance(c, dict)
        assert set(c.keys()) == {"label", "reference", "claim"}
        assert all(isinstance(c[k], str) and c[k] for k in c)


def test_rank_guides_emits_citations():
    result = rank_guides(
        guides=[{"protospacer": "ATGCATGCATGCATGCATGG"}],
        reference="ATGCATGCATGCATGCATGGAGG" + "C" * 50,
        nuclease="cas9",
    )
    _assert_well_formed_citations(result["citations"])


def test_predict_offtargets_emits_citations():
    result = predict_offtargets(
        protospacer="ATGCATGCATGCATGCATGG",
        reference="ATGCATGCATGCATGCATGGAGG" + "C" * 50,
        nuclease="cas9",
    )
    _assert_well_formed_citations(result["citations"])


def test_predict_editing_efficiency_emits_citations():
    result = predict_editing_efficiency(
        protospacer="ATGCATGCATGCATGCATGG",
        pam="AGG",
        nuclease="cas9",
        delivery="rnp",
        outcome="nhej",
    )
    _assert_well_formed_citations(result["citations"])


def test_colony_calculator_emits_citations_for_preset():
    result = colony_calculator(preset="hdr_mammalian", desired_clones=1)
    _assert_well_formed_citations(result["citations"])


def test_colony_calculator_emits_citations_for_custom_efficiency():
    result = colony_calculator(editing_efficiency=0.5, desired_clones=1)
    _assert_well_formed_citations(result["citations"])


def test_interpret_ice_tide_cites_ice_paper_for_ice():
    result = interpret_ice_tide(editing_pct=80, r_squared=0.95, tool="ice")
    _assert_well_formed_citations(result["citations"])
    labels = " ".join(c["label"] for c in result["citations"])
    assert "Hsiau" in labels


def test_interpret_ice_tide_cites_tide_paper_for_tide():
    result = interpret_ice_tide(editing_pct=80, r_squared=0.95, tool="tide")
    _assert_well_formed_citations(result["citations"])
    labels = " ".join(c["label"] for c in result["citations"])
    assert "Brinkman" in labels


def test_verify_edit_emits_citations():
    reference = (
        "CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTT"
        "TAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
    )
    result = verify_edit(
        protospacer="TCAGAAACCTGCCAGTTTGC",
        reference=reference,
        nuclease="cas9",
    )
    _assert_well_formed_citations(result["citations"])
    labels = " ".join(c["label"] for c in result["citations"])
    # Cas9 mechanism + primer Tm + ICE/TIDE protocol papers should all appear
    assert "Jinek" in labels
    assert "Wallace" in labels


