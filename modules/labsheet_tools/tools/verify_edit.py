from __future__ import annotations

from typing import Optional

from modules.crispr_tools.tools.citations import cites, format_citations

_SUPPORTED_NUCLEASES = {"cas9", "cas12a"}


def _reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[b] for b in reversed(seq))


def _find_protospacer_in_reference(protospacer: str, reference: str) -> tuple[int, str]:
    pos = reference.find(protospacer)
    if pos != -1:
        return pos, "+"

    # guide targets the minus strand — protospacer appears as its RC in the forward sequence
    rc = _reverse_complement(protospacer)
    pos = reference.find(rc)
    if pos != -1:
        return pos, "-"

    raise ValueError(
        f"Protospacer not found in reference sequence on either strand. "
        "Ensure the protospacer was designed from this reference."
    )


def _design_sequencing_primer(
    reference: str,
    cut_position: int,
    direction: str,
    primer_len: int = 20,
    offset: int = 150,
) -> tuple[str, int]:
    if direction == "forward":
        start = max(0, cut_position - offset)
        end = start + primer_len
        if end > len(reference):
            end = len(reference)
            start = max(0, end - primer_len)
        primer = reference[start:end]
        return primer, start
    else:
        start = min(len(reference), cut_position + offset)
        end = start + primer_len
        if end > len(reference):
            end = len(reference)
        fwd_seq = reference[start:end]
        primer = _reverse_complement(fwd_seq)
        return primer, start


def _calculate_tm_wallace(primer: str) -> float:
    """
    Estimate primer melting temperature using the Wallace rule:
        Tm = 2*(A+T) + 4*(G+C)  [degrees Celsius]

    Source: Wallace et al. 1979, Nucleic Acids Res 6(11):3543-3557.

    Valid for short oligos (~14-20 bp) at standard salt conditions. For longer
    primers this tends to overestimate Tm, but it's the standard quick-check
    used in undergraduate molecular biology labs and is good enough to flag
    mismatched primer pairs.
    """
    at = sum(1 for b in primer if b in "AT")
    gc = sum(1 for b in primer if b in "GC")
    return 2 * at + 4 * gc


def _evaluate_primers(fwd_primer: str, rev_primer: str) -> tuple[float, float, list[str]]:
    """
    Calculate Tm for each primer and return warnings for common PCR failure modes:
      - Tm outside 50-65 degC range (too unstable or too sticky)
      - |Tm_fwd - Tm_rev| > 5 degC (one primer will dominate; asymmetric amplification)
      - GC content outside 40-60% (hard to get clean amplification)
      - Poly-N runs of 4+ (slippage, secondary structure)
    Sources: Dieffenbach et al. 1993, PCR Methods Appl 3(3):S30-37;
             Rychlik 1995, Mol Biotechnol 3(2):129-134.
    """
    tm_fwd = _calculate_tm_wallace(fwd_primer)
    tm_rev = _calculate_tm_wallace(rev_primer)
    warnings: list[str] = []

    for label, primer, tm in [("forward", fwd_primer, tm_fwd), ("reverse", rev_primer, tm_rev)]:
        if tm < 50:
            warnings.append(
                f"{label} primer Tm is {tm:.1f}degC (< 50degC): may not anneal reliably. "
                "Consider a longer primer or a more GC-rich region."
            )
        elif tm > 65:
            warnings.append(
                f"{label} primer Tm is {tm:.1f}degC (> 65degC): may cause non-specific binding. "
                "Consider a shorter primer or a less GC-rich region."
            )

        gc_frac = sum(1 for b in primer if b in "GC") / len(primer) if primer else 0
        if gc_frac < 0.40:
            warnings.append(f"{label} primer GC content {gc_frac:.0%} is below 40% — weak binding.")
        elif gc_frac > 0.60:
            warnings.append(f"{label} primer GC content {gc_frac:.0%} is above 60% — risk of secondary structure.")

        for base in "ATGC":
            if base * 4 in primer:
                warnings.append(
                    f"{label} primer contains a poly-{base} run of 4+: risk of slippage/mispriming."
                )
                break

    if abs(tm_fwd - tm_rev) > 5:
        warnings.append(
            f"Primer Tm mismatch: forward {tm_fwd:.1f}degC vs reverse {tm_rev:.1f}degC "
            f"(delta {abs(tm_fwd - tm_rev):.1f}degC > 5degC). One primer will dominate at "
            "the annealing step, leading to asymmetric PCR. Redesign so Tms are within 5degC."
        )

    return tm_fwd, tm_rev, warnings


def _validate_pam(
    reference: str, ps_position: int, strand: str, protospacer_len: int, nuclease: str
) -> str:
    """Validate PAM and return its sequence. Raises ValueError if PAM is absent or wrong."""
    if nuclease == "cas9":
        # NGG PAM sits immediately 3' of the protospacer (Jinek et al. 2012)
        if strand == "+":
            pam_start = ps_position + protospacer_len
            if pam_start + 3 > len(reference):
                raise ValueError(
                    "Protospacer is at the very end of the reference — no room for a PAM. "
                    "Provide a longer reference sequence."
                )
            pam_sequence = reference[pam_start : pam_start + 3]
            if pam_sequence[1] != "G" or pam_sequence[2] != "G":
                raise ValueError(
                    f"Expected NGG PAM after protospacer but found '{pam_sequence}'. "
                    "Ensure the protospacer was designed with an NGG PAM."
                )
        else:
            # on the minus strand, CCN on the forward strand (= RC of NGG) sits
            # immediately before the RC protospacer
            pam_start_fwd = ps_position - 3
            if pam_start_fwd < 0:
                raise ValueError(
                    "Protospacer is at the very start of the reference — no room for a PAM. "
                    "Provide a longer reference sequence."
                )
            pam_fwd = reference[pam_start_fwd : ps_position]
            pam_sequence = _reverse_complement(pam_fwd)
        return pam_sequence

    else:  # cas12a
        # TTTV PAM (V = A/C/G) sits immediately 5' of the protospacer on the target strand.
        # (Zetsche et al. 2015, Cell 163(3):759-771; Fonfara et al. 2016, Nature)
        if strand == "+":
            if ps_position < 4:
                raise ValueError(
                    "Protospacer is at the very start of the reference — no room for a TTTV PAM. "
                    "Provide a longer reference sequence."
                )
            pam_fwd = reference[ps_position - 4 : ps_position]
            if pam_fwd[:3] != "TTT" or pam_fwd[3] not in "ACG":
                raise ValueError(
                    f"Expected TTTV PAM before protospacer but found '{pam_fwd}'. "
                    "Ensure the protospacer was designed with a TTTV PAM for Cas12a."
                )
            return pam_fwd
        else:
            # On the minus strand the TTTV PAM appears as BAAA (B = C/G/T = RC of V)
            # on the forward strand immediately after the RC protospacer.
            pam_end = ps_position + protospacer_len + 4
            if pam_end > len(reference):
                raise ValueError(
                    "Protospacer is at the very end of the reference — no room for a TTTV PAM. "
                    "Provide a longer reference sequence."
                )
            pam_fwd = reference[ps_position + protospacer_len : ps_position + protospacer_len + 4]
            pam_sequence = _reverse_complement(pam_fwd)
            if pam_sequence[:3] != "TTT" or pam_sequence[3] not in "ACG":
                raise ValueError(
                    f"Expected TTTV PAM on minus strand but found '{pam_sequence}' "
                    f"(RC of '{pam_fwd}'). "
                    "Ensure the protospacer was designed with a TTTV PAM for Cas12a."
                )
            return pam_sequence


def _compute_cut_position(
    ps_position: int, strand: str, protospacer_len: int, nuclease: str
) -> int:
    if nuclease == "cas9":
        # SpCas9 cuts between nt 17 and 18, i.e. 3bp upstream of the NGG PAM.
        # (Jinek et al. 2012, Science 337:816-821)
        if strand == "+":
            return ps_position + 17
        else:
            return ps_position + protospacer_len - 17
    else:  # cas12a
        # LbCas12a makes a staggered cut: non-template strand between positions 18-19
        # (counting from the PAM-proximal end), template strand between 23-24,
        # producing a 5-nt 5' overhang. We report the non-template strand cut.
        # (Zetsche et al. 2015, Cell; Fonfara et al. 2016, Nature 532:517-521)
        if strand == "+":
            return ps_position + 18
        else:
            return ps_position + protospacer_len - 18


class VerifyEdit:
    """
    Description:
        After the CRISPR experiment is done, this tool helps the user figure out
        where Cas9 or Cas12a cut and gives them everything they need to verify
        whether the edit actually worked.

        For SpCas9: cuts between nt 17 and 18 of the protospacer (3bp upstream
        of the NGG PAM), producing a blunt double-strand break. PAM is 3' of the
        protospacer. (Jinek et al. 2012, Science 337:816-821)

        For LbCas12a: makes a staggered cut between positions 18-19 on the
        non-template strand and 23-24 on the template strand (5-nt 5' overhang).
        PAM is TTTV (V = A/C/G) and sits 5' of the 23bp protospacer.
        (Zetsche et al. 2015, Cell 163(3):759-771; Fonfara et al. 2016, Nature)

        Given the protospacer and reference, this tool calculates the exact cut
        site, designs two sequencing primers flanking it (~150bp each side), and
        returns a step-by-step ICE/TIDE protocol.

        ICE (Synthego) and TIDE (Brinkman et al. 2014) both work by detecting the
        mixed/noisy Sanger signal that appears downstream of the cut when edited
        and unedited alleles are sequenced together.

        Reference input can be a resource name like "pBR322", a raw sequence,
        FASTA, or GenBank — the framework resolves it automatically.

        Supports circular references via is_circular=True so the protospacer
        can be found even when it spans the origin of a plasmid.

    Input:
        protospacer (str): the DNA protospacer (no PAM) used in the edit.
                           20bp for Cas9, 23bp for Cas12a.
                           e.g. "TCAGAAACCTGCCAGTTTGC" (Cas9)
        reference (str):   the original unedited reference sequence. accepts
                           resource name, raw string, FASTA, or GenBank.
        nuclease (str):    "cas9" (default) or "cas12a". Determines PAM
                           identity and cut site calculation.
        primer_offset (int): how far from the cut site to place the primers.
                             default 150bp.
        primer_len (int): length of the sequencing primers. default 20bp.
        is_circular (bool): if True, wraps the reference before searching so
                            protospacers spanning the origin are found. default False.

    Output:
        dict with these keys:
            - protospacer: the input protospacer
            - nuclease: which nuclease was used
            - strand: "+" forward or "-" reverse
            - protospacer_position: where the protospacer starts in the reference (0-indexed)
            - pam_sequence: the PAM found adjacent to the protospacer
            - cut_position: where the nuclease cuts (0-indexed, non-template strand for Cas12a)
            - forward_primer: sequence of the upstream sequencing primer
            - forward_primer_position: where that primer starts
            - forward_primer_tm: estimated Tm in degC (Wallace rule)
            - reverse_primer: sequence of the downstream sequencing primer
            - reverse_primer_position: where that primer starts
            - reverse_primer_tm: estimated Tm in degC (Wallace rule)
            - tm_difference: |Tm_fwd - Tm_rev| in degC (good primer pairs have <5)
            - primer_warnings: list of QC warnings (Tm out of range, Tm mismatch,
                               GC content, poly-N runs). Empty list = primers are clean.
            - amplicon_sequence: the reference sequence between the two primers
            - amplicon_length: how long that amplicon is in bp
            - cut_offset_in_amplicon: where the cut falls within the amplicon
            - interpretation_guide: step-by-step ICE/TIDE instructions

    Tests:
        - Case:
            Input: protospacer="TCAGAAACCTGCCAGTTTGC", reference="CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC", nuclease="cas9"
            Expected Output: cut_position=32, pam_sequence="TGG", strand="+"
            Description: standard Cas9 forward-strand case.
        - Case:
            Input: protospacer="GCAAACTGGCAGGTTTCTGA", reference="CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACG...", nuclease="cas9"
            Expected Output: strand="-"
            Description: reverse-complement Cas9 protospacer — tool finds it on minus strand.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGATGAT", reference="TTTAATGATGATGATGATGATGATGATAAAAAAAAAAAA", nuclease="cas12a"
            Expected Output: cut_position=22, pam_sequence="TTTA", strand="+"
            Description: Cas12a forward-strand case; TTTA PAM before 23bp protospacer, cut at pos 4+18=22.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGATGAT", reference="AAAAAAAAAAAATCATCATCATCATCATCATCATTAAAAAAAAAAA", nuclease="cas12a"
            Expected Output: strand="-"
            Description: Cas12a minus-strand case; TTTA PAM on minus strand (appears as TAAA on forward).
        - Case:
            Input: protospacer="AAAAAAAAAAAAAAAAAAAA", reference="ATGCATGCATGC", nuclease="cas9"
            Expected Exception: ValueError
            Description: protospacer not in reference → error.
        - Case:
            Input: protospacer="", reference="ATGCATGCATGC"
            Expected Exception: ValueError
            Description: empty protospacer → error.
        - Case:
            Input: protospacer="TCAGAAACCTGCCAGTTTGC", reference=""
            Expected Exception: ValueError
            Description: empty reference → error.
        - Case:
            Input: protospacer="TCAGAAACCTGCCAGTTTGC", reference="ATGCATGCATGC", nuclease="talenuclease"
            Expected Exception: ValueError
            Description: unsupported nuclease → error.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        protospacer: str,
        reference: str,
        nuclease: str = "cas9",
        primer_offset: int = 150,
        primer_len: int = 20,
        is_circular: bool = False,
    ) -> dict:

        protospacer = protospacer.upper().strip()
        reference = reference.upper().strip()
        nuclease = nuclease.lower().strip()

        if nuclease not in _SUPPORTED_NUCLEASES:
            raise ValueError(
                f"Unsupported nuclease '{nuclease}'. Choose 'cas9' or 'cas12a'."
            )
        if not protospacer:
            raise ValueError("Protospacer must not be empty.")
        if not reference:
            raise ValueError("Reference sequence must not be empty.")

        invalid_ps = [b for b in set(protospacer) if b not in "ATGC"]
        if invalid_ps:
            raise ValueError(
                f"Invalid base(s) in protospacer: {sorted(invalid_ps)}. "
                "Only standard bases A, T, G, C are accepted."
            )

        protospacer_len = len(protospacer)
        ref_len = len(reference)

        # for circular references, wrap enough bases to catch sites spanning the origin
        pam_len = 3 if nuclease == "cas9" else 4
        search_ref = reference + reference[: protospacer_len + pam_len] if is_circular else reference
        ps_position, strand = _find_protospacer_in_reference(protospacer, search_ref)
        ps_position = ps_position % ref_len

        pam_sequence = _validate_pam(reference, ps_position, strand, protospacer_len, nuclease)
        cut_position = _compute_cut_position(ps_position, strand, protospacer_len, nuclease)
        cut_position = max(0, min(cut_position, ref_len - 1))

        fwd_primer, fwd_primer_pos = _design_sequencing_primer(
            reference, cut_position, "forward", primer_len=primer_len, offset=primer_offset
        )
        rev_primer, rev_primer_pos = _design_sequencing_primer(
            reference, cut_position, "reverse", primer_len=primer_len, offset=primer_offset
        )

        tm_fwd, tm_rev, primer_warnings = _evaluate_primers(fwd_primer, rev_primer)

        amplicon_start = fwd_primer_pos
        amplicon_end = min(rev_primer_pos + primer_len, ref_len)
        amplicon_sequence = reference[amplicon_start:amplicon_end]
        amplicon_length = len(amplicon_sequence)
        cut_offset_in_amplicon = cut_position - amplicon_start

        if nuclease == "cas9":
            cut_description = (
                f"SpCas9 makes a blunt cut between nt 17 and 18 (3bp upstream of the NGG PAM). "
                f"Cut site is {cut_offset_in_amplicon} bp from the amplicon start "
                f"(absolute position {cut_position} in reference)."
            )
        else:
            cut_description = (
                f"LbCas12a makes a staggered cut: non-template strand between positions 18-19, "
                f"template strand between positions 23-24, producing a 5-nt 5' overhang. "
                f"Reported cut (non-template strand) is {cut_offset_in_amplicon} bp from "
                f"the amplicon start (absolute position {cut_position} in reference)."
            )

        interpretation_guide = (
            f"Edit Verification Protocol (ICE/TIDE) — {nuclease.upper()}:\n"
            f"\n"
            f"1. PCR: Amplify the target locus using the forward primer "
            f"(pos {fwd_primer_pos}) and reverse primer (pos {rev_primer_pos}). "
            f"Expected amplicon size: {amplicon_length} bp.\n"
            f"\n"
            f"2. Cut site: {cut_description}\n"
            f"\n"
            f"3. Sanger sequencing: Submit the PCR amplicon for Sanger sequencing "
            f"using either the forward or reverse sequencing primer.\n"
            f"\n"
            f"4. ICE analysis (Synthego): Upload the Sanger .ab1 trace file and "
            f"the reference amplicon sequence. ICE will decompose the mixed trace "
            f"signal downstream of the cut site into indel alleles and report "
            f"editing efficiency (KO score) and the indel spectrum.\n"
            f"\n"
            f"5. TIDE analysis (Brinkman et al. 2014): Upload the edited .ab1 trace "
            f"and an unedited control trace. TIDE fits a linear model to the "
            f"sequence downstream of the cut site to quantify indel frequencies.\n"
            f"\n"
            f"6. Expected result for successful editing: Trace signal becomes "
            f"mixed or noisy at position ~{cut_offset_in_amplicon} bp in the amplicon. "
            f"ICE/TIDE should report indel efficiency > 0% if editing occurred. "
            f"The most common indels are +1 insertions and small deletions (1-10 bp) "
            f"centered on the cut site.\n"
            f"\n"
            f"7. Negative control: Use an unedited sample to confirm clean trace "
            f"signal. No mixed peaks should appear at the cut site."
        )

        # Cite the cut-site mechanism papers (which justify cut_position),
        # the primer-design rules (Wallace Tm + Dieffenbach QC), and the
        # ICE/TIDE protocol papers (which justify the interpretation_guide).
        if nuclease == "cas9":
            citation_keys = ["jinek_2012"]
        else:
            citation_keys = ["zetsche_2015", "fonfara_2016"]
        citation_keys.extend([
            "wallace_1979",
            "dieffenbach_1993",
            "hsiau_2019_ice",
            "brinkman_2014_tide",
        ])

        return {
            "protospacer": protospacer,
            "nuclease": nuclease,
            "strand": strand,
            "protospacer_position": ps_position,
            "pam_sequence": pam_sequence,
            "cut_position": cut_position,
            "forward_primer": fwd_primer,
            "forward_primer_position": fwd_primer_pos,
            "forward_primer_tm": round(tm_fwd, 1),
            "reverse_primer": rev_primer,
            "reverse_primer_position": rev_primer_pos,
            "reverse_primer_tm": round(tm_rev, 1),
            "tm_difference": round(abs(tm_fwd - tm_rev), 1),
            "primer_warnings": primer_warnings,
            "amplicon_sequence": amplicon_sequence,
            "amplicon_length": amplicon_length,
            "cut_offset_in_amplicon": cut_offset_in_amplicon,
            "interpretation_guide": interpretation_guide,
            "citations": format_citations(cites(*citation_keys)),
        }


_instance = VerifyEdit()
_instance.initiate()
verify_edit = _instance.run
