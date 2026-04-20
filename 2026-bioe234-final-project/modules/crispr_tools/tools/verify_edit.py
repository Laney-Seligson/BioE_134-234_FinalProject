from __future__ import annotations

from typing import Optional


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


class VerifyEdit:
    """
    Description:
        After the CRISPR experiment is done, this tool helps the user figure out
        where Cas9 cut and gives them everything they need to verify whether the
        edit actually worked.

        Cas9 always cuts at the same predictable position: between nucleotides 17
        and 18 of the protospacer (3bp upstream of the PAM). Given the protospacer
        and reference sequence, this tool calculates the exact cut site, designs
        two sequencing primers flanking it (~150bp each side), and returns a
        step-by-step ICE/TIDE protocol.

        ICE (Synthego) and TIDE (Brinkman et al. 2014) both work by detecting the
        mixed/noisy Sanger signal that appears downstream of the cut when edited
        and unedited alleles are sequenced together.

        Reference input can be a resource name like "pBR322", a raw sequence,
        FASTA, or GenBank — the framework resolves it automatically.

        Supports circular references via is_circular=True so the protospacer
        can be found even when it spans the origin of a plasmid.

    Input:
        protospacer (str): the 20bp DNA protospacer (no PAM) used in the edit.
                           e.g. "TCAGAAACCTGCCAGTTTGC"
        reference (str):   the original unedited reference sequence. accepts
                           resource name, raw string, FASTA, or GenBank.
        primer_offset (int): how far from the cut site to place the primers.
                             default 150bp.
        primer_len (int): length of the sequencing primers. default 20bp.
        is_circular (bool): if True, wraps the reference before searching so
                            protospacers spanning the origin are found. default False.

    Output:
        dict with these keys:
            - protospacer: the input protospacer
            - strand: "+" forward or "-" reverse
            - protospacer_position: where the protospacer starts in the reference (0-indexed)
            - pam_sequence: the 3bp PAM found after the protospacer (should be NGG)
            - cut_position: where Cas9 cuts (between nt 17-18 of protospacer, 0-indexed)
            - forward_primer: sequence of the upstream sequencing primer
            - forward_primer_position: where that primer starts
            - reverse_primer: sequence of the downstream sequencing primer
            - reverse_primer_position: where that primer starts
            - amplicon_sequence: the reference sequence between the two primers
            - amplicon_length: how long that amplicon is in bp
            - cut_offset_in_amplicon: where the cut falls within the amplicon
            - interpretation_guide: step-by-step ICE/TIDE instructions

    Tests:
        - Case:
            Input: protospacer="TCAGAAACCTGCCAGTTTGC", reference="CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
            Expected Output: cut_position=32, pam_sequence="TGG", strand="+"
            Description: standard forward-strand case.
        - Case:
            Input: protospacer="GCAAACTGGCAGGTTTCTGA", reference="CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACG..."
            Expected Output: strand="-"
            Description: reverse-complement protospacer — tool should find it on the minus strand.
        - Case:
            Input: protospacer="AAAAAAAAAAAAAAAAAAAA", reference="ATGCATGCATGC"
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
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        protospacer: str,
        reference: str,
        primer_offset: int = 150,
        primer_len: int = 20,
        is_circular: bool = False,
    ) -> dict:

        protospacer = protospacer.upper().strip()
        reference = reference.upper().strip()

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

        # for circular references, wrap enough bases to catch sites spanning the origin
        search_ref = reference + reference[:len(protospacer) + 3] if is_circular else reference
        ps_position, strand = _find_protospacer_in_reference(protospacer, search_ref)

        # keep reported position within the original reference length
        ref_len = len(reference)
        ps_position = ps_position % ref_len

        # verify PAM; on the minus strand the NGG appears as CCN before the RC protospacer
        if strand == "+":
            pam_start = ps_position + len(protospacer)
            if pam_start + 3 > ref_len:
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
            # on the minus strand, CCN (= RC of NGG) sits just left of the RC protospacer
            pam_start_fwd = ps_position - 3
            if pam_start_fwd < 0:
                raise ValueError(
                    "Protospacer is at the very start of the reference — no room for a PAM. "
                    "Provide a longer reference sequence."
                )
            pam_fwd = reference[pam_start_fwd : ps_position]
            pam_sequence = _reverse_complement(pam_fwd)

        # Cas9 cuts between nt 17 and 18, i.e. 3bp upstream of the PAM (well established)
        if strand == "+":
            cut_position = ps_position + 17
        else:
            cut_position = ps_position + len(protospacer) - 17

        cut_position = max(0, min(cut_position, ref_len - 1))

        fwd_primer, fwd_primer_pos = _design_sequencing_primer(
            reference, cut_position, "forward", primer_len=primer_len, offset=primer_offset
        )
        rev_primer, rev_primer_pos = _design_sequencing_primer(
            reference, cut_position, "reverse", primer_len=primer_len, offset=primer_offset
        )

        amplicon_start = fwd_primer_pos
        amplicon_end = min(rev_primer_pos + primer_len, ref_len)
        amplicon_sequence = reference[amplicon_start:amplicon_end]
        amplicon_length = len(amplicon_sequence)
        cut_offset_in_amplicon = cut_position - amplicon_start

        interpretation_guide = (
            f"Edit Verification Protocol (ICE/TIDE):\n"
            f"\n"
            f"1. PCR: Amplify the target locus using the forward primer "
            f"(pos {fwd_primer_pos}) and reverse primer (pos {rev_primer_pos}). "
            f"Expected amplicon size: {amplicon_length} bp.\n"
            f"\n"
            f"2. Sanger sequencing: Submit the PCR amplicon for Sanger sequencing "
            f"using either the forward or reverse sequencing primer. The cut site "
            f"is {cut_offset_in_amplicon} bp from the amplicon start "
            f"(absolute position {cut_position} in reference).\n"
            f"\n"
            f"3. ICE analysis (Synthego): Upload the Sanger .ab1 trace file and "
            f"the reference amplicon sequence. ICE will decompose the mixed trace "
            f"signal downstream of the cut site into indel alleles and report "
            f"editing efficiency (KO score) and the indel spectrum.\n"
            f"\n"
            f"4. TIDE analysis (Brinkman et al. 2014): Upload the edited .ab1 trace "
            f"and an unedited control trace. TIDE fits a linear model to the "
            f"sequence downstream of the cut site to quantify indel frequencies.\n"
            f"\n"
            f"5. Expected result for successful editing: Trace signal becomes "
            f"mixed or noisy at position ~{cut_offset_in_amplicon} bp in the amplicon. "
            f"ICE/TIDE should report indel efficiency > 0% if editing occurred. "
            f"The most common indels are +1 insertions and small deletions (1-10 bp) "
            f"centered on the cut site.\n"
            f"\n"
            f"6. Negative control: Use an unedited sample to confirm clean trace "
            f"signal. No mixed peaks should appear at the cut site."
        )

        return {
            "protospacer": protospacer,
            "strand": strand,
            "protospacer_position": ps_position,
            "pam_sequence": pam_sequence,
            "cut_position": cut_position,
            "forward_primer": fwd_primer,
            "forward_primer_position": fwd_primer_pos,
            "reverse_primer": rev_primer,
            "reverse_primer_position": rev_primer_pos,
            "amplicon_sequence": amplicon_sequence,
            "amplicon_length": amplicon_length,
            "cut_offset_in_amplicon": cut_offset_in_amplicon,
            "interpretation_guide": interpretation_guide,
        }


_instance = VerifyEdit()
_instance.initiate()
verify_edit = _instance.run
