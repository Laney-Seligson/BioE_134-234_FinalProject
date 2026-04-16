from __future__ import annotations

from typing import Optional


# same helper as in predict_offtargets.
# keeping it as a module-level function so the class stays clean.
def _reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[b] for b in reversed(seq))


def _find_protospacer_in_reference(protospacer: str, reference: str) -> tuple[int, str]:
    # try forward strand first — just a simple string search
    pos = reference.find(protospacer)
    if pos != -1:
        return pos, "+"

    # if not found, try the reverse complement of the protospacer.
    # this handles the case where the guide targets the minus strand —
    # the protospacer itself would appear as its RC in the forward sequence.
    rc = _reverse_complement(protospacer)
    pos = reference.find(rc)
    if pos != -1:
        return pos, "-"

    # if still not found, something is wrong — guide and reference don't match
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
    # for the forward primer: go upstream (left) of the cut site by ~offset bp,
    # then take primer_len bases going right. this gives a primer that reads
    # toward the cut site during sequencing.
    if direction == "forward":
        start = max(0, cut_position - offset)  # don't go below position 0
        end = start + primer_len
        # handle edge case if we're near the end of the sequence
        if end > len(reference):
            end = len(reference)
            start = max(0, end - primer_len)
        primer = reference[start:end]
        return primer, start

    else:  # reverse primer
        # go downstream (right) of the cut site by ~offset bp, take primer_len
        # bases, then reverse complement — so it reads back toward the cut site
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
        where Cas9 cut and gives you everything they need to verify whether the
        edit actually worked.

        The idea: Cas9 always cuts at the same predictable position — between
        nucleotides 17 and 18 of the protospacer, counting from the PAM-distal
        end. That's 3bp upstream of the PAM. So if I know the protospacer and
        the reference sequence, I can calculate the exact cut site.

        Then I design two sequencing primers flanking that cut site (~150bp on
        each side). You PCR-amplify that region from your edited cells, send it
        for Sanger sequencing, and upload the trace to ICE (Synthego) or TIDE
        (Brinkman et al. 2014). Those tools look for the characteristic "noisy"
        signal that appears downstream of the cut when a mix of edited and
        unedited cells is sequenced together.

        Steps this function does:
          1. finds the protospacer in the reference (forward or reverse strand)
          2. checks the NGG PAM is actually there after it
          3. calculates the cut position (between nt 17 and 18)
          4. designs a forward sequencing primer ~150bp upstream
          5. designs a reverse sequencing primer ~150bp downstream
          6. extracts the expected amplicon sequence between those primers
          7. returns a step-by-step ICE/TIDE protocol with the exact coordinates

        Reference input can be a resource name like "pBR322", a raw sequence,
        FASTA, or GenBank — the framework resolves it automatically.

    Input:
        protospacer (str): the 20bp DNA protospacer (no PAM) used in the edit.
                           e.g. "TCAGAAACCTGCCAGTTTGC"
        reference (str):   the original unedited reference sequence — used to
                           find the cut site and design primers. accepts resource
                           name, raw string, FASTA, or GenBank.
        primer_offset (int): how far from the cut site to place the primers.
                             default 150bp. if your reference is short (like a
                             test sequence), reduce this to like 10 or 20.
        primer_len (int): length of the sequencing primers. default 20bp.

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
            Expected Output: cut_position is an int >= 0, pam_sequence == "TGG"
            Description: standard BioE134 test sequence — should find protospacer at position 15, PAM TGG, cut at 32.
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
        # nothing to set up here — all the logic lives in run()
        pass

    def run(
        self,
        protospacer: str,
        reference: str,
        primer_offset: int = 150,
        primer_len: int = 20,
    ) -> dict:

        # clean up input
        protospacer = protospacer.upper().strip()
        reference = reference.upper().strip()

        if not protospacer:
            raise ValueError("Protospacer must not be empty.")
        if not reference:
            raise ValueError("Reference sequence must not be empty.")

        # protospacer should only have standard bases
        invalid_ps = [b for b in set(protospacer) if b not in "ATGC"]
        if invalid_ps:
            raise ValueError(
                f"Invalid base(s) in protospacer: {sorted(invalid_ps)}. "
                "Only standard bases A, T, G, C are accepted."
            )

        # step 1: find where the protospacer is in the reference sequence.
        # it could be on the forward strand or the reverse strand —
        # _find_protospacer_in_reference checks both and tells us which.
        ps_position, strand = _find_protospacer_in_reference(protospacer, reference)

        # step 2: verify the PAM is actually there.
        # for the forward strand: PAM is the 3 bases immediately after the protospacer.
        # for the reverse strand: it's more confusing — on the forward reference the
        # PAM would appear as CCN (the RC of NGG) in the 3 bases BEFORE the RC protospacer.
        if strand == "+":
            pam_start = ps_position + len(protospacer)
            if pam_start + 3 > len(reference):
                raise ValueError(
                    "Protospacer is at the very end of the reference — no room for a PAM. "
                    "Provide a longer reference sequence."
                )
            pam_sequence = reference[pam_start : pam_start + 3]
            # PAM must be xGG (NGG pattern — any base then two G's)
            if pam_sequence[1] != "G" or pam_sequence[2] != "G":
                raise ValueError(
                    f"Expected NGG PAM after protospacer but found '{pam_sequence}'. "
                    "Ensure the protospacer was designed with an NGG PAM."
                )
        else:
            # on the minus strand, the PAM sits just to the left of the RC protospacer
            # on the forward strand. it looks like CCN there, but converting it to RC
            # gives us the NGG representation.
            pam_start_fwd = ps_position - 3
            if pam_start_fwd < 0:
                raise ValueError(
                    "Protospacer is at the very start of the reference — no room for a PAM. "
                    "Provide a longer reference sequence."
                )
            pam_fwd = reference[pam_start_fwd : ps_position]
            pam_sequence = _reverse_complement(pam_fwd)  # convert CCN → NGG view

        # step 3: calculate the cut position.
        # Cas9 creates a blunt-end cut between positions 17 and 18 of the protospacer,
        # which is 3bp upstream of the PAM. this is well established in the literature.
        # on the forward strand: cut = start of protospacer + 17
        # on the reverse strand: the math flips because everything is RC'd
        if strand == "+":
            cut_position = ps_position + 17
        else:
            cut_position = ps_position + len(protospacer) - 17

        # clamp to valid range just in case (shouldn't happen with real sequences)
        cut_position = max(0, min(cut_position, len(reference) - 1))

        # step 4: design the two sequencing primers flanking the cut site.
        # forward primer goes upstream (left), reverse primer goes downstream (right).
        # the offset controls how far from the cut site they start —
        # 150bp is standard for Sanger sequencing window, but for short test sequences
        # you'd pass a smaller number.
        fwd_primer, fwd_primer_pos = _design_sequencing_primer(
            reference, cut_position, "forward", primer_len=primer_len, offset=primer_offset
        )
        rev_primer, rev_primer_pos = _design_sequencing_primer(
            reference, cut_position, "reverse", primer_len=primer_len, offset=primer_offset
        )

        # step 5: extract the amplicon — everything between the two primers.
        # also calculate where the cut falls inside the amplicon (needed for ICE/TIDE).
        amplicon_start = fwd_primer_pos
        amplicon_end = min(rev_primer_pos + primer_len, len(reference))
        amplicon_sequence = reference[amplicon_start:amplicon_end]
        amplicon_length = len(amplicon_sequence)
        cut_offset_in_amplicon = cut_position - amplicon_start

        # step 6: build the ICE/TIDE protocol with all the coordinates filled in.
        # ICE = Synthego's web tool, TIDE = the original Brinkman et al. 2014 method.
        # both work by looking at how the Sanger trace gets "noisy" at the cut site
        # when you have a mix of edited and unedited alleles.
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


# module-level alias — same pattern used in all the other tools in this project.
# means you can do either:
#   from modules.crispr_tools.tools.verify_edit import verify_edit
#   result = verify_edit(protospacer, reference)
# or use the class directly if you need more control.
_instance = VerifyEdit()
_instance.initiate()
verify_edit = _instance.run
