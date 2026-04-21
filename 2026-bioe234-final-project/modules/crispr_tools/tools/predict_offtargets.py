from __future__ import annotations

from typing import Optional


def _reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[b] for b in reversed(seq))


class PredictOfftargets:
    """
    Description:
        Scans a reference DNA sequence for potential CRISPR off-target sites.
        Takes the 20bp protospacer from the gRNA design step and slides it
        across every position of the reference on both strands counting
        how many bases don't match at each position.

        Any site with mismatches at or below the threshold (default <= 3) gets
        flagged. Logic: if guide is similar enough to another region of the genome,
        it can cut there too, which we do not want.

        The scoring uses seed-region logic — the "seed region" is positions 1-12
        counting from the PAM end of the guide (Hsu et al. 2013). Mismatches
        there are much more dangerous than mismatches at the far end, because
        that's where Cas9 first contacts the DNA.

        Also checks for an NGG PAM after each candidate site, since without a
        PAM Cas9 basically can't cut.

        The reference can be passed as a resource name like "pBR322", a raw
        sequence, FASTA, or GenBank — the framework handles the conversion.

        Supports circular references (e.g. plasmids) via is_circular=True,
        which wraps the sequence before scanning so sites spanning the origin
        are not missed.

    Input:
        protospacer (str): the 20bp DNA protospacer sequence WITHOUT the PAM.
                           e.g. "TCAGAAACCTGCCAGTTTGC"
                           only A, T, G, C — no ambiguous bases here.
        reference (str):   the DNA sequence to scan for off-target sites.
                           can be a resource name, raw string, FASTA, or GenBank.
        max_mismatches (int): how many mismatches to still flag as a potential
                              off-target. default is 3.
        is_circular (bool): if True, wraps the reference before scanning so
                            sites spanning the origin are captured. default False.

    Output:
        dict with keys:
            - protospacer: the input protospacer (echoed back for reference)
            - reference_length: how long the reference was
            - strands_scanned: always 2 (forward + reverse complement)
            - sites_evaluated: total windows checked across both strands
            - offtarget_sites: list of dicts, one per flagged site, each with:
                - position: where in the reference (0-indexed, forward strand)
                - strand: "+" forward or "-" reverse
                - sequence: the 20bp window at that position
                - mismatches: how many bases differ from the protospacer
                - mismatch_positions: which positions (1 = PAM-proximal end)
                - seed_mismatches: mismatches in the dangerous positions 1-12
                - has_pam: True if NGG follows this site
                - risk: "HIGH", "MEDIUM", or "LOW"
            - high_risk_count: how many HIGH risk sites were found
            - specificity_summary: one sentence summary of the result

    Tests:
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="ATGATGATGATGATGATGATAGG"
            Expected Output: offtarget_sites has at least one entry with mismatches=0, has_pam=True
            Description: exact match + NGG PAM on forward strand should be flagged HIGH risk.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="CCCCCCCCCCCCCCCCCCCCAGG"
            Expected Output: offtarget_sites is empty or all entries have mismatches > 0
            Description: no match → empty list or only high-mismatch sites.
        - Case:
            Input: protospacer="", reference="ATGATGATG"
            Expected Exception: ValueError
            Description: empty protospacer should raise an error.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference=""
            Expected Exception: ValueError
            Description: empty reference should raise an error.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="ATGATGATGATGATGATGATAGG", max_mismatches=0
            Expected Output: only the exact-match site is returned
            Description: max_mismatches=0 means only perfect matches come through.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="CCCCCCCCCCCCCCCCATGA", max_mismatches=3
            Expected Output: site at position 19 has has_pam=False
            Description: protospacer at last valid window with no room for PAM.
    """

    _seed_region_length: int

    def initiate(self) -> None:
        # positions 1-12 from the PAM end are the seed region per Hsu et al. 2013
        self._seed_region_length = 12

    def run(
        self,
        protospacer: str,
        reference: str,
        max_mismatches: int = 3,
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

        invalid_ref = [b for b in set(reference) if b not in "ATGCN"]
        if invalid_ref:
            raise ValueError(
                f"Invalid base(s) in reference: {sorted(invalid_ref)}. "
                "Only standard bases A, T, G, C (and N) are accepted."
            )

        guide_len = len(protospacer)
        ref_len = len(reference)

        # for circular references, append enough bases to catch sites spanning the origin
        wrap_len = guide_len + 3  # protospacer + PAM
        scan_ref = reference + reference[:wrap_len] if is_circular else reference

        offtarget_sites = []
        sites_evaluated = 0

        for strand, seq in [("+", scan_ref), ("-", _reverse_complement(scan_ref))]:
            for i in range(len(seq) - guide_len + 1):
                window = seq[i : i + guide_len]

                if "N" in window:
                    continue

                sites_evaluated += 1

                mismatches = sum(1 for a, b in zip(protospacer, window) if a != b)

                if mismatches > max_mismatches:
                    continue

                # position numbering: 1 = PAM-proximal end (most dangerous)
                mismatch_positions = []
                for j in range(guide_len):
                    pam_proximal_idx = guide_len - 1 - j
                    if protospacer[pam_proximal_idx] != window[pam_proximal_idx]:
                        mismatch_positions.append(j + 1)

                seed_mismatches = sum(
                    1 for p in mismatch_positions if p <= self._seed_region_length
                )

                pam_start = i + guide_len
                has_pam = False
                if pam_start + 2 < len(seq):
                    pam_candidate = seq[pam_start : pam_start + 3]
                    if len(pam_candidate) == 3 and pam_candidate[1] == "G" and pam_candidate[2] == "G":
                        has_pam = True

                if mismatches == 0:
                    risk = "HIGH"
                elif seed_mismatches == 0 and has_pam:
                    risk = "HIGH"
                elif seed_mismatches <= 1 and has_pam:
                    risk = "MEDIUM"
                elif mismatches <= 2 and has_pam:
                    risk = "MEDIUM"
                else:
                    risk = "LOW"

                # convert back to forward-strand coordinates in the original reference;
                # for circular sequences, positions >= ref_len wrap back to 0-based origin
                if strand == "+":
                    reported_position = i % ref_len
                else:
                    reported_position = (ref_len - i - guide_len) % ref_len

                offtarget_sites.append({
                    "position": reported_position,
                    "strand": strand,
                    "sequence": window,
                    "mismatches": mismatches,
                    "mismatch_positions": mismatch_positions,
                    "seed_mismatches": seed_mismatches,
                    "has_pam": has_pam,
                    "risk": risk,
                })

        offtarget_sites.sort(key=lambda s: (s["mismatches"], -s["seed_mismatches"]))

        high_risk_count = sum(1 for s in offtarget_sites if s["risk"] == "HIGH")

        if not offtarget_sites:
            specificity_summary = (
                f"No off-target sites found within {max_mismatches} mismatches, "
                "the guide appears highly specific for this reference."
            )
        elif high_risk_count == 1 and offtarget_sites[0]["mismatches"] == 0:
            specificity_summary = (
                f"1 on-target site found (0 mismatches). "
                f"{len(offtarget_sites) - 1} additional site(s) within "
                f"{max_mismatches} mismatches — review before proceeding."
            )
        else:
            specificity_summary = (
                f"{len(offtarget_sites)} potential off-target site(s) found "
                f"within {max_mismatches} mismatches ({high_risk_count} HIGH risk). "
                "Consider redesigning the guide or choosing a more specific protospacer."
            )

        return {
            "protospacer": protospacer,
            "reference_length": ref_len,
            "strands_scanned": 2,
            "sites_evaluated": sites_evaluated,
            "offtarget_sites": offtarget_sites,
            "high_risk_count": high_risk_count,
            "specificity_summary": specificity_summary,
        }


_instance = PredictOfftargets()
_instance.initiate()
predict_offtargets = _instance.run
