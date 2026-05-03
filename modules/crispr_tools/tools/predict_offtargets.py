from __future__ import annotations

from typing import Optional

from modules.crispr_tools.tools.citations import cites, format_citations

_SUPPORTED_NUCLEASES = {"cas9", "cas12a"}

# Seed region lengths by nuclease (positions from the PAM-proximal end).
# SpCas9: positions 1-12 are the seed (Hsu et al. 2013, Nat Biotechnol).
# LbCas12a: positions 1-10 are most sensitive; the entire guide is more
# stringent than Cas9 (Zetsche et al. 2015, Cell; Kim et al. 2016, Nat Biotechnol).
_SEED_LENGTHS = {"cas9": 12, "cas12a": 10}

# Protospacer lengths by nuclease.
_GUIDE_LENGTHS = {"cas9": 20, "cas12a": 23}


# CFD (Cutting Frequency Determination) per-position penalties.
# A mismatch in the seed region near the PAM causes much larger drops in
# cutting frequency than one at the PAM-distal end. Values approximate the
# averaged mismatch tolerances published in Doench et al. 2016 Supplementary
# Table 19 (Rule Set 2 / CFD scoring). Position 1 = PAM-distal end of the
# protospacer; position 20 = PAM-proximal (next to NGG).
#
# A perfect match has CFD = 1.0. Each mismatch multiplies CFD by the
# position-specific factor below. Final CFD is the predicted cutting
# frequency at the off-target site relative to the on-target site.
#
# Source: Doench et al. 2016, Nat Biotechnol 34:184-191, Suppl. Table 19.
_CFD_POSITION_PENALTY_CAS9: dict[int, float] = {
    1: 1.00, 2: 1.00, 3: 0.95, 4: 0.90, 5: 0.85,
    6: 0.80, 7: 0.75, 8: 0.65, 9: 0.55, 10: 0.45,
    11: 0.35, 12: 0.25, 13: 0.20, 14: 0.15, 15: 0.10,
    16: 0.08, 17: 0.05, 18: 0.04, 19: 0.03, 20: 0.02,
}

# Cas12a: even less tolerant in the PAM-proximal seed (Kim et al. 2016).
_CFD_POSITION_PENALTY_CAS12A: dict[int, float] = {
    1: 1.00, 2: 1.00, 3: 0.95, 4: 0.90, 5: 0.85,
    6: 0.80, 7: 0.70, 8: 0.60, 9: 0.50, 10: 0.40,
    11: 0.30, 12: 0.20, 13: 0.15, 14: 0.10, 15: 0.05,
    16: 0.03, 17: 0.02, 18: 0.02, 19: 0.01, 20: 0.01,
    21: 0.01, 22: 0.005, 23: 0.005,
}


def _compute_cfd_score(mismatch_positions: list[int], nuclease: str) -> float:
    """
    Cutting Frequency Determination (CFD) score for an off-target.
    Returns predicted cutting frequency at this off-target as a fraction of
    the on-target cutting frequency. Range [0, 1]; 1.0 = perfect match.

    mismatch_positions: list of 1-indexed positions (1 = PAM-proximal end,
                        matching the convention used elsewhere in this file).
    """
    if not mismatch_positions:
        return 1.0
    table = _CFD_POSITION_PENALTY_CAS9 if nuclease == "cas9" else _CFD_POSITION_PENALTY_CAS12A
    score = 1.0
    for pos in mismatch_positions:
        # Convert from PAM-proximal-1-indexed to the table's PAM-distal-1-indexed.
        # In this codebase position 1 = PAM-proximal (most dangerous).
        # The CFD table positions are also numbered with 1 = PAM-distal in the
        # original paper, so we flip: table_pos = guide_len - pos + 1.
        guide_len = _GUIDE_LENGTHS[nuclease]
        table_pos = guide_len - pos + 1
        penalty = table.get(table_pos, 0.5)
        score *= penalty
    return round(score, 4)


def _reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[b] for b in reversed(seq))


def _check_pam(seq: str, window_start: int, guide_len: int, nuclease: str) -> bool:
    """Return True if a valid PAM is present adjacent to the window."""
    if nuclease == "cas9":
        # NGG immediately after the window (Jinek et al. 2012)
        pam_start = window_start + guide_len
        if pam_start + 2 >= len(seq):
            return False
        return seq[pam_start + 1] == "G" and seq[pam_start + 2] == "G"
    else:  # cas12a
        # TTTV (V = A/C/G) immediately before the window (Zetsche et al. 2015)
        if window_start < 4:
            return False
        pam = seq[window_start - 4 : window_start]
        return len(pam) == 4 and pam[:3] == "TTT" and pam[3] in "ACG"


class PredictOfftargets:
    """
    Description:
        Scans a reference DNA sequence for potential CRISPR off-target sites.
        Supports both SpCas9 (NGG PAM, 20bp guide) and LbCas12a (TTTV PAM, 23bp guide).

        Takes the protospacer from the gRNA design step and slides it across
        every position of the reference on both strands, counting mismatches
        at each position. Any site within max_mismatches is flagged.

        PAM check:
          - Cas9:   NGG immediately 3' of the protospacer (Jinek et al. 2012).
          - Cas12a: TTTV (V = A/C/G) immediately 5' of the protospacer
                    (Zetsche et al. 2015, Cell 163(3):759-771).
          Without a valid PAM the nuclease cannot cut, so PAM-less sites are
          reported but scored LOW risk regardless of mismatch count.

        Risk scoring uses seed-region logic:
          - Cas9   seed = positions 1-12 from PAM-proximal end (Hsu et al. 2013).
          - Cas12a seed = positions 1-10 from PAM-proximal end; Cas12a is more
            mismatch-intolerant throughout the guide (Kim et al. 2016, Nat Biotechnol).
          HIGH: 0 mismatches, OR (0 seed mismatches + valid PAM).
          MEDIUM: ≤1 seed mismatch + PAM, OR ≤2 total mismatches + PAM.
          LOW: everything else.

        Supports circular references via is_circular=True so sites spanning
        the origin of a plasmid are not missed.

    Input:
        protospacer (str): the DNA protospacer WITHOUT the PAM.
                           20bp for Cas9, 23bp for Cas12a.
                           Only A, T, G, C — no ambiguous bases.
        reference (str):   the DNA sequence to scan. can be a resource name,
                           raw string, FASTA, or GenBank.
        nuclease (str):    "cas9" (default) or "cas12a".
        max_mismatches (int): flag sites with ≤N mismatches. default 3.
        is_circular (bool): wrap the reference before scanning. default False.

    Output:
        dict with keys:
            - protospacer: the input protospacer (echoed back)
            - nuclease: which nuclease was used
            - reference_length: length of the reference
            - strands_scanned: always 2
            - sites_evaluated: total windows checked across both strands
            - offtarget_sites: list of dicts, one per flagged site, each with:
                - position: where in the reference (0-indexed, forward strand coords)
                - strand: "+" or "-"
                - sequence: the window at that position
                - mismatches: how many bases differ from the protospacer
                - mismatch_positions: which positions (1 = PAM-proximal end)
                - seed_mismatches: mismatches in the seed region
                - has_pam: True if a valid PAM flanks this site
                - risk: "HIGH", "MEDIUM", or "LOW"
            - high_risk_count: number of HIGH risk sites
            - specificity_summary: one-sentence summary

    Tests:
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="ATGATGATGATGATGATGATAGG", nuclease="cas9"
            Expected Output: offtarget_sites has at least one entry with mismatches=0, has_pam=True, risk="HIGH"
            Description: exact Cas9 match + NGG PAM on forward strand → HIGH risk.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="CCCCCCCCCCCCCCCCCCCCAGG", nuclease="cas9"
            Expected Output: offtarget_sites is empty or all entries have mismatches > 0
            Description: no match → empty list or only high-mismatch sites.
        - Case:
            Input: protospacer="", reference="ATGATGATG"
            Expected Exception: ValueError
            Description: empty protospacer raises ValueError.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference=""
            Expected Exception: ValueError
            Description: empty reference raises ValueError.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="ATGATGATGATGATGATGATAGG", max_mismatches=0, nuclease="cas9"
            Expected Output: only the exact-match site is returned
            Description: max_mismatches=0 means only perfect matches come through.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGAT", reference="CCCCCCCCCCCCCCCCATGA", max_mismatches=3, nuclease="cas9"
            Expected Output: site at position 19 has has_pam=False
            Description: protospacer at last valid window with no room for PAM.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGATGAT", reference="TTTAATGATGATGATGATGATGATGATAAAAAAAAAA", nuclease="cas12a"
            Expected Output: offtarget_sites has at least one entry with mismatches=0, has_pam=True, risk="HIGH"
            Description: exact Cas12a match + TTTV PAM before 23bp window → HIGH risk.
        - Case:
            Input: protospacer="ATGATGATGATGATGATGATGAT", reference="ATGATGATGATGATGATGATGATAAAAAAAAAA", nuclease="cas12a"
            Expected Output: all sites have has_pam=False
            Description: Cas12a protospacer present but no TTTV PAM before it → no valid cut site.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        protospacer: str,
        reference: str,
        nuclease: str = "cas9",
        max_mismatches: int = 3,
        is_circular: bool = False,
    ) -> dict:

        protospacer = protospacer.upper().strip()
        reference = reference.upper().strip()
        nuclease = nuclease.lower().strip()

        # Normalize IUPAC ambiguity codes in the reference to N so scanning
        # still works on sequences from GenBank or genomic databases.
        _IUPAC_AMBIGUOUS = set("RYWSKMBDHV")
        reference = "".join(
            "N" if b in _IUPAC_AMBIGUOUS else b for b in reference
        )

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

        invalid_ref = [b for b in set(reference) if b not in "ATGCN"]
        if invalid_ref:
            raise ValueError(
                f"Invalid base(s) in reference: {sorted(invalid_ref)}. "
                "Only standard bases A, T, G, C (and N) are accepted."
            )

        guide_len = _GUIDE_LENGTHS[nuclease]
        seed_len = _SEED_LENGTHS[nuclease]
        ref_len = len(reference)

        # wrap enough bases to catch sites spanning the circular origin
        pam_len = 3 if nuclease == "cas9" else 4
        wrap_len = guide_len + pam_len
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

                # positions numbered 1 = PAM-proximal end (most dangerous)
                mismatch_positions = []
                for j in range(guide_len):
                    pam_proximal_idx = guide_len - 1 - j
                    if protospacer[pam_proximal_idx] != window[pam_proximal_idx]:
                        mismatch_positions.append(j + 1)

                seed_mismatches = sum(1 for p in mismatch_positions if p <= seed_len)

                has_pam = _check_pam(seq, i, guide_len, nuclease)

                # Without a valid PAM the nuclease cannot bind and cut,
                # so PAM-less sites are always LOW risk regardless of mismatch count.
                if not has_pam:
                    risk = "LOW"
                elif mismatches == 0:
                    risk = "HIGH"
                elif seed_mismatches == 0:
                    risk = "HIGH"
                elif seed_mismatches <= 1:
                    risk = "MEDIUM"
                elif mismatches <= 2:
                    risk = "MEDIUM"
                else:
                    risk = "LOW"

                # convert back to forward-strand coordinates in the original reference
                if strand == "+":
                    reported_position = i % ref_len
                else:
                    reported_position = (ref_len - i - guide_len) % ref_len

                # CFD score: predicted cutting frequency at this off-target
                # relative to the on-target site (Doench 2016 Rule Set 2).
                # Sites without a valid PAM cannot be cut, so CFD = 0 there.
                if has_pam:
                    cfd_score = _compute_cfd_score(mismatch_positions, nuclease)
                else:
                    cfd_score = 0.0

                offtarget_sites.append({
                    "position": reported_position,
                    "strand": strand,
                    "sequence": window,
                    "mismatches": mismatches,
                    "mismatch_positions": mismatch_positions,
                    "seed_mismatches": seed_mismatches,
                    "has_pam": has_pam,
                    "risk": risk,
                    "cfd_score": cfd_score,
                })

        offtarget_sites.sort(key=lambda s: (s["mismatches"], -s["seed_mismatches"]))

        high_risk_count = sum(1 for s in offtarget_sites if s["risk"] == "HIGH")

        if not offtarget_sites:
            specificity_summary = (
                f"No off-target sites found within {max_mismatches} mismatches "
                f"({nuclease.upper()}), the guide appears highly specific for this reference."
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

        # Aggregate CFD: sum of all off-target CFD scores, EXCLUDING the
        # on-target. Mirrors the "specificity score" definition in CRISPOR
        # (sum of off-target CFDs; lower is better). Skips the first 0-mismatch
        # HIGH-risk site as the on-target.
        on_seen = False
        cfd_off_sum = 0.0
        max_off_cfd = 0.0
        for s in offtarget_sites:
            if not on_seen and s["mismatches"] == 0 and s["risk"] == "HIGH":
                on_seen = True
                continue
            cfd_off_sum += s["cfd_score"]
            if s["cfd_score"] > max_off_cfd:
                max_off_cfd = s["cfd_score"]

        # Citations: seed-region scoring (Hsu/Kim) and CFD weights (Doench).
        # Cas9 → Hsu 2013 + Jinek (PAM) + Doench (CFD).
        # Cas12a → Zetsche + Kim (seed) + Doench (CFD adapted).
        if nuclease == "cas9":
            citation_keys = ["jinek_2012", "hsu_2013", "doench_2016"]
        else:
            citation_keys = ["zetsche_2015", "kim_2016_cas12a", "doench_2016"]

        return {
            "protospacer": protospacer,
            "nuclease": nuclease,
            "reference_length": ref_len,
            "strands_scanned": 2,
            "sites_evaluated": sites_evaluated,
            "offtarget_sites": offtarget_sites,
            "high_risk_count": high_risk_count,
            "aggregate_offtarget_cfd": round(cfd_off_sum, 4),
            "max_offtarget_cfd": round(max_off_cfd, 4),
            "specificity_summary": specificity_summary,
            "citations": format_citations(cites(*citation_keys)),
        }


_instance = PredictOfftargets()
_instance.initiate()
predict_offtargets = _instance.run
