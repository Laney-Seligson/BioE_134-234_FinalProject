class CasSelector:
    """
    Description:
        Counts the available PAM sites in a DNA sequence and recommends
        whether to use Cas9 or Cas12a for CRISPR editing.

        Rather than using a fixed GC% threshold, this tool directly counts
        how many target sites each nuclease can reach in the sequence:

            Cas9 (SpCas9) requires an NGG PAM site (G-rich, 3' of protospacer).
            Cas12a (LbCas12a) requires a TTTV PAM site (T-rich, 5' of protospacer,
            where V = A, C, or G).

        The nuclease with more available PAM sites in the sequence is recommended.
        This approach is grounded in the literature: PAM site density directly
        determines how many genomic loci each system can target (Zetsche et al.
        2015, Nature Biotechnology; Addgene CRISPR Guide).

        Both strands are scanned since either strand can serve as the target.
        Ties default to Cas9 as the more broadly validated system.

        The framework resolves the input before calling run(): a GenBank file,
        FASTA string, resource name (e.g. "pBR322"), or a raw sequence string
        are all accepted and converted to a clean uppercase sequence automatically.

    Input:
        seq (str): DNA sequence. Accepts a resource name, a raw sequence
                   string, a FASTA-formatted string, or a GenBank-formatted
                   string. The framework resolves the format automatically.

    Output:
        dict: A dictionary with keys:
            - ngg_count (int): Number of NGG PAM sites found (both strands).
            - tttv_count (int): Number of TTTV PAM sites found (both strands).
            - recommendation (str): "Cas9" or "Cas12a".
            - rationale (str): One-sentence explanation of the recommendation.

    Tests:
        - Case:
            Input: seq="GCGCGCGCGCGG"
            Expected Output: {"recommendation": "Cas9"}
            Description: GC-rich sequence has more NGG PAM sites.
        - Case:
            Input: seq="ATTTAATTTAATTTC"
            Expected Output: {"recommendation": "Cas12a"}
            Description: AT-rich sequence has more TTTV PAM sites.
        - Case:
            Input: seq=""
            Expected Exception: ValueError
            Description: Empty sequence raises ValueError.
    """

    def initiate(self) -> None:
        pass

    def _reverse_complement(self, seq: str) -> str:
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join(complement.get(b, "N") for b in reversed(seq))

    def run(self, seq: str) -> dict:
        """Count NGG and TTTV PAM sites on both strands and recommend a Cas system."""
        seq = seq.upper()

        if not seq:
            raise ValueError("Sequence must not be empty.")

        if len(seq) < 4:
            raise ValueError("Sequence is too short to contain any PAM sites.")

        rc = self._reverse_complement(seq)

        # Count NGG PAM sites on both strands (Cas9)
        ngg_count = sum(
            1 for i in range(len(seq) - 2)
            if seq[i + 1] == "G" and seq[i + 2] == "G"
        ) + sum(
            1 for i in range(len(rc) - 2)
            if rc[i + 1] == "G" and rc[i + 2] == "G"
        )

        # Count TTTV PAM sites on both strands (LbCas12a, V = A/C/G)
        tttv_count = sum(
            1 for i in range(len(seq) - 3)
            if seq[i] == "T" and seq[i + 1] == "T" and seq[i + 2] == "T"
            and seq[i + 3] in "ACG"
        ) + sum(
            1 for i in range(len(rc) - 3)
            if rc[i] == "T" and rc[i + 1] == "T" and rc[i + 2] == "T"
            and rc[i + 3] in "ACG"
        )

        if ngg_count >= tttv_count:
            recommendation = "Cas9"
            rationale = (
                f"Found {ngg_count} NGG PAM sites vs {tttv_count} TTTV PAM sites — "
                "Cas9 has more targetable sites in this sequence."
            )
        else:
            recommendation = "Cas12a"
            rationale = (
                f"Found {tttv_count} TTTV PAM sites vs {ngg_count} NGG PAM sites — "
                "Cas12a (LbCas12a) has more targetable sites in this sequence."
            )

        return {
            "ngg_count": ngg_count,
            "tttv_count": tttv_count,
            "recommendation": recommendation,
            "rationale": rationale,
        }


_instance = CasSelector()
_instance.initiate()
cas_selector = _instance.run
