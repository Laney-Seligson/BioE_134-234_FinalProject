class DesignCas12aCrrna:
    """
    Description:
        Designs an LbCas12a crRNA for a given DNA target sequence.

        Scans a DNA target sequence for all TTTV PAM sites (TTTA, TTTC, or
        TTTG), takes the 23 nt immediately downstream of each as the
        protospacer, prepends the LbCas12a direct repeat, and converts T -> U
        to produce the crRNA. Returns up to 10 candidates in order of
        appearance.

    Input:
        seq (str): DNA target sequence (A, T, G, C only).

    Output:
        list: Up to 10 dictionaries, each with keys:
            - crrna_sequence (str): Full RNA crRNA (direct repeat + spacer).
            - protospacer (str): The 23 nt DNA protospacer.
            - pam_site (str): The 4 bp TTTV PAM immediately before the protospacer.

    Tests:
        - Case:
            Input: seq="CCCTAGATGCCTTTTAGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
            Expected Output: protospacer is 23 bp
            Description: Standard test with a TTTV PAM present.
        - Case:
            Input: seq=""
            Expected Exception: ValueError
            Description: Empty sequence raises ValueError.
        - Case:
            Input: seq="ATGATGATGATGATGATGATGATG"
            Expected Exception: ValueError
            Description: No TTTV PAM site found raises ValueError.
    """

    _direct_repeat: str

    def initiate(self) -> None:
        self._direct_repeat = "AATTTCTACTAAGTGTAGAT"

    def run(self, seq: str) -> list:
        seq = seq.upper()

        if not seq:
            raise ValueError("Sequence must not be empty.")

        invalid = [b for b in set(seq) if b not in "ATGC"]
        if invalid:
            raise ValueError(f"Invalid base(s) in sequence: {sorted(invalid)}.")

        results = []

        for i in range(len(seq) - 26):
            if seq[i:i+3] == "TTT" and seq[i+3] in "ACG":
                pam = seq[i:i+4]
                protospacer = seq[i+4:i+27]
                crrna_rna = (self._direct_repeat + protospacer).replace("T", "U")
                results.append({
                    "crrna_sequence": crrna_rna,
                    "protospacer": protospacer,
                    "pam_site": pam,
                })

        if not results:
            raise ValueError("No TTTV PAM site found in sequence.")

        return results[:10]


_instance = DesignCas12aCrrna()
_instance.initiate()
design_cas12a_crrna = _instance.run
