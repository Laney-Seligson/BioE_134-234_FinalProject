class DesignCas9Grna:
    """
    Description:
        Designs a Cas9 gRNA for a given DNA target sequence.

        Finds the first NGG PAM, takes the 20 nt immediately upstream as the
        protospacer, appends the tracrRNA scaffold, and converts T -> U to
        produce the gRNA.

    Input:
        seq (str): DNA target sequence (A, T, G, C only).

    Output:
        dict: A dictionary with keys:
            - grna_sequence (str): Full RNA gRNA (protospacer + scaffold).
            - protospacer (str): The 20 nt DNA protospacer.
            - pam_site (str): The 3 bp NGG PAM immediately after the protospacer.

    Tests:
        - Case:
            Input: seq="CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
            Expected Output: protospacer is 20 bp
            Description: Standard test with an NGG PAM present.
        - Case:
            Input: seq=""
            Expected Exception: ValueError
            Description: Empty sequence raises ValueError.
        - Case:
            Input: seq="ATGATGATGATGATGATGATGATG"
            Expected Exception: ValueError
            Description: No NGG PAM site found raises ValueError.
    """

    _tracrRNA: str

    def initiate(self) -> None:
        self._tracrRNA = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

    def run(self, seq: str) -> dict:
        seq = seq.upper()

        if not seq:
            raise ValueError("Sequence must not be empty.")

        invalid = [b for b in set(seq) if b not in "ATGC"]
        if invalid:
            raise ValueError(f"Invalid base(s) in sequence: {sorted(invalid)}.")

        for i in range(20, len(seq) - 2):
            if seq[i+1] == "G" and seq[i+2] == "G":
                pam = seq[i:i+3]
                protospacer = seq[i-20:i]
                grna_rna = (protospacer + self._tracrRNA).replace("T", "U")
                return {
                    "grna_sequence": grna_rna,
                    "protospacer": protospacer,
                    "pam_site": pam,
                }

        raise ValueError("No NGG PAM site found in sequence.")


_instance = DesignCas9Grna()
_instance.initiate()
design_cas9_grna = _instance.run
