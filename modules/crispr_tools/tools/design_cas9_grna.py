class DesignCas9Grna:
    """
    Description:
        Designs a Cas9 gRNA for a given DNA target sequence.

        Scans a DNA target sequence for all NGG PAM sites, takes the 20 nt
        immediately upstream of each as the protospacer, appends the tracrRNA
        scaffold, and converts T -> U to produce the gRNA. Returns up to 10
        candidates in order of appearance.

    Citations: 
        - Claim:SpCas9 uses a NGG PAM - Jinek et al., 2012, Science 337(6096):816-821 (original Cas9 gRNA design).
        - Claim:Guide Targets the 20nt sequence upstream of NGG - Cong et al., 2013, Science
        - Claim:TracrRNA scaffold sequence - Jinek et al., 2012, Science 337(6096):816-821 (original Cas9 gRNA design).
    Input:
        seq (str): DNA target sequence (A, T, G, C only).

    Output:
        list: Up to 10 dictionaries, each with keys:
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

    _sgRNA_scaffold: str

    def initiate(self) -> None:
        self._sgRNA_scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

    def run(self, seq: str) -> list:
        seq = seq.upper()

        if not seq:
            raise ValueError("Sequence must not be empty.")

        invalid = [b for b in set(seq) if b not in "ATGC"]
        if invalid:
            raise ValueError(f"Invalid base(s) in sequence: {sorted(invalid)}.")

        results = []

        for i in range(20, len(seq) - 2):
            if seq[i+1:i+3] == "GG":
                pam = seq[i:i+3]
                protospacer = seq[i-20:i]
                grna_rna = (protospacer + self._sgRNA_scaffold).replace("T", "U")
                results.append({
                    "grna_sequence": grna_rna,
                    "protospacer": protospacer,
                    "pam_site": pam,
                })

        if not results:
            raise ValueError("No NGG PAM site found in sequence.")

        return results[:10]


_instance = DesignCas9Grna()
_instance.initiate()
design_cas9_grna = _instance.run
