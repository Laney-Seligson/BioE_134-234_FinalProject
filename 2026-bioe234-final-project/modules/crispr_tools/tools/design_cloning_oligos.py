class DesignCloningOligos:
    """
    Description:
        Designs the top and bottom strand oligos needed to clone a CRISPR
        guide RNA protospacer into a linearized expression vector.

        After a protospacer is selected (e.g. by design_cas9_grna or
        design_cas12a_crrna), it must be inserted into a vector as a short
        double-stranded DNA. This is done by ordering two complementary oligos,
        annealing them, and ligating them into the restriction-digested vector.

        The core operation is vector-agnostic:
            top oligo    = top_overhang + protospacer
            bottom oligo = bottom_overhang + reverse_complement(protospacer)

        The overhangs are determined by the restriction enzyme used to
        linearize the vector, not by the Cas system. Default overhangs
        (CACC / AAAC) match pX330 and similar BbsI-digested Cas9 vectors
        (Ran et al. 2013, Nature Protocols 8, 2281-2308.
        DOI: 10.1038/nprot.2013.143). For other vectors, supply the
        appropriate overhangs as arguments.

        If the vector uses a U6 promoter (common in Cas9 vectors), RNA
        Pol III transcription is most efficient when the first transcribed
        base is G. If the protospacer does not start with G, this tool
        automatically prepends one and notes it in the output so the user
        knows the synthesized oligo differs slightly from the protospacer.

    Input:
        protospacer (str): The DNA protospacer sequence (20 bp for Cas9,
                           23 bp for Cas12a, or any length). Raw DNA string
                           only — no resource names or file formats.
        top_overhang (str): 5' overhang added to the top strand oligo.
                            Default is "CACC" (pX330 / BbsI).
        bottom_overhang (str): 5' overhang added to the bottom strand oligo.
                               Default is "AAAC" (pX330 / BbsI).

    Output:
        dict: A dictionary with keys:
            - top_oligo (str): Full top strand oligo sequence (5' to 3').
            - bottom_oligo (str): Full bottom strand oligo sequence (5' to 3').
            - g_prepended (bool): True if a G was added to the protospacer
                                  for U6 promoter compatibility.
            - notes (str): Plain-language description of what was done.

    Tests:
        - Case:
            Input: protospacer="TCAGAAACCTGCCAGTTTGC", top_overhang="CACC", bottom_overhang="AAAC"
            Expected Output: {"top_oligo": "CACCTCAGAAACCTGCCAGTTTGC", "g_prepended": False}
            Description: Protospacer starts with T, no G prepended. Standard Cas9 overhangs.
        - Case:
            Input: protospacer="GCAGAAACCTGCCAGTTTGC", top_overhang="CACC", bottom_overhang="AAAC"
            Expected Output: {"top_oligo": "CACCGCAGAAACCTGCCAGTTTGC", "g_prepended": False}
            Description: Protospacer starts with G, no G prepended needed.
        - Case:
            Input: protospacer="ACAGAAACCTGCCAGTTTGC", top_overhang="CACC", bottom_overhang="AAAC"
            Expected Output: {"top_oligo": "CACCGACAGAAACCTGCCAGTTTGC", "g_prepended": True}
            Description: Protospacer starts with A, G prepended for U6 compatibility.
        - Case:
            Input: protospacer=""
            Expected Exception: ValueError
            Description: Empty protospacer raises ValueError.
        - Case:
            Input: protospacer="ATGX"
            Expected Exception: ValueError
            Description: Non-DNA character raises ValueError.
    """

    _complement: dict

    def initiate(self) -> None:
        # Complement table for reverse complement — from BioE134 homework
        self._complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def _reverse_complement(self, seq: str) -> str:
        """Return the reverse complement of a DNA sequence."""
        return "".join(self._complement[b] for b in reversed(seq))

    def run(
        self,
        protospacer: str,
        top_overhang: str = "CACC",
        bottom_overhang: str = "AAAC",
    ) -> dict:
        """Design top and bottom cloning oligos for a CRISPR protospacer."""
        protospacer = protospacer.upper()
        top_overhang = top_overhang.upper()
        bottom_overhang = bottom_overhang.upper()

        if not protospacer:
            raise ValueError("Protospacer must not be empty.")

        invalid = [b for b in set(protospacer) if b not in "ATGC"]
        if invalid:
            raise ValueError(
                f"Invalid base(s) in protospacer: {sorted(invalid)}. "
                "Only standard bases A, T, G, C are accepted."
            )

        # Prepend G if protospacer does not start with G, for U6 promoter
        # compatibility (Ran et al. 2013). Only applied when using default
        # overhangs — if the user supplies custom overhangs they control this.
        g_prepended = False
        if top_overhang == "CACC" and protospacer[0] != "G":
            protospacer = "G" + protospacer
            g_prepended = True

        rc = self._reverse_complement(protospacer)

        top_oligo    = top_overhang + protospacer
        bottom_oligo = bottom_overhang + rc

        if g_prepended:
            notes = (
                "A G was prepended to the protospacer for U6 promoter "
                "compatibility (Ran et al. 2013). The synthesized oligo is "
                "one base longer than the original protospacer."
            )
        else:
            notes = (
                "Protospacer starts with G — no modification needed for "
                "U6 promoter compatibility."
            )

        return {
            "top_oligo": top_oligo,
            "bottom_oligo": bottom_oligo,
            "g_prepended": g_prepended,
            "notes": notes,
        }


_instance = DesignCloningOligos()
_instance.initiate()
design_cloning_oligos = _instance.run
