import re as _re


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
        if _re.search(r'\.\.\.\s*\([\d,]+\s*chars?\)', seq):
            raise ValueError(
                "seq is a truncated display string, not a DNA sequence. "
                "Pass the full sequence from sequence_info['sequence'] as returned "
                "by crispr_fetch_target_sequence."
            )

        seq = seq.upper()

        if not seq:
            raise ValueError("Sequence must not be empty.")

        invalid = [b for b in set(seq) if b not in "ATGC"]
        if invalid:
            raise ValueError(f"Invalid base(s) in sequence: {sorted(invalid)}.")

        all_guides = []
        for i in range(20, len(seq) - 2):
            if seq[i+1:i+3] == "GG":
                protospacer = seq[i-20:i]
                gc = sum(1 for b in protospacer if b in "GC") / 20
                all_guides.append({
                    "grna_sequence": (protospacer + self._sgRNA_scaffold).replace("T", "U"),
                    "protospacer": protospacer,
                    "pam_site": seq[i:i+3],
                    "_gc": gc,
                    "_pos": i,
                })

        if not all_guides:
            raise ValueError("No NGG PAM site found in sequence.")

        # Match cas_selector's validity filter (GC 30-70%, consistent range).
        # Sample evenly across the sequence so the 10 returned guides cover the
        # full gene rather than clustering at the first NGG-dense region.
        preferred = [g for g in all_guides if 0.30 <= g["_gc"] <= 0.70]
        pool = preferred if len(preferred) >= 10 else all_guides

        if len(pool) <= 10:
            selected = pool
        else:
            bucket_size = len(pool) / 10
            selected = [pool[int(i * bucket_size)] for i in range(10)]

        return [{k: v for k, v in g.items() if not k.startswith("_")} for g in selected]


_instance = DesignCas9Grna()
_instance.initiate()
design_cas9_grna = _instance.run
