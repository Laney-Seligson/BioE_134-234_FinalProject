import re as _re


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

        for i in range(len(seq) - 26):
            if seq[i:i+3] == "TTT" and seq[i+3] in "ACG":
                protospacer = seq[i+4:i+27]
                gc = sum(1 for b in protospacer if b in "GC") / 23
                all_guides.append({
                    "crrna_sequence": (self._direct_repeat + protospacer).replace("T", "U"),
                    "protospacer": protospacer,
                    "pam_site": seq[i:i+4],
                    "_gc": gc,
                    "_pos": i,
                })

        if not all_guides:
            raise ValueError("No TTTV PAM site found in sequence.")

        # Match cas_selector's validity filter (GC 30-70%).
        # Sample evenly across the sequence so candidates span the full gene.
        preferred = [g for g in all_guides if 0.30 <= g["_gc"] <= 0.70]
        pool = preferred if len(preferred) >= 10 else all_guides

        if len(pool) <= 10:
            selected = pool
        else:
            bucket_size = len(pool) / 10
            selected = [pool[int(i * bucket_size)] for i in range(10)]

        return [{k: v for k, v in g.items() if not k.startswith("_")} for g in selected]


_instance = DesignCas12aCrrna()
_instance.initiate()
design_cas12a_crrna = _instance.run
