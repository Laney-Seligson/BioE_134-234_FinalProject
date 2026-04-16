class DesignCas9Grna:
    """
    Description:
        Designs an optimal Cas9 guide RNA (gRNA) for a given DNA target sequence.

        Finds every valid NGG PAM site in the sequence, extracts the 20 bp
        protospacer upstream of each, and scores each candidate using three
        empirically validated rules from Doench et al. 2016 (Nature Biotechnology
        34, 184-191. DOI: 10.1038/nbt.3437):

            1. GC content of the protospacer should be between 40% and 70%.
               Outside this range, on-target efficiency drops significantly.

            2. No run of 4 or more consecutive T bases (poly-T). Such runs
               terminate RNA Polymerase III transcription prematurely, producing
               a truncated, non-functional gRNA.

            3. The final base of the protospacer (position 20, immediately
               adjacent to the PAM) should be G. Doench et al. identified
               a strong positional preference for G at this position.

        Each rule contributes 1 point to an efficiency score (0-3). The
        candidate with the highest score is selected. If multiple candidates
        tie, the first occurrence is returned. The selected protospacer is
        then appended to the standard tracrRNA scaffold and converted from
        DNA to RNA (T -> U).

        The framework resolves the input before calling run(): a GenBank file,
        FASTA string, resource name (e.g. "pBR322"), or a raw sequence string
        are all accepted and converted to a clean uppercase sequence automatically.

    Input:
        seq (str): DNA target sequence. Accepts a resource name, a raw sequence
                   string, a FASTA-formatted string, or a GenBank-formatted
                   string. Only standard bases A, T, G, C are accepted.

    Output:
        dict: A dictionary with keys:
            - grna_sequence (str): Full RNA gRNA sequence (protospacer + scaffold).
            - protospacer (str): The selected 20 bp DNA protospacer sequence.
            - pam_site (str): The 3 bp PAM sequence (NGG) immediately following the protospacer.
            - gc_fraction (float): GC fraction of the protospacer (0.0 to 1.0).
            - efficiency_score (int): Doench rule score from 0 (worst) to 3 (best).
            - warnings (list): List of rule violations for the selected guide, if any.
            - candidates_evaluated (int): Total number of PAM sites evaluated.

    Tests:
        - Case:
            Input: seq="CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
            Expected Output: protospacer="TCAGAAACCTGCCAGTTTGC", efficiency_score=3
            Description: Standard test case from BioE134 homework. Best guide scores 3/3.
        - Case:
            Input: seq="ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGG"
            Expected Output: protospacer is 20 bp, efficiency_score >= 0
            Description: Minimal sequence with one PAM site at the end.
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
        # tracrRNA scaffold sequence — from BioE134 CRISPR homework
        self._tracrRNA = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"

    def run(self, seq: str) -> dict:
        """Find the best Cas9 protospacer in seq and return the full gRNA."""
        seq = seq.upper()  # from BioE134 homework: uppercase input before processing

        if not seq:
            raise ValueError("Sequence must not be empty.")

        # Validation — from BioE134 homework: validate_DNA logic
        invalid = [b for b in set(seq) if b not in "ATGC"]
        if invalid:
            raise ValueError(
                f"Invalid base(s) in sequence: {sorted(invalid)}. "
                "Only standard bases A, T, G, C are accepted."
            )

        # Find all NGG PAM sites and score each candidate protospacer.
        # PAM must start at index >= 20 so the 20 bp protospacer fits upstream.
        # Searching for all sites (not just the first) is original to this tool —
        # the BioE134 homework only finds the first PAM.
        candidates = []
        for i in range(20, len(seq) - 2):  # from BioE134 homework: find_PAM logic
            if seq[i + 1] == "G" and seq[i + 2] == "G":
                pam = seq[i:i + 3]
                protospacer = seq[i - 20:i]  # from BioE134 homework: find_protospacer logic

                # --- Doench et al. 2016 scoring rules (original to this tool) ---
                gc = sum(1 for b in protospacer if b in "GC")
                gc_fraction = gc / 20

                gc_ok     = 0.4 <= gc_fraction <= 0.7   # rule 1: GC content 40-70%
                no_poly_t = "TTTT" not in protospacer    # rule 2: no poly-T run
                good_end  = protospacer[-1] == "G"       # rule 3: position 20 = G

                score = int(gc_ok) + int(no_poly_t) + int(good_end)

                candidates.append({
                    "protospacer": protospacer,
                    "pam": pam,
                    "gc_fraction": round(gc_fraction, 4),
                    "gc_ok": gc_ok,
                    "no_poly_t": no_poly_t,
                    "good_end": good_end,
                    "score": score,
                })

        if not candidates:
            raise ValueError("No valid NGG PAM site found in sequence.")

        # Select the highest scoring candidate — original to this tool
        best = max(candidates, key=lambda c: c["score"])

        # Build warnings for any failed rules
        warnings = []
        if not best["gc_ok"]:
            warnings.append(
                f"Protospacer GC content is {best['gc_fraction']:.1%} — "
                "optimal range is 40-70% (Doench et al. 2016)."
            )
        if not best["no_poly_t"]:
            warnings.append(
                "Protospacer contains a poly-T run (TTTT+) which may cause "
                "premature RNA Pol III termination (Doench et al. 2016)."
            )
        if not best["good_end"]:
            warnings.append(
                "Protospacer position 20 is not G — efficiency may be reduced "
                "(Doench et al. 2016)."
            )

        # Combine protospacer + scaffold — from BioE134 homework: design_cas9_DNA logic
        grna_dna = best["protospacer"] + self._tracrRNA

        # Convert T -> U to produce RNA — from BioE134 homework: design_cas9_gRNA logic
        grna_rna = grna_dna.replace("T", "U")

        return {
            "grna_sequence": grna_rna,
            "protospacer": best["protospacer"],
            "pam_site": best["pam"],
            "gc_fraction": best["gc_fraction"],
            "efficiency_score": best["score"],
            "warnings": warnings,
            "candidates_evaluated": len(candidates),
        }


_instance = DesignCas9Grna()
_instance.initiate()
design_cas9_grna = _instance.run
