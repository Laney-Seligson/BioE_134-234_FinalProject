class DesignCas12aCrrna:
    """
    Description:
        Designs an optimal FnCas12a CRISPR RNA (crRNA) for a given DNA target sequence.

        Finds every valid TTTV PAM site in the sequence (where V = A, C, or G),
        extracts the 23 bp protospacer immediately DOWNSTREAM of each PAM, and
        scores each candidate using efficiency rules drawn from two papers:

            Zetsche et al. 2015 (Nature Biotechnology 33, 551-557.
            DOI: 10.1038/nbt.3347) — established the FnCas12a PAM (TTTV),
            protospacer orientation (downstream of PAM), and crRNA structure
            (direct repeat + 23 bp spacer).

            Pausch et al. 2020 (Nature Communications 11, 2708.
            DOI: 10.1038/s41467-020-16669-9) — identified sequence-level
            efficiency rules for Cas12a:

            1. GC content of the protospacer should be between 30% and 70%.
               Cas12a tolerates a broader GC range than Cas9 because it operates
               in AT-rich genomes, but extremes still reduce efficiency.

            2. No run of 4 or more consecutive T bases (poly-T). Such runs
               can interfere with crRNA transcription and processing.

            3. The first base of the protospacer (immediately after the PAM)
               should be C or G. Cleavage kinetics are improved when the
               PAM-proximal base is a strong (C/G) base pair.

        Each rule contributes 1 point to an efficiency score (0-3). The
        candidate with the highest score is selected. If multiple candidates
        tie, the first occurrence is returned. The selected protospacer is
        prepended with the FnCas12a direct repeat and converted from DNA to
        RNA (T -> U).

        Key difference from Cas9: the protospacer lies DOWNSTREAM of the PAM
        (not upstream), and there is no tracrRNA — Cas12a processes its own
        single-guide crRNA from a direct repeat scaffold.

        The framework resolves the input before calling run(): a GenBank file,
        FASTA string, resource name (e.g. "pBR322"), or a raw sequence string
        are all accepted and converted to a clean uppercase sequence automatically.

    Input:
        seq (str): DNA target sequence. Accepts a resource name, a raw sequence
                   string, a FASTA-formatted string, or a GenBank-formatted
                   string. Only standard bases A, T, G, C are accepted.

    Output:
        dict: A dictionary with keys:
            - crrna_sequence (str): Full RNA crRNA sequence (direct repeat + spacer).
            - protospacer (str): The selected 23 bp DNA protospacer sequence.
            - pam_site (str): The 4 bp TTTV PAM sequence immediately preceding the protospacer.
            - gc_fraction (float): GC fraction of the protospacer (0.0 to 1.0).
            - efficiency_score (int): Rule-based score from 0 (worst) to 3 (best).
            - warnings (list): List of rule violations for the selected guide, if any.
            - candidates_evaluated (int): Total number of TTTV PAM sites evaluated.

    Tests:
        - Case:
            Input: seq="CCCTAGATGCCTTTTAGCTCAGAAACCTGCCAGTTTGCTGGCACGTTTTTTTCTTTTGTCTTTAGTTCTCACGTTTGTCATACTTGACAACGCTTCTTTAACCAAATATAATTGTTC"
            Expected Output: protospacer is 23 bp, efficiency_score >= 0
            Description: Standard test case from BioE134 homework adapted for Cas12a.
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
    _valid_v_bases: str

    def initiate(self) -> None:
        # FnCas12a direct repeat sequence — from Zetsche et al. 2015
        self._direct_repeat = "AATTTCTACTGTTGTAGAT"
        # V in TTTV = any base except T
        self._valid_v_bases = "ACG"

    def run(self, seq: str) -> dict:
        """Find the best FnCas12a protospacer in seq and return the full crRNA."""
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

        # Find all TTTV PAM sites and score each candidate protospacer.
        # PAM = 4 bp (TTTV), protospacer = 23 bp DOWNSTREAM of PAM.
        # Searching for all sites (not just the first) is original to this tool —
        # the BioE134 homework only finds the first PAM.
        candidates = []
        for i in range(len(seq) - 26):  # need 4 bp PAM + 23 bp protospacer
            if (seq[i]     == "T" and
                    seq[i + 1] == "T" and
                    seq[i + 2] == "T" and
                    seq[i + 3] in self._valid_v_bases):  # from BioE134 homework: find_PAM logic adapted for TTTV

                pam = seq[i:i + 4]
                protospacer = seq[i + 4:i + 27]  # from BioE134 homework: find_protospacer logic, downstream for Cas12a

                # --- Zetsche 2015 / Pausch 2020 scoring rules (original to this tool) ---
                gc = sum(1 for b in protospacer if b in "GC")
                gc_fraction = gc / 23

                gc_ok      = 0.3 <= gc_fraction <= 0.7  # rule 1: GC content 30-70% (Pausch 2020)
                no_poly_t  = "TTTT" not in protospacer   # rule 2: no poly-T run (Pausch 2020)
                good_start = protospacer[0] in "CG"      # rule 3: first base = C or G (Pausch 2020)

                score = int(gc_ok) + int(no_poly_t) + int(good_start)

                candidates.append({
                    "protospacer": protospacer,
                    "pam": pam,
                    "gc_fraction": round(gc_fraction, 4),
                    "gc_ok": gc_ok,
                    "no_poly_t": no_poly_t,
                    "good_start": good_start,
                    "score": score,
                })

        if not candidates:
            raise ValueError("No valid TTTV PAM site found in sequence.")

        # Select the highest scoring candidate — original to this tool
        best = max(candidates, key=lambda c: c["score"])

        # Build warnings for any failed rules
        warnings = []
        if not best["gc_ok"]:
            warnings.append(
                f"Protospacer GC content is {best['gc_fraction']:.1%} — "
                "optimal range is 30-70% for Cas12a (Pausch et al. 2020)."
            )
        if not best["no_poly_t"]:
            warnings.append(
                "Protospacer contains a poly-T run (TTTT+) which may interfere "
                "with crRNA transcription or processing (Pausch et al. 2020)."
            )
        if not best["good_start"]:
            warnings.append(
                "First protospacer base is not C or G — cleavage kinetics may "
                "be reduced at the PAM-proximal position (Pausch et al. 2020)."
            )

        # Build full crRNA: direct repeat + protospacer spacer — from BioE134 homework: design_FnCas12a_crRNA logic
        crrna_dna = self._direct_repeat + best["protospacer"]

        # Convert T -> U to produce RNA — from BioE134 homework: design_cas9_gRNA T->U logic
        crrna_rna = crrna_dna.replace("T", "U")

        return {
            "crrna_sequence": crrna_rna,
            "protospacer": best["protospacer"],
            "pam_site": best["pam"],
            "gc_fraction": best["gc_fraction"],
            "efficiency_score": best["score"],
            "warnings": warnings,
            "candidates_evaluated": len(candidates),
        }


_instance = DesignCas12aCrrna()
_instance.initiate()
design_cas12a_crrna = _instance.run
