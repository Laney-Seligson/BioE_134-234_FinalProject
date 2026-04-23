class CasSelector:
    """
    Description:
        Recommends Cas9 (SpCas9) or Cas12a (LbCas12a) for a CRISPR experiment.

        Decision logic (three steps):

        Step 1 -- Valid guide counts.
              For each PAM site on both strands, the adjacent spacer is
              evaluated for quality:
                - PAM: NGG (3') for SpCas9; TTTV (5', V = A/C/G) for LbCas12a.
                - Spacer quality: GC content 30-70%; no TTTT run (Pol III
                  termination); no homopolymer run >= 5 bases.
              Spacers that fail any filter are excluded.

        Step 2 -- Is the difference meaningful?
              If one nuclease has zero valid guides and the other does not,
              the feasible nuclease wins outright. If both have valid guides,
              the nuclease with more wins only when its count is at least
              margin_threshold times higher (default 1.5x). Below that ratio
              the counts are treated as a tie and Step 3 decides.
              (Zetsche et al., 2015, Cell 163(3):759-771)

        Step 3 -- Secondary criteria (used when guide counts are tied or
              both zero):
              a. Multiplexing (num_targets >= 2) -> Cas12a.
                 Cas12a processes a crRNA array from a single transcript,
                 avoiding a separate transcription unit per guide.
                 (Zetsche et al., 2017, Nat Biotechnol 35(1):31-34)
              b. High specificity required -> Cas12a.
                 Cas12a is more mismatch-intolerant than SpCas9, with
                 off-target burden comparable to engineered high-fidelity
                 Cas9 variants.
                 (Kleinstiver et al., 2016, Nat Biotechnol 34:869-874)
              c. Otherwise -> Cas9 (well-established default).

        AT-richness note: when gc_content < 0.45, AT-richness is reported
        in the rationale as supporting context for a Cas12a recommendation
        but does not override the guide-count or secondary-criteria signals.
        (IDT Alt-R Cas12a documentation; Zetsche et al., 2015)

        Repair strategy is noted in the rationale but does not change the
        recommendation -- both nucleases support HDR and NHEJ.
        (Kim et al., 2016, Nat Biotechnol 34:863-868)

        The framework resolves the input before calling run(): a GenBank file,
        FASTA string, resource name (e.g. "pBR322"), or a raw sequence string
        are all accepted and converted to a clean uppercase sequence automatically.

    Input:
        seq (str): DNA sequence to analyze. Accepts a resource name, a raw
                   sequence string, a FASTA-formatted string, or a GenBank-
                   formatted string. The framework resolves the format
                   automatically.
        repair_template (bool): True if a homology-directed repair (HDR) donor
                   template will be provided; False for NHEJ-based knockout.
                   Informs the rationale but does not change the recommendation.
        num_targets (int): Number of distinct genomic loci to edit. Values >= 2
                   indicate multiplexing and favor Cas12a in Step 3. Defaults to 1.
        high_specificity (bool): True if minimizing off-target edits is a
                   priority (e.g., therapeutic context). Favors Cas12a in
                   Step 3. Defaults to False.
        system (str or None): If provided, must be "Cas9" or "Cas12a". Bypasses
                   all analysis and returns the specified system directly.
                   Defaults to None (auto-select).
        cas12a_spacer_len (int): Length of the Cas12a spacer to evaluate for
                   quality filtering. Must be 20-24. Defaults to 23 (LbCas12a
                   standard).
        margin_threshold (float): Minimum ratio by which one nuclease must
                   exceed the other's valid guide count for the count difference
                   to drive the recommendation (Step 2). Must be >= 1.0.
                   Default 1.5 means one nuclease needs at least 50% more valid
                   guides (e.g. 40 vs 36 is a tie; 40 vs 20 is a clear win).

    Output:
        dict: A dictionary with keys:
            - ngg_count (int): Number of NGG PAM sites found (both strands).
            - tttv_count (int): Number of TTTV PAM sites found (both strands).
            - cas9_valid_guides (int): NGG sites with spacers passing quality filters.
            - cas12a_valid_guides (int): TTTV sites with spacers passing quality filters.
            - gc_content (float): GC fraction of the sequence (0.0-1.0).
            - recommendation (str): "Cas9" or "Cas12a".
            - rationale (str): Explanation of the recommendation and any
                               relevant secondary considerations.

    Tests:
        - Case:
            Input: seq="GCGCGCGCGCGG", repair_template=False, num_targets=1, high_specificity=False
            Expected Output: {"recommendation": "Cas9"}
            Description: GC-rich sequence has many more NGG PAMs than TTTV PAMs; Cas9 wins on guide count.
        - Case:
            Input: seq="ATTTAATTTAATTTC", repair_template=False, num_targets=1, high_specificity=False
            Expected Output: {"recommendation": "Cas12a"}
            Description: AT-rich sequence; Cas12a wins on guide count and AT-richness is noted as support.
        - Case:
            Input: seq="GCGCGCGCGCGG", repair_template=True, num_targets=1, high_specificity=False
            Expected Output: {"recommendation": "Cas9"}
            Description: HDR alone does not change the recommendation.
        - Case:
            Input: seq="GCGCGCGCGCGG", repair_template=False, num_targets=3, high_specificity=False
            Expected Output: {"recommendation": "Cas12a"}
            Description: Guide counts favor Cas9 but are overridden by multiplexing in Step 3 only
                         when they fall within the margin. With default margin 1.5, if Cas9 clearly
                         dominates, Cas9 still wins; otherwise multiplexing tips to Cas12a.
        - Case:
            Input: seq="GCGCGCGCGCGG", repair_template=False, num_targets=1, high_specificity=True
            Expected Output: {"recommendation": "Cas12a"}
            Description: High specificity wins in Step 3 when guide counts are within margin.
        - Case:
            Input: seq=""
            Expected Exception: ValueError
            Description: Empty sequence raises ValueError.
        - Case:
            Input: seq="GCG", repair_template=False, num_targets=0, high_specificity=False
            Expected Exception: ValueError
            Description: num_targets must be at least 1.
        - Case:
            Input: seq="GCGCGCGCGCGG", repair_template=False, num_targets=1, high_specificity=False, margin_threshold=0.5
            Expected Exception: ValueError
            Description: margin_threshold must be >= 1.0.
    """

    def initiate(self) -> None:
        pass

    def _reverse_complement(self, seq: str) -> str:
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join(complement[b] for b in reversed(seq))

    def _has_simple_hairpin(self, spacer: str, k: int = 4) -> bool:
        """
        Simple hairpin proxy:
        checks for short reverse-complement matches within the spacer.
        This approximates potential secondary structure (hairpins).
        """
        comp = {"A": "T", "T": "A", "G": "C", "C": "G"}

        def rev_comp(s):
            return "".join(comp[b] for b in reversed(s))

        for i in range(len(spacer) - k + 1):
            sub = spacer[i:i + k]
            rc = rev_comp(sub)
            for j in range(len(spacer) - k + 1):
                if j != i and spacer[j:j + k] == rc:
                    return True
        return False

    def _spacer_passes_quality(self, spacer: str, return_reason: bool = False):
        gc = (spacer.count("G") + spacer.count("C")) / len(spacer)
        if not (0.30 <= gc <= 0.70):
            return (False, "bad GC") if return_reason else False
        if "TTTT" in spacer:
            return (False, "TTTT run") if return_reason else False
        if any(base * 5 in spacer for base in "ACGT"):
            return (False, "homopolymer") if return_reason else False
        if self._has_simple_hairpin(spacer):
            return (False, "hairpin") if return_reason else False
        return (True, None) if return_reason else True

    def run(
        self,
        seq: str,
        repair_template: bool,
        num_targets: int = 1,
        high_specificity: bool = False,
        system: str = None,
        cas12a_spacer_len: int = 23,
        margin_threshold: float = 1.5,
        debug: bool = False,
    ) -> dict:
        if system is not None:
            if system not in ("Cas9", "Cas12a"):
                raise ValueError("system must be 'Cas9' or 'Cas12a'.")
            return {
                "ngg_count": None,
                "tttv_count": None,
                "cas9_valid_guides": None,
                "cas12a_valid_guides": None,
                "gc_content": None,
                "recommendation": system,
                "rationale": f"User specified {system} directly.",
            }

        seq = seq.upper()

        if not seq:
            raise ValueError("Sequence must not be empty.")
        if any(b not in "ATGC" for b in seq):
            raise ValueError("Sequence contains invalid bases.")
        if len(seq) < 4:
            raise ValueError("Sequence is too short to contain any PAM sites.")
        if num_targets < 1:
            raise ValueError("num_targets must be at least 1.")
        if not (20 <= cas12a_spacer_len <= 24):
            raise ValueError("cas12a_spacer_len must be between 20 and 24.")
        if margin_threshold < 1.0:
            raise ValueError("margin_threshold must be >= 1.0.")

        rc = self._reverse_complement(seq)
        gc_content = (seq.count("G") + seq.count("C")) / len(seq)
        at_rich = gc_content < 0.45
        debug_logs = [] if debug else None

        ngg_count = 0
        tttv_count = 0
        cas9_valid_guides = 0
        cas12a_valid_guides = 0

        for strand in (seq, rc):
            for i in range(len(strand) - 2):
                if strand[i + 1] == "G" and strand[i + 2] == "G":
                    ngg_count += 1
                    if i >= 20:
                        spacer = strand[i - 20:i]
                        if debug:
                            ok, reason = self._spacer_passes_quality(spacer, return_reason=True)
                            if ok:
                                cas9_valid_guides += 1
                            else:
                                debug_logs.append(f"Cas9 reject @ {i}: {reason} | {spacer}")
                        else:
                            if self._spacer_passes_quality(spacer):
                                cas9_valid_guides += 1
            for i in range(len(strand) - 3):
                if (strand[i] == "T" and strand[i + 1] == "T"
                        and strand[i + 2] == "T" and strand[i + 3] in "ACG"):
                    tttv_count += 1
                    if i + 4 + cas12a_spacer_len <= len(strand):
                        spacer = strand[i + 4:i + 4 + cas12a_spacer_len]
                        if debug:
                            ok, reason = self._spacer_passes_quality(spacer, return_reason=True)
                            if ok:
                                cas12a_valid_guides += 1
                            else:
                                debug_logs.append(f"Cas12a reject @ {i}: {reason} | {spacer}")
                        else:
                            if self._spacer_passes_quality(spacer):
                                cas12a_valid_guides += 1

        # --- Step 2: guide count signal ---
        guide_decision = None
        guide_rationale = ""

        if cas9_valid_guides == 0 and cas12a_valid_guides == 0:
            guide_rationale = (
                f"No valid guides found for either nuclease "
                f"({ngg_count} NGG PAMs, {tttv_count} TTTV PAMs; GC {gc_content:.0%})."
            )
        elif cas9_valid_guides == 0:
            guide_decision = "Cas12a"
            guide_rationale = (
                f"Cas9 has no valid guides in this sequence; "
                f"Cas12a has {cas12a_valid_guides} ({tttv_count} TTTV PAMs; GC {gc_content:.0%})."
            )
        elif cas12a_valid_guides == 0:
            guide_decision = "Cas9"
            guide_rationale = (
                f"Cas12a has no valid guides in this sequence; "
                f"Cas9 has {cas9_valid_guides} ({ngg_count} NGG PAMs; GC {gc_content:.0%})."
            )
        else:
            ratio = max(cas9_valid_guides, cas12a_valid_guides) / min(cas9_valid_guides, cas12a_valid_guides)
            if ratio >= margin_threshold:
                if cas9_valid_guides > cas12a_valid_guides:
                    guide_decision = "Cas9"
                    guide_rationale = (
                        f"Guide count favors Cas9: {cas9_valid_guides} valid Cas9 guides vs "
                        f"{cas12a_valid_guides} valid Cas12a guides "
                        f"(ratio {ratio:.1f}x >= threshold {margin_threshold}x; "
                        f"GC {gc_content:.0%})."
                    )
                else:
                    guide_decision = "Cas12a"
                    guide_rationale = (
                        f"Guide count favors Cas12a: {cas12a_valid_guides} valid Cas12a guides vs "
                        f"{cas9_valid_guides} valid Cas9 guides "
                        f"(ratio {ratio:.1f}x >= threshold {margin_threshold}x; "
                        f"GC {gc_content:.0%})."
                    )
            else:
                guide_rationale = (
                    f"Guide counts are within margin: {cas9_valid_guides} valid Cas9 guides vs "
                    f"{cas12a_valid_guides} valid Cas12a guides "
                    f"(ratio {ratio:.1f}x < threshold {margin_threshold}x; "
                    f"GC {gc_content:.0%}); using secondary criteria."
                )

        # --- Step 3: secondary criteria (when guide counts are tied or both zero) ---
        if guide_decision is not None:
            recommendation = guide_decision
            rationale = guide_rationale
        elif num_targets >= 2:
            recommendation = "Cas12a"
            rationale = (
                guide_rationale + f" Multiplexing ({num_targets} targets): Cas12a processes "
                "a crRNA array from a single transcript without a separate transcription "
                "unit per guide (Zetsche et al., 2017, Nat Biotechnol)."
            )
        elif high_specificity:
            recommendation = "Cas12a"
            rationale = (
                guide_rationale + " High specificity required: Cas12a is more mismatch-"
                "intolerant than SpCas9, with off-target burden comparable to engineered "
                "high-fidelity Cas9 variants (Kleinstiver et al., 2016, Nat Biotechnol)."
            )
        else:
            recommendation = "Cas9"
            rationale = guide_rationale + " Defaulting to Cas9 as the well-established standard."

        # --- AT-richness: supporting context only ---
        if at_rich and recommendation == "Cas12a":
            rationale += (
                f" AT-rich sequence (GC {gc_content:.0%}) further supports Cas12a: "
                "its TTTV PAM is more prevalent in AT-rich regions "
                "(Zetsche et al., 2015, Cell 163(3):759-771)."
            )

        # --- Secondary notes on repair strategy ---
        notes = []

        if not repair_template:
            notes.append(
                "NHEJ/knockout: Cas9's blunt DSB and mature knockout protocols make "
                "it the standard choice for indel-based disruption"
            )

        if repair_template:
            notes.append(
                "Templated edit: for a precise SNV or short sequence rewrite, consider "
                "base or prime editing (Cas9-based) as a first option. For a donor "
                "knock-in where staggered overhangs are advantageous, Cas12a is worth "
                "testing in parallel"
            )

        if notes:
            rationale += " | " + "; ".join(notes) + "."

        if debug and debug_logs:
            print("\n--- Guide Rejection Debug ---")
            for log in debug_logs[:20]:
                print(log)
            if len(debug_logs) > 20:
                print(f"... ({len(debug_logs) - 20} more)")

        return {
            "ngg_count": ngg_count,
            "tttv_count": tttv_count,
            "cas9_valid_guides": cas9_valid_guides,
            "cas12a_valid_guides": cas12a_valid_guides,
            "gc_content": round(gc_content, 3),
            "recommendation": recommendation,
            "rationale": rationale,
        }


_instance = CasSelector()
_instance.initiate()
cas_selector = _instance.run
