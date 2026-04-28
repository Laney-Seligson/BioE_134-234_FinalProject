from __future__ import annotations

from typing import Optional

from modules.crispr_tools.tools.citations import cites, format_citations


# Quality and efficiency thresholds drawn from published guidance and the
# Synthego ICE / Brinkman TIDE manuals.
#   - R^2 >= 0.9: high-quality fit (Brinkman et al. 2014, Nucleic Acids Res
#     42(22):e168). Below 0.9 the indel deconvolution is unreliable and the
#     reported percentages should not be trusted.
#   - KO score / indel %: experiment-specific, but standard cutoffs are:
#       <10%   = failed / no detectable editing
#       10-30% = marginal; consider re-transfecting or re-picking colonies
#       30-70% = good for a population
#       >70%   = excellent, comparable to RNP delivery benchmarks
#     (Synthego ICE documentation; Hsiau et al. 2019 bioRxiv 251082)
_R2_GOOD = 0.90
_R2_MARGINAL = 0.80

_EDITING_FAILED = 10.0
_EDITING_MARGINAL = 30.0
_EDITING_GOOD = 70.0


def _classify_efficiency(editing_pct: float) -> str:
    if editing_pct < _EDITING_FAILED:
        return "FAILED"
    if editing_pct < _EDITING_MARGINAL:
        return "MARGINAL"
    if editing_pct < _EDITING_GOOD:
        return "GOOD"
    return "EXCELLENT"


def _classify_fit_quality(r_squared: float) -> str:
    if r_squared >= _R2_GOOD:
        return "HIGH"
    if r_squared >= _R2_MARGINAL:
        return "MODERATE"
    return "LOW"


class InterpretIceTide:
    """
    Description:
        Parses ICE (Synthego) or TIDE (Brinkman et al. 2014) output and
        translates the raw numbers into a plain-English interpretation
        with actionable next steps. Closes the loop on the verify_edit
        protocol: verify_edit tells the user how to set up ICE/TIDE,
        and this tool tells them what to do with the result.

        Inputs accepted:
          - editing_pct:  the indel/KO percentage reported by ICE or TIDE
          - r_squared:    the fit quality (R^2) reported by the tool
          - tool:         "ice" or "tide" (affects naming in the report only)
          - indel_distribution (optional): dict mapping indel size (e.g.
            "+1", "-3", "0") to percentage. If supplied, the most common
            indel and its frequency are highlighted.
          - sample_id (optional): label for the report.

        Decision logic:
          1. If R^2 < 0.80, the result is unreliable regardless of the
             editing percentage. Flag for re-sequencing.
          2. Otherwise classify editing percentage:
             <10% FAILED, 10-30% MARGINAL, 30-70% GOOD, >=70% EXCELLENT.
          3. If indel_distribution is supplied, identify the dominant indel
             and warn if 0 (unedited) is the largest fraction even when
             editing_pct is non-trivial (suggests background noise).

        Sources:
          - Brinkman et al. 2014, Nucleic Acids Res 42(22):e168 (TIDE)
          - Hsiau et al. 2019, bioRxiv 251082 (ICE)
          - Synthego ICE Analysis documentation

    Input:
        editing_pct (float): percent of reads with indels (KO score for ICE,
                             total efficiency for TIDE). Range 0-100.
        r_squared (float):   R^2 of the deconvolution fit. Range 0-1.
        tool (str):          "ice" (default) or "tide". Cosmetic.
        indel_distribution (dict, optional): {"+1": 45.2, "-3": 12.1, ...}
        sample_id (str, optional): label included in the report.

    Output:
        dict with keys:
            - sample_id
            - tool
            - editing_pct
            - r_squared
            - efficiency_classification: FAILED / MARGINAL / GOOD / EXCELLENT
            - fit_quality: HIGH / MODERATE / LOW
            - is_reliable: bool (False if R^2 < 0.80)
            - dominant_indel: most common indel from indel_distribution, or None
            - dominant_indel_pct: its percentage, or None
            - warnings: list of human-readable warnings
            - next_steps: list of recommended actions
            - summary: one-paragraph plain-English overall verdict

    Tests:
        - Case:
            Input: editing_pct=85, r_squared=0.95, tool="ice"
            Expected Output: efficiency_classification == "EXCELLENT", is_reliable == True
            Description: high editing + good fit -> success.
        - Case:
            Input: editing_pct=5, r_squared=0.95, tool="ice"
            Expected Output: efficiency_classification == "FAILED"
            Description: low editing despite good fit -> failed experiment.
        - Case:
            Input: editing_pct=85, r_squared=0.5, tool="ice"
            Expected Output: is_reliable == False, warning about R^2
            Description: high editing but bad fit -> result not trustworthy.
        - Case:
            Input: editing_pct=-5, r_squared=0.9
            Expected Exception: ValueError
            Description: editing_pct out of range raises.
        - Case:
            Input: editing_pct=50, r_squared=1.5
            Expected Exception: ValueError
            Description: r_squared out of range raises.
        - Case:
            Input: editing_pct=50, r_squared=0.9, indel_distribution={"+1": 30, "0": 50, "-3": 20}
            Expected Output: warning that unedited reads dominate
            Description: dominant indel is 0 -> background noise warning.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        editing_pct: float,
        r_squared: float,
        tool: str = "ice",
        indel_distribution: Optional[dict] = None,
        sample_id: Optional[str] = None,
    ) -> dict:

        tool = tool.lower().strip()
        if tool not in {"ice", "tide"}:
            raise ValueError(f"tool must be 'ice' or 'tide', got '{tool}'.")

        if not (0 <= editing_pct <= 100):
            raise ValueError(
                f"editing_pct must be between 0 and 100, got {editing_pct}."
            )
        if not (0 <= r_squared <= 1):
            raise ValueError(
                f"r_squared must be between 0 and 1, got {r_squared}."
            )

        efficiency_class = _classify_efficiency(editing_pct)
        fit_class = _classify_fit_quality(r_squared)
        is_reliable = r_squared >= _R2_MARGINAL

        warnings: list[str] = []
        next_steps: list[str] = []

        if not is_reliable:
            warnings.append(
                f"R^2 = {r_squared:.2f} is below 0.80 — the indel deconvolution "
                f"is unreliable. The reported {editing_pct:.1f}% editing should not "
                "be trusted. Re-sequence with a cleaner trace."
            )
            next_steps.append(
                "Re-PCR and re-sequence the sample. Check primer specificity, "
                "trace quality (SNR > 20), and that the unedited control trace "
                "is clean."
            )
        elif fit_class == "MODERATE":
            warnings.append(
                f"R^2 = {r_squared:.2f} is acceptable but not ideal (< 0.90). "
                "Treat the editing percentage as approximate."
            )

        # Indel distribution analysis
        dominant_indel = None
        dominant_pct = None
        if indel_distribution:
            cleaned = {k: float(v) for k, v in indel_distribution.items()}
            dominant_indel = max(cleaned, key=cleaned.get)
            dominant_pct = cleaned[dominant_indel]

            # If the "0" (unedited) bin is the largest but editing_pct claims
            # significant editing, that's contradictory and worth flagging.
            unedited_pct = cleaned.get("0", 0.0)
            if dominant_indel == "0" and editing_pct > _EDITING_FAILED:
                warnings.append(
                    f"Unedited reads ('0') are the dominant fraction at "
                    f"{unedited_pct:.1f}% even though overall editing is "
                    f"reported at {editing_pct:.1f}%. The non-zero indel "
                    "channels may be background noise — verify with an "
                    "unedited control trace."
                )

        # Efficiency-driven next steps (only meaningful if reliable)
        if is_reliable:
            if efficiency_class == "FAILED":
                next_steps.append(
                    "Editing failed (< 10%). Possible causes: (1) gRNA didn't "
                    "load — check delivery method, (2) target site is "
                    "inaccessible (chromatin), (3) PAM was misidentified. "
                    "Re-check guide design with predict_offtargets and "
                    "confirm the protospacer matches the reference."
                )
            elif efficiency_class == "MARGINAL":
                next_steps.append(
                    "Editing is marginal (10-30%). Consider re-transfecting "
                    "with a higher dose, switching to RNP delivery, or "
                    "picking more colonies (use colony_calculator to estimate)."
                )
            elif efficiency_class == "GOOD":
                next_steps.append(
                    "Editing is good (30-70%). Pick colonies and screen by "
                    "Sanger sequencing of single clones to recover homozygous "
                    "edited lines."
                )
            else:  # EXCELLENT
                next_steps.append(
                    "Editing is excellent (>= 70%). Proceed to clone isolation; "
                    "even a small number of picked colonies should yield "
                    "homozygous edits."
                )

        # Build summary paragraph
        sample_label = f"Sample '{sample_id}': " if sample_id else ""
        verdict_text = {
            "EXCELLENT": "highly successful editing",
            "GOOD": "successful editing",
            "MARGINAL": "marginal editing",
            "FAILED": "failed editing",
        }[efficiency_class]

        if not is_reliable:
            summary = (
                f"{sample_label}{tool.upper()} reports {editing_pct:.1f}% editing "
                f"with R^2 = {r_squared:.2f}. The fit quality is too low to "
                "trust this result — re-sequence before drawing conclusions."
            )
        else:
            summary = (
                f"{sample_label}{tool.upper()} reports {editing_pct:.1f}% editing "
                f"with R^2 = {r_squared:.2f} ({fit_class.lower()} fit quality). "
                f"This indicates {verdict_text} ({efficiency_class.lower()})."
            )
            if dominant_indel is not None:
                summary += (
                    f" Dominant indel: {dominant_indel} at {dominant_pct:.1f}%."
                )

        # Cite the paper that defines the analysis tool the user ran.
        if tool == "ice":
            citation_keys = ["hsiau_2019_ice"]
        else:
            citation_keys = ["brinkman_2014_tide"]

        return {
            "sample_id": sample_id,
            "tool": tool,
            "editing_pct": editing_pct,
            "r_squared": r_squared,
            "efficiency_classification": efficiency_class,
            "fit_quality": fit_class,
            "is_reliable": is_reliable,
            "dominant_indel": dominant_indel,
            "dominant_indel_pct": dominant_pct,
            "warnings": warnings,
            "next_steps": next_steps,
            "summary": summary,
            "citations": format_citations(cites(*citation_keys)),
        }


_instance = InterpretIceTide()
_instance.initiate()
interpret_ice_tide = _instance.run
