from __future__ import annotations

from typing import Optional

from modules.crispr_tools.tools.predict_offtargets import predict_offtargets
from modules.crispr_tools.tools._citations import cites, format_citations

_SUPPORTED_NUCLEASES = {"cas9", "cas12a"}


def _score_efficiency(protospacer: str, nuclease: str) -> tuple[int, dict]:
    """
    Score a protospacer on predicted on-target cutting efficiency.

    Cas9 rules (Doench et al. 2016, Nat Biotechnol 34:184-191):
      - GC content 40-70% (binding energy sweet spot)               +1
      - No TTTT run (Pol III terminator would truncate the sgRNA)   +1
      - G at PAM-proximal end (position 20; the base adjacent       +1
        to the NGG PAM correlates with higher cleavage activity)
      Max score: 3.

    Cas12a rules (Kim et al. 2018; Zetsche et al. 2015):
      - GC content 40-70%                                           +1
      - No TTTT run in the spacer (4+ T's look like another PAM)    +1
      Max score: 2. (Fewer PAM-adjacent preference rules are well-
      established for Cas12a than for Cas9.)
    """
    gc_count = sum(1 for b in protospacer if b in "GC")
    gc_content = gc_count / len(protospacer) if protospacer else 0
    gc_ok = 0.40 <= gc_content <= 0.70

    no_polyt = "TTTT" not in protospacer

    score = 0
    details = {
        "gc_content": round(gc_content, 3),
        "gc_content_ok": gc_ok,
        "no_polyt_run": no_polyt,
    }

    if gc_ok:
        score += 1
    if no_polyt:
        score += 1

    if nuclease == "cas9":
        # position 20 (PAM-proximal) = last base of the protospacer
        g_at_pam_proximal = len(protospacer) > 0 and protospacer[-1] == "G"
        details["g_at_pam_proximal"] = g_at_pam_proximal
        if g_at_pam_proximal:
            score += 1

    return score, details


def _score_specificity(
    protospacer: str, reference: str, nuclease: str, max_mismatches: int
) -> tuple[int, dict]:
    """
    Score a protospacer on off-target specificity by calling predict_offtargets.

    Rules:
      - No HIGH risk off-target sites (on-target counts as HIGH but
        is excluded from this tally)                                 +1
      - ≤1 MEDIUM risk off-target site                                +1
      - Total off-target sites ≤5                                     +1
    Max score: 3.
    """
    report = predict_offtargets(
        protospacer=protospacer,
        reference=reference,
        nuclease=nuclease,
        max_mismatches=max_mismatches,
    )

    # on-target = first 0-mismatch HIGH-risk site; exclude it from off-target counts
    sites = report["offtarget_sites"]
    off_sites = []
    on_target_seen = False
    for s in sites:
        if not on_target_seen and s["mismatches"] == 0 and s["risk"] == "HIGH":
            on_target_seen = True
            continue
        off_sites.append(s)

    high = sum(1 for s in off_sites if s["risk"] == "HIGH")
    medium = sum(1 for s in off_sites if s["risk"] == "MEDIUM")
    total = len(off_sites)

    score = 0
    if high == 0:
        score += 1
    if medium <= 1:
        score += 1
    if total <= 5:
        score += 1

    details = {
        "high_risk_offtargets": high,
        "medium_risk_offtargets": medium,
        "total_offtargets": total,
        "on_target_found": on_target_seen,
    }
    return score, details


class RankGuides:
    """
    Description:
        Ranks a list of gRNA candidates (from design_cas9_grna or
        design_cas12a_crrna) by combined on-target efficiency and
        off-target specificity scoring. Replaces the default "pick
        guides[0]" behavior in run_full_crispr_workflow with an
        evidence-based guide selection step.

        Scoring has two components:

        1. Efficiency (on-target cutting likelihood)
           - GC content 40-70%                           +1
           - No TTTT run (Pol III terminator)            +1
           - G at PAM-proximal end (Cas9 only)           +1
           Sources: Doench et al. 2016 (Nat Biotechnol);
                    Briner et al. 2014 (Mol Cell);
                    Zetsche et al. 2015 (Cell).

        2. Specificity (off-target safety) — uses predict_offtargets
           - No HIGH risk off-target sites                +1
           - ≤1 MEDIUM risk off-target site               +1
           - Total off-target sites ≤5                    +1
           Sources: Hsu et al. 2013 (Nat Biotechnol);
                    Fu et al. 2013 (Nat Biotechnol).

        Total score = efficiency + specificity. Max 6 for Cas9, 5 for Cas12a.
        Guides are returned sorted by total score (desc), then by
        efficiency (desc), then by total off-targets (asc).

    Input:
        guides (list[dict]): list of guide dicts as returned by
                             design_cas9_grna or design_cas12a_crrna.
                             Each must have a "protospacer" key.
        reference (str):     genome/plasmid sequence to scan for off-targets.
                             Resource name, raw string, FASTA, or GenBank.
        nuclease (str):      "cas9" (default) or "cas12a".
        max_mismatches (int): mismatch threshold passed through to
                              predict_offtargets. Default 3.

    Output:
        dict with keys:
            - ranked_guides: input guides with added efficiency_score,
                             specificity_score, total_score, and detail dicts;
                             sorted best-first.
            - best_guide: the top-ranked guide (ranked_guides[0]).
            - scoring_rationale: one-paragraph plain-English summary
                                 of why the best guide won.

    Tests:
        - Case:
            Input: guides=[{"protospacer": "ATGATGATGATGATGATGAG"}],
                   reference="ATGATGATGATGATGATGAGAGG"+"C"*50, nuclease="cas9"
            Expected Output: ranked_guides[0]["efficiency_score"] == 3
            Description: GC=50%, no TTTT, G at pos 20 → max efficiency.
        - Case:
            Input: guides=[], reference="AAA", nuclease="cas9"
            Expected Exception: ValueError
            Description: empty guide list → error.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        guides: list,
        reference: str,
        nuclease: str = "cas9",
        max_mismatches: int = 3,
    ) -> dict:

        nuclease = nuclease.lower().strip()
        if nuclease not in _SUPPORTED_NUCLEASES:
            raise ValueError(
                f"Unsupported nuclease '{nuclease}'. Choose 'cas9' or 'cas12a'."
            )
        if not guides:
            raise ValueError("guides must not be empty.")
        if not reference or not reference.strip():
            raise ValueError("reference must not be empty.")

        scored = []
        for guide in guides:
            protospacer = guide.get("protospacer")
            if not protospacer:
                raise ValueError("Each guide must have a 'protospacer' key.")

            eff_score, eff_details = _score_efficiency(protospacer, nuclease)
            spec_score, spec_details = _score_specificity(
                protospacer, reference, nuclease, max_mismatches
            )

            scored.append({
                **guide,
                "efficiency_score": eff_score,
                "efficiency_details": eff_details,
                "specificity_score": spec_score,
                "specificity_details": spec_details,
                "total_score": eff_score + spec_score,
            })

        # primary: total score desc; secondary: efficiency desc; tertiary: total off-targets asc
        scored.sort(
            key=lambda g: (
                -g["total_score"],
                -g["efficiency_score"],
                g["specificity_details"]["total_offtargets"],
            )
        )

        best = scored[0]
        rationale = (
            f"Best guide '{best['protospacer']}' scored {best['total_score']} "
            f"(efficiency {best['efficiency_score']}, specificity {best['specificity_score']}). "
            f"GC content: {best['efficiency_details']['gc_content']:.0%}. "
            f"Off-target sites: {best['specificity_details']['total_offtargets']} total "
            f"({best['specificity_details']['high_risk_offtargets']} HIGH, "
            f"{best['specificity_details']['medium_risk_offtargets']} MEDIUM)."
        )

        # Citations are based on which scoring rules were applied:
        # - Doench 2016 for the on-target efficiency rules
        # - Hsu 2013 for off-target seed-region scoring (via predict_offtargets)
        # - Cas12a-specific work when nuclease=cas12a
        citation_keys = ["doench_2016", "hsu_2013"]
        if nuclease == "cas12a":
            citation_keys.extend(["zetsche_2015", "kim_2016_cas12a", "kim_2018_cas12a"])

        return {
            "ranked_guides": scored,
            "best_guide": best,
            "scoring_rationale": rationale,
            "citations": format_citations(cites(*citation_keys)),
        }


_instance = RankGuides()
_instance.initiate()
rank_guides = _instance.run
