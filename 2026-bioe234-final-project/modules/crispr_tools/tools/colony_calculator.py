from __future__ import annotations

import math


def _binomial_at_least_k(n: int, p: float, k: int) -> float:
    """
    Probability of at least k successes in n Bernoulli(p) trials.
    Computed as 1 - sum_{i=0}^{k-1} C(n,i) p^i (1-p)^(n-i).
    """
    if k <= 0:
        return 1.0
    if n < k:
        return 0.0
    cumulative = 0.0
    for i in range(k):
        cumulative += math.comb(n, i) * (p ** i) * ((1 - p) ** (n - i))
    return 1.0 - cumulative


# Typical editing efficiency ranges by experimental context, drawn from
# published benchmarks. Numbers are rough order-of-magnitude defaults
# users can override.
#   - Cas9 + RNP delivery to mammalian cells: ~50-80% (Kim et al. 2014;
#     Lin et al. 2014, eLife 3:e04766)
#   - Cas9 plasmid transfection of mammalian cells: ~10-30%
#   - Cas9 in E. coli (lambda Red recombineering or pCRISPR): ~30-70%
#     (Jiang et al. 2013, Nat Biotechnol 31:233-239)
#   - HDR (precise edits): ~1-10% in most cell types
#     (Paquet et al. 2016, Nature 533:125-129)
_EFFICIENCY_PRESETS = {
    "cas9_rnp_mammalian": 0.65,
    "cas9_plasmid_mammalian": 0.20,
    "cas9_ecoli": 0.50,
    "cas12a_ecoli": 0.40,
    "hdr_mammalian": 0.05,
    "hdr_ecoli": 0.10,
}


class ColonyCalculator:
    """
    Description:
        Given an expected editing efficiency, computes how many colonies (or
        clones) to screen to be confident of recovering a desired number of
        successfully edited clones. Uses the binomial distribution: each
        colony is an independent Bernoulli trial with success probability
        equal to the editing efficiency.

        Why this matters: a common newcomer mistake is to pick 3-5 colonies
        when efficiency is low (e.g. HDR at 5%), find none edited, and
        conclude the experiment failed when in fact insufficient colonies
        were screened. With p = 0.05 you need ~59 colonies for a 95% chance
        of finding even one edit.

        Math:
            For "at least k edited clones with confidence C":
                find smallest n such that P(X >= k) >= C
                where X ~ Binomial(n, p)

            For the special case k=1:
                n >= log(1 - C) / log(1 - p)

        Sources:
          - Kim et al. 2014, Genome Res 24:1012-1019 (Cas9 RNP efficiency)
          - Paquet et al. 2016, Nature 533:125-129 (HDR efficiency benchmarks)
          - Jiang et al. 2013, Nat Biotechnol 31:233-239 (Cas9 in bacteria)

    Input:
        editing_efficiency (float): expected fraction of colonies that will
            be successfully edited. Range (0, 1]. Required UNLESS preset is given.
        desired_clones (int): how many edited clones you want to recover.
            Default 1.
        confidence (float): probability of recovering at least desired_clones.
            Default 0.95.
        preset (str, optional): use a published efficiency benchmark instead
            of supplying editing_efficiency. One of:
            'cas9_rnp_mammalian' (0.65), 'cas9_plasmid_mammalian' (0.20),
            'cas9_ecoli' (0.50), 'cas12a_ecoli' (0.40),
            'hdr_mammalian' (0.05), 'hdr_ecoli' (0.10).
        max_colonies (int): safety cap on the search. Default 10000.

    Output:
        dict with keys:
            - editing_efficiency: efficiency used in the calculation
            - desired_clones: how many edits we want
            - confidence: target confidence level
            - colonies_to_pick: minimum n satisfying P(X >= k) >= C
            - expected_edits: n * p (mean of the binomial)
            - probability_at_chosen_n: actual P(X >= k) at the recommended n
            - safety_margin_recommendation: a 1.5x bump for real-world losses
              (failed PCR, contaminated wells, sequencing fails)
            - recommendation: one-paragraph plain-English advice

    Tests:
        - Case:
            Input: editing_efficiency=0.5, desired_clones=1, confidence=0.95
            Expected Output: colonies_to_pick == 5
            Description: P(X>=1) = 1-0.5^n >= 0.95 -> n=5 (1-0.5^5=0.969).
        - Case:
            Input: editing_efficiency=0.05, desired_clones=1, confidence=0.95
            Expected Output: colonies_to_pick == 59
            Description: low-efficiency HDR case; 59 colonies for 95% confidence.
        - Case:
            Input: editing_efficiency=0, desired_clones=1
            Expected Exception: ValueError
            Description: zero efficiency means no number of colonies suffices.
        - Case:
            Input: editing_efficiency=1.5
            Expected Exception: ValueError
            Description: efficiency must be in (0, 1].
        - Case:
            Input: preset="hdr_mammalian", desired_clones=3
            Expected Output: colonies_to_pick > 60
            Description: 3 HDR clones at 5% efficiency requires many colonies.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        editing_efficiency: float = None,
        desired_clones: int = 1,
        confidence: float = 0.95,
        preset: str = None,
        max_colonies: int = 10000,
    ) -> dict:

        if preset is not None:
            preset = preset.lower().strip()
            if preset not in _EFFICIENCY_PRESETS:
                raise ValueError(
                    f"Unknown preset '{preset}'. Available: {sorted(_EFFICIENCY_PRESETS)}."
                )
            editing_efficiency = _EFFICIENCY_PRESETS[preset]

        if editing_efficiency is None:
            raise ValueError(
                "Must provide either editing_efficiency or preset."
            )
        if not (0 < editing_efficiency <= 1):
            raise ValueError(
                f"editing_efficiency must be in (0, 1], got {editing_efficiency}."
            )
        if desired_clones < 1:
            raise ValueError(f"desired_clones must be >= 1, got {desired_clones}.")
        if not (0 < confidence < 1):
            raise ValueError(f"confidence must be in (0, 1), got {confidence}.")

        # find smallest n such that P(X >= desired_clones) >= confidence
        n = desired_clones
        while n <= max_colonies:
            prob = _binomial_at_least_k(n, editing_efficiency, desired_clones)
            if prob >= confidence:
                break
            n += 1
        else:
            raise ValueError(
                f"Could not reach confidence {confidence} for {desired_clones} edits "
                f"at efficiency {editing_efficiency} within {max_colonies} colonies. "
                "Check inputs or raise max_colonies."
            )

        actual_prob = _binomial_at_least_k(n, editing_efficiency, desired_clones)
        expected_edits = round(n * editing_efficiency, 2)
        # 1.5x bump accounts for typical real-world losses: failed PCR,
        # well contamination, sequencing dropouts. Empirical rule of thumb.
        safety_n = math.ceil(n * 1.5)

        if editing_efficiency < 0.10:
            advice = (
                "Low editing efficiency — screening burden is high. Consider "
                "improving delivery (RNP > plasmid), adding selection markers, "
                "or using a co-selection strategy."
            )
        elif editing_efficiency < 0.30:
            advice = (
                "Moderate editing efficiency. Pick the recommended number plus "
                "a safety margin to account for failed PCRs and contamination."
            )
        else:
            advice = (
                "Good editing efficiency — relatively few colonies needed. "
                "A standard 24-well plate is usually more than sufficient."
            )

        recommendation = (
            f"Pick {n} colonies to be {confidence:.0%} confident of recovering "
            f"at least {desired_clones} edited clone(s) at {editing_efficiency:.0%} "
            f"editing efficiency. Recommended safety total: {safety_n} colonies "
            f"(1.5x bump for failed PCR/contamination losses). {advice}"
        )

        return {
            "editing_efficiency": editing_efficiency,
            "desired_clones": desired_clones,
            "confidence": confidence,
            "colonies_to_pick": n,
            "expected_edits": expected_edits,
            "probability_at_chosen_n": round(actual_prob, 4),
            "safety_margin_recommendation": safety_n,
            "recommendation": recommendation,
        }


_instance = ColonyCalculator()
_instance.initiate()
colony_calculator = _instance.run
