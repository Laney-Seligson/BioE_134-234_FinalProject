from __future__ import annotations

# Lookup table for common model organisms.
# Each entry: nuclease, delivery method, and a brief rationale.
# Sources: Anzalone et al. 2020 (Nature Biotechnology), Mout et al. 2017 (ACS Nano),
#          Concordet & Haeussler 2018 (NAR), Addgene CRISPR guide.
_SPECIES_DB: dict[str, dict] = {
    "human": {
        "nuclease": "Cas9",
        "delivery": "RNP electroporation or lentiviral vector",
        "rationale": (
            "Human cells are well-characterized for SpCas9; RNP electroporation "
            "minimizes off-targets by limiting Cas9 exposure time, while lentiviral "
            "delivery suits stable integration workflows."
        ),
    },
    "mouse": {
        "nuclease": "Cas9",
        "delivery": "RNP microinjection (zygote) or AAV (somatic)",
        "rationale": (
            "Mouse zygote microinjection with Cas9 RNP is the gold standard for "
            "germline edits; AAV serotypes (AAV9, AAVrh10) are preferred for "
            "in vivo somatic delivery."
        ),
    },
    "zebrafish": {
        "nuclease": "Cas9",
        "delivery": "RNP microinjection into 1-cell embryo",
        "rationale": (
            "Zebrafish embryos are highly amenable to direct RNP injection; "
            "the large, transparent embryo allows visual confirmation of injection "
            "success. Cas9 is preferred due to abundant NGG PAM sites in the "
            "zebrafish GC-balanced genome."
        ),
    },
    "drosophila": {
        "nuclease": "Cas9",
        "delivery": "Plasmid injection or transgenic Cas9 line",
        "rationale": (
            "Fly genetics favor stable transgenic Cas9/gRNA lines crossed to "
            "target strains; direct embryo injection with plasmid is also common "
            "for one-off experiments."
        ),
    },
    "yeast": {
        "nuclease": "Cas9",
        "delivery": "Plasmid transformation (lithium acetate)",
        "rationale": (
            "S. cerevisiae has highly efficient homologous recombination, making "
            "Cas9 + donor template plasmid transformation the simplest approach. "
            "The AT-rich genome still has sufficient NGG sites for Cas9."
        ),
    },
    "plant": {
        "nuclease": "Cas12a",
        "delivery": "Agrobacterium-mediated transformation or biolistics",
        "rationale": (
            "Plant genomes are often AT-rich, favoring Cas12a (TTTV PAM). "
            "Agrobacterium T-DNA delivery is most common for dicots; "
            "biolistics (gene gun) is preferred for monocots like maize and wheat."
        ),
    },
    "e_coli": {
        "nuclease": "Cas9",
        "delivery": "Plasmid electroporation",
        "rationale": (
            "E. coli is straightforward: transform Cas9 + gRNA plasmid by "
            "electroporation. Lambda Red recombineering is often combined with "
            "Cas9 to improve editing efficiency."
        ),
    },
    "rat": {
        "nuclease": "Cas9",
        "delivery": "RNP microinjection (zygote) or AAV (somatic)",
        "rationale": (
            "Similar to mouse; RNP microinjection into rat zygotes is the standard "
            "germline approach. AAV8 and AAV9 are effective for liver and CNS "
            "somatic delivery respectively."
        ),
    },
}

# Normalize species name input to match database keys
_ALIASES: dict[str, str] = {
    "homo sapiens": "human",
    "mus musculus": "mouse",
    "danio rerio": "zebrafish",
    "drosophila melanogaster": "drosophila",
    "fly": "drosophila",
    "fruit fly": "drosophila",
    "saccharomyces cerevisiae": "yeast",
    "s. cerevisiae": "yeast",
    "arabidopsis": "plant",
    "arabidopsis thaliana": "plant",
    "escherichia coli": "e_coli",
    "ecoli": "e_coli",
    "rattus norvegicus": "rat",
}


class SpeciesRecommender:
    """
    Description:
        Given a target species, recommends which CRISPR nuclease to use (Cas9 vs
        Cas12a) and the most appropriate delivery method for that organism. This
        extends the sequence-based cas_selector tool by adding species-aware
        biological context: delivery constraints, genome characteristics, and
        established experimental protocols.

        Covers 8 common model organisms: human, mouse, zebrafish, Drosophila,
        yeast, plant, E. coli, and rat. Accepts common aliases and Latin names.

    Input:
        species (str): common name or Latin name of the target organism.
                       e.g. "mouse", "Mus musculus", "zebrafish", "plant"

    Output:
        dict with keys:
            - species: normalized species name
            - nuclease: recommended CRISPR nuclease ("Cas9" or "Cas12a")
            - delivery: recommended delivery method
            - rationale: biological justification for the recommendation
            - supported: True if the species is in the database, False otherwise
            - note: present only when species is unsupported; suggests a fallback

    Tests:
        - Case:
            Input: species="mouse"
            Expected Output: nuclease="Cas9", "RNP" in delivery
            Description: mouse → Cas9 + RNP microinjection.
        - Case:
            Input: species="plant"
            Expected Output: nuclease="Cas12a"
            Description: AT-rich plant genomes → Cas12a.
        - Case:
            Input: species="Mus musculus"
            Expected Output: nuclease="Cas9"
            Description: Latin name alias resolves to mouse.
        - Case:
            Input: species=""
            Expected Exception: ValueError
            Description: empty species raises ValueError.
        - Case:
            Input: species="tardigrade"
            Expected Output: supported=False
            Description: unknown species returns unsupported=False with a fallback note.
    """

    def initiate(self) -> None:
        pass

    def run(self, species: str) -> dict:
        if not species or not species.strip():
            raise ValueError("Species name must not be empty.")

        normalized = species.strip().lower().replace("-", "_")

        # resolve alias first, then look up directly
        key = _ALIASES.get(normalized, normalized)
        entry = _SPECIES_DB.get(key)

        if entry is None:
            return {
                "species": species.strip(),
                "nuclease": None,
                "delivery": None,
                "rationale": None,
                "supported": False,
                "note": (
                    f"'{species}' is not in the species database. "
                    "Use crispr_cas_selector with the target sequence to get a "
                    "nuclease recommendation based on GC content, then consult "
                    "the literature for delivery methods specific to your organism."
                ),
            }

        return {
            "species": key,
            "nuclease": entry["nuclease"],
            "delivery": entry["delivery"],
            "rationale": entry["rationale"],
            "supported": True,
        }


_instance = SpeciesRecommender()
_instance.initiate()
species_recommender = _instance.run
