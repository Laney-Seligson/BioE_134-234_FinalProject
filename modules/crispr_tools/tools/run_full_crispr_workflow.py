from __future__ import annotations

from typing import Optional


# ---------------------------------------------------------------------------
# Organism → vector recommendations (with citations)
# ---------------------------------------------------------------------------

_VECTOR_RECOMMENDATIONS: dict[str, list[dict]] = {
    "mammalian": [
        {
            "vector_key": "px330",
            "vector_name": "pX330",
            "use_case": "Transient transfection of mammalian cells (MEFs, HEK293, primary cells). All-in-one SpCas9 + sgRNA; fastest path from spacer to edit.",
            "delivery": "Lipofection or electroporation",
            "citation": "Cong et al. Science 2013. doi:10.1126/science.1231143",
        },
        {
            "vector_key": "lenticrispr_v2",
            "vector_name": "lentiCRISPR v2",
            "use_case": "Stable Cas9 + sgRNA integration in mammalian cell lines via lentiviral transduction; puromycin selection.",
            "delivery": "Lentiviral transduction",
            "citation": "Sanjana et al. Nat Methods 2014. doi:10.1038/nmeth.3047",
        },
        {
            "vector_key": "px458_gibson",
            "vector_name": "pX458",
            "use_case": "SpCas9 + GFP reporter; FACS-based enrichment of transfected cells improves editing efficiency.",
            "delivery": "Lipofection or electroporation",
            "citation": "Ran et al. Nat Protoc 2013. doi:10.1038/nprot.2013.143",
        },
    ],
    "zebrafish": [
        {
            "vector_key": "pdr274",
            "vector_name": "pDR274",
            "use_case": "T7-driven sgRNA for zebrafish embryo microinjection. Co-inject with Cas9 mRNA or protein into 1-cell embryos.",
            "delivery": "Microinjection into 1-cell embryo",
            "citation": "Hwang et al. PLoS One 2013. doi:10.1371/journal.pone.0068708; Varshney et al. Nat Protoc 2016. doi:10.1038/nprot.2016.099",
        },
    ],
    "ecoli": [
        {
            "vector_key": "pcrispr",
            "vector_name": "pCRISPR::rpsL",
            "use_case": "E. coli genome editing via two-plasmid system (pCas9 + pCRISPR). rpsL counter-selection improves recombinant recovery.",
            "delivery": "Electroporation",
            "citation": "Jiang et al. Nat Biotechnol 2013. doi:10.1038/nbt.2508",
        },
        {
            "vector_key": "ptargetf",
            "vector_name": "pTargetF",
            "use_case": "E. coli sgRNA delivery; pairs with pCas9-CR4 for chromosomal editing.",
            "delivery": "Electroporation",
            "citation": "Jiang et al. Nat Biotechnol 2015. doi:10.1038/nbt.3234",
        },
    ],
    "yeast": [
        {
            "vector_key": "pml104",
            "vector_name": "pML104",
            "use_case": "All-in-one S. cerevisiae CRISPR vector. SpCas9 under TEF1 promoter; sgRNA under SNR52 RNA Pol III promoter. Episomal 2μ origin; URA3 selection.",
            "delivery": "Electroporation / lithium acetate transformation",
            "citation": "Laughery et al. Yeast 2015. doi:10.1002/yea.3080",
        },
    ],
    "plant": [
        {
            "vector_key": "pkse401",
            "vector_name": "pKSE401",
            "use_case": "Binary vector for Agrobacterium-mediated CRISPR in Arabidopsis and other dicots. Expresses SpCas9 under 35S and sgRNA under AtU6 promoter.",
            "delivery": "Agrobacterium tumefaciens (floral dip or tissue culture)",
            "citation": "Xing et al. BMC Plant Biol 2014. doi:10.1186/s12870-014-0327-y",
        },
        {
            "vector_key": "phee401e",
            "vector_name": "pHEE401E",
            "use_case": "High-efficiency Arabidopsis CRISPR using egg cell-specific promoter for Cas9; improves heritable editing rates in T1 plants.",
            "delivery": "Agrobacterium tumefaciens (floral dip)",
            "citation": "Wang et al. Plant Cell 2015. doi:10.1105/tpc.15.00454",
        },
        {
            "vector_key": "pcbc_dt1t2",
            "vector_name": "pCBC-DT1T2",
            "use_case": "Dual-sgRNA binary vector for multiplex editing in Arabidopsis; Golden Gate assembly of two sgRNAs into a single T-DNA.",
            "delivery": "Agrobacterium tumefaciens (floral dip)",
            "citation": "Xing et al. BMC Plant Biol 2014. doi:10.1186/s12870-014-0327-y",
        },
    ],
}

_MAMMALIAN_KEYWORDS = {
    "human", "homo sapiens", "mouse", "mus musculus", "rat", "rattus",
    "hamster", "monkey", "primate", "mammal", "hela", "hek", "cho",
    "293", "jurkat", "macaque",
}
_ZEBRAFISH_KEYWORDS = {"zebrafish", "danio rerio", "danio"}
_ECOLI_KEYWORDS = {"escherichia", "e. coli", "ecoli", "e.coli"}
_PLANT_KEYWORDS = {
    "arabidopsis", "arabidopsis thaliana", "plant", "tobacco", "nicotiana",
    "maize", "zea mays", "rice", "oryza sativa", "tomato", "solanum lycopersicum",
    "wheat", "triticum", "soybean", "glycine max", "poplar", "populus",
}
_YEAST_KEYWORDS = {"saccharomyces", "s. cerevisiae", "yeast", "cerevisiae"}


def _classify_organism(organism: str) -> str:
    o = organism.lower()
    if any(kw in o for kw in _ECOLI_KEYWORDS):
        return "ecoli"
    if any(kw in o for kw in _ZEBRAFISH_KEYWORDS):
        return "zebrafish"
    if any(kw in o for kw in _PLANT_KEYWORDS):
        return "plant"
    if any(kw in o for kw in _YEAST_KEYWORDS):
        return "yeast"
    if any(kw in o for kw in _MAMMALIAN_KEYWORDS):
        return "mammalian"
    return "mammalian"  # default fallback for unlisted organisms

from modules.construction_file_tools.tools.create_construction_file import CreateConstructionFile
from modules.construction_file_tools.tools.validate_construction_file import ValidateConstructionFile
from modules.crispr_tools.tools.design_cas9_grna import DesignCas9Grna
from modules.crispr_tools.tools.cas_selector import CasSelector
from modules.crispr_tools.tools.design_cas12a_crrna import DesignCas12aCrrna
from modules.crispr_tools.tools.design_cloning_oligos import CRISPRCloningDesigner
from modules.crispr_tools.tools.fetch_target_sequence import FetchTargetSequence
from modules.crispr_tools.tools.rank_guides import rank_guides


_CONSTRUCTION_INPUT_FIELDS = {
    "construct_name",
    "assembly_strategy",
    "backbone_name",
    "backbone_sequence",
    "insert_name",
    "insert_sequence",
    "insert_forward_primer_name",
    "insert_forward_primer_sequence",
    "insert_reverse_primer_name",
    "insert_reverse_primer_sequence",
    "vector_forward_primer_name",
    "vector_forward_primer_sequence",
    "vector_reverse_primer_name",
    "vector_reverse_primer_sequence",
    "top_oligo_name",
    "top_oligo_sequence",
    "bottom_oligo_name",
    "bottom_oligo_sequence",
    "enzyme",
    "cell_strain",
    "selection",
    "temperature_c",
    "notes",
}


class RunFullCrisprWorkflow:
    """
    Description:
        End-to-end CRISPR pipeline: fetches target sequence, designs and ranks
        guides, selects a protospacer, designs cloning oligos, and builds and
        validates a cloning record.

    Input:
        query (str): Gene name, local resource key, or raw DNA sequence.
        organism (str): Target organism. Default: "Escherichia coli".
        vector (str): Cloning vector key (e.g. "px330"). Prompted if absent.
        construct_name (str): Optional name for the output construct.
        guide_index (int): Rank index of guide to use. Default: 0.
        validate_strict (bool): Strict validation mode. Default: False.
        force_vector (bool): Skip Cas-system compatibility check. Default: False.

    Output:
        dict: status="ready" with sequence_info, guides, selected_guide,
              cloning, and validation; or status="needs_user_input".

    Tests:
        - Case:
            Input: query="rpsL", organism="Escherichia coli", vector="pcrispr"
            Expected Output: status == "ready", "selected_guide" in result
            Description: Full workflow with a valid gene and vector.
        - Case:
            Input: query="rpsL", organism="Escherichia coli", vector=""
            Expected Output: status == "needs_user_input", missing_fields contains "vector"
            Description: Missing vector returns a structured prompt.
        - Case:
            Input: query=""
            Expected Exception: ValueError
            Description: Empty query raises ValueError.
    """

    def initiate(self) -> None:
        self.sequence_fetcher = FetchTargetSequence()
        self.sequence_fetcher.initiate()

        self.cas9_designer = DesignCas9Grna()
        self.cas9_designer.initiate()

        self.cas12a_designer = DesignCas12aCrrna()
        self.cas12a_designer.initiate()

        self.oligo_designer = CRISPRCloningDesigner()
        self.oligo_designer.initiate()

        self.construction_builder = CreateConstructionFile()
        self.construction_builder.initiate()

        self.construction_validator = ValidateConstructionFile()
        self.construction_validator.initiate()

        self.cas_selector = CasSelector()
        self.cas_selector.initiate()

    def run(
        self,
        query: str,
        organism: str = "Escherichia coli",
        vector: str = "",
        construct_name: Optional[str] = None,
        guide_index: int = 0,
        validate_strict: bool = False,
        force_vector: bool = False,
        confirmed: bool = False,
    ) -> dict:
        workflow_trace: list[dict] = []

        def add_trace(step: int, tool: str, inputs: dict, output_summary: dict) -> None:
            workflow_trace.append(
                {
                    "step": step,
                    "tool": tool,
                    "inputs": inputs,
                    "output_summary": output_summary,
                }
            )

        if not query or not query.strip():
            raise ValueError("query must not be empty.")
        if not vector or not vector.strip():
            category = _classify_organism(organism)
            recs = _VECTOR_RECOMMENDATIONS.get(category, _VECTOR_RECOMMENDATIONS["mammalian"])
            return {
                "status": "needs_user_input",
                "missing_fields": ["vector"],
                "organism": organism,
                "organism_category": category,
                "questions": [
                    f"Which vector would you like to use for {organism}? "
                    f"See vector_recommendations for literature-backed options."
                ],
                "vector_recommendations": recs,
                "workflow_trace": workflow_trace,
            }
        if not isinstance(guide_index, int) or guide_index < 0:
            raise ValueError("guide_index must be a non-negative integer.")

        sequence_info = self.sequence_fetcher.run(query=query, organism=organism)
        add_trace(
            1,
            "crispr_fetch_target_sequence",
            {"query": query, "organism": organism},
            {
                "source": sequence_info.get("source"),
                "resource": sequence_info.get("resource"),
                "organism": sequence_info.get("organism"),
                "length": sequence_info.get("length"),
                "ncbi_gene_id": sequence_info.get("ncbi_gene_id"),
                "ncbi_accession": sequence_info.get("ncbi_accession"),
            },
        )
        seq = sequence_info["sequence"]
        spec = self.oligo_designer.resolve_vector(vector)

        if not confirmed:
            method = spec.cloning_method if spec else "custom"
            enzyme = spec.enzyme if spec else "N/A"
            strain = spec.cell_strain if spec else "N/A"
            selection = spec.selection if spec else "N/A"
            vector_name = spec.name if spec else vector
            return {
                "status": "needs_user_input",
                "missing_fields": ["confirmed"],
                "questions": [
                    f"Please confirm you want to proceed with the following before running the full workflow:\n"
                    f"  Vector: {vector_name}\n"
                    f"  Cloning method: {method}\n"
                    f"  Enzyme: {enzyme}\n"
                    f"  Cell strain: {strain}\n"
                    f"  Selection: {selection}\n"
                    f"Reply with confirmed=true to continue, or choose a different vector."
                ],
                "vector": vector_name,
                "cloning_method": method,
                "enzyme": enzyme,
                "cell_strain": strain,
                "selection": selection,
                "workflow_trace": workflow_trace,
            }

        if spec is None:
            raise ValueError(
            "Custom vectors are not supported in the full workflow yet. "
            "Please use a known vector preset or provide enzyme/overhang handling."
        )

        cas_recommendation = self.cas_selector.run(
            seq=seq,
            repair_template=False,
        )
        add_trace(
            2,
            "crispr_cas_selector",
            {"sequence_length": len(seq), "repair_template": False},
            {
                "recommendation": cas_recommendation.get("recommendation"),
                "cas9_valid_guides": cas_recommendation.get("cas9_valid_guides"),
                "cas12a_valid_guides": cas_recommendation.get("cas12a_valid_guides"),
                "gc_content": cas_recommendation.get("gc_content"),
            },
        )

        def normalize_system(system: str) -> str:
            return {"SpCas9": "Cas9", "Cas9": "Cas9", "Cas12a": "Cas12a"}.get(system, system)

        recommended_system = normalize_system(cas_recommendation["recommendation"])
        vector_system = normalize_system(spec.nuclease_system)

        if recommended_system != vector_system and not force_vector:
            return {
                "status": "needs_user_input",
                "sequence_info": sequence_info,
                "cas_recommendation": cas_recommendation,
                "selected_vector": spec.name,
                "vector_system": vector_system,
                "recommended_system": recommended_system,
                "questions": [
                    f"The recommended system for this target is {recommended_system}, "
                    f"but the selected vector '{spec.name}' is meant for {vector_system}. "
                    "Do you still want to use this vector?"
                ],
                "options": [
                    f"Yes, continue with {spec.name} ({vector_system})",
                    f"No, switch to a {recommended_system}-compatible vector"
                ],
                "workflow_trace": workflow_trace,
            }
                    
        if spec.nuclease_system == "SpCas9":
            guides = self.cas9_designer.run(seq)
            guide_tool = "crispr_design_cas9_grna"
        elif spec.nuclease_system == "Cas12a":
            guides = self.cas12a_designer.run(seq)
            guide_tool = "crispr_design_cas12a_crrna"
        else:
            raise ValueError(f"Unsupported nuclease system: {spec.nuclease_system}")
        add_trace(
            3,
            guide_tool,
            {"sequence_length": len(seq)},
            {"guide_count": len(guides), "nuclease_system": spec.nuclease_system},
        )
        nuclease_key = "cas12a" if spec.nuclease_system == "Cas12a" else "cas9"
        ranking = rank_guides(guides=guides, reference=seq, nuclease=nuclease_key)
        ranked_guides    = ranking["ranked_guides"]
        scoring_rationale = ranking["scoring_rationale"]
        add_trace(
            4,
            "crispr_rank_guides",
            {"guide_count": len(guides), "reference_length": len(seq), "nuclease": nuclease_key},
            {
                "ranked_guide_count": len(ranked_guides),
                "best_guide": ranking.get("best_guide", {}).get("protospacer"),
                "best_score": ranking.get("best_guide", {}).get("total_score"),
            },
        )

        if guide_index >= len(ranked_guides):
            raise ValueError(
                f"guide_index {guide_index} is out of range for {len(ranked_guides)} designed guides."
            )

        selected_guide = ranked_guides[guide_index]
        protospacer = selected_guide["protospacer"]

        cloning = self.oligo_designer.run(
            vector=vector,
            protospacer=protospacer,
            construct_name=construct_name,
        )
        add_trace(
            5,
            "crispr_design_cloning_oligos",
            {"vector": vector, "protospacer": protospacer, "construct_name": construct_name},
            {
                "status": cloning.get("status"),
                "vector": cloning.get("vector"),
                "cloning_method": cloning.get("cloning_method"),
                "enzyme": cloning.get("enzyme"),
                "top_oligo": cloning.get("top_oligo"),
                "bottom_oligo": cloning.get("bottom_oligo"),
                "end_structure": cloning.get("end_structure"),
            },
        )
        if cloning.get("status") != "ready":
            return {
                "status": cloning.get("status", "needs_user_input"),
                "workflow_trace": workflow_trace,
                "sequence_info": sequence_info,
                "guides": ranked_guides,
                "selected_guide": selected_guide,
                "scoring_rationale": scoring_rationale,
                "protospacer": protospacer,
                "cloning": cloning,
            }

        construction_inputs = cloning["construction_file_inputs"]
        construction_payload = {
            key: value
            for key, value in construction_inputs.items()
            if key in _CONSTRUCTION_INPUT_FIELDS
        }
        construction = self.construction_builder.run(
            input_mode="sequence_build",
            **construction_payload,
        )
        add_trace(
            6,
            "create_construction_file",
            {"input_mode": "sequence_build", "assembly_strategy": construction_payload.get("assembly_strategy")},
            {
                "construct_name": construction.get("construct_name"),
                "assembly_strategy": construction.get("assembly_strategy"),
                "file_name": construction.get("file_name"),
            },
        )
        validation = self.construction_validator.run(
            strict=validate_strict,
            **construction_payload,
        )
        add_trace(
            7,
            "validate_construction_file",
            {"assembly_strategy": construction_payload.get("assembly_strategy"), "strict": validate_strict},
            {
                "summary": validation.get("summary"),
                "is_valid": validation.get("is_valid"),
                "detail_count": len(validation.get("details", [])),
            },
        )

        return {
            "status": "ready",
            "workflow_trace": workflow_trace,
            "sequence_info": sequence_info,
            "guides": ranked_guides,
            "selected_guide": selected_guide,
            "scoring_rationale": scoring_rationale,
            "protospacer": protospacer,
            "cloning": cloning,
            "construction_file_inputs": construction_inputs,
            "construction": construction,
            "validation": validation,
        }


_instance = RunFullCrisprWorkflow()
_instance.initiate()
run_full_crispr_workflow = _instance.run
