from __future__ import annotations

from platform import system
from typing import Optional

from modules.crispr_tools.tools.create_construction_file import CreateConstructionFile
from modules.crispr_tools.tools.construction_file_validation import ValidateConstructionFile
from modules.crispr_tools.tools.design_cas9_grna import DesignCas9Grna
from modules.crispr_tools.tools.crispr_cas_selector import CasSelector
from modules.crispr_tools.tools.design_cas12a_grna import DesignCas12aGrna
from modules.crispr_tools.tools.design_cloning_oligos import CRISPRCloningDesigner
from modules.crispr_tools.tools.fetch_target_sequence import FetchTargetSequence


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
    "enzyme",
    "cell_strain",
    "selection",
    "temperature_c",
    "notes",
}


class RunFullCrisprWorkflow:
    """
    Resolve a target sequence, design Cas9 guides, select one protospacer,
    design cloning oligos, build a construction file, and validate it.
    """

    def initiate(self) -> None:
        self.sequence_fetcher = FetchTargetSequence()
        self.sequence_fetcher.initiate()

        self.cas9_designer = DesignCas9Grna()
        self.cas9_designer.initiate()

        self.cas12a_designer = DesignCas12aGrna()
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
        force_vector : bool = False,
    ) -> dict:
        if not query or not query.strip():
            raise ValueError("query must not be empty.")
        if not vector or not vector.strip():
            raise ValueError("vector must not be empty.")
        if not isinstance(guide_index, int) or guide_index < 0:
            raise ValueError("guide_index must be a non-negative integer.")

        sequence_info = self.sequence_fetcher.run(query=query, organism=organism)
        seq = sequence_info["sequence"]
        spec = self.oligo_designer.resolve_vector(vector)

        if spec is None:
            raise ValueError(
            "Custom vectors are not supported in the full workflow yet. "
            "Please use a known vector preset or provide enzyme/overhang handling."
        )

        cas_recommendation = self.cas_selector.run(
            seq=seq,
            repair_template=False,
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
            }
                    
        if spec.nuclease_system == "SpCas9":
            guides = self.cas9_designer.run(seq)
        elif spec.nuclease_system == "Cas12a":
            guides = self.cas12a_designer.run(seq)
        else:
            raise ValueError(f"Unsupported nuclease system: {spec.nuclease_system}")
        if guide_index >= len(guides):
            raise ValueError(
                f"guide_index {guide_index} is out of range for {len(guides)} designed guides."
            )

        selected_guide = guides[guide_index]
        protospacer = selected_guide["protospacer"]

        cloning = self.oligo_designer.run(
            vector=vector,
            protospacer=protospacer,
            construct_name=construct_name,
        )
        if cloning.get("status") != "ready":
            return {
                "status": cloning.get("status", "needs_user_input"),
                "sequence_info": sequence_info,
                "guides": guides,
                "selected_guide": selected_guide,
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
        validation = self.construction_validator.run(
            strict=validate_strict,
            **construction_payload,
        )

        return {
            "status": "ready",
            "sequence_info": sequence_info,
            "guides": guides,
            "selected_guide": selected_guide,
            "protospacer": protospacer,
            "cloning": cloning,
            "construction_file_inputs": construction_inputs,
            "construction": construction,
            "validation": validation,
        }


_instance = RunFullCrisprWorkflow()
_instance.initiate()
run_full_crispr_workflow = _instance.run
