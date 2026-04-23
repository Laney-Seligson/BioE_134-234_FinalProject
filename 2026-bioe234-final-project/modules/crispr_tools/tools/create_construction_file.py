import json
from pathlib import Path
from typing import Any


class CreateConstructionFile:
    """
    Generate one of three outputs depending on input_mode:

    1) sequence_build
       Create a sequence-based construction file for Golden Gate, Gibson,
       or direct synthesis workflows.

    2) paper_info
       Create a compact paper-derived information resource. This stores only
       important workflow facts from a paper and does NOT create shorthand.

    3) paper_shorthand
       Create a shorthand construction/workflow file from paper-derived facts.
       This mode does not require DNA sequences and must not invent them.

    Notes:
    - sequence_build preserves the original sequence-centric behavior.
    - paper_info is a compact metadata/resource object only.
    - paper_shorthand is a workflow abstraction layer intended for literature.
    """

    def initiate(self) -> None:
        self.allowed_input_modes = {"sequence_build", "paper_info", "paper_shorthand"}
        self.allowed_strategies = {"GoldenGate", "Gibson", "DirectSynthesis"}
        self.allowed_paper_assembly_methods = {
            "GoldenGate",
            "Gibson",
            "Gateway",
            "DirectSynthesis",
            "Unknown",
        }
        self.allowed_part_types = {"oligo", "primer", "dsdna", "plasmid", "fragment"}
        self.allowed_step_types = {"PCR", "GoldenGate", "Gibson", "DirectSynthesis", "Transform"}
        self.allowed_shorthand_ops = {
            "clone",
            "assemble",
            "transform",
            "pcr",
            "sequence",
            "screen",
        }

    def run(
        self,
        input_mode: str = "sequence_build",

        # -------- sequence_build fields --------
        construct_name: str = "",
        assembly_strategy: str = "",
        backbone_name: str = "",
        backbone_sequence: str = "",
        insert_name: str = "",
        insert_sequence: str = "",
        insert_forward_primer_name: str = "",
        insert_forward_primer_sequence: str = "",
        insert_reverse_primer_name: str = "",
        insert_reverse_primer_sequence: str = "",
        vector_forward_primer_name: str = "",
        vector_forward_primer_sequence: str = "",
        vector_reverse_primer_name: str = "",
        vector_reverse_primer_sequence: str = "",
        enzyme: str = "",
        cell_strain: str = "",
        selection: str = "",
        temperature_c: int = 37,
        notes: str = "",

        # -------- paper_info / paper_shorthand shared fields --------
        paper_id: str = "",
        title: str = "",
        source_pdf: str = "",
        organism: str = "",
        system: str = "",
        targets: Any = None,
        vectors: Any = None,
        enzymes: Any = None,
        assembly_method: str = "",
        delivery_method: Any = None,
        validation_methods: Any = None,
        key_constraints: Any = None,
        paper_notes: Any = None,
        paper_info_json: Any = None,
    ) -> dict:

        input_mode = self._normalize_input_mode(input_mode)
        if input_mode not in self.allowed_input_modes:
            raise ValueError(f"input_mode must be one of {sorted(self.allowed_input_modes)}.")

        if input_mode == "paper_info":
            paper_resource = self._build_paper_info_resource(
                paper_id=paper_id,
                title=title,
                source_pdf=source_pdf,
                organism=organism,
                system=system,
                targets=targets,
                vectors=vectors,
                enzymes=enzymes,
                assembly_method=assembly_method,
                delivery_method=delivery_method,
                validation_methods=validation_methods,
                key_constraints=key_constraints,
                paper_notes=paper_notes,
                paper_info_json=paper_info_json,
            )
            rendered_text = self._render_paper_info_resource(paper_resource)
            file_stub = paper_resource["paper_id"] or "paper_info"
            return {
                "mode": "paper_info",
                "resource_type": paper_resource["resource_type"],
                "paper_id": paper_resource["paper_id"],
                "title": paper_resource["title"],
                "file_name": f"{file_stub}_important_info.json",
                "paper_important_info": paper_resource,
                "text": rendered_text,
            }

        if input_mode == "paper_shorthand":
            shorthand_result = self._build_paper_shorthand(
                paper_id=paper_id,
                title=title,
                source_pdf=source_pdf,
                organism=organism,
                system=system,
                targets=targets,
                vectors=vectors,
                enzymes=enzymes,
                assembly_method=assembly_method,
                delivery_method=delivery_method,
                validation_methods=validation_methods,
                key_constraints=key_constraints,
                paper_notes=paper_notes,
                paper_info_json=paper_info_json,
            )
            return shorthand_result

        # -------- sequence_build mode --------
        self._require_nonempty_string(construct_name, "construct_name")
        self._require_nonempty_string(assembly_strategy, "assembly_strategy")
        self._require_nonempty_string(backbone_name, "backbone_name")
        self._require_nonempty_string(backbone_sequence, "backbone_sequence")
        self._require_nonempty_string(insert_name, "insert_name")
        self._require_nonempty_string(insert_sequence, "insert_sequence")

        assembly_strategy = self._normalize_assembly_strategy(assembly_strategy)
        if assembly_strategy not in self.allowed_strategies:
            raise ValueError(f"assembly_strategy must be one of {sorted(self.allowed_strategies)}.")

        backbone_sequence = self._resolve_sequence_input(backbone_sequence, "backbone_sequence")
        insert_sequence = self._resolve_sequence_input(insert_sequence, "insert_sequence")

        self._validate_user_inputs(
            assembly_strategy=assembly_strategy,
            insert_forward_primer_name=insert_forward_primer_name,
            insert_forward_primer_sequence=insert_forward_primer_sequence,
            insert_reverse_primer_name=insert_reverse_primer_name,
            insert_reverse_primer_sequence=insert_reverse_primer_sequence,
            vector_forward_primer_name=vector_forward_primer_name,
            vector_forward_primer_sequence=vector_forward_primer_sequence,
            vector_reverse_primer_name=vector_reverse_primer_name,
            vector_reverse_primer_sequence=vector_reverse_primer_sequence,
            enzyme=enzyme,
        )

        parts = self._build_parts(
            backbone_name=backbone_name,
            backbone_sequence=backbone_sequence,
            insert_name=insert_name,
            insert_sequence=insert_sequence,
            insert_forward_primer_name=insert_forward_primer_name,
            insert_forward_primer_sequence=insert_forward_primer_sequence,
            insert_reverse_primer_name=insert_reverse_primer_name,
            insert_reverse_primer_sequence=insert_reverse_primer_sequence,
            vector_forward_primer_name=vector_forward_primer_name,
            vector_forward_primer_sequence=vector_forward_primer_sequence,
            vector_reverse_primer_name=vector_reverse_primer_name,
            vector_reverse_primer_sequence=vector_reverse_primer_sequence,
        )
        validated_parts = self._validate_parts(parts)

        operations = self._build_operations(
            construct_name=construct_name,
            assembly_strategy=assembly_strategy,
            backbone_name=backbone_name,
            insert_name=insert_name,
            insert_forward_primer_name=insert_forward_primer_name,
            insert_reverse_primer_name=insert_reverse_primer_name,
            vector_forward_primer_name=vector_forward_primer_name,
            vector_reverse_primer_name=vector_reverse_primer_name,
            enzyme=enzyme,
            cell_strain=cell_strain,
            selection=selection,
            temperature_c=temperature_c,
        )
        validated_operations = self._validate_operations(operations, validated_parts)
        construction_file_txt = self._render_construction_file(validated_parts, validated_operations)

        structured_construction_file = {
            "construct_name": construct_name,
            "assembly_strategy": assembly_strategy,
            "parts": validated_parts,
            "operations": validated_operations,
            "notes": notes,
        }

        return {
            "mode": "sequence_build",
            "construct_name": construct_name,
            "assembly_strategy": assembly_strategy,
            "file_name": f"{construct_name}_construction.txt",
            "structured_construction_file": structured_construction_file,
            "construction_file_txt": construction_file_txt,
            "text": construction_file_txt,
        }

    # ------------------------------------------------------------------
    # Shared utilities
    # ------------------------------------------------------------------

    def _normalize_input_mode(self, mode: str) -> str:
        if not isinstance(mode, str) or not mode.strip():
            raise ValueError("input_mode must be a non-empty string.")

        normalized = mode.strip().lower().replace("-", "_").replace(" ", "_")
        mapping = {
            "sequencebuild": "sequence_build",
            "sequence_build": "sequence_build",
            "paperinfo": "paper_info",
            "paper_info": "paper_info",
            "papershorthand": "paper_shorthand",
            "paper_shorthand": "paper_shorthand",
        }
        return mapping.get(normalized, mode.strip())

    def _require_nonempty_string(self, value: str, field_name: str) -> None:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{field_name} must be a non-empty string.")

    def _coerce_list(self, value: Any) -> list[str]:
        if value is None:
            return []
        if isinstance(value, list):
            items = value
        elif isinstance(value, str):
            if not value.strip():
                return []
            items = value.split(",")
        else:
            raise ValueError("List-like paper fields must be provided as a list or comma-separated string.")

        cleaned: list[str] = []
        for item in items:
            if isinstance(item, str):
                token = item.strip()
                if token:
                    cleaned.append(token)
            else:
                raise ValueError("All items in list-like paper fields must be strings.")
        return cleaned

    def _sanitize_name(self, text: str, fallback: str) -> str:
        if not isinstance(text, str):
            return fallback
        cleaned = []
        for char in text.strip():
            if char.isalnum() or char in {"_", "-"}:
                cleaned.append(char)
            elif char.isspace() or char in {"/", ":", ".", ",", ";", "(", ")"}:
                cleaned.append("_")
        out = "".join(cleaned).strip("_")
        return out or fallback

    def _normalize_sequence(self, sequence: str) -> str:
        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError("sequence must be a non-empty string.")

        raw = sequence.strip()

        cleaned = []
        for char in raw.upper():
            if char.isalpha():
                cleaned.append(char)
        cleaned = "".join(cleaned)

        if not cleaned:
            raise ValueError("sequence became empty after normalization.")

        invalid = set(cleaned) - set("ACGTN")
        if invalid:
            raise ValueError(
                "Sequence contains invalid DNA characters after resolution. "
                f"Invalid characters: {sorted(invalid)}"
            )
        return cleaned

    def _normalize_assembly_strategy(self, strategy: str) -> str:
        if not isinstance(strategy, str) or not strategy.strip():
            raise ValueError("assembly_strategy must be a non-empty string.")

        s = strategy.strip().lower().replace("_", "").replace(" ", "")
        mapping = {
            "goldengate": "GoldenGate",
            "gibson": "Gibson",
            "directsynthesis": "DirectSynthesis",
        }
        return mapping.get(s, strategy.strip())

    def _normalize_paper_assembly_method(self, value: str) -> str:
        if not isinstance(value, str) or not value.strip():
            return "Unknown"

        s = value.strip().lower().replace("_", "").replace(" ", "")
        mapping = {
            "goldengate": "GoldenGate",
            "gibson": "Gibson",
            "gateway": "Gateway",
            "directsynthesis": "DirectSynthesis",
            "unknown": "Unknown",
        }
        normalized = mapping.get(s, value.strip())
        if normalized not in self.allowed_paper_assembly_methods:
            return value.strip()
        return normalized

    def _paper_info_dir(self):
        from pathlib import Path
        return Path(__file__).resolve().parents[1] / "data" / "paper_info"

    def _normalize_paper_id(self, paper_id: str) -> str:
        value = paper_id.strip()
        prefix = "resource://paper_info/"
        if value.startswith(prefix):
            value = value[len(prefix):]
        if value.endswith(".json"):
            value = value[:-5]
        return value

    def _resolve_paper_info(
        self,
        paper_id: str,
        paper_info_json: Any,
    ) -> dict:
        if isinstance(paper_info_json, dict) and paper_info_json:
            data = dict(paper_info_json)
            if data.get("resource_type") == "paper_important_info_v1":
                return data

        if isinstance(paper_info_json, str) and paper_info_json.strip():
            try:
                data = json.loads(paper_info_json)
            except Exception:
                data = None
            if isinstance(data, dict) and data.get("resource_type") == "paper_important_info_v1":
                return data

        normalized_id = self._normalize_paper_id(paper_id) if isinstance(paper_id, str) and paper_id.strip() else ""
        if not normalized_id:
            return {}

        path = self._paper_info_dir() / f"{normalized_id}.json"
        if not path.exists():
            return {}

        try:
            import json
            data = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}

        if isinstance(data, dict) and data.get("resource_type") == "paper_important_info_v1":
            return data
        return {}


    # ------------------------------------------------------------------
    # sequence_build helpers
    # ------------------------------------------------------------------

    def _validate_user_inputs(
        self,
        assembly_strategy: str,
        insert_forward_primer_name: str,
        insert_forward_primer_sequence: str,
        insert_reverse_primer_name: str,
        insert_reverse_primer_sequence: str,
        vector_forward_primer_name: str,
        vector_forward_primer_sequence: str,
        vector_reverse_primer_name: str,
        vector_reverse_primer_sequence: str,
        enzyme: str,
    ) -> None:
        if assembly_strategy in {"GoldenGate", "Gibson"}:
            required_pairs = [
                ("insert_forward_primer", insert_forward_primer_name, insert_forward_primer_sequence),
                ("insert_reverse_primer", insert_reverse_primer_name, insert_reverse_primer_sequence),
                ("vector_forward_primer", vector_forward_primer_name, vector_forward_primer_sequence),
                ("vector_reverse_primer", vector_reverse_primer_name, vector_reverse_primer_sequence),
            ]

            missing_fields = []
            for label, primer_name, primer_seq in required_pairs:
                if not primer_name.strip():
                    missing_fields.append(f"{label}_name")
                if not primer_seq.strip():
                    missing_fields.append(f"{label}_sequence")

            if assembly_strategy == "GoldenGate" and not enzyme.strip():
                missing_fields.append("enzyme")

            if missing_fields:
                raise ValueError(
                    f"Missing required fields for {assembly_strategy}: {', '.join(missing_fields)}."
                )

        if assembly_strategy == "GoldenGate" and not enzyme.strip():
            raise ValueError("enzyme is required for GoldenGate workflows.")

    def _build_parts(
        self,
        backbone_name: str,
        backbone_sequence: str,
        insert_name: str,
        insert_sequence: str,
        insert_forward_primer_name: str,
        insert_forward_primer_sequence: str,
        insert_reverse_primer_name: str,
        insert_reverse_primer_sequence: str,
        vector_forward_primer_name: str,
        vector_forward_primer_sequence: str,
        vector_reverse_primer_name: str,
        vector_reverse_primer_sequence: str,
    ) -> list:
        parts = [
            {
                "part_type": "plasmid",
                "name": backbone_name,
                "sequence": self._normalize_sequence(backbone_sequence),
                "description": "Backbone plasmid",
            },
            {
                "part_type": "dsdna",
                "name": insert_name,
                "sequence": self._normalize_sequence(insert_sequence),
                "description": "Insert sequence",
            },
        ]

        primer_entries = [
            ("oligo", insert_forward_primer_name, insert_forward_primer_sequence, "Insert forward primer"),
            ("oligo", insert_reverse_primer_name, insert_reverse_primer_sequence, "Insert reverse primer"),
            ("oligo", vector_forward_primer_name, vector_forward_primer_sequence, "Vector forward primer"),
            ("oligo", vector_reverse_primer_name, vector_reverse_primer_sequence, "Vector reverse primer"),
        ]

        for part_type, name, sequence, description in primer_entries:
            if name.strip() and sequence.strip():
                parts.append(
                    {
                        "part_type": part_type,
                        "name": name,
                        "sequence": self._normalize_sequence(sequence),
                        "description": description,
                    }
                )

        return parts

    def _resolve_sequence_input(self, value: str, field_name: str) -> str:
        """
        Resolve a sequence input that may be a raw DNA string, a known resource
        name like 'pET28a', an MCP-style resource URI, or a local FASTA/GenBank file.
        Returns a raw DNA sequence string.
        """
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{field_name} must be a non-empty string.")

        raw = value.strip()

        # Already looks like DNA; let _normalize_sequence clean it later.
        dna_chars = set("ACGTNacgtn")
        compact = "".join(ch for ch in raw if ch.isalpha())
        if compact and set(compact) <= dna_chars:
            return raw

        # MCP-style resource URI.
        if raw.startswith("resource://"):
            resource_name = raw.split("/")[-1]
            resolved = self._load_sequence_resource_by_name(resource_name)
            if resolved:
                return resolved
            raise ValueError(
                f"{field_name} resource '{raw}' could not be resolved to a sequence."
            )

        # Plain resource names like pET28a, pBR322, pUC19.
        resolved = self._load_sequence_resource_by_name(raw)
        if resolved:
            return resolved

        # Local file path.
        maybe_path = Path(raw)
        if maybe_path.exists() and maybe_path.is_file():
            suffix = maybe_path.suffix.lower()
            if suffix in {".fa", ".fasta", ".fna"}:
                return self._read_fasta_file(maybe_path)
            if suffix in {".gb", ".gbk", ".genbank"}:
                return self._read_genbank_file(maybe_path)

        raise ValueError(
            f"{field_name} does not look like a valid DNA sequence and could not be "
            f"resolved as a known resource or file: '{raw}'"
        )

    def _load_sequence_resource_by_name(self, name: str) -> str | None:
        """Resolve known sequence resources from module data directories."""
        normalized = name.strip()

        candidate_paths = [
            Path("modules/crispr_tools/data") / f"{normalized}.gb",
            Path("modules/crispr_tools/data") / f"{normalized}.gbk",
            Path("modules/crispr_tools/data") / f"{normalized}.fasta",
            Path("modules/crispr_tools/data") / f"{normalized}.fa",
            Path("modules/seq_basics/data") / f"{normalized}.gb",
            Path("modules/seq_basics/data") / f"{normalized}.gbk",
            Path("modules/seq_basics/data") / f"{normalized}.fasta",
            Path("modules/seq_basics/data") / f"{normalized}.fa",
        ]

        for path in candidate_paths:
            if path.exists():
                suffix = path.suffix.lower()
                if suffix in {".fa", ".fasta", ".fna"}:
                    return self._read_fasta_file(path)
                if suffix in {".gb", ".gbk", ".genbank"}:
                    return self._read_genbank_file(path)
        return None

    def _read_fasta_file(self, path: Path) -> str:
        lines = path.read_text(encoding="utf-8").splitlines()
        seq_lines = [line.strip() for line in lines if line.strip() and not line.startswith(">")]
        seq = "".join(seq_lines)
        return self._normalize_sequence(seq)

    def _read_genbank_file(self, path: Path) -> str:
        text = path.read_text(encoding="utf-8")
        if "ORIGIN" not in text:
            raise ValueError(f"GenBank file '{path}' is missing ORIGIN section.")

        origin = text.split("ORIGIN", 1)[1]
        origin = origin.split("//", 1)[0]
        letters_only = "".join(ch for ch in origin if ch.isalpha())
        return self._normalize_sequence(letters_only)

    def _build_operations(
        self,
        construct_name: str,
        assembly_strategy: str,
        backbone_name: str,
        insert_name: str,
        insert_forward_primer_name: str,
        insert_reverse_primer_name: str,
        vector_forward_primer_name: str,
        vector_reverse_primer_name: str,
        enzyme: str,
        cell_strain: str,
        selection: str,
        temperature_c: int,
    ) -> list:
        operations = []

        if assembly_strategy in {"GoldenGate", "Gibson"}:
            insert_pcr_product = f"{insert_name}_pcr"
            vector_pcr_product = f"{backbone_name}_pcr"

            operations.append(
                {
                    "step_number": 1,
                    "step_type": "PCR",
                    "inputs": [insert_forward_primer_name, insert_reverse_primer_name, insert_name],
                    "parameters": {
                        "forward_primer": insert_forward_primer_name,
                        "reverse_primer": insert_reverse_primer_name,
                        "template": insert_name,
                    },
                    "output": insert_pcr_product,
                }
            )

            operations.append(
                {
                    "step_number": 2,
                    "step_type": "PCR",
                    "inputs": [vector_forward_primer_name, vector_reverse_primer_name, backbone_name],
                    "parameters": {
                        "forward_primer": vector_forward_primer_name,
                        "reverse_primer": vector_reverse_primer_name,
                        "template": backbone_name,
                    },
                    "output": vector_pcr_product,
                }
            )

            if assembly_strategy == "GoldenGate":
                operations.append(
                    {
                        "step_number": 3,
                        "step_type": "GoldenGate",
                        "inputs": [vector_pcr_product, insert_pcr_product],
                        "parameters": {"enzyme": enzyme},
                        "output": construct_name,
                    }
                )
            else:
                operations.append(
                    {
                        "step_number": 3,
                        "step_type": "Gibson",
                        "inputs": [vector_pcr_product, insert_pcr_product],
                        "parameters": {"reagent": "GibsonMix", "overlap_bp": 20},
                        "output": construct_name,
                    }
                )

        elif assembly_strategy == "DirectSynthesis":
            operations.append(
                {
                    "step_number": 1,
                    "step_type": "DirectSynthesis",
                    "inputs": [insert_name],
                    "parameters": {},
                    "output": construct_name,
                }
            )

        if cell_strain.strip() and selection.strip():
            operations.append(
                {
                    "step_number": len(operations) + 1,
                    "step_type": "Transform",
                    "inputs": [construct_name],
                    "parameters": {
                        "cells": cell_strain,
                        "selection": selection,
                        "temperature_c": temperature_c,
                    },
                    "output": f"{construct_name}_e",
                }
            )

        return operations

    def _normalize_sequence(self, sequence: str) -> str:
        if not isinstance(sequence, str) or not sequence.strip():
            raise ValueError("sequence must be a non-empty string.")

        raw = sequence.strip()

        # Reject obvious unresolved placeholders / resource names
        # instead of silently converting them into fake DNA.
        if raw.startswith("resource://"):
            raise ValueError(
                f"sequence '{raw}' was not resolved before normalization."
            )

        cleaned = []
        for char in raw.upper():
            if char.isalpha():
                cleaned.append(char)

        cleaned = "".join(cleaned)

        if not cleaned:
            raise ValueError("sequence became empty after normalization.")

        # Auto-convert RNA to DNA (U → T) — construction files are always DNA
        cleaned = cleaned.replace("U", "T")

        invalid = set(cleaned) - set("ACGTN")
        if invalid:
            raise ValueError(
                "sequence contains invalid DNA characters or appears to be an unresolved "
                f"resource/placeholder: {sorted(invalid)}"
            )

        return cleaned

    def _validate_parts(self, parts: list) -> list:
        validated = []
        seen_names = set()

        for part in parts:
            if not isinstance(part, dict):
                raise ValueError("Each part must be a dictionary.")

            for field in ("part_type", "name", "sequence"):
                if field not in part:
                    raise ValueError(f"Each part must include '{field}'.")

            part_type = part["part_type"]
            name = part["name"]
            sequence = part["sequence"]
            description = part.get("description", "")

            if part_type not in self.allowed_part_types:
                raise ValueError(f"Invalid part_type '{part_type}'.")
            if name in seen_names:
                raise ValueError(f"Duplicate part name '{name}'.")
            seen_names.add(name)

            validated.append(
                {
                    "part_type": part_type,
                    "name": name,
                    "sequence": self._normalize_sequence(sequence),
                    "description": description,
                }
            )

        return validated

    def _validate_step_specific_fields(
        self,
        step_number: int,
        step_type: str,
        inputs: list,
        parameters: dict,
    ) -> None:
        if step_type == "PCR":
            required = {"forward_primer", "reverse_primer", "template"}
            missing = [field for field in required if field not in parameters]
            if missing:
                raise ValueError(f"PCR step {step_number} missing required parameter(s): {missing}.")

        elif step_type == "GoldenGate":
            if len(inputs) < 2:
                raise ValueError(f"GoldenGate step {step_number} requires at least 2 inputs.")
            if "enzyme" not in parameters:
                raise ValueError(f"GoldenGate step {step_number} requires 'enzyme'.")

        elif step_type == "Gibson":
            if len(inputs) < 2:
                raise ValueError(f"Gibson step {step_number} requires at least 2 inputs.")
            if "overlap_bp" not in parameters and "overlap_notes" not in parameters:
                raise ValueError(
                    f"Gibson step {step_number} requires 'overlap_bp' or 'overlap_notes'."
                )

        elif step_type == "DirectSynthesis":
            if len(inputs) != 1:
                raise ValueError(f"DirectSynthesis step {step_number} should have exactly 1 input.")

        elif step_type == "Transform":
            if len(inputs) != 1:
                raise ValueError(f"Transform step {step_number} should have exactly 1 input.")
            required = {"cells", "selection"}
            missing = [field for field in required if field not in parameters]
            if missing:
                raise ValueError(
                    f"Transform step {step_number} missing required parameter(s): {missing}."
                )

    def _validate_operations(self, operations: list, parts: list) -> list:
        validated = []
        part_names = {part["name"] for part in parts}
        produced_outputs = set()
        seen_steps = set()

        sorted_ops = sorted(operations, key=lambda op: op.get("step_number", 0))

        for op in sorted_ops:
            if not isinstance(op, dict):
                raise ValueError("Each operation must be a dictionary.")

            for field in ("step_number", "step_type", "inputs", "output"):
                if field not in op:
                    raise ValueError(f"Each operation must include '{field}'.")

            step_number = op["step_number"]
            step_type = op["step_type"]
            inputs = op["inputs"]
            output = op["output"]
            parameters = op.get("parameters", {})
            description = op.get("description", "")

            if not isinstance(step_number, int):
                raise ValueError("step_number must be an integer.")
            if step_number in seen_steps:
                raise ValueError(f"Duplicate step_number '{step_number}'.")
            seen_steps.add(step_number)

            if step_type not in self.allowed_step_types:
                raise ValueError(f"Invalid step_type '{step_type}'.")
            if not isinstance(inputs, list):
                raise ValueError("inputs must be a list.")
            if not isinstance(output, str) or not output.strip():
                raise ValueError("output must be a non-empty string.")
            if not isinstance(parameters, dict):
                raise ValueError("parameters must be a dictionary.")

            available_names = part_names | produced_outputs
            for item in inputs:
                if item not in available_names:
                    raise ValueError(
                        f"Operation step {step_number} references unknown input '{item}'."
                    )

            self._validate_step_specific_fields(step_number, step_type, inputs, parameters)

            validated.append(
                {
                    "step_number": step_number,
                    "step_type": step_type,
                    "inputs": inputs,
                    "parameters": parameters,
                    "output": output,
                    "description": description,
                }
            )
            produced_outputs.add(output)

        return validated

    def _render_construction_file(self, parts: list, operations: list) -> str:
        lines = []

        def _shorten(seq: str, max_len: int = 40) -> str:
            if len(seq) <= max_len:
                return seq
            return f"{seq[:20]}...{seq[-20:]} ({len(seq)} bp)"

        for op in operations:
            step_type = op["step_type"]
            inputs = op["inputs"]
            parameters = op["parameters"]
            output = op["output"]

            if step_type == "PCR":
                lines.append(
                    f"{'PCR':<14}"
                    f"{parameters.get('forward_primer', ''):<14}"
                    f"{parameters.get('reverse_primer', ''):<14}"
                    f"{parameters.get('template', ''):<18}"
                    f"{output}"
                )
            elif step_type == "GoldenGate":
                lines.append(
                    f"{'GoldenGate':<14}"
                    f"{inputs[0] if len(inputs) > 0 else '':<14}"
                    f"{inputs[1] if len(inputs) > 1 else '':<14}"
                    f"{parameters.get('enzyme', ''):<18}"
                    f"{output}"
                )
            elif step_type == "Gibson":
                reagent = parameters.get("reagent", "GibsonMix")
                lines.append(
                    f"{'Gibson':<14}"
                    f"{inputs[0] if len(inputs) > 0 else '':<14}"
                    f"{inputs[1] if len(inputs) > 1 else '':<14}"
                    f"{reagent:<18}"
                    f"{output}"
                )
            elif step_type == "DirectSynthesis":
                lines.append(f"{'DirectSynthesis':<18}{inputs[0] if inputs else '':<16}{output}")
            elif step_type == "Transform":
                lines.append(
                    f"{'Transform':<14}"
                    f"{inputs[0] if inputs else '':<22}"
                    f"{parameters.get('cells', ''):<12}"
                    f"{parameters.get('selection', ''):<10}"
                    f"{str(parameters.get('temperature_c', '')):<8}"
                    f"{output}"
                )

        lines.append("")

        for part in parts:
            seq_display = part["sequence"] if part["part_type"] == "oligo" else _shorten(part["sequence"])
            lines.append(f"{part['part_type']:<14}{part['name']:<14}{seq_display}")

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # paper_info helpers
    # ------------------------------------------------------------------

    def _build_paper_info_resource(
        self,
        paper_id: str,
        title: str,
        source_pdf: str,
        organism: str,
        system: str,
        targets: Any,
        vectors: Any,
        enzymes: Any,
        assembly_method: str,
        delivery_method: Any,
        validation_methods: Any,
        key_constraints: Any,
        paper_notes: Any,
        paper_info_json: Any = None,
    ) -> dict:
        self._require_nonempty_string(paper_id, "paper_id")

        resolved_info = self._resolve_paper_info(paper_id, paper_info_json)
        if resolved_info:
            title = title or resolved_info.get("title", "")
            source_pdf = source_pdf or resolved_info.get("source_pdf", "")
            organism = organism or resolved_info.get("organism", "")
            system = system or resolved_info.get("system", "")
            targets = targets if self._coerce_list(targets) else resolved_info.get("targets", [])
            vectors = vectors if self._coerce_list(vectors) else resolved_info.get("vectors", [])
            enzymes = enzymes if self._coerce_list(enzymes) else resolved_info.get("enzymes", [])
            assembly_method = assembly_method or resolved_info.get("assembly_method", "")
            delivery_method = delivery_method if self._coerce_list(delivery_method) else resolved_info.get("delivery_method", [])
            validation_methods = validation_methods if self._coerce_list(validation_methods) else resolved_info.get("validation_methods", [])
            key_constraints = key_constraints if self._coerce_list(key_constraints) else resolved_info.get("key_constraints", [])
            paper_notes = paper_notes if self._coerce_list(paper_notes) else resolved_info.get("notes", [])

        title = title.strip() if isinstance(title, str) else ""
        if not title:
            title = self._normalize_paper_id(paper_id)

        resource = {
            "resource_type": "paper_important_info_v1",
            "paper_id": paper_id.strip(),
            "title": title.strip(),
            "source_pdf": source_pdf.strip() if isinstance(source_pdf, str) else "",
            "organism": organism.strip() if isinstance(organism, str) else "",
            "system": system.strip() if isinstance(system, str) else "",
            "targets": self._coerce_list(targets),
            "vectors": self._coerce_list(vectors),
            "enzymes": self._coerce_list(enzymes),
            "assembly_method": self._normalize_paper_assembly_method(assembly_method),
            "delivery_method": self._coerce_list(delivery_method),
            "validation_methods": self._coerce_list(validation_methods),
            "key_constraints": self._coerce_list(key_constraints),
            "notes": self._coerce_list(paper_notes),
        }
        return resource

    def _render_paper_info_resource(self, resource: dict) -> str:
        return json.dumps(resource, indent=2)

    # ------------------------------------------------------------------
    # paper_shorthand helpers
    # ------------------------------------------------------------------

    def _build_paper_shorthand(
        self,
        paper_id: str,
        title: str = "",
        source_pdf: str = "",
        organism: str = "",
        system: str = "",
        targets: Any = None,
        vectors: Any = None,
        enzymes: Any = None,
        assembly_method: str = "",
        delivery_method: Any = None,
        validation_methods: Any = None,
        key_constraints: Any = None,
        paper_notes: Any = None,
        paper_info_json: Any = None,
    ) -> dict:
        self._require_nonempty_string(paper_id, "paper_id")

        resolved_info = self._resolve_paper_info(paper_id, paper_info_json)
        if resolved_info:
            title = title or resolved_info.get("title", "")
            source_pdf = source_pdf or resolved_info.get("source_pdf", "")
            organism = organism or resolved_info.get("organism", "")
            system = system or resolved_info.get("system", "")
            targets = targets if self._coerce_list(targets) else resolved_info.get("targets", [])
            vectors = vectors if self._coerce_list(vectors) else resolved_info.get("vectors", [])
            enzymes = enzymes if self._coerce_list(enzymes) else resolved_info.get("enzymes", [])
            assembly_method = assembly_method or resolved_info.get("assembly_method", "")
            delivery_method = delivery_method if self._coerce_list(delivery_method) else resolved_info.get("delivery_method", [])
            validation_methods = validation_methods if self._coerce_list(validation_methods) else resolved_info.get("validation_methods", [])
            key_constraints = key_constraints if self._coerce_list(key_constraints) else resolved_info.get("key_constraints", [])
            paper_notes = paper_notes if self._coerce_list(paper_notes) else resolved_info.get("notes", [])

        title = title.strip() if isinstance(title, str) else ""
        if not title:
            title = self._normalize_paper_id(paper_id)

        organism = organism.strip() if isinstance(organism, str) else ""
        system = system.strip() if isinstance(system, str) else ""
        source_pdf = source_pdf.strip() if isinstance(source_pdf, str) else ""
        targets_list = self._coerce_list(targets)
        vectors_list = self._coerce_list(vectors)
        enzymes_list = self._coerce_list(enzymes)
        delivery_list = self._coerce_list(delivery_method)
        validation_list = self._coerce_list(validation_methods)
        constraints_list = self._coerce_list(key_constraints)
        notes_list = self._coerce_list(paper_notes)
        assembly_method = self._normalize_paper_assembly_method(assembly_method)

        declarations = self._infer_shorthand_declarations(vectors_list, enzymes_list, organism, delivery_list)
        steps = self._infer_shorthand_steps(
            system=system,
            organism=organism,
            targets=targets_list,
            vectors=vectors_list,
            enzymes=enzymes_list,
            assembly_method=assembly_method,
            delivery_method=delivery_list,
            validation_methods=validation_list,
        )

        shorthand_txt = self._render_paper_shorthand(
            title=title,
            paper_id=paper_id,
            source_pdf=source_pdf,
            organism=organism,
            system=system,
            declarations=declarations,
            steps=steps,
            constraints=constraints_list,
            notes=notes_list,
        )

        return {
            "mode": "paper_shorthand",
            "resource_type": "paper_shorthand_v1",
            "paper_id": paper_id,
            "title": title,
            "source_pdf": source_pdf,
            "organism": organism,
            "system": system,
            "assembly_method": assembly_method,
            "declarations": declarations,
            "steps": steps,
            "key_constraints": constraints_list,
            "notes": notes_list,
            "file_name": f"{self._sanitize_name(paper_id, 'paper')}_shorthand.txt",
            "paper_shorthand_txt": shorthand_txt,
            "text": shorthand_txt,
        }

    def _infer_shorthand_declarations(
        self,
        vectors: list[str],
        enzymes: list[str],
        organism: str,
        delivery_method: list[str],
    ) -> list[str]:
        declarations: list[str] = []
        seen = set()

        for vector in vectors:
            line = f"plasmid {self._sanitize_name(vector, 'vector')}"
            if line not in seen:
                declarations.append(line)
                seen.add(line)

        for enzyme in enzymes:
            line = f"enzyme {self._sanitize_name(enzyme, 'enzyme')}"
            if line not in seen:
                declarations.append(line)
                seen.add(line)

        org_name = self._sanitize_name(organism, "target_cells") if organism else ""
        if org_name:
            line = f"cell {org_name}"
            if line not in seen:
                declarations.append(line)
                seen.add(line)

        if any("agrobacterium" in method.lower() for method in delivery_method):
            line = "host Agrobacterium"
            if line not in seen:
                declarations.append(line)
                seen.add(line)

        return declarations

    def _infer_shorthand_steps(
        self,
        system: str,
        organism: str,
        targets: list[str],
        vectors: list[str],
        enzymes: list[str],
        assembly_method: str,
        delivery_method: list[str],
        validation_methods: list[str],
    ) -> list[str]:
        steps: list[str] = []

        sanitized_vectors = [self._sanitize_name(v, "vector") for v in vectors]
        sanitized_targets = [self._sanitize_name(t, "target") for t in targets]
        sanitized_enzymes = [self._sanitize_name(e, "enzyme") for e in enzymes]
        sanitized_organism = self._sanitize_name(organism, "target_cells") if organism else "target_cells"

        cas_vector = None
        guide_vector = None
        for vector in sanitized_vectors:
            lower = vector.lower()
            if cas_vector is None and "cas" in lower:
                cas_vector = vector
            if guide_vector is None and ("grna" in lower or "sgrna" in lower or "tracr" in lower or "cr" in lower):
                guide_vector = vector

        if cas_vector is None and sanitized_vectors:
            cas_vector = sanitized_vectors[0]
        if guide_vector is None and len(sanitized_vectors) > 1:
            guide_vector = sanitized_vectors[1]
        elif guide_vector is None and sanitized_vectors:
            guide_vector = sanitized_vectors[0]

        chosen_enzyme = sanitized_enzymes[0] if sanitized_enzymes else "enzyme"

        is_crispr = "crispr" in system.lower() if system else False
        uses_agro = any("agrobacterium" in method.lower() for method in delivery_method)
        uses_biolistic = any(
            "bombard" in method.lower() or "biolistic" in method.lower() or "particle" in method.lower()
            for method in delivery_method
        )

        if is_crispr and guide_vector:
            steps.append(f"clone spacerPair {guide_vector} {chosen_enzyme} grna_entry")
            if cas_vector:
                steps.append(f"assemble grna_entry {cas_vector} {assembly_method} crispr_plasmid")
            else:
                steps.append(f"assemble grna_entry backbone {assembly_method} crispr_plasmid")
        elif len(sanitized_vectors) >= 2:
            steps.append(f"assemble {sanitized_vectors[1]} {sanitized_vectors[0]} {assembly_method} assembled_construct")
        elif sanitized_vectors:
            steps.append(f"assemble insert {sanitized_vectors[0]} {assembly_method} assembled_construct")
        else:
            raise ValueError("paper_shorthand requires at least one vector or structured paper info with vectors.")

        if uses_agro:
            steps.append("transform crispr_plasmid Agrobacterium electroporation agro_strain")
            steps.append(f"transform agro_strain {sanitized_organism} infection edited_cells")
        elif uses_biolistic:
            steps.append(f"transform crispr_plasmid {sanitized_organism} bombardment edited_cells")
        elif delivery_method:
            generic_method = self._sanitize_name(delivery_method[0], "delivery")
            steps.append(f"transform crispr_plasmid {sanitized_organism} {generic_method} edited_cells")

        if sanitized_targets:
            first_target = sanitized_targets[0]
            if any("pcr" in method.lower() for method in validation_methods):
                steps.append(f"pcr {first_target}_F {first_target}_R {first_target} amplicon")
                amplicon_name = "amplicon"
            else:
                amplicon_name = first_target
        else:
            amplicon_name = "amplicon"

        if any("sanger" in method.lower() or "sequencing" in method.lower() or "sequence" in method.lower() for method in validation_methods):
            steps.append(f"sequence {amplicon_name} Sanger mutation_call")

        if any("phenotype" in method.lower() or "reporter" in method.lower() or "gus" in method.lower() for method in validation_methods):
            steps.append("screen edited_cells phenotype validated_line")

        return steps

    def _render_paper_shorthand(
        self,
        title: str,
        paper_id: str,
        source_pdf: str,
        organism: str,
        system: str,
        declarations: list[str],
        steps: list[str],
        constraints: list[str],
        notes: list[str],
    ) -> str:
        lines = [
            f"# Paper shorthand for: {title}",
            f"# paper_id: {paper_id}",
        ]
        if source_pdf:
            lines.append(f"# source_pdf: {source_pdf}")
        if organism:
            lines.append(f"# organism: {organism}")
        if system:
            lines.append(f"# system: {system}")
        for constraint in constraints:
            lines.append(f"# constraint: {constraint}")
        for note in notes:
            lines.append(f"# note: {note}")

        if declarations:
            lines.append("")
            lines.extend(declarations)

        if steps:
            lines.append("")
            lines.extend(steps)

        return "\n".join(lines)


_instance = CreateConstructionFile()
_instance.initiate()
create_construction_file = _instance.run


def prompt_optional(prompt_text: str) -> str:
    return input(prompt_text).strip()


def prompt_required(prompt_text: str) -> str:
    value = input(prompt_text).strip()
    while not value:
        print("This field is required.")
        value = input(prompt_text).strip()
    return value


def prompt_int(prompt_text: str, default: int = 37) -> int:
    value = input(f"{prompt_text} [{default}]: ").strip()
    if not value:
        return default
    while True:
        try:
            return int(value)
        except ValueError:
            value = input("Please enter an integer: ").strip()


def prompt_csv(prompt_text: str) -> str:
    return input(prompt_text).strip()


def main() -> None:
    print("=== Construction File Generator ===")
    print("Modes: sequence_build, paper_info, paper_shorthand")
    print("Leave optional fields blank if they are not needed.\n")

    input_mode = prompt_required("Input mode (sequence_build, paper_info, paper_shorthand): ")
    normalized_mode = input_mode.strip().lower().replace("-", "_").replace(" ", "_")

    try:
        if normalized_mode == "sequence_build":
            construct_name = prompt_required("Construct name: ")
            assembly_strategy = prompt_required(
                "Assembly strategy (GoldenGate, Gibson, DirectSynthesis): "
            )
            backbone_name = prompt_required("Backbone name: ")
            backbone_sequence = prompt_required("Backbone sequence: ")
            insert_name = prompt_required("Insert name: ")
            insert_sequence = prompt_required("Insert sequence: ")

            insert_forward_primer_name = ""
            insert_forward_primer_sequence = ""
            insert_reverse_primer_name = ""
            insert_reverse_primer_sequence = ""
            vector_forward_primer_name = ""
            vector_forward_primer_sequence = ""
            vector_reverse_primer_name = ""
            vector_reverse_primer_sequence = ""
            enzyme = ""

            normalized_strategy = assembly_strategy.strip().lower().replace("_", "").replace(" ", "")
            if normalized_strategy in {"goldengate", "gibson"}:
                print("\n--- Primer information required for this strategy ---")
                insert_forward_primer_name = prompt_required("Insert forward primer name: ")
                insert_forward_primer_sequence = prompt_required("Insert forward primer sequence: ")
                insert_reverse_primer_name = prompt_required("Insert reverse primer name: ")
                insert_reverse_primer_sequence = prompt_required("Insert reverse primer sequence: ")
                vector_forward_primer_name = prompt_required("Vector forward primer name: ")
                vector_forward_primer_sequence = prompt_required("Vector forward primer sequence: ")
                vector_reverse_primer_name = prompt_required("Vector reverse primer name: ")
                vector_reverse_primer_sequence = prompt_required("Vector reverse primer sequence: ")

            if normalized_strategy == "goldengate":
                enzyme = prompt_required("Assembly enzyme: ")

            print("\n--- Optional transformation information ---")
            cell_strain = prompt_optional("Cell strain: ")
            selection = prompt_optional("Selection marker: ")
            temperature_c = prompt_int("Transformation temperature (C)", default=37)
            notes = prompt_optional("Notes: ")

            result = create_construction_file(
                input_mode="sequence_build",
                construct_name=construct_name,
                assembly_strategy=assembly_strategy,
                backbone_name=backbone_name,
                backbone_sequence=backbone_sequence,
                insert_name=insert_name,
                insert_sequence=insert_sequence,
                insert_forward_primer_name=insert_forward_primer_name,
                insert_forward_primer_sequence=insert_forward_primer_sequence,
                insert_reverse_primer_name=insert_reverse_primer_name,
                insert_reverse_primer_sequence=insert_reverse_primer_sequence,
                vector_forward_primer_name=vector_forward_primer_name,
                vector_forward_primer_sequence=vector_forward_primer_sequence,
                vector_reverse_primer_name=vector_reverse_primer_name,
                vector_reverse_primer_sequence=vector_reverse_primer_sequence,
                enzyme=enzyme,
                cell_strain=cell_strain,
                selection=selection,
                temperature_c=temperature_c,
                notes=notes,
            )

            print("\n=== Construction file generated ===")
            print(f"File name: {result['file_name']}\n")
            print(result["construction_file_txt"])

            save_choice = input("\nSave to file? (y/n): ").strip().lower()
            if save_choice == "y":
                with open(result["file_name"], "w") as f:
                    f.write(result["construction_file_txt"])
                print(f"Saved to {result['file_name']}")

        elif normalized_mode == "paper_info":
            print("\n--- Paper important information mode ---")
            paper_id = prompt_required("Paper ID: ")
            title = prompt_optional("Title (optional if paper_info JSON exists): ")
            source_pdf = prompt_optional("Source PDF: ")
            organism = prompt_optional("Organism: ")
            system = prompt_optional("System: ")
            targets = prompt_csv("Targets (comma-separated): ")
            vectors = prompt_csv("Vectors (comma-separated): ")
            enzymes = prompt_csv("Enzymes (comma-separated): ")
            assembly_method = prompt_optional("Assembly method: ")
            delivery_method = prompt_csv("Delivery method(s) (comma-separated): ")
            validation_methods = prompt_csv("Validation method(s) (comma-separated): ")
            key_constraints = prompt_csv("Key constraints (comma-separated): ")
            paper_notes = prompt_csv("Notes (comma-separated): ")

            result = create_construction_file(
                input_mode="paper_info",
                paper_id=paper_id,
                title=title,
                source_pdf=source_pdf,
                organism=organism,
                system=system,
                targets=targets,
                vectors=vectors,
                enzymes=enzymes,
                assembly_method=assembly_method,
                delivery_method=delivery_method,
                validation_methods=validation_methods,
                key_constraints=key_constraints,
                paper_notes=paper_notes,
            )

            print("\n=== Paper information resource generated ===\n")
            print(result["text"])

            save_choice = input("\nSave paper info to JSON file? (y/n): ").strip().lower()
            if save_choice == "y":
                with open(result["file_name"], "w") as f:
                    json.dump(result["paper_important_info"], f, indent=2)
                print(f"Saved to {result['file_name']}")

        elif normalized_mode == "paper_shorthand":
            print("\n--- Paper shorthand mode ---")
            paper_id = prompt_required("Paper ID: ")
            title = prompt_optional("Title (optional if paper_info JSON exists): ")
            source_pdf = prompt_optional("Source PDF: ")
            organism = prompt_optional("Organism: ")
            system = prompt_optional("System: ")
            targets = prompt_csv("Targets (comma-separated): ")
            vectors = prompt_csv("Vectors (comma-separated): ")
            enzymes = prompt_csv("Enzymes (comma-separated): ")
            assembly_method = prompt_optional("Assembly method: ")
            delivery_method = prompt_csv("Delivery method(s) (comma-separated): ")
            validation_methods = prompt_csv("Validation method(s) (comma-separated): ")
            key_constraints = prompt_csv("Key constraints (comma-separated): ")
            paper_notes = prompt_csv("Notes (comma-separated): ")

            result = create_construction_file(
                input_mode="paper_shorthand",
                paper_id=paper_id,
                title=title,
                source_pdf=source_pdf,
                organism=organism,
                system=system,
                targets=targets,
                vectors=vectors,
                enzymes=enzymes,
                assembly_method=assembly_method,
                delivery_method=delivery_method,
                validation_methods=validation_methods,
                key_constraints=key_constraints,
                paper_notes=paper_notes,
            )

            print("\n=== Paper shorthand generated ===\n")
            print(result["paper_shorthand_txt"])

            save_choice = input("\nSave shorthand to file? (y/n): ").strip().lower()
            if save_choice == "y":
                with open(result["file_name"], "w") as f:
                    f.write(result["paper_shorthand_txt"])
                print(f"Saved to {result['file_name']}")

        else:
            raise ValueError("Unknown input mode.")

    except ValueError as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main()
