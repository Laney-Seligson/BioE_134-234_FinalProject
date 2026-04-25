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
        self.allowed_strategies = {"GoldenGate", "Gibson", "DirectSynthesis", "TypeIISOligoCloning", "RestrictionLigation"}
        self.allowed_paper_assembly_methods = {
            "GoldenGate",
            "Gibson",
            "Gateway",
            "DirectSynthesis",
            "Unknown",
        }
        self.allowed_part_types = {"oligo", "primer", "dsdna", "plasmid", "fragment"}
        self.allowed_step_types = {"PCR", "GoldenGate", "Gibson", "DirectSynthesis", "TypeIISOligoCloning", "Transform", "RestrictionLigation"}
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
        top_oligo_name: str = "",
        top_oligo_sequence: str = "",
        bottom_oligo_name: str = "",
        bottom_oligo_sequence: str = "",
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

        assembly_strategy = self._normalize_assembly_strategy(assembly_strategy)
        if assembly_strategy not in self.allowed_strategies:
            raise ValueError(f"assembly_strategy must be one of {sorted(self.allowed_strategies)}.")

        if assembly_strategy == "TypeIISOligoCloning":
            # This workflow clones annealed top/bottom oligos directly into a
            # Type IIS-digested backbone, so it does not require an insert PCR
            # template or insert_sequence.
            if not insert_name.strip():
                insert_name = "guide_oligo_pair"
            insert_sequence = insert_sequence or ""
        else:
            self._require_nonempty_string(insert_name, "insert_name")
            self._require_nonempty_string(insert_sequence, "insert_sequence")

        backbone_sequence = self._resolve_sequence_input(backbone_sequence, "backbone_sequence")
        if insert_sequence.strip():
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
            top_oligo_name=top_oligo_name,
            top_oligo_sequence=top_oligo_sequence,
            bottom_oligo_name=bottom_oligo_name,
            bottom_oligo_sequence=bottom_oligo_sequence,
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
            top_oligo_name=top_oligo_name,
            top_oligo_sequence=top_oligo_sequence,
            bottom_oligo_name=bottom_oligo_name,
            bottom_oligo_sequence=bottom_oligo_sequence,
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
            top_oligo_name=top_oligo_name,
            bottom_oligo_name=bottom_oligo_name,
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
            "gibsonassembly": "Gibson",
            "directsynthesis": "DirectSynthesis",
            "typeiis": "TypeIISOligoCloning",
            "typeiiscloning": "TypeIISOligoCloning",
            "typeiisoligocloning": "TypeIISOligoCloning",
            "restrictionligation": "RestrictionLigation",
            "restrictionfragmentligation": "RestrictionLigation",
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
        paper_info_json: Any = None,
    ) -> dict:
        """
        Resolve curated paper information from either:
        - a dict passed directly as paper_info_json,
        - a JSON string passed as paper_info_json,
        - a file in modules/crispr_tools/data/paper_info/ named <paper_id>.json,
        - or a v2 file named <paper_id>_v2.json.

        Supports both paper_important_info_v1 and paper_important_info_v2.
        This is intentionally done inside Create_Construction_File so Gemini does
        not need to pass nested fields like oligo_information between tool calls.
        """
        allowed_types = {"paper_important_info_v1", "paper_important_info_v2"}

        def _valid(data: Any) -> dict:
            if not isinstance(data, dict):
                return {}
            if data.get("resource_type") not in allowed_types:
                return {}
            data = dict(data)
            data.setdefault("paper_id", paper_id or data.get("paper_id", ""))
            data.setdefault("title", data.get("paper_id", ""))
            data.setdefault("source_pdf", "")
            data.setdefault("organism", "")
            data.setdefault("system", "")
            data.setdefault("targets", [])
            data.setdefault("vectors", [])
            data.setdefault("enzymes", [])
            data.setdefault("assembly_method", "")
            data.setdefault("delivery_method", [])
            data.setdefault("validation_methods", [])
            data.setdefault("key_constraints", [])
            data.setdefault("notes", [])
            data.setdefault("oligo_information", None)
            return data

        # 1. Direct dict input.
        direct = _valid(paper_info_json)
        if direct:
            return direct

        # 2. JSON string input.
        if isinstance(paper_info_json, str) and paper_info_json.strip():
            try:
                direct = _valid(json.loads(paper_info_json))
                if direct:
                    return direct
            except Exception:
                pass

        # 3. Local curated paper-info file lookup.
        normalized_id = self._normalize_paper_id(paper_id) if isinstance(paper_id, str) and paper_id.strip() else ""
        if not normalized_id:
            return {}

        paper_dir = self._paper_info_dir()
        candidate_paths = [
            paper_dir / f"{normalized_id}.json",
            paper_dir / f"{normalized_id}_v2.json",
        ]

        # If the user passed an id ending in _v2, also try the base name.
        if normalized_id.endswith("_v2"):
            base_id = normalized_id[:-3]
            candidate_paths.extend([
                paper_dir / f"{base_id}.json",
                paper_dir / f"{base_id}_v2.json",
            ])

        for path in candidate_paths:
            if not path.exists():
                continue
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except Exception:
                continue
            resolved = _valid(data)
            if resolved:
                return resolved

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
        top_oligo_name: str = "",
        top_oligo_sequence: str = "",
        bottom_oligo_name: str = "",
        bottom_oligo_sequence: str = "",
        enzyme: str = "",
    ) -> None:
        if assembly_strategy == "GoldenGate":
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

            if not enzyme.strip():
                missing_fields.append("enzyme")

            if missing_fields:
                raise ValueError(
                    f"Missing required fields for GoldenGate: {', '.join(missing_fields)}."
                )

        if assembly_strategy == "Gibson":
            missing_fields = []
            if not insert_forward_primer_name.strip():
                missing_fields.append("insert_forward_primer_name")
            if not insert_forward_primer_sequence.strip():
                missing_fields.append("insert_forward_primer_sequence")
            if not insert_reverse_primer_name.strip():
                missing_fields.append("insert_reverse_primer_name")
            if not insert_reverse_primer_sequence.strip():
                missing_fields.append("insert_reverse_primer_sequence")
            if missing_fields:
                raise ValueError(
                    f"Missing required fields for Gibson: {', '.join(missing_fields)}."
                )

        if assembly_strategy == "TypeIISOligoCloning":
            missing_fields = []
            if not top_oligo_name.strip():
                missing_fields.append("top_oligo_name")
            if not top_oligo_sequence.strip():
                missing_fields.append("top_oligo_sequence")
            if not bottom_oligo_name.strip():
                missing_fields.append("bottom_oligo_name")
            if not bottom_oligo_sequence.strip():
                missing_fields.append("bottom_oligo_sequence")
            if not enzyme.strip():
                missing_fields.append("enzyme")
            if missing_fields:
                raise ValueError(
                    "Missing required fields for TypeIISOligoCloning: "
                    + ", ".join(missing_fields)
                    + "."
                )

        if assembly_strategy == "RestrictionLigation":
            missing_fields = []
            if not insert_forward_primer_name.strip():
                missing_fields.append("insert_forward_primer_name")
            if not insert_forward_primer_sequence.strip():
                missing_fields.append("insert_forward_primer_sequence")
            if not insert_reverse_primer_name.strip():
                missing_fields.append("insert_reverse_primer_name")
            if not insert_reverse_primer_sequence.strip():
                missing_fields.append("insert_reverse_primer_sequence")
            if not enzyme.strip():
                missing_fields.append("enzyme")
            if missing_fields:
                raise ValueError(
                    "Missing required fields for RestrictionLigation: "
                    + ", ".join(missing_fields)
                    + "."
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
        top_oligo_name: str = "",
        top_oligo_sequence: str = "",
        bottom_oligo_name: str = "",
        bottom_oligo_sequence: str = "",
    ) -> list:
        parts = [
            {
                "part_type": "plasmid",
                "name": backbone_name,
                "sequence": self._normalize_sequence(backbone_sequence),
                "description": "Backbone plasmid",
            }
        ]

        if insert_sequence.strip():
            parts.append(
                {
                    "part_type": "dsdna",
                    "name": insert_name,
                    "sequence": self._normalize_sequence(insert_sequence),
                    "description": "Insert sequence",
                }
            )

        primer_entries = [
            ("oligo", insert_forward_primer_name, insert_forward_primer_sequence, "Insert forward primer"),
            ("oligo", insert_reverse_primer_name, insert_reverse_primer_sequence, "Insert reverse primer"),
            ("oligo", vector_forward_primer_name, vector_forward_primer_sequence, "Vector forward primer"),
            ("oligo", vector_reverse_primer_name, vector_reverse_primer_sequence, "Vector reverse primer"),
            ("oligo", top_oligo_name, top_oligo_sequence, "Top oligo for Type IIS oligo cloning"),
            ("oligo", bottom_oligo_name, bottom_oligo_sequence, "Bottom oligo for Type IIS oligo cloning"),
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
        top_oligo_name: str,
        bottom_oligo_name: str,
        enzyme: str,
        cell_strain: str,
        selection: str,
        temperature_c: int,
    ) -> list:
        operations = []

        if assembly_strategy in {"GoldenGate", "Gibson"}:
            insert_pcr_product = f"{insert_name}_pcr"
            vector_pcr_product = f"{backbone_name}_pcr"
            has_vector_primers = (
                vector_forward_primer_name.strip() and vector_reverse_primer_name.strip()
            )

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

            if has_vector_primers:
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

            assembly_inputs = (
                [vector_pcr_product, insert_pcr_product]
                if has_vector_primers
                else [backbone_name, insert_pcr_product]
            )

            if assembly_strategy == "GoldenGate":
                operations.append(
                    {
                        "step_number": len(operations) + 1,
                        "step_type": "GoldenGate",
                        "inputs": assembly_inputs,
                        "parameters": {"enzyme": enzyme},
                        "output": construct_name,
                    }
                )
            else:
                operations.append(
                    {
                        "step_number": len(operations) + 1,
                        "step_type": "Gibson",
                        "inputs": assembly_inputs,
                        "parameters": {"reagent": "GibsonMix", "overlap_bp": 20},
                        "output": construct_name,
                    }
                )

        elif assembly_strategy == "RestrictionLigation":
            insert_pcr_product = f"{insert_name}_pcr"

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
                    "step_type": "RestrictionLigation",
                    "inputs": [backbone_name, insert_pcr_product],
                    "parameters": {"enzyme": enzyme},
                    "output": construct_name,
                }
            )

        elif assembly_strategy == "TypeIISOligoCloning":
            operations.append(
                {
                    "step_number": 1,
                    "step_type": "TypeIISOligoCloning",
                    "inputs": [top_oligo_name, bottom_oligo_name, backbone_name],
                    "parameters": {
                        "enzyme": enzyme,
                        "overhangs": "vector-specific (e.g. CACC/AAAC)",
                    },
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

        elif step_type == "TypeIISOligoCloning":
            if len(inputs) != 3:
                raise ValueError(
                    f"TypeIISOligoCloning step {step_number} should have top oligo, bottom oligo, and backbone inputs."
                )
            if "enzyme" not in parameters or not str(parameters.get("enzyme", "")).strip():
                raise ValueError(f"TypeIISOligoCloning step {step_number} requires 'enzyme'.")

        elif step_type == "RestrictionLigation":
            if len(inputs) != 2:
                raise ValueError(f"RestrictionLigation step {step_number} should have vector and insert inputs.")
            if "enzyme" not in parameters or not str(parameters.get("enzyme", "")).strip():
                raise ValueError(f"RestrictionLigation step {step_number} requires 'enzyme'.")

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
            elif step_type == "TypeIISOligoCloning":
                lines.append(
                    f"{'TypeIIS':<14}"
                    f"{inputs[0] if len(inputs) > 0 else '':<14}"
                    f"{inputs[1] if len(inputs) > 1 else '':<14}"
                    f"{inputs[2] if len(inputs) > 2 else '':<18}"
                    f"{parameters.get('enzyme', ''):<10}"
                    f"{output}"
                )
            elif step_type == "RestrictionLigation":
                lines.append(
                    f"{'Ligate':<14}"
                    f"{inputs[0] if len(inputs) > 0 else '':<14}"
                    f"{inputs[1] if len(inputs) > 1 else '':<14}"
                    f"{parameters.get('enzyme', ''):<18}"
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

            # v2/protocol-style paper records include oligo_information. These are
            # reusable protocol/template papers rather than single construct papers,
            # so they need protocol shorthand instead of normal plasmid assembly shorthand.
            if resolved_info.get("oligo_information"):
                if title:
                    resolved_info["title"] = title
                if source_pdf:
                    resolved_info["source_pdf"] = source_pdf
                if organism:
                    resolved_info["organism"] = organism
                if system:
                    resolved_info["system"] = system
                if self._coerce_list(targets):
                    resolved_info["targets"] = self._coerce_list(targets)
                if self._coerce_list(vectors):
                    resolved_info["vectors"] = self._coerce_list(vectors)
                if self._coerce_list(enzymes):
                    resolved_info["enzymes"] = self._coerce_list(enzymes)
                if assembly_method:
                    resolved_info["assembly_method"] = assembly_method
                if self._coerce_list(delivery_method):
                    resolved_info["delivery_method"] = self._coerce_list(delivery_method)
                if self._coerce_list(validation_methods):
                    resolved_info["validation_methods"] = self._coerce_list(validation_methods)
                if self._coerce_list(key_constraints):
                    resolved_info["key_constraints"] = self._coerce_list(key_constraints)
                if self._coerce_list(paper_notes):
                    resolved_info["notes"] = self._coerce_list(paper_notes)
                return self._build_protocol_template_shorthand(resolved_info)

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

            # v2/protocol-style paper records include oligo_information. These are
            # reusable protocol/template papers rather than single construct papers,
            # so they need protocol shorthand instead of normal plasmid assembly shorthand.
            if resolved_info.get("oligo_information"):
                if title:
                    resolved_info["title"] = title
                if source_pdf:
                    resolved_info["source_pdf"] = source_pdf
                if organism:
                    resolved_info["organism"] = organism
                if system:
                    resolved_info["system"] = system
                if self._coerce_list(targets):
                    resolved_info["targets"] = self._coerce_list(targets)
                if self._coerce_list(vectors):
                    resolved_info["vectors"] = self._coerce_list(vectors)
                if self._coerce_list(enzymes):
                    resolved_info["enzymes"] = self._coerce_list(enzymes)
                if assembly_method:
                    resolved_info["assembly_method"] = assembly_method
                if self._coerce_list(delivery_method):
                    resolved_info["delivery_method"] = self._coerce_list(delivery_method)
                if self._coerce_list(validation_methods):
                    resolved_info["validation_methods"] = self._coerce_list(validation_methods)
                if self._coerce_list(key_constraints):
                    resolved_info["key_constraints"] = self._coerce_list(key_constraints)
                if self._coerce_list(paper_notes):
                    resolved_info["notes"] = self._coerce_list(paper_notes)
                return self._build_protocol_template_shorthand(resolved_info)

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

    def _build_protocol_template_shorthand(self, paper_info: dict) -> dict:
        """Build shorthand for protocol-style v2 paper records with oligo_information."""
        if not isinstance(paper_info, dict):
            raise ValueError("paper_info must be a dictionary for protocol shorthand generation.")

        paper_id = str(paper_info.get("paper_id", "paper")).strip() or "paper"
        title = str(paper_info.get("title", paper_id)).strip() or paper_id
        source_pdf = str(paper_info.get("source_pdf", "")).strip()
        organism = str(paper_info.get("organism", "organism")).strip() or "organism"
        system = str(paper_info.get("system", "system")).strip() or "system"

        vectors = self._coerce_list(paper_info.get("vectors", []))
        enzymes = self._coerce_list(paper_info.get("enzymes", []))
        delivery_methods = self._coerce_list(paper_info.get("delivery_method", []))
        validation_methods = self._coerce_list(paper_info.get("validation_methods", []))
        constraints = self._coerce_list(paper_info.get("key_constraints", []))
        notes = self._coerce_list(paper_info.get("notes", []))

        oligo_info = paper_info.get("oligo_information", {}) or {}
        if not isinstance(oligo_info, dict):
            oligo_info = {}
        explicit_sequences = oligo_info.get("explicit_sequences", []) or []
        design_rules = oligo_info.get("design_rules", []) or []
        interpretation_notes = self._coerce_list(oligo_info.get("interpretation_notes", []))

        declarations = []
        seen = set()

        def add_declaration(line: str) -> None:
            if line not in seen:
                declarations.append(line)
                seen.add(line)

        for vector in vectors:
            add_declaration(f"plasmid {self._sanitize_name(vector, 'vector')}")
        for enzyme in enzymes:
            add_declaration(f"enzyme {self._sanitize_name(enzyme, 'enzyme')}")
        add_declaration(f"cell {self._sanitize_name(organism, 'target_cells')}")

        steps = []
        for item in explicit_sequences:
            if not isinstance(item, dict):
                continue
            name = self._sanitize_name(str(item.get("name", "oligo_template")), "oligo_template")
            category = self._sanitize_name(str(item.get("category", "template")), "template")
            steps.append(f"oligo_template {name} {category}")

        enzyme_tokens = [self._sanitize_name(e, "enzyme") for e in enzymes]
        vector_tokens = [self._sanitize_name(v, "vector") for v in vectors]
        bbsI_present = any(e.lower() == "bbsi" for e in enzyme_tokens)
        if bbsI_present and vector_tokens:
            guide_vector = None
            for vector in vector_tokens:
                lower = vector.lower()
                if "px330" in lower or "sgrna" in lower or "grna" in lower:
                    guide_vector = vector
                    break
            guide_vector = guide_vector or vector_tokens[0]
            steps.append(f"digest {guide_vector} BbsI guide_vector_digest")
            steps.append("ligate guide_oligo_pair guide_vector_digest sgRNA_construct")

        explicit_text = json.dumps(explicit_sequences).lower()
        has_t7 = "t7" in explicit_text or any("t7" in e.lower() for e in enzymes)
        if has_t7:
            steps.append("ivt sgRNA_template T7 sgRNA")

        delivery = self._sanitize_name(delivery_methods[0], "delivery") if delivery_methods else "delivery"
        target_cell = self._sanitize_name(organism, "target_cells")
        if "microinjection" in delivery.lower() or "injection" in delivery.lower():
            steps.append(f"deliver Cas9 sgRNA {target_cell}_zygote {delivery} edited_zygote")
        elif "electroporation" in delivery.lower():
            steps.append(f"deliver Cas9 sgRNA {target_cell}_zygote electroporation edited_zygote")
        else:
            steps.append(f"deliver Cas9 sgRNA {target_cell} {delivery} edited_cells")

        validation_lower = " ".join(validation_methods).lower()
        if "pcr" in validation_lower:
            steps.append("pcr genotype_F genotype_R edited_sample genotype_amplicon")
        if "t7e1" in validation_lower or "survey" in validation_lower or "mismatch" in validation_lower:
            steps.append("screen genotype_amplicon mismatch_assay indel_screen")
        if "rflp" in validation_lower:
            steps.append("screen genotype_amplicon RFLP knockin_screen")
        if "sequence" in validation_lower or "sequencing" in validation_lower:
            steps.append("sequence genotype_amplicon Sanger confirmed_edit")

        design_rule_comments = []
        for rule in design_rules:
            if not isinstance(rule, dict):
                continue
            rule_name = rule.get("name", "design_rule")
            description = rule.get("description", "")
            if description:
                design_rule_comments.append(f"# design_rule {rule_name}: {description}")

        lines = [f"# Paper shorthand for: {title}", f"# paper_id: {paper_id}"]
        if source_pdf:
            lines.append(f"# source_pdf: {source_pdf}")
        lines.extend([f"# organism: {organism}", f"# system: {system}", "# protocol_style: true"])
        for constraint in constraints:
            lines.append(f"# constraint: {constraint}")
        for note in notes:
            lines.append(f"# note: {note}")
        for note in interpretation_notes:
            lines.append(f"# note: {note}")

        if design_rule_comments:
            lines.append("")
            lines.extend(design_rule_comments)
        if declarations:
            lines.append("")
            lines.extend(declarations)
        if steps:
            lines.append("")
            lines.extend(steps)
        shorthand_txt = "\n".join(lines)

        return {
            "mode": "paper_shorthand",
            "resource_type": "paper_shorthand_v1",
            "paper_id": paper_id,
            "title": title,
            "source_pdf": source_pdf,
            "organism": organism,
            "system": system,
            "protocol_style": True,
            "declarations": declarations,
            "steps": steps,
            "key_constraints": constraints,
            "notes": notes + interpretation_notes,
            "oligo_information": oligo_info,
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
                "Assembly strategy (GoldenGate, Gibson, DirectSynthesis, TypeIISOligoCloning): "
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
            top_oligo_name = ""
            top_oligo_sequence = ""
            bottom_oligo_name = ""
            bottom_oligo_sequence = ""
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

            if normalized_strategy in {"typeiis", "typeiiscloning", "typeiisoligocloning"}:
                print("\n--- Type IIS oligo cloning information required ---")
                top_oligo_name = prompt_required("Top oligo name: ")
                top_oligo_sequence = prompt_required("Top oligo sequence: ")
                bottom_oligo_name = prompt_required("Bottom oligo name: ")
                bottom_oligo_sequence = prompt_required("Bottom oligo sequence: ")
                enzyme = prompt_required("Type IIS enzyme (e.g. BbsI, BsaI): ")

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
                top_oligo_name=top_oligo_name,
                top_oligo_sequence=top_oligo_sequence,
                bottom_oligo_name=bottom_oligo_name,
                bottom_oligo_sequence=bottom_oligo_sequence,
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
