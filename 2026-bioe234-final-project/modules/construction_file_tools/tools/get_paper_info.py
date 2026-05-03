
import json
from pathlib import Path
from typing import Any


class GetPaperInfo:
    """
    Description:
        Load a curated paper information JSON record from the construction_file_tools data directory.
        Supports both paper_important_info_v1 and paper_important_info_v2 records and
        normalizes the result so downstream tools can use a consistent structure.

    Input:
        paper_id (str): Identifier of the paper JSON file to load, without the .json suffix.

    Output:
        dict: Normalized paper information record.

    Tests:
        - Case:
            Input: paper_id="miao_2013_targeted_mutagenesis_rice"
            Expected Output: A dictionary with resource_type, paper_id, title, organism, system.
            Description: Loads a v1 paper record.
        - Case:
            Input: paper_id="hall_2018_genome_editing_mice_crispr_cas9"
            Expected Output: A dictionary with resource_type, paper_id, title, organism, system, oligo_information.
            Description: Loads a v2 paper record.
        - Case:
            Input: paper_id="missing_paper"
            Expected Output: ValueError
            Description: Missing paper file should raise a clear error.
    """

    ALLOWED_RESOURCE_TYPES = {
        "paper_important_info_v1",
        "paper_important_info_v2",
    }

    def initiate(self) -> None:
        pass

    def _candidate_paths(self, paper_id: str) -> list[Path]:
        base_dirs = [
            Path("modules/construction_file_tools/data/paper_info"),
            Path("modules/construction_file_tools/data"),
        ]
        candidates: list[Path] = []
        for base in base_dirs:
            candidates.append(base / f"{paper_id}.json")
            # helpful fallback for files stored with an explicit v2 suffix
            candidates.append(base / f"{paper_id}_v2.json")
        return candidates

    def _load_json(self, paper_id: str) -> dict[str, Any]:
        for path in self._candidate_paths(paper_id):
            if path.exists() and path.is_file():
                text = path.read_text(encoding="utf-8")
                data = json.loads(text)
                if not isinstance(data, dict):
                    raise ValueError("Paper info JSON must decode to an object/dictionary.")
                return data
        searched = ", ".join(str(p) for p in self._candidate_paths(paper_id))
        raise ValueError(
            f"Could not find paper info JSON for paper_id='{paper_id}'. "
            f"Searched: {searched}"
        )

    def _normalize_v2(self, data: dict[str, Any]) -> dict[str, Any]:
        # Preserve v2-only content while ensuring v1-style fields always exist.
        normalized = dict(data)
        normalized.setdefault("targets", [])
        normalized.setdefault("vectors", [])
        normalized.setdefault("enzymes", [])
        normalized.setdefault("assembly_method", "Unknown")
        normalized.setdefault("delivery_method", [])
        normalized.setdefault("validation_methods", [])
        normalized.setdefault("key_constraints", [])
        normalized.setdefault("notes", [])
        normalized.setdefault("oligo_information", {
            "explicit_sequences": [],
            "design_rules": [],
            "interpretation_notes": [],
        })
        return normalized

    def _normalize_v1(self, data: dict[str, Any]) -> dict[str, Any]:
        normalized = dict(data)
        normalized.setdefault("targets", [])
        normalized.setdefault("vectors", [])
        normalized.setdefault("enzymes", [])
        normalized.setdefault("assembly_method", "Unknown")
        normalized.setdefault("delivery_method", [])
        normalized.setdefault("validation_methods", [])
        normalized.setdefault("key_constraints", [])
        normalized.setdefault("notes", [])
        return normalized

    def run(self, paper_id: str) -> dict[str, Any]:
        if not isinstance(paper_id, str) or not paper_id.strip():
            raise ValueError("paper_id must be a non-empty string.")

        paper_id = paper_id.strip()
        data = self._load_json(paper_id)

        resource_type = data.get("resource_type")
        if resource_type not in self.ALLOWED_RESOURCE_TYPES:
            raise ValueError(
                "paper info JSON must include "
                "resource_type='paper_important_info_v1' or "
                "resource_type='paper_important_info_v2'."
            )

        if resource_type == "paper_important_info_v2":
            normalized = self._normalize_v2(data)
        else:
            normalized = self._normalize_v1(data)

        # Ensure paper_id is always present and matches the requested ID if omitted.
        normalized.setdefault("paper_id", paper_id)

        # Useful minimal required fields for downstream tools.
        title = normalized.get("title", "")
        organism = normalized.get("organism", "")
        system = normalized.get("system", "")

        if not isinstance(title, str):
            raise ValueError("paper info JSON field 'title' must be a string.")
        if organism and not isinstance(organism, str):
            raise ValueError("paper info JSON field 'organism' must be a string.")
        if system and not isinstance(system, str):
            raise ValueError("paper info JSON field 'system' must be a string.")

        return normalized


_instance = GetPaperInfo()
_instance.initiate()
get_paper_info = _instance.run
