import json
from pathlib import Path
from typing import Any


class GetPaperInfo:
    """
    Description:
        Load a curated paper information JSON file from the module data directory.
        This tool is intended for non-sequence literature metadata such as paper-
        derived workflow facts that are later consumed by other tools like
        create_construction_file in paper_shorthand mode.

    Input:
        paper_id (str): Stable paper identifier, such as
            'miao_2013_targeted_mutagenesis_rice'.

    Output:
        dict: Parsed paper information resource.
    """

    def initiate(self) -> None:
        # Assumes this file will live at modules/crispr_tools/tools/get_paper_info.py
        # and curated JSON files will live at modules/crispr_tools/data/paper_info/*.json
        self.paper_dir = (
            Path(__file__).resolve().parents[1] / "data" / "paper_info"
        )

    def run(self, paper_id: str) -> dict:
        self._require_nonempty_string(paper_id, "paper_id")
        normalized_id = self._normalize_paper_id(paper_id)
        file_path = self.paper_dir / f"{normalized_id}.json"

        if not file_path.exists():
            available = self._available_paper_ids()
            raise ValueError(
                f"Paper info file not found for paper_id='{normalized_id}'. "
                f"Expected file: {file_path}. Available paper_ids: {available}"
            )

        try:
            with open(file_path, "r", encoding="utf-8") as f:
                data = json.load(f)
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in {file_path}: {e}") from e

        if not isinstance(data, dict):
            raise ValueError("Paper info JSON must decode to an object/dictionary.")

        if data.get("resource_type") != "paper_important_info_v1":
            raise ValueError(
                "paper info JSON must include resource_type='paper_important_info_v1'."
            )

        if data.get("paper_id") != normalized_id:
            raise ValueError(
                f"paper_id mismatch: file requested '{normalized_id}' but JSON contains "
                f"'{data.get('paper_id', '')}'."
            )

        return data

    def _require_nonempty_string(self, value: str, field_name: str) -> None:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{field_name} must be a non-empty string.")

    def _normalize_paper_id(self, paper_id: str) -> str:
        value = paper_id.strip()
        prefix = "resource://paper_info/"
        if value.startswith(prefix):
            value = value[len(prefix):]
        if value.endswith(".json"):
            value = value[:-5]
        return value

    def _available_paper_ids(self) -> list[str]:
        if not self.paper_dir.exists():
            return []
        return sorted(p.stem for p in self.paper_dir.glob("*.json"))


_instance = GetPaperInfo()
_instance.initiate()
get_paper_info = _instance.run
