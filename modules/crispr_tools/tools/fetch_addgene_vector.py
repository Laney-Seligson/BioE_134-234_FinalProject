from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

_ADDGENE_API_BASE = "https://api.developers.addgene.org"


class FetchAddgeneVector:
    """
    Description:
        Fetches plasmid metadata from the Addgene Developers API by plasmid ID.
        Used as a fallback in the cloning oligo designer when a vector is not in
        the local preset registry. Requires an Addgene developer API token set
        via ADDGENE_API_KEY env var or passed directly.

    Input:
        addgene_id (int | str): Numeric Addgene plasmid ID (e.g. 67639 or
            "addgene:67639").
        api_key (str): Addgene developer token. Falls back to ADDGENE_API_KEY.
        include_sequences (bool): If True, fetches full sequences (requires
            catalog:retrieve-with-sequences scope). Default: False.

    Output:
        dict: status="ready" with keys: name, addgene_id, url, clone_method_raw,
              enzyme, promoter, backbone, bacterial_resistance, resistance_markers,
              growth_strain, description, article_doi; or status="needs_user_input"
              if the API key is missing.

    Tests:
        - Case:
            Input: addgene_id=67639, api_key="" (ADDGENE_API_KEY not set)
            Expected Output: status == "needs_user_input", "api_key" in missing_fields
            Description: Missing API key returns structured prompt.
        - Case:
            Input: addgene_id="not_a_number"
            Expected Exception: ValueError
            Description: Non-numeric ID raises ValueError.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        addgene_id: "int | str",
        api_key: str = "",
        include_sequences: bool = False,
    ) -> dict:
        raw = str(addgene_id).strip()
        if raw.lower().startswith("addgene:"):
            raw = raw[8:].strip()
        if not raw.isdigit():
            raise ValueError(
                f"addgene_id must be a numeric Addgene plasmid ID, got '{addgene_id}'."
            )
        plasmid_id = int(raw)

        key = api_key.strip() or os.environ.get("ADDGENE_API_KEY", "").strip()
        if not key:
            return {
                "status": "needs_user_input",
                "missing_fields": ["api_key"],
                "questions": [
                    "An Addgene developer API token is required. "
                    "Set the ADDGENE_API_KEY environment variable or pass api_key directly. "
                    "Register at https://www.addgene.org/developer to get a token."
                ],
            }

        endpoint = (
            f"/catalog/plasmid-with-sequences/{plasmid_id}/"
            if include_sequences
            else f"/catalog/plasmid/{plasmid_id}/"
        )
        req = urllib.request.Request(
            _ADDGENE_API_BASE + endpoint,
            headers={"Authorization": f"Token {key}"},
        )
        try:
            with urllib.request.urlopen(req, timeout=10) as resp:
                data = json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise ValueError(f"Addgene plasmid #{plasmid_id} not found.")
            if e.code == 401:
                raise ValueError(
                    "Addgene API authentication failed — check your API key."
                )
            raise ValueError(f"Addgene API error {e.code}: {e.reason}")
        except urllib.error.URLError as e:
            raise ValueError(
                f"Network error fetching Addgene #{plasmid_id}: {e.reason}"
            )

        return self._parse(data)

    def _parse(self, data: dict) -> dict:
        inserts = data.get("inserts") or []
        ins_cloning = (inserts[0].get("cloning") or {}) if inserts else {}
        cloning_info = data.get("cloning") or {}
        article = data.get("article") or {}

        enzyme = (
            ins_cloning.get("cloning_site_5")
            or ins_cloning.get("cloning_site_3")
            or ""
        )

        return {
            "status": "ready",
            "addgene_id": data.get("id"),
            "name": data.get("name", ""),
            "url": data.get("url", f"https://www.addgene.org/{data.get('id', '')}/"),
            "description": data.get("description", ""),
            "clone_method_raw": ins_cloning.get("clone_method", ""),
            "enzyme": enzyme,
            "promoter": cloning_info.get("promoter", ""),
            "backbone": cloning_info.get("backbone", ""),
            "vector_types": cloning_info.get("vector_types", []),
            "bacterial_resistance": data.get("bacterial_resistance", ""),
            "resistance_markers": data.get("resistance_markers", []),
            "growth_strain": data.get("growth_strain", ""),
            "growth_temp": data.get("growth_temp"),
            "growth_notes": data.get("growth_notes", ""),
            "article_doi": article.get("doi", ""),
            "article_pmid": article.get("pubmed_id"),
            **({"sequences": data["sequences"]} if "sequences" in data else {}),
        }


_instance = FetchAddgeneVector()
_instance.initiate()
fetch_addgene_vector = _instance.run
