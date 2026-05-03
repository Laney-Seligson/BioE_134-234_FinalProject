"""
semantic_wrapper.py

Semantic wrapper for biology concept discovery using:
- OLS4 for GO / ontology term search

Pipeline:
Natural language query -> parsed query -> OLS4 GO / ontology terms

This tool is intentionally separated from gene retrieval.
Gene retrieval should be handled by the annotation tool:
go_term_gene_lookup
"""

from __future__ import annotations

import re
import time
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Optional

import requests


OLS4_SEARCH_URL = "https://www.ebi.ac.uk/ols4/api/search"


@dataclass
class ParsedQuery:
    raw_query: str
    organism: Optional[str]
    keywords: List[str]
    ontology_terms: List[str]


@dataclass
class GOTerm:
    go_id: str
    label: str
    definition: Optional[str] = None


@dataclass
class Result:
    parsed_query: ParsedQuery
    go_terms: List[GOTerm]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "parsed_query": asdict(self.parsed_query),
            "go_terms": [asdict(x) for x in self.go_terms],
        }


def get_json(
    url: str,
    params: Optional[Dict[str, Any]] = None,
    timeout: tuple[int, int] = (10, 45),
    max_retries: int = 3,
) -> Dict[str, Any]:
    last_error = None

    for attempt in range(max_retries):
        try:
            response = requests.get(
                url,
                params=params,
                timeout=timeout,
                headers={"User-Agent": "semantic-wrapper/1.0"},
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as exc:
            last_error = exc
            print(
                f"[semantic_wrapper] attempt {attempt + 1}/{max_retries} "
                f"failed for {url}: {exc}"
            )
            time.sleep(1.5 * (attempt + 1))

    print(f"[semantic_wrapper] FINAL FAIL for {url}: {last_error}")
    return {}


ORGANISMS = {
    "yeast": "Saccharomyces cerevisiae",
    "saccharomyces cerevisiae": "Saccharomyces cerevisiae",
    "human": "Homo sapiens",
    "homo sapiens": "Homo sapiens",
    "mouse": "Mus musculus",
    "mus musculus": "Mus musculus",
    "ecoli": "Escherichia coli",
    "e coli": "Escherichia coli",
    "e. coli": "Escherichia coli",
    "escherichia coli": "Escherichia coli",
}


def parse_query(q: str) -> ParsedQuery:
    q_lower = q.lower()

    organism = None
    for key, value in ORGANISMS.items():
        if key in q_lower:
            organism = value
            break

    tokens = re.findall(r"[a-zA-Z]+", q_lower)

    stop = {
        "find",
        "genes",
        "gene",
        "related",
        "to",
        "in",
        "of",
        "and",
        "for",
        "associated",
        "involved",
        "with",
        "show",
        "me",
        "the",
        "return",
        "most",
        "relevant",
        "matches",
        "biological",
        "functions",
    }

    keywords = [t for t in tokens if t not in stop]

    ontology_terms: List[str] = []

    if "oxidative" in q_lower and "stress" in q_lower:
        ontology_terms.append("response to oxidative stress")
    elif "dna" in q_lower and "repair" in q_lower:
        ontology_terms.append("DNA repair")
    elif "immune" in q_lower and "response" in q_lower:
        ontology_terms.append("immune response")
    elif "cell" in q_lower and "cycle" in q_lower:
        ontology_terms.append("cell cycle")
    elif keywords:
        ontology_terms.append(" ".join(keywords[:4]))

    return ParsedQuery(
        raw_query=q,
        organism=organism,
        keywords=keywords,
        ontology_terms=ontology_terms,
    )


def search_go(term: str, rows: int = 5) -> List[GOTerm]:
    data = get_json(
        OLS4_SEARCH_URL,
        params={
            "q": term,
            "ontology": "go",
            "rows": rows,
        },
    )

    docs = data.get("response", {}).get("docs", [])
    results: List[GOTerm] = []
    seen = set()

    for doc in docs:
        go_id = doc.get("obo_id")
        label = doc.get("label")
        description = doc.get("description")

        if isinstance(description, list):
            definition = description[0] if description else None
        else:
            definition = description

        if go_id and label and go_id not in seen:
            seen.add(go_id)
            results.append(
                GOTerm(
                    go_id=go_id,
                    label=label,
                    definition=definition,
                )
            )

    return results


class SemanticGeneWrapper:
    def run(self, query: str) -> Result:
        parsed = parse_query(query)

        go_terms: List[GOTerm] = []
        for term in parsed.ontology_terms:
            go_terms.extend(search_go(term))
            time.sleep(0.2)

        return Result(
            parsed_query=parsed,
            go_terms=go_terms,
        )