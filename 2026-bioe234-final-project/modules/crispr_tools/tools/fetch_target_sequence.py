"""
Fetch Target Sequence
=====================
Resolves a gene name, local reference key, or raw DNA sequence into a clean
DNA string that can be passed directly to design_cas9_grna, design_cas12a_crrna,
or crispr_cas_selector.

Resolution order
----------------
1. Raw DNA string (only A/T/G/C)   → returned as-is after validation
2. Known local resource name        → resolved from bundled file via MCP registry
3. Gene name + organism             → fetched from NCBI Entrez (live API call)

The NCBI path uses three Entrez calls:
  esearch  → gene ID from gene name + organism
  elink    → linked nucleotide record ID
  efetch   → FASTA for that nucleotide record

No API key is required for low-volume use, but NCBI rate-limits unauthenticated
requests to 3/second. Set NCBI_API_KEY in your .env for higher limits.
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path
from typing import Optional

import requests

from modules.seq_basics._plumbing.resolve import resolve_to_seq, list_resources


_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_VALID_DNA   = set("ATGC")


class FetchTargetSequence:

    def initiate(self) -> None:
        pass

    def run(
        self,
        query: str,
        organism: str = "Escherichia coli",
    ) -> dict:
        """
        Resolve a query to a clean DNA sequence.

        Parameters
        ----------
        query : str
            One of:
            - A raw DNA string (A/T/G/C only) — returned immediately.
            - A local resource name (e.g. "ecoli_rpsl", "pBR322") — resolved
              from the bundled MCP resource registry.
            - A gene name (e.g. "rpsL", "lacZ", "recA") — fetched from NCBI.
        organism : str
            Organism name used for NCBI gene search. Ignored when query is a
            raw sequence or local resource name. Default: "Escherichia coli".

        Returns
        -------
        dict with keys:
            sequence  (str)  — clean uppercase DNA string
            source    (str)  — "raw_input" | "local_resource" | "ncbi"
            resource  (str)  — resource name or gene name used
            organism  (str)  — organism used (NCBI path only)
            length    (int)  — sequence length in bp
            note      (str)  — human-readable provenance note
        """
        query = query.strip()
        if not query:
            raise ValueError("query must not be empty.")

        # ── 1. Raw DNA sequence ───────────────────────────────────────────────
        cleaned = re.sub(r"\s", "", query).upper()
        if set(cleaned) <= _VALID_DNA and len(cleaned) >= 10:
            return {
                "sequence": cleaned,
                "source": "raw_input",
                "resource": "raw_input",
                "organism": "",
                "length": len(cleaned),
                "note": "Input was a raw DNA sequence — used as-is.",
            }

        # ── 2. Local resource (MCP registry) ─────────────────────────────────
        resources = list_resources()
        # Accept bare filename stem too (e.g. "gene.fna" → "gene")
        resource_key = query
        if resource_key not in resources:
            # try matching by stem
            for k in resources:
                if Path(k).stem.lower() == query.lower() or k.lower() == query.lower():
                    resource_key = k
                    break

        if resource_key in resources:
            seq = resolve_to_seq(resource_key)
            return {
                "sequence": seq,
                "source": "local_resource",
                "resource": resource_key,
                "organism": "",
                "length": len(seq),
                "note": f"Resolved from local bundled resource '{resource_key}'.",
            }

        # ── 3. NCBI Entrez fetch ──────────────────────────────────────────────
        seq, nuc_id, accession = _fetch_from_ncbi(query, organism)
        return {
            "sequence": seq,
            "source": "ncbi",
            "resource": query,
            "organism": organism,
            "ncbi_gene_id": nuc_id,
            "ncbi_accession": accession,
            "length": len(seq),
            "note": (
                f"Fetched from NCBI Entrez: gene '{query}' in '{organism}'. "
                f"Nucleotide record: {accession}. "
                "Verify the record matches your intended target before use."
            ),
        }


# ---------------------------------------------------------------------------
# NCBI helpers
# ---------------------------------------------------------------------------

def _ncbi_get(url: str, params: dict) -> dict:
    """GET an Entrez endpoint with rate-limit back-off."""
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        params["api_key"] = api_key

    for attempt in range(3):
        resp = requests.get(url, params=params, timeout=15)
        if resp.status_code == 429:
            time.sleep(2 ** attempt)
            continue
        resp.raise_for_status()
        return resp
    raise RuntimeError("NCBI API returned 429 Too Many Requests after retries.")


def _fetch_from_ncbi(gene_name: str, organism: str) -> tuple[str, str, str]:
    """
    Return (sequence, gene_id, accession) for gene_name in organism.

    Uses three Entrez calls:
      esearch   → gene ID from gene name + organism
      esummary  → chromosomal coordinates (accession, start, stop, strand)
      efetch    → FASTA for just that gene region

    Fetching only the gene region (via seq_start/seq_stop) avoids pulling the
    entire chromosome, which could be millions of bp.
    """
    # Step 1: gene ID
    search_resp = _ncbi_get(
        f"{_ENTREZ_BASE}/esearch.fcgi",
        {"db": "gene", "term": f"{gene_name}[Gene] AND {organism}[Organism]",
         "retmode": "json"},
    )
    gene_ids = search_resp.json()["esearchresult"]["idlist"]
    if not gene_ids:
        raise ValueError(
            f"Gene '{gene_name}' not found in NCBI for organism '{organism}'. "
            "Check the spelling or try a more specific organism name."
        )
    gene_id = gene_ids[0]

    # Step 2: chromosomal coordinates from gene summary
    summary_resp = _ncbi_get(
        f"{_ENTREZ_BASE}/esummary.fcgi",
        {"db": "gene", "id": gene_id, "retmode": "json"},
    )
    doc = summary_resp.json()["result"][gene_id]
    genomic_info = doc.get("genomicinfo", [])
    if not genomic_info:
        raise ValueError(
            f"No genomic coordinates found for gene '{gene_name}' (gene ID {gene_id}). "
            "The gene record exists in NCBI but has no mapped chromosomal location."
        )

    loc        = genomic_info[0]
    accession  = loc["chraccver"]           # e.g. "NC_000913.3"
    chr_start  = int(loc["chrstart"])       # 0-based; may be > chrstop on minus strand
    chr_stop   = int(loc["chrstop"])

    # NCBI genomicinfo uses chrstart > chrstop to indicate minus strand
    if chr_start > chr_stop:
        seq_start, seq_stop, strand = chr_stop + 1, chr_start + 1, 2
    else:
        seq_start, seq_stop, strand = chr_start + 1, chr_stop + 1, 1

    # Step 3: fetch just the gene region as FASTA
    fasta_resp = _ncbi_get(
        f"{_ENTREZ_BASE}/efetch.fcgi",
        {
            "db":       "nuccore",
            "id":       accession,
            "rettype":  "fasta",
            "retmode":  "text",
            "seq_start": seq_start,
            "seq_stop":  seq_stop,
            "strand":    strand,
        },
    )
    fasta_text = fasta_resp.text.strip()
    if not fasta_text.startswith(">"):
        raise ValueError(
            f"NCBI returned unexpected content for {accession}:{seq_start}-{seq_stop}."
        )

    lines  = fasta_text.splitlines()
    raw_seq = "".join(lines[1:])
    seq     = re.sub(r"[^ATGCNatgcn]", "", raw_seq).upper()

    if not seq:
        raise ValueError(
            f"Fetched FASTA for {accession}:{seq_start}-{seq_stop} contained no sequence."
        )

    return seq, gene_id, accession


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

_instance = FetchTargetSequence()
_instance.initiate()
fetch_target_sequence = _instance.run
