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
   - target_type="genomic_locus"   → chromosomal coordinates via esearch/esummary/efetch
   - target_type="cds"             → annotated CDS feature via elink/efetch (GenBank)

Citations:

National Center for Biotechnology Information. (n.d.). Entrez programming utilities help (E-utilities). U.S. National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK25501/

National Center for Biotechnology Information. (n.d.). Gene database help. U.S. National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK3841/

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
from Bio import Entrez, SeqIO

from modules.seq_basics._plumbing.resolve import resolve_to_seq, list_resources
from modules.crispr_tools.tools._utils import normalize_organism, VALID_DNA


_ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class FetchTargetSequence:
    """
    Description:
        Resolves a gene name, local resource key, or raw DNA string to a
        clean uppercase DNA sequence. For gene-name queries, supports two
        target types:
          - "genomic_locus" (default): full gene region from chromosomal
            coordinates, including introns and flanking sequence.
          - "cds": annotated coding sequence only, extracted from the top
            RefSeq mRNA record. Falls back to the full mRNA record if no
            CDS feature is annotated.

    Input:
        query (str): Raw DNA sequence, local resource name, or gene name.
        organism (str): Organism for NCBI lookup. Default: "Escherichia coli".
        target_type (str): "genomic_locus" or "cds". Default: "genomic_locus".

    Output:
        dict: Keys: sequence (str), source (str), resource (str),
              organism (str), target_type (str), length (int), note (str).
              CDS results also include product (str).

    Tests:
        - Case:
            Input: query="ATGCATGCATGC"
            Expected Output: source == "raw_input", sequence == "ATGCATGCATGC"
            Description: Raw DNA is returned as-is.
        - Case:
            Input: query="rpsL", organism="Escherichia coli", target_type="genomic_locus"
            Expected Output: source == "ncbi", target_type == "genomic_locus", length > 0
            Description: Gene name triggers NCBI genomic locus fetch.
        - Case:
            Input: query="lacZ", organism="Escherichia coli", target_type="cds"
            Expected Output: source == "ncbi", target_type == "cds", length > 0
            Description: Gene name triggers NCBI CDS fetch.
        - Case:
            Input: query=""
            Expected Exception: ValueError
            Description: Empty query raises ValueError.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        query: str,
        organism: str = "Escherichia coli",
        target_type: str = "genomic_locus",
    ) -> dict:
        query = query.strip()
        if not query:
            raise ValueError("query must not be empty.")

        # ── 1. Raw DNA sequence ───────────────────────────────────────────────
        cleaned = re.sub(r"\s", "", query).upper()
        if set(cleaned) <= VALID_DNA and len(cleaned) >= 10:
            return {
                "sequence": cleaned,
                "source": "raw_input",
                "resource": "raw_input",
                "organism": "",
                "target_type": target_type,
                "length": len(cleaned),
                "note": "Input was a raw DNA sequence — used as-is.",
            }

        # ── 2. Local resource (MCP registry) ─────────────────────────────────
        resources = list_resources()
        resource_key = query
        if resource_key not in resources:
            for k in resources:
                if Path(k).stem.lower() == query.lower() or k.lower() == query.lower():
                    resource_key = k
                    break

        if resource_key in resources:
            seq = resolve_to_seq(resource_key)
            note = f"Resolved from local bundled resource '{resource_key}'."
            default_organism = "Escherichia coli"
            if organism and organism.strip().lower() != default_organism.lower():
                note += (
                    f" WARNING: organism '{organism}' was specified but ignored — "
                    f"'{resource_key}' is a local resource, not a gene lookup. "
                    "To fetch a gene from this organism, use the gene name instead."
                )
            return {
                "sequence": seq,
                "source": "local_resource",
                "resource": resource_key,
                "organism": "",
                "target_type": target_type,
                "length": len(seq),
                "note": note,
            }

        # ── 3. NCBI Entrez fetch ──────────────────────────────────────────────
        organism = normalize_organism(organism)
        _configure_biopython_entrez()

        if target_type == "cds":
            seq, gene_id, accession, product = _fetch_cds_from_ncbi(query, organism)
            return {
                "sequence": seq,
                "source": "ncbi",
                "resource": query,
                "organism": organism,
                "target_type": "cds",
                "ncbi_gene_id": gene_id,
                "ncbi_accession": accession,
                "product": product,
                "length": len(seq),
                "note": (
                    f"Fetched CDS from NCBI Entrez: gene '{query}' in '{organism}'. "
                    f"Accession: {accession}. "
                    "Verify the record matches your intended target before use."
                ),
            }
        else:
            seq, gene_id, accession = _fetch_locus_from_ncbi(query, organism)
            return {
                "sequence": seq,
                "source": "ncbi",
                "resource": query,
                "organism": organism,
                "target_type": "genomic_locus",
                "ncbi_gene_id": gene_id,
                "ncbi_accession": accession,
                "length": len(seq),
                "note": (
                    f"Fetched genomic locus from NCBI Entrez: gene '{query}' in '{organism}'. "
                    f"Nucleotide record: {accession}. "
                    "Verify the record matches your intended target before use."
                ),
            }


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _configure_biopython_entrez() -> None:
    Entrez.email = os.environ.get("NCBI_EMAIL", "user@example.com")
    api_key = os.environ.get("NCBI_API_KEY")
    if api_key:
        Entrez.api_key = api_key


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


def _esearch_gene_id(gene_name: str, organism: str) -> str:
    """Return the top NCBI Gene ID for gene_name in organism."""
    search_strategies = [
        f"{gene_name}[Gene] AND {organism}[Organism]",
        f"{gene_name}[Gene/Protein Name] AND {organism}[Organism]",
        f"{gene_name}[All Fields] AND {organism}[Organism]",
    ]
    for term in search_strategies:
        resp = _ncbi_get(
            f"{_ENTREZ_BASE}/esearch.fcgi",
            {"db": "gene", "term": term, "retmode": "json"},
        )
        ids = resp.json()["esearchresult"]["idlist"]
        if ids:
            return ids[0]
    raise ValueError(
        f"Gene '{gene_name}' not found in NCBI for organism '{organism}'. "
        "Check the spelling or try a more specific organism name."
    )


# ---------------------------------------------------------------------------
# Genomic locus path (chromosomal coordinates → FASTA slice)
# ---------------------------------------------------------------------------

def _fetch_locus_from_ncbi(gene_name: str, organism: str) -> tuple[str, str, str]:
    """Return (sequence, gene_id, accession) for the full genomic gene region."""
    gene_id = _esearch_gene_id(gene_name, organism)

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

    loc       = genomic_info[0]
    accession = loc["chraccver"]
    chr_start = int(loc["chrstart"])
    chr_stop  = int(loc["chrstop"])

    if chr_start > chr_stop:
        seq_start, seq_stop, strand = chr_stop + 1, chr_start + 1, 2
    else:
        seq_start, seq_stop, strand = chr_start + 1, chr_stop + 1, 1

    fasta_resp = _ncbi_get(
        f"{_ENTREZ_BASE}/efetch.fcgi",
        {
            "db":        "nuccore",
            "id":        accession,
            "rettype":   "fasta",
            "retmode":   "text",
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

    lines   = fasta_text.splitlines()
    raw_seq = "".join(lines[1:])
    seq     = re.sub(r"[^ATGCNatgcn]", "", raw_seq).upper()

    if not seq:
        raise ValueError(
            f"Fetched FASTA for {accession}:{seq_start}-{seq_stop} contained no sequence."
        )

    return seq, gene_id, accession


# ---------------------------------------------------------------------------
# CDS path (mRNA elink → GenBank CDS feature extraction)
# ---------------------------------------------------------------------------

def _fetch_cds_from_ncbi(gene_name: str, organism: str) -> tuple[str, str, str, str]:
    """Return (cds_sequence, gene_id, accession, product) for the annotated CDS."""
    gene_id = _esearch_gene_id(gene_name, organism)

    # elink: Gene → RefSeq mRNA accessions
    handle = Entrez.elink(
        dbfrom="gene", db="nucleotide", id=gene_id,
        linkname="gene_nuccore_refseqrna",
    )
    link_result = Entrez.read(handle)
    handle.close()

    accession_ids: list[str] = []
    for link_set in link_result:
        for link in link_set.get("LinkSetDb", []):
            for item in link.get("Link", []):
                accession_ids.append(item["Id"])

    if not accession_ids:
        # fallback: nucleotide search for mRNA records
        handle = Entrez.esearch(
            db="nucleotide",
            term=f"{gene_name}[Gene Name] AND {organism}[Organism] AND mRNA[Filter]",
            retmax=3,
        )
        nt_result = Entrez.read(handle)
        handle.close()
        accession_ids = nt_result.get("IdList", [])

    if not accession_ids:
        raise ValueError(
            f"No mRNA records found for gene '{gene_name}' in '{organism}'. "
            "The gene may lack a RefSeq mRNA record."
        )

    handle = Entrez.efetch(
        db="nucleotide", id=accession_ids[0], rettype="gb", retmode="text"
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()
    accession = record.id

    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = str(feature.extract(record.seq)).upper()
            product = feature.qualifiers.get("product", ["unknown"])[0]
            return cds_seq, gene_id, accession, product

    # no CDS annotation — return the full mRNA record sequence
    full_seq = str(record.seq).upper()
    return full_seq, gene_id, accession, "unknown"


# ---------------------------------------------------------------------------
# Module-level singleton
# ---------------------------------------------------------------------------

_instance = FetchTargetSequence()
_instance.initiate()
fetch_target_sequence = _instance.run
