from __future__ import annotations

from typing import Optional
from urllib.error import HTTPError, URLError

from Bio import Entrez, SeqIO


_ORGANISM_ALIASES = {
    "ecoli": "Escherichia coli",
    "e.coli": "Escherichia coli",
    "e. coli": "Escherichia coli",
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "yeast": "Saccharomyces cerevisiae",
    "arabidopsis": "Arabidopsis thaliana",
    "zebrafish": "Danio rerio",
    "fly": "Drosophila melanogaster",
    "worm": "Caenorhabditis elegans",
    "rat": "Rattus norvegicus",
    "tobacco": "Nicotiana tabacum",
    "rice": "Oryza sativa",
    "corn": "Zea mays",
    "maize": "Zea mays",
}


class LookupGeneSequence:
    """
    Description:
        Fetches the coding DNA sequence (CDS) for a gene from NCBI Entrez
        given a gene name and organism. Queries the NCBI Gene database for
        the gene ID, then retrieves the mRNA record and extracts the CDS.

        Common organism names are accepted (e.g. "E. coli", "mouse",
        "Arabidopsis") in addition to full scientific names. Returns the
        CDS nucleotide sequence plus key metadata (gene ID, accession,
        protein product name).

        Uses the nucleotide database fallback if no annotated CDS is found
        in the top mRNA hit, returning the full record sequence instead.

        An NCBI email is required (passed once at initiate time or as a
        parameter) so NCBI can contact you if your queries cause issues.
        A free NCBI API key increases the rate limit from 3 to 10
        requests/second, but is not required.

    Input:
        gene_name (str): gene symbol or name, e.g. "cscB", "cas9", "lacZ"
        organism (str): organism name or common alias, e.g. "E. coli",
                        "mouse", "Arabidopsis thaliana"
        email (str): email address for NCBI Entrez (required by NCBI policy).
                     defaults to the email set during initiate().
        api_key (str): optional NCBI API key to increase rate limit.
        max_results (int): how many Gene DB hits to try before giving up.
                           default 5.

    Output:
        dict with keys:
            - gene_name: the queried gene name (echoed)
            - organism: the resolved scientific name
            - gene_id: NCBI Gene ID (string)
            - accession: the mRNA/nucleotide accession used
            - sequence: the CDS nucleotide sequence (5' to 3')
            - sequence_length: length in bp
            - product: protein product name from the CDS annotation
            - source: "CDS annotation" or "full record"
            - summary: one-sentence description

    Tests:
        - Case:
            Input: gene_name="lacZ", organism="E. coli", email="test@example.com"
            Expected Output: "sequence" is a non-empty string of ATGC; "gene_id" is non-empty
            Description: classic E. coli gene should be retrievable from NCBI.
        - Case:
            Input: gene_name="GAPDH", organism="human", email="test@example.com"
            Expected Output: "sequence_length" > 100; "organism" == "Homo sapiens"
            Description: common alias "human" resolves to Homo sapiens.
        - Case:
            Input: gene_name="", organism="E. coli", email="test@example.com"
            Expected Exception: ValueError
            Description: empty gene name raises ValueError.
        - Case:
            Input: gene_name="lacZ", organism="", email="test@example.com"
            Expected Exception: ValueError
            Description: empty organism raises ValueError.
        - Case:
            Input: gene_name="xyzzy_fake_gene_9999", organism="E. coli", email="test@example.com"
            Expected Exception: ValueError
            Description: gene not found in NCBI raises ValueError with a clear message.
    """

    _default_email: str
    _default_api_key: Optional[str]

    def initiate(
        self,
        email: str = "user@example.com",
        api_key: Optional[str] = None,
    ) -> None:
        self._default_email = email
        self._default_api_key = api_key

    def run(
        self,
        gene_name: str,
        organism: str,
        email: Optional[str] = None,
        api_key: Optional[str] = None,
        max_results: int = 5,
    ) -> dict:
        gene_name = gene_name.strip()
        organism = organism.strip()

        if not gene_name:
            raise ValueError("gene_name must not be empty.")
        if not organism:
            raise ValueError("organism must not be empty.")

        resolved_organism = _ORGANISM_ALIASES.get(organism.lower(), organism)

        Entrez.email = email or self._default_email
        if api_key or self._default_api_key:
            Entrez.api_key = api_key or self._default_api_key

        # search NCBI Gene for the gene + organism
        query = f"{gene_name}[Gene Name] AND {resolved_organism}[Organism]"
        try:
            handle = Entrez.esearch(db="gene", term=query, retmax=max_results)
            search_result = Entrez.read(handle)
            handle.close()
        except (HTTPError, URLError) as exc:
            raise ValueError(f"NCBI network error during gene search: {exc}") from exc

        gene_ids = search_result.get("IdList", [])
        if not gene_ids:
            raise ValueError(
                f"No gene found for '{gene_name}' in '{resolved_organism}'. "
                "Check the gene symbol and organism name."
            )

        # try each gene ID until we get a usable mRNA/CDS
        for gene_id in gene_ids:
            result = self._fetch_cds_for_gene(gene_id, gene_name, resolved_organism)
            if result is not None:
                return result

        raise ValueError(
            f"Found {len(gene_ids)} gene ID(s) for '{gene_name}' in "
            f"'{resolved_organism}' but could not retrieve a CDS sequence. "
            "The gene may lack a RefSeq mRNA record."
        )

    def _fetch_cds_for_gene(
        self, gene_id: str, gene_name: str, organism: str
    ) -> Optional[dict]:
        # get the Gene record to find linked mRNA accessions
        try:
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
            handle.close()

            # use elink to get nucleotide records linked to this gene
            handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id, linkname="gene_nuccore_refseqrna")
            link_result = Entrez.read(handle)
            handle.close()
        except (HTTPError, URLError):
            return None

        accession_ids = []
        for link_set in link_result:
            for link in link_set.get("LinkSetDb", []):
                for link_item in link.get("Link", []):
                    accession_ids.append(link_item["Id"])

        if not accession_ids:
            # fallback: search nucleotide directly
            query = f"{gene_name}[Gene Name] AND {organism}[Organism] AND mRNA[Filter]"
            try:
                handle = Entrez.esearch(db="nucleotide", term=query, retmax=3)
                nt_result = Entrez.read(handle)
                handle.close()
                accession_ids = nt_result.get("IdList", [])
            except (HTTPError, URLError):
                return None

        if not accession_ids:
            return None

        # fetch the top accession and extract CDS
        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession_ids[0],
                rettype="gb",
                retmode="text",
            )
            record = SeqIO.read(handle, "genbank")
            handle.close()
        except (HTTPError, URLError, ValueError):
            return None

        accession = record.id

        # look for an annotated CDS feature
        for feature in record.features:
            if feature.type == "CDS":
                cds_seq = str(feature.extract(record.seq))
                product = feature.qualifiers.get("product", ["unknown"])[0]
                return {
                    "gene_name": gene_name,
                    "organism": organism,
                    "gene_id": gene_id,
                    "accession": accession,
                    "sequence": cds_seq.upper(),
                    "sequence_length": len(cds_seq),
                    "product": product,
                    "source": "CDS annotation",
                    "summary": (
                        f"Retrieved {len(cds_seq)} bp CDS for {gene_name} "
                        f"({product}) from {organism} via accession {accession}."
                    ),
                }

        # no CDS annotation — return the full record sequence as fallback
        full_seq = str(record.seq).upper()
        return {
            "gene_name": gene_name,
            "organism": organism,
            "gene_id": gene_id,
            "accession": accession,
            "sequence": full_seq,
            "sequence_length": len(full_seq),
            "product": "unknown",
            "source": "full record",
            "summary": (
                f"No CDS annotation found; returning full record ({len(full_seq)} bp) "
                f"for {gene_name} in {organism} via accession {accession}."
            ),
        }


_instance = LookupGeneSequence()
_instance.initiate()
lookup_gene_sequence = _instance.run
