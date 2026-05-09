### Karina:

The main goal of my individual scope is to bridge natural-language biological questions with structured biological knowledge sources such as:
- OLS4 Ontology Lookup Service
- NCBI Gene / NCBI E-utilities

The implemented tools enable users to:
- convert biological questions into ontology searches,
- retrieve relevant GO terms and ontology definitions,
- identify genes associated with biological processes,
- and retrieve genomic locus or sequence information.

The repository contains four primary tools:
1. **Semantic Tool**
   - Performs ontology-driven semantic search using biological natural-language queries.
2. **Annotation Tool**
   - Retrieves genes associated with Gene Ontology (GO) terms.
3. **Locus Tool** *(earlier prototype; not used in the final pipeline)*
   - Retrieves genomic coordinate and chromosome locus information.
4. **Sequence Tool** *(earlier prototype; not used in the final pipeline)*
   - Retrieves nucleotide sequence records and FASTA information.

Each tool follows the MCP architecture introduced in the course:
- backend Python implementation,
- MCP wrapper metadata,
- structured JSON outputs,
- pytest validation,
- and example prompts for LLM interaction.

<b>1. Semantic Tool</b>

- What it does:
  - Uses natural-language biology queries to identify biologically relevant ontology terms.
  - Converts human-readable biological questions into structured ontology searches.
  - Example:
    - “oxidative stress in yeast”
    - “human genes involved in blood disorders”
    
- Code files:
  - `modules/semantic_tools/semantic_wrapper.py`
    - Main backend workflow for semantic searching.
    - Parses the raw user query.
    - Detects organism names from common aliases, such as:
      - “yeast” → `Saccharomyces cerevisiae`
      - “human” → `Homo sapiens`
    - Extracts biological keywords from the query.
    - Sends ontology search requests to OLS4.
    - Parses ontology API responses.
    - Selects relevant GO terms.
    - Formats ontology information into structured JSON.
    - Includes basic API error handling.
    - Returns:
      - parsed query
      - organism
      - ontology terms
      - GO IDs
      - GO labels
      - ontology definitions
  - `modules/semantic_tools/tools/semantic_gene_search.py`
    - MCP wrapper that exposes the semantic search backend as an MCP tool using the Function Object Pattern with `initiate()` and `run()`.
  - `modules/semantic_tools/tools/semantic_gene_search.json`
    - MCP metadata definition file describing the tool schema and registration information.
- Pytests (`tests/unit/test_semantic_wrapper.py`):
  - Tests organism extraction:
    - “yeast” → `Saccharomyces cerevisiae`
  - Tests ontology term parsing:
    - “oxidative stress” → “response to oxidative stress”
  - Tests GO term retrieval:
    - verifies `GO:0006979` is returned
  - Tests structured dictionary output from the semantic wrapper
- Sample Prompt:
  - “What GO ID is related to oxidative stress in yeast?”
- Output:
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: What GO id related to oxidative stress in yeast?

  [Tool call] → semantic_gene_search
  {
    "query": "oxidative stress response yeast"
  }
  [Tool result] ← semantic_gene_search:
  {
    "parsed_query": {
      "raw_query": "oxidative stress response yeast",
      "organism": "Saccharomyces cerevisiae",
      "keywords": [
        "oxidative",
        "stress",
        "response",
        "yeast"
      ],
      "ontology_terms": [
        "response to oxidative stress"
      ]
    },
    "go_terms": [
      {
        "go_id": "GO:0006979",
        "label": "response to oxidative stress",
        "definition": "Any process that results in a change in state or activity of a cell or an organism (in terms of movement, secretion, enz\u2026 (345 chars)"
      },
      {
        "go_id": "GO:0034599",
        "label": "cellular response to oxidative stress",
        "definition": "Any process that results in a change in state or activity of a cell (in terms of movement, secretion, enzyme production,\u2026 (330 chars)"
      },
      {
        "go_id": "GO:1902882",
        "label": "regulation of response to oxidative stress",
        "definition": "Any process that modulates the frequency, rate or extent of response to oxidative stress."
      },
      {
        "go_id": "GO:1900407",
        "label": "regulation of cellular response to oxidative stress",
        "definition": "Any process that modulates the frequency, rate or extent of cellular response to oxidative stress."
      },
      {
        "go_id": "GO:1902883",
        "label": "negative regulation of response to oxidative stress",
        "definition": "Any process that stops, prevents or reduces the frequency, rate or extent of response to oxidative stress."
      }
    ]
  }

  Gemini: The Gene Ontology (GO) term for "response to oxidative stress" is **GO:0006979**.
  Other related terms include:
  *   **GO:0034599**: Cellular response to oxidative stress
  *   **GO:1902882**: Regulation of response to oxidative stress
  *   **GO:1900407**: Regulation of cellular response to oxidative stress
  ```
  </details>

<b>2. Annotation Tool</b>

- What it does:
  - This tool connects biological processes or ontology terms to possible genes.
  - Retrieves genes associated with a Gene Ontology (GO) term.
  - Uses the GO term label to search the NCBI Gene database.
  - Can optionally filter results by organism.
  - Returns structured gene annotations including:
    - gene ID
    - gene symbol
    - gene description
    - organism
- Code files:
  - `modules/annotation_tools/go_term_to_genes.py`
    - Main backend logic for annotation-based gene retrieval.
    - Defines the `GeneHit` and `LookupResult` dataclasses used for structured outputs.
    - Uses the NCBI E-utilities API to retrieve genes associated with a GO term.
    - `get_json()`:
      - sends HTTP requests, applies timeout/retry handling, returns parsed JSON responses
    - `search_gene_ids()`:
      - builds NCBI Gene search queries using: GO label, optional organism filter
      - retrieves matching NCBI Gene IDs using `esearch.fcgi`
    - `summarize_gene_ids()`:
      - retrieves gene summaries using `esummary.fcgi`
      - extracts gene symbol, description, and organism
      - converts records into `GeneHit` objects
    - `run()`:
      - performs the full workflow: search gene IDs, summarize gene results and return structured lookup output

  - `modules/annotation_tools/tools/go_term_gene_lookup.py`
    - MCP wrapper that exposes the annotation backend as an MCP tool using the Function Object Pattern with `initiate()` and `run()`.
  - `modules/annotation_tools/tools/go_term_gene_lookup.json`
    - MCP metadata definition file describing the tool schema and registration information.
- Pytests (`tests/unit/test_go_term_to_genes.py`):
  - Tests structured output shape:
    - verifies `go_id` is `GO:0006979`
    - verifies organism is `Saccharomyces cerevisiae`
    - verifies `genes` is a list
    - verifies at least one gene is returned
  - Tests plausible oxidative stress gene retrieval:
    - checks whether `YAP1` or `SOD1` appears in the returned gene symbols
- Sample Prompt:
  - “Find yeast genes associated with GO:0006979 response to oxidative stress.”
- Output:
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Find yeast genes associated with GO:0006979 response to oxidative stress.

  [Tool call] → go_term_gene_lookup
  {
    "organism": "Saccharomyces cerevisiae",
    "go_label": "response to oxidative stress",
    "go_id": "GO:0006979"
  }
  [Tool result] ← go_term_gene_lookup:
  {
    "go_id": "GO:0006979",
    "go_label": "response to oxidative stress",
    "organism": "Saccharomyces cerevisiae",
    "genes": [
      {
        "gene_id": "851752",
        "symbol": "SUP35",
        "description": "translation termination factor GTPase eRF3",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "856831",
        "symbol": "RAD51",
        "description": "recombinase RAD51",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "856709",
        "symbol": "GCN4",
        "description": "amino acid starvation-responsive transcription factor GCN4",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "852876",
        "symbol": "PMA1",
        "description": "H(+)-exporting P2-type ATPase PMA1",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "852088",
        "symbol": "SNF1",
        "description": "AMP-activated serine/threonine-protein kinase catalytic subunit SNF1",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "852433",
        "symbol": "MEC1",
        "description": "protein kinase MEC1",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "850803",
        "symbol": "HOG1",
        "description": "mitogen-activated protein kinase HOG1",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "853568",
        "symbol": "SOD1",
        "description": "superoxide dismutase SOD1",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "855005",
        "symbol": "YAP1",
        "description": "DNA-binding transcription factor YAP1",
        "organism": "Saccharomyces cerevisiae S288C"
      },
      {
        "gene_id": "856425",
        "symbol": "SLT2",
        "description": "mitogen-activated serine/threonine-protein kinase SLT2",
        "organism": "Saccharomyces cerevisiae S288C"
      }
    ]
  }

  Gemini: The following *Saccharomyces cerevisiae* genes are associated with the GO term **GO:0006979 (response to oxidative stress)**:

  | Gene Symbol | Description |
  | :--- | :--- |
  | **SUP35** | Translation termination factor GTPase eRF3 |
  | **RAD51** | Recombinase RAD51 |
  | **GCN4** | Amino acid starvation-responsive transcription factor GCN4 |
  | **PMA1** | H(+)-exporting P2-type ATPase PMA1 |
  | **SNF1** | AMP-activated serine/threonine-protein kinase catalytic subunit SNF1 |
  | **MEC1** | Protein kinase MEC1 |
  | **HOG1** | Mitogen-activated protein kinase HOG1 |
  | **SOD1** | Superoxide dismutase SOD1 |
  | **YAP1** | DNA-binding transcription factor YAP1 |
  | **SLT2** | Mitogen-activated serine/threonine-protein kinase SLT2 |
  ```
<br>
<b> 3. Locus Tool (Not used in pipeline)</b>

- What it does:
  - Retrieves genomic locus information for genes.
  - Can return chromosome accession, genomic coordinates, strand direction, and optional FASTA sequences.
  - This tool was part of an earlier version of the pipeline and is not integrated with the main pipeline.
- Code files:
  - `modules/locus_tools/gene_locus_lookup.py`
    - Main backend logic for genomic locus retrieval.
    - Uses NCBI Gene and nucleotide records to retrieve genomic coordinate information.
    - Returns structured locus metadata and optional FASTA sequences.
  - `modules/locus_tools/tools/gene_locus_lookup_tool.py`
    - MCP wrapper that exposes the locus backend as an MCP tool.
  - `modules/locus_tools/tools/gene_locus_lookup_tool.json`
    - MCP metadata definition file describing tool schema and registration information.
- Pytests (`tests/unit/gene_locus_lookup.py`)::
  - Tests locus coordinate retrieval.
  - Tests FASTA sequence retrieval.
  - Tests structured dictionary output.
- Sample Prompt:
  - “Find the genomic locus for YAP1 in yeast.”
- Output:
  <details>
  <summary>Click to expand details</summary>

  ```json
    You: Find the genomic locus for YAP1 in yeast.

    [Tool call] → gene_locus_lookup_tool
    {
      "organism": "Saccharomyces cerevisiae",
      "gene_symbol": "YAP1"
    }
    [Tool result] ← gene_locus_lookup_tool:
    {
      "query_gene_symbol": "YAP1",
      "query_gene_id": null,
      "organism": "Saccharomyces cerevisiae",
      "resolved_gene_id": "855005",
      "resolved_symbol": "YAP1",
      "gene_description": "DNA-binding transcription factor YAP1",
      "loci": [
        {
          "chr_accession": "NC_001145.3",
          "chr_loc": "XIII",
          "start_0_based": 253847,
          "stop_0_based": 255799,
          "start_1_based": 253848,
          "stop_1_based": 255800,
          "strand": "plus",
          "exon_count": 1,
          "fasta": ">NC_001145.3:253848-255800 Saccharomyces cerevisiae S288C chromosome XIII, complete sequence\nATGAGTGTGTCTACCGCCAAGAGGTCG\u2026 (2075 chars)"
        }
      ]
    }

    Gemini: The *YAP1* gene in *Saccharomyces cerevisiae* is located on **Chromosome XIII**.

    Here are the details:
    *   **Gene Symbol:** *YAP1*
    *   **Chromosome:** XIII
    *   **Genomic Coordinates (1-based):** 253,848 – 255,800
    *   **Strand:** Plus (+)
    *   **NCBI Gene ID:** 855005
  ```
  </details>

<br>
<b>4. Sequence Tool (Not used in pipeline)</b>

- What it does:
  - Retrieves nucleotide sequence records associated with genes.
  - Can optionally return FASTA sequences.
  - This tool was used in an earlier sequence retrieval workflow before the current CRISPR pipeline integration.

- Code files:
  - `modules/sequence_tools/gene_sequence_lookup.py`
    - Main backend logic for nucleotide sequence retrieval.
    - Uses NCBI Gene and linked nucleotide records to retrieve accession and sequence information.
  - `modules/sequence_tools/tools/gene_sequence_lookup_tool.py`
    - MCP wrapper that exposes the sequence lookup backend as an MCP tool.
  - `modules/sequence_tools/tools/gene_sequence_lookup_tool.json`
    - MCP metadata definition file describing tool schema and registration information.
- Pytests:
  - Tests nucleotide record retrieval.
  - Tests FASTA sequence retrieval.
  - Tests structured output formatting.
- Sample Prompt:
  - “Retrieve nucleotide sequence records for YAP1 in yeast.”
- Output:
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Retrieve nucleotide sequence records for YAP1 in yeast.

  [Tool call] → gene_sequence_lookup_tool
  {
    "gene_symbol": "YAP1",
    "organism": "Saccharomyces cerevisiae"
  }
  [Tool result] ← gene_sequence_lookup_tool:
  {
    "query_gene_symbol": "YAP1",
    "query_gene_id": null,
    "organism": "Saccharomyces cerevisiae",
    "resolved_gene_id": "855005",
    "resolved_symbol": "YAP1",
    "gene_description": "DNA-binding transcription factor YAP1",
    "nucleotide_records": [
      {
        "nucleotide_id": "330443688",
        "accession": "NC_001145",
        "title": "Saccharomyces cerevisiae S288C chromosome XIII, complete sequence",
        "fasta": ">NC_001145.3 Saccharomyces cerevisiae S288C chromosome XIII, complete sequence\nCCACACACACACCACACCCACACCACACCCACACCCACACC\u2026 (937718 chars)"
      },
      {
        "nucleotide_id": "329138949",
        "accession": "BK006946",
        "title": "TPA_inf: Saccharomyces cerevisiae S288C chromosome XIII, complete sequence",
        "fasta": ">BK006946.2 TPA_inf: Saccharomyces cerevisiae S288C chromosome XIII, complete sequence\nCCACACACACACCACACCCACACCACACCCACA\u2026 (937726 chars)"
      },
      {
        "nucleotide_id": "296147112",
        "accession": "NM_001182362",
        "title": "Saccharomyces cerevisiae S288C DNA-binding transcription factor YAP1 (YAP1), partial mRNA",
        "fasta": ">NM_001182362.1 Saccharomyces cerevisiae S288C DNA-binding transcription factor YAP1 (YAP1), partial mRNA\nATGAGTGTGTCTAC\u2026 (2088 chars)"
      }
    ]
  }

  Gemini: I have retrieved the nucleotide sequence information for the *YAP1* gene in *Saccharomyces cerevisiae*. 

  The **NM_001182362.1** record provides the partial mRNA sequence for *YAP1* (DNA-binding transcription factor), which is useful for tasks like primer design. You can also find the gene within the context of the full chromosome XIII sequences (e.g., **NC_001145**).

  *   **NCBI Gene Symbol:** *YAP1*
  *   **Accession (mRNA):** NM_001182362.1
  *   **Description:** DNA-binding transcription factor YAP1

  If you need to perform specific analyses, such as guide RNA design or sequence translation for a specific region, please let me know!
  ```
<br>

<b> Pytest Overview </b>
The pytest suite validates that each tool can successfully run its expected workflow and return structured biological outputs. Since most of these tools depend on external APIs such as OLS4 and NCBI E-utilities, the tests function mainly as integration and smoke tests rather than isolated unit tests. They verify that the API-based workflows are still reachable, that the returned data has the expected structure, and that biologically plausible results are produced.

The semantic tool test checks that a natural-language query such as `"oxidative stress in yeast"` is correctly parsed into the organism `Saccharomyces cerevisiae` and mapped to the expected Gene Ontology term `GO:0006979`.

The annotation tool tests check that the GO term lookup returns a structured result containing the correct GO ID, organism, and a non-empty gene list. They also verify that plausible oxidative stress-related yeast genes, such as `YAP1` or `SOD1`, appear in the returned results.

The locus and sequence lookup tests check that the tools can resolve the yeast gene `YAP1`, return valid NCBI-linked records, and provide structured outputs such as genomic coordinates, strand direction, nucleotide records, and optional FASTA sequences.

The pipeline smoke test runs a multi-step workflow across the semantic, annotation, and locus tools. It confirms that the tools can work together by first identifying a GO term from a natural-language query, then retrieving associated yeast genes, and finally resolving locus information for `YAP1`.

Because these tests rely on external databases, failures may sometimes reflect API downtime, rate limits, or changes in database search results rather than errors in the local code. However, they are useful for confirming that the full MCP tool workflows remain functional end-to-end.

<b>Citation:</b>

  - James McLaughlin, Josh Lagrimas, Haider Iqbal, Helen Parkinson, Henriette Harmse, OLS4: a new Ontology Lookup Service for a growing interdisciplinary knowledge ecosystem, Bioinformatics, Volume 41, Issue 5, May 2025, btaf279, https://doi.org/10.1093/bioinformatics/btaf279
  - Sayers E. A General Introduction to the E-utilities. 2009 May 26 [Updated 2022 Nov 17]. In: Entrez® Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25497/