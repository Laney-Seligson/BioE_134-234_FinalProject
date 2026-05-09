# BioE134/234 CRISPR Pipeline Project
Names: Emory Elizabeth Adelman, Jillian Ho, Karina Jin, Laney Seligson

Presentation slide: [click this](https://docs.google.com/presentation/d/1eYHWdArwtcOk_REBWtsuq_iwyRDxQa1vGYRxENukCQE/edit?usp=sharing)


---
## Table of Contents
1. [Overview](#1-overview)
2. [Project Structure](#2-project-structure)
3. [Demo](#3-demo-on-whole-pipeline)
4. [Setup](#4-setup)
5. [Individual Scope](#5-individual-scope)
    - [Karina](#karina)
    - [Emory](#emory)
    - [Laney](#laney)
    - [Jillian](#jillian)

## 1. Overview

Designing CRISPR experiments typically requires coordinating multiple tools, manual sequence handling, and careful validation, making the process time-consuming and error-prone. This project is a modular pipeline designed to streamline and automate the full CRISPR and cloning workflow, taking a user from <b>an initial idea to a lab-ready experimental plan</b>. Instead of handling each step manually, the system organizes the process into connected components that guide users through sequence selection, design, and validation in a structured way. It includes built-in gene and sequence lookup functionality, allowing users to start from a gene name, organism, or identifier and automatically retrieve the relevant DNA sequence for downstream design. The pipeline supports upstream and downstream design. Including but not limited to, identifying target sequences, designing and evaluating guide RNAs, oligo design, construction file creation, construction file validation, and lab sheet generation. Overall, the goal of this project is to create a more efficient, reproducible, and reliable approach to experimental design in synthetic biology by combining structured workflows with intelligent automation.

---

## 2. Project structure
![Pipeline Flowchart](readme_appendix/pipeline_flowchart.png)

---

## 3. Demo on whole pipeline
https://youtu.be/QgIpJ8BjLJQ 

---

## 4. Setup
You need to have at least Python 3.10 for this repos to work.

<b>1. Create a virtual environment </b>
Open a terminal in VS Code (`Terminal >> New Terminal`):
```bash
python -m venv .venv
source .venv/bin/activate        # Mac / Linux
# .venv\Scripts\activate         # Windows
```
You should see `(.venv)` at the start of your terminal line. Then install dependencies:
```bash
pip install -r requirements.txt
```

<b>2. Get your Gemini API key </b>
Get your gemini api key from **[https://aistudio.google.com/api-keys](https://aistudio.google.com/api-keys)**
In the project root folder, create a file named exactly **`.env`**.
Inside the file, paste the following.
```
GEMINI_API_KEY="paste_your_key_here"
```
<b>Never upload `.env` to GitHub. Ensure `.env` is listed in `.gitignore`. 
Run this in your project folder terminal: </b>
```
echo ".env" >> .gitignore
```

<b>3. Run the client </b>

```bash
python client_gemini.py
```
Once you see the 
```
Type a request. Ctrl-C to quit.
You: 
```
The client is ready and you can start entering prompts.

## 5. Individual Scope

### Karina:
<div style="margin-left: 20px;">
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
<details>
<summary>Click to expand details</summary>

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
</details>

<br>
<b>4. Sequence Tool (Not used in pipeline)</b>
<details>
<summary>Click to expand sample output</summary>

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
</details>
<br>
<b>Citation:</b>

  - James McLaughlin, Josh Lagrimas, Haider Iqbal, Helen Parkinson, Henriette Harmse, OLS4: a new Ontology Lookup Service for a growing interdisciplinary knowledge ecosystem, Bioinformatics, Volume 41, Issue 5, May 2025, btaf279, https://doi.org/10.1093/bioinformatics/btaf279
  - Sayers E. A General Introduction to the E-utilities. 2009 May 26 [Updated 2022 Nov 17]. In: Entrez® Programming Utilities Help [Internet]. Bethesda (MD): National Center for Biotechnology Information (US); 2010-. Available from: https://www.ncbi.nlm.nih.gov/books/NBK25497/
</div>


![CRISPR](https://img.shields.io/badge/CRISPR-Workflow-F527F2?style=for-the-badge)


## Emory Adelman — CRISPR Toolchain

My tools support CRISPR education and wet-lab protocol design. Every design decision — vector selection, guide scoring, primer logic — is backed by literature citations. The workflow begins when the user provides a gene name: the tool fetches the genomic locus or coding sequence from NCBI, selects a Cas nuclease using a heuristic I developed (the nuclease with more high-quality PAM-compatible guides wins), designs and ranks up to 10 guide sequences, then passes the best guide and all cloning variables to teammates' tools for construction file and lab sheet generation.
      



## 1. [run_full_crispr_workflow](modules/crispr_tools/tools/run_full_crispr_workflow.py)

- What it does: 
  - NEEDS gene and Organism name to successfully call tool 
  - Fetches a target sequence, genomic locus or cds from the NCBI GenBank file for the gene chosen 
  - design guide or crRNA sequence using the 10 found in the fetched sequence, meant to cover the gene broadly from PAM-based windows across the whole sequence 
  - rank guide and crispr RNAs with [rank_guides.py](modules/crispr_tools/tools/rank_guides.py), authored by Jillian and I, and find the single best one to use out of the 10 
  - off-target prediction calculated with [predict_offtargets.py](modules/crispr_tools/tools/predict_offtargets.py) 
  - design cloning oligos 
  - prepares construction file inputs 

The MCP will then stop and ask if the user wants to generate a construction file or lab sheet, then Laney and Jillian's tools are called.

**Example shown in demo video!, this output is very long** 

> **Benchmark alignment — Hsu et al. 2013:** The demo targets EMX1 in *Homo sapiens* with pX330, one of the canonical SpCas9 benchmark loci from [Hsu, P.D. et al. 2013, *Nat Biotechnol*](https://doi.org/10.1038/nbt.2647). Alignment points: NGG PAM targeting (SpCas9), mammalian guide generation via U6-driven sgRNA scaffold, pX330 BbsI CACC/AAAC sticky-end overhangs, and zero predicted off-target sites (low-off-target prioritization).

Possible Unhappy Paths:
- A user may say "Run full CRISPR workflow on genes that cause cystic fibrosis in Homo sapiens." If only a disease name is given, semantic tools are called to identify a target gene before the workflow proceeds.
- After gene resolution, the organism is occasionally re-prompted if it was inferred rather than explicitly stated.
- Vectors are not always known; the MCP supplies organism-appropriate recommendations or waits for user input before calling `design_cloning_oligos`.
- Construction file regions that cannot be determined are printed as incomplete or `"N"` — users should review those fields before submitting to wet lab.
 
- MCP Wrapper: [run_full_crispr_workflow.json](modules/crispr_tools/tools/run_full_crispr_workflow.json)
- Pytests ([tests/unit/test_run_full_crispr_workflow.py](tests/unit/test_run_full_crispr_workflow.py)): guide selection, vector prompting, workflow confirmation gate, upstream-selected gene confirmation, empty query error, per-guide score fields
- Sample Prompt: “Design a CRISPR edit targeting lacZ in E.Coli using pTargetF.”
- Citations:
  - Ran, F. Ann, Patrick D. Hsu, Jason Wright, Vineeta Agarwala, David A. Scott, and Feng Zhang. 2013. “Genome Engineering Using the CRISPR-Cas9 System.” Nature Protocols 8 (11): 2281–2308. https://doi.org/10.1038/nprot.2013.143. *(pX330 BbsI cloning protocol; basis for mammalian guide-insertion workflow)*
  - Jiang, Yihui, Bin Chen, Caifeng Duan, Bing Sun, Jian Yang, and Sheng Yang. 2015. “Multigene Editing in the Escherichia coli Genome via the CRISPR-Cas9 System.” Applied and Environmental Microbiology 81 (7): 2506–2514. https://doi.org/10.1128/AEM.04023-14. *(pTargetF/pCas9-CR4 two-plasmid E. coli editing system; basis for bacterial guide-insertion workflow)*
  - Hsu, Patrick D., David A. Scott, Joshua A. Weinstein, F. Ann Ran, Silvana Konermann, Vineeta Agarwala, Yinqing Li, et al. 2013. “DNA Targeting Specificity of RNA-Guided Cas9 Nucleases.” Nature Biotechnology 31 (9): 827–832. https://doi.org/10.1038/nbt.2647. *(canonical benchmark loci including EMX1; SpCas9 off-target specificity framework that validates this workflow's NGG PAM and low-off-target prioritization)*

## 2. [fetch_target_sequence.py](modules/crispr_tools/tools/fetch_target_sequence.py)
- What it does: 
  - Resolves a gene name or accession to a clean DNA sequence using NCBI Entrez tools
  - Closely related to Karina's gene lookup tool integrated directly into `run_full_crispr_workflow` for semi- seamless end-to-end operation
  - Uses `_utils.normalize_organism()` to canonicalize common organism aliases (e.g. "e. coli" → "Escherichia coli") before any NCBI call

  Contains 'fetch_locus_from_ncbi' which uses esearch for the gene ID,esummary for the chromosomal coordinates, and efetch to get a fasta slice from designated region of the chromosome 

  Also,'fetch_cds_from_ncbi' function makes it so we have a esearch for gene ID, elink for RefSeq mRNA accession IDs, and efetch to extract CDS from a GenBank record.  
  
  
  Future work:

  If the organism is unspecified or the gene name is ambiguous, the tool raises an error and defers to Karina's semantic search tools. A planned improvement is to replace this tool entirely with Karina's `gene_sequence_lookup_tool` and `gene_locus_*` scripts, which can infer more context. The key difference is that the current tool is stricter — it raises an error rather than returning a possibly incorrect gene — while Karina's tools can deduce more from natural language input.
  
  - MCP Wrapper: [fetch_target_sequence.json](modules/crispr_tools/tools/fetch_target_sequence.json)
  - Pytests ([tests/unit/test_fetch_target_sequence.py](tests/unit/test_fetch_target_sequence.py)): raw DNA passthrough, uppercase normalization, invalid DNA error, empty query error, required output keys
---
  - Sample Happy Prompt: “Fetch the sequence for sma-2, sma-3 and sma-4 in C.elegans.” 
  
  Notice that all genes are formatted simply and that the organism is stated. This is usually the way this tool works best. 

  - Output:

![Sample Output](https://img.shields.io/badge/expand-sample_output-F527F2?style=flat-square)
<details>
<summary>Click to expand sample output</summary>
        
      ```json
      You: Fetch sequences to edit sma-2, sma-3, and sma-4 in c.elegans             
      [Tool call] → crispr_fetch_target_sequence
      {
        "organism": "Caenorhabditis elegans",
        "query": "sma-2"
      }
      [Tool result] ← crispr_fetch_target_sequence:
      {
        "sequence": "GGACACAAGTTGACTGTCGCCGGGTGGAAACACATCTGTGAAGATGATTGCACAATTTGGAGAGACGGAGAGTTGAAAACAATGAGTAGACAGTTAGACGGTTGATGATGGGCTGGTGAA\u2026 (7410 chars)",
        "source": "ncbi",
        "resource": "sma-2",
        "organism": "Caenorhabditis elegans",
        "target_type": "genomic_locus",
        "ncbi_gene_id": "176229",
        "ncbi_accession": "NC_003281.10",
        "length": 7410,
        "note": "Fetched genomic locus from NCBI Entrez: gene 'sma-2' in 'Caenorhabditis elegans'. Nucleotide record: NC_003281.10. Verif\u2026 (173 chars)"
      }
      
      [Tool call] → crispr_fetch_target_sequence
      {
        "organism": "Caenorhabditis elegans",
        "query": "sma-3"
      }
      [Tool result] ← crispr_fetch_target_sequence:
      {
        "sequence": "TCAAAATGAACGGATTACTGCATATGCATGGTCCAGCTGTCAAAAAGCTTTTAGGTTGGAAAATTGGAGAAGACGAAGAAAAATGGTGTGAGAAGGCTGTTGAGGCATTGGTGAAAAAAT\u2026 (2995 chars)",
        "source": "ncbi",
        "resource": "sma-3",
        "organism": "Caenorhabditis elegans",
        "target_type": "genomic_locus",
        "ncbi_gene_id": "175955",
        "ncbi_accession": "NC_003281.10",
        "length": 2995,
        "note": "Fetched genomic locus from NCBI Entrez: gene 'sma-3' in 'Caenorhabditis elegans'. Nucleotide record: NC_003281.10. Verif\u2026 (173 chars)"
      }
      
      [Tool call] → crispr_fetch_target_sequence
      {
        "query": "sma-4",
        "organism": "Caenorhabditis elegans"
      }
      [Tool result] ← crispr_fetch_target_sequence:
      {
        "sequence": "TTCTCGAAAAACGTCTGAATTGGTCGAATAATGTTTCATCCAGGTATGACATCTCAGCCGTCTACTTCTAACCAAATGTATTATGATCCTTTGTATGGAGCAGAGCAAATTGTACAGGTG\u2026 (3966 chars)",
        "source": "ncbi",
        "resource": "sma-4",
        "organism": "Caenorhabditis elegans",
        "target_type": "genomic_locus",
        "ncbi_gene_id": "175815",
        "ncbi_accession": "NC_003281.10",
        "length": 3966,
        "note": "Fetched genomic locus from NCBI Entrez: gene 'sma-4' in 'Caenorhabditis elegans'. Nucleotide record: NC_003281.10. Verif\u2026 (173 chars)"
      }
      
      Gemini: I have fetched the genomic locus sequences for your target genes in *C. elegans*: **sma-2**, **sma-3**, and **sma-4**.
      
      I am ready to help you proceed with the CRISPR design workflow for these targets. Because you have multiple target genes, I will follow the multi-target workflow:
      
      1.  **System Selection**: I will run the `crispr_cas_selector` on each sequence to determine the recommended nuclease (Cas9 vs. Cas12a).
      2.  **Guide Design**: I will proceed with guide design based on that recommendation.
      3.  **Aggregation**: I will then aggregate these results to provide a unified strategy.
      ```
</details>
  
  
  - Citations: 
    - National Center for Biotechnology Information. n.d. "Entrez Programming Utilities Help (E-utilities)." U.S National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK25501/.
  
    - National Center for Biotechnology Information. n.d. "Gene Database Help." U.S. National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK3841/.
  
## 3. [cas_selector.py](modules/crispr_tools/tools/cas_selector.py)
  - What it does: 
  - Recommends Cas9 or Cas12a based on guide availability, multiplexing needs, and specificity requirements by scanning the full genomic locus (or CDS if specified) from `fetch_target_sequence`
  - For multiplexing: results across multiple genomic loci are aggregated by majority vote
  
  > **Note on AT-rich genomes:** Cas12a may provide broader targeting coverage in AT-rich genomes because TTTV PAMs can become more abundant than NGG PAMs — but this is a heuristic based on PAM density, not a universal biological rule. The selector uses guide count ratios, not GC alone, to make the recommendation (Zetsche et al. 2015; Kim et al. 2016).
  
  Possible Unhappy Paths: 
  If a person is not specific with the gene name, semantic tools will be called along with fetch sequence. The user will then have to go through a series of questions to find genes of interest. 
  
  - MCP Wrapper: [cas_selector.json](modules/crispr_tools/tools/cas_selector.json)
  - Pytests ([tests/unit/test_cas_selector.py](tests/unit/test_cas_selector.py)): GC-rich → Cas9, AT-rich → Cas12a (PAM-density heuristic), multiplexing override, empty sequence error
---
  - Sample Happy Prompt: “Which Cas nuclease should I use for sma-2, sma-3, and sma-4 in c.elegans? 
  - Output:

![Sample Output](https://img.shields.io/badge/expand-sample_output-F527F2?style=flat-square)
<details>
<summary>Click to expand sample output</summary>

        ```json 
        You: What Cas system would I use to edit sma-2, sma-3 and sma-4 in C.elegans all at once?
        
        [Gemini busy (503). Retrying in 2s...]
        
        [Gemini busy (503). Retrying in 4s...]
        
        [Tool call] → crispr_cas_selector
        {
          "seq": "GGACACAAGTTGACTGTCGCCGGGTGGAAACACATCTGTGAAGATGATTGCACAATTTGGAGAGACGGAGAGTTGAAAACAATGAGTAGACAGTTAGACGGTTGATGATGGGCTGGTGAAGCTGAGAACAATTTATGAGATACCAGGCACTTTAGAGAGGCTACAGAGGTAGAAAGAAAAAAAGGTAAAGGAGTGCGCATTTACAGGCTATTTTCAATGTTTCTTATTTATAGAATCTTCTTGTGATCTATTAGTTTTTACAATATTTCAAGAAGATTAGAATCATGGTAAGTTTTTATTTGTATTTAAGAATAAACATTCAAAAAAAATTCGAAATTTTGTAATTCAAATTTTTTAAAGAAAAAAAAAAATAGAATTCAAATTGTTGAACTCAATGAGTCTTTAAAATTTTGAAAAAAATACGAAAATCATTTTTAAAAATTTCAAAAAATGTAGCTACATTTTTTATGTTTTCCTATAAAGAATTTGTTTAAAAATTTCCCAAAATATTGTTGTGAAATAAGAAACTGAAATATCTTCTCGATATTGTCTAACCACTACCAGTTTGTAAAAATTATTATTACTTCCTCAAAACCAGGGCCTTGAGTCCTGGAACGTGTGTAATTTTTCAGCTACAGTACCCCGTGAACACCATAGACTTTGAAATTCTCAGCTTTCAAAGTTTGAACCTTCCGTTTTACATATTTTAATATTTGACGTAGTTTAAACATTTTTCCAAAATATGTGATTTCTACAATAAAATTCTTGAAAAAATCGAGACAAGAGAGACAAATAGTAATAGAGAGTACGGTAGCTGAAAAATTTCATTCGTATCAGAACCCAATTCCCAAGAATTCTTATTTTTCGAAAATAGATAAAATTGTATAACCATTTTTTATTTCACAACAAAAAACCACGAACCTTTCAGATAAACTTCGACGGAATAAAAAAGATCACGGAGAGACTGAAATGGAAACAAGGTAAGTCTCATTTTTGTCCAAAAGACCAATAGGGTGGGGAAATATTACAATCAACGGCGGAAATCTAGAATCGAAACTGAAATCGGATTTCTGACGGAAATAAGTATTTTACAGGAGACGAGGATGAGAATTGGGCAAAAAAGGCGATCGACAACTTGATGAAGAAGCTGATAAAGCATAACAAGCAGGCGTTGGAAAATTTGGAATTCGCGTTGAGGGTTCATGTGTTTTGTGTCTTTTGAATCTTTCAATTTGATCTTAAGATTTCTTTTTTATGAATAGAAATGAAATAATCTAGCAATAGACTACACATACTTTAAACTCGTTAAAAATCAATAAAAACTACATATATCCGATATTTCAATCGTAGAATAACCTCAAATCAATAACTATATCATTTCAGTGCCAAGGACAACAAAAAACCGAATGCGTGACAATTCCGCGGTCACTAGATGGACGACTTCAAATTTCGCATCGAAAGGTGATTTTAGTATAGAATAAGAGAATAATGGGGAGAAGAAGTAAACAAATTCAAAATGAAGGGGTGCCCTATTTGTTTATGAACATCAACTAGAAACAGAGACAGACTTACTTATTCGATTAATCGACCGACTAGATTTTGAAAAATGATATTAAATTTTTAAATTGTAAAATTGTGTCAACTAGTCGTTTGAACAATTTCAAAAATTGATTTGATCTTTTTTGATTTTATGGGCAAAAGTCCACGTGTCAGTCCCATTTGCATTTTTAATTTATAAATTTGCATAATTTCTACAGTACATCCTAACGATATTAGAAAAAACATAGATTTTTCTGCTAATTTTTTGATATTTTCACAAAAATTACACCAGCCAAAACTTAACCTAGATTATGAAAAATCATATTAAATTTATGAAATTTAAAAGTTTAAAAAATTTATAAACGTTAGCTGATCTTTTGTTCGATTTTCAATTTTTTTTGAAATGAAGCTTCTTCGTGAAATGCTTGTTTGAACAGCTGGAGATTATCATTTTTCACAGAAAATTCAAATAATCCACTCTGAGAAAAAATGTTAGGAACAATTCCAAATCCTATGTACTAAAAATTCAGGCACTTCCTCATGTGATATATTGCCGTGTCTATCGTTGGCCAGATCTTCAATCACATCATGAATTGAAAGCAATTGAAGATTGCCGATTCTGCTACGAATCAGGGCAAAAGGATATTTGCATTAATCCATATCATTACAAAAGAGTTCATGCAACTGGAGTACTTCCACCAGTTCTTGTACCTCGTTATTCCGAGAAACCACCACAAGAAGTGCCCCCAACTCTGGCGAAATTTCAGTTGATGGAGATGAGCGGATCAAGAATGCCGCAGAATGTCAATATGGCAAACGTTAATTTTACGTAAGTGTTTAAAAAATATATGGTATTACTGATTGAAGAAAAATTGAAATATCTGGTATTTATGAGATTATCGGGAACAAAAAAACAAAAAAGCGAAGAATAAAATAAGAAAACGAACACAGTAGGTAAGGAACGGAACTGGAAAACGGATGGGTTGATTTTGAGAAACCCGAAAACTTAAAAAAAAAACGGGGAAATGAGGGGAAAAACAGGTAACAAAGCAGCCTATCAACTATTCAAAAATTCTTTAAATTAATACATATCGAGCATTTCAGCTGCACGTTTCTTAGCATTCATCAATGGAATCATTCTCTTCAACATGAGAAGTTCTCGTGGTGGGGCAAATAGAGCTGCATTTCGAATTGATTCCATTGCTCGCTCGTAATTTTCTCTGAAAGGACTTGGATTTTTTAGTGGAACTCTGAAACTTTTTTACTTACTTCTTGATCATCGGCTCTGGCCCGTTAGTTCCCAAAAGATGAAGAACGTGAGCCGAGAAATTGAATGGATCAACTGCAATTGCTGCGTTGGCACTTTCAGTTGTCAATCTGAGTCGTTCAGCTTCTGGTTTGTGCACGTGGAGATCACTGGAAAATTAATTGTATTTAACGTAGAACTTAGAAGCCGATGGGAATTTCAATTTAAAATTCGCAGGCTGGTCCCTAATGACGTGATTGAACAAAAAATTCGCACTTTTTGATTTTTCAACTTTCCGATAGCTTAAATTTTTGCCATTTTCCTTCCTTTTAAAAGAAAATTATACTTCTTCTTTTTTCAACAACAGAAACAAAATTTATTTCAAAAATCGAATTTTATGAGTATATTCTTCGATTCTTTTTTTTTTTTTTTTGAAAAAAATCGAAAAGTCAAAAAAAATTTTTAATTTTCGATTTTTTGTTTTATTTTTTTTTTCATAAAATCTCGGCAACTGGTACCAGTCGGCGAATTTCAGAAAATTAAAATCTAAACTTACGCAGTCTGCAAATGCCTACTTGGCACCAAATATCTCGTAGTAATCCTTGCAGCCTCTTCGAGCAACTTCTTCGAGCTCTCCTCATTATTAGTCATCATATATACTTTGGACATAAGAACCAGATTGAAATGAGTAGATCTACTCTCTGCGAGTTCTTCCATTTCTCGGATTCGTCGATGAACATCAGCGTAGTCTTCGCTGGTACACAACTTCAGGAATGTAACAATCATGGAAAGATGTGGAGTGAGAGTTCCGTTCTTTTCAGCTTCTGTGAGAAGTCTATTGGCATTGTTCTGGAAAATTAAATAATTATAAATAAGAAAATTATGAATTTAAATTTACTGTATCCGCCATCATTTGTTGCTCATGATTCGTCACTGCTGTATCGAGAACACATCCATTATACTGAAACCATATTTTTAGAGAATAAAAATGATTTAATCGAATTACCACAAAAATAGATGCAAACGAATAAAGTTTGAAAAGTCTTTCCTTATCGTCTTTTGGAACTGATTCCAGGAATCTGAAAATTAGAAAAAATTAGACAAGACCTGGAAATATTTTGGAAAATCAAGAAAAAAAACTAACTTATCCAAATCTGCCATAGTATACATTAGCTCAATCTTGAGAATATCGAAAGCATCTCTCATTAGTTTTTTGCGTTCCTCATCCTCTTCTGATCCAATCAGTTCGTTAATTTTTGCTCCAACTTTATCGATATATCTCGTTGCCTCGTTAGGATCAGTACTGTAGCTTAGCATTTTTCCACGCAAGATCATTGCCGGAGGATAGCTGTTTTCCTCAGTAAGCAGATCAATCAGTTCTTCGTATTTTCCTTTTTTGAACTTTTCCAGTGCAAGTTGGTATTGTCTGAAATTTATAATTATTTTTACTGTATATCTAGAATTAGAGCTTACTTCTGGTCTTCACTGAGCACTTCTGGTTCGAATCTGGTGAGATCTGTTTTGAAAAGATCAAGAATATCCCAAGTATCACACCACTGCTGAATTCTGAACAAAGCTACTGGCTGTGGCTTAACATCTTCGATTTTGAGTGGTTTTCCAATGTGACAGTCGGACATCTAGAATAAAAAGTTGAGTTTCTTTAAGAAATTGTTAGGCTCTTACGGGTTCTATAGCAAGAAGCTTTCCAAAGAATTCGAAATTCGTCTTATTCAGGCTAGAATCCAAAGTGAACGCGGCAAGAAGATCTAAAAACAGGAATTGATACTTTTCTCGTAACATCTTCAGAAACTCACACGCCAGAGCATCCTGTTTCTTTCCGACGTCATTAAGTGCTTTTGCTGCTCTCAAGTATGCTTTGGTGTATTTTTTGTCTACTTTCAATGCAGCCATACAGTCGTTGAGGATATCAAAAGGGCTGTGTCCAACCTTTTCACGACAGGCTGCACGATTTTGGTAGAGCATTGCGACAATTTGATCAGAAGAGTTTGGTCCGGCCTTTTCGACTCCTTTTGTGAATGCTTCTAGTGCGGAATCGTATTGTTTCTCTTTAAACTTCAAATTTCCCAATGCTTTTATTTCCTCCAAGTCCCTTTCCAGCTGAAATTAGTTTTGAAATTGATATCATTGTACACTAGAAACTCACGTCTGAACCTCCAAAACTCTTATATACCAAGTATCCTACACCTGCAACTGTCGCTGCAGCTGCCACACCGATCAGCACTTTTTTCGTTTGGTCGCTGATTCCAGTCGTTTCGACCATTATATCTACAAAATAATAAATAGATAAATGATATAAGTTAAGATTAAGTTAGAAAAATTGAGAATTCAAAGAGAATGAAGAAAAATTGGTAAAATGATGTAAAATTTCAGTGATAAGAAAAGTAGGCAATATTTACGATAAACAGCATGCGGTGGAATTCCGGAGTATCGTCGCAAGTTGTAAATAAAATTGGATCCAAAATAACTTTGAGAGGGAAAAGGTTACGAAAAGAAAAAAAGGTTGAAATATTATCAATATATTTTCTAATTAGGGATTTAATAAATAATTTTCTAGAGCCAACCAATTTCATCAATACAATCCAAATGGAATTGAGGAAATGGATACTTCTCAGAAATTCGACATTCCACCAGGAGTTCCAACATGTCTCGTGCCATTCGATAAAGTTTGGGAAGAACAATTCTGGGCAACTGTTAGCTATTATGAATTGAATACAAGAGTCGGAGAACAGGTTAAAGTGTCTTCAACAACTATAACTATCGATGGATTCACTGATCCATGCATCAATGGAAGCAAAATTAGCCTGGGACTGTTTTCGAATGTAAATCGAAATGCAACTATTGAAAATACAAGAAGACATATTGGAAACGGTAAGAGCACCCACGGCCTCGCCACGAACGAATTACTTTCAAATCAAGGCCGAGCAAAACACGCGTGCGCCTTCAAGAAGTACAGTAACCTGATTTTCGGGAGATTGGGTACCGTACTTCGGCGCATAATTTTTGCAAGACTATGACTAAAAAATTGAAAATTTTCAAACTTTTTGTAATTCCAGGTGTCAAGCTAACCTACGTACGTAGCAATGGCTCACTCTTTGCTCAATGTGAAAGTGATTCAGCTATATTCGTACAATCTAGTAATTGTAATTATATCAATGGATTTCACTCAACAACTGTTGTTAAAATTGCAAATAAGTGTTCACTGAAGATTTTCGATATGGAGATTTTTAGACAGGTAAGCTAATTCACGTAATTGCAATGTTTGAGTCAAGAAAAGGTCTCGTCACGCACAAAATAACTTAAAATTGCAGTCTCACAGCAATATGTGTACCTTTAAGCAGTACGGTAGCCTATTCTTATACTTTTTGAAAAAAACTGGGTTACTGTACTTTTTAAAGGTACACACATTAAAACTGGCGAAGCCATGTTCTTGGGCAAAAAATAGTAACTGAAACTATTAGGTCTCCAAATAAGTTCCGGGTCAAAAATCATAACTTTGTTCACTGCGTATCGATTTTTATGAAACTGTGGGAATTTATGTTATCAACCATGATCTTTCATTTGACAATAGTCACAAAATTTTTTTGCCGTCCGAAGTGCCCTAACTCGGAGCCAAATTTTTCAGGCATTTTTCAGATCTCGCTTCTTTTACGTTTTGATTTGAGGTTCGTGTGCGGATTTAGCTTTGTTTAGTACACAATGTAAGAAAACAAGAAAAGTTTGGAAAAAATCCGTTAAAAAAAAAATTTTTTTTGTCGGTCGTCAAAAAATGTTCAAAAAAATTTTTGTCCAAAATTCTTGATTTTTCATACAAAAATGATGTAACCATGTGCAAAATAAATAGTTTGCACATGGTTACATCATATTTGTATAAAAAATCAAGAATTTTCGACAAAAACTTTTTTGAGCATTTTTTGACGACCGACAAAAAAATTTTTTTTTGGACGGATTTTTTTCTAAACTTTTCTTGTTTTCTTACATTGTGTACTAAACAAAGCTAAATCCGCACACAAACCTCAAATCAAAGCGTAAAAGAAGCGAGATCTGAAAAATGCCTGAAAACATTGGCTCCGAGTTAGGGCACTTCGGACGGCCAAAAAATTTTGTGACTATTGTCAAATGAAAGATCATGGTTGATAACATAAATTCCCACAGTTTCATAAAAATCGATACGCAGTGAACAAAGTTATGATTTTTGACCCGGAACTTATTTGGAGACCTAATAAGAAAAACTAAATTTTTCAGCTTCTCGAGGACTGCAGTCGCCGTGGATTCGACGCCTCATTTGATCTACAAAAAATGACATTTATTCGTATGTCATTTGTCAAAGGTTGGGGAGCTGAATATCAACGACAAGATGTCACTTCAACACCTTGTTGGATTGAAATTCATCTTCATGCTCCATTGGCTGTAAGTTCTGAAATTCTTGTTTTCTTTAGAAATAATACAATTTTCAGTGGCTCGACCGAGTTTTGTCTACAATGGGTCCAACTCCACGGCCAATTTCTTCAATTTCTTGAATCATTTTTATTATTTTTCAAATAACACTACCTTGTTGAATTGGAATCAGAAATAGAATAATAAAACAAAAAAAAAACAAACAAAGA",
          "repair_template": false,
          "organism": "Caenorhabditis elegans"
        }
        [Tool result] ← crispr_cas_selector:
        {
          "ngg_count": 422,
          "tttv_count": 489,
          "cas9_valid_guides": 143,
          "cas12a_valid_guides": 102,
          "gc_content": 0.34,
          "recommendation": "Cas12a",
          "rationale": "Guide counts are within margin: 143 valid Cas9 guides vs 102 valid Cas12a guides (ratio 1.4x < threshold 3.0x; GC 34%); \u2026 (738 chars)"
        }
        
        [Tool call] → crispr_cas_selector
        {
          "seq": "TCAAAATGAACGGATTACTGCATATGCATGGTCCAGCTGTCAAAAAGCTTTTAGGTTGGAAAATTGGAGAAGACGAAGAAAAATGGTGTGAGAAGGCTGTTGAGGCATTGGTGAAAAAATTAAAAAAGAAAAATAATGGATGTGGTTAGTCAAATATATTTTTTAAAGGAATTTATAATAATTTAGGGAGGGACATTTCTGAAAAAGGGTCTCGACGCGCTTATGAAAGTTTCAAAACTATTTTTGTGATTACCTTTAAATTTACAGTAGAAATATTTGAGCCTCCTTCCTCTTAGTCTAATTTTTCTTACGATTTCTTATGATTTTCATTTTTTTTAACTTGCAAAAAACTATAAAAACAGTAAGAAGAAAGCACCAAAAAACAAACTTAACATTTTCGAATGTATTCGAAAATACTATTGGCGCGGCGAGACCCATCCGAATGTTAAGCATATTGTTACGGATACAATATAATATTATGATTAACATTTTCCGATGACCCTTTTCGATTCATGAATCAAAGTGCACCTTATTTATAAGTTAGATCATAAATATTAAAAAAAAAGAAATTCAACCGATCCATTTCAAATTCTAAATTCCAGGAACCCTCGAAGATCTGGAATGTGTTCTAGCAAACCCATGTACAAACTCAAGATGCATTACGATAGCAAAAAGTTTAGACGGTCGGTTACAGGTATGTGGCGCCGTTGTGGTTGCTTAACGTAAATTTGTTCGCTAATCAGACAACTGTGATAACTTTCAAAAATAACACTTTGGTTTGAAATAACGACTTTTTGACGCCAAACTTCCGGCTTTTCAGGTATCCCATAAGAAAGGTCTTCCACATGTGATTTACTGTCGAGTATGGCGGTGGCCGGATATAAGTTCTCCTCATGTAAGTTAATTGGATATATTAATTTGAAATTTTTGAATGTTCATTTCAGGAATTACGTTCAATTGATACATGCTCCTACCCATATGAATCTAGTTCCAAAACCATGTACATCTGCATTAATCCATATCACTATCAGAGGTTTGTTCTGAAAGTTGTTTAAAATTGGAATACATTTCAAAAATGAATAAAATTCCAGATTATCACGACCACAAGGCTTAAATTCTTCAATGCCATCACCACAACCAATATCATCTCCAAACAGTAAGCGTTGGATCTAAAAATTGTTTAAACATATCAGAATCTTCTATTAAAAAAAAACAACTAAACTCTCTTTAAAATTCTCTAGGGAAATTAAATGAAAAATACATAAAATTGCAGCTATCTGGCAGTCTAGTGGAAGCTCAACGGCATCTTGTGCAAGCTCGCCGAGCCCGTCGGTTTTTTCAGAAGATGGAGGAGAAGTGGTTAGTTTTTTTGCAATACACATTTTTGCAAAAATTATAATTTTCAGCAAGTTCACCAGCGACCTCCACCGTTCCGTCATCCCAAATCCTGGGCCCAAATTACGTATTTTGAGTTGAATTCGAGAGTTGGAGAAGTTTTTAAAGTGAGGTTTTGATTTTTCAAAACAAAACAAAATTGGAGTTTCAGCTAGTGAATCTGAGCATTACAGTCGACGGATATACTAATCCATCAAATTCAAATACAAGGTAAGATTTTGAAAGACATTCGAAAGTAGAAATAAATTGATTTTAAAAGTGAAAAAACCAAAAAATAATTAGTTTTTTTCGATTTCGGTCGGACTTCTAAAAAATGTTCTTACAAAATATAAATATTATTAGTTTGATTTAAATTTCCTCCATTAATACGGTATAATACATATCTTCTCTATAATATTATCATAGTTATTAAAAAGTGGACACTTCGAGTTAATTAACGTGAAAATTGCTGATTCGAAGCTTTAACTTCGAGAAGTTCTGAATAATAATAATTCTTACCGGTTTCTTTATTAAAATAATCCGGAAATTTCAGAATATGTCTTGGACAATTAACCAATGTGAATCGAAATGGAACAATTGAGAACACACGGATGCATATTGGAAAAGGAATTCAGCTTGATAATAAAGAAGATCAAATGCATATTATGATTACAAATAACTCTGATATGCCAGTTTTTGTGCAATCCAAAAATACGAATCTCATGATGAACATGCCACTGGTCAAAGTCTGCCGAATACCACCGCACAGTCAATTATGCGTTTTTGAGTTTAATTTATTTTTTCAAGTGAGGTTTTTAAAATTTAAAATTATTTTCATGCTATTCTTCAACAATTGTCAAAAATTAAGATGTTGGAACAGTCGTGCAACGATAGCGATGGATTGAATGAACTCTCCAAACACTGTTTTATTCGTATTTCTTTTGTAAAGGTTTGTTCTAAATGTTTGAATTCCTTGTAGATTTATACTTGTTAAATTTATCGTGAGATTTTCAAATGAAACTTCGATAAACAAATTTTCAGGGCTGGGGAGAAGATTATCCACGTCAAGACGTCACATCCACGCCTTGCTGGCTTGAACTTCGGCTCAATGTGCCACTGGCTGTAAGTTTTTTTGAGAAGGAATCGTATCTGGAAATTTATCAAATTATCAAAAATATTCTATAAAATTTTTGTACATTTTGTAAATATCATCCAAAATTAAGCCATATTTCTATTTTGGGTACTTTTCCGATTTCTGCCCCCCCAATTTTTTAAGTTCAATATTTCAATTTCAATTCTCACGTGAATGTTTATTTCAATACTAAATAACAAAAAAGTTAACGTTTCATTAATGAATATATCGAAAACATTTTGGGGGCAGAAATCGGAAAAGTACCCTATTTTGTTTGAGTTTTCTGTGAAAATGAAATATTTTCAGTATATCGATCAAAAAATGAAGCAAACTCCAAGGACGAACTTGATGGAACCAAATTCAATGACATAATGATTTGTTAAAAACGCTTTCCTTCTCGCAACCACCTGAATTCAGAATTTACCTACTCTAATGATGATACTTGTAACATGTTTTACATTTGCTGGTTCTTTTTAATTATGAAGATTT",
          "repair_template": false,
          "organism": "Caenorhabditis elegans"
        }
        [Tool result] ← crispr_cas_selector:
        {
          "ngg_count": 183,
          "tttv_count": 207,
          "cas9_valid_guides": 48,
          "cas12a_valid_guides": 26,
          "gc_content": 0.323,
          "recommendation": "Cas12a",
          "rationale": "Guide counts are within margin: 48 valid Cas9 guides vs 26 valid Cas12a guides (ratio 1.8x < threshold 3.0x; GC 32%); us\u2026 (736 chars)"
        }
        
        [Tool call] → crispr_cas_selector
        {
          "repair_template": false,
          "seq": "TTCTCGAAAAACGTCTGAATTGGTCGAATAATGTTTCATCCAGGTATGACATCTCAGCCGTCTACTTCTAACCAAATGTATTATGATCCTTTGTATGGAGCAGAGCAAATTGTACAGGTGAGCACTATATTGTAGGAGACAGTTCAAGAGAGAGCGCGAGAGAATGAGAGAATATGAACATTTTTGGGGAACACTTATCTATTATAGTAGTCTAGTAATTGCTTCTGATTCTTAAAATACAAATAAAACTGTCTTCGTGCTTTGCATTTAGAAACTGAAGTACAAAATTCAGACAATATTAATTTTAACTATTTCAGTGCAATCCAATGGATTACCACCAGGCAAATATTTTATGCGGAATGCAGTATTTCAACAACTCGCATAATCGATATCCGGTTGGTTCAGTGATATTTGAAAACTAATTCTCTCAAAAATAAAATGGACCTCATTTTACTATTTTTCTGCATATGTAACGAAAAAATCAATAAACTTATATCCCAGTGACCAAATCGATAAGCTTTTCAAAATTTATTAAATCAAAAAATTTCAGTTATTACCTCAAATGCCACCACAGTTCACAAATGATCATCCATATGATTTTCCGAACGTGCCGACAATAAGTACTCTTGACGAAGCGTCTTCTTTTAATGGATTTTTGATACCTTCTCAACCATCTTCATATAACAATAATAACATTTCTTGTGTATTTACGCCTACACCATGTACATCATCACAAGCTTCATCCCAACCTCCACCTACACCTACAGTCAATCCAACACCGATTCCACCGAATGCTGGAGCAGTTCTCACAACTGCAATGGATTCATGTCAGCAGATTTCACATGTTTTACAATGTTATCAGCAGGTAGAAGGGGTCTTTGAAAAATAGAAAAATTTTAAGTAAAAATATAATGTTTTAAAATACTTGTTTTCCTATTAAAATACATTTCAAGTTGTATGGATTTTAAGTTTTCAAAAAAAATTGCTGAGAAAGTTGCATATTCATATCTGCTTGGGATTTATGAAAATTTTCGAACCAAAAATTTAGATTAAACTGAAACTTCTCTACATTTGAAAATAACGTACTAGACTTGAAATTTTTATTTGATTTCCTCTGATGAAATCATTGAATAAAACTCTGAAAATTTCAGAGCTCCAGAATTTTTTCTTTCTTGTAAACTTTAACATAAAACTTTAGGAAAGAAAAAACTCTGGAATTTTTACAAAATAGTTCTGAAACTACTTTTCCACTTTTATTTTTCAGAACAATTTTGAAAAAGGTATTCCAAAACTTGTCAAGCAAAACAAAAGAAGTTTTAATAATTAAAAAAAACTCCAAACTACAAAACAAACATTGTGGGTACGGCAGGCGCCAAGATAGCTGTGGAATTTGAAAACTGGAGCATTGAAAGAGAATTAAACAATTTACACGTGTTTCAGGGTGGAGAGGACAGTGATTTTGTACGGAAAGCCATCGAATCTCTTGTCAAAAAATTGAAAGACAAACGAATTGAATTGGATGCATTGATTACAGCAGTTACTTCAAATGGAAAACAACCAACTGGATGTGTTACTATTCAGGTTTGAAATAATACAACCTTTTGCTTCTGTACAGTAACCTAAAATTTTTTATTTTATGCTTGGAATTTTCGAAAATTTCAAAATTCCAGATATAAATTTCAGCAGTTTATGGGTATTTCAAAGAGAAAAACATGTTGCTTGCCAACTTTCTAAAATAACCAACAAAATTTTCAAATTTTTTTTCGTGAATTTGAAGTAAAATAATAATGTCAAAGTTTTTTCCGGTTTGAACATTTTAAATGGATTTTTGTTAGTACTTCCGATTCTTCAATTTTAAATACAGTTATCATCTTTAAAAGTTTAGAGATCCTTGGATGGACGCCTTCAAGTTGCTGGACGGAAAGGTGTTCCACATGTAGTTTACGCAAGAATATGGAGATGGCCAAAAGTTAGCAAAAATGAACTTGTGGTACATGAAAATGTTTGAACTCGGAAATAAACTGAACTATTCCAGAAATTAGTTCAATGTCAAACTTCATCAGATCATCCTGATAATATTTGCATTAATCCTTATCACTACGAAAGAGTTGTCTCAAATCGGATAACAAGTGCAGATCAATCACTTCACGTGGAGAATTCACCTATGAAGAGCGAATATTTAGGAGATGCAGGAGTTATTGATAGCTGCAGTGATTGGCCGAACACACCTCCTGATAACAATTTTAATGGTAAGAGTTGAACTCCAAAACTGTAAGTAGAGGTGGCTGCTCTCTCTCTCTGACTTTTATGCCTGCCTACGTACCTTCTAATACTTATTTGTTTGATATGGATGTTTAGTGAAGATAAAGGGTAGATAGAGGCATTTCTCATCTGCCCAAGATGAGCATGAATATATTTAATACAAAATCAACACTGAGAATTTTAGAGACCGATTTTAAATGTGACCCAATTTTTTTCAGGAGGATTTGCACCAGATCAACCTCAGCTAGTCACACCGATTATTTCTGATATTCCGATAGATCTCAATCAAATATATGTTCCAACACCTCCACAATTACTTGATAATTGGTGTTCAATCATTTATTATGAACTGGATACACCCATTGGTGAAACCTTTAAGGTATGTTTTTCTATGAAATCTGATGACTATTCATTCATGGTGCAAATCGCCTAGAAATTTTTGTGAAAGAGCGTGATTCTATCAAAAATGTGCAATGGCGTTATTAAATAATTTTAATTTTCTAAATTTTAGGTCAGCGCCCGCGATCATGGAAAAGTAATAGTAGATGGTGGAATGGATCCACATGGAGAAAATGAAGGAAGGTTATGTCTTGGAGCACTTTCAAATGTTCATCGAACAGAAGCTAGTGAAAAAGCAAGGTAATTTTAATTCAGGATCATGAAGAAAAATATAACACTGTTTTGCAGAATCCATATTGGTCGTGGTGTTGAATTAACCGCTCACGCCGATGGCAACATCTCTATCACTTCAAATTGTAAAATTTTCGTCCGATCCGGGTACCTAGACTACACTCATGGATCCGAGTACAGTAGTAAAGCACATCGATTCACTCCGAATGAATCTTCATTTACAGTATGTTCCCATTCAGGGAGTTTTACAAATGATAGAATGTTACAGGTATTCGACATCAGATGGGCCTACATGCAGATGCTCCGACGTTCTCGATCAAGCAATGAAGCAGTACGAGCACAGGCGGCAGCTGTTGCTGGATATGCACCGGTTAGTAGTAATTTGCAAGCAAAAAAACGACTCGACGCACCATATGGTGCATCATATGATCAGTGTACTGATTTTCTGCACGGACCCAATGCCATAGTTCACTTAAAAATTATTAGTGAAACTGACTAAATTTGCGACCAGGGTGATGTATTATTTCTGAGAATCAAAGATTCAATTTACCCTGATGATGCTTCGGGCATCAGCTAGCTAATTGTACACTTTTAATTCACAGACTCAAAGCATCATATATCACGAACACGAACAAGCAAAGCCTGATACCCAAGTTTAAAAAATATTTAAAAATTCAGATGTCCGTGATGCCAGCGATAATGCCTGACAGTGGAGTTGACAGGATGCGTCGTGACTTCTGTACAATTGCGATTTCATTTGTAAAAGCTTGGGGTGATGTCTACCAGAGGAAAACTATTAAGGTATGTCGACTACTCTATAAGAACTAGTCACCAATGACTTTTTTCCAGGAAACCCCATGTTGGATTGAAGTAACTTTGCACCGTCCTCTACAAATACTGGACCAGCTTCTCAAAAATAGTTCTCAATTTGGAAGTTCCTAACACTTATTTTTTATTTTTCTTCAATAATGATATCTACCTGTATTTACGTACCCCAGCAGTTTCTCATATTCATGTGCGTTGTAAAATCTCCTTTTTTGAATTTTTTCACTTTGGAAATAAAGTTTTTAGCAT",
          "organism": "Caenorhabditis elegans"
        }
        [Tool result] ← crispr_cas_selector:
        {
          "ngg_count": 244,
          "tttv_count": 244,
          "cas9_valid_guides": 96,
          "cas12a_valid_guides": 61,
          "gc_content": 0.347,
          "recommendation": "Cas12a",
          "rationale": "Guide counts are within margin: 96 valid Cas9 guides vs 61 valid Cas12a guides (ratio 1.6x < threshold 3.0x; GC 35%); us\u2026 (736 chars)"
        }
        
        Gemini: For editing all three genes (*sma-2*, *sma-3*, and *sma-4*) simultaneously in *C. elegans*, **Cas12a** is the recommended system.
        
        ### Recommendation Rationale
        I analyzed the genomic sequences for each target individually, and in all three cases, Cas12a was recommended by the selector tool based on the following:
        
        *   **Sequence Composition:** *C. elegans* genomes and these genes specifically have high AT-content (GC content ~32–35%). Cas12a recognizes a **TTTV** PAM, which is significantly more abundant and easier to target in AT-rich regions compared to Cas9's **NGG** PAM.
        *   **Multiplexing Efficiency:** Cas12a may offer broader multiplexing coverage via crRNA arrays (Zetsche et al. 2017). Since you intend to edit three genes simultaneously, Cas12a allows you to use a single crRNA array to express all guides from one transcript — a practical advantage for multi-locus editing, though efficiency can vary by context.
        
        ### Per-Gene Breakdown
        
        | Gene | Recommendation | Rationale Highlights |
        | :--- | :--- | :--- |
        | **sma-2** | **Cas12a** | AT-rich sequence (34% GC); higher abundance of TTTV PAMs (489 vs 422 NGG). |
        | **sma-3** | **Cas12a** | AT-rich sequence (32% GC); higher abundance of TTTV PAMs (207 vs 183 NGG). |
        | **sma-4** | **Cas12a** | AT-rich sequence (35% GC); TTTV/NGG PAMs well-distributed for Cas12a targeting. |
        
        ### Next Steps
        Since the selector recommends Cas12a for this multi-gene edit, would you like me to proceed with:
        
        1.  **Designing Cas12a crRNA guides** for each of these targets?
        2.  **Ranking these guides** based on efficiency and specificity?
        3.  **Assisting with the cloning strategy** for a Cas12a multiplexing vector?
        
        ```
</details>
  
  - Citations: 
    - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ian M. Slaymaker, Kira S. Makarova, Patrick Essletzbichler, Sara E. Volz, et al. 2015. "Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System." Cell 163 (3): 759–771. https://doi.org/10.1016/j.cell.2015.09.038.
  
    - Kim, Daesik, Byungkuk Min, Jungeun Kim, Seokjoong Kim, and Jin-Soo Kim. 2016. "Genome-Wide Analysis Reveals Specificities of Cpf1 Endonucleases in Human Cells." Nature Biotechnology 34 (8): 863–868. https://doi.org/10.1038/nbt.3609.
  
    - Zetsche, Bernd, Matthias Heidenreich, Prarthana Mohanraju, Ines Fedorova, Jeroen Kneppers, Ellen M. DeGennaro, Naomi Winblad, et al. 2017. "Multiplex Gene Editing by CRISPR-Cpf1 Using a Single crRNA Array." Nature Biotechnology 35 (1): 31–34. https://doi.org/10.1038/nbt.3737.
  
  
  
## 4. [design_cas9_grna](modules/crispr_tools/tools/design_cas9_grna.py)
  - What it does: 
    - Enumerates SpCas9 guide RNA candidates from NGG PAM sites across the sequence; candidate scoring applies GC content (40–65% optimum) and poly-T checks as heuristics derived from Doench et al. 2016 Rule Set 2 — downstream `rank_guides` applies weighted heuristics for final ranking (not a reproduction of the trained Doench ML model)
    - It does this with a prompt that is directly about designing guides if nuclease, genes and organism are known. 
    - For multiplexing, say on 3 different genes, each gene has to go through 1 individual workflow. 
  
    Possible Unhappy Paths: 
    If genes and organism are known: [fetch_target_sequence.py](modules/crispr_tools/tools/fetch_target_sequence.py) and `cas_selector` are used before designing guides. 
    If only the cas nuclease is known then we have to start from the semantic tools that belong to Karina. 
  
  - MCP Wrapper: [design_cas9_grna.json](modules/crispr_tools/tools/design_cas9_grna.json)
  - Pytests ([tests/unit/test_design_cas9_grna.py](tests/unit/test_design_cas9_grna.py)): 20 bp protospacer, NGG PAM, max 10 guides, no-PAM error, empty input error
---
  - Sample Happy Prompt: “Design Cas9 guides for the lacZ locus.” 
  
  - Output:

![Sample Output](https://img.shields.io/badge/expand-sample_output-F527F2?style=flat-square)
<details>
<summary>Click to expand sample output</summary>
  
    ```json
    After calling previous tools →
    [Tool result] ← crispr_design_cas9_grna:
    {
      "guides": [
        {
          "grna_sequence": "GACCAUGAUUACGGAUUCACGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "protospacer": "GACCATGATTACGGATTCAC",
          "pam_site": "TGG",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.45,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "protospacer": "CGTCGTGACTGGGAAAACCC",
          "pam_site": "TGG",
          "grna_sequence": "CGUCGUGACUGGGAAAACCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.6,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "grna_sequence": "CACAUCCCCCUUUCGCCAGCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "pam_site": "TGG",
          "protospacer": "CACATCCCCCTTTCGCCAGC",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.65,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "protospacer": "CAGCTGGCGTAATAGCGAAG",
          "pam_site": "AGG",
          "grna_sequence": "CAGCUGGCGUAAUAGCGAAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.55,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "protospacer": "CAACAGTTGCGCAGCCTGAA",
          "pam_site": "TGG",
          "grna_sequence": "CAACAGUUGCGCAGCCUGAAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.55,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "pam_site": "TGG",
          "protospacer": "TGCGCAGCCTGAATGGCGAA",
          "grna_sequence": "UGCGCAGCCUGAAUGGCGAAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.6,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "grna_sequence": "AUGGCGAAUGGCGCUUUGCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "pam_site": "TGG",
          "protospacer": "ATGGCGAATGGCGCTTTGCC",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.6,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "grna_sequence": "AUGGCGCUUUGCCUGGUUUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "protospacer": "ATGGCGCTTTGCCTGGTTTC",
          "pam_site": "CGG",
          "efficiency_score": 2,
          "efficiency_details": {
            "gc_content": 0.55,
            "gc_content_ok": true,
            "no_polyt_run": true
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 3
        },
        {
          "grna_sequence": "UCGUUUUACAACGUCGUGACGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "protospacer": "TCGTTTTACAACGTCGTGAC",
          "pam_site": "TGG",
          "efficiency_score": 1,
          "efficiency_details": {
            "gc_content": 0.45,
            "gc_content_ok": true,
            "no_polyt_run": false
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 2
        },
        {
          "grna_sequence": "CGUUUUACAACGUCGUGACUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
          "protospacer": "CGTTTTACAACGTCGTGACT",
          "pam_site": "GGG",
          "efficiency_score": 1,
          "efficiency_details": {
            "gc_content": 0.45,
            "gc_content_ok": true,
            "no_polyt_run": false
          },
          "specificity_score": 1,
          "specificity_details": {
            "high_risk_offtargets": 0,
            "medium_risk_offtargets": 0,
            "total_offtargets": 0,
            "on_target_found": true
          },
          "total_score": 2
        }
      ],
      "best_guide": {
        "grna_sequence": "GACCAUGAUUACGGAUUCACGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "GACCATGATTACGGATTCAC",
        "pam_site": "TGG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.45,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 3
      },
      "scoring_rationale": "Best guide 'GACCATGATTACGGATTCAC' scored 3 (efficiency 2, specificity 1). GC content: 45%. Off-target sites: 0 total (0 \u2026 (136 chars)",
    }
    Gemini: I have designed and ranked CRISPR-Cas9 gRNAs targeting the *lacZ* gene in *E. coli*.
    The top-ranked guide is:
    **Protospacer:** `GACCATGATTACGGATTCAC`
    **PAM:** `TGG`
    **Rationale:** This guide has a GC content of 45%, contains no poly-T runs (which would cause premature termination), and has zero predicted off-target sites in the *E. coli* genome, making it highly specific.
    Would you like to design cloning oligos for a specific vector (e.g., pTargetF, or a custom one) using this guide, or would you like me to check the editing efficiency prediction?
    ```
</details>
  
  - Citations: 
    - Jinek, Martin, Krzysztof Chylinski, Ines Fonfara, Michael Hauer, Jennifer A. Doudna, and Emmanuelle Charpentier. 2012. "A Programmable Dual-RNA–Guided DNA Endonuclease in Adaptive Bacterial Immunity." Science 337 (6096): 816–821. https://doi.org/10.1126/science.1225829.
  
    - Hsu, Patrick D., David A. Scott, Joshua A. Weinstein, F. Ann Ran, Silvana Konermann, Vineeta Agarwala, Yinqing Li, et al. 2013. "DNA Targeting Specificity of RNA-Guided Cas9 Nucleases." Nature Biotechnology 31 (9): 827–832. https://doi.org/10.1038/nbt.2647.
  
    - Doench, John G., Nicolo Fusi, Meagan Sullender, Mudra Hegde, Emma W. Vaimberg, Katherine F. Donovan, Ian Smith, et al. 2016. "Optimized sgRNA Design to Maximize Activity and Minimize Off-Target Effects of CRISPR-Cas9." Nature Biotechnology 34 (2): 184–191. https://doi.org/10.1038/nbt.3437. *(GC content optimum 40–65%, poly-T termination penalty, PAM-proximal position weights used in candidate scoring)*
  
## 5. [design_cas12a_crrna.py](modules/crispr_tools/tools/design_cas12a_crrna.py)
  - What it does: 
    - Scans a sequence for TTTV PAM sites and returns up to 10 LbCas12a crRNA candidates
    - downstream `rank_guides` applies weighted heuristics for final ranking (not a reproduction of the trained Doench ML model)
    - For multiplexing across multiple genes, each gene goes through a separate workflow.
  
    Possible Unhappy Paths: 
    If genes and organism are known: [fetch_target_sequence.py](modules/crispr_tools/tools/fetch_target_sequence.py) and `cas_selector` are used before designing guides. 
    If only the cas nuclease is known then we have to start from the semantic tools that belong to Karina. 
     
  
  - MCP Wrapper: [design_cas12a_crrna.json](modules/crispr_tools/tools/design_cas12a_crrna.json)
  - Pytests ([tests/unit/test_design_cas12a_crrna.py](tests/unit/test_design_cas12a_crrna.py)): 23 bp protospacer, TTTV PAM, max 10 guides, no-PAM error, empty input error
---
  - Sample Happy Prompt: “Give me Cas12a guides for the araB locus.”
  - Output:

![Sample Output](https://img.shields.io/badge/expand-sample_output-F527F2?style=flat-square)
<details>
<summary>Click to expand sample output</summary>
  
        ```json
        After calling previous tools →
        [Tool result] ← crispr_design_cas12a_crrna:
        {
          "result": [
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUGCAGUGAUUCUGUGCGAGCUUUG",
              "protospacer": "GCAGTGATTCTGTGCGAGCTTTG",
              "pam_site": "TTTG"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUGCGGUGGACUGCGCUACCGGUGA",
              "protospacer": "GCGGTGGACTGCGCTACCGGTGA",
              "pam_site": "TTTG"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUGAUGCCCCGAAUAACCAGUUCC",
              "protospacer": "TGATGCCCCGAATAACCAGTTCC",
              "pam_site": "TTTG"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUCCGAAAACCCGAACGCGAUGUUC",
              "protospacer": "CCGAAAACCCGAACGCGATGTTC",
              "pam_site": "TTTG"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUGCCACGCGCCGGGCAACGUUGA",
              "protospacer": "TGCCACGCGCCGGGCAACGTTGA",
              "pam_site": "TTTG"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUUCCAGCGAAUGGUUCUGGGCAA",
              "protospacer": "TTCCAGCGAATGGTTCTGGGCAA",
              "pam_site": "TTTA"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUCGGUACCACCCGCCCGCAGGAUA",
              "protospacer": "CGGTACCACCCGCCCGCAGGATA",
              "pam_site": "TTTC"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUUUGAUGAGCUGGACCCGAUCCU",
              "protospacer": "TTTGATGAGCTGGACCCGATCCT",
              "pam_site": "TTTC"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUAUGAGCUGGACCCGAUCCUCAAU",
              "protospacer": "ATGAGCTGGACCCGATCCTCAAT",
              "pam_site": "TTTG"
            },
            {
              "crrna_sequence": "AAUUUCUACUAAGUGUAGAUCCUUCCCCGCUGUUCACUGACAC",
              "protospacer": "CCTTCCCCGCTGTTCACTGACAC",
              "pam_site": "TTTG"
            }
          ]
        }
        Gemini: I have designed several Cas12a crRNA candidates targeting the *araB* gene in *E. coli*. These were identified by scanning for the TTTV PAM sites (TTTA, TTTC, TTTG) required for LbCas12a.
        | Protospacer (23 nt) | PAM | crRNA Candidate |
        | :--- | :--- | :--- |
        | `GCAGTGATTCTGTGCGAGCTTTG` | TTTG | `AAUUUCUACUAAGUGUAGAUGCAGUGAUUCUGUGCGAGCUUUG` |
        | `GCGGTGGACTGCGCTACCGGTGA` | TTTG | `AAUUUCUACUAAGUGUAGAUGCGGUGGACUGCGCUACCGGUGA` |
        | `TGATGCCCCGAATAACCAGTTCC` | TTTG | `AAUUUCUACUAAGUGUAGAUUGAUGCCCCGAAUAACCAGUUCC` |
        | `CCGAAAACCCGAACGCGATGTTC` | TTTG | `AAUUUCUACUAAGUGUAGAUCCGAAAACCCGAACGCGAUGUUC` |
        | `TGCCACGCGCCGGGCAACGTTGA` | TTTG | `AAUUUCUACUAAGUGUAGAUUGCCACGCGCCGGGCAACGUUGA` |
        | `TTCCAGCGAATGGTTCTGGGCAA` | TTTA | `AAUUUCUACUAAGUGUAGAUUUCCAGCGAAUGGUUCUGGGCAA` |
        | `CGGTACCACCCGCCCGCAGGATA` | TTTC | `AAUUUCUACUAAGUGUAGAUCGGUACCACCCGCCCGCAGGAUA` |
        | `TTTGATGAGCTGGACCCGATCCT` | TTTC | `AAUUUCUACUAAGUGUAGAUUUUGAUGAGCUGGACCCGAUCCU` |
        | `ATGAGCTGGACCCGATCCTCAAT` | TTTG | `AAUUUCUACUAAGUGUAGAUAUGAGCUGGACCCGAUCCUCAAU` |
        | `CCTTCCCCGCTGTTCACTGACAC` | TTTG | `AAUUUCUACUAAGUGUAGAUCCUUCCCCGCUGUUCACUGACAC` |
        These crRNAs include the LbCas12a direct repeat sequence (`AAUUUCUACUAAGUGUAGAU`) at the 5' end. 
        If you would like to proceed, I can help you select the best candidates based on specificity, design cloning oligos for your vector, or provide an editing efficiency prediction for a specific guide. How would          you like to continue?
        ```
</details>
    
  - Citations:
    - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ian M. Slaymaker, Kira S. Makarova, Patrick Essletzbichler, Sara E. Volz, et al. 2015. "Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System." Cell 163 (3): 759–771. https://doi.org/10.1016/j.cell.2015.09.038.
    - Fonfara, Ines, Hagen Richter, Majda Bratovič, Anaïs Le Rhun, and Emmanuelle Charpentier. 2016. "The CRISPR-Associated DNA-Cleaving Enzyme Cpf1 Also Processes Precursor CRISPR RNA." Nature 532 (7600): 517–521. https://doi.org/10.1038/nature17945.
    - Kim, Hoon Cho, Sunghyeok Ye, Daesik Kim, Young-Chang Cho, Ju-Hyun Lee, and Jin-Soo Kim. 2018. "In Vivo High-Throughput Profiling of CRISPR-Cpf1 Activity." Nature Methods 15 (1): 10–12. https://doi.org/10.1038/nbt.4061. *(Position-specific nucleotide preferences for Cas12a on-target activity; informs crRNA candidate scoring)*
    - Hu, Xianrong, Kejian Wang, Yiping Tong, et al. 2020. "Technological Advances in CRISPR-Cas-Based Plant Genome Editing." The Crop Journal 8 (3): 403–407. https://doi.org/10.1016/j.cj.2019.06.007. *(Cas12a crRNA array cloning strategy using BbsI/BsmBI into pML104/pML107 vectors)*
  
## 6. [design_cloning_oligos.py](modules/crispr_tools/tools/design_cloning_oligos.py)

- What it does:
    - Takes the vector and guide information from the guide design tools and
    - Designs annealed oligos or PCR primers into a preset vector 
    - TypeIIS (BbsI / BsmBI / BsaI annealed-oligo sticky-end ligation) — Type IIS enzymes cut outside their recognition sequence, leaving vector-specific 4-nt sticky ends (e.g. CACC/AAAC for BbsI in pX330). These overhangs are determined by the vector backbone, not user design (Ran et al. 2013; Addgene Zhang Lab cloning protocol, Sept 2015).
    - RestrictionLigation
    - Gibson
    - Golden Gate cloning (Engler et al. 2008)
  
  > **Vector provenance:** Preset vectors (pX330, pML104, pML107, pDD162, pRC11, BPK4446) were derived from published vector maps and Addgene repository records. pX330 overhangs and protocol follow Ran et al. 2013. pML104/pML107 follow Addgene depositor records. Cas12a-specific vectors follow Hu et al. 2020.
  
  Possible Unhappy Paths: 
    When left without vector, the MCP will ask for information such as 
  
    Future Work/ Not Super established:  
    - Add more vectors to Dictionary presets 
    - **Addgene API fallback:** preset vectors are handled locally; if `vector` is a numeric Addgene plasmid ID (e.g. `"42230"` or `"addgene:67639"`), the tool attempts to fetch plasmid metadata from the Addgene Developers API. Requires `ADDGENE_API_KEY` in `.env`. Dynamic results include citations in the same `{ "label", "reference", "claim" }` format as other tools. If no API key is available, use a preset vector or `vector="custom"`. This feature does not have and will not have its own JSON wrapper — it is handled entirely within `design_cloning_oligos`.
    - if organism/vector info is missing → returns `needs_user_input` with human-readable questions instead of an error
  
  - MCP Wrapper: [design_cloning_oligos.json](modules/crispr_tools/tools/design_cloning_oligos.json)
  - Pytests ([tests/unit/test_design_cloning_oligos.py](tests/unit/test_design_cloning_oligos.py)): pCRISPR E. coli produces ready oligos; pml104 enzyme is BclI-SwaI, top overhang is GATC, bottom overhang is blank, top/bottom oligo sequences verified; pml107 present with LEU2 selection; organism mismatch → needs_user_input
---
  - Sample Happy Prompt: “Design oligos to clone this protospacer into px330.” 
  - Output:

![Sample Output](https://img.shields.io/badge/expand-sample_output-F527F2?style=flat-square)
<details>
<summary>Click to expand sample output</summary>
    
        ```json
        [Tool result] ← crispr_design_cloning_oligos:
        {
          "status": "ready",
          "cloning_method": "TypeIISOligoCloning",
          "vector": "pX330",
          "enzyme": "BbsI",
          "top_overhang": "CACC",
          "bottom_overhang": "AAAC",
          "top_oligo_name": "px330_guide_top",
          "bottom_oligo_name": "px330_guide_bottom",
          "top_oligo": "CACCGCAGTGATTCTGTGCGAGCTTTG",
          "bottom_oligo": "AAACCAAAGCTCGCACAGAATCACTGC",
          "g_prepended": false,
          "final_protospacer": "GCAGTGATTCTGTGCGAGCTTTG",
          "construction_file_inputs": {
            "construct_name": "px330_araB_edit",
            "cloning_method": "TypeIISOligoCloning",
            "assembly_strategy": "TypeIISOligoCloning",
            "backbone_name": "px330",
            "backbone_sequence": "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAAC... (8484 bp)",
            "insert_name": "px330_annealed_guide_insert",
            "insert_sequence": "GCAGTGATTCTGTGCGAGCTTTG",
            "insert_forward_primer_name": "",
            "insert_forward_primer_sequence": "",
            "insert_reverse_primer_name": "",
            "insert_reverse_primer_sequence": "",
            "top_oligo_name": "px330_guide_top",
            "top_oligo_sequence": "CACCGCAGTGATTCTGTGCGAGCTTTG",
            "bottom_oligo_name": "px330_guide_bottom",
            "bottom_oligo_sequence": "AAACCAAAGCTCGCACAGAATCACTGC",
            "vector_forward_primer_name": "",
            "vector_forward_primer_sequence": "",
            "vector_reverse_primer_name": "",
            "vector_reverse_primer_sequence": "",
            "enzyme": "BbsI",
            "cell_strain": "Any mammalian",
            "selection": "Amp",
            "temperature_c": 37,
            "notes": "Vector 'pX330' \u2014 BbsI (TypeIIS) cloning. All-in-one SpCas9 + sgRNA vector.  BbsI digest leaves CACC/AAAC overhangs.  5\u2032G\u2026 (347 chars)"
          }
        }
        [Gemini busy (503). Retrying in 2s...]
        Gemini: I have designed the cloning oligos for your *araB* protospacer (`GCAGTGATTCTGTGCGAGCTTTG`) to be cloned into the pX330 vector using BbsI-mediated annealed-oligo cloning. 
        ### Cloning Oligos
        *   **Top Oligo (`px330_guide_top`):** `CACCGCAGTGATTCTGTGCGAGCTTTG`
        *   **Bottom Oligo (`px330_guide_bottom`):** `AAACCAAAGCTCGCACAGAATCACTGC`
        **Notes:**
        - A **G** was not prepended because the +1 position of your protospacer is already a **G**, which is ideal for U6 promoter transcription efficiency in pX330.
        - The top and bottom strands include the necessary `CACC` and `AAAC` overhangs for BbsI-digested pX330.
        All CRISPR design and validation steps are complete. Would you like me to generate:
        **(a) a construction file** — a structured record of the cloning workflow,
        **(b) a lab sheet** — bench-ready step-by-step protocol (requires a construction file),
        **(c) both**, or
        **(d) neither**?
        ```
</details>
  
  - Citations: 
    - Ran, F. Ann, Patrick D. Hsu, Jason Wright, Vineeta Agarwala, David A. Scott, and Feng Zhang. 2013. "Genome Engineering Using the CRISPR-Cas9 System." Nature Protocols 8 (11): 2281–2308. https://doi.org/10.1038/nprot.2013.143. *(pX330 BbsI cloning protocol; CACC/AAAC overhang logic for annealed-oligo ligation)*
  
    - Jiang, Yihui, Bin Chen, Caifeng Duan, Bing Sun, Jian Yang, and Sheng Yang. 2015. "Multigene Editing in the Escherichia coli Genome via the CRISPR-Cas9 System." Applied and Environmental Microbiology 81 (7): 2506–2514. https://doi.org/10.1128/AEM.04023-14. *(pTargetF/pCas9-CR4 two-plasmid E. coli editing system; SpeI restriction ligation)*
  
    - Sanjana, Neville E., Ophir Shalem, and Feng Zhang. 2014. "Improved Vectors and Genome-Wide Libraries for CRISPR Screening." Nature Methods 11 (8): 783–784. https://doi.org/10.1038/nmeth.3047.
  
    - Jiang, Wenyan, David Bikard, David Cox, Feng Zhang, and Luciano A. Marraffini. 2013. "RNA-Guided Editing of Bacterial Genomes Using CRISPR-Cas Systems." Nature Biotechnology 31 (3): 233–239. https://doi.org/10.1038/nbt.2508.
  
    - Gibson, Daniel G., Lei Young, Ray-Yuan Chuang, J. Craig Venter, Clyde A. Hutchison III, and Hamilton O. Smith. 2009. "Enzymatic Assembly of DNA Molecules up to Several Hundred Kilobases." Nature Methods 6 (5): 343–345. https://doi.org/10.1038/nmeth.1318.
  
    - Engler, Carola, Ramona Kandzia, and Sylvestre Marillonnet. 2008. "A One Pot, One Step, Precision Cloning Method with High Throughput Capability." PLOS ONE 3 (11): e3647. https://doi.org/10.1371/journal.pone.0003647. *(Golden Gate Assembly; BsaI and Type IIS enzyme logic)*
  
    - Laughery, Marian F., Taylor Hunter, Angela Brown, Jake Hoopes, Tiffany Ostbye, Trevor Shumaker, and John J. Wyrick. 2015. "New Vectors for Simple and Streamlined CRISPR-Cas9 Editing in *Saccharomyces cerevisiae*." Yeast 32 (12): 711–720. https://doi.org/10.1002/yea.3098.
  
    - Addgene. "Cloning of Oligos for sgRNA/shRNA." Zhang Lab protocol, September 2015. *(Standard BbsI/BsmBI annealed-oligo ligation workflow for CRISPR guide cloning; basis for TypeIISOligoCloning overhang conventions)*
  
    - Addgene. "pX330-U6-Chimeric_BB-CBh-hSpCas9 (Plasmid #42230)." https://www.addgene.org/42230/
  
    - NovoPro Bioscience. "pX330-U6-Chimeric_BB-CBh-hSpCas9 vector (Cat. No. V005940)." https://www.novoprolabs.com/vector/Vgy3dima
  
    - Hu, Xianrong, Kejian Wang, Yiping Tong, et al. 2020. "Technological Advances in CRISPR-Cas-Based Plant Genome Editing." The Crop Journal 8 (3): 403–407. https://doi.org/10.1016/j.cj.2019.06.007. *(Cas12a crRNA cloning into pML104/pML107 vectors using BbsI/BsmBI)*
  
    - New England Biolabs. "Which Restriction Enzymes Are Used in Golden Gate Assembly?" NEB FAQ. https://www.neb.com/en-us/faqs/which-restriction-enzymes-are-used-in-golden-gate-assembly. *(BbsI, BsmBI, BsaI cut outside recognition site leaving 4-nt overhangs; basis for TypeIIS enzyme selection)*
  
  ## 7.[_utils.py](modules/crispr_tools/tools/_utils.py)( NOT callable by MCP, only for support, no test files )
  
  - What it does:
    
    - Centralizes logic that was previously duplicated or missing across multiple CRISPR tools
    - `normalize_organism(name)` — converts common organism aliases ("e. coli", "human", "yeast", "c. elegans", etc.) to their canonical scientific names before any NCBI Entrez call. This fixed a bug where `fetch_target_sequence` failed for the same organism/gene combinations that `semantic_gene_search` handled correctly: the fetch tool was passing aliases like "e. coli" directly to NCBI, which rejects them.
    - `VALID_DNA` — shared set of valid DNA bases (`ATGC`), previously duplicated in [fetch_target_sequence.py](modules/crispr_tools/tools/fetch_target_sequence.py) and [run_full_crispr_workflow.py](modules/crispr_tools/tools/run_full_crispr_workflow.py)
    - `reverse_complement(seq)` — shared DNA reverse-complement function, previously duplicated as a module-level function in [design_cloning_oligos.py](modules/crispr_tools/tools/design_cloning_oligos.py) and as an instance method in [cas_selector.py](modules/crispr_tools/tools/cas_selector.py)
    
    There are definitely more scripts that could utilize a centralized _utils function. Can be a future goal. 
  
  - Files that now import from [_utils.py](modules/crispr_tools/tools/_utils.py):
    
    - [fetch_target_sequence.py](modules/crispr_tools/tools/fetch_target_sequence.py) — now normalizes organism before NCBI lookup (bug fix)
    - [lookup_gene_sequence.py](modules/crispr_tools/tools/lookup_gene_sequence.py) — replaced local `_ORGANISM_ALIASES` dict with shared `normalize_organism()`
    - [design_cloning_oligos.py](modules/crispr_tools/tools/design_cloning_oligos.py) — replaced local `_reverse_complement` and `_COMPLEMENT`; uses `normalize_organism()` in vector compatibility check
    - [cas_selector.py](modules/crispr_tools/tools/cas_selector.py) — replaced inline `_reverse_complement` instance method
    - [run_full_crispr_workflow.py](modules/crispr_tools/tools/run_full_crispr_workflow.py) — replaced local `_VALID_DNA` definition
### Laney:
<div style="margin-left: 20px;">
<b> 1. create_construction_file.py:</b> 
create_construction_file.py has two separate abilities:
<div style="margin-left: 20px;">
a) Can create a full construction file for a user who wants one for their specific inputs, and can take information from prior conversation in the same chat.

  Takes information already provided by the user in the conversation or prompts the user for missing information. Generates a construction file easily viewable inline to the user in the following format:

```
PCR         primerF         primerR         insert_template        insert_pcr
PCR         vectorF         vectorR         backbone_plasmid       vector_pcr
assembly_method   vector_pcr      insert_pcr      enzyme_or_reagent      assembled_plasmid
Transform   assembled_plasmid   competent_cells   antibiotic   temperature   final_construct

dsdna       insert_template    INSERT_DNA_SEQUENCE
plasmid     backbone_plasmid   BACKBONE_SEQUENCE
oligo       primerF        SEQUENCE_FORWARD_PRIMER
oligo       primerR        SEQUENCE_REVERSE_PRIMER
oligo       vectorF        SEQUENCE_VECTOR_FORWARD
oligo       vectorR        SEQUENCE_VECTOR_REVERSE
```

The user can provide any information they want, except that, due to API model limitations, backbone sequences are too long for the model to process. Therefore, we have curated a list of plasmids that the user can ask about, and the model will provide information on the ones they may select. The user then just needs to mention the plasmid they would like to use.

Example:
>Note: Break up the entire prompt into multiple prompts in order not to overwhelm the API tokens. Also, it may call errors as it realizes you did not give it all the needed information in one single prompt, but the API should continue to prompt you for the missing information.
```
You: I want to make a construction file. I want to do Golden Gate, my construct name is pET28a_REP24, my insert name is REP24, my insert sequence is atgaaaaatgttttaatggttactacttctcatgatgttatgggtaattctaatgaaaaaactggtttatggttatctgaattaactcatccttattattctattattgataaaaatattaatattgatattgtttctattatgggtggtgaaattcctattgatcctaattctgttgctcaagaagattattataatgataaatttttagctgatgataatttaaaaaatattatgaaaaattctacttctttacgtgatgttaatattaaagaatatgatgctattatttttgctggtggtcatggtactatgtgggattttcctaataatgctaatattcattctaaagttttagatatttatgctaaaaatggtgttattggtgctatttgtcatggtgttgctgctttaattaatgttaaagataataatggtcaaaatattattcgtgataaagaagttactggtttttctaataatgaagaaaaaattgttggtttaactgatgttgttcctttttctttagaagattctttagttgaagctggtgctaaatattcttctgcttctgaatggcaatcttatgttaaatctgattctaaaattattactgctcaaaatcctcaatctgctactgattttgctaaagctattaaacaatctttatttaat
Gemini: Should ask for the missing information
You: my backbone name is pET28a, and I want to use the plasmid pET28a for the backbone sequence, my insert forward primer name is repF, my insert forward primer sequence is ccataGGTCTCaATGAAAAATGTTTTAATGGTTACTA, my insert reverse primer name is repR, my insert reverse primer sequence is cagatGGTCTCaCGAGATTAAATAAAGATTGTTTAAT, my vector forward primer name is vecF
Gemini: Should ask for the missing information
You: my vector reverse primer name is vecR, my vector forward primer sequence is ccataGGTCTCaCTCGAGCACCACCACCACCACCACT, my vector reverse primer sequence is cagatGGTCTCaTCATGGTATATCTCCTTCTTAAAGT, my restriction enzyme is BsaI. Output my construction file
Gemini: Should output the construction file
```

b) Create a shorthand construction file from a paper
  The user can inquire about more information about what papers are available to choose from. They can then ask for a shorthand construction file to be created about that paper. This tool is then able to generate a shorthand construction file with the paper’s details.

Example:
```
You: What papers are available for me to look at?
Gemini: list of available papers. May ask you if you want to create a shorthand construction file
You: Yes, I want a shorthand construction file for Miao 2013
Gemini: Should output the shorthand construction file
```
</div>

<b>2. Create_construction_file.json:</b>
  This file is the C9 JSON wrapper for the construction file generator. It defines how the MCP framework and Gemini should understand and call the construction tool, including its name, description, inputs, outputs, examples, and execution details. The wrapper tells the model that this tool can generate a full sequence-based construction file, create a structured paper-information record, or generate a shorthand workflow from paper-derived metadata, depending on the selected input mode.

<b>3. validate_construction_file.py:</b>
This file contains the Python implementation of the construction file validation tool. Its job is to check whether a proposed cloning workflow is biologically consistent. It validates things like whether primers anneal correctly to the intended templates, whether they are oriented properly, whether a plausible amplicon can be formed, and whether the overall workflow structure is biologically reasonable. It can retrieve information from the construction file previously generated. This tool is useful after generating a construction file because it provides a second layer of error checking before the workflow is trusted for downstream use.

Supported step types and what is validated:
- **PCR**: primer annealing, orientation, and predicted amplicon length
- **GoldenGate**: BsaI overhang compatibility for circular assembly
- **Gibson**: overlap length and sequence compatibility at both junctions

Example prompt (happy path):
>Note: Have the API create a construction file first. For this example, let’s assume there are 2 PCR steps and a Gibson step in the construction file. In this case, the construction file is valid.
```
You: Validate this
Gemini: PCR step 1 passed, PCR step 2 passed, Gibson passed. Overall: Pass
```
Example prompt (sad path):
>Note: Have the API create a construction file first. For this example, let’s assume there are 2 PCR steps and a goldengate step in the construction file. In this case, the goldengate overhangs do not overlap.
```
You: Validate this
Gemini: PCR step 1 passed, PCR step 2 passed, goldengate failed due to not aligned overhangs. Overall: Fail
```

<b>4. validate_construction_file.json:</b>
This file is the JSON wrapper for the construction file validator. It describes the validator to the MCP framework by specifying the tool metadata, accepted inputs, expected outputs, and how the Python implementation should be executed. Its main purpose is to let Gemini recognize when the user is asking to check, verify, or debug a cloning design and then correctly route that request to the validation logic.

<b>5. Get_paper_info.py:</b>
This file contains the Python implementation of the paper information loader. Instead of forcing the model to guess details from a paper title alone, this tool reads a curated JSON record from the module’s data directory and returns structured metadata about a paper, such as organism, system, vectors, enzymes, assembly method, delivery method, validation methods, and major constraints. This makes the paper-based workflow much more reliable because downstream tools can use structured paper information rather than hallucinating experimental details.

<b>6. Get_paper_info.json:</b>
This file is the JSON wrapper for the paper information loader. It defines how Gemini can call the tool by providing a paper_id, and it tells the framework where the corresponding Python implementation lives. The purpose of this wrapper is to expose curated paper metadata as a callable MCP tool so that other tools, such as shorthand workflow generation, can use literature-derived information in a structured and reproducible way.

<b>Citations:</b>
  - Hall, Bradford, Andrew Cho, Advait Limaye, Kyoungin Cho, Jaspal Khillan, and Ashok B. Kulkarni. 2018. “Genome Editing in Mice Using CRISPR/Cas9 Technology.” Current Protocols in Cell Biology 81 (1): e57. https://doi.org/10.1002/cpcb.57.
  - Miao, Jin, Dongshu Guo, Jinzhe Zhang, Qingpei Huang, Genji Qin, Xin Zhang, Jianmin Wan, Hongya Gu, and Li-Jia Qu. 2013. “Targeted Mutagenesis in Rice Using CRISPR-Cas System.” Cell Research 23 (10): 1233–36. https://doi.org/10.1038/cr.2013.123.
  - Hall, Bradford, Andrew Cho, Advait Limaye, Kyoungin Cho, Jaspal Khillan, and Ashok B. Kulkarni. 2018. “Genome Editing in Mice Using CRISPR/Cas9 Technology.” Current Protocols in Cell Biology 81 (1): e57. https://doi.org/10.1002/cpcb.57.
  - “Super Simulator - SynBio Project Tutorials.” n.d. https://ucb-bioe-anderson-lab.github.io/cloning-tutorials/tools/supersimulator/#step-1-paste-your-construction-file.
  - New England Biolabs. n.d. “Golden Gate Assembly Domestication Tutorial.” https://www.neb.com/en-us/applications/cloning-and-synthetic-biology/dna-assembly-and-cloning/golden-gate-assembly?srsltid=AfmBOopg3N59s77LuEIArVnb6sf9Mxs3QQIvBPdP0-s_u6WLgFnhUjxh.
  - New England Biolabs. n.d. “Introduction to Gibson AssemblyTM.” https://www.neb.com/en-us/applications/cloning-and-synthetic-biology/dna-assembly-and-cloning/gibson-assembly?srsltid=AfmBOoojQV5Vac9DwTV5iHmNu7VOOlzIHbsOp5k2wHydVIZj5Fqrerij.
</div>

### Jillian:
<div style="margin-left: 20px;">
<b>1. predict_editing_efficiency</b>

- What it does: Predicts on-target editing efficiency BEFORE the experiment is run, using a simplified linear approximation of Doench 2016 Rule Set 2. Position-specific nucleotide weights along the protospacer, PAM context (NGG > NAG,NGGT > NGGA), GC content optimum, and poly-T penalty. Adjusts for delivery method (RNP / plasmid / lentivirus / AAV / electroporation) and editing outcome (NHEJ / HDR / base / prime editing). Returns predicted % efficiency with a ±15% confidence range, feature breakdown, and red-flag warnings.Supports both Cas9 (20 bp protospacer) and Cas12a (23 bp).

- Sample prompt: "Predict the editing efficiency for protospacer ACAGTTGGTATATTAGGAGG with NGG PAM at TGG, delivered via electroporation for an NHEJ knockout."

- Output:
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Predict the editing efficiency for protospacer ACAGTTGGTATATTAGGAGG with NGG PAM at TGG, delivered via electroporation for an NHEJ knockout.
  [Tool call] → crispr_predict_editing_efficiency
  {
    "delivery": "electroporation",
    "nuclease": "cas9",
    "protospacer": "ACAGTTGGTATATTAGGAGG",
    "outcome": "nhej",
    "pam": "TGG"
  }
  [Tool result] ← crispr_predict_editing_efficiency:
  {
    "protospacer": "ACAGTTGGTATATTAGGAGG",
    "nuclease": "cas9",
    "raw_score": 10.8,
    "on_target_efficiency_pct": 81.4,
    "confidence_range": [
      66.4,
      96.4
    ],
    "delivery": "electroporation",
    "outcome": "nhej",
    "delivery_multiplier": 0.85,
    "outcome_multiplier": 1.0,
    "feature_contributions": {
      "position_specific_nucleotides": 4.8,
      "pam_context": 5.0,
      "gc_content": 1.0,
      "polyt_penalty": 0.0
    },
    "interpretation": "Predicted NHEJ knockout efficiency: 81.4% (range 66.4-96.4%) for CAS9 delivered via electroporation. This is high \u2014 this\u2026 (232 chars)",
    "warnings": [],
    "citations": [
      {
        "label": "Doench et al. 2016, Nat Biotechnol 34:184-191",
        "reference": "https://doi.org/10.1038/nbt.3437",
        "claim": "Rule Set 2: position-specific nucleotide weights, GC content optimum 40-65%, NGGT > NGGA PAM context, and CFD off-target\u2026 (147 chars)"
      },
      {
        "label": "Hsu et al. 2013, Nat Biotechnol 31:827-832",
        "reference": "https://doi.org/10.1038/nbt.2647",
        "claim": "SpCas9 seed region is positions 1-12 from the PAM-proximal end; mismatches in the seed reduce cutting much more than PAM\u2026 (189 chars)"
      }
    ]
  }
  Gemini: The predicted NHEJ knockout efficiency for your protospacer **ACAGTTGGTATATTAGGAGG** with a **TGG** PAM, delivered via **electroporation**, is **81.4%** (confidence range: 66.4–96.4%).
  This indicates a high predicted editing efficiency. The prediction is based on the Doench 2016 Rule Set 2 model, which considers the protospacer nucleotide sequence, PAM context, and GC content, adjusted for your specific delivery method and expected outcome. No specific warnings (such as poly-T runs or weak PAM contexts) were detected for this guide.
  ```
  </details>


- Citation:
  - Doench, J. G., et al. 2016. "Optimized sgRNA design to maximize activity and minimize off-target effects of CRISPR-Cas9." Nat Biotechnol 34(2):184-191. https://doi.org/10.1038/nbt.3437
  - Hsu, P. D., et al. 2013. "DNA targeting specificity of RNA-guided Cas9 nucleases." Nat Biotechnol 31(9):827-832. https://doi.org/10.1038/nbt.2647
  - Kim, H. K., et al. 2018. "Deep learning improves prediction of CRISPR-Cpf1 guide RNA activity." Nat Biotechnol 36(3):239-241. https://doi.org/10.1038/nbt.4061
  - Lin, S., et al. 2014. "Enhanced homology-directed human genome engineering by controlled timing of CRISPR/Cas9 delivery." eLife 3:e04766. https://doi.org/10.7554/eLife.04766

<b> 2. crispr_predict_offtargets </b>
- What it does: Slides the protospacer across both strands of a reference DNA sequence, counts mismatches at each position, checks for a valid PAM (NGG for Cas9, TTTV for Cas12a), and assigns each candidate site a HIGH/MEDIUM/LOW risk score using seed-region logic (Hsu 2013 for Cas9, Kim 2016 for Cas12a). Reports an aggregate CFD-style specificity score. Handles circular references, IUPAC ambiguity codes, and both nucleases.
- Sample prompt: "Check if guide TCAGAAACCTGCCAGTTTGC has any off-target sites in pBR322 with up to 2 mismatches."
- Output:
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Check if guide TCAGAAACCTGCCAGTTTGC has any off-target sites in pBR322 with up to 2 mismatches.
  [Tool call] → crispr_predict_offtargets
  {
    "max_mismatches": 2,
    "reference": "pBR322",
    "protospacer": "TCAGAAACCTGCCAGTTTGC"
  }
  [Tool result] ← crispr_predict_offtargets:
  {
    "protospacer": "TCAGAAACCTGCCAGTTTGC",
    "nuclease": "cas9",
    "reference_length": 4361,
    "strands_scanned": 2,
    "sites_evaluated": 8684,
    "offtarget_sites": [],
    "high_risk_count": 0,
    "aggregate_offtarget_cfd": 0.0,
    "max_offtarget_cfd": 0.0,
    "specificity_summary": "No off-target sites found within 2 mismatches (CAS9), the guide appears highly specific for this reference.",
    "citations": [
      {
        "label": "Jinek et al. 2012, Science 337:816-821",
        "reference": "https://doi.org/10.1126/science.1225829",
        "claim": "SpCas9 cuts 3bp upstream of the NGG PAM producing a blunt double-strand break."
      },
      {
        "label": "Hsu et al. 2013, Nat Biotechnol 31:827-832",
        "reference": "https://doi.org/10.1038/nbt.2647",
        "claim": "SpCas9 seed region is positions 1-12 from the PAM-proximal end; mismatches in the seed reduce cutting much more than PAM\u2026 (189 chars)"
      },
      {
        "label": "Doench et al. 2016, Nat Biotechnol 34:184-191",
        "reference": "https://doi.org/10.1038/nbt.3437",
        "claim": "Rule Set 2: position-specific nucleotide weights, GC content optimum 40-65%, NGGT > NGGA PAM context, and CFD off-target\u2026 (147 chars)"
      }
    ]
  }
  Gemini: The guide **TCAGAAACCTGCCAGTTTGC** is highly specific for the **pBR322** reference sequence. 
  My analysis evaluated 8,684 potential binding sites across both strands and found **no off-target sites** with 2 or fewer mismatches. The aggregate off-target CFD score is 0, confirming that this guide should not cause unintended cutting elsewhere in the pBR322 plasmid.
  ```

- Citations:
  - Hsu, P. D., et al. 2013. "DNA targeting specificity of RNA-guided Cas9 nucleases." Nat Biotechnol 31(9):827-832. https://doi.org/10.1038/nbt.2647
  - Fu, Y., et al. 2013. "High-frequency off-target mutagenesis induced by CRISPR-Cas nucleases in human cells." Nat Biotechnol 31(9):822-826. https://doi.org/10.1038/nbt.2623
  - Zetsche, B., et al. 2015. "Cpf1 is a single RNA-guided endonuclease of a class 2 CRISPR-Cas system." Cell 163(3):759-771. https://doi.org/10.1016/j.cell.2015.09.038
  - Kim, D., et al. 2016. "Genome-wide analysis reveals specificities of Cpf1 endonucleases in human cells." Nat Biotechnol 34(8):863-868. https://doi.org/10.1038/nbt.3609
  - Jinek, M., et al. 2012. "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity." Science 337:816-821. https://doi.org/10.1126/science.1225829

<b> 3. generate_lab_sheet </b>
- What it does: Converts a structured construction file into a printable, bench-ready lab protocol with reagent volumes, thermocycler programs, and citations on every step. Renders PCR (Q5), Golden Gate, Gibson, Type IIS oligo cloning, RestrictionLigation, Transform, Pick, Miniprep, plasmid Sanger sequencing, optional CRISPR-specific delivery (RNP / electroporation / lipofection / IVT sgRNA), and a separate post-edit GENOMIC verification section that runs verify_edit on the genomic reference and prepares the ICE/TIDE submission protocol. Auto-extracts the protospacer from the construction record's parts list, then chains into predict_editing_efficiency to size the colony plan with a guide-specific number instead of a preset. Supports custom storage labels via location_overrides (boxA → "Bench-3 freezer / Rack 2", etc.) so lab sheets match the user's actual ELN/freezer system.
- Sample prompt: "Run the full CRISPR workflow for ADE2 in Saccharomyces cerevisiae using pml104 and generate a lab sheet."
- Output: 
  > The upstream workflow tools run first, then generate_lab_sheet ties everything together
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Run the full CRISPR workflow for ADE2 in Saccharomyces cerevisiae using pml104 and generate a lab sheet.
  [Tool call] → crispr_run_full_workflow
  {
    "query": "ADE2",
    "organism": "Saccharomyces cerevisiae",
    "vector": "pml104",
    "construct_name": "ADE2_edit_pml104"
  }
  [Tool result] ← crispr_run_full_workflow:
  {
    "status": "needs_user_input",
    "missing_fields": [
      "confirmed"
    ],
    "sequence_info": {
      "sequence": "ATGGATTCTAGAACAGTTGGTATATTAGGAGGGGGACAATTGGGACGTATGATTGTTGAGGCAGCAAACAGGCTCAACATTAAGACGGTAATACTAGATGCTGAAAATTCTCCTGCCAAA\u2026 (1716 chars)",
      "source": "ncbi",
      "resource": "ADE2",
      "organism": "Saccharomyces cerevisiae",
      "ncbi_gene_id": "854295",
      "ncbi_accession": "NC_001147.6",
      "length": 1716,
      "note": "Fetched from NCBI Entrez: gene 'ADE2' in 'Saccharomyces cerevisiae'. Nucleotide record: NC_001147.6. Verify the record m\u2026 (159 chars)"
    },
    "cas_recommendation": {
      "ngg_count": 142,
      "tttv_count": 74,
      "cas9_valid_guides": 68,
      "cas12a_valid_guides": 24,
      "gc_content": 0.402,
      "recommendation": "Cas9",
      "rationale": "Guide count favors Cas9: 68 valid Cas9 guides vs 24 valid Cas12a guides (ratio 2.8x >= threshold 1.5x; GC 40%). | NHEJ/k\u2026 (231 chars)"
    },
    "recommended_system": "Cas9",
    "selected_vector": "pML104",
    "vector_system": "Cas9",
    "cloning_method": "TypeIISOligoCloning",
    "enzyme": "BclI-SwaI",
    "cell_strain": "Saccharomyces cerevisiae (lithium acetate transformation)",
    "selection": "URA3 (yeast) / Ampicillin (bacteria)",
    "questions": [
      "The sequence-based Cas selector recommends Cas9 (68 Cas9 guides; 24 Cas12a guides). The selected vector 'pML104' is comp\u2026 (248 chars)"
    ],
    "options": [
      "Yes, continue with confirmed=true",
      "No, stop here or choose a different vector"
    ],
    "workflow_trace": [
      {
        "step": 1,
        "tool": "crispr_fetch_target_sequence",
        "inputs": {
          "query": "ADE2",
          "organism": "Saccharomyces cerevisiae"
        },
        "output_summary": {
          "source": "ncbi",
          "resource": "ADE2",
          "organism": "Saccharomyces cerevisiae",
          "length": 1716,
          "ncbi_gene_id": "854295",
          "ncbi_accession": "NC_001147.6"
        }
      },
      {
        "step": 2,
        "tool": "crispr_cas_selector",
        "inputs": {
          "sequence_length": 1716,
          "repair_template": false
        },
        "output_summary": {
          "recommendation": "Cas9",
          "cas9_valid_guides": 68,
          "cas12a_valid_guides": 24,
          "gc_content": 0.402
        }
      }
    ]
  }
  [Tool call] → crispr_run_full_workflow
  {
    "confirmed": true,
    "query": "ADE2",
    "vector": "pml104",
    "construct_name": "ADE2_edit_pml104",
    "organism": "Saccharomyces cerevisiae"
  }
  [Tool result] ← crispr_run_full_workflow:
  {
    "status": "ready",
    "workflow_trace": [
      {
        "step": 1,
        "tool": "crispr_fetch_target_sequence",
        "inputs": {
          "query": "ADE2",
          "organism": "Saccharomyces cerevisiae"
        },
        "output_summary": {
          "source": "ncbi",
          "resource": "ADE2",
          "organism": "Saccharomyces cerevisiae",
          "length": 1716,
          "ncbi_gene_id": "854295",
          "ncbi_accession": "NC_001147.6"
        }
      },
      {
        "step": 2,
        "tool": "crispr_cas_selector",
        "inputs": {
          "sequence_length": 1716,
          "repair_template": false
        },
        "output_summary": {
          "recommendation": "Cas9",
          "cas9_valid_guides": 68,
          "cas12a_valid_guides": 24,
          "gc_content": 0.402
        }
      },
      {
        "step": 3,
        "tool": "crispr_design_cas9_grna",
        "inputs": {
          "sequence_length": 1716
        },
        "output_summary": {
          "guide_count": 10,
          "nuclease_system": "SpCas9"
        }
      },
      {
        "step": 4,
        "tool": "crispr_rank_guides",
        "inputs": {
          "guide_count": 10,
          "reference_length": 1716,
          "nuclease": "cas9"
        },
        "output_summary": {
          "ranked_guide_count": 10,
          "best_guide": "ACAGTTGGTATATTAGGAGG",
          "best_score": 3
        }
      },
      {
        "step": 5,
        "tool": "crispr_design_cloning_oligos",
        "inputs": {
          "vector": "pml104",
          "protospacer": "ACAGTTGGTATATTAGGAGG",
          "construct_name": "ADE2_edit_pml104"
        },
        "output_summary": {
          "status": "ready",
          "vector": "pML104",
          "cloning_method": "TypeIISOligoCloning",
          "enzyme": "BclI-SwaI",
          "top_oligo": "GATCACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
          "bottom_oligo": "CTAGCTCTAAAACCCTCCTAATATACCAACTGT",
          "end_structure": "BclI-compatible 5' GATC overhang and SwaI blunt end"
        }
      },
      {
        "step": 6,
        "tool": "create_construction_file",
        "inputs": {
          "input_mode": "sequence_build",
          "assembly_strategy": "TypeIISOligoCloning"
        },
        "output_summary": {
          "construct_name": "ADE2_edit_pml104",
          "assembly_strategy": "TypeIISOligoCloning",
          "file_name": "ADE2_edit_pml104_construction.txt"
        }
      },
      {
        "step": 7,
        "tool": "validate_construction_file",
        "inputs": {
          "assembly_strategy": "TypeIISOligoCloning",
          "strict": false
        },
        "output_summary": {
          "summary": "PASS",
          "is_valid": true,
          "detail_count": 2
        }
      }
    ],
    "sequence_info": {
      "sequence": "ATGGATTCTAGAACAGTTGGTATATTAGGAGGGGGACAATTGGGACGTATGATTGTTGAGGCAGCAAACAGGCTCAACATTAAGACGGTAATACTAGATGCTGAAAATTCTCCTGCCAAA\u2026 (1716 chars)",
      "source": "ncbi",
      "resource": "ADE2",
      "organism": "Saccharomyces cerevisiae",
      "ncbi_gene_id": "854295",
      "ncbi_accession": "NC_001147.6",
      "length": 1716,
      "note": "Fetched from NCBI Entrez: gene 'ADE2' in 'Saccharomyces cerevisiae'. Nucleotide record: NC_001147.6. Verify the record m\u2026 (159 chars)"
    },
    "guides": [
      {
        "grna_sequence": "ACAGUUGGUAUAUUAGGAGGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "ACAGTTGGTATATTAGGAGG",
        "pam_site": "GGG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.4,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 3
      },
      {
        "grna_sequence": "UAUAUUAGGAGGGGGACAAUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "TATATTAGGAGGGGGACAAT",
        "pam_site": "TGG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.4,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 3
      },
      {
        "grna_sequence": "AUAUUAGGAGGGGGACAAUUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "ATATTAGGAGGGGGACAATT",
        "pam_site": "GGG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.4,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 3
      },
      {
        "grna_sequence": "AUUGGGACGUAUGAUUGUUGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "ATTGGGACGTATGATTGTTG",
        "pam_site": "AGG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.4,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 3
      },
      {
        "grna_sequence": "UGAUUGUUGAGGCAGCAAACGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "TGATTGTTGAGGCAGCAAAC",
        "pam_site": "AGG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.45,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 3
      },
      {
        "grna_sequence": "UCUAGAACAGUUGGUAUAUUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "TCTAGAACAGTTGGTATATT",
        "pam_site": "AGG",
        "efficiency_score": 1,
        "efficiency_details": {
          "gc_content": 0.3,
          "gc_content_ok": false,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 2
      },
      {
        "grna_sequence": "AGAACAGUUGGUAUAUUAGGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "AGAACAGTTGGTATATTAGG",
        "pam_site": "AGG",
        "efficiency_score": 1,
        "efficiency_details": {
          "gc_content": 0.35,
          "gc_content_ok": false,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 2
      },
      {
        "grna_sequence": "GAACAGUUGGUAUAUUAGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "GAACAGTTGGTATATTAGGA",
        "pam_site": "GGG",
        "efficiency_score": 1,
        "efficiency_details": {
          "gc_content": 0.35,
          "gc_content_ok": false,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 2
      },
      {
        "grna_sequence": "AACAGUUGGUAUAUUAGGAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "AACAGTTGGTATATTAGGAG",
        "pam_site": "GGG",
        "efficiency_score": 1,
        "efficiency_details": {
          "gc_content": 0.35,
          "gc_content_ok": false,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 2
      },
      {
        "grna_sequence": "AAACAGGCUCAACAUUAAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
        "protospacer": "AAACAGGCTCAACATTAAGA",
        "pam_site": "CGG",
        "efficiency_score": 1,
        "efficiency_details": {
          "gc_content": 0.35,
          "gc_content_ok": false,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": true
        },
        "total_score": 2
      }
    ],
    "selected_guide": {
      "grna_sequence": "ACAGUUGGUAUAUUAGGAGGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "ACAGTTGGTATATTAGGAGG",
      "pam_site": "GGG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.4,
        "gc_content_ok": true,
        "no_polyt_run": true
      },
      "specificity_score": 1,
      "specificity_details": {
        "high_risk_offtargets": 0,
        "medium_risk_offtargets": 0,
        "total_offtargets": 0,
        "on_target_found": true
      },
      "total_score": 3
    },
    "scoring_rationale": "Best guide 'ACAGTTGGTATATTAGGAGG' scored 3/3 (efficiency 2/2, specificity 1/1). GC content: 40%. Off-target sites: 0 tot\u2026 (142 chars)",
    "protospacer": "ACAGTTGGTATATTAGGAGG",
    "cloning": {
      "status": "ready",
      "cloning_method": "TypeIISOligoCloning",
      "vector": "pML104",
      "enzyme": "BclI-SwaI",
      "citations": [
        {
          "label": "Laughery et al. Yeast 2015",
          "reference": "https://doi.org/10.1002/yea.3098",
          "claim": "pML104/pML107 all-in-one S. cerevisiae CRISPR vectors; BclI-SwaI guide spacer cloning"
        },
        {
          "label": "Addgene #67638",
          "reference": "https://www.addgene.org/67638/",
          "claim": "pML104 plasmid repository record; pRSII426 backbone; URA3 selection; BclI-SwaI guide cloning"
        },
        {
          "label": "Ran et al. Nat Protoc 2013",
          "reference": "https://doi.org/10.1038/nprot.2013.143",
          "claim": "General annealed-oligo Type IIS guide cloning protocol"
        }
      ],
      "top_overhang": "GATC",
      "bottom_overhang": "",
      "top_oligo_name": "pml104_guide_top",
      "bottom_oligo_name": "pml104_guide_bottom",
      "top_oligo": "GATCACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
      "bottom_oligo": "CTAGCTCTAAAACCCTCCTAATATACCAACTGT",
      "end_structure": "BclI-compatible 5' GATC overhang and SwaI blunt end",
      "g_prepended": false,
      "final_protospacer": "ACAGTTGGTATATTAGGAGG",
      "construction_file_inputs": {
        "construct_name": "ADE2_edit_pml104",
        "cloning_method": "TypeIISOligoCloning",
        "assembly_strategy": "TypeIISOligoCloning",
        "backbone_name": "pML104",
        "backbone_sequence": "N",
        "insert_name": "pml104_annealed_guide_insert",
        "insert_sequence": "ACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
        "insert_forward_primer_name": "",
        "insert_forward_primer_sequence": "",
        "insert_reverse_primer_name": "",
        "insert_reverse_primer_sequence": "",
        "top_oligo_name": "pml104_guide_top",
        "top_oligo_sequence": "GATCACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
        "bottom_oligo_name": "pml104_guide_bottom",
        "bottom_oligo_sequence": "CTAGCTCTAAAACCCTCCTAATATACCAACTGT",
        "vector_forward_primer_name": "",
        "vector_forward_primer_sequence": "",
        "vector_reverse_primer_name": "",
        "vector_reverse_primer_sequence": "",
        "enzyme": "BclI-SwaI",
        "cell_strain": "Saccharomyces cerevisiae (lithium acetate transformation)",
        "selection": "URA3 (yeast) / Ampicillin (bacteria)",
        "temperature_c": 37,
        "notes": "Vector 'pML104' \u2014 BclI-SwaI (TypeIIS) cloning. All-in-one S. cerevisiae CRISPR vector (pRSII426 backbone, 2\u03bc, high copy)\u2026 (693 chars)"
      }
    },
    "construction_file_inputs": {
      "construct_name": "ADE2_edit_pml104",
      "cloning_method": "TypeIISOligoCloning",
      "assembly_strategy": "TypeIISOligoCloning",
      "backbone_name": "pML104",
      "backbone_sequence": "N",
      "insert_name": "pml104_annealed_guide_insert",
      "insert_sequence": "ACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
      "insert_forward_primer_name": "",
      "insert_forward_primer_sequence": "",
      "insert_reverse_primer_name": "",
      "insert_reverse_primer_sequence": "",
      "top_oligo_name": "pml104_guide_top",
      "top_oligo_sequence": "GATCACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
      "bottom_oligo_name": "pml104_guide_bottom",
      "bottom_oligo_sequence": "CTAGCTCTAAAACCCTCCTAATATACCAACTGT",
      "vector_forward_primer_name": "",
      "vector_forward_primer_sequence": "",
      "vector_reverse_primer_name": "",
      "vector_reverse_primer_sequence": "",
      "enzyme": "BclI-SwaI",
      "cell_strain": "Saccharomyces cerevisiae (lithium acetate transformation)",
      "selection": "URA3 (yeast) / Ampicillin (bacteria)",
      "temperature_c": 37,
      "notes": "Vector 'pML104' \u2014 BclI-SwaI (TypeIIS) cloning. All-in-one S. cerevisiae CRISPR vector (pRSII426 backbone, 2\u03bc, high copy)\u2026 (693 chars)"
    },
    "construction": {
      "mode": "sequence_build",
      "construct_name": "ADE2_edit_pml104",
      "assembly_strategy": "TypeIISOligoCloning",
      "file_name": "ADE2_edit_pml104_construction.txt",
      "structured_construction_file": {
        "construct_name": "ADE2_edit_pml104",
        "assembly_strategy": "TypeIISOligoCloning",
        "parts": [
          {
            "part_type": "plasmid",
            "name": "pML104",
            "sequence": "N",
            "description": "Backbone plasmid"
          },
          {
            "part_type": "dsdna",
            "name": "pml104_annealed_guide_insert",
            "sequence": "ACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
            "description": "Insert sequence"
          },
          {
            "part_type": "oligo",
            "name": "pml104_guide_top",
            "sequence": "GATCACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG",
            "description": "Top oligo for Type IIS oligo cloning"
          },
          {
            "part_type": "oligo",
            "name": "pml104_guide_bottom",
            "sequence": "CTAGCTCTAAAACCCTCCTAATATACCAACTGT",
            "description": "Bottom oligo for Type IIS oligo cloning"
          }
        ],
        "operations": [
          {
            "step_number": 1,
            "step_type": "TypeIISOligoCloning",
            "inputs": [
              "pml104_guide_top",
              "pml104_guide_bottom",
              "pML104"
            ],
            "parameters": {
              "enzyme": "BclI-SwaI",
              "overhangs": "vector-specific (e.g. CACC/AAAC)"
            },
            "output": "ADE2_edit_pml104",
            "description": ""
          },
          {
            "step_number": 2,
            "step_type": "Transform",
            "inputs": [
              "ADE2_edit_pml104"
            ],
            "parameters": {
              "cells": "Saccharomyces cerevisiae (lithium acetate transformation)",
              "selection": "URA3 (yeast) / Ampicillin (bacteria)",
              "temperature_c": 37
            },
            "output": "ADE2_edit_pml104_e",
            "description": ""
          }
        ],
        "notes": "Vector 'pML104' \u2014 BclI-SwaI (TypeIIS) cloning. All-in-one S. cerevisiae CRISPR vector (pRSII426 backbone, 2\u03bc, high copy)\u2026 (693 chars)"
      },
      "construction_file_txt": "TypeIIS       pml104_guide_toppml104_guide_bottompML104            BclI-SwaI ADE2_edit_pml104\nTransform     ADE2_edit_pm\u2026 (491 chars)",
      "text": "TypeIIS       pml104_guide_toppml104_guide_bottompML104            BclI-SwaI ADE2_edit_pml104\nTransform     ADE2_edit_pm\u2026 (491 chars)"
    },
    "validation": {
      "is_valid": true,
      "summary": "PASS",
      "details": [
        {
          "step_number": 1,
          "step_type": "TypeIISOligoCloning",
          "output_name": "ADE2_edit_pml104",
          "is_valid": null,
          "message": "TypeIISOligoCloning biological validation is not implemented in version 1."
        },
        {
          "step_number": 2,
          "step_type": "Transform",
          "output_name": "ADE2_edit_pml104_e",
          "is_valid": null,
          "message": "Transform biological validation is not implemented in version 1."
        }
      ],
      "readable_summary": "=== Construction File Validation Report ===\nConstruct: ADE2_edit_pml104\nAssembly strategy: TypeIISOligoCloning\nOverall r\u2026 (543 chars)"
    }
  }
  [Tool call] → generate_lab_sheet
  {
    "nuclease": "cas9",
    "construction_record": {
      "construct_name": "ADE2_edit_pml104",
      "operations": [
        {
          "step_type": "TypeIISOligoCloning",
          "description": "",
          "parameters": {
            "overhangs": "vector-specific (e.g. CACC/AAAC)",
            "enzyme": "BclI-SwaI"
          },
          "output": "ADE2_edit_pml104",
          "inputs": [
            "pml104_guide_top",
            "pml104_guide_bottom",
            "pML104"
          ],
          "step_number": 1
        },
        {
          "step_type": "Transform",
          "description": "",
          "parameters": {
            "cells": "Saccharomyces cerevisiae (lithium acetate transformation)",
            "selection": "URA3 (yeast) / Ampicillin (bacteria)",
            "temperature_c": 37
          },
          "inputs": [
            "ADE2_edit_pml104"
          ],
          "output": "ADE2_edit_pml104_e",
          "step_number": 2
        }
      ],
      "assembly_strategy": "TypeIISOligoCloning",
      "notes": "Vector 'pML104' \u2014 BclI-SwaI (TypeIIS) cloning. All-in-one S. cerevisiae CRISPR vector (pRSII426 backbone, 2\u03bc, high copy). Cas9 under pTDH3 (GAP) promoter; sgRNA under pSNR52 RNA Pol III promoter. Guide spacer cloned by annealing two oligos into BclI/SwaI-digested vector; the top oligo has a 5\u2032 GATC BclI-compatible overhang and the opposite end is blunt. IMPORTANT: BclI is Dam-methylation-sensitive \u2014 purify plasmid from a Dam- E. coli strain (e.g. SCS110) before digestion or BclI will not cut. Include GTTTTAGAGCTAG after the 20 nt guide to rebuild the 5\u2032 sgRNA segment. No 5\u2032G added. This insert is an sgRNA protospacer oligo (annealed oligo cloning), not a gene or multi-fragment insert.",
      "parts": [
        {
          "name": "pML104",
          "description": "Backbone plasmid",
          "part_type": "plasmid",
          "sequence": "N"
        },
        {
          "description": "Insert sequence",
          "name": "pml104_annealed_guide_insert",
          "part_type": "dsdna",
          "sequence": "ACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG"
        },
        {
          "part_type": "oligo",
          "description": "Top oligo for Type IIS oligo cloning",
          "name": "pml104_guide_top",
          "sequence": "GATCACAGTTGGTATATTAGGAGGGTTTTAGAGCTAG"
        },
        {
          "description": "Bottom oligo for Type IIS oligo cloning",
          "name": "pml104_guide_bottom",
          "part_type": "oligo",
          "sequence": "CTAGCTCTAAAACCCTCCTAATATACCAACTGT"
        }
      ]
    }
  }
  [Tool result] ← generate_lab_sheet:
  {
    "construct_name": "ADE2_edit_pml104",
    "assembly_strategy": "TypeIISOligoCloning",
    "lab_sheet_text": "A: TypeIISOligoCloning\n\nprotocol: Type IIS Oligo Cloning (anneal + digest-ligate, 10 uL)\nenzyme: BclI-SwaI\noverhangs: ve\u2026 (4441 chars)",
    "lab_sheet_tsv": "section\tlabel\tprimer1\tprimer1_location\tprimer2\tprimer2_location\ttemplate\ttemplate_location\tfragments\tproduct\tproduct_loc\u2026 (1441 chars)",
    "step_count": 6,
    "protocol_sources": [
      {
        "step": "Type IIS Oligo Cloning (anneal + digest-ligate, 10 uL)",
        "sources": [
          {
            "label": "Engler et al. 2008, PLoS ONE 3(11):e3647",
            "reference": "https://doi.org/10.1371/journal.pone.0003647",
            "claim": "Golden Gate cloning: type IIs restriction enzyme (e.g. BsaI) and T4 DNA ligase in a single reaction, with thermocycling \u2026 (155 chars)"
          },
          {
            "label": "NEB Golden Gate Assembly Kit (BsaI-HFv2), product E1601",
            "reference": "https://www.neb.com/en-us/products/e1601-neb-golden-gate-assembly-kit-bsai-hfv2",
            "claim": "NEB recommended Golden Gate reaction: 1 uL T4 ligase buffer, 0.5 uL T4 ligase, 0.5 uL BsaI-HFv2, ~1 uL DNA mix, ddH2O to\u2026 (142 chars)"
          }
        ]
      },
      {
        "step": "Chemical Transformation (heat shock)",
        "sources": [
          {
            "label": "Inoue et al. 1990, Gene 96(1):23-28",
            "reference": "https://doi.org/10.1016/0378-1119(90)90336-p",
            "claim": "Chemically competent E. coli prepared with the Inoue method achieve >1e8 cfu/ug; standard heat-shock protocol: 30 min on\u2026 (177 chars)"
          },
          {
            "label": "Sambrook & Russell 2001, Molecular Cloning: A Laboratory Manual, 3rd ed.",
            "reference": "https://www.cshlpress.com/default.tpl?action=full&--eqskudatarq=525",
            "claim": "Standard transformation, plating, and antibiotic selection protocols; rescue/recovery in SOC or 2YT for 45-60 min before\u2026 (142 chars)"
          }
        ]
      },
      {
        "step": "Plasmid Miniprep (Zymo silica column)",
        "sources": [
          {
            "label": "Zymo Research Plasmid Miniprep-Classic, product D4015",
            "reference": "https://www.zymoresearch.com/products/zyppy-plasmid-miniprep-kit",
            "claim": "Silica-column plasmid miniprep: lyse 600 uL overnight culture, neutralize, bind to column, wash, elute in 30-50 uL eluti\u2026 (139 chars)"
          }
        ]
      },
      {
        "step": "Sanger Sequencing Premix Submission",
        "sources": [
          {
            "label": "UC Berkeley DNA Sequencing Facility submission guidelines",
            "reference": "https://mcb.berkeley.edu/barker/dnaseq/",
            "claim": "Standard Sanger premix submission: ~10 uL final volume containing 200-500 ng plasmid DNA and 3.2 pmol primer (typical: 4\u2026 (171 chars)"
          }
        ]
      }
    ],
    "citations": [
      {
        "label": "Engler et al. 2008, PLoS ONE 3(11):e3647",
        "reference": "https://doi.org/10.1371/journal.pone.0003647",
        "claim": "Golden Gate cloning: type IIs restriction enzyme (e.g. BsaI) and T4 DNA ligase in a single reaction, with thermocycling \u2026 (155 chars)"
      },
      {
        "label": "NEB Golden Gate Assembly Kit (BsaI-HFv2), product E1601",
        "reference": "https://www.neb.com/en-us/products/e1601-neb-golden-gate-assembly-kit-bsai-hfv2",
        "claim": "NEB recommended Golden Gate reaction: 1 uL T4 ligase buffer, 0.5 uL T4 ligase, 0.5 uL BsaI-HFv2, ~1 uL DNA mix, ddH2O to\u2026 (142 chars)"
      },
      {
        "label": "Inoue et al. 1990, Gene 96(1):23-28",
        "reference": "https://doi.org/10.1016/0378-1119(90)90336-p",
        "claim": "Chemically competent E. coli prepared with the Inoue method achieve >1e8 cfu/ug; standard heat-shock protocol: 30 min on\u2026 (177 chars)"
      },
      {
        "label": "Sambrook & Russell 2001, Molecular Cloning: A Laboratory Manual, 3rd ed.",
        "reference": "https://www.cshlpress.com/default.tpl?action=full&--eqskudatarq=525",
        "claim": "Standard transformation, plating, and antibiotic selection protocols; rescue/recovery in SOC or 2YT for 45-60 min before\u2026 (142 chars)"
      },
      {
        "label": "Zymo Research Plasmid Miniprep-Classic, product D4015",
        "reference": "https://www.zymoresearch.com/products/zyppy-plasmid-miniprep-kit",
        "claim": "Silica-column plasmid miniprep: lyse 600 uL overnight culture, neutralize, bind to column, wash, elute in 30-50 uL eluti\u2026 (139 chars)"
      },
      {
        "label": "UC Berkeley DNA Sequencing Facility submission guidelines",
        "reference": "https://mcb.berkeley.edu/barker/dnaseq/",
        "claim": "Standard Sanger premix submission: ~10 uL final volume containing 200-500 ng plasmid DNA and 3.2 pmol primer (typical: 4\u2026 (171 chars)"
      },
      {
        "label": "Kim et al. 2014, Genome Res 24:1012-1019",
        "reference": "https://doi.org/10.1101/gr.171322.113",
        "claim": "Cas9 RNP delivery to mammalian cells achieves ~50-80% editing efficiency, higher than plasmid transfection."
      },
      {
        "label": "Paquet et al. 2016, Nature 533:125-129",
        "reference": "https://doi.org/10.1038/nature17664",
        "claim": "HDR knockin efficiency is typically ~5-10% in mammalian cells, ~1/10 of NHEJ cutting efficiency."
      },
      {
        "label": "Jiang et al. 2013, Nat Biotechnol 31:233-239",
        "reference": "https://doi.org/10.1038/nbt.2508",
        "claim": "CRISPR-Cas9 editing in E. coli typically reaches 30-70% efficiency with pCRISPR/lambda Red recombineering systems."
      }
    ],
    "verify_before_use": "Volumes and conditions in this lab sheet are transcribed from canonical sources (NEB, Zymo, original method papers) into\u2026 (359 chars)",
    "colony_plan": [
      {
        "transform_index": 0,
        "colonies_to_pick": 5,
        "labels": [
          "A1A",
          "A1B",
          "A1C",
          "A1D",
          "A1E"
        ],
        "efficiency_used": 0.527,
        "preset": null,
        "rationale": "Pick 5 colonies to be 95% confident of recovering at least 1 edited clone(s) at 53% editing efficiency. Recommended safe\u2026 (305 chars)"
      }
    ],
    "protocols_io_links": [
      {
        "protocol": "Type IIS Oligo Cloning (anneal + digest-ligate, 10 uL)",
        "search_url": "https://www.protocols.io/search?q=Type+IIS+Oligo+Cloning+%28anneal+%2B+digest-ligate%2C+10+uL%29"
      },
      {
        "protocol": "Chemical Transformation (heat shock)",
        "search_url": "https://www.protocols.io/search?q=Chemical+Transformation+%28heat+shock%29"
      },
      {
        "protocol": "Plasmid Miniprep (Zymo silica column)",
        "search_url": "https://www.protocols.io/search?q=Plasmid+Miniprep+%28Zymo+silica+column%29"
      },
      {
        "protocol": "Sanger Sequencing Premix Submission",
        "search_url": "https://www.protocols.io/search?q=Sanger+Sequencing+Premix+Submission"
      }
    ],
    "predicted_editing_efficiency": {
      "on_target_efficiency_pct": 52.7,
      "confidence_range": [
        37.7,
        67.7
      ],
      "interpretation": "Predicted NHEJ knockout efficiency: 52.7% (range 37.7-67.7%) for CAS9 delivered via plasmid. This is moderate \u2014 usable, \u2026 (228 chars)",
      "warnings": [],
      "delivery": "plasmid",
      "outcome": "nhej"
    },
    "verify_edit_summary": null
  }
  [Gemini busy (503). Retrying in 2s...]
  Gemini: The full CRISPR workflow for targeting *ADE2* in *Saccharomyces cerevisiae* with the pML104 vector is complete.
  ### 1. Design Summary
  *   **Target Gene:** *ADE2* (NCBI Accession NC_001147.6)
  *   **System:** SpCas9
  *   **Vector:** pML104 (BclI-SwaI Type IIS cloning)
  *   **Best Guide:** `ACAGTTGGTATATTAGGAGG`
  *   **Predicted Editing Efficiency:** 52.7% (NHEJ knockout, Cas9 plasmid)
  *   **Off-target Assessment:** 0 total predicted off-targets (0 HIGH, 0 MEDIUM risk).
  ### 2. Cloning Protocol (Lab Sheet)
  Below is your bench-ready lab protocol.
  ```text
  A: TypeIISOligoCloning
  protocol: Type IIS Oligo Cloning (anneal + digest-ligate, 10 uL)
  enzyme: BclI-SwaI
  overhangs: vector-specific (e.g. CACC/AAAC)
  reaction:
  1 uL 100 uM top oligo
  1 uL 100 uM bottom oligo
  1 uL T4 DNA Ligase Buffer (10X) — also serves as anneal buffer
  1 uL Vector (50 ng, undigested)
  0.5 uL Type IIS enzyme (BbsI-HF / BsmBI-v2)
  0.5 uL T4 DNA Ligase
  5 uL ddH2O
  program: Phase 1 - Anneal oligos: combine 1 uL top + 1 uL bottom + 1 uL T4 ligase buffer + 7 uL ddH2O in a PCR tube. Heat 95C 5 min, ramp -0.1C/sec to 25C (or boil and slow-cool in a beaker). Dilute 1:200 in ddH2O before use.
  Phase 2 - Digest-ligate: in a fresh tube combine the recipe above using 1 uL of the diluted oligo duplex in place of the anneal mix. Run [37C 5 min, 16C 5 min] x 30; 60C 5 min; 80C 5 min; 4C hold.
  source:
  dna                  location
  pml104_guide_top     oligos1/A1
  pml104_guide_bottom  oligos1/B1
  pML104               boxA/A1
  samples:
  label  fragments                                             product
  A      pml104_guide_top+pml104_guide_bottom duplex / pML104  ADE2_edit_pml104
  destination: thermocycler1A
  notes:
  Use BbsI-HF (NEB R3539) for pX330/pSpCas9-style Cas9 vectors and BsmBI-v2 (NEB R0739) for Lenti-Guide-Puro / pCRISPR. Vector overhangs (e.g. CACC/AAAC) come from design_cloning_oligos. Resulting plasmid carries the protospacer in-frame with the sgRNA scaffold.
  A: Transform
  samples:
  label  product           strain                                                    antibiotic                            incubate
  A1     ADE2_edit_pml104  Saccharomyces cerevisiae (lithium acetate transformation)  URA3 (yeast) / Ampicillin (bacteria)  37°C
  A: Pick
  samples:
  source  product           strain                                                    antibiotic                            incubate  number  labels
  A1      ADE2_edit_pml104  Saccharomyces cerevisiae (lithium acetate transformation)  URA3 (yeast) / Ampicillin (bacteria)  37°C      5       A1A, A1B, A1C, A1D, A1E
  colonies to pick: 5 (efficiency 53%)
  rationale: Pick 5 colonies to be 95% confident of recovering at least 1 edited clone(s) at 53% editing efficiency. Recommended safety total: 8 colonies.
  A: Miniprep
  (See lab sheet above for details)
  ```
  ### 3. Next Steps
  *   **Customization:** The lab sheet uses placeholder storage labels (e.g., `boxA/A1`). If you want to use real lab freezer/rack labels, provide me with a mapping (e.g., `{'boxA': 'Freezer-1 Rack-2', 'oligos1': 'Oligo-Box-A'}`) and I will regenerate the sheet.
  *   **Genomic Verification:** Once you have your edited colonies, do not rely solely on the plasmid sequence. You should PCR-amplify the genomic *ADE2* locus and Sanger sequence it. 
  *   **ICE/TIDE Analysis:** Once you have your Sanger results, paste the trace data and ICE/TIDE metrics (editing %, R²) back to me for interpretation.
  Would you like to design sequencing primers for the *ADE2* locus verification now?
  ```

- Citations:
  - Engler, C., Kandzia, R., Marillonnet, S. 2008. "A One Pot, One Step, Precision Cloning Method with High Throughput Capability." PLoS ONE 3(11):e3647. https://doi.org/10.1371/journal.pone.0003647
  - Gibson, D. G., et al. 2009. "Enzymatic assembly of DNA molecules up to several hundred kilobases." Nat Methods 6:343-345. https://doi.org/10.1038/nmeth.1318
  - Inoue, H., et al. 1990. "High efficiency transformation of Escherichia coli with plasmids." Gene 96(1):23-28. https://doi.org/10.1016/0378-1119(90)90336-p
  - Sambrook, J., Russell, D. W. 2001. Molecular Cloning: A Laboratory Manual, 3rd ed. Cold Spring Harbor Laboratory Press.
  - NEB Q5 High-Fidelity DNA Polymerase, product M0491. https://www.neb.com/en-us/products/m0491-q5-high-fidelity-dna-polymerase
  - NEB Golden Gate Assembly Kit (BsaI-HFv2), product E1601. https://www.neb.com/en-us/products/e1601-neb-golden-gate-assembly-kit-bsai-hfv2
  - NEB Gibson Assembly Master Mix, product E2611. https://www.neb.com/en-us/products/e2611-gibson-assembly-master-mix
  - Zymo Research Plasmid Miniprep-Classic, product D4015. https://www.zymoresearch.com/products/zyppy-plasmid-miniprep-kit
  - Laughery, M. F., et al. 2015. "New vectors for simple and streamlined CRISPR-Cas9 genome editing in Saccharomyces cerevisiae." Yeast 32(12):711-720. https://doi.org/10.1002/yea.3098

<b> 4. colony_calculator </b>
- What it does: Computes how many colonies to pick to recover a desired number of edited clones at a target confidence level, given the editing efficiency. Uses the binomial distribution: each colony is an independent Bernoulli trial. Includes 7 published-benchmark presets (cas9_yeast, cas9_ecoli,
cas9_plasmid_mammalian, cas9_rnp_mammalian, cas12a_ecoli, hdr_mammalian, hdr_ecoli). Returns the minimum colony count, the actual probability achieved at that count, a 1.5x safety-margin recommendation, and branch-specific advice (low-efficiency screens get a "consider improving delivery" hint).
- Sample prompt: "Given my 45% editing efficiency, how many colonies should I pick to recover 2 homozygous clones with 99% confidence?"
- Output:
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Given my 45% editing efficiency, how many colonies should I pick to recover 2 homozygous clones with 99% confidence?
  [Tool call] → crispr_colony_calculator
  {
    "desired_clones": 2,
    "confidence": 0.99,
    "editing_efficiency": 0.45
  }
  [Tool result] ← crispr_colony_calculator:
  {
    "editing_efficiency": 0.45,
    "desired_clones": 2,
    "confidence": 0.99,
    "colonies_to_pick": 12,
    "expected_edits": 5.4,
    "probability_at_chosen_n": 0.9917,
    "safety_margin_recommendation": 18,
    "recommendation": "Pick 12 colonies to be 99% confident of recovering at least 2 edited clone(s) at 45% editing efficiency. Recommended saf\u2026 (307 chars)",
    "suggested_next_prompts": [
      "Design Sanger primers to verify the edit at the genomic locus.",
      "Once I have my Sanger results back, interpret an ICE result of 45% editing with R-squared 0.93.",
      "What if I want 3 edited clones instead of 2? Recompute."
    ],
    "citations": [
      {
        "label": "Kim et al. 2014, Genome Res 24:1012-1019",
        "reference": "https://doi.org/10.1101/gr.171322.113",
        "claim": "Cas9 RNP delivery to mammalian cells achieves ~50-80% editing efficiency, higher than plasmid transfection."
      },
      {
        "label": "Paquet et al. 2016, Nature 533:125-129",
        "reference": "https://doi.org/10.1038/nature17664",
        "claim": "HDR knockin efficiency is typically ~5-10% in mammalian cells, ~1/10 of NHEJ cutting efficiency."
      },
      {
        "label": "Jiang et al. 2013, Nat Biotechnol 31:233-239",
        "reference": "https://doi.org/10.1038/nbt.2508",
        "claim": "CRISPR-Cas9 editing in E. coli typically reaches 30-70% efficiency with pCRISPR/lambda Red recombineering systems."
      }
    ]
  }
  Gemini: To recover 2 edited clones with **99% confidence** given an editing efficiency of **45%**, you should pick **12 colonies**.
  To account for potential real-world losses—such as failed PCR, sequencing dropouts, or contamination—I recommend a safety margin, bringing your total to **18 colonies**.
  Since you have a relatively high editing efficiency, this is a very manageable number of colonies; a single 24-well plate will be more than sufficient for your screening.
  **Suggested Next Steps:**
  * [Design Sanger primers to verify the edit at the genomic locus.]
  * [Once I have my Sanger results back, interpret an ICE result of 45% editing with R-squared 0.93.]
  * [What if I want 3 edited clones instead of 2? Recompute.]
  ```
- Citations:
  - Kim, S., et al. 2014. "Highly efficient RNA-guided genome editing in human cells via delivery of purified Cas9 ribonucleoproteins." Genome Res 24:1012-1019. https://doi.org/10.1101/gr.171322.113
  - Paquet, D., et al. 2016. "Efficient introduction of specific homozygous and heterozygous mutations using CRISPR/Cas9." Nature 533:125-129. https://doi.org/10.1038/nature17664
  - Jiang, W., et al. 2013. "RNA-guided editing of bacterial genomes using CRISPR-Cas systems." Nat Biotechnol 31:233-239. https://doi.org/10.1038/nbt.2508
  - Laughery, M. F., et al. 2015. "New vectors for simple and streamlined CRISPR-Cas9 genome editing in Saccharomyces cerevisiae." Yeast 32(12):711-720. https://doi.org/10.1002/yea.3098

<b> 5. verify_edit </b>
- What it does: Designs Sanger sequencing primers flanking a Cas9 or Cas12a cut site at the genomic locus. Computes the cut position (3 bp upstream of NGG for Cas9), the expected amplicon size, and the cut offset within the amplicon. Validates each primer for Tm, GC content, and poly-N runs and warns if the pair has a |ΔTm| > 5°C. Returns a complete ICE/TIDE submission protocol with all coordinates pre-filled. The output ties directly into interpret_ice_tide once Sanger results come back.
- Sample prompt: "Design verification primers for protospacer TCAGAAACCTGCCAGTTTGC with the reference sequence CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTT."
- Output: 
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Design verification primers for protospacer TCAGAAACCTGCCAGTTTGC with the reference sequence CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTT...
  [Tool call] → crispr_verify_edit
  {
    "protospacer": "TCAGAAACCTGCCAGTTTGC",
    "reference": "CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTT"
  }
  [Tool result] ← crispr_verify_edit:
  {
    "protospacer": "TCAGAAACCTGCCAGTTTGC",
    "nuclease": "cas9",
    "strand": "+",
    "protospacer_position": 15,
    "pam_sequence": "TGG",
    "cut_position": 32,
    "forward_primer": "CCCTAGATGCCTGGCTCAGA",
    "forward_primer_position": 0,
    "forward_primer_tm": 64,
    "reverse_primer": "",
    "reverse_primer_position": 44,
    "reverse_primer_tm": 0,
    "tm_difference": 64,
    "primer_warnings": [
      "reverse primer Tm is 0.0degC (< 50degC): may not anneal reliably. Consider a longer primer or a more GC-rich region.",
      "reverse primer GC content 0% is below 40% \u2014 weak binding.",
      "Primer Tm mismatch: forward 64.0degC vs reverse 0.0degC (delta 64.0degC > 5degC). One primer will dominate at the anneal\u2026 (190 chars)"
    ],
    "amplicon_sequence": "CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGCTGGCACGTT",
    "amplicon_length": 44,
    "cut_offset_in_amplicon": 32,
    "interpretation_guide": "Edit Verification Protocol (ICE/TIDE) \u2014 CAS9:\n\n1. PCR: Amplify the target locus using the forward primer (pos 0) and rev\u2026 (1321 chars)",
    "suggested_next_prompts": [
      "Generate a lab sheet that includes this genomic verification step.",
      "Once I have my Sanger trace back, interpret an ICE result of 45% editing with R-squared 0.93.",
      "How many colonies should I pick to recover an edited clone at typical Cas9 efficiency?"
    ],
    "citations": [
      {
        "label": "Jinek et al. 2012, Science 337:816-821",
        "reference": "https://doi.org/10.1126/science.1225829",
        "claim": "SpCas9 cuts 3bp upstream of the NGG PAM producing a blunt double-strand break."
      },
      {
        "label": "Wallace et al. 1979, Nucleic Acids Res 6(11):3543-3557",
        "reference": "https://doi.org/10.1093/nar/6.11.3543",
        "claim": "Primer melting temperature for short oligos approximated by Tm = 2*(A+T) + 4*(G+C) degC."
      },
      {
        "label": "Dieffenbach et al. 1993, PCR Methods Appl 3(3):S30-37",
        "reference": "https://doi.org/10.1101/gr.3.3.S30",
        "claim": "PCR primer design guidelines: Tm 50-65 degC, GC 40-60%, no poly-N runs, primer pair |delta-Tm| < 5 degC."
      },
      {
        "label": "Hsiau et al. 2019, bioRxiv 251082",
        "reference": "https://doi.org/10.1101/251082",
        "claim": "ICE (Synthego) deconvolves Sanger trace mixed peaks into indel alleles and reports a KO score; R^2 below 0.80 indicates \u2026 (142 chars)"
      },
      {
        "label": "Brinkman et al. 2014, Nucleic Acids Res 42(22):e168",
        "reference": "https://doi.org/10.1093/nar/gku936",
        "claim": "TIDE quantifies CRISPR indel frequencies by linear decomposition of mixed Sanger trace signal downstream of the cut site\u2026 (121 chars)"
      }
    ]
  }
  Gemini: I have designed the verification primers for your protospacer. Please note that the reverse primer design currently has warnings due to the very short length of the reference sequence provided (44 bp), which is limiting primer design options.
  ### Edit Verification Protocol (ICE/TIDE) — CAS9
  1. **PCR**: Amplify the target locus using the forward primer (`CCCTAGATGCCTGGCTCAGA`) and the target region. **Note:** Due to the short reference sequence, the reverse primer could not be optimally designed; you should design a longer reverse primer flanking this 44bp region to ensure successful PCR.
  2. **Cut site**: SpCas9 cuts 32 bp from the amplicon start.
  3. **Sanger Sequencing**: Submit the PCR product for Sanger sequencing.
  4. **ICE/TIDE Analysis**: Upload the trace to ICE (Synthego) or TIDE, specifying the cut position at **32 bp** into your amplicon.
  5. **Expected Result**: Indel signal becomes noisy/mixed starting at position ~32 bp.
  **Important Warnings:**
  * The reverse primer failed to meet standard quality metrics (Tm 0.0°C). You will need to use a longer template (at least 200–300 bp of the surrounding genome) to allow for the design of a valid, high-quality reverse primer pair.
  **Suggested Next Steps:**
  * [Generate a lab sheet that includes this genomic verification step.]
  * [Once I have my Sanger trace back, interpret an ICE result of 45% editing with R-squared 0.93.]
  * [How many colonies should I pick to recover an edited clone at typical Cas9 efficiency?]
  ```
- Citation:
  - Jinek, M., et al. 2012. "A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity." Science 337:816-821. https://doi.org/10.1126/science.1225829
  - Wallace, R. B., et al. 1979. "Hybridization of synthetic oligodeoxyribonucleotides to phi chi 174 DNA: the effect of single base pair mismatch." Nucleic Acids Res 6(11):3543-3557. https://doi.org/10.1093/nar/6.11.3543
  - Dieffenbach, C. W., et al. 1993. "General concepts for PCR primer design." PCR Methods Appl 3(3):S30-37. https://doi.org/10.1101/gr.3.3.S30
  - Hsiau, T., et al. 2019. "Inference of CRISPR Edits from Sanger Trace Data." bioRxiv 251082. https://doi.org/10.1101/251082
  - Brinkman, E. K., et al. 2014. "Easy quantitative assessment of genome editing by sequence trace decomposition." Nucleic Acids Res 42(22):e168. https://doi.org/10.1093/nar/gku936

<b> 6. interpret_ice_tide </b>
- What it does: Classifies the output of an ICE (Synthego) or TIDE (Brinkman 2014) Sanger trace analysis. Returns an editing-efficiency band (FAILED < 10%, MARGINAL 10-30%, GOOD 30-70%, EXCELLENT ≥ 70%), an R²-based fit-quality classification (LOW < 0.80, MODERATE 0.80-0.89, HIGH ≥ 0.90), a reliability flag (False if R² < 0.80, then result should not be trusted regardless of editing percentage), warnings, and concrete next-step recommendations. If indel_distribution is supplied the dominant indel is highlighted and contradictions with editing_pct are flagged.
- Sample prompt: "My ICE score is 45% with R² = 0.93, dominant indel +1 at 38%."
- Output: 
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: My ICE score is 45% with R² = 0.93, dominant indel +1 at 38%.
  [Tool call] → crispr_interpret_ice_tide
  {
    "editing_pct": 45,
    "r_squared": 0.93,
    "indel_distribution": {
      "\"+1\"": 38
    },
    "tool": "ice"
  }
  [Tool result] ← crispr_interpret_ice_tide:
  {
    "sample_id": null,
    "tool": "ice",
    "editing_pct": 45.0,
    "r_squared": 0.93,
    "efficiency_classification": "GOOD",
    "fit_quality": "HIGH",
    "is_reliable": true,
    "dominant_indel": "\"+1\"",
    "dominant_indel_pct": 38.0,
    "warnings": [],
    "next_steps": [
      "Editing is good (30-70%). Pick colonies and screen by Sanger sequencing of single clones to recover homozygous edited li\u2026 (124 chars)"
    ],
    "suggested_next_prompts": [
      "How many colonies should I pick at 45% editing efficiency to recover one homozygous clone?",
      "Design Sanger primers to verify single clones at the genomic locus.",
      "I picked clones \u2014 interpret two more ICE results: A1A is 95% R^2=0.97, A1B is 12% R^2=0.85."
    ],
    "summary": "ICE reports 45.0% editing with R^2 = 0.93 (high fit quality). This indicates successful editing (good). Dominant indel: \u2026 (134 chars)",
    "citations": [
      {
        "label": "Hsiau et al. 2019, bioRxiv 251082",
        "reference": "https://doi.org/10.1101/251082",
        "claim": "ICE (Synthego) deconvolves Sanger trace mixed peaks into indel alleles and reports a KO score; R^2 below 0.80 indicates \u2026 (142 chars)"
      }
    ]
  }
  Gemini: Your ICE analysis results indicate a successful editing experiment:
  *   **Editing Efficiency (KO score):** 45.0%
  *   **Fit Quality (R²):** 0.93 (High)
  *   **Verdict:** **GOOD** — your experiment worked well.
  *   **Dominant Indel:** +1 insertion at 38%.
  Because the fit quality is high (R² ≥ 0.90), these results are reliable. 
  **Next Steps:**
  Since your efficiency is in the 30–70% range, you should now proceed to pick individual colonies and screen them by Sanger sequencing to isolate homozygous edited lines.
  **You could ask next:**
  * [How many colonies should I pick at 45% editing efficiency to recover one homozygous clone?]
  * [Design Sanger primers to verify single clones at the genomic locus.]
  * [I picked clones — interpret two more ICE results: A1A is 95% R^2=0.97, A1B is 12% R^2=0.85.]
  ```
- Citations:
  - Hsiau, T., et al. 2019. "Inference of CRISPR Edits from Sanger Trace Data." bioRxiv 251082. https://doi.org/10.1101/251082
  - Brinkman, E. K., et al. 2014. "Easy quantitative assessment of genome editing by sequence trace decomposition." Nucleic Acids Res 42(22):e168. https://doi.org/10.1093/nar/gku936 

<b> 7. Shared infrastructure I built that the team also uses:</b>
- citations.py
  - repo-wide citation registry. Every claim a CRISPR or labsheet tool returns links back to a paper or manufacturer manual through this single source of truth. Imported by predict_offtargets, predict_editing_efficiency, verify_edit, interpret_ice_tide, colony_calculator, and generate_lab_sheet.
- _protocols.py
  - auditable wet-lab protocol registry (Q5 PCR, Golden Gate, Gibson, Type IIS oligo cloning, RestrictionLigation, chemical transformation, Zymo miniprep, Sanger sequencing, RNP assembly, electroporation, lipofection, IVT sgRNA). Reagent volumes are transcribed from manufacturer manuals or original method papers, never invented. Each entry carries a source citation key so generate_lab_sheet can ship a fully auditable protocol_sources field.
</div>
