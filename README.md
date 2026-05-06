# BioE134/234 CRISPR Pipeline Project
Names: Emory Elizabeth Adelman, Jillian Ho, Karina Jin, Laney Seligson

Presentation slide: [click this](https://docs.google.com/presentation/d/1eYHWdArwtcOk_REBWtsuq_iwyRDxQa1vGYRxENukCQE/edit?usp=sharing)


---
## Table of Contents
1. [Overview](#1-overview)
2. [Project Structure](#2-project-structure)
3. [How the Pipeline Works](#3-how-the-pipeline-works)
4. [Setup](#4-setup)
5. [Individual Scope](#5-individual-scope)
    - [Karina](#karina)
    - [Emory](#emory)
    - [Laney](#laney)
    - [Jillian](#jillian)
6. [Demo](#6-demo-for-whole-crispr-pipeline)
7. [Citation](#7-citation)

## 1. Overview

Designing CRISPR experiments typically requires coordinating multiple tools, manual sequence handling, and careful validation, making the process time-consuming and error-prone. This project is a modular pipeline designed to streamline and automate the full CRISPR and cloning workflow, taking a user from <b>an initial idea to a lab-ready experimental plan</b>. Instead of handling each step manually, the system organizes the process into connected components that guide users through sequence selection, design, and validation in a structured way. It includes built-in gene and sequence lookup functionality, allowing users to start from a gene name, organism, or identifier and automatically retrieve the relevant DNA sequence for downstream design. The pipeline supports upstream and downstream design. Including but not limited to, identifying target sequences, designing and evaluating guide RNAs, oligo design, construction file creation, construction file validation, and lab sheet generation. Overall, the goal of this project is to create a more efficient, reproducible, and reliable approach to experimental design in synthetic biology by combining structured workflows with intelligent automation.

---

## 2. Project structure
![Pipeline Flowchart](readme_appendix/pipeline_flowchart.png)

---

## 3. How the pipeline works

```
Tentative leave for more use
```

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

### Emory: 
<div style="margin-left: 20px;">
<b>1. run_full_crispr_workflow </b>

- What it does: 
  - Fetch a target sequence from the NCBI GenBank file for the gene chosen 
  - design guide RNA using the first 10 found in the fetched sequence
  - rank guide RNAs with rank_guides, authored by Jillian, and find the single best one to use out of the 10 
  - design cloning oligos 
  - build and validate construction file scripts authored by Laney
 
- MCP Wrapper: `run_full_crispr_workflow.json`
- Pytests (`tests/unit/test_run_full_crispr_workflow.py`): guide selection, vector prompting, workflow confirmation gate, upstream-selected gene confirmation, empty query error, per-guide score fields
- Sample Prompt: “Design a CRISPR edit targeting lacZ in E.Coli using pTargetF.”

<b>2. fetch_target_sequence</b>
- What it does: 
  - Resolves a gene name or accession to a clean DNA sequence using NCBI Entrez tools 

- MCP Wrapper: `fetch_target_sequence.json`
- Pytests (`tests/unit/test_fetch_target_sequence.py`): raw DNA passthrough, uppercase normalization, invalid DNA error, empty query error, required output keys
- Sample Prompt: “Fetch the sequence for recA in E.Coli.” 
- Output: 
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  You: Fetch the sequence for recA in E.Coli
  [Tool call] → crispr_fetch_target_sequence
  {
    "organism": "Escherichia coli",
    "query": "recA"
  }
  [Tool result] ← crispr_fetch_target_sequence:
  {
    "sequence": "ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACC\u2026 (1062 chars)",
    "source": "ncbi",
    "resource": "recA",
    "organism": "Escherichia coli",
    "ncbi_gene_id": "947170",
    "ncbi_accession": "NC_000913.3",
    "length": 1062,
    "note": "Fetched from NCBI Entrez: gene 'recA' in 'Escherichia coli'. Nucleotide record: NC_000913.3. Verify the record matches y\u2026 (151 chars)"
  }
  Gemini: The DNA sequence for the `recA` gene in *Escherichia coli* (NCBI Gene ID: 947170, NC_000913.3) has been retrieved:
  `ATGGCTATCGACGAAAACAAACAGAAAGCGTTGGCGGCAGCACTGGGCCAGATTGAGAAACAATTTGGTAAAGGCTCCATCATGCGCCTGGGTGAAGACCGTTCCATGGATGTGGAAACCATCTCTACCGGTTCGCTTTCACTGGATATCGCGCTTGGGGCAGGTGGTCTGCCGATGGGCCGTATCGTCGAAATCTACGGACCGGAATCTTCCGGTAAAACCACGCTGACGCTGCAGGTGATCGCCGCAGCGCAGCGTGAAGGTAAAACCTGTGCGTTTATCGATGCTGAACACGCGCTGGACCCAATCTACGCACGTAAACTGGGCGTCGATATCGACAACCTGCTGTGCTCCCAGCCGGACACCGGCGAGCAGGCACTGGAAATCTGTGACGCCCTGGCGCGTTCTGGCGCAGTAGACGTTATCGTCGTTGACTCCGTGGCGGCACTGACGCCGAAAGCGGAAATCGAAGGCGAAATCGGCGACTCTCACATGGGCCTTGCGGCACGTATGATGAGCCAGGCGATGCGTAAGCTGGCGGGTAACCTGAAGCAGTCCAACACGCTGCTGATCTTCATCAACCAGATCCGTATGAAAATTGGTGTGATGTTCGGTAACCCGGAAACCACTACCGGTGGTAACGCGCTGAAATTCTACGCCTCTGTTCGTCTCGACATCCGTCGTATCGGCGCGGTGAAAGAGGGCGAAAACGTGGTGGGTAGCGAAACCCGCGTGAAAGTGGTGAAGAACAAAATCGCTGCGCCGTTTAAACAGGCTGAATTCCAGATCCTCTACGGCGAAGGTATCAACTTCTACGGCGAACTGGTTGACCTGGGCGTAAAAGAGAAGCTGATCGAGAAAGCAGGCGCGTGGTACAGCTACAAAGGTGAGAAGATCGGTCAGGGTAAAGCGAATGCGACTGCCTGGCTGAAAGATAACCCGGAAACCGCGAAAGAGATCGAGAAGAAAGTACGTGAGTTGCTGCTGAGCAACCCGAACTCAACGCCGGATTTCTCTGTAGATGATAGCGAAGGCGTAGCAGAAACTAACGAAGATTTTTAA`
  ```
  </details>

- Citations: 
  - National Center for Biotechnology Information. n.d. "Entrez Programming Utilities Help (E-utilities)." U.S National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK25501/.

  - National Center for Biotechnology Information. n.d. "Gene Database Help." U.S. National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK3841/.

<b> 3. design_cas9_grna </b>
- What it does: 
  - Scans a sequence for NGG PAM sites and returns up to 10 SpCas9 guide RNA candidates 

- MCP Wrapper: `design_cas9_grna.json`
- Pytests (`tests/unit/test_design_cas9_grna.py`): 20 bp protospacer, NGG PAM, max 10 guides, no-PAM error, empty input error
- Sample Prompt: “Design Cas9 guides for the lacZ locus.” 
- Output: 
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

<b> 4. design_cas12a_crrna </b>
- What it does: 
  - Scans a sequence for TTTV PAM sites and returns up to 10 LbCas12a crRNA candidates 

- MCP Wrapper: `design_cas12a_crrna.json`
- Pytests (`tests/unit/test_design_cas12a_crrna.py`): 23 bp protospacer, TTTV PAM, max 10 guides, no-PAM error, empty input error
- Sample Prompt: “Give me Cas12a guides for the araB locus.”
- Output: 
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
  If you would like to proceed, I can help you select the best candidates based on specificity, design cloning oligos for your vector, or provide an editing efficiency prediction for a specific guide. How would you like to continue?
  ```
  </details>
- Citations:
  - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ian M. Slaymaker, Kira S. Makarova, Patrick Essletzbichler, Sara E. Volz, et al. 2015. "Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System." Cell 163 (3): 759–771. https://doi.org/10.1016/j.cell.2015.09.038.
  - Fonfara, Ines, Hagen Richter, Majda Bratovič, Anaïs Le Rhun, and Emmanuelle Charpentier. 2016. "The CRISPR-Associated DNA-Cleaving Enzyme Cpf1 Also Processes Precursor CRISPR RNA." Nature 532 (7600): 517–521. https://doi.org/10.1038/nature17945.

<b> 5. cas_selector</b>
- What it does: 
  - Recommends Cas9 or Cas12a based on guide availability, multiplexing needs, and specificity requirements 
- MCP Wrapper: `cas_selector.json`
- Pytests (`tests/unit/test_cas_selector.py`): GC-rich → Cas9, AT-rich → Cas12a, multiplexing override, empty sequence error
- Sample Prompt: “Which Cas nuclease should I use for 3 targets in E.Coli?” 
- Citations: 
  - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ian M. Slaymaker, Kira S. Makarova, Patrick Essletzbichler, Sara E. Volz, et al. 2015. "Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System." Cell 163 (3): 759–771. https://doi.org/10.1016/j.cell.2015.09.038.

  - Kim, Daesik, Byungkuk Min, Jungeun Kim, Seokjoong Kim, and Jin-Soo Kim. 2016. "Genome-Wide Analysis Reveals Specificities of Cpf1 Endonucleases in Human Cells." Nature Biotechnology 34 (8): 863–868. https://doi.org/10.1038/nbt.3609.

  - Zetsche, Bernd, Matthias Heidenreich, Prarthana Mohanraju, Ines Fedorova, Jeroen Kneppers, Ellen M. DeGennaro, Naomi Winblad, et al. 2017. "Multiplex Gene Editing by CRISPR-Cpf1 Using a Single crRNA Array." Nature Biotechnology 35 (1): 31–34. https://doi.org/10.1038/nbt.3737.


<b> 6. design_cloning_oligos </b>
- What it does: 
  - Designs annealed oligos or PCR primers into a preset vector 
    - TypeIIS (BbsI / BsmBI / BsaI annealed-oligo sticky-end ligation)
    - RestrictionLigation
    - Gibson
  - Golden Gate cloning 
  - **Addgene API fallback:** preset vectors are handled locally; if `vector` is a numeric Addgene plasmid ID (e.g. `"42230"` or `"addgene:67639"`), the tool attempts to fetch plasmid metadata from the Addgene Developers API. Requires `ADDGENE_API_KEY` in `.env`. Dynamic results include citations in the same `{ "label", "reference", "claim" }` format as other tools. If no API key is available, use a preset vector or `vector="custom"`. This feature does not have and will not have its own JSON wrapper — it is handled entirely within `design_cloning_oligos`.
  - **pX330 backbone:** `px330` now uses a local 8484 bp GenBank backbone resource instead of the previous `"N"` placeholder. Addgene #42230 is cited for the plasmid record and NovoPro V005940 is cited as the public GenBank sequence source.
  - if organism/vector info is missing → returns `needs_user_input` with human-readable questions instead of an error

- MCP Wrapper: `design_cloning_oligos.json`
- Pytests (`tests/unit/test_design_cloning_oligos.py`): pCRISPR E. coli produces ready oligos; pml104 enzyme is BclI-SwaI, top overhang is GATC, bottom overhang is blank, top/bottom oligo sequences verified; pml107 present with LEU2 selection; organism mismatch → needs_user_input
- Sample Prompt: “Design oligos to clone this protospacer into px330.” 
- Output: 
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
  - Ran, F. Ann, Patrick D. Hsu, Jason Wright, Vineeta Agarwala, David A. Scott, and Feng Zhang. 2013. "Genome Engineering Using the CRISPR-Cas9 System." Nature Protocols 8 (11): 2281–2308. https://doi.org/10.1038/nprot.2013.143.

  - Sanjana, Neville E., Ophir Shalem, and Feng Zhang. 2014. "Improved Vectors and Genome-Wide Libraries for CRISPR Screening." Nature Methods 11 (8): 783–784. https://doi.org/10.1038/nmeth.3047.

  - Jiang, Wenyan, David Bikard, David Cox, Feng Zhang, and Luciano A. Marraffini. 2013. "RNA-Guided Editing of Bacterial Genomes Using CRISPR-Cas Systems." Nature Biotechnology 31 (3): 233–239. https://doi.org/10.1038/nbt.2508.

  - Gibson, Daniel G., Lei Young, Ray-Yuan Chuang, J. Craig Venter, Clyde A. Hutchison III, and Hamilton O. Smith. 2009. "Enzymatic Assembly of DNA Molecules up to Several Hundred Kilobases." Nature Methods 6 (5): 343–345. https://doi.org/10.1038/nmeth.1318.

  - Engler, Carola, Ramona Kandzia, and Sylvestre Marillonnet. 2008. "A One Pot, One Step, Precision Cloning Method with High Throughput Capability." PLOS ONE 3 (11): e3647. https://doi.org/10.1371/journal.pone.0003647.

  - Laughery, Marian F., Taylor Hunter, Angela Brown, Jake Hoopes, Tiffany Ostbye, Trevor Shumaker, and John J. Wyrick. 2015. "New Vectors for Simple and Streamlined CRISPR-Cas9 Editing in *Saccharomyces cerevisiae*." Yeast 32 (12): 711–720. https://doi.org/10.1002/yea.3098.

  - Addgene. "pX330-U6-Chimeric_BB-CBh-hSpCas9 (Plasmid #42230)." https://www.addgene.org/42230/

  - NovoPro Bioscience. "pX330-U6-Chimeric_BB-CBh-hSpCas9 vector (Cat. No. V005940)." https://www.novoprolabs.com/vector/Vgy3dima
</div>

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
</div>

### Jillian: 
<div style="margin-left: 20px;">
<b>1. predict_offtargets </b>
- What it does:
  - Scans a reference DNA sequence on both strands for possible CRISPR off-target sites matching a guide within a chosen mismatch threshold
  - Supports Cas9 NGG PAMs and Cas12a TTTV PAMs
  - Scores each candidate site as HIGH, MEDIUM, or LOW risk using PAM presence, total mismatches, seed-region mismatches, and a simplified CFD-style score
  - Supports circular references so plasmid sites that span the sequence origin are not missed

- MCP Wrapper: `predict_offtargets.json`
- Sample Prompt: "Check if the guide TCAGAAACCTGCCAGTTTGC has any off-target sites in pBR322."
- Output:
  ```json
  [Tool result] ← crispr_predict_offtargets:
  {
    "protospacer": "TCAGAAACCTGCCAGTTTGC",
    "nuclease": "cas9",
    "reference_length": 4361,
    "strands_scanned": ["+", "-"],
    "sites_evaluated": 8684,
    "offtarget_sites": [],
    "high_risk_count": 0,
    "aggregate_offtarget_cfd": 0,
    "max_offtarget_cfd": 0,
    "specificity_summary": "No off-target sites detected within 3 mismatches for cas9."
  }
  ```

- Citations:
  - Hsu, Patrick D., David A. Scott, Joshua A. Weinstein, F. Ann Ran, Silvana Konermann, Vineeta Agarwala, Yinqing Li, et al. 2013. "DNA Targeting Specificity of RNA-Guided Cas9 Nucleases." Nature Biotechnology 31 (9): 827-832. https://doi.org/10.1038/nbt.2647.

  - Doench, John G., Nicolo Fusi, Meagan Sullender, Mudra Hegde, Emma W. B. Hegde, et al. 2016. "Optimized sgRNA Design to Maximize Activity and Minimize Off-Target Effects of CRISPR-Cas9." Nature Biotechnology 34 (2): 184-191. https://doi.org/10.1038/nbt.3437.

<b> 2. rank_guides </b>
- What it does: 
  - Scores guide candidates 
    - GC content
    - poly-T runs
    - off-target risk
  - returns a sorted list with rationale — scoring denominator is now explicit (X/3, efficiency X/2, specificity X/1) to prevent misreading the scale
- MCP Wrapper: `rank_guides.json`
- Pytests (`tests/unit/test_rank_guides.py`): sorted output, poly-T penalty, GC scoring, citations present, empty list error
- Sample Prompt: “Rank these guides against the lacZ sequence.”
- Output: 
  <details>
  <summary>Click to expand sample output</summary>

  ```json
  [Tool result] ← crispr_rank_guides:
  {
    "ranked_guides": [
      {
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUGCAGUGAUUCUGUGCGAGCUUUG",
        "protospacer": "GCAGTGATTCTGTGCGAGCTTTG",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.522,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUGCGGUGGACUGCGCUACCGGUGA",
        "protospacer": "GCGGTGGACTGCGCTACCGGTGA",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.696,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUUGAUGCCCCGAAUAACCAGUUCC",
        "protospacer": "TGATGCCCCGAATAACCAGTTCC",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.522,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUCCGAAAACCCGAACGCGAUGUUC",
        "protospacer": "CCGAAAACCCGAACGCGATGTTC",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.565,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "protospacer": "TGCCACGCGCCGGGCAACGTTGA",
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUUGCCACGCGCCGGGCAACGUUGA",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.696,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "pam_site": "TTTA",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUUUCCAGCGAAUGGUUCUGGGCAA",
        "protospacer": "TTCCAGCGAATGGTTCTGGGCAA",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.522,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "protospacer": "CGGTACCACCCGCCCGCAGGATA",
        "pam_site": "TTTC",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUCGGUACCACCCGCCCGCAGGAUA",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.696,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "protospacer": "TTTGATGAGCTGGACCCGATCCT",
        "pam_site": "TTTC",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUUUUGAUGAGCUGGACCCGAUCCU",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.522,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUAUGAGCUGGACCCGAUCCUCAAU",
        "protospacer": "ATGAGCTGGACCCGATCCTCAAT",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.522,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      },
      {
        "protospacer": "CCTTCCCCGCTGTTCACTGACAC",
        "pam_site": "TTTG",
        "grna_sequence": "AAUUUCUACUAAGUGUAGAUCCUUCCCCGCUGUUCACUGACAC",
        "efficiency_score": 2,
        "efficiency_details": {
          "gc_content": 0.609,
          "gc_content_ok": true,
          "no_polyt_run": true
        },
        "specificity_score": 1,
        "specificity_details": {
          "high_risk_offtargets": 0,
          "medium_risk_offtargets": 0,
          "total_offtargets": 0,
          "on_target_found": false
        },
        "total_score": 3
      }
    ],
    "best_guide": {
      "pam_site": "TTTG",
      "grna_sequence": "AAUUUCUACUAAGUGUAGAUGCAGUGAUUCUGUGCGAGCUUUG",
      "protospacer": "GCAGTGATTCTGTGCGAGCTTTG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.522,
        "gc_content_ok": true,
        "no_polyt_run": true
      },
      "specificity_score": 1,
      "specificity_details": {
        "high_risk_offtargets": 0,
        "medium_risk_offtargets": 0,
        "total_offtargets": 0,
        "on_target_found": false
      },
      "total_score": 3
    },
    "scoring_rationale": "Best guide 'GCAGTGATTCTGTGCGAGCTTTG' scored 3/3 (efficiency 2/2, specificity 1/1). GC content: 52%. Off-target sites: 0 total \u2026 (139 chars)",
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
      },
      {
        "label": "Zetsche et al. 2015, Cell 163(3):759-771",
        "reference": "https://doi.org/10.1016/j.cell.2015.09.038",
        "claim": "LbCas12a recognizes a TTTV (V=A/C/G) PAM 5' of a 23bp protospacer and produces a staggered cut with a 5-nt 5' overhang."
      },
      {
        "label": "Kim et al. 2016, Nat Biotechnol 34:863-868",
        "reference": "https://doi.org/10.1038/nbt.3609",
        "claim": "Cas12a seed region is positions 1-10 from the PAM-proximal end; Cas12a is more mismatch-intolerant than Cas9 across the \u2026 (131 chars)"
      },
      {
        "label": "Kim et al. 2018, Nat Biotechnol 36:239-241",
        "reference": "https://doi.org/10.1038/nbt.4061",
        "claim": "Position-specific nucleotide preferences for Cas12a on-target activity (DeepCpf1 model)."
      }
    ]
  }
  Gemini: I have ranked the Cas12a crRNA candidates for the *araB* gene based on their on-target efficiency and off-target specificity against the *lacZ* reference sequence.
  The top-ranked guide is:
  **Protospacer:** `GCAGTGATTCTGTGCGAGCTTTG`
  **PAM:** `TTTG`
  **Rationale:** This guide has a GC content of 52%, contains no poly-T runs, and has zero predicted off-target sites in the *lacZ* reference sequence, making it the most specific and likely efficient candidate.
  Would you like to design cloning oligos for a specific vector using this guide, or would you like to proceed with a different step?
  ```
</details>

- Citations: 
  - Doench, John G., Nicolo Fusi, Meagan Sullender, Mudra Hegde, Emma W. Vaimberg, Katherine F. Donovan, Ian Smith, et al. 2016. "Optimized sgRNA Design to Maximize Activity and Minimize Off-Target Effects of CRISPR-Cas9." Nature Biotechnology 34 (2): 184–191. https://doi.org/10.1038/nbt.3437.

  - Hsu, Patrick D., David A. Scott, Joshua A. Weinstein, F. Ann Ran, Silvana Konermann, Vineeta Agarwala, Yinqing Li, et al. 2013. "DNA Targeting Specificity of RNA-Guided Cas9 Nucleases." Nature Biotechnology 31 (9): 827–832. https://doi.org/10.1038/nbt.2647.

  - Kim, Heon Seok, Yongmuk Kim, Kunwoo Lee, Sung Il Cho, and Jin-Soo Kim. 2018. "Deep Learning Improves Prediction of CRISPR-Cpf1 Guide RNA Activity." Nature Biotechnology 36 (3): 239–241. https://doi.org/10.1038/nbt.4061.

  - Bogenhagen, Daniel F., and Donald D. Brown. 1981. "Nucleotide Sequences in Xenopus 5S DNA Required for Transcription Termination." Cell 24 (1): 261–270. https://doi.org/10.1016/0092-8674(81)90522-5.
</div>

## 6. Demo for whole CRISPR pipeline
[Insert video link]
