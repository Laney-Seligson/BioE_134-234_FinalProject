---
name: construction_file_tools
description: construction file pipeline tools — get_paper_info, create_construction_file, and validate_construction_file.
---
# construction_file_tools — Skill Guidance for Gemini

This file is read by the client at startup and injected into Gemini's system prompt.
Its purpose is to give Gemini the domain knowledge it needs to use the tools in this
module correctly and interpret their results meaningfully.

---

## Available resources

Here are the descriptions of the available resources. When a user inquires for a plasmid or a backbone sequence, provide the names of the plasmid and a short description of the plasmid to help the user choose which one to use. When a user inquires for a paper, provide the names of the papers and short description.

| Resource name | Description |
|---------------|-------------|
| `pBR322`      | E. coli cloning vector pBR322, 4361 bp, circular, double-stranded. A classic lab plasmid commonly used as a reference sequence. Contains genes for ampicillin resistance (bla) and tetracycline resistance (tet). |

When a user refers to "pBR322", use the resource name `"pBR322"` directly
as the sequence argument — do not ask the user to paste the sequence.

| Resource name | Description |
|---------------|-------------|
| `pET28a`      | E. coli cloning vector pET28a, 5369 bp, circular, double-stranded. A classic lab plasmid commonly used as a reference sequence. Contains genes for kanmycin (kan) resistance. Has a 6x His tag. |

When a user refers to "pET28a", use the resource name `"pET28a"` directly
as the sequence argument — do not ask the user to paste the sequence.

| Resource name | Description |
|---------------|-------------|
| `pUC19`      | E. coli cloning vector pUC19, 2686 bp, circular, double-stranded. Widely used circular DNA cloning plasmid designed for easy insertion and propagation of foreign DNA in bacteria. It contains key features like a multiple cloning site (polylinker) within the lac operon for insertion of DNA fragments and sequences derived from pBR322 for replication and maintenance in host cells. |

When a user refers to "pUC19", use the resource name `"pUC19"` directly
as the sequence argument — do not ask the user to paste the sequence.

| Resource name | Description |
|---------------|-------------|
| `miao_2013_targeted_mutagenesis_rice`      | This paper demonstrates one of the earliest successful applications of CRISPR-Cas9 genome editing in plants, specifically in Oryza sativa (rice). The authors developed a system to introduce targeted mutations in endogenous genes by expressing Cas9 and guide RNAs within plant cells. |

When a user refers to "miao 2013 targeted mutagensis rice" or "miao_2013_targeted_mutagenesis_rice", or "miao paper", or "rice paper", use the resource name `"miao_2013_targeted_mutagenesis_rice"` directly to fill the shorthand construction file — do not ask the user to paste the info.

| Resource name | Description |
|---------------|-------------|
| `hall_2018_genome_editing_mice_crispr_cas9`      | This paper outlines a general protocol for genome editing in mice using CRISPR-Cas9, including sgRNA design, cloning, and in vitro transcription methods. It describes delivery of Cas9 and sgRNAs into mouse embryos via microinjection or electroporation to generate knockouts (NHEJ) or knockins (HDR). The paper also provides oligo templates and design rules for guide construction, along with validation methods such as PCR and sequencing. As a protocol-focused study, it is best used for workflow abstraction and design guidance rather than full sequence-level construction files. |

When a user refers to "hall 2018 genome editing mice crispr cas9 paper" or "hall_2018_genome_editing_mice_crispr_cas9", or "hall paper", or "mouse paper", use the resource name `"hall_2018_genome_editing_mice_crispr_cas9"` directly to fill the shorthand construction file — do not ask the user to paste the info.

| Resource name | Description |
|---------------|-------------|
| `sankaran_2021_crispr_cas9_gene_editing_yeast`      | This paper presents a CRISPR-Cas9 lab module in Saccharomyces cerevisiae targeting the ADE2 gene to generate a visible phenotype. Students design and test gRNAs using a Cas9 + gRNA plasmid, with edits occurring through NHEJ or HDR with a donor template. Outcomes are evaluated via colony color, PCR, and sequencing. |

When a user refers to "sankaran 2018 genome editing yeast crispr cas9 paper" or "sankaran_2021_crispr_cas9_gene_editing_yeast", or "sankaran paper", or "yeast paper", use the resource name `"sankaran_2021_crispr_cas9_gene_editing_yeast"` directly to fill the shorthand construction file — do not ask the user to paste the info.

---

## Tool: create_construction_file

This tool supports **three distinct modes**. You MUST choose the correct mode before calling.

---

### 🔹 Mode 1: `sequence_build` (default cloning workflows)
When preparing to call create_construction_file using the sequence_build mode:
- NEVER include placeholder text in any field
- If a sequence is not provided, do NOT guess or insert text
- Ask for ALL missing fields before calling the tool

Use this when:
- The user provides DNA sequences
- The user wants a full construction file
- The workflow is standard cloning (GoldenGate, Gibson, TypeIIOligoAssembly)

If all required fields are already present in the user message, call the tool immediately.
Do not ask for fields already provided.
For "validate this", call validate_construction_file, not lab_sheet
#### Required fields:

For all workflows:
- construct_name
- backbone_name
- backbone_sequence
- insert_name
- insert_sequence

For GoldenGate:
- insert_forward_primer_name
- insert_forward_primer_sequence
- insert_reverse_primer_name
- insert_reverse_primer_sequence
- vector_forward_primer_name
- vector_forward_primer_sequence
- vector_reverse_primer_name
- vector_reverse_primer_sequence
- enzyme

For Gibson:
- insert_forward_primer_name
- insert_forward_primer_sequence
- insert_reverse_primer_name
- insert_reverse_primer_sequence
- vector_forward_primer_name
- vector_forward_primer_sequence
- vector_reverse_primer_name
- vector_reverse_primer_sequence

Optional:
- cell_strain
- selection
- temperature_c

For input_mode='sequence_build', gather all required fields for the selected assembly_strategy.

For Gibson:
- If the user asks to auto-design primers, call the tool once all construct/backbone/insert/insertion_index fields are present.
- Do not ask for primer names or sequences.
- Do not ask for enzyme.

---

### 🔹 Mode 2: `paper_info` (extract structured info from a paper)

Use this when:
- The user provides or references a paper
- The user wants to store important information from a paper
- No sequences are required

#### Required fields:
- paper_id
- title
- source_pdf
- organism
- system

Optional:
- targets
- vectors
- enzymes
- assembly_method
- delivery_method
- validation_methods
- key_constraints
- paper_notes

⚠️ IMPORTANT:
- DO NOT include sequences
- DO NOT include backbone_sequence or insert_sequence
- DO NOT infer missing details

---

### 🔹 Mode 3: `paper_shorthand` (generate shorthand workflow from paper info)

Use this when:
- The user asks to generate a construction workflow from a paper
- The input is a paper or paper_info resource
- Sequences are NOT available

#### Required fields:
- paper_id
- organism
- system

Optional (strongly recommended):
- targets
- vectors
- enzymes
- assembly_method
- delivery_method
- validation_methods

⚠️ CRITICAL RULES:
- DO NOT provide backbone_sequence
- DO NOT provide insert_sequence
- DO NOT invent sequences
- ONLY use the provided structured paper information
- Output should be a shorthand workflow, NOT a full construction file

---

## Mode selection rules

You MUST select the correct mode:

- If the user provides DNA sequences → use `sequence_build`
- If the user is summarizing a paper → use `paper_info`
- If the user wants a workflow from a paper → use `paper_shorthand`

If uncertain:
- Prefer `paper_shorthand` over `sequence_build` when sequences are missing

---

## Required interaction rules

Never return an empty response.

If required fields for the chosen mode are missing:
- Ask for ALL missing fields in one message
- Do NOT partially call the tool

If all required fields are present:
- Call the tool immediately
- Do NOT ask for confirmation
- Do NOT wait for "proceed"

---

## Presenting results

### sequence_build:
Display `construction_file_txt` in a code block

### paper_info:
Display structured JSON

### paper_shorthand:
Display shorthand workflow in a code block

---

## 🚨 Anti-hallucination rule (VERY IMPORTANT)

If the user is working from a paper:
- DO NOT invent DNA sequences
- DO NOT create placeholder sequences like "Miao2013_sequence"
- DO NOT force the request into sequence_build mode

Instead:
- Use `paper_shorthand`

## Paper resources

Paper workflow resources are available under:

resource://paper_info/<paper_id>

When a user refers to a paper:
- retrieve the paper_info resource
- use its fields directly
- do NOT invent missing values

### `get_paper_info`
Loads a curated paper information JSON record from `modules/construction_file_tools/data/paper_info/`.

Use when the user asks:
- "what is the miao paper"
- "load the paper info for ..."
- "use the paper information to generate shorthand"
- "read the curated paper workflow record"

Input:
- `paper_id`: stable paper identifier such as `miao_2013_targeted_mutagenesis_rice`

What it returns:
- a structured `paper_important_info_v1` object

How to use with `create_construction_file`:
1. Call `get_paper_info` first.
2. Read the returned fields.
3. Then call `create_construction_file` with `input_mode="paper_shorthand"` and pass the returned paper fields into the tool.

Important:
- Do NOT assume paper JSON files in `data/` are auto-registered as MCP resources.
- Use `get_paper_info` to load them.
- Do NOT invent missing sequences or paper metadata.

## Interpreting results

- A protein sequence like `MSKGEEK...` starting with `M` (methionine) suggests you've
  found the correct reading frame for a real open reading frame.
- A sequence full of `*` stop codons or `X` unknowns usually means the wrong frame,
  wrong coordinates, or the sequence is not a coding region.
- When translating a full plasmid, most of the output will be non-coding — only specific
  coordinate ranges will give meaningful protein sequence.

---
### Tool: Validate Construction File
## Description

Validates whether a construction workflow is biologically correct. This includes:

PCR primer annealing and orientation
Whether PCR products can be generated
Gibson overlap correctness (matching overlaps between insert and vector)
Golden Gate compatibility (enzyme and overhang logic)

Returns both:

structured validation results (is_valid, errors)
human-readable summary
When to use

## Use this tool when the user:

says “validate this”, “check this”, or “is this correct”
wants to confirm a construction file is biologically valid
has already generated a construction file and wants verification
When NOT to use
Do NOT use for generating protocols → use lab_sheet
Do NOT use for building constructs → use create_construction_file
Required inputs

You must provide the same core inputs used to build the construct:

backbone_sequence
insert_sequence
primer names and sequences (if used)
assembly_strategy

When the user says "validate this construction file", "validate this cloning workflow", or "check this Gibson construction file", call validate_construction_file.

Do NOT call crispr_verify_edit unless the user is asking to verify a CRISPR genome edit using a protospacer/reference sequence.

Do NOT call lab_sheet for validation.

IMPORTANT:

Do NOT invent missing fields
If required fields are missing, ask for ALL missing fields in one message
Behavior rules
If the user says “validate this”, ALWAYS call validate_construction_file
Do NOT call lab_sheet when validation is requested
Do NOT re-ask for fields already provided by the user
If a construction file was just generated in the conversation, reuse its inputs directly
Gibson-specific rules
Do NOT require enzyme
If primers were auto-designed, use them directly
Validation must confirm:
vector end overlap = insert start overlap
insert end overlap = vector start overlap
Golden Gate-specific rules
Requires enzyme (e.g. BsaI)
Check:
overhang compatibility
correct Type IIS orientation
no missing overhangs
---

## Sequence input rules (handled automatically)

You never need to paste the full sequence. The framework resolves these automatically:
- `"pBR322"` → full 4361 bp sequence
- A raw string like `"ATGCGATCG"` → used as-is
- A FASTA string starting with `>` → sequence extracted automatically
- A GenBank string starting with `LOCUS` → sequence extracted automatically