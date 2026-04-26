---
name: crispr_tools
description: CRISPR pipeline tools — guide RNA design, off-target prediction, cloning oligo design, edit verification, and ICE/TIDE interpretation.
---

# crispr_tools — Skill Guidance for Gemini

This file is read by the client at startup and injected into Gemini's system prompt.
Its purpose is to give Gemini the domain knowledge it needs to use the tools in this
module correctly and interpret their results meaningfully.

---

## What this module does

The `crispr_tools` module provides fundamental tools to go through the crispr pipeline.

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


---

## Tools and when to use them

### `dna_reverse_complement`
Returns the reverse complement of a DNA or RNA sequence.

Use when the user asks for:
- "reverse complement of X"
- "complement of the bottom strand"
- "what does the antisense strand look like"
- "flip the sequence"

The result is the same length as the input. Uppercase output.

### `dna_translate`
Translates a DNA coding sequence to a protein sequence using the standard genetic code.

Use when the user asks to:
- "translate", "get the protein", "what protein does this encode"
- work with a specific reading frame (1, 2, or 3)
- translate a specific region using `start` / `end` coordinates (0-indexed, end is exclusive)

**Frame guidance:**
- Frame 1 — start reading from the first base (default)
- Frame 2 — skip 1 base, then read triplets
- Frame 3 — skip 2 bases, then read triplets

**Stop codons** appear as `*` in the output. **Unrecognised codons** appear as `X`.

**Coordinate example:** "translate bases 100 to 200" → `start=100, end=200`
"translate the first 60bp" → `start=0, end=60` (or omit start, set `end=60`)

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
- The workflow is standard cloning (GoldenGate, Gibson, DirectSynthesis)

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

### Required behavior after create_construction_file

After successfully calling `create_construction_file` and presenting the result, **always** follow up with:

> "Would you like me to also generate a bench-ready lab sheet for this construct? It will give you step-by-step instructions for PCR, assembly, transformation, and sequencing."

If the user says yes (or anything like "sure", "yes please", "go ahead"), immediately call `crispr_lab_sheet` with the `structured_construction_file` returned by `create_construction_file`.

**Important:** Only offer the construction file and lab sheet after the full CRISPR validation workflow is complete (guide design, off-target check, and cloning oligo design). Do not offer either output mid-workflow. If the full workflow offer (see below) already asked about both construction file and lab sheet together and the user chose "construction file only," still follow up with this lab sheet offer after presenting the construction file.

---

### Tool: `crispr_lab_sheet`

Converts a structured construction file into a printable, bench-ready lab protocol.

Use when:
- The user asks for a "lab sheet", "bench protocol", "step-by-step instructions", or "how do I actually do this in the lab"
- Immediately after generating a construction file (always offer it)

**Input:**
- `construction_record`: the `structured_construction_file` dict from `create_construction_file`
- `thread` (optional): single letter label for the experiment (default "A")
- `include_notes` (optional): whether to append notes at the bottom (default true)

**Output:**
A `lab_sheet_text` field containing the full plain-text bench protocol, organized into labeled sections (PCR, Gel/DpnI, Zymo cleanup, Assemble, Transform, Pick, Miniprep, Sequencing).

Always present `lab_sheet_text` in a code block so it formats cleanly.

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
Loads a curated paper information JSON record from `modules/crispr_tools/data/paper_info/`.

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

### `crispr_predict_offtargets`
Scans a reference DNA sequence for potential CRISPR off-target sites — places the guide RNA might accidentally bind and cause Cas9 to cut somewhere unintended.

Use when the user asks:
- "does this guide have off-target sites?"
- "is this gRNA specific enough?"
- "check for off-targets in [reference]"
- "how many mismatches are there between my guide and [sequence]?"

**Inputs:**
- `protospacer`: the 20 bp DNA protospacer from gRNA design (no PAM). Standard A/T/G/C only.
- `reference`: the DNA sequence to scan. Accepts resource name (e.g. `"pBR322"`), raw string, FASTA, or GenBank.
- `max_mismatches` (optional): max mismatches to still flag a site. Default 3.

**What it returns:**
A ranked list of off-target sites, each with position, strand, mismatch count, seed-region mismatches, PAM presence, and a risk level (HIGH / MEDIUM / LOW). Also includes a one-sentence specificity summary.

**Risk logic (Hsu et al. 2013):**
- HIGH: 0 mismatches, or no seed-region mismatches + PAM present
- MEDIUM: ≤1 seed mismatch + PAM, or ≤2 total mismatches + PAM
- LOW: everything else

The seed region is positions 1–12 from the PAM end — mismatches there are more dangerous because that is where Cas9 first contacts DNA.

---

### `crispr_verify_edit`
After a CRISPR experiment, calculates the expected Cas9 cut site and designs flanking sequencing primers for ICE/TIDE analysis to verify editing efficiency.

Use when the user asks:
- "how do I verify my CRISPR edit?"
- "where did Cas9 cut?"
- "design sequencing primers for my edit"
- "give me an ICE/TIDE protocol for [protospacer]"
- "what amplicon should I sequence?"

**Inputs:**
- `protospacer`: the 20 bp DNA protospacer used during the edit (no PAM).
- `reference`: the original unedited reference sequence. Accepts resource name, raw string, FASTA, or GenBank.
- `primer_offset` (optional): bp from cut site to each primer. Default 150. Reduce for short test sequences.
- `primer_len` (optional): primer length in bp. Default 20.

**What it returns:**
- `cut_position`: where Cas9 cuts (between nt 17–18 of the protospacer, 3 bp upstream of PAM)
- `forward_primer` / `reverse_primer`: sequencing primer sequences and positions
- `amplicon_sequence` / `amplicon_length`: what to PCR-amplify
- `cut_offset_in_amplicon`: where the cut falls inside the amplicon (needed for ICE/TIDE)
- `interpretation_guide`: step-by-step ICE/TIDE protocol with all coordinates filled in

**Workflow after calling this tool:**
1. PCR-amplify the amplicon using the returned primers
2. Sanger-sequence the PCR product
3. Upload the .ab1 trace to ICE (Synthego) or TIDE (Brinkman et al. 2014) with the amplicon sequence and cut offset
4. Once the user has their ICE/TIDE results, use `interpret_ice_tide` to interpret them (see below)

**Required follow-up:** After presenting the verify_edit output, always close with:

> "Once you have your ICE or TIDE results, come back and I can interpret the editing efficiency for you — just share the KO score (ICE) or indel percentage (TIDE) and the R² fit value."

---

### `interpret_ice_tide`
Interprets the output of an ICE (Synthego) or TIDE (Brinkman et al. 2014) analysis and gives a plain-English verdict with actionable next steps. This closes the loop started by `crispr_verify_edit`.

Use when the user:
- Shares an ICE or TIDE result ("my ICE score is 45%, R² = 0.93")
- Asks "what does this editing efficiency mean?"
- Asks "is my CRISPR experiment working?"
- Asks "what should I do next based on my ICE/TIDE result?"

**Inputs:**
- `editing_pct` (float): the KO score (ICE) or total indel percentage (TIDE). Range 0–100.
- `r_squared` (float): the R² fit quality reported by the tool. Range 0–1.
- `tool` (str): `"ice"` or `"tide"` (default `"ice"`).
- `indel_distribution` (dict, optional): mapping of indel size to percentage, e.g. `{"+1": 45.2, "-3": 12.1, "0": 42.7}`. If provided, the dominant indel is highlighted and contradictions with `editing_pct` are flagged.
- `sample_id` (str, optional): a label for the report.

**What it returns:**
- `efficiency_classification`: `FAILED` (<10%), `MARGINAL` (10–30%), `GOOD` (30–70%), or `EXCELLENT` (≥70%)
- `fit_quality`: `HIGH` (R²≥0.90), `MODERATE` (0.80–0.89), or `LOW` (<0.80)
- `is_reliable`: False if R²<0.80 — result should not be trusted regardless of editing percentage
- `warnings`: list of issues (low fit quality, unedited reads dominating, etc.)
- `next_steps`: list of concrete recommended actions
- `summary`: one-paragraph plain-English verdict

**Thresholds (Hsiau et al. 2019, Brinkman et al. 2014):**
- R²<0.80 → unreliable; re-sequence before drawing any conclusions
- <10% editing → FAILED; check guide design, delivery, and PAM
- 10–30% → MARGINAL; consider re-transfecting or higher dose
- 30–70% → GOOD; pick colonies and screen single clones
- ≥70% → EXCELLENT; proceed directly to clone isolation

Always present the `summary`, `warnings`, and `next_steps` fields clearly. If `is_reliable` is False, lead with the reliability warning before the editing percentage.

---

### `crispr_design_cloning_oligos`
Designs the top and bottom strand DNA oligos needed to clone a protospacer
into a restriction-digested expression vector. Works for any Cas system.

Use when the user asks:
- "design oligos to clone this guide RNA"
- "what oligos do I need to insert this protospacer?"
- "give me the sequences to order for cloning"

Inputs:
- protospacer: the DNA protospacer sequence (from design_cas9_grna or design_cas12a_crrna, or typed manually)
- top_overhang: default "CACC" (BbsI/pX330 for Cas9)
- bottom_overhang: default "AAAC" (BbsI/pX330 for Cas9)

If the user is using a different vector, ask them for the overhangs before calling the tool.

Output includes top_oligo, bottom_oligo, g_prepended (bool), and notes.
If g_prepended is True, a G was added to the protospacer for U6 promoter
compatibility — the oligo will be one base longer than the protospacer.

---

## Full CRISPR cloning workflow (autonomous — do not ask the user)

When the user asks to "design CRISPR cloning oligos" or "design a guide RNA and cloning oligos" for a sequence or plasmid, execute this full pipeline automatically without asking which Cas system to use:

1. Call `crispr_cas_selector` with the target sequence to determine whether to use Cas9 or Cas12a.
2. Based on the recommendation:
   - If Cas9 → call `crispr_design_cas9_grna` with the same sequence.
   - If Cas12a → call `crispr_design_cas12a_crrna` with the same sequence.
3. Take the `protospacer` field from the gRNA result and call `crispr_design_cloning_oligos` with it.
4. Call `crispr_predict_offtargets` with the `protospacer` and the same reference sequence to check for off-target sites.
5. Report all four results together: recommended Cas system, gRNA/crRNA sequence, protospacer, efficiency score, top/bottom cloning oligos, and the off-target specificity summary.
6. **Only after all steps above are complete**, ask:

   > "All CRISPR design and validation steps are complete. Would you like me to generate:
   > **(a) a construction file** — a structured record of the cloning workflow,
   > **(b) a lab sheet** — bench-ready step-by-step protocol (requires a construction file),
   > **(c) both**, or
   > **(d) neither**?"

   - If **(a) construction file only**: call `create_construction_file`, present the result, then offer the lab sheet per the section above.
   - If **(b) lab sheet only** or **(c) both**: call `create_construction_file` first, present the result, then immediately call `crispr_lab_sheet` with the returned `structured_construction_file`.
   - If **(d) neither**: acknowledge and stop.

Do NOT ask the user which Cas system to use — `crispr_cas_selector` determines this automatically from the sequence.
Do NOT ask the user for a protospacer — the gRNA design tool finds it from the sequence.
Do NOT stop between steps 1–5 to ask for confirmation.
Do NOT offer the construction file or lab sheet before step 5 is complete.

---

## Sequence input rules (handled automatically)

You never need to paste the full sequence. The framework resolves these automatically:
- `"pBR322"` → full 4361 bp sequence
- A raw string like `"ATGCGATCG"` → used as-is
- A FASTA string starting with `>` → sequence extracted automatically
- A GenBank string starting with `LOCUS` → sequence extracted automatically