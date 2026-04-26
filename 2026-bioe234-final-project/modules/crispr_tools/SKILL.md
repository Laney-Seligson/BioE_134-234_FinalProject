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

**Required behavior after `crispr_lab_sheet`:** Once the lab sheet is presented, always ask:

> "Would you like a colony picking estimate? After transformation, you'll need to screen enough colonies to be confident of recovering an edited clone — the number depends on your expected editing efficiency and delivery method. I can calculate how many to pick."

If the user says yes (or anything like "sure", "yes", "go ahead"):
- Call `crispr_colony_calculator` with the appropriate preset (infer from context — see tool section below)
- Present the result clearly: efficiency used, colonies to pick, safety margin, and the plain-English recommendation

Then close with:

> "Once you've picked colonies, run a diagnostic PCR, and sent samples for Sanger sequencing, come back and I'll design your sequencing primers and help you interpret your ICE/TIDE results."

If the user says no or skips: close with just the ICE/TIDE handoff prompt above.

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

### `crispr_colony_calculator`
Estimates how many colonies to pick to be 95% confident of recovering at least one edited clone, using the binomial distribution and published efficiency benchmarks. Call this **automatically after every lab sheet** — it bridges the cloning protocol to the colony screening step.

Use when:
- After presenting a lab sheet (always — do not wait for the user to ask)
- The user asks "how many colonies should I pick?"
- The user asks "how confident are we this will work?" or "what are the odds of success?"

**Inputs:**
- `preset` (str, preferred): one of —
  - `cas9_plasmid_mammalian` — 20% (plasmid transfection, e.g. pX330)
  - `cas9_rnp_mammalian` — 65% (RNP / electroporation; Kim et al. 2014)
  - `cas9_ecoli` — 50% (E. coli; Jiang et al. 2013)
  - `cas12a_ecoli` — 40%
  - `hdr_mammalian` — 5% (precise knock-in via HDR; Paquet et al. 2016)
  - `hdr_ecoli` — 10%
- `editing_efficiency` (float 0–1): use instead of preset if the user provides their own estimate
- `desired_clones` (int, default 1): how many edited clones they want to recover
- `confidence` (float, default 0.95): probability threshold

**How to choose preset — infer from context:**
- pX330 or any plasmid vector → `cas9_plasmid_mammalian`
- User mentions RNP, electroporation, nucleofection → `cas9_rnp_mammalian`
- Bacterial / E. coli target → `cas9_ecoli`
- Precise knock-in / HDR → `hdr_mammalian` or `hdr_ecoli`
- Unclear → ask "are you delivering by plasmid transfection or RNP?"

**Preferred: use `predict_editing_efficiency` instead of a preset when available.**
If `crispr_predict_editing_efficiency` has already been called for this guide, pass its `on_target_efficiency_pct` (divided by 100) as `editing_efficiency` directly — this is more accurate than a generic preset. If not yet called, fall back to the preset table above.

**Fallback — connecting guide quality score to preset range:**
The `total_score` from `crispr_rank_guides` (max 6 for Cas9) indicates where within the preset range efficiency will likely fall:
- Score 5–6 → upper end of preset range (e.g. ~25–30% for plasmid mammalian)
- Score 3–4 → mid-range (~15–20%)
- Score ≤2 → lower end; consider redesigning the guide before running the experiment

Always mention the guide score or efficiency estimate alongside the colony count so the user understands it is a conditional prediction.

**What it returns:**
- `colonies_to_pick`: minimum colonies for the desired confidence
- `safety_margin_recommendation`: 1.5× bump for real-world losses (failed PCR, contamination)
- `recommendation`: plain-English advice

---

### `crispr_predict_editing_efficiency`
Predicts on-target editing efficiency **before the experiment** using a simplified Doench 2016 Rule Set 2 model. Returns a predicted % efficiency and a ±15% confidence range. Use this to give the user a guide-specific estimate — more accurate than the generic delivery presets in `colony_calculator`.

Use when:
- The user asks "how well will this guide work?" or "what efficiency do we expect?"
- During the full design workflow — call after `crispr_rank_guides` to give a pre-experiment efficiency prediction for the best guide
- When offering colony picking guidance — call this first to get a guide-specific efficiency, then pass it to `colony_calculator`

**Inputs:**
- `protospacer` (str): 20 bp (Cas9) or 23 bp (Cas12a). No PAM.
- `pam` (str): the PAM adjacent to the protospacer (e.g. `"AGG"` for Cas9, `"TTTA"` for Cas12a).
- `nuclease` (str): `"cas9"` (default) or `"cas12a"`.
- `downstream_3nt` (str, optional): 3 nt immediately after the PAM — used for NGGT vs NGGA PAM preference (Cas9 only).
- `delivery` (str, default `"plasmid"`): one of `"rnp"`, `"plasmid"`, `"lentivirus"`, `"aav"`, `"electroporation"`.
- `outcome` (str, default `"nhej"`): one of `"nhej"` (knockout), `"hdr"` (knockin), `"base_edit"`, `"prime_edit"`.

**What it returns:**
- `on_target_efficiency_pct`: predicted % editing after delivery/outcome adjustment
- `confidence_range`: [low, high] as ±15 percentage points
- `feature_contributions`: breakdown of what helped/hurt the score (position weights, PAM context, GC, poly-T)
- `interpretation`: plain-English verdict (high / moderate / low / very low)
- `warnings`: red flags like poly-T runs, weak PAM, GC out of range

**How efficiency feeds into colony picking:**
Pass `on_target_efficiency_pct / 100` as `editing_efficiency` to `colony_calculator` for a guide-specific colony count rather than a generic preset.

**Caveat:** This is a simplified linear approximation of Doench Rule Set 2, not the full Azimuth/CRISPOR model. Treat predictions as estimates, not guarantees — especially for HDR where cell type and donor design dominate.

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
A ranked list of off-target sites, each with position, strand, mismatch count, seed-region mismatches, PAM presence, risk level (HIGH / MEDIUM / LOW), and a `cfd_score`. Two aggregate fields: `aggregate_offtarget_cfd` (sum of all off-target CFD scores; lower = more specific guide) and `max_offtarget_cfd` (worst single off-target; above 0.1 warrants concern). Also includes a one-sentence specificity summary.

**Risk logic (Hsu et al. 2013):**
- HIGH: 0 mismatches, or no seed-region mismatches + PAM present
- MEDIUM: ≤1 seed mismatch + PAM, or ≤2 total mismatches + PAM
- LOW: everything else

The seed region is positions 1–12 from the PAM end — mismatches there are more dangerous because that is where Cas9 first contacts DNA.

**CFD score (Doench 2016):** Each off-target's `cfd_score` is the predicted cutting frequency relative to the on-target (0 = no cutting, 1 = same as on-target). Sites without a PAM have CFD = 0. Use `max_offtarget_cfd` as a quick specificity flag: >0.1 means at least one off-target could be cut at >10% the on-target rate.

---

### `crispr_verify_edit`
Calculates the expected Cas9 cut site and designs flanking sequencing primers for ICE/TIDE analysis. Call this **before the experiment** — the primers need to be ordered alongside the cloning oligos so the user has them ready when it's time to sequence.

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

Note: `crispr_verify_edit` designs the **sequencing primers** — this is prep work done before the experiment so the primers are ready to order. The actual interpretation of results is done later with `interpret_ice_tide` after sequencing data is in hand. Do NOT call `crispr_verify_edit` during the automated design workflow (steps 1–7 above) — it is offered as a standalone step when the user asks to set up their verification protocol.

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
3. Call `crispr_rank_guides` with the full `guides` list returned in step 2 and the same reference sequence. This scores all candidates on efficiency (GC content, no poly-T, PAM-proximal G) and specificity (off-target risk), then selects the best one.
4. Call `crispr_predict_editing_efficiency` with `best_guide.protospacer` and its PAM to get a guide-specific predicted efficiency % (adjust `delivery` and `outcome` from context if known, otherwise use defaults: `delivery="plasmid"`, `outcome="nhej"`).
5. Take `best_guide.protospacer` from the ranking result and call `crispr_design_cloning_oligos` with it.
6. Call `crispr_predict_offtargets` with `best_guide.protospacer` and the same reference sequence to get the full off-target site list including CFD scores.
7. Report all results together. **Always include:**
   - The `scoring_rationale` from `crispr_rank_guides` (why this guide was selected)
   - The `interpretation` and `on_target_efficiency_pct` from `crispr_predict_editing_efficiency` (pre-experiment efficiency estimate)
   - The `max_offtarget_cfd` from `crispr_predict_offtargets` (specificity flag — flag if >0.1)
7. **Only after all steps above are complete**, ask:

   > "All CRISPR design and validation steps are complete. Would you like me to generate:
   > **(a) a construction file** — a structured record of the cloning workflow,
   > **(b) a lab sheet** — bench-ready step-by-step protocol (requires a construction file),
   > **(c) both**, or
   > **(d) neither**?"

   - If **(a) construction file only**: call `create_construction_file`, present the result, then offer the lab sheet per the section above.
   - If **(b) lab sheet only** or **(c) both**: call `create_construction_file` first, present the result, then immediately call `crispr_lab_sheet` with the returned `structured_construction_file`. After presenting the lab sheet, follow the required behavior described in the `crispr_lab_sheet` section above (offer colony picking, then close with the sequencing/ICE-TIDE handoff).
   - If **(d) neither**: acknowledge and stop.

Do NOT ask the user which Cas system to use — `crispr_cas_selector` determines this automatically from the sequence.
Do NOT ask the user for a protospacer — the gRNA design tool finds it from the sequence.
Do NOT stop between steps 1–7 to ask for confirmation.
Do NOT offer the construction file or lab sheet before step 7 is complete.
Do NOT call `crispr_verify_edit` during this workflow — it is a post-experiment tool used only after the physical experiment has been run.
Do NOT offer ICE/TIDE interpretation during this workflow — it belongs in a completely separate conversation after the user has sequencing results in hand.

---

## Sequence input rules (handled automatically)

You never need to paste the full sequence. The framework resolves these automatically:
- `"pBR322"` → full 4361 bp sequence
- A raw string like `"ATGCGATCG"` → used as-is
- A FASTA string starting with `>` → sequence extracted automatically
- A GenBank string starting with `LOCUS` → sequence extracted automatically