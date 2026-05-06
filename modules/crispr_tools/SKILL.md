---
name: crispr_tools
description: CRISPR pipeline tools — guide RNA design, off-target prediction, cloning oligo design, and editing efficiency prediction.
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
- MEDIUM: <=1 seed mismatch + PAM, or <=2 total mismatches + PAM
- LOW: everything else

The seed region is positions 1-12 from the PAM end — mismatches there are more dangerous because that is where Cas9 first contacts DNA.

**CFD score (Doench 2016):** Each off-target's `cfd_score` is the predicted cutting frequency relative to the on-target (0 = no cutting, 1 = same as on-target). Sites without a PAM have CFD = 0. Use `max_offtarget_cfd` as a quick specificity flag: >0.1 means at least one off-target could be cut at >10% the on-target rate.

---

### `crispr_predict_editing_efficiency`
Predicts on-target editing efficiency **before the experiment** using a simplified Doench 2016 Rule Set 2 model. Returns a predicted % efficiency and a +-15% confidence range. Use this to give the user a guide-specific estimate — more accurate than the generic delivery presets in `colony_calculator`.

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
- `confidence_range`: [low, high] as +-15 percentage points
- `feature_contributions`: breakdown of what helped/hurt the score (position weights, PAM context, GC, poly-T)
- `interpretation`: plain-English verdict (high / moderate / low / very low)
- `warnings`: red flags like poly-T runs, weak PAM, GC out of range

**How efficiency feeds into colony picking:**
Pass `on_target_efficiency_pct / 100` as `editing_efficiency` to `colony_calculator` for a guide-specific colony count rather than a generic preset.

**Caveat:** This is a simplified linear approximation of Doench Rule Set 2, not the full Azimuth/CRISPOR model. Treat predictions as estimates, not guarantees — especially for HDR where cell type and donor design dominate.

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

## When to use `crispr_run_full_workflow` vs individual tools

**Use `crispr_run_full_workflow` ONLY when the user's request is explicitly a one-shot full workflow** — meaning they name a gene, an organism, AND a vector in the same message (e.g. "Design a CRISPR edit targeting lacZ in E. coli using pTargetF").

**Use the individual tools (steps below) for ALL other requests**, including:
- "Design Cas9 guides for X" → fetch → design_cas9_grna → rank_guides
- "Fetch the sequence for X" → fetch only
- "Design guides for X" → fetch → cas_selector → design_cas9_grna or design_cas12a_crrna → rank_guides
- "Rank these guides" → rank_guides only
- Any request that does not name a specific vector upfront

Calling `crispr_run_full_workflow` when the user only asked for guides hides the intermediate results (guide list, ranking rationale) that the user needs to see.

**When the user names a gene vaguely or asks for "the most researched" / "a common" gene:** search for candidates using `semantic_gene_search` or `go_term_gene_lookup`, then present the options to the user and ask them which gene to target. Do NOT pick one automatically.

**When `crispr_run_full_workflow` returns `status: "needs_user_input"`:** present the `questions` field verbatim to the user and the `vector_recommendations` list (name + use_case for each). Wait for the user to choose. Do NOT select a vector from the recommendations yourself — always ask the user to choose.

**The workflow result does NOT include the raw target sequence** — it is stripped to avoid truncation artifacts. If a downstream tool (e.g. `crispr_verify_edit`) needs the reference sequence, call `crispr_fetch_target_sequence` again with the same query and organism from `sequence_info`.

---

## Full CRISPR cloning workflow (autonomous — do not ask the user)

When the user asks to "design CRISPR cloning oligos" or "design a guide RNA and cloning oligos" for a sequence or plasmid, execute this full pipeline automatically without asking which Cas system to use:

1. Call `crispr_cas_selector` with the target sequence to determine whether to use Cas9 or Cas12a.
2. Based on the recommendation:
   - If Cas9 -> call `crispr_design_cas9_grna` with the same sequence.
   - If Cas12a -> call `crispr_design_cas12a_crrna` with the same sequence.
3. Call `crispr_rank_guides` with the full `guides` list returned in step 2 and the same reference sequence. This scores all candidates on efficiency (GC content, no poly-T, PAM-proximal G) and specificity (off-target risk), then selects the best one.
4. Call `crispr_predict_editing_efficiency` with `best_guide.protospacer` and its PAM to get a guide-specific predicted efficiency % (adjust `delivery` and `outcome` from context if known, otherwise use defaults: `delivery="plasmid"`, `outcome="nhej"`).
5. Take `best_guide.protospacer` from the ranking result and call `crispr_design_cloning_oligos` with it.
6. Call `crispr_predict_offtargets` with `best_guide.protospacer` and the same reference sequence to get the full off-target site list including CFD scores.
7. Report all results together. **Always include:**
   - The `scoring_rationale` from `crispr_rank_guides` (why this guide was selected)
   - The `interpretation` and `on_target_efficiency_pct` from `crispr_predict_editing_efficiency` (pre-experiment efficiency estimate)
   - The `max_offtarget_cfd` from `crispr_predict_offtargets` (specificity flag — flag if >0.1)
8. **Only after all steps above are complete**, ask:

   > "All CRISPR design and validation steps are complete. Would you like me to generate:
   > **(a) a construction file** — a structured record of the cloning workflow,
   > **(b) a lab sheet** — bench-ready step-by-step protocol (requires a construction file),
   > **(c) both**, or
   > **(d) neither**?"

   - If **(a) construction file only**: call `create_construction_file`, present the result, then offer the lab sheet.
   - If **(b) lab sheet only** or **(c) both**: call `create_construction_file` first, present the result, then immediately call `crispr_lab_sheet`. After presenting the lab sheet, offer colony picking (see `labsheet_tools` SKILL for colony calculator guidance), then close with the sequencing/ICE-TIDE handoff.
   - If **(d) neither**: acknowledge and stop.

Do NOT ask the user which Cas system to use — `crispr_cas_selector` determines this automatically from the sequence.
Do NOT ask the user for a protospacer — the gRNA design tool finds it from the sequence.
Do NOT stop between steps 1-7 to ask for confirmation.
Do NOT offer the construction file or lab sheet before step 7 is complete.
Do NOT call `crispr_verify_edit` during this workflow — offer it separately after the construction file/lab sheet are done so the user can order sequencing primers alongside cloning oligos.
Do NOT offer ICE/TIDE interpretation during this workflow — it belongs after the user has sequencing results in hand.

---

## Sequence tool priority

For any CRISPR-related request, always use `crispr_fetch_target_sequence` to resolve a gene name to a DNA sequence — it returns a clean, guide-design-ready sequence. `gene_sequence_lookup_tool` and `lookup_gene_sequence` return full genome FASTA records that are too large for guide design and will crash the pipeline.

---

## Sequence input rules (handled automatically)

You never need to paste the full sequence. The framework resolves these automatically:
- `"pBR322"` -> full 4361 bp sequence
- A raw string like `"ATGCGATCG"` -> used as-is
- A FASTA string starting with `>` -> sequence extracted automatically
- A GenBank string starting with `LOCUS` -> sequence extracted automatically
