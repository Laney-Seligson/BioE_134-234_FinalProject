# crispr_tools — Skill Guidance for Gemini

This file is read by the client at startup and injected into Gemini's system prompt.
Its purpose is to give Gemini the domain knowledge it needs to use the tools in this
module correctly and interpret their results meaningfully.

---

## What this module does

The `crispr_tools` module provides fundamental tools to go through the CRISPR pipeline,
including generating and validating cloning construction workflows.

---

## Available resources

Here are the descriptions of the available resources. When a user inquires for a plasmid
or a backbone sequence, provide the names of the plasmid and a short description to help
the user choose.

| Resource name | Description                                                                                                                           |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `pBR322`      | E. coli cloning vector pBR322, 4361 bp, circular, double-stranded. Contains ampicillin (bla) and tetracycline (tet) resistance genes. |
| `pET28a`      | E. coli cloning vector pET28a, 5369 bp, circular, double-stranded. Contains kanamycin resistance and a 6x His tag.                    |

When a user refers to these plasmids, use the resource name directly (e.g., `"pET28a"`).
Do NOT ask the user to paste sequences.

---

## Tools and when to use them

### `dna_reverse_complement`

Returns the reverse complement of a DNA or RNA sequence.

Use when the user asks for:

* reverse complement
* antisense strand
* flipping a sequence

---

### `dna_translate`

Translates DNA into protein.

Use when the user asks to:

* translate a sequence
* get a protein sequence
* evaluate reading frames

---

## Tool: create_construction_file

This tool generates a cloning construction workflow.

---

### Required inputs

For all workflows:

* construct_name
* host_organism
* backbone_name
* backbone_sequence
* insert_name
* insert_sequence

For GoldenGate:

* insert_forward_primer_name
* insert_forward_primer_sequence
* insert_reverse_primer_name
* insert_reverse_primer_sequence
* vector_forward_primer_name
* vector_forward_primer_sequence
* vector_reverse_primer_name
* vector_reverse_primer_sequence
* enzyme

Optional:

* cell_strain
* selection
* temperature_c


## 🚨 CRITICAL EXECUTION RULES (MANDATORY)

The assistant must ALWAYS take exactly ONE action per turn:

1. Ask for missing required fields
2. Call `create_construction_file`

The assistant must NEVER:

* Return an empty response
* Return a placeholder response
* Wait for the user to say "proceed"
* Delay execution after all fields are present


## Required interaction behavior

If any required fields are missing:

* Ask for ALL missing fields in ONE message
* Do NOT ask multiple follow-ups
* Do NOT call the tool yet


## Required execution behavior

If ALL required fields are present:

* Call `create_construction_file` IMMEDIATELY
* Do NOT ask for confirmation
* Do NOT wait for another message
* Do NOT return text instead of calling the tool

---

## Fallback rule (VERY IMPORTANT)

If unsure whether all required fields are present:

* ASSUME they are present
* CALL the tool anyway

It is ALWAYS better to attempt a tool call than to return an empty response.

---

## Mandatory example behavior

This is REQUIRED behavior:

User provides missing fields:

Assistant MUST:

* Recognize all required fields are now present
* Immediately call `create_construction_file`
* Return the result

Assistant MUST NOT:

* Return an empty response
* Ask for confirmation
* Wait for "proceed"
* Skip the tool call

Failure to call the tool when inputs are complete is incorrect behavior.

---

## Presenting results

After calling `create_construction_file`:

* If `construction_file_txt` is returned:

  * Display it directly in a code block

* Do NOT show raw JSON

* Do NOT show MCP output format

Correct format:

```
PCR           ...
PCR           ...
GoldenGate    ...

plasmid       ...
dsdna         ...
oligo         ...
```

---

## Tool: validate_construction_file

Validates whether a construction workflow is biologically correct.

---

### When to use

Use when the user asks:

* "is this valid?"
* "does this PCR work?"
* "check my cloning workflow"
* "verify this assembly"

Also use after generating a construction file IF the user asks for confirmation.

### `crispr_cas_selector`
Analyzes the GC and AT content of a DNA sequence and recommends Cas9 or Cas12a.

Use when the user asks:
- "which CRISPR system should I use for this organism?"
- "should I use Cas9 or Cas12a?"
- "what is the GC content of this genome?"
- "is this sequence GC-rich or AT-rich?"

Output includes gc_fraction, at_fraction, recommendation ("Cas9" or "Cas12a"),
and a one-sentence rationale. GC >= 50% → Cas9 (NGG PAM abundant);
GC < 50% → Cas12a (TTTV PAM abundant).

---

### `crispr_design_cas9_grna`
Scans a DNA sequence for all NGG PAM sites, scores each candidate protospacer
using Doench et al. 2016 efficiency rules, and returns the optimal Cas9 gRNA.

Use when the user asks:
- "design a Cas9 gRNA for this sequence"
- "what is the best guide RNA for Cas9?"
- "find a Cas9 target site in pBR322"

Scoring rules (each worth 1 point, max score 3):
- GC content of protospacer between 40-70%
- No poly-T run (TTTT+) in protospacer
- Final base of protospacer (position 20) is G

Output includes the full RNA gRNA sequence, protospacer, PAM site,
efficiency_score (0-3), and any warnings.

---

### `crispr_design_cas12a_crrna`
Scans a DNA sequence for all TTTV PAM sites, scores each candidate protospacer
using Zetsche et al. 2015 and Pausch et al. 2020 efficiency rules, and returns
the optimal FnCas12a crRNA.

Use when the user asks:
- "design a Cas12a crRNA for this sequence"
- "what is the best guide RNA for Cas12a?"
- "find a Cas12a target site in pBR322"

Key difference from Cas9: PAM is TTTV (not NGG), protospacer is 23 bp and lies
DOWNSTREAM of the PAM (not upstream), and the scaffold is a direct repeat
(not tracrRNA).

Scoring rules (each worth 1 point, max score 3):
- GC content of protospacer between 30-70%
- No poly-T run (TTTT+) in protospacer
- First base of protospacer (PAM-proximal) is C or G

Output includes the full RNA crRNA sequence, protospacer, PAM site,
efficiency_score (0-3), and any warnings.

---

## Required execution rules (validation)

If user requests validation:

* Call `validate_construction_file` immediately

Do NOT:

* Return an empty response
* Delay execution
* Ask unnecessary follow-up questions

---

## What validation checks (Version 1)

* PCR primer annealing
* Primer orientation
* Valid amplicon formation

---

## Not yet implemented

* GoldenGate validation (future versions)
* Gibson overlaps
* enzyme cut simulation

These steps appear as `[SKIP]`

---

## Interpreting validation results

* `[PASS]` → valid
* `[FAIL]` → incorrect, must fix
* `[SKIP]` → not implemented

---

## Sequence input rules

Sequences are automatically resolved:

* `"pET28a"` → full plasmid sequence
* raw strings → used directly
* FASTA / GenBank → parsed automatically

No need to paste full sequences.
