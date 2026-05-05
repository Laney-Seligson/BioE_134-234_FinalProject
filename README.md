# BioE 134/234 Final Project — CRISPR Biodesign Automation Pipeline

A Python-based MCP (Model Context Protocol) server that exposes CRISPR design and cloning tools to an AI assistant (Gemini). Each tool follows the Function Object Pattern (`initiate` / `run`) and is auto-registered as a callable MCP endpoint.

---

## Team Pipeline Overview

```
DNA sequence input
      ↓
Sequence Analyser (auto-annotate, GC content, PAM sites)
      ↓
gRNA Design ──────────────── Off-target Prediction
(PAM scanning, Doench/        (seed-region mismatch scoring,
 Zetsche scoring)              risk ranking)
      ↓
Sequence-based Cas Selector
(Cas-system choice comes from target DNA)
      ↓
Construction File ─────────── LabPlanner Output
(PCR · digest · ligation       (LabSheet + LabPacket)
 steps as JSON)
      ↓
Edit Verification
(cut site · sequencing primers · ICE/TIDE protocol)
```

---

## Implemented Tools

| Tool | File | Author |
|------|------|--------|
| `dna_reverse_complement` | `seq_basics/tools/reverse_complement.py` | |
| `dna_translate` | `seq_basics/tools/translate.py` | |
| `crispr_cas_selector` | `crispr_tools/tools/cas_selector.py` | |
| `crispr_design_cas9_grna` | `crispr_tools/tools/design_cas9_grna.py` | |
| `crispr_design_cas12a_crrna` | `crispr_tools/tools/design_cas12a_crrna.py` | |
| `crispr_predict_offtargets` | `crispr_tools/tools/predict_offtargets.py` | Jillian Ho |
| `crispr_verify_edit` | `crispr_tools/tools/verify_edit.py` | Jillian Ho |
| `create_construction_file` | `crispr_tools/tools/create_construction_file.py` | |
| `validate_construction_file` | `crispr_tools/tools/construction_file_validation.py` | |

---

## Individual Scope — Jillian Ho

**Tools:** `crispr_predict_offtargets` and `crispr_verify_edit`

### crispr_predict_offtargets

**Problem it solves:** After designing a guide RNA, there is no guarantee it only cuts at the intended site. Any genomic region sufficiently similar to the protospacer could be cut too, causing unintended edits. This tool scans a reference sequence on both strands and flags every potential off-target site before the experiment is run.

**How it works:**
- Slides a 20 bp window across every position of the reference (forward + reverse complement strand)
- Counts mismatches between each window and the protospacer
- Applies **seed-region logic** (Hsu et al. 2013): mismatches in positions 1–12 from the PAM end are weighted more heavily because that's where Cas9 initially contacts DNA
- Checks for an NGG PAM immediately after each candidate window
- Assigns a risk level — **HIGH** (0 mismatches, or no seed mismatches + PAM present), **MEDIUM** (1–2 seed mismatches + PAM), **LOW** (everything else)
- Returns all sites ranked by risk, plus a one-sentence specificity summary

**Inputs:** protospacer (20 bp, no PAM), reference sequence (resource name or raw string), optional `max_mismatches` (default 3)

**Output:** ranked list of off-target sites with position, strand, mismatch count, seed mismatch count, PAM presence, and risk level

**MCP wrapper:** `predict_offtargets.json` — `seq_params: ["reference"]` causes the framework to auto-resolve resource names like `pBR322` before the tool runs

---

### crispr_verify_edit

**Problem it solves:** After a CRISPR experiment, researchers need to confirm the edit actually happened. This tool takes the protospacer and the original reference sequence and produces everything needed to run ICE or TIDE analysis on Sanger sequencing data.

**How it works:**
1. Finds the protospacer in the reference on the forward or reverse strand
2. Verifies the NGG PAM is present after it
3. Calculates the Cas9 cut position — always between nucleotides 17 and 18 of the protospacer (3 bp upstream of the PAM)
4. Designs a forward sequencing primer ~150 bp upstream of the cut site
5. Designs a reverse sequencing primer ~150 bp downstream
6. Extracts the expected amplicon sequence between the two primers
7. Returns a step-by-step ICE/TIDE protocol with all coordinates filled in

**Inputs:** protospacer (20 bp, no PAM), reference sequence (resource name or raw string), optional `primer_offset` (default 150 bp), optional `primer_len` (default 20 bp)

**Output:** cut position, PAM sequence, forward/reverse primer sequences and positions, amplicon sequence and length, cut offset within amplicon, full ICE/TIDE interpretation guide

**MCP wrapper:** `verify_edit.json` — `seq_params: ["reference"]` auto-resolves resource names

---

## Test Prompts

See `modules/crispr_tools/tools/prompts.json` for the full list. Prompts covering Jillian's tools:

- `"Check if the guide TCAGAAACCTGCCAGTTTGC has any off-target sites in pBR322."` → `crispr_predict_offtargets`
- `"Are there any off-target sites for ATGATGATGATGATGATGAT in the sequence ATGATGATGATGATGATGATAGG with up to 2 mismatches?"` → `crispr_predict_offtargets`
- `"I edited pBR322 using the guide TCAGAAACCTGCCAGTTTGC. How do I verify the edit worked?"` → `crispr_verify_edit`
- `"Design an ICE/TIDE verification plan for protospacer TCAGAAACCTGCCAGTTTGC in the sequence CCCTAGATGCCTGGCTCAGAAACCTGCCAGTTTGC..."` → `crispr_verify_edit`

---

## Pytests

```bash
cd 2026-bioe234-final-project
pytest tests/test_tools.py -vv
```

Tests covering Jillian's tools (`tests/test_tools.py`):

| Test | What it checks |
|------|---------------|
| `test_predict_offtargets_exact_match_high_risk` | Exact match + NGG PAM → HIGH risk |
| `test_predict_offtargets_no_match_returns_empty` | No similar sites → no results within threshold |
| `test_predict_offtargets_empty_protospacer_raises` | Empty protospacer → ValueError |
| `test_predict_offtargets_empty_reference_raises` | Empty reference → ValueError |
| `test_predict_offtargets_max_mismatches_zero` | `max_mismatches=0` → only perfect matches returned |
| `test_verify_edit_standard_case` | Correct cut site (pos 32), PAM (TGG), strand (+) |
| `test_verify_edit_protospacer_not_in_reference_raises` | Protospacer not in reference → ValueError |
| `test_verify_edit_empty_protospacer_raises` | Empty protospacer → ValueError |
| `test_verify_edit_empty_reference_raises` | Empty reference → ValueError |

---

## Running the MCP Server

```bash
cd 2026-bioe234-final-project
python server.py
```

## Running the Gemini Client

```bash
cd 2026-bioe234-final-project
python client_gemini.py
```
