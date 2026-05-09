![CRISPR](https://img.shields.io/badge/CRISPR-Workflow-F527F2?style=for-the-badge)


## Emory Adelman — CRISPR Toolchain

> **DISCLAIMER**  
> Citations here do not cover every piece of material in the script. The README focuses primarily on benchmark workflows and representative end-user examples, while the scripts themselves provide more granular justification for individual implementation decisions and helper logic.

My tools support CRISPR education and wet-lab protocol design. Every design decision — vector selection, guide scoring, primer logic — is backed by literature citations. The workflow begins when the user provides a gene name: the tool fetches the genomic locus or coding sequence from NCBI, selects a Cas nuclease using a heuristic I developed (the nuclease with more high-quality PAM-compatible guides wins), designs and ranks up to 10 guide sequences, then passes the best guide and all cloning variables to teammates' tools for construction file and lab sheet generation.
      



## 1. [run_full_crispr_workflow](modules/crispr_tools/tools/run_full_crispr_workflow.py)

- What it does: 
  - NEEDS gene and Organism name to successfully call tool 
  - Fetches a target sequence, genomic locus or cds from the NCBI GenBank file for the gene chosen 
  - design guide or crRNA sequence using the 10 found in the fetched sequence, meant to cover the gene broadly from PAM-based windows across the whole sequence 
  - rank guide and crispr RNAs with [rank_guides.py](modules/crispr_tools/tools/rank_guides.py), authored by Jillian and I, and find the single best one to use out of the 10 
  - off-target prediction is calculated with [predict_offtargets.py](modules/crispr_tools/tools/predict_offtargets.py) 
  - design cloning oligos 
  - prepares construction file inputs 

The MCP will then stop and ask if the user wants to generate a construction file or lab sheet, then Laney and Jillian's tools are called. TypeIIS oligo cloning steps (BbsI/BsmBI annealed-oligo ligation into digested vector backbones) follow the standard protocols used in `create_construction_file` and `validate_construction_file` (Addgene/Zhang Lab 2015; NEB FAQ; Hu et al. 2020).

![Sample Output](https://img.shields.io/badge/expand-sample_output-F527F2?style=flat-square)
<details>
<summary>Click to expand sample output</summary>

```
You: Run a CRISPR edit workflow for EMX1 in homo sapiens using px330 

[Tool call] → crispr_run_full_workflow
{
  "vector": "px330",
  "organism": "homo sapiens",
  "query": "EMX1"
}
[Tool result] ← crispr_run_full_workflow:
{
  "status": "needs_user_input",
  "missing_fields": [
    "target_type"
  ],
  "query": "EMX1",
  "organism": "homo sapiens",
  "questions": [
    "Should I fetch the full genomic locus or just the coding sequence (CDS) for 'EMX1' in homo sapiens?"
  ],
  "options": [
    "genomic_locus \u2014 full gene region including introns and flanking sequence (recommended for most CRISPR experiments)",
    "cds \u2014 coding sequence only, mRNA-derived; use when targeting transcribed exons specifically"
  ],
  "continue_with": {
    "query": "EMX1",
    "organism": "homo sapiens",
    "vector": "px330",
    "confirmed": false,
    "source_query": "",
    "gene_confirmed": false,
    "target_type": "<genomic_locus or cds>"
  },
  "workflow_trace": []
}

Gemini: Please confirm whether I should fetch the full genomic locus or just the coding sequence (CDS) for EMX1 in homo sapiens:

*   **genomic_locus**: full gene region including introns and flanking sequence (recommended for most CRISPR experiments)
*   **cds**: coding sequence only, mRNA-derived (use when targeting transcribed exons specifically)

You: genomic locus 

[Tool call] → crispr_run_full_workflow
{
  "vector": "px330",
  "target_type": "genomic_locus",
  "confirmed": true,
  "query": "EMX1",
  "organism": "homo sapiens",
  "gene_confirmed": true
}
[Tool result] ← crispr_run_full_workflow:
{
  "status": "needs_user_input",
  "missing_fields": [
    "generate_docs"
  ],
  "sequence_info": {
    "source": "ncbi",
    "resource": "EMX1",
    "organism": "Homo sapiens",
    "target_type": "genomic_locus",
    "ncbi_gene_id": "2016",
    "ncbi_accession": "NC_000002.12",
    "length": 18703,
    "note": "Fetched genomic locus from NCBI Entrez: gene 'EMX1' in 'Homo sapiens'. Nucleotide record: NC_000002.12. Verify the recor\u2026 (162 chars)"
  },
  "guides": [
    {
      "grna_sequence": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GCTCTCCAAACCACGGCTCC",
      "pam_site": "CGG",
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
      "grna_sequence": "UUUAGAAGUUACUGCCGGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TTTAGAAGTTACTGCCGGGA",
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
      "grna_sequence": "GUCCAGGACUCCGAAGGUCGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTCCAGGACTCCGAAGGTCG",
      "pam_site": "GGG",
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
      "grna_sequence": "ACUACGUGGUGGGCGCCGAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "ACTACGTGGTGGGCGCCGAG",
      "pam_site": "CGG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.7,
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
      "grna_sequence": "CACUUAACCUCUCUGAGCUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CACTTAACCTCTCTGAGCTC",
      "pam_site": "TGG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.5,
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
      "grna_sequence": "GGAAAGAAGGGCAGGGUUUGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAAAGAAGGGCAGGGTTTG",
      "pam_site": "AGG",
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
      "grna_sequence": "GUGUAGGGGCAGCACCUCCUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTGTAGGGGCAGCACCTCCT",
      "pam_site": "GGG",
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
      "grna_sequence": "CAGCUGGUCAGAGGGGACCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CAGCTGGTCAGAGGGGACCC",
      "pam_site": "CGG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.7,
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
      "grna_sequence": "UACAUAGAUGUUUCCAGAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TACATAGATGTTTCCAGAGA",
      "pam_site": "TGG",
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
      "grna_sequence": "GGAGAAUGGAGACUUUUGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAGAATGGAGACTTTTGGA",
      "pam_site": "AGG",
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
  "selected_guide": {
    "grna_sequence": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
    "protospacer": "GCTCTCCAAACCACGGCTCC",
    "pam_site": "CGG",
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
  "scoring_rationale": "Best guide 'GCTCTCCAAACCACGGCTCC' scored 3/3 (efficiency 2/2, specificity 1/1). GC content: 65%. Off-target sites: 0 tot\u2026 (142 chars)",
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "cloning": {
    "status": "ready",
    "cloning_method": "TypeIISOligoCloning",
    "vector": "pX330",
    "enzyme": "BbsI",
    "top_overhang": "CACC",
    "bottom_overhang": "AAAC",
    "top_oligo_name": "px330_guide_top",
    "bottom_oligo_name": "px330_guide_bottom",
    "top_oligo": "CACCGCTCTCCAAACCACGGCTCC",
    "bottom_oligo": "AAACGGAGCCGTGGTTTGGAGAGC",
    "end_structure": "two sticky-end overhangs",
    "g_prepended": false,
    "final_protospacer": "GCTCTCCAAACCACGGCTCC",
    "construction_file_inputs": {
      "construct_name": "px330_guide_construct",
      "assembly_strategy": "TypeIISOligoCloning",
      "backbone_name": "px330",
      "backbone_sequence": "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGA\u2026 (8484 chars)",
      "insert_name": "px330_annealed_guide_insert",
      "insert_sequence": "GCTCTCCAAACCACGGCTCC",
      "insert_forward_primer_name": "",
      "insert_forward_primer_sequence": "",
      "insert_reverse_primer_name": "",
      "insert_reverse_primer_sequence": "",
      "top_oligo_name": "px330_guide_top",
      "top_oligo_sequence": "CACCGCTCTCCAAACCACGGCTCC",
      "bottom_oligo_name": "px330_guide_bottom",
      "bottom_oligo_sequence": "AAACGGAGCCGTGGTTTGGAGAGC",
      "top_overhang": "CACC",
      "bottom_overhang": "AAAC",
      "vector_forward_primer_name": "",
      "vector_forward_primer_sequence": "",
      "vector_reverse_primer_name": "",
      "vector_reverse_primer_sequence": "",
      "enzyme": "BbsI",
      "cell_strain": "Any mammalian",
      "selection": "Amp",
      "temperature_c": 37,
      "notes": "Vector 'pX330' \u2014 BbsI (TypeIIS) cloning. All-in-one SpCas9 + sgRNA vector.  BbsI digest leaves CACC/AAAC overhangs.  5\u2032G\u2026 (603 chars)"
    }
  },
  "questions": [
    "Oligo design is complete. Would you like to also generate a construction file and lab sheet, or return the oligo results\u2026 (126 chars)"
  ],
  "options": [
    "Yes, generate both construction file and lab sheet (generate_docs=true)",
    "No, return oligo results only (generate_docs=false)"
  ],
  "continue_with": {
    "query": "EMX1",
    "organism": "homo sapiens",
    "vector": "px330",
    "confirmed": true,
    "source_query": "",
    "target_type": "genomic_locus",
    "generate_docs": true,
    "gene_confirmed": true
  },
  "workflow_trace": [
    {
      "step": 1,
      "tool": "crispr_fetch_target_sequence",
      "inputs": {
        "query": "EMX1",
        "organism": "homo sapiens",
        "target_type": "genomic_locus"
      },
      "output_summary": {
        "source": "ncbi",
        "resource": "EMX1",
        "organism": "Homo sapiens",
        "length": 18703,
        "target_type": "genomic_locus",
        "ncbi_gene_id": "2016",
        "ncbi_accession": "NC_000002.12"
      }
    },
    {
      "step": 2,
      "tool": "crispr_cas_selector",
      "inputs": {
        "sequence_length": 18703,
        "repair_template": false
      },
      "output_summary": {
        "recommendation": "Cas9",
        "cas9_valid_guides": 1499,
        "cas12a_valid_guides": 172,
        "gc_content": 0.551
      }
    },
    {
      "step": 3,
      "tool": "crispr_design_cas9_grna",
      "inputs": {
        "sequence_length": 18703
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
        "reference_length": 18703,
        "nuclease": "cas9"
      },
      "output_summary": {
        "ranked_guide_count": 10,
        "best_guide": "GCTCTCCAAACCACGGCTCC",
        "best_score": 3
      }
    },
    {
      "step": 5,
      "tool": "crispr_design_cloning_oligos",
      "inputs": {
        "vector": "px330",
        "protospacer": "GCTCTCCAAACCACGGCTCC",
        "construct_name": null
      },
      "output_summary": {
        "status": "ready",
        "vector": "pX330",
        "cloning_method": "TypeIISOligoCloning",
        "enzyme": "BbsI",
        "top_oligo": "CACCGCTCTCCAAACCACGGCTCC",
        "bottom_oligo": "AAACGGAGCCGTGGTTTGGAGAGC",
        "end_structure": "two sticky-end overhangs"
      }
    }
  ]
}

[Tool call] → crispr_predict_editing_efficiency
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "pam": "CGG"
}
[Tool result] ← crispr_predict_editing_efficiency:
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "nuclease": "cas9",
  "raw_score": 5.5,
  "on_target_efficiency_pct": 40.2,
  "confidence_range": [
    25.2,
    55.2
  ],
  "delivery": "plasmid",
  "outcome": "nhej",
  "delivery_multiplier": 0.55,
  "outcome_multiplier": 1.0,
  "feature_contributions": {
    "position_specific_nucleotides": -0.5,
    "pam_context": 5.0,
    "gc_content": 1.0,
    "polyt_penalty": 0.0
  },
  "interpretation": "Predicted NHEJ knockout efficiency: 40.2% (range 25.2-55.2%) for CAS9 delivered via plasmid. This is moderate \u2014 usable, \u2026 (228 chars)",
  "warnings": []
}

[Tool call] → crispr_predict_offtargets
{
  "reference": "NC_000002.12",
  "protospacer": "GCTCTCCAAACCACGGCTCC"
}
[05/08/26 17:28:10] Error calling tool 'crispr_predict_offtargets'              
                    ╭─────────── Traceback (most recent call last) ────────────╮
                    │ /opt/anaconda3/lib/python3.12/site-packages/fastmcp/serv │
                    │ er/server.py:987 in call_tool                            │
                    │                                                          │
                    │ /opt/anaconda3/lib/python3.12/site-packages/fastmcp/tool │
                    │ s/tool.py:354 in _run                                    │
                    │                                                          │
                    │                 ... 7 frames hidden ...                  │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /seq_basics/_plumbing/resolve.py:84 in resolve_to_seq    │
                    │                                                          │
                    │    81 │   │   │   return _clean_sequence(remainder)      │
                    │    82 │                                                  │
                    │    83 │   # 4. Assume it's a raw sequence string - clean │
                    │ ❱  84 │   return _clean_sequence(input_value)            │
                    │    85                                                    │
                    │    86                                                    │
                    │    87 def _parse_file(path: Path) -> str:                │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /seq_basics/_plumbing/resolve.py:152 in _clean_sequence  │
                    │                                                          │
                    │   149 │   # Validate characters                          │
                    │   150 │   invalid = set(cleaned) - VALID_SEQUENCE_CHARS  │
                    │   151 │   if invalid:                                    │
                    │ ❱ 152 │   │   raise ValueError(                          │
                    │   153 │   │   │   │   │   │   f"Invalid sequence charact │
                    │   154 │   │   │   │   │   │   "Valid DNA/RNA/IUPAC chara │
                    │   155 )                                                  │
                    ╰──────────────────────────────────────────────────────────╯
                    ValueError: Invalid sequence characters detected: ['.',     
                    '_']. Valid DNA/RNA/IUPAC characters are: ATUCGRSYKWMN.     
[Tool result] ← crispr_predict_offtargets:
{
  "error": "Error calling tool 'crispr_predict_offtargets': Invalid sequence characters detected: ['.', '_']. Valid DNA/RNA/IUPAC ch\u2026 (147 chars)"
}

[Tool call] → crispr_predict_offtargets
{
  "reference": "EMX1",
  "protospacer": "GCTCTCCAAACCACGGCTCC"
}
[05/08/26 17:28:22] Error calling tool 'crispr_predict_offtargets'              
                    ╭─────────── Traceback (most recent call last) ────────────╮
                    │ /opt/anaconda3/lib/python3.12/site-packages/fastmcp/serv │
                    │ er/server.py:987 in call_tool                            │
                    │                                                          │
                    │ /opt/anaconda3/lib/python3.12/site-packages/fastmcp/tool │
                    │ s/tool.py:354 in _run                                    │
                    │                                                          │
                    │                 ... 7 frames hidden ...                  │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /seq_basics/_plumbing/resolve.py:84 in resolve_to_seq    │
                    │                                                          │
                    │    81 │   │   │   return _clean_sequence(remainder)      │
                    │    82 │                                                  │
                    │    83 │   # 4. Assume it's a raw sequence string - clean │
                    │ ❱  84 │   return _clean_sequence(input_value)            │
                    │    85                                                    │
                    │    86                                                    │
                    │    87 def _parse_file(path: Path) -> str:                │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /seq_basics/_plumbing/resolve.py:152 in _clean_sequence  │
                    │                                                          │
                    │   149 │   # Validate characters                          │
                    │   150 │   invalid = set(cleaned) - VALID_SEQUENCE_CHARS  │
                    │   151 │   if invalid:                                    │
                    │ ❱ 152 │   │   raise ValueError(                          │
                    │   153 │   │   │   │   │   │   f"Invalid sequence charact │
                    │   154 │   │   │   │   │   │   "Valid DNA/RNA/IUPAC chara │
                    │   155 )                                                  │
                    ╰──────────────────────────────────────────────────────────╯
                    ValueError: Invalid sequence characters detected: ['E',     
                    'X']. Valid DNA/RNA/IUPAC characters are: ATUCGRSYKWMN.     
[Tool result] ← crispr_predict_offtargets:
{
  "error": "Error calling tool 'crispr_predict_offtargets': Invalid sequence characters detected: ['E', 'X']. Valid DNA/RNA/IUPAC ch\u2026 (147 chars)"
}

[Tool call] → crispr_predict_offtargets
{
  "reference": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
  "protospacer": "GCTCTCCAAACCACGGCTCC"
}
[05/08/26 17:31:08] Error calling tool 'crispr_predict_offtargets'              
                    ╭─────────── Traceback (most recent call last) ────────────╮
                    │ /opt/anaconda3/lib/python3.12/site-packages/fastmcp/serv │
                    │ er/server.py:987 in call_tool                            │
                    │                                                          │
                    │ /opt/anaconda3/lib/python3.12/site-packages/fastmcp/tool │
                    │ s/tool.py:354 in _run                                    │
                    │                                                          │
                    │                 ... 6 frames hidden ...                  │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /seq_basics/_plumbing/register.py:354 in wrapped         │
                    │                                                          │
                    │   351 │   │   # Drop any kwargs the underlying function  │
                    │   352 │   │   if _func_params is not None:               │
                    │   353 │   │   │   kwargs = {k: v for k, v in kwargs.item │
                    │ ❱ 354 │   │   return func(**kwargs)                      │
                    │   355 │                                                  │
                    │   356 │   wrapped.__name__ = mcp_name                    │
                    │   357 │   wrapped.__doc__  = description                 │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /crispr_tools/tools/predict_offtargets.py:227 in run     │
                    │                                                          │
                    │   224 │   │                                              │
                    │   225 │   │   invalid_ref = [b for b in set(reference) i │
                    │   226 │   │   if invalid_ref:                            │
                    │ ❱ 227 │   │   │   raise ValueError(                      │
                    │   228 │   │   │   │   f"Invalid base(s) in reference: {s │
                    │   229 │   │   │   │   "Only standard bases A, T, G, C (a │
                    │   230 │   │   │   )                                      │
                    ╰──────────────────────────────────────────────────────────╯
                    ValueError: Invalid base(s) in reference: ['U']. Only       
                    standard bases A, T, G, C (and N) are accepted.             
[Tool result] ← crispr_predict_offtargets:
{
  "error": "Error calling tool 'crispr_predict_offtargets': Invalid base(s) in reference: ['U']. Only standard bases A, T, G, C (and\u2026 (137 chars)"
}

[Tool call] → crispr_predict_offtargets
{
  "reference": "GCTCTCCAAACCACGGCTCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
  "protospacer": "GCTCTCCAAACCACGGCTCC"
}
[Tool result] ← crispr_predict_offtargets:
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "nuclease": "cas9",
  "reference_length": 96,
  "strands_scanned": 2,
  "sites_evaluated": 154,
  "offtarget_sites": [
    {
      "position": 0,
      "strand": "+",
      "sequence": "GCTCTCCAAACCACGGCTCC",
      "mismatches": 0,
      "mismatch_positions": [],
      "seed_mismatches": 0,
      "has_pam": false,
      "risk": "LOW",
      "cfd_score": 0.0
    }
  ],
  "high_risk_count": 0,
  "aggregate_offtarget_cfd": 0.0,
  "max_offtarget_cfd": 0.0,
  "specificity_summary": "1 potential off-target site(s) found within 3 mismatches (0 HIGH risk). Consider redesigning the guide or choosing a mor\u2026 (143 chars)"
}

Gemini: The CRISPR design for targeting the **EMX1** gene in *Homo sapiens* using the **pX330** vector is complete.

### Guide RNA Selection
*   **Best Guide:** `GCTCTCCAAACCACGGCTCC`
*   **Scoring Rationale:** This guide was selected because it achieved a total score of 3/3 (Efficiency 2/2, Specificity 1/1). It has 65% GC content and 0 predicted off-target sites (0 HIGH or MEDIUM risk).

### Pre-experiment Efficiency Prediction
*   **Predicted NHEJ Knockout Efficiency:** 40.2% (±15% confidence interval: 25.2–55.2%)
*   **Interpretation:** Moderate efficiency. This guide is usable, but you should expect to screen a slightly larger number of colonies to recover an edited clone.

### Off-target Analysis
*   **Specificity Flag:** The max predicted off-target CFD score is 0.0, indicating high specificity for the target site.

### Cloning Oligos
To clone this protospacer into the pX330 vector (using the BbsI restriction site), order the following annealed oligos:
*   **Top:** `CACCGCTCTCCAAACCACGGCTCC`
*   **Bottom:** `AAACGGAGCCGTGGTTTGGAGAGC`

***

All CRISPR design and validation steps are complete. Would you like me to generate:
(a) **a construction file** — a structured record of the cloning workflow,
(b) **a lab sheet** — bench-ready step-by-step protocol (requires a construction file),
(c) **both**, or
(d) **neither**?

You: 
```

</details>

> **What happened?**
> - `predict_editing_efficiency` was called autonomously by Gemini as an extra validation step — it is **not** invoked by `run_full_crispr_workflow.py` directly and is not part of the core pipeline.
> - The three `crispr_predict_offtargets` errors occurred because Gemini initially passed `"NC_000002.12"` (an NCBI accession) and `"EMX1"` (a gene name) as the `reference` parameter. `predict_offtargets` requires a raw DNA sequence; the `predict_offtargets.json` tool description was updated after this session to explicitly forbid accessions/gene names and direct callers to use `crispr_fetch_target_sequence` first.

> **Benchmark alignment — Hsu et al. 2013:** The sample above targets EMX1 in *Homo sapiens* with pX330, one of the canonical SpCas9 benchmark loci from [Hsu, P.D. et al. 2013, *Nat Biotechnol*](https://doi.org/10.1038/nbt.2647). Alignment points: NGG PAM targeting (SpCas9), mammalian guide generation via U6-driven sgRNA scaffold, pX330 BbsI CACC/AAAC sticky-end overhangs, and zero predicted off-target sites (low-off-target prioritization).

Possible Unhappy Paths:
- Predict editing efficiency and predict offtargets could be called in pipeline before correct gene is used
- A user may say "Run full CRISPR workflow on genes that cause cystic fibrosis in Homo sapiens." If only a disease name is given, semantic tools are called to identify a target gene before the workflow proceeds.
- After gene resolution, the organism is occasionally re-prompted if it was inferred rather than explicitly stated.
- Vectors are not always known; the MCP supplies organism-appropriate recommendations or waits for user input before calling `design_cloning_oligos`.
- Construction file regions that cannot be determined are printed as incomplete or `"N"` — users should review those fields before submitting to wet lab.
 
- MCP Wrapper: [run_full_crispr_workflow.json](modules/crispr_tools/tools/run_full_crispr_workflow.json)
- Pytests ([tests/unit/test_run_full_crispr_workflow.py](tests/unit/test_run_full_crispr_workflow.py)): guide selection, vector prompting, workflow confirmation gate, upstream-selected gene confirmation, empty query error, per-guide score fields
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

> **What happened ?**
> - `fetch_target_sequence.py` normalized the organism alias ("c.elegans" → "Caenorhabditis elegans") before querying NCBI.
> - The tool then used Entrez `esearch`, `esummary`, and `efetch` operations to retrieve genomic loci for sma-2, sma-3, and sma-4.
> - The returned accession and genomic coordinate information were formatted into downstream-compatible sequence records for CRISPR guide discovery.

> **Benchmark alignment — NCBI Entrez E-utilities:**  
> This workflow follows standard NCBI Entrez gene-resolution practices by converting canonical gene symbols into genomic loci suitable for downstream genome engineering workflows.

- Citations:
  - National Center for Biotechnology Information. n.d. "Entrez Programming Utilities Help (E-utilities)." U.S National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK25501/. *(basis for esearch/esummary/elink/efetch workflow used for automated sequence retrieval and accession resolution)*

  - National Center for Biotechnology Information. n.d. "Gene Database Help." U.S. National Library of Medicine. https://www.ncbi.nlm.nih.gov/books/NBK3841/. *(gene locus and accession-resolution framework used before downstream CRISPR guide discovery)*
  
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

> **What happened ?**
> - `crispr_cas_selector` evaluated both NGG-compatible SpCas9 windows and TTTV-compatible LbCas12a windows across all supplied loci.
> - The workflow aggregated guide statistics across sma-2, sma-3, and sma-4 before generating a single nuclease recommendation.
> - Multiplexing support was treated as a tie-breaker when guide counts between nucleases fell within the configured threshold margin.


> **Benchmark alignment — Zetsche et al. 2015 / 2017:**  
> This workflow mirrors canonical Cas12a multiplexing logic by comparing TTTV PAM availability against NGG PAM availability across multiple AT-rich loci in *Caenorhabditis elegans*. Alignment points include LbCas12a PAM recognition, crRNA-array compatibility, and multiplex editing rationale.
  
- Citations:
  - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ian M. Slaymaker, Kira S. Makarova, Patrick Essletzbichler, et al. 2015. “Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System.” Cell 163 (3): 759–771. https://doi.org/10.1016/j.cell.2015.09.038. *(LbCas12a TTTV PAM recognition, crRNA architecture, and Cas12a targeting rationale used by the selector)*

  - Kim, Daesik, Jinu Kim, Jeongbin Hur, Kwang-Hyun Been, Sung Roh Yoon, and Jin-Soo Kim. 2016. “Genome-wide Analysis Reveals Specificities of Cpf1 Endonucleases in Human Cells.” Nature Biotechnology 34 (8): 863–868. https://doi.org/10.1038/nbt.3609. *(Cas12a specificity behavior and PAM-density heuristic supporting Cas9-vs-Cas12a recommendation logic)*

  - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ian M. Slaymaker, Maximilian V. Makarova, Patrick Essletzbichler, et al. 2017. “Multiplex Gene Editing by CRISPR-Cpf1 Using a Single crRNA Array.” Nature Biotechnology 35 (1): 31–34. https://doi.org/10.1038/nbt.3737. *(multiplex Cas12a crRNA-array framework and multi-locus editing rationale used for aggregated recommendations)*

  - Mali, Prashant, Luhan Yang, Kevin M. Esvelt, John Aach, Marc Guell, James E. DiCarlo, Julie E. Norville, and George M. Church. 2013. “RNA-Guided Human Genome Engineering via Cas9.” Science 339 (6121): 823–826. https://doi.org/10.1126/science.1232033. *(U6 promoter-driven sgRNA expression constraints supporting poly-T filtering and practical guide-design heuristics)*

  - Ghorbal, Mehdi, Michelle Gorman, Christopher R. Macpherson, Rui Martins, Anja Scherf, and Jake Baum. 2014. “Genome Editing in the Human Malaria Parasite Plasmodium falciparum Using the CRISPR-Cas9 System.” Nature Methods 11 (8): 819–821. https://doi.org/10.1038/nmeth.3063. *(extreme AT-rich genome editing context motivating Cas12a-favoring AT-richness heuristics)*

- Minkah, Nkuli K., et al. 2018. “CRISPR-Cas9-Mediated Genome Editing of Plasmodium falciparum.” Cell Host & Microbe. *(AT-rich genome editing validation supporting TTTV PAM abundance rationale in highly AT-biased loci)*
  
  
  
## 4. [design_cas9_grna](modules/crispr_tools/tools/design_cas9_grna.py) with [rank_guides.py](modules/crispr_tools/tools/rank_guides.py)
  - What it does: 
    - Enumerates SpCas9 guide RNA candidates from NGG PAM sites across the sequence; candidate scoring applies GC content (40–65% optimum) and poly-T checks as heuristics derived from Doench et al. 2016 Rule Set 2 — downstream [rank_guides.py](modules/crispr_tools/tools/rank_guides.py) applies weighted heuristics for final ranking (not a reproduction of the trained Doench ML model)
     - Shared downstream scoring layer: candidate guides generated here are passed into [rank_guides.py](modules/crispr_tools/tools/rank_guides.py), where efficiency heuristics (GC content, poly-T filtering, homopolymer penalties) and specificity metrics from `predict_offtargets.py` are combined into a final ranked score used by downstream cloning workflows.
    - Infrastructure integration: although [rank_guides.py](modules/crispr_tools/tools/rank_guides.py) is not exposed as a standalone end-user MCP workflow, it includes its own MCP JSON wrapper [rank_guides.json](modules/crispr_tools/tools/rank_guides.json) and dedicated pytest suite (`tests/unit/test_rank_guides.py`) for validation of guide ranking, scoring consistency, off-target aggregation, and best-guide selection behavior.
    - Downstream [rank_guides.py](modules/crispr_tools/tools/rank_guides.py) applies weighted heuristics for final ranking (not a reproduction of the trained Doench ML model)
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

You: Design guide RNA for EXM1 in Homo Sapiens  

[Tool call] → crispr_fetch_target_sequence
{
  "query": "EXM1",
  "organism": "Homo sapiens"
}
[05/08/26 18:05:11] Error calling tool 'crispr_fetch_target_sequence'           
                    ╭─────────── Traceback (most recent call last) ────────────╮
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/.venv/l │
                    │ ib/python3.12/site-packages/fastmcp/tools/tool_manager.p │
                    │ y:160 in call_tool                                       │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/.venv/l │
                    │ ib/python3.12/site-packages/fastmcp/tools/tool.py:381 in │
                    │ run                                                      │
                    │                                                          │
                    │                 ... 3 frames hidden ...                  │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /crispr_tools/tools/fetch_target_sequence.py:241 in      │
                    │ _fetch_locus_from_ncbi                                   │
                    │                                                          │
                    │   238                                                    │
                    │   239 def _fetch_locus_from_ncbi(gene_name: str, organis │
                    │   240 │   """Return (sequence, gene_id, accession) for t │
                    │ ❱ 241 │   gene_id = _esearch_gene_id(gene_name, organism │
                    │   242 │                                                  │
                    │   243 │   summary_resp = _ncbi_get(                      │
                    │   244 │   │   f"{_ENTREZ_BASE}/esummary.fcgi",           │
                    │                                                          │
                    │ /Users/emory/Documents/BioE_134-234_FinalProject/modules │
                    │ /crispr_tools/tools/fetch_target_sequence.py:229 in      │
                    │ _esearch_gene_id                                         │
                    │                                                          │
                    │   226 │   │   ids = resp.json()["esearchresult"]["idlist │
                    │   227 │   │   if ids:                                    │
                    │   228 │   │   │   return ids[0]                          │
                    │ ❱ 229 │   raise ValueError(                              │
                    │   230 │   │   f"Gene '{gene_name}' not found in NCBI for │
                    │   231 │   │   "Check the spelling or try a more specific │
                    │   232 │   )                                              │
                    ╰──────────────────────────────────────────────────────────╯
                    ValueError: Gene 'EXM1' not found in NCBI for organism 'Homo
                    sapiens'. Check the spelling or try a more specific organism
                    name.                                                       
[Tool result] ← crispr_fetch_target_sequence:
{
  "error": "Error calling tool 'crispr_fetch_target_sequence': Gene 'EXM1' not found in NCBI for organism 'Homo sapiens'. Check the \u2026 (166 chars)"
}

[Tool call] → semantic_gene_search
{
  "query": "EXM1",
  "organism": "Homo sapiens"
}
[Tool result] ← semantic_gene_search:
{
  "parsed_query": {
    "raw_query": "EXM1",
    "organism": null,
    "keywords": [
      "exm"
    ],
    "ontology_terms": [
      "exm"
    ]
  },
  "go_terms": []
}

[Tool call] → semantic_gene_search
{
  "organism": "Homo sapiens",
  "query": "EMX1"
}
[Tool result] ← semantic_gene_search:
{
  "parsed_query": {
    "raw_query": "EMX1",
    "organism": null,
    "keywords": [
      "emx"
    ],
    "ontology_terms": [
      "emx"
    ]
  },
  "go_terms": []
}

[Tool call] → crispr_fetch_target_sequence
{
  "organism": "Homo sapiens",
  "query": "EMX1"
}
[Tool result] ← crispr_fetch_target_sequence:
{
  "sequence": "GAGCGCTCTCCAAACCACGGCTCCCGGCGCTCAGGCGGGCCGCTGCCAAGACCCGGCCTGGAGTCCCCGCAGAGTTGCGCGGCGCACGGACCCCGTGGCCTTGGGGCGTCAGGAGGCCCA\u2026 (18703 chars)",
  "source": "ncbi",
  "resource": "EMX1",
  "organism": "Homo sapiens",
  "target_type": "genomic_locus",
  "ncbi_gene_id": "2016",
  "ncbi_accession": "NC_000002.12",
  "length": 18703,
  "note": "Fetched genomic locus from NCBI Entrez: gene 'EMX1' in 'Homo sapiens'. Nucleotide record: NC_000002.12. Verify the recor\u2026 (162 chars)"
}

[Tool call] → crispr_design_cas9_grna
{
  "seq": "GAGCGCTCTCCAAACCACGGCTCCCGGCGCTCAGGCGGGCCGCTGCCAAGACCCGGCCTGGAGTCCCCGCAGAGTTGCGCGGCGCACGGACCCCGTGGCCTTGGGGCGTCAGGAGGCCCA\u2026 (18703 chars)"
}
[Tool result] ← crispr_design_cas9_grna:
{
  "result": [
    {
      "grna_sequence": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GCTCTCCAAACCACGGCTCC",
      "pam_site": "CGG"
    },
    {
      "grna_sequence": "UUUAGAAGUUACUGCCGGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TTTAGAAGTTACTGCCGGGA",
      "pam_site": "AGG"
    },
    {
      "grna_sequence": "GUCCAGGACUCCGAAGGUCGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTCCAGGACTCCGAAGGTCG",
      "pam_site": "GGG"
    },
    {
      "grna_sequence": "UACAUAGAUGUUUCCAGAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TACATAGATGTTTCCAGAGA",
      "pam_site": "TGG"
    },
    {
      "grna_sequence": "ACUACGUGGUGGGCGCCGAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "ACTACGTGGTGGGCGCCGAG",
      "pam_site": "CGG"
    },
    {
      "grna_sequence": "CACUUAACCUCUCUGAGCUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CACTTAACCTCTCTGAGCTC",
      "pam_site": "TGG"
    },
    {
      "grna_sequence": "GGAAAGAAGGGCAGGGUUUGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAAAGAAGGGCAGGGTTTG",
      "pam_site": "AGG"
    },
    {
      "grna_sequence": "GGAGAAUGGAGACUUUUGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAGAATGGAGACTTTTGGA",
      "pam_site": "AGG"
    },
    {
      "grna_sequence": "GUGUAGGGGCAGCACCUCCUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTGTAGGGGCAGCACCTCCT",
      "pam_site": "GGG"
    },
    {
      "grna_sequence": "CAGCUGGUCAGAGGGGACCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CAGCTGGTCAGAGGGGACCC",
      "pam_site": "CGG"
    }
  ]
}

[Tool call] → crispr_rank_guides
{
  "nuclease": "cas9",
  "reference": "GAGCGCTCTCCAAACCACGGCTCCCGGCGCTCAGGCGGGCCGCTGCCAAGACCCGGCCTGGAGTCCCCGCAGAGTTGCGCGGCGCACGGACCCCGTGGCCTTGGGGCGTCAGGAGGCCCA\u2026 (18703 chars)",
  "guides": [
    {
      "pam_site": "CGG",
      "grna_sequence": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GCTCTCCAAACCACGGCTCC"
    },
    {
      "pam_site": "AGG",
      "grna_sequence": "UUUAGAAGUUACUGCCGGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TTTAGAAGTTACTGCCGGGA"
    },
    {
      "grna_sequence": "GUCCAGGACUCCGAAGGUCGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTCCAGGACTCCGAAGGTCG",
      "pam_site": "GGG"
    },
    {
      "grna_sequence": "UACAUAGAUGUUUCCAGAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TACATAGATGTTTCCAGAGA",
      "pam_site": "TGG"
    },
    {
      "grna_sequence": "ACUACGUGGUGGGCGCCGAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "ACTACGTGGTGGGCGCCGAG",
      "pam_site": "CGG"
    },
    {
      "grna_sequence": "CACUUAACCUCUCUGAGCUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CACTTAACCTCTCTGAGCTC",
      "pam_site": "TGG"
    },
    {
      "grna_sequence": "GGAAAGAAGGGCAGGGUUUGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAAAGAAGGGCAGGGTTTG",
      "pam_site": "AGG"
    },
    {
      "pam_site": "AGG",
      "grna_sequence": "GGAGAAUGGAGACUUUUGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAGAATGGAGACTTTTGGA"
    },
    {
      "pam_site": "GGG",
      "grna_sequence": "GUGUAGGGGCAGCACCUCCUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTGTAGGGGCAGCACCTCCT"
    },
    {
      "pam_site": "CGG",
      "grna_sequence": "CAGCUGGUCAGAGGGGACCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CAGCTGGTCAGAGGGGACCC"
    }
  ]
}
[Tool result] ← crispr_rank_guides:
{
  "ranked_guides": [
    {
      "pam_site": "CGG",
      "grna_sequence": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GCTCTCCAAACCACGGCTCC",
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
      "pam_site": "AGG",
      "grna_sequence": "UUUAGAAGUUACUGCCGGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TTTAGAAGTTACTGCCGGGA",
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
      "grna_sequence": "GUCCAGGACUCCGAAGGUCGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTCCAGGACTCCGAAGGTCG",
      "pam_site": "GGG",
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
      "grna_sequence": "ACUACGUGGUGGGCGCCGAGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "ACTACGTGGTGGGCGCCGAG",
      "pam_site": "CGG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.7,
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
      "grna_sequence": "CACUUAACCUCUCUGAGCUCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CACTTAACCTCTCTGAGCTC",
      "pam_site": "TGG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.5,
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
      "grna_sequence": "GGAAAGAAGGGCAGGGUUUGGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAAAGAAGGGCAGGGTTTG",
      "pam_site": "AGG",
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
      "pam_site": "GGG",
      "grna_sequence": "GUGUAGGGGCAGCACCUCCUGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GTGTAGGGGCAGCACCTCCT",
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
      "pam_site": "CGG",
      "grna_sequence": "CAGCUGGUCAGAGGGGACCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "CAGCTGGTCAGAGGGGACCC",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.7,
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
      "grna_sequence": "UACAUAGAUGUUUCCAGAGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "TACATAGATGTTTCCAGAGA",
      "pam_site": "TGG",
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
      "pam_site": "AGG",
      "grna_sequence": "GGAGAAUGGAGACUUUUGGAGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
      "protospacer": "GGAGAATGGAGACTTTTGGA",
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
    "pam_site": "CGG",
    "grna_sequence": "GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC",
    "protospacer": "GCTCTCCAAACCACGGCTCC",
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
  "scoring_rationale": "Best guide 'GCTCTCCAAACCACGGCTCC' scored 3/3 (efficiency 2/2, specificity 1/1). GC content: 65%. Off-target sites: 0 tot\u2026 (142 chars)"
}

Gemini: The CRISPR guide RNA design for targeting the **EMX1** gene in *Homo sapiens* has been completed.

### Best Guide RNA Candidate
The top-ranked candidate was selected based on a combined score of on-target efficiency (GC content, absence of poly-T runs, and PAM-proximal G) and off-target safety (using the CFD scoring system to ensure zero predicted off-target sites).

*   **Protospacer:** `GCTCTCCAAACCACGGCTCC`
*   **PAM:** `CGG`
*   **gRNA Sequence:** `GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC`
*   **Rationale:** This guide achieved a perfect score of 3/3, demonstrating optimal GC content (65%) and no predicted off-target binding sites.

---

### Next Steps
Would you like me to proceed with one of the following?

1.  **Generate a construction file** — to create a structured cloning workflow record.
2.  **Generate a lab sheet** — for a bench-ready, step-by-step experimental protocol (requires a construction file first).
3.  **Run off-target validation** — for a detailed CFD analysis.
4.  **Predict editing efficiency** — for a pre-experiment performance estimate.

Please let me know how you would like to proceed.

You: 

```
</details>

> **What happened?**
> EMX1 was not recognizable in NCBI, therefore semantic tools were used instead. This will be diagnosed with future work. 

> **Benchmark alignment — Hsu et al. 2013 / Ran et al. 2013:**  
> The sample above targets EMX1 in *Homo sapiens*, one of the canonical SpCas9 benchmark loci used in early mammalian CRISPR-Cas9 specificity studies. Alignment points include NGG PAM scanning (SpCas9), 20 nt protospacer extraction upstream of the PAM, and generation of a U6-compatible sgRNA scaffold architecture consistent with pX330 mammalian guide-expression systems.
  
  - Citations: 
    - Jinek, Martin, Krzysztof Chylinski, Ines Fonfara, Michael Hauer, Jennifer A. Doudna, and Emmanuelle Charpentier. 2012. "A Programmable Dual-RNA–Guided DNA Endonuclease in Adaptive Bacterial Immunity." Science 337 (6096): 816–821. https://doi.org/10.1126/science.1225829.
    *(programmable Cas9 and sgRNA architecture)*
  
    - Hsu, Patrick D., David A. Scott, Joshua A. Weinstein, F. Ann Ran, Silvana Konermann, Vineeta Agarwala, Yinqing Li, et al. 2013. "DNA Targeting Specificity of RNA-Guided Cas9 Nucleases." Nature Biotechnology 31 (9): 827–832. https://doi.org/10.1038/nbt.2647.
    *(EMX1 benchmark loci and SpCas9 NGG specificity framework)*

    - Ran, F. Ann, Patrick D. Hsu, Jason Wright, Vineeta Agarwala, David A. Scott, and Feng Zhang. 2013. “Genome Engineering Using the CRISPR-Cas9 System.” Nature Protocols 8 (11): 2281–2308. https://doi.org/10.1038/nprot.2013.143. *(pX330 BbsI cloning protocol; basis for mammalian guide-insertion workflow)*
  
    - Doench, John G., Nicolo Fusi, Meagan Sullender, Mudra Hegde, Emma W. Vaimberg, Katherine F. Donovan, Ian Smith, et al. 2016. "Optimized sgRNA Design to Maximize Activity and Minimize Off-Target Effects of CRISPR-Cas9." Nature Biotechnology 34 (2): 184–191. https://doi.org/10.1038/nbt.3437. *(GC content optimum 40–65%, poly-T termination penalty, PAM-proximal position weights used in candidate scoring)*
      
    - Cong, Le, F. Ann Ran, David Cox, Shuailiang Lin, Robert Barretto, Naomi Habib, Patrick D. Hsu, et al. 2013. “Multiplex Genome Engineering Using CRISPR/Cas Systems.” Science 339 (6121): 819–823. https://doi.org/10.1126/science.1231143. *(mammalian SpCas9 guide design using 20 nt target sequences adjacent to NGG PAMs)*
  
## 5. [design_cas12a_crrna.py](modules/crispr_tools/tools/design_cas12a_crrna.py) with [rank_guides.py](modules/crispr_tools/tools/rank_guides.py)
  - What it does: 
    - Scans a sequence for TTTV PAM sites and returns up to 10 LbCas12a crRNA candidates
    - Shared downstream scoring layer: candidate guides generated here are passed into [rank_guides.py](modules/crispr_tools/tools/rank_guides.py), where efficiency heuristics (GC content, poly-T filtering, homopolymer penalties) and specificity metrics from `predict_offtargets.py` are combined into a final ranked score used by downstream cloning workflows.
    - Infrastructure integration: although [rank_guides.py](modules/crispr_tools/tools/rank_guides.py) is not exposed as a standalone end-user MCP workflow, it includes its own MCP JSON wrapper [rank_guides.json](modules/crispr_tools/tools/rank_guides.json) and dedicated pytest suite (`tests/unit/test_rank_guides.py`) for validation of guide ranking, scoring consistency, off-target aggregation, and best-guide selection behavior.
    - Downstream [rank_guides.py](modules/crispr_tools/tools/rank_guides.py) applies weighted heuristics for final ranking (not a reproduction of the trained Doench ML model)
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
     You: Design Cas12a crRNAs for sma-2 in C.elegans and rank 

[Tool call] → crispr_design_cas12a_crrna
{
  "organism": "C. elegans",
  "seq": "GGACACAAGTTGACTGTCGCCGGGTGGAAACACATCTGTGAAGATGATTGCACAATTTGGAGAGACGGAGAGTTGAAAACAATGAGTAGACAGTTAGACGGTTGATGATGGGCTGGTGAA\u2026 (7410 chars)"
}
[Tool result] ← crispr_design_cas12a_crrna:
{
  "result": [
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUGAGAGACGGAGAGUUGAAAACAA",
      "protospacer": "GAGAGACGGAGAGTTGAAAACAA",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUCCAAAAGACCAAUAGGGUGGGG",
      "protospacer": "TCCAAAAGACCAATAGGGTGGGG",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUACAAAAAUUACACCAGCCAAAAC",
      "protospacer": "ACAAAAATTACACCAGCCAAAAC",
      "pam_site": "TTTC"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUGCACGUGGAGAUCACUGGAAAA",
      "protospacer": "TGCACGTGGAGATCACTGGAAAA",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUUGCUCAUGAUUCGUCACUGCUG",
      "protospacer": "TTGCTCATGATTCGTCACTGCTG",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUAAAAGAUCAAGAAUAUCCCAAGU",
      "protospacer": "AAAAGATCAAGAATATCCCAAGT",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUGUAGAGCAUUGCGACAAUUUGAU",
      "protospacer": "GTAGAGCATTGCGACAATTTGAT",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUAGAGGGAAAAGGUUACGAAAAGA",
      "protospacer": "AGAGGGAAAAGGTTACGAAAAGA",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUAAAAAAACUGGGUUACUGUACUU",
      "protospacer": "AAAAAAACTGGGTTACTGTACTT",
      "pam_site": "TTTG"
    },
    {
      "crrna_sequence": "AAUUUCUACUAAGUGUAGAUUCGGUCGUCAAAAAAUGUUCAAA",
      "protospacer": "TCGGTCGTCAAAAAATGTTCAAA",
      "pam_site": "TTTG"
    }
  ]
}

[Tool call] → crispr_rank_guides
{
  "reference": "GGACACAAGTTGACTGTCGCCGGGTGGAAACACATCTGTGAAGATGATTGCACAATTTGGAGAGACGGAGAGTTGAAAACAATGAGTAGACAGTTAGACGGTTGATGATGGGCTGGTGAA\u2026 (7410 chars)",
  "nuclease": "cas12a",
  "guides": [
    "GAGAGACGGAGAGTTGAAAACAA",
    "TCCAAAAGACCAATAGGGTGGGG",
    "ACAAAAATTACACCAGCCAAAAC",
    "TGCACGTGGAGATCACTGGAAAA",
    "TTGCTCATGATTCGTCACTGCTG",
    "AAAAGATCAAGAATATCCCAAGT",
    "GTAGAGCATTGCGACAATTTGAT",
    "AGAGGGAAAAGGTTACGAAAAGA",
    "AAAAAAACTGGGTTACTGTACTT",
    "TCGGTCGTCAAAAAATGTTCAAA"
  ]
}
[Tool result] ← crispr_rank_guides:
{
  "ranked_guides": [
    {
      "protospacer": "GAGAGACGGAGAGTTGAAAACAA",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.435,
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
      "protospacer": "TCCAAAAGACCAATAGGGTGGGG",
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
        "on_target_found": true
      },
      "total_score": 3
    },
    {
      "protospacer": "TGCACGTGGAGATCACTGGAAAA",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.478,
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
      "protospacer": "TTGCTCATGATTCGTCACTGCTG",
      "efficiency_score": 2,
      "efficiency_details": {
        "gc_content": 0.478,
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
      "protospacer": "ACAAAAATTACACCAGCCAAAAC",
      "efficiency_score": 1,
      "efficiency_details": {
        "gc_content": 0.348,
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
      "protospacer": "AAAAGATCAAGAATATCCCAAGT",
      "efficiency_score": 1,
      "efficiency_details": {
        "gc_content": 0.304,
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
      "protospacer": "GTAGAGCATTGCGACAATTTGAT",
      "efficiency_score": 1,
      "efficiency_details": {
        "gc_content": 0.391,
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
      "protospacer": "AGAGGGAAAAGGTTACGAAAAGA",
      "efficiency_score": 1,
      "efficiency_details": {
        "gc_content": 0.391,
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
      "protospacer": "AAAAAAACTGGGTTACTGTACTT",
      "efficiency_score": 1,
      "efficiency_details": {
        "gc_content": 0.304,
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
      "protospacer": "TCGGTCGTCAAAAAATGTTCAAA",
      "efficiency_score": 1,
      "efficiency_details": {
        "gc_content": 0.348,
        "gc_content_ok": false,
        "no_polyt_run": true
      },
      "specificity_score": 1,
      "specificity_details": {
        "high_risk_offtargets": 0,
        "medium_risk_offtargets": 1,
        "total_offtargets": 1,
        "on_target_found": true
      },
      "total_score": 2
    }
  ],
  "best_guide": {
    "protospacer": "GAGAGACGGAGAGTTGAAAACAA",
    "efficiency_score": 2,
    "efficiency_details": {
      "gc_content": 0.435,
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
  "scoring_rationale": "Best guide 'GAGAGACGGAGAGTTGAAAACAA' scored 3/3 (efficiency 2/2, specificity 1/1). GC content: 44%. Off-target sites: 0 \u2026 (145 chars)"
}

Gemini: The Cas12a crRNA design for the *sma-2* gene in *C. elegans* is complete. The top-ranked crRNA was selected based on a combination of high on-target efficiency scores (GC content optimization and lack of poly-T runs) and specificity (zero predicted off-target binding sites).

### Top-Ranked crRNA
*   **Protospacer:** `GAGAGACGGAGAGTTGAAAACAA`
*   **PAM Site:** `TTTG` (5' to the protospacer)
*   **Total Score:** 3/3
*   **Rationale:** This guide has a GC content of 44%, which is within the optimal range, and contains no poly-T sequence runs that could cause premature termination of Pol III-driven transcription. Most importantly, it shows no predicted off-target binding sites in the *C. elegans* genome.

### Summary of Ranked Candidates
| Protospacer | GC% | Total Score | Specificity Notes |
| :--- | :--- | :--- | :--- |
| **GAGAGACGGAGAGTTGAAAACAA** | 44% | 3 | 0 total off-targets |
| **TCCAAAAGACCAATAGGGTGGGG** | 52% | 3 | 0 total off-targets |
| **TGCACGTGGAGATCACTGGAAAA** | 48% | 3 | 0 total off-targets |
| **TTGCTCATGATTCGTCACTGCTG** | 48% | 3 | 0 total off-targets |
| **ACAAAAATTACACCAGCCAAAAC** | 35% | 2 | Low GC% |
| **AAAAGATCAAGAATATCCCAAGT** | 30% | 2 | Low GC% |
| **GTAGAGCATTGCGACAATTTGAT** | 39% | 2 | Low GC% |
| **AGAGGGAAAAGGTTACGAAAAGA** | 39% | 2 | Low GC% |
| **AAAAAAACTGGGTTACTGTACTT** | 30% | 2 | Low GC% |
| **TCGGTCGTCAAAAAATGTTCAAA** | 35% | 2 | 1 medium risk off-target |

All designs were scanned for TTTV PAM sites.

**Next Steps:**
1.  Would you like me to predict the on-target editing efficiency for this top guide?
2.  Would you like me to design the cloning oligos needed to clone this crRNA into an expression vector?
3.  Would you like me to design flanking Sanger sequencing primers to verify the edit later?

You: 

```
</details>

> **What Happened?**
> - `design_cas12a_crrna.py` scanned the fetched genomic locus for valid TTTV PAM sites compatible with LbCas12a.
> - The tool generated up to 10 candidate crRNAs by extracting the 23 nt protospacer immediately downstream of each TTTV PAM and prepending the Cas12a direct repeat sequence.
> - These candidate crRNAs were then passed into `rank_guides.py`, where GC-content heuristics, poly-T filtering, homopolymer penalties, and off-target metrics from `predict_offtargets.py` were aggregated into a final ranked score.
> - The highest-scoring crRNA was selected as the recommended downstream cloning and editing candidate.

> **Benchmark alignment — Zetsche et al. 2015 / Kim et al. 2016:**
> The workflow follows canonical LbCas12a targeting logic by scanning for TTTV PAMs and generating staggered-cut-compatible crRNAs. The ranking stage prioritizes low predicted off-target activity and biologically reasonable spacer composition, consistent with early Cas12a specificity and multiplexing studies.
    
  - Citations:
  - Zetsche, Bernd, Jonathan S. Gootenberg, Omar O. Abudayyeh, Ioanna M. Slaymaker, Kira S. Makarova, Patrick Essletzbichler, Sara E. Volz, et al. 2015. “Cpf1 Is a Single RNA-Guided Endonuclease of a Class 2 CRISPR-Cas System.” Cell 163 (3): 759–771. https://doi.org/10.1016/j.cell.2015.09.038. *(canonical Cas12a/Cpf1 TTTV PAM targeting and crRNA architecture framework used in this workflow)*

  - Kim, Hui K., Minju Song, Ji Hyun Lee, Hyunho Ko, Youngbin Kim, Seung Hwan Jang, and Jin-Soo Kim. 2016. “In Vivo High-Throughput Profiling of CRISPR–Cpf1 Activity.” Nature Methods 14 (2): 153–159. https://doi.org/10.1038/nmeth.4104. *(Cas12a activity and targeting-efficiency benchmark supporting spacer-selection heuristics and PAM-based targeting logic)*
  
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
        Gemini: The CRISPR guide RNA design for targeting the **EMX1** gene in *Homo sapiens* has been completed.

### Best Guide RNA Candidate
The top-ranked candidate was selected based on a combined score of on-target efficiency (GC content, absence of poly-T runs, and PAM-proximal G) and off-target safety (using the CFD scoring system to ensure zero predicted off-target sites).

*   **Protospacer:** `GCTCTCCAAACCACGGCTCC`
*   **PAM:** `CGG`
*   **gRNA Sequence:** `GCUCUCCAAACCACGGCUCCGUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC`
*   **Rationale:** This guide achieved a perfect score of 3/3, demonstrating optimal GC content (65%) and no predicted off-target binding sites.

---

### Next Steps
Would you like me to proceed with one of the following?

1.  **Generate a construction file** — to create a structured cloning workflow record.
2.  **Generate a lab sheet** — for a bench-ready, step-by-step experimental protocol (requires a construction file first).
3.  **Run off-target validation** — for a detailed CFD analysis.
4.  **Predict editing efficiency** — for a pre-experiment performance estimate.

Please let me know how you would like to proceed.

You: design cloning oligos 

[Tool call] → crispr_predict_editing_efficiency
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "pam": "CGG"
}
[Tool result] ← crispr_predict_editing_efficiency:
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "nuclease": "cas9",
  "raw_score": 5.5,
  "on_target_efficiency_pct": 40.2,
  "confidence_range": [
    25.2,
    55.2
  ],
  "delivery": "plasmid",
  "outcome": "nhej",
  "delivery_multiplier": 0.55,
  "outcome_multiplier": 1.0,
  "feature_contributions": {
    "position_specific_nucleotides": -0.5,
    "pam_context": 5.0,
    "gc_content": 1.0,
    "polyt_penalty": 0.0
  },
  "interpretation": "Predicted NHEJ knockout efficiency: 40.2% (range 25.2-55.2%) for CAS9 delivered via plasmid. This is moderate \u2014 usable, \u2026 (228 chars)",
  "warnings": []
}

[Tool call] → crispr_predict_offtargets
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "reference": "GAGCGCTCTCCAAACCACGGCTCCCGGCGCTCAGGCGGGCCGCTGCCAAGACCCGGCCTGGAGTCCCCGCAGAGTTGCGCGGCGCACGGACCCCGTGGCCTTGGGGCGTCAGGAGGCCCA\u2026 (18703 chars)"
}
[Tool result] ← crispr_predict_offtargets:
{
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "nuclease": "cas9",
  "reference_length": 18703,
  "strands_scanned": 2,
  "sites_evaluated": 37368,
  "offtarget_sites": [
    {
      "position": 4,
      "strand": "+",
      "sequence": "GCTCTCCAAACCACGGCTCC",
      "mismatches": 0,
      "mismatch_positions": [],
      "seed_mismatches": 0,
      "has_pam": true,
      "risk": "HIGH",
      "cfd_score": 1.0
    }
  ],
  "high_risk_count": 1,
  "aggregate_offtarget_cfd": 0.0,
  "max_offtarget_cfd": 0.0,
  "specificity_summary": "1 on-target site found (0 mismatches). 0 additional site(s) within 3 mismatches \u2014 review before proceeding."
}

[Tool call] → crispr_design_cloning_oligos
{
  "bottom_overhang": "AAAC",
  "prepend_g": true,
  "construct_name": "pX330-EMX1-guide",
  "protospacer": "GCTCTCCAAACCACGGCTCC",
  "top_overhang": "CACC",
  "vector": "px330",
  "enzyme": "BbsI"
}
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
  "top_oligo": "CACCGCTCTCCAAACCACGGCTCC",
  "bottom_oligo": "AAACGGAGCCGTGGTTTGGAGAGC",
  "end_structure": "two sticky-end overhangs",
  "g_prepended": false,
  "final_protospacer": "GCTCTCCAAACCACGGCTCC",
  "construction_file_inputs": {
    "construct_name": "pX330-EMX1-guide",
    "assembly_strategy": "TypeIISOligoCloning",
    "backbone_name": "px330",
    "backbone_sequence": "GAGGGCCTATTTCCCATGATTCCTTCATATTTGCATATACGATACAAGGCTGTTAGAGAGATAATTGGAATTAATTTGACTGTAAACACAAAGATATTAGTACAAAATACGTGACGTAGA\u2026 (8484 chars)",
    "insert_name": "px330_annealed_guide_insert",
    "insert_sequence": "GCTCTCCAAACCACGGCTCC",
    "insert_forward_primer_name": "",
    "insert_forward_primer_sequence": "",
    "insert_reverse_primer_name": "",
    "insert_reverse_primer_sequence": "",
    "top_oligo_name": "px330_guide_top",
    "top_oligo_sequence": "CACCGCTCTCCAAACCACGGCTCC",
    "bottom_oligo_name": "px330_guide_bottom",
    "bottom_oligo_sequence": "AAACGGAGCCGTGGTTTGGAGAGC",
    "top_overhang": "CACC",
    "bottom_overhang": "AAAC",
    "vector_forward_primer_name": "",
    "vector_forward_primer_sequence": "",
    "vector_reverse_primer_name": "",
    "vector_reverse_primer_sequence": "",
    "enzyme": "BbsI",
    "cell_strain": "Any mammalian",
    "selection": "Amp",
    "temperature_c": 37,
    "notes": "Vector 'pX330' \u2014 BbsI (TypeIIS) cloning. All-in-one SpCas9 + sgRNA vector.  BbsI digest leaves CACC/AAAC overhangs.  5\u2032G\u2026 (603 chars)"
  }
}

Gemini: The cloning oligos for the EMX1-targeting guide RNA (`GCTCTCCAAACCACGGCTCC`) into the **pX330** vector have been designed. This workflow uses BbsI-mediated annealed-oligo cloning.

### Designed Cloning Oligos
You will need to order these two synthetic DNA oligonucleotides:

*   **Top Oligo (`px330_guide_top`):** `5′-CACCGCTCTCCAAACCACGGCTCC-3′`
*   **Bottom Oligo (`px330_guide_bottom`):** `5′-AAACGGAGCCGTGGTTTGGAGAGC-3′`

---

### Construction File
This structured record captures the design parameters for your cloning workflow:

```text
Construct: pX330-EMX1-guide
Vector: px330 (Addgene #42230)
Method: TypeIISOligoCloning (BbsI)
Insert: GCTCTCCAAACCACGGCTCC

--- Cloning Setup ---
Enzyme: BbsI
Selection: Ampicillin
Host: Mammalian cells

--- Oligos ---
Top Oligo: CACCGCTCTCCAAACCACGGCTCC
Bottom Oligo: AAACGGAGCCGTGGTTTGGAGAGC

--- Notes ---
Vector 'pX330' — BbsI (TypeIIS) cloning. All-in-one SpCas9 + sgRNA vector.
BbsI digest leaves CACC/AAAC overhangs.
```

---

### Next Steps
All CRISPR design and validation steps are complete. Would you like me to generate:

**(a) a lab sheet** — a bench-ready, step-by-step protocol (which includes a colony picking estimate),
**(b) sequencing primer design** — to order flanking primers for post-experiment ICE/TIDE verification,
**(c) both**, or
**(d) neither**?

You: 

> **Notes on session above**
> - `crispr_predict_editing_efficiency` and `crispr_predict_offtargets` were not called directly by `run_full_crispr_workflow.py`; they were autonomously selected by Gemini as downstream validation tools after guide selection completed.
> - These validation tools were authored by other project contributors and are integrated into the MCP ecosystem as optional downstream analysis modules.
> - `predict_editing_efficiency` was used to estimate pre-experimental NHEJ editing performance for the selected EMX1 SpCas9 guide.
> - `predict_offtargets` was used to perform mismatch-based specificity scanning against the fetched EMX1 genomic locus sequence.
> - `design_cloning_oligos` was then called because the user explicitly requested guide-cloning preparation for the pX330 workflow.
</details>

> **Benchmark alignment — Hsu et al. 2013 / Ran et al. 2013:**  
> This workflow targets EMX1 in *Homo sapiens* using the pX330 mammalian SpCas9 system, matching one of the canonical benchmark loci from early CRISPR-Cas9 specificity studies. Alignment points include NGG PAM targeting, U6-driven sgRNA architecture, pX330 BbsI-mediated Type IIS cloning with CACC/AAAC overhangs, and low-off-target guide prioritization.


  
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
  
  ## 7. [_utils.py](modules/crispr_tools/tools/_utils.py) ( NOT callable by MCP, only for support, no test files ) 
  
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
