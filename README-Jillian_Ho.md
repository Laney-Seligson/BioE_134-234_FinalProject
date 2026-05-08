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
- What it does: Converts a structured construction file into a printable, bench-ready lab packet matching the UC Berkeley LabPlanner format with reagent volumes, thermocycler programs, and citations on every step. Renders one section per physical bench session — PCR (with auto-mastermix for 4+ parallel samples, 10% dead volume, shared-reagent detection across reactions), Gel + DpnI, Zymo cleanup, Digest, Ligate (or Assemble for Golden Gate / Gibson / Type IIS oligo cloning), Transform (with auto-inferred rescue step: any antibiotic except β-lactams (Amp/Carb) and yeast auxotrophic markers (URA3/LEU2/HIS3/TRP1/...) requires a 60-min recovery before plating), Plate, Pick, Miniprep, and plasmid Sanger sequencing (instructions match the professor's LabPlanner format verbatim — 487 µL ddH2O / 13.3 µL of 100 µM oligo dilution recipe; 237 Stanley Hall pickup; "save the stock at..." reminder; embeds the protospacer as the expected read content so the user knows what to look for in the trace). Auto-extracts the protospacer from the construction record's parts list, then chains into predict_editing_efficiency to size the colony plan with a guide-specific number instead of a preset. Optional CRISPR-specific delivery section (RNP / electroporation / lipofection / IVT sgRNA) and a separate post-edit GENOMIC verification section that runs verify_edit on the genomic reference and prepares the ICE/TIDE submission protocol. Outputs both human-readable text and a flat-row TSV with 17 columns suitable for spreadsheet import. Supports custom storage labels via location_overrides (boxA → "Bench-3 freezer / Rack 2", etc.) so lab sheets match the user's actual ELN/freezer system.
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