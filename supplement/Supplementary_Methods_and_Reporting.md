# Supplementary Methods and Reporting

## S1. Protocol Summary
This study implemented a reproducible, script-based literature mining protocol for EV cargo/signaling comparison across neurodegeneration, tumor metastasis, and cardiac repair. All code and outputs are located at https://github.com/MichaLie/ev-negative-space-framework.

Run date: February 15, 2026  
Date window: January 1, 2015 to February 15, 2026  
Data source: PubMed (NCBI E-utilities)

## S2. Query Design
Core EV clause: `(extracellular vesicle OR exosome OR exosomes)`

Context clauses were strict and domain-targeted:
- Neurodegeneration: neurodegenerative/Alzheimer/Parkinson/ALS/FTD/Huntington terms.
- Tumor metastasis: cancer/tumor AND metastatic/dissemination terms.
- Cardiac repair: myocardial/cardiac injury terms AND repair/regeneration/cardioprotection terms.

Cargo axes:
- proteins/proteomics
- miRNA

The translational confidence layer was intentionally restricted to RNA and protein/marker cargo. Lipid/sphingolipid terms were retained only for background pathway mapping in the initial atlas.

Pathway axes:
- mTOR
- Wnt
- Notch
- NF-kB
- Complement
- Sphingolipid/Ceramide
- Autophagy
- Hypoxia/HIF-1
- Integrin/Src
- TGF-beta

Full query registry is provided in:  
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/query_registry.csv`

## S3. Screening Rules
For each pathway-context cell, up to 40 records were retrieved by relevance and annotated with metadata and publication types.

Inclusion logic:
- `include_primary_high`: publication type consistent with primary research, relevance score >= 3, explicit context hit, explicit pathway hit.
- `include_contextual`: relevance score >= 2 with explicit context hit.
- otherwise `exclude_low_relevance`.

Fields retained: PMID, DOI, title, journal, publication date, publication type, abstract, relevance components.

Screening outputs:
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/screening_table_pathway_context_top40.csv`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/screening_summary.json`

## S4. Quantitative Outputs
Primary matrices:
- Cargo-context counts: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/pubmed_cargo_context_counts_strict_2015_2026.csv`
- Pathway-context counts: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/pubmed_pathway_context_counts_strict_2015_2026.csv`
- Year trends: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/pubmed_context_year_counts_2010_2026_strict.csv`

Derived prioritization:
- Negative-space index: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/negative_space_priority_index.csv`
- Phase-I roadmap: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/cross_field_phase1_roadmap.csv`
- Confidence-weighted roadmap: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/negative_space_with_cargo_confidence.csv`
- Robustness hardening manifest: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/robustness_hardening_manifest.json`

## S5. Figures
Main figures:
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure1_cargo_context_bars.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure2_annual_growth_contexts.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure3_pathway_heatmap.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure4_directionality_matrix_final.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure5_negative_space_priority.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure6_cargo_confidence_scores.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/figures/figure7_negative_space_vs_cargo_confidence.png`

Supplementary figures:
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/figureS1_screening_flow.png`
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/figureS2_mechanistic_confidence_final.png`

## S6. Reproducibility Manifest
Machine-readable manifest:  
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/methods_manifest.json`

RNA/protein confidence manifest:
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/cargo_confidence_manifest.json`

ClinicalTrials.gov snapshot (accessed February 15, 2026):
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/clinicaltrials_ev_snapshot.json`

## S7. Limitations of Literature-Mining Approach
- Bibliometric density is a proxy for research focus, not biological effect size.
- Indexing and keyword overlap can introduce residual off-topic retrieval despite strict context logic.
- Full-text adjudication is required for effect-size synthesis and formal causal claims.
- RNA/protein confidence scores are panel-based and depend on database completeness (miRBase/UniProt/ExoCarta) and query semantics.

## S8. Full-Text Adjudication (Top-Priority Set)
We added a full-text layer for the top-priority paper set:
- Priority set: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/fulltext_request_top30.csv`
- Coverage report: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/fulltext_request_top30_coverage.csv`
- Final adjudication set: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/fulltext_final_adjudication_set.csv`
- Structured adjudication table: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/fulltext_mechanistic_adjudication_final.csv`
- Directionality matrix: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/directionality_matrix_final.csv`
- NA-aware directionality means and n-by-cell source data: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/directionality_proxy_scores_with_counts.csv`
- Confidence summary by pathway: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/mechanistic_confidence_by_pathway_final.csv`
- Manual extraction seed table: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/effect_size_causal_matrix_seed.csv`

Coverage achieved: 42/42 full texts in the finalized adjudication set across 5 pathway axes (Integrin/Src n=11, Autophagy n=10, TGF-beta n=8, Hypoxia/HIF-1 n=7, Wnt n=6). The initial 32-paper set (4 pathways) was expanded to 42 papers by adding 7 Hypoxia/HIF-1 papers (4 donor-context tumor metastasis, 3 target-context neurodegeneration) and 3 Wnt target-context cardiac repair papers. Pre-specified replacement papers were added for the unavailable DOI `10.1016/j.cell.2025.01.025` (PMID 39938515). The gap/replacement rationale is documented in:
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/fulltext_gap_and_replacement_plan.md`.

Neurodegeneration context audit file (substring-vs-boundary screening check):
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/neuro_context_false_positives_wordboundary.csv`

## S9. RNA/Protein Cargo-Confidence Pipeline
Pipeline script:
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/scripts/cargo_confidence_rna_protein_pipeline.py`

Panel definition:
- 31 prespecified candidates (12 RNA, 19 proteins)
- Pathway assignment: Integrin/Src, TGF-beta, Wnt, mTOR, NF-kB, Complement, Autophagy, Hypoxia/HIF-1

Primary data sources:
- ExoCarta cargo archives (miRNA and protein/mRNA tables)
- miRBase mature FASTA (RNA identifier validation)
- UniProt reviewed human entries (protein identifier validation)
- PubMed E-utilities (context, biomarker, mechanistic query strata)
- Final full-text adjudication table for title-level support feature

Scoring equation (0-100 scale):
- Identifier validity: 0.25
- ExoCarta evidence: 0.25
- Context breadth / literature breadth: 0.25
- Mechanistic evidence: 0.15
- Biomarker evidence: 0.10

Tier thresholds:
- High: score >= 70
- Moderate: score 40-69.9
- Low: score < 40

Outputs:
- Candidate table: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/cargo_confidence_rna_protein.csv`
- Pathway summary: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/pathway_cargo_confidence_summary.csv`
- Query log: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/cargo_confidence_query_log.csv`
- Confidence-weighted negative-space table: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/negative_space_with_cargo_confidence.csv`

## S10. Query/Formulation Robustness Hardening (A + C + F)
Pipeline script:
`https://github.com/MichaLie/ev-negative-space-framework/blob/main/scripts/robustness_hardening_analysis.py`

Run date: February 15, 2026  
Window: January 1, 2015 to February 15, 2026

### A) Matched-stringency sensitivity across all three contexts
Three aligned query tiers were run for neurodegeneration, tumor metastasis, and cardiac repair:
- strict
- relaxed
- broadest

Outputs:
- Tiered counts by pathway/context: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/sensitivity_matched_query_counts_by_tier.csv`
- Tiered priority/rank table: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/sensitivity_matched_query_priority_ranks.csv`
- Rank stability summary (top3 + Spearman vs strict): `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/sensitivity_matched_query_rank_stability.csv`

### C) Priority-index formulation stability
Alternative formulations were benchmarked across:
- Smoothing constants: 0.5, 1.0, 2.0
- Score families:
  - ratio only
  - ratio × log(total evidence)
  - ratio × log(donor-context evidence)
  - ratio × maturity fraction (baseline family)

Outputs:
- Full formulation matrix + per-formulation ranks: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/priority_index_alternative_formulations.csv`
- Spearman-vs-baseline summary with top3 per formulation: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/priority_index_rank_stability_summary.csv`

### F) Integrin/Src partition analysis
Under strict context clauses, the combined Integrin/Src query was partitioned into:
- integrin-only (integrin NOT SRC)
- Src-only (SRC NOT integrin)
- both (integrin AND SRC)

Output:
- Context-level partition counts, percentages, and full query strings: `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/integrin_src_partition_counts.csv`

## S11. Experimental Starter-Kit Table (B)
To convert top-ranked negative-space signals into falsifiable first-pass experiments, we created a structured starter-kit table for the top five strict-priority axes (Integrin/Src, Autophagy, TGF-beta, Wnt, Hypoxia/HIF-1). Each row contains:
- transfer direction
- donor/target anchor PMIDs
- starter EV source class
- recipient model stack
- primary efficacy endpoints
- mandatory safety gates
- minimum EV characterization
- explicit falsification criteria

Output:
- `https://github.com/MichaLie/ev-negative-space-framework/blob/main/supplement/experimental_starter_kit_top5_axes.csv`
