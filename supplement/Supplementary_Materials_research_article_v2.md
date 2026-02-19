# Supplementary Materials

## Supplementary Methods

### Literature-mining framework
The study used a reproducible PubMed E-utilities pipeline with a date window from January 1, 2015 to February 15, 2026. The core extracellular vesicle query was crossed with three context clauses (neurodegeneration, tumor metastasis, and cardiac repair), two cargo axes (proteins/proteomics and miRNA), and ten pathway axes.

### Screening and adjudication
For each pathway-context cell, up to 40 records were retrieved by relevance and screened using prespecified context/pathway matching and publication-type filters. A full-text adjudication layer evaluated 42 papers across five top-priority axes with structured extraction fields and a strict mechanistic confidence rubric.

### Priority framework and robustness
Pathway prioritization used a smoothed asymmetry ratio multiplied by the mature-context fraction (>=100 records). Robustness analyses included cardiac-only relaxation sensitivity, matched-stringency tiers across all contexts, threshold variation, alternative index formulations and smoothing constants, and Integrin/Src partitioning.

### Cargo-confidence and enrichment consistency check
A prespecified 31-candidate RNA/protein panel was scored using identifier validity, EV database support, context breadth, mechanistic support, and biomarker enrichment. Enrichr over-representation analysis confirmed internal consistency of the protein panel's pathway annotations.

## Supplementary Results Summary

### Full-text adjudication
The adjudication set included 42/42 available full texts across Integrin/Src, Autophagy, TGF-beta, Hypoxia/HIF-1, and Wnt axes. Confidence distribution was 38 High and 4 Moderate under the strict rubric.

### Matched-stringency analysis
Matched strict/relaxed/broadest context tiers showed moderate rank reordering versus strict (Spearman 0.6121 and 0.6364), with stable overlap of the strict high-priority core.

### Index formulation stability
Across 12 formulation/smoothing variants, Integrin/Src remained rank #1 in all formulations.

### Integrin/Src partition
Integrin-only records represented the majority component of the combined Integrin/Src signal in all three contexts.

## Supplementary Figures

**Supplementary Figure S1.** Screening workflow summary from retrieval through deduplication and relevance classification.

**Supplementary Figure S2.** Mechanistic confidence distribution in the 42-paper full-text adjudication set.

## Supplementary Tables and Files

*Note: Table numbering is non-consecutive (e.g., S4 to S8) to preserve cross-reference stability with the analysis pipeline outputs. Gaps reflect tables removed during iterative refinement.*

- Supplementary Table S2: Query registry
- Supplementary Table S3: Screening table (top-40 per pathway-context cell)
- Supplementary Table S4: Cross-field phase-I roadmap
- Supplementary Table S8: Full-text mechanistic adjudication table
- Supplementary Table S9: Cargo-confidence query log
- Supplementary Table S11: Cardiac-only query relaxation sensitivity
- Supplementary Table S12: Maturity-threshold sensitivity
- Supplementary Table S13: Enrichr pathway enrichment results
- Supplementary Table S15: Matched-stringency ranking outputs
- Supplementary Table S16: Alternative index formulation stability outputs
- Supplementary Table S17: Integrin/Src partition outputs
- Supplementary Table S18: Experimental starter-kit for top five strict-priority axes
- Supplementary File S10: Cargo-confidence manifest
- Supplementary File S19: Robustness hardening manifest
- Supplementary File S20: Directionality means and per-cell counts (NA-aware matrix source data)
- Supplementary File S21: Neurodegeneration context false-positive audit under word-boundary matching

## Data and Code Availability

All analysis scripts and output artifacts used in this manuscript are publicly available at https://github.com/MichaLie/ev-negative-space-framework.
