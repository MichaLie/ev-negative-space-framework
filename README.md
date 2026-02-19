# Negative-Space Framework for Translational EV Research

**A Reproducible Framework for Mapping Translational Negative Space in Extracellular Vesicle Research**

Liegertová M, Malý J (2026). *Submitted to Patterns (Cell Press).*

## Overview

This repository contains all analysis scripts, data outputs, and supplementary materials for a quantitative evidence-synthesis framework that maps cross-domain pathway asymmetries in extracellular vesicle (EV) research across neurodegeneration, tumor metastasis, and cardiac repair.

The framework integrates five modular layers:
1. **Bibliometric asymmetry mapping** — PubMed-based pathway-level evidence counts
2. **Full-text mechanistic adjudication** — 42 papers, strict 8-criterion rubric
3. **Cargo-confidence scoring** — 31 RNA/protein candidates scored against miRBase, UniProt, ExoCarta
4. **Pathway enrichment consistency check** — Enrichr over-representation analysis
5. **Clinical-trial landscape survey** — ClinicalTrials.gov snapshot

## Repository Structure

```
scripts/                           # Analysis pipelines
  systematic_ev_negative_space_pipeline.py   # Main bibliometric + screening pipeline
  cargo_confidence_rna_protein_pipeline.py   # RNA/protein cargo scoring
  robustness_hardening_analysis.py           # Matched-stringency, formulation & partition analyses
  sensitivity_threshold_and_cardiac_relaxation.py  # Threshold sweep, continuous variant & cardiac relaxation

figures/                           # Main manuscript figures (1-7)

supplement/                        # Supplementary figures, tables, and methods
  Supplementary_Materials_research_article_v2.md
  Supplementary_Methods_and_Reporting.md
  figureS1_screening_flow.png
  figureS2_mechanistic_confidence_final.png
  query_registry.csv               # Complete query strings (Table S2)
  screening_table_pathway_context_top40.csv  # Full screening set (Table S3)
  experimental_starter_kit_top5_axes.csv     # Starter-kit table (Table S18)
  ... (additional CSV/JSON outputs)

*.csv / *.json                     # Root-level data outputs
  negative_space_priority_index.csv          # Priority rankings (Table S4)
  cargo_confidence_rna_protein.csv           # Cargo scores (Table S9)
  fulltext_mechanistic_adjudication_final.csv # Adjudication data (Table S8)
  clinicaltrials_ev_snapshot.json            # Clinical trials data
  ...
```

## Note on Figure Filenames

The analysis scripts (`scripts/systematic_ev_negative_space_pipeline.py`) generate intermediate figure files with working names (e.g., `figure1_pathway_heatmap_strict.png`). For the final manuscript, figures were renumbered and renamed to match the submission order (e.g., `figure1_cargo_context_bars.png`). The repository `figures/` directory contains the manuscript-ready versions. The script output filenames reflect the development-phase naming convention.

## Reproducibility

All queries, scoring weights, screening rules, and decision thresholds are fully specified in:
- `supplement/methods_manifest.json` — machine-readable protocol manifest
- `supplement/query_registry.csv` — complete PubMed query strings
- `supplement/cargo_confidence_manifest.json` — cargo scoring parameters

The pipeline is designed for temporal updating: re-running scripts with a new date window produces an updated atlas without methodological changes.

## Data Sources

- **PubMed/MEDLINE** (NCBI E-utilities) — bibliometric queries
- **ExoCarta** — EV cargo database
- **Vesiclepedia v5.1** — independent EV cargo concordance
- **miRBase** — miRNA identifier validation
- **UniProt** — protein identifier validation
- **Enrichr** — pathway enrichment analysis
- **ClinicalTrials.gov** — clinical trial landscape

## Key Results

- **Integrin/Src** ranks #1 across all 12 tested index formulations (median Spearman = 0.75)
- **Cardiac repair** is consistently the most evidence-sparse domain for pathways mature in cancer and neurodegeneration
- **Context-dependent functional inversion** (e.g., miR-21 pro-metastatic in tumors, cardioprotective after infarction) demands safety-gated cross-domain transfer

## License

Code: MIT License
Data and figures: CC-BY 4.0

## Contact

Michaela Liegertová — michaela.liegertova@ujep.cz
Centre for Nanomaterials and Biotechnology, Faculty of Science, Jan Evangelista Purkyne University in Usti nad Labem, Czech Republic
