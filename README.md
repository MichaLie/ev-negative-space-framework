# EV Negative-Space Framework

This repository contains the code, frozen outputs, figures, and machine-readable supplementary artifacts supporting the manuscript:

`Negative-Space Mapping in Extracellular Vesicle Signaling: Cardiac Repair as a Cross-Disease Translational Priority`

It is organized as a public reproducibility package for inspection, reuse, and extension.

## Overview

The framework combines five linked layers:

1. PubMed-based pathway and cargo landscape mapping across neurodegeneration, tumor metastasis, and cardiac repair
2. Rule-based screening and quality control of pathway-context records
3. Full-text mechanistic adjudication of the highest-priority axes
4. RNA/protein cargo-confidence scoring from structured external evidence sources
5. Translational outputs including the ranked roadmap, clinical-trial gap snapshot, and experimental starter kit

## Repository layout

```text
scripts/
  systematic_ev_negative_space_pipeline.py
  cargo_confidence_rna_protein_pipeline.py
  robustness_hardening_analysis.py
  sensitivity_threshold_and_cardiac_relaxation.py
  supplementary_strengthening.py
  render_figure1_workflow.py

figures/
  figure1_workflow.png
  figure1_workflow.svg
  figure2_cargo_context_bars.png
  figure3_annual_growth_contexts.png
  figure4_pathway_heatmap.png
  figure5_directionality_matrix_final.png
  figure6_negative_space_priority.png
  figure7_cargo_confidence_scores.png
  figure8_negative_space_vs_cargo_confidence.png

supplement/
  query_registry.csv
  screening_table_pathway_context_top40.csv
  cargo_confidence_query_log.csv
  methods_manifest.json
  cargo_confidence_manifest.json
  robustness_hardening_manifest.json
  ev_nomenclature_sensitivity_table.csv
  cargo_weight_sensitivity_summary.csv
  mechanistic_rigor_summary.csv
  figureS1_screening_flow.png
  figureS2_mechanistic_confidence_final.png
  ... additional machine-readable supplementary outputs

root-level outputs/
  pubmed_pathway_context_counts_strict_2015_2026.csv
  pubmed_context_year_counts_2010_2026_strict.csv
  pubmed_cargo_context_counts_strict_2015_2026.csv
  negative_space_priority_index.csv
  cargo_confidence_rna_protein.csv
  pathway_enrichment_results.csv
  clinicaltrials_ev_snapshot.json
  ... additional final analysis outputs
```

## Main figures

The `figures/` directory contains all eight main figures from the study:

| # | File | Content |
|---|---|---|
| 1 | `figure1_workflow.png` | Workflow schematic of the five-layer framework |
| 2 | `figure2_cargo_context_bars.png` | Cargo-class literature density by disease context |
| 3 | `figure3_annual_growth_contexts.png` | Annual EV literature growth, 2010–2026 |
| 4 | `figure4_pathway_heatmap.png` | Pathway-level EV evidence density heatmap |
| 5 | `figure5_directionality_matrix_final.png` | Directionality signal by pathway and context |
| 6 | `figure6_negative_space_priority.png` | Negative-space priority index by EV pathway |
| 7 | `figure7_cargo_confidence_scores.png` | RNA/protein cargo-confidence scores for 31 candidates |
| 8 | `figure8_negative_space_vs_cargo_confidence.png` | Pathway negative-space priority versus mean cargo confidence |

The workflow schematic is reproducible via `scripts/render_figure1_workflow.py`.

## Dependencies

The scripts use a small Python stack:

- `pandas`
- `requests`
- `matplotlib`
- `seaborn`
- `scipy`

Install with:

```bash
pip install -r requirements.txt
```

## Reproducibility notes

- The frozen outputs in this repository correspond to the manuscript analysis window ending on `2026-02-15`.
- Some scripts call live external services such as PubMed E-utilities, Enrichr, and ClinicalTrials.gov-related APIs or snapshots. Later reruns are therefore date-sensitive.
- The CSV and JSON outputs in this repository are the canonical frozen artifacts for the reported analysis.

## Suggested reading order

1. `README.md`
2. `scripts/systematic_ev_negative_space_pipeline.py`
3. `scripts/render_figure1_workflow.py`
4. `supplement/query_registry.csv`
5. `negative_space_priority_index.csv`
6. `supplement/ev_nomenclature_sensitivity_table.csv`
7. `supplement/cargo_weight_sensitivity_summary.csv`
8. `supplement/mechanistic_rigor_summary.csv`

## Contact

Michaela Liegertova  
Centre for Nanomaterials and Biotechnology, Faculty of Science, Jan Evangelista Purkyne University in Usti nad Labem  
`michaela.liegertova@ujep.cz`
