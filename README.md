# Negative-Space Framework for Translational EV Research

This repository contains the analysis scripts, frozen outputs, and figures supporting the manuscript:

`A Reproducible Framework for Mapping Translational Evidence Gaps in Extracellular Vesicle Research`

It is intended as the public computational artifact behind the paper and is organized for inspection, reuse, and reproducibility.

## Scope

Included:

- analysis scripts used to generate the bibliometric, cargo-confidence, robustness, and JEV-strengthening outputs
- frozen CSV/JSON outputs used by the manuscript and supplementary tables
- final main figures and supplementary figures
- machine-readable manifests and query registries needed for inspection and reruns

Deliberately excluded:

- manuscript text and journal submission files
- cover letter and TOC submission assets
- narrative supplementary prose
- internal planning/request files that are not part of the scientific artifact
- submission-packaging scripts and assets

## Repository layout

```text
scripts/
  systematic_ev_negative_space_pipeline.py
  cargo_confidence_rna_protein_pipeline.py
  robustness_hardening_analysis.py
  sensitivity_threshold_and_cardiac_relaxation.py
  jev_submission_strengthening.py

figures/
  figure1_cargo_context_bars.png
  figure2_annual_growth_contexts.png
  figure3_pathway_heatmap.png
  figure4_directionality_matrix_final.png
  figure5_negative_space_priority.png
  figure6_cargo_confidence_scores.png
  figure7_negative_space_vs_cargo_confidence.png

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

The path layout intentionally stays close to the original working repository so the scripts can still locate the frozen outputs they expect.

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
- Some scripts call live external services such as PubMed E-utilities, Enrichr, and ClinicalTrials.gov-related APIs or snapshots. Those reruns are therefore date-sensitive.
- The frozen CSV/JSON outputs in this repository are the canonical paper-supporting artifacts; rerunning online queries later may not reproduce the exact same counts.

## Suggested reading order

1. `README.md`
2. `scripts/systematic_ev_negative_space_pipeline.py`
3. `supplement/query_registry.csv`
4. `negative_space_priority_index.csv`
5. `scripts/jev_submission_strengthening.py`
6. `supplement/ev_nomenclature_sensitivity_table.csv`
7. `supplement/cargo_weight_sensitivity_summary.csv`
8. `supplement/mechanistic_rigor_summary.csv`

## Contact

Michaela Liegertova  
Centre for Nanomaterials and Biotechnology, Faculty of Science, Jan Evangelista Purkyne University in Usti nad Labem  
`michaela.liegertova@ujep.cz`
