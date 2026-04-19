# EV Negative-Space Framework

Public data-availability package for a reproducible extracellular vesicle negative-space analysis across neurodegeneration, tumor metastasis, and cardiac repair.

This repository contains frozen outputs, analysis scripts, figures, and machine-readable supplementary files.

## Layout

- `scripts/`: analysis and figure-generation scripts
- `figures/`: main figures
- `supplement/`: machine-readable supplementary tables, manifests, and supplementary figures
- root-level `.csv` and `.json` files: frozen analysis outputs

## Key files

- `negative_space_priority_index.csv`: pathway-level priority scores
- `cargo_confidence_rna_protein.csv`: candidate-level cargo-confidence scores
- `negative_space_with_cargo_confidence.csv`: integrated priority-confidence view
- `supplement/query_registry.csv`: query definitions and registry notes
- `clinicaltrials_ev_snapshot.json`: frozen ClinicalTrials.gov snapshot

## Reproducibility

- The frozen outputs correspond to the analysis window ending on `2026-02-15`.
- Some scripts call live external services such as PubMed E-utilities, Enrichr, and ClinicalTrials.gov. Reruns are therefore date-sensitive.
- The CSV and JSON files in this repository are the canonical machine-readable artifacts.

## Dependencies

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

## Contact

`michaela.liegertova@ujep.cz`
