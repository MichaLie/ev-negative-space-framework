#!/usr/bin/env python3
"""Generate supplementary robustness analyses and package-ready tables.

This script adds three manuscript-supporting analyses:
1. EV nomenclature sensitivity across alternative EV clauses.
2. Cargo-confidence weight sensitivity across plausible weighting schemes.
3. MISEV-aligned adjudication rigor summary from the full-text layer.

Outputs are written both to the repository supplement/ directory and, when
available, into sibling submission package directories.
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import requests
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
SUPP = ROOT / "supplement"
SUPP.mkdir(parents=True, exist_ok=True)

PACKAGE_ROOT = ROOT.parent
PACKAGE_TABLES = PACKAGE_ROOT / "05_SUPPLEMENTARY_TABLES"

RUN_DATE = "2026-03-09"
MINDATE = "2015/01/01"
MAXDATE = "2026/02/15"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
SMOOTHING = 1.0

CONTEXTS = {
    "neurodegeneration": '(neurodegeneration OR neurodegenerative OR Alzheimer OR Parkinson OR "amyotrophic lateral sclerosis" OR ALS OR frontotemporal dementia OR Huntington)',
    "tumor_metastasis": "((cancer OR tumor OR tumour) AND (metastasis OR metastatic OR premetastatic OR organotropism OR dissemination))",
    "cardiac_repair": "((myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury) AND (repair OR regeneration OR cardioprotection OR post-infarction))",
}

PATHWAYS = {
    "mTOR": "(mTOR OR PI3K OR AKT)",
    "Wnt": "(Wnt OR beta-catenin OR CTNNB1)",
    "Notch": "(Notch OR NOTCH1)",
    "NF-kB": '("NF-kappa B" OR NF-kB OR NFKB1)',
    "Complement": "(complement OR C3 OR C5)",
    "Sphingolipid/Ceramide": "(sphingolipid OR ceramide OR sphingomyelinase OR S1P)",
    "Autophagy": "(autophagy OR mitophagy OR lysosome)",
    "Hypoxia/HIF-1": "(hypoxia OR HIF-1 OR HIF1A)",
    "Integrin/Src": "(integrin OR SRC)",
    "TGF-beta": '("TGF-beta" OR TGFB1)',
}

EV_CLAUSES = {
    "current_baseline": "(extracellular vesicle OR exosome OR exosomes)",
    "legacy_exosome": "(exosome OR exosomes)",
    "umbrella_ev": '("extracellular vesicle" OR "extracellular vesicles" OR "small extracellular vesicle" OR "small extracellular vesicles")',
    "expanded_ev": '("extracellular vesicle" OR "extracellular vesicles" OR exosome OR exosomes OR microvesicle OR microvesicles OR ectosome OR ectosomes OR "small extracellular vesicle" OR "small extracellular vesicles")',
}

CARGO_WEIGHTS = {
    "baseline": {"identifier": 0.25, "exocarta": 0.25, "breadth": 0.25, "mechanistic": 0.15, "biomarker": 0.10},
    "equal": {"identifier": 0.20, "exocarta": 0.20, "breadth": 0.20, "mechanistic": 0.20, "biomarker": 0.20},
    "exocarta_heavy": {"identifier": 0.20, "exocarta": 0.35, "breadth": 0.20, "mechanistic": 0.15, "biomarker": 0.10},
    "breadth_heavy": {"identifier": 0.20, "exocarta": 0.20, "breadth": 0.35, "mechanistic": 0.15, "biomarker": 0.10},
    "mechanistic_heavy": {"identifier": 0.20, "exocarta": 0.20, "breadth": 0.20, "mechanistic": 0.30, "biomarker": 0.10},
    "biomarker_heavy": {"identifier": 0.20, "exocarta": 0.20, "breadth": 0.20, "mechanistic": 0.15, "biomarker": 0.25},
}

CARGO_COMPONENTS = {
    "identifier": "identifier_score",
    "exocarta": "exocarta_score",
    "breadth": "breadth_score",
    "mechanistic": "mechanistic_score",
    "biomarker": "biomarker_norm",
}

RIGOR_COLUMNS = [
    "in_vivo",
    "in_vitro",
    "human_samples",
    "causal_perturbation",
    "ev_isolation",
    "ev_markers",
    "context_endpoint_match",
]


def esearch_count(term: str) -> int:
    response = requests.get(
        BASE,
        params={
            "db": "pubmed",
            "term": term,
            "retmode": "json",
            "rettype": "count",
            "mindate": MINDATE,
            "maxdate": MAXDATE,
            "datetype": "pdat",
        },
        timeout=30,
    )
    response.raise_for_status()
    return int(response.json()["esearchresult"]["count"])


def compute_priority_table(counts_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pathway, grp in counts_df.groupby("pathway"):
        vals = grp.set_index("context")["count"]
        max_count = int(vals.max())
        min_count = int(vals.min())
        contexts_ge100 = int((vals >= 100).sum())
        ratio = (max_count + SMOOTHING) / (min_count + SMOOTHING)
        priority = ratio * (contexts_ge100 / 3)
        rows.append(
            {
                "pathway": pathway,
                "max_count": max_count,
                "min_count": min_count,
                "contexts_ge100": contexts_ge100,
                "smoothed_ratio": ratio,
                "priority_index": priority,
            }
        )
    out = pd.DataFrame(rows).sort_values(
        ["priority_index", "pathway"], ascending=[False, True]
    ).reset_index(drop=True)
    out["rank"] = range(1, len(out) + 1)
    return out


def build_nomenclature_sensitivity() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base_counts = pd.read_csv(ROOT / "pubmed_pathway_context_counts_strict_2015_2026.csv")
    base_counts = base_counts[["pathway", "context", "count"]].copy()
    base_counts["ev_clause_label"] = "current_baseline"
    base_counts["ev_clause"] = EV_CLAUSES["current_baseline"]

    count_frames = [base_counts]
    for label, clause in EV_CLAUSES.items():
        if label == "current_baseline":
            continue
        rows = []
        for context, context_clause in CONTEXTS.items():
            for pathway, pathway_clause in PATHWAYS.items():
                query = f"{clause} AND {context_clause} AND {pathway_clause}"
                rows.append(
                    {
                        "ev_clause_label": label,
                        "ev_clause": clause,
                        "context": context,
                        "pathway": pathway,
                        "count": esearch_count(query),
                        "query": query,
                    }
                )
                time.sleep(0.12)
        count_frames.append(pd.DataFrame(rows))

    counts = pd.concat(count_frames, ignore_index=True)

    priority_frames = []
    for label, grp in counts.groupby("ev_clause_label"):
        priority = compute_priority_table(grp[["pathway", "context", "count"]])
        priority["ev_clause_label"] = label
        priority_frames.append(priority)
    priority = pd.concat(priority_frames, ignore_index=True)

    rank_maps = {
        label: grp.set_index("pathway")["rank"].sort_index()
        for label, grp in priority.groupby("ev_clause_label")
    }
    base_rank = rank_maps["current_baseline"]
    base_top5 = set(
        priority[priority["ev_clause_label"] == "current_baseline"]
        .sort_values("rank")
        .head(5)["pathway"]
    )

    summary_rows = []
    for label, grp in priority.groupby("ev_clause_label"):
        ranks = grp.set_index("pathway")["rank"].sort_index()
        top5 = grp.sort_values("rank").head(5)["pathway"].tolist()
        top2 = grp.sort_values("rank").head(2)["pathway"].tolist()
        summary_rows.append(
            {
                "ev_clause_label": label,
                "spearman_vs_current_baseline": round(
                    float(spearmanr(base_rank, ranks).statistic), 4
                ),
                "top5_overlap_vs_current_baseline": len(base_top5 & set(top5)),
                "top2_pathways": "; ".join(top2),
                "top5_pathways": "; ".join(top5),
            }
        )
    summary = pd.DataFrame(summary_rows).sort_values("ev_clause_label")

    wide_counts = counts.pivot_table(
        index=["ev_clause_label", "pathway"],
        columns="context",
        values="count",
        aggfunc="first",
    ).reset_index()
    table = wide_counts.merge(
        priority[
            [
                "ev_clause_label",
                "pathway",
                "contexts_ge100",
                "smoothed_ratio",
                "priority_index",
                "rank",
            ]
        ],
        on=["ev_clause_label", "pathway"],
        how="left",
    ).sort_values(["ev_clause_label", "rank", "pathway"])
    return counts, table, summary


def build_cargo_weight_sensitivity() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(ROOT / "cargo_confidence_rna_protein.csv")
    long_rows = []
    summary_rows = []

    base = None
    base_top10 = set()

    for scheme, weights in CARGO_WEIGHTS.items():
        score = sum(df[CARGO_COMPONENTS[k]] * v for k, v in weights.items()) * 100
        ranked = df[["cargo_type", "candidate", "pathway"]].copy()
        ranked["scheme"] = scheme
        ranked["score"] = score.round(6)
        ranked["rank"] = ranked["score"].rank(method="min", ascending=False).astype(int)
        ranked = ranked.sort_values(["rank", "candidate"]).reset_index(drop=True)
        long_rows.append(ranked)

        if scheme == "baseline":
            base = ranked.set_index("candidate")["rank"].sort_index()
            base_top10 = set(ranked.head(10)["candidate"])
            summary_rows.append(
                {
                    "scheme": scheme,
                    "spearman_vs_baseline": 1.0,
                    "top10_overlap_vs_baseline": 10,
                    "top5_candidates": "; ".join(ranked.head(5)["candidate"]),
                }
            )
        else:
            cur = ranked.set_index("candidate")["rank"].sort_index()
            summary_rows.append(
                {
                    "scheme": scheme,
                    "spearman_vs_baseline": round(
                        float(spearmanr(base, cur).statistic), 4
                    ),
                    "top10_overlap_vs_baseline": len(base_top10 & set(ranked.head(10)["candidate"])),
                    "top5_candidates": "; ".join(ranked.head(5)["candidate"]),
                }
            )

    long_df = pd.concat(long_rows, ignore_index=True)
    summary_df = pd.DataFrame(summary_rows).sort_values("scheme")
    return long_df, summary_df


def build_rigor_summary() -> pd.DataFrame:
    df = pd.read_csv(ROOT / "fulltext_mechanistic_adjudication_final.csv")
    frames = []

    overall = {
        "summary_level": "overall",
        "pathway": "all",
        "n_papers": len(df),
        "mean_mechanistic_score_strict": round(df["mechanistic_score_strict"].mean(), 3),
        "high_confidence_count": int((df["mechanistic_confidence_strict"] == "High").sum()),
        "moderate_confidence_count": int((df["mechanistic_confidence_strict"] == "Moderate").sum()),
    }
    for col in RIGOR_COLUMNS:
        overall[f"{col}_count"] = int(df[col].sum())
        overall[f"{col}_pct"] = round(float(df[col].mean() * 100), 1)
    frames.append(pd.DataFrame([overall]))

    by_pathway = (
        df.groupby("pathway")
        .agg(
            n_papers=("pathway", "size"),
            mean_mechanistic_score_strict=("mechanistic_score_strict", "mean"),
            high_confidence_count=("mechanistic_confidence_strict", lambda s: int((s == "High").sum())),
            moderate_confidence_count=("mechanistic_confidence_strict", lambda s: int((s == "Moderate").sum())),
            **{f"{col}_count": (col, "sum") for col in RIGOR_COLUMNS},
            **{f"{col}_pct": (col, "mean") for col in RIGOR_COLUMNS},
        )
        .reset_index()
    )
    by_pathway["summary_level"] = "by_pathway"
    by_pathway["mean_mechanistic_score_strict"] = by_pathway["mean_mechanistic_score_strict"].round(3)
    for col in RIGOR_COLUMNS:
        by_pathway[f"{col}_pct"] = (by_pathway[f"{col}_pct"] * 100).round(1)

    frames.append(by_pathway)
    return pd.concat(frames, ignore_index=True)


def write_outputs() -> None:
    counts, nomenclature_table, nomenclature_summary = build_nomenclature_sensitivity()
    cargo_long, cargo_summary = build_cargo_weight_sensitivity()
    rigor_summary = build_rigor_summary()

    counts.to_csv(SUPP / "ev_nomenclature_sensitivity_counts.csv", index=False)
    nomenclature_table.to_csv(SUPP / "ev_nomenclature_sensitivity_table.csv", index=False)
    nomenclature_summary.to_csv(SUPP / "ev_nomenclature_rank_stability.csv", index=False)
    cargo_long.to_csv(SUPP / "cargo_weight_sensitivity_candidates.csv", index=False)
    cargo_summary.to_csv(SUPP / "cargo_weight_sensitivity_summary.csv", index=False)
    rigor_summary.to_csv(SUPP / "mechanistic_rigor_summary.csv", index=False)

    if PACKAGE_TABLES.exists():
        nomenclature_table.to_csv(
            PACKAGE_TABLES / "TableS5_ev_nomenclature_sensitivity.csv", index=False
        )
        rigor_summary.to_csv(
            PACKAGE_TABLES / "TableS6_mechanistic_rigor_summary.csv", index=False
        )
        cargo_summary.to_csv(
            PACKAGE_TABLES / "TableS14_cargo_weight_sensitivity.csv", index=False
        )

    print("Wrote:")
    print(" - supplement/ev_nomenclature_sensitivity_counts.csv")
    print(" - supplement/ev_nomenclature_sensitivity_table.csv")
    print(" - supplement/ev_nomenclature_rank_stability.csv")
    print(" - supplement/cargo_weight_sensitivity_candidates.csv")
    print(" - supplement/cargo_weight_sensitivity_summary.csv")
    print(" - supplement/mechanistic_rigor_summary.csv")
    if PACKAGE_TABLES.exists():
        print(" - package TableS5, TableS6, TableS14")


if __name__ == "__main__":
    write_outputs()
