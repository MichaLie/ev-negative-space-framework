#!/usr/bin/env python3
"""Generate robustness artifacts for A + C + F hardening bundle.

A) Matched-stringency context sensitivity across strict/relaxed/broadest tiers.
C) Priority-index rank stability across alternative formulations and smoothing.
F) Integrin/Src partition analysis under strict context queries.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Dict, List

import pandas as pd
import requests

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MINDATE = "2015/01/01"
MAXDATE = "2026/02/15"
RUN_DATE = "2026-02-15"
SLEEP = 0.34

ROOT = Path(__file__).resolve().parent.parent
SUPP = ROOT / "supplement"
SUPP.mkdir(parents=True, exist_ok=True)

EV_CLAUSE = "(extracellular vesicle OR exosome OR exosomes)"

CONTEXT_QUERY_TIERS = {
    "strict": {
        "neurodegeneration": "(neurodegeneration OR neurodegenerative OR Alzheimer OR Parkinson OR \"amyotrophic lateral sclerosis\" OR ALS OR frontotemporal dementia OR Huntington)",
        "tumor_metastasis": "((cancer OR tumor OR tumour) AND (metastasis OR metastatic OR premetastatic OR organotropism OR dissemination))",
        "cardiac_repair": "((myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury) AND (repair OR regeneration OR cardioprotection OR post-infarction))",
    },
    "relaxed": {
        "neurodegeneration": "(neurodegeneration OR neurodegenerative OR Alzheimer OR Parkinson OR \"amyotrophic lateral sclerosis\" OR ALS OR dementia OR neuroinflammation OR neurotoxicity)",
        "tumor_metastasis": "(cancer OR tumor OR tumour OR neoplasm OR malignancy)",
        "cardiac_repair": "(myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury)",
    },
    "broadest": {
        "neurodegeneration": "(neurologic OR neurological OR \"nervous system\" OR brain OR CNS)",
        "tumor_metastasis": "(cancer OR tumor OR tumour OR neoplasm OR malignancy OR carcinoma OR oncology)",
        "cardiac_repair": "(cardiac OR heart OR myocardial OR cardiovascular OR cardiomyopathy)",
    },
}

PATHWAYS = {
    "mTOR": "(mTOR OR PI3K OR AKT)",
    "Wnt": "(Wnt OR beta-catenin OR CTNNB1)",
    "Notch": "(Notch OR NOTCH1)",
    "NF-kB": "(\"NF-kappa B\" OR NF-kB OR NFKB1)",
    "Complement": "(complement OR C3 OR C5)",
    "Sphingolipid/Ceramide": "(sphingolipid OR ceramide OR sphingomyelinase OR S1P)",
    "Autophagy": "(autophagy OR mitophagy OR lysosome)",
    "Hypoxia/HIF-1": "(hypoxia OR HIF-1 OR HIF1A)",
    "Integrin/Src": "(integrin OR SRC)",
    "TGF-beta": "(\"TGF-beta\" OR TGFB1)",
}

INTEGRIN_TERM = "(integrin)"
SRC_TERM = "(SRC)"


def _req(endpoint: str, params: Dict[str, str]) -> Dict:
    r = requests.get(f"{BASE}/{endpoint}", params=params, timeout=45)
    r.raise_for_status()
    return r.json()


def esearch_count(term: str) -> int:
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "mindate": MINDATE,
        "maxdate": MAXDATE,
        "datetype": "pdat",
        "retmax": "0",
    }
    d = _req("esearch.fcgi", params)
    return int(d["esearchresult"].get("count", "0"))


def compute_priority_from_counts(
    counts: Dict[str, int], smoothing: float, maturity_threshold: int
) -> Dict[str, float]:
    max_count = max(counts.values())
    min_count = min(counts.values())
    ratio = (max_count + smoothing) / (min_count + smoothing)
    ge = sum(1 for v in counts.values() if v >= maturity_threshold)
    return {
        "max_to_min_ratio": ratio,
        "contexts_ge_threshold": ge,
        "priority_index": ratio * (ge / 3.0),
    }


def build_matched_stringency_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows: List[Dict] = []
    for tier, tier_cfg in CONTEXT_QUERY_TIERS.items():
        for context, context_clause in tier_cfg.items():
            for pathway, pathway_clause in PATHWAYS.items():
                term = f"{EV_CLAUSE} AND {context_clause} AND {pathway_clause}"
                count = esearch_count(term)
                rows.append(
                    {
                        "run_date": RUN_DATE,
                        "mindate": MINDATE,
                        "maxdate": MAXDATE,
                        "tier": tier,
                        "context": context,
                        "pathway": pathway,
                        "count": count,
                        "query": term,
                    }
                )
                time.sleep(SLEEP)

    counts_df = pd.DataFrame(rows)
    counts_df.to_csv(SUPP / "sensitivity_matched_query_counts_by_tier.csv", index=False)

    priority_rows: List[Dict] = []
    for tier, tier_df in counts_df.groupby("tier"):
        mat = tier_df.pivot(index="pathway", columns="context", values="count")
        for pathway, row in mat.iterrows():
            counts = {
                "neurodegeneration": int(row["neurodegeneration"]),
                "tumor_metastasis": int(row["tumor_metastasis"]),
                "cardiac_repair": int(row["cardiac_repair"]),
            }
            met = compute_priority_from_counts(counts, smoothing=1.0, maturity_threshold=100)
            priority_rows.append(
                {
                    "tier": tier,
                    "pathway": pathway,
                    "neurodegeneration_count": counts["neurodegeneration"],
                    "tumor_metastasis_count": counts["tumor_metastasis"],
                    "cardiac_repair_count": counts["cardiac_repair"],
                    "max_to_min_ratio": round(met["max_to_min_ratio"], 4),
                    "contexts_ge100": int(met["contexts_ge_threshold"]),
                    "priority_index": round(met["priority_index"], 4),
                }
            )

    priority_df = pd.DataFrame(priority_rows)
    priority_df["rank"] = (
        priority_df.groupby("tier")["priority_index"]
        .rank(method="first", ascending=False)
        .astype(int)
    )
    priority_df = priority_df.sort_values(["tier", "rank", "pathway"])
    priority_df.to_csv(SUPP / "sensitivity_matched_query_priority_ranks.csv", index=False)

    rank_pivot = priority_df.pivot(index="pathway", columns="tier", values="rank")
    rank_corr = rank_pivot.corr(method="spearman")

    stability_rows: List[Dict] = []
    for tier in ["strict", "relaxed", "broadest"]:
        tdf = priority_df[priority_df["tier"] == tier].sort_values("rank")
        top3 = tdf.head(3)["pathway"].tolist()
        stability_rows.append(
            {
                "tier": tier,
                "top1": top3[0],
                "top2": top3[1],
                "top3": top3[2],
                "spearman_vs_strict": round(float(rank_corr.loc["strict", tier]), 4),
            }
        )
    stability_df = pd.DataFrame(stability_rows)
    stability_df.to_csv(SUPP / "sensitivity_matched_query_rank_stability.csv", index=False)

    return counts_df, priority_df, stability_df


def build_index_stability_tables(strict_priority_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = strict_priority_df.copy()

    alt_rows: List[Dict] = []
    for _, r in base.iterrows():
        counts = {
            "neurodegeneration": int(r["neurodegeneration_count"]),
            "tumor_metastasis": int(r["tumor_metastasis_count"]),
            "cardiac_repair": int(r["cardiac_repair_count"]),
        }
        max_count = max(counts.values())
        min_count = min(counts.values())
        total_count = sum(counts.values())
        pathway = r["pathway"]

        out: Dict[str, float | str] = {
            "pathway": pathway,
            "max_count": max_count,
            "min_count": min_count,
            "total_count": total_count,
        }

        for s in (0.5, 1.0, 2.0):
            ratio = (max_count + s) / (min_count + s)
            ge100 = sum(1 for v in counts.values() if v >= 100)
            out[f"ratio_only_s{s}"] = ratio
            out[f"ratio_x_logtotal_s{s}"] = ratio * math.log1p(total_count)
            out[f"ratio_x_donor_s{s}"] = ratio * math.log1p(max_count)
            out[f"ratio_x_maturity100_s{s}"] = ratio * (ge100 / 3.0)
        alt_rows.append(out)

    alt_df = pd.DataFrame(alt_rows)

    score_cols = [c for c in alt_df.columns if c.startswith("ratio_")]
    for c in score_cols:
        alt_df[f"rank__{c}"] = alt_df[c].rank(method="first", ascending=False).astype(int)
    alt_df = alt_df.sort_values("rank__ratio_x_maturity100_s1.0")
    alt_df.to_csv(SUPP / "priority_index_alternative_formulations.csv", index=False)

    baseline_rank = alt_df.set_index("pathway")["rank__ratio_x_maturity100_s1.0"]

    summary_rows: List[Dict] = []
    for c in score_cols:
        rank_col = f"rank__{c}"
        rk = alt_df.set_index("pathway")[rank_col]
        top3 = alt_df.sort_values(rank_col).head(3)["pathway"].tolist()
        summary_rows.append(
            {
                "formulation": c,
                "spearman_vs_baseline": round(float(rk.corr(baseline_rank, method="spearman")), 4),
                "top1": top3[0],
                "top2": top3[1],
                "top3": top3[2],
            }
        )
    summary_df = pd.DataFrame(summary_rows).sort_values("spearman_vs_baseline", ascending=False)
    summary_df.to_csv(SUPP / "priority_index_rank_stability_summary.csv", index=False)

    return alt_df, summary_df


def build_integrin_src_partition() -> pd.DataFrame:
    rows: List[Dict] = []
    strict_contexts = CONTEXT_QUERY_TIERS["strict"]
    for context, context_clause in strict_contexts.items():
        q_combined = f"{EV_CLAUSE} AND {context_clause} AND ({INTEGRIN_TERM} OR {SRC_TERM})"
        q_integrin_only = f"{EV_CLAUSE} AND {context_clause} AND ({INTEGRIN_TERM} NOT {SRC_TERM})"
        q_src_only = f"{EV_CLAUSE} AND {context_clause} AND ({SRC_TERM} NOT {INTEGRIN_TERM})"
        q_both = f"{EV_CLAUSE} AND {context_clause} AND ({INTEGRIN_TERM} AND {SRC_TERM})"

        combined = esearch_count(q_combined)
        time.sleep(SLEEP)
        integrin_only = esearch_count(q_integrin_only)
        time.sleep(SLEEP)
        src_only = esearch_count(q_src_only)
        time.sleep(SLEEP)
        both = esearch_count(q_both)
        time.sleep(SLEEP)

        rows.append(
            {
                "run_date": RUN_DATE,
                "context": context,
                "combined_integrin_or_src": combined,
                "integrin_only_not_src": integrin_only,
                "src_only_not_integrin": src_only,
                "integrin_and_src": both,
                "integrin_only_pct_of_combined": round((integrin_only / combined) * 100, 2) if combined else 0.0,
                "src_only_pct_of_combined": round((src_only / combined) * 100, 2) if combined else 0.0,
                "both_pct_of_combined": round((both / combined) * 100, 2) if combined else 0.0,
                "query_combined": q_combined,
                "query_integrin_only": q_integrin_only,
                "query_src_only": q_src_only,
                "query_both": q_both,
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(SUPP / "integrin_src_partition_counts.csv", index=False)
    return df


def main() -> None:
    counts_df, priority_df, matched_stability_df = build_matched_stringency_tables()
    strict_priority = priority_df[priority_df["tier"] == "strict"].copy()
    _, index_stability_df = build_index_stability_tables(strict_priority)
    partition_df = build_integrin_src_partition()

    manifest = {
        "run_date": RUN_DATE,
        "window": {"mindate": MINDATE, "maxdate": MAXDATE},
        "outputs": [
            "supplement/sensitivity_matched_query_counts_by_tier.csv",
            "supplement/sensitivity_matched_query_priority_ranks.csv",
            "supplement/sensitivity_matched_query_rank_stability.csv",
            "supplement/priority_index_alternative_formulations.csv",
            "supplement/priority_index_rank_stability_summary.csv",
            "supplement/integrin_src_partition_counts.csv",
        ],
        "matched_query_tiers": list(CONTEXT_QUERY_TIERS.keys()),
        "index_formulations_summary_rows": int(len(index_stability_df)),
        "integrin_src_contexts": int(len(partition_df)),
    }
    (SUPP / "robustness_hardening_manifest.json").write_text(json.dumps(manifest, indent=2))

    print("Saved robustness artifacts:")
    for p in manifest["outputs"]:
        print("-", p)
    print("- supplement/robustness_hardening_manifest.json")
    print("\nMatched-query top3 by tier:")
    print(matched_stability_df.to_string(index=False))
    print("\nIndex stability (top rows):")
    print(index_stability_df.head(8).to_string(index=False))
    print("\nIntegrin/Src partition:")
    print(
        partition_df[
            [
                "context",
                "combined_integrin_or_src",
                "integrin_only_not_src",
                "src_only_not_integrin",
                "integrin_and_src",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
