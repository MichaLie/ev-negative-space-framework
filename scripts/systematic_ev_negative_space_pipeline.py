#!/usr/bin/env python3
"""Systematic EV negative-space mining pipeline.

Builds reproducible evidence tables for:
- EV cargo class counts across contexts
- EV pathway x context evidence density
- PMID-level screening set with publication type + relevance scoring
- Negative-space prioritization index and roadmap
- Figures for manuscript and supplement

Run date anchored to analysis day unless overridden.
"""

from __future__ import annotations

import json
import math
import re
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MINDATE = "2015/01/01"
MAXDATE = "2026/02/15"
RUN_DATE = "2026-02-15"
RETMAX_SCREEN = 40
SLEEP = 0.34

ROOT = Path(__file__).resolve().parent.parent  # repo root: Contro_CC/
SUPP = ROOT / "supplement"
FIG = ROOT
SUPP.mkdir(parents=True, exist_ok=True)

EV_CLAUSE = "(extracellular vesicle OR exosome OR exosomes)"

CONTEXTS = {
    "neurodegeneration": {
        "query": "(neurodegeneration OR neurodegenerative OR Alzheimer OR Parkinson OR \"amyotrophic lateral sclerosis\" OR ALS OR frontotemporal dementia OR Huntington)",
        "keywords": [
            "neurodegeneration",
            "neurodegenerative",
            "alzheimer",
            "parkinson",
            "als",
            "dementia",
            "neuronal",
            "brain",
        ],
    },
    "tumor_metastasis": {
        "query": "((cancer OR tumor OR tumour) AND (metastasis OR metastatic OR premetastatic OR organotropism OR dissemination))",
        "keywords": [
            "cancer",
            "tumor",
            "tumour",
            "metastasis",
            "metastatic",
            "pre-metastatic",
            "premetastatic",
            "organotrop",
        ],
    },
    "cardiac_repair": {
        "query": "((myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury) AND (repair OR regeneration OR cardioprotection OR post-infarction))",
        "keywords": [
            "cardiac",
            "heart",
            "myocardial",
            "infarction",
            "cardiomyocyte",
            "angiogenesis",
            "fibrosis",
        ],
    },
}

CARGO = {
    "proteins/proteomics": {
        "query": "(protein OR proteomics)",
        "keywords": ["protein", "proteomic", "proteome"],
    },
    "miRNA": {
        "query": "(miRNA OR microRNA)",
        "keywords": ["mir", "microrna"],
    },
    "lipids/lipidomics": {
        "query": "(lipid OR lipidomics OR sphingolipid OR ceramide)",
        "keywords": ["lipid", "lipidomic", "ceramide", "sphingolipid", "s1p"],
    },
}

PATHWAYS = {
    "mTOR": {
        "query": "(mTOR OR PI3K OR AKT)",
        "keywords": ["mtor", "pi3k", "akt"],
    },
    "Wnt": {
        "query": "(Wnt OR beta-catenin OR CTNNB1)",
        "keywords": ["wnt", "beta-catenin", "ctnnb1"],
    },
    "Notch": {
        "query": "(Notch OR NOTCH1)",
        "keywords": ["notch", "notch1"],
    },
    "NF-kB": {
        "query": "(\"NF-kappa B\" OR NF-kB OR NFKB1)",
        "keywords": ["nf-kappa", "nf-kb", "nfkb"],
    },
    "Complement": {
        "query": "(complement OR C3 OR C5)",
        "keywords": ["complement", "c3", "c5", "c3ar"],
    },
    "Sphingolipid/Ceramide": {
        "query": "(sphingolipid OR ceramide OR sphingomyelinase OR S1P)",
        "keywords": ["sphingolipid", "ceramide", "sphingomyelin", "s1p"],
    },
    "Autophagy": {
        "query": "(autophagy OR mitophagy OR lysosome)",
        "keywords": ["autophagy", "mitophagy", "lysosome", "autophagic"],
    },
    "Hypoxia/HIF-1": {
        "query": "(hypoxia OR HIF-1 OR HIF1A)",
        "keywords": ["hypoxia", "hif-1", "hif1a"],
    },
    "Integrin/Src": {
        "query": "(integrin OR SRC)",
        "keywords": ["integrin", "itg", "src"],
    },
    "TGF-beta": {
        "query": "(\"TGF-beta\" OR TGFB1)",
        "keywords": ["tgf", "tgfb"],
    },
}

MECH_KEYWORDS = [
    "promot",
    "inhibit",
    "mediate",
    "transfer",
    "deliver",
    "uptake",
    "regulat",
    "activate",
    "suppress",
    "therapeut",
    "metast",
    "repair",
]

PRIMARY_EXCLUDE_TYPES = {
    "Review",
    "Systematic Review",
    "Meta-Analysis",
    "Editorial",
    "Comment",
    "Letter",
    "News",
}

ROADMAP_NOTES = {
    "TGF-beta": {
        "translational_hypothesis": "Leverage metastasis-mature EV TGF-beta control to tune post-MI fibrosis and scar quality while preserving required wound-healing signaling.",
        "safety_gate": "Gate pro-fibrotic drift, EMT signatures, and impaired scar integrity.",
        "starter_models": "Cardiac fibroblast activation assays; mouse MI with scar histology, collagen organization, and ventricular function readouts.",
    },
    "Integrin/Src": {
        "translational_hypothesis": "Repurpose tumor organotropism logic (integrin-coded tropism and integrin-Src signaling) to improve EV homing and retention in injured myocardium.",
        "safety_gate": "Gate pro-migratory/invasive signaling, aberrant angiogenesis, and oncogenic cargo transfer.",
        "starter_models": "iPSC-cardiomyocyte/fibroblast/endothelial co-cultures plus mouse MI biodistribution and fibrosis/arrhythmia monitoring.",
    },
    "Wnt": {
        "translational_hypothesis": "Translate controlled Wnt/beta-catenin EV signaling to improve angiogenesis and cardiomyocyte survival after myocardial injury.",
        "safety_gate": "Gate uncontrolled proliferation, arrhythmogenic remodeling, and fibrosis.",
        "starter_models": "Wnt-reporter cardiac/endothelial assays and MI or ischemia-reperfusion models with ECG and vascular/fibrosis endpoints.",
    },
    "Hypoxia/HIF-1": {
        "translational_hypothesis": "Transfer hypoxia/HIF-1 EV programming from metastasis to neurodegeneration for neurovascular protection and ferroptosis resistance.",
        "safety_gate": "Pro-angiogenic mispatterning and hypoxia-driven oncogenic response screening.",
        "starter_models": "Neuron-glia oxidative-stress and ferroptosis models plus AD/PD in vivo systems with BBB-integrity and behavior readouts.",
    },
    "NF-kB": {
        "translational_hypothesis": "Translate EV-mediated NF-kB tuning to optimize post-MI inflammatory resolution and repair timing.",
        "safety_gate": "Gate immune suppression, delayed healing, and off-target systemic inflammation.",
        "starter_models": "Macrophage polarization and endothelial activation assays plus MI immunophenotyping time-course in vivo.",
    },
    "mTOR": {
        "translational_hypothesis": "Evaluate whether EV-mediated PI3K-AKT-mTOR signaling can be harnessed for cardioprotection without inducing maladaptive hypertrophy.",
        "safety_gate": "Gate pathological hypertrophy, metabolic derangement, and pro-growth overactivation.",
        "starter_models": "iPSC-cardiomyocyte stress and hypertrophy assays plus MI/ischemia-reperfusion long-term remodeling follow-up.",
    },
    "Notch": {
        "translational_hypothesis": "Assess EV-mediated Notch/Jagged signaling to enhance reparative angiogenesis and cardiomyocyte survival after injury.",
        "safety_gate": "Gate aberrant vascular patterning and fibrotic remodeling.",
        "starter_models": "Endothelial Notch reporter assays with in vivo MI vascular-density and fibrosis endpoints.",
    },
    "Complement": {
        "translational_hypothesis": "Test whether EV-associated complement modulation can reduce reperfusion injury in cardiac repair without compromising host defense.",
        "safety_gate": "Gate infection susceptibility, off-target complement depletion, and systemic inflammatory liabilities.",
        "starter_models": "Cardiomyocyte/endothelial injury models with complement readouts plus in vivo ischemia-reperfusion infarct and biomarker tracking.",
    },
    "Sphingolipid/Ceramide": {
        "translational_hypothesis": "Explore EV sphingolipid/ceramide programs to improve cardiomyocyte stress resistance and apoptosis control in cardiac repair.",
        "safety_gate": "Gate pro-apoptotic ceramide accumulation, inflammatory lipid signaling, and systemic lipid toxicity.",
        "starter_models": "Cardiomyocyte ceramide flux and apoptosis assays plus in vivo ischemia-reperfusion cell-death and lipidomics checks.",
    },
    "Autophagy": {
        "translational_hypothesis": "Transfer EV autophagy/mitophagy modulation strategies from neurodegeneration to improve cardiomyocyte stress tolerance after ischemia-reperfusion.",
        "safety_gate": "Excess autophagic flux and mitochondrial dysfunction.",
        "starter_models": "Hypoxia/reoxygenation cardiomyocyte assays and MI/ischemia-reperfusion models with LC3/p62 flux and function endpoints.",
    },
}


def _req(endpoint: str, params: Dict[str, str], expect_json: bool = True):
    url = f"{BASE}/{endpoint}"
    r = requests.get(url, params=params, timeout=40)
    r.raise_for_status()
    return r.json() if expect_json else r.text


def esearch(term: str, retmax: int = 0, sort: str = "relevance") -> Tuple[int, List[str]]:
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "mindate": MINDATE,
        "maxdate": MAXDATE,
        "datetype": "pdat",
        "retmax": str(retmax),
        "sort": sort,
    }
    data = _req("esearch.fcgi", params)
    d = data["esearchresult"]
    return int(d.get("count", "0")), d.get("idlist", [])


def esummary(pmids: List[str]) -> Dict[str, Dict[str, str]]:
    if not pmids:
        return {}
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "json"}
    data = _req("esummary.fcgi", params)
    result = data.get("result", {})
    out = {}
    for pmid in pmids:
        rec = result.get(pmid, {})
        doi = ""
        for aid in rec.get("articleids", []):
            if aid.get("idtype") == "doi":
                doi = aid.get("value", "")
                break
        out[pmid] = {
            "title": rec.get("title", ""),
            "pubdate": rec.get("pubdate", ""),
            "journal": rec.get("fulljournalname", ""),
            "doi": doi,
        }
    return out


def efetch_details(pmids: List[str]) -> Dict[str, Dict[str, str]]:
    if not pmids:
        return {}
    params = {"db": "pubmed", "id": ",".join(pmids), "retmode": "xml"}
    xml_text = _req("efetch.fcgi", params, expect_json=False)
    root = ET.fromstring(xml_text)
    out = {}

    for article in root.findall(".//PubmedArticle"):
        pmid_node = article.find(".//PMID")
        if pmid_node is None or not pmid_node.text:
            continue
        pmid = pmid_node.text.strip()

        abs_parts = []
        for abst in article.findall(".//Abstract/AbstractText"):
            txt = (abst.text or "").strip()
            if txt:
                abs_parts.append(txt)
        abstract = " ".join(abs_parts)

        pub_types = []
        for ptype in article.findall(".//PublicationTypeList/PublicationType"):
            txt = (ptype.text or "").strip()
            if txt:
                pub_types.append(txt)

        out[pmid] = {
            "abstract": abstract,
            "publication_types": "; ".join(sorted(set(pub_types))),
        }
    return out


def compute_relevance(text: str, context_kw: List[str], pathway_kw: List[str]) -> Tuple[int, bool, bool, bool, bool]:
    t = text.lower()
    ev_hit = "exosome" in t or "extracellular vesicle" in t or bool(re.search(r'\bev\b', t))
    context_hit = any(re.search(rf'\b{re.escape(k)}\b', t) for k in context_kw)
    pathway_hit = any(re.search(rf'\b{re.escape(k)}\b', t) for k in pathway_kw)
    mech_hit = any(k in t for k in MECH_KEYWORDS)
    score = int(ev_hit) + int(context_hit) + int(pathway_hit) + int(mech_hit)
    return score, ev_hit, context_hit, pathway_hit, mech_hit


def is_primary(publication_types: str) -> bool:
    if not publication_types:
        return True
    pts = {p.strip() for p in publication_types.split(";") if p.strip()}
    return len(pts.intersection(PRIMARY_EXCLUDE_TYPES)) == 0


def attention_level(count: int) -> str:
    if count >= 500:
        return "High"
    if count >= 100:
        return "Moderate"
    return "Low"


def mine_counts_and_screening():
    # Cargo x context counts
    cargo_rows = []
    query_registry = []
    for cargo_name, cargo_cfg in CARGO.items():
        for ctx_name, ctx_cfg in CONTEXTS.items():
            term = f"{EV_CLAUSE} AND {ctx_cfg['query']} AND {cargo_cfg['query']}"
            count, _ = esearch(term, retmax=0)
            cargo_rows.append(
                {
                    "run_date": RUN_DATE,
                    "mindate": MINDATE,
                    "maxdate": MAXDATE,
                    "context": ctx_name,
                    "cargo_axis": cargo_name,
                    "count": count,
                    "query": term,
                }
            )
            query_registry.append(
                {
                    "query_group": "cargo_context",
                    "label_1": cargo_name,
                    "label_2": ctx_name,
                    "query": term,
                }
            )
            time.sleep(SLEEP)

    cargo_df = pd.DataFrame(cargo_rows)
    cargo_df.to_csv(ROOT / "pubmed_cargo_context_counts_strict_2015_2026.csv", index=False)

    # Pathway x context counts + screening pulls
    path_rows = []
    screening_rows = []
    seen_pmids = set()

    for pathway, p_cfg in PATHWAYS.items():
        for ctx_name, ctx_cfg in CONTEXTS.items():
            term = f"{EV_CLAUSE} AND {ctx_cfg['query']} AND {p_cfg['query']}"
            count, pmids = esearch(term, retmax=RETMAX_SCREEN)
            path_rows.append(
                {
                    "run_date": RUN_DATE,
                    "mindate": MINDATE,
                    "maxdate": MAXDATE,
                    "pathway": pathway,
                    "context": ctx_name,
                    "count": count,
                    "query": term,
                }
            )
            query_registry.append(
                {
                    "query_group": "pathway_context",
                    "label_1": pathway,
                    "label_2": ctx_name,
                    "query": term,
                }
            )

            if pmids:
                meta = esummary(pmids)
                details = efetch_details(pmids)
                for rank, pmid in enumerate(pmids, start=1):
                    m = meta.get(pmid, {})
                    d = details.get(pmid, {})
                    title = m.get("title", "")
                    abstract = d.get("abstract", "")
                    text_blob = f"{title} {abstract}"
                    score, ev_hit, context_hit, pathway_hit, mech_hit = compute_relevance(
                        text_blob,
                        context_kw=ctx_cfg["keywords"],
                        pathway_kw=p_cfg["keywords"],
                    )
                    primary_flag = is_primary(d.get("publication_types", ""))
                    if score >= 3 and primary_flag and context_hit and pathway_hit:
                        decision = "include_primary_high"
                    elif score >= 2 and context_hit:
                        decision = "include_contextual"
                    else:
                        decision = "exclude_low_relevance"

                    screening_rows.append(
                        {
                            "pathway": pathway,
                            "context": ctx_name,
                            "rank": rank,
                            "total_hits_for_cell": count,
                            "pmid": pmid,
                            "title": title,
                            "pubdate": m.get("pubdate", ""),
                            "journal": m.get("journal", ""),
                            "doi": m.get("doi", ""),
                            "publication_types": d.get("publication_types", ""),
                            "is_primary": primary_flag,
                            "relevance_score": score,
                            "ev_hit": ev_hit,
                            "context_hit": context_hit,
                            "pathway_hit": pathway_hit,
                            "mechanistic_hit": mech_hit,
                            "screen_decision": decision,
                            "abstract": abstract,
                        }
                    )
                    seen_pmids.add(pmid)
            time.sleep(SLEEP)

    path_df = pd.DataFrame(path_rows)
    path_df.to_csv(ROOT / "pubmed_pathway_context_counts_strict_2015_2026.csv", index=False)

    screening_df = pd.DataFrame(screening_rows)
    screening_df.to_csv(SUPP / "screening_table_pathway_context_top40.csv", index=False)

    qr = pd.DataFrame(query_registry)
    qr.to_csv(SUPP / "query_registry.csv", index=False)

    # Screening summary (PRISMA-like)
    if not screening_df.empty:
        summary = {
            "records_identified_topN": int(len(screening_df)),
            "unique_pmids_screened": int(screening_df["pmid"].nunique()),
            "included_primary_high": int((screening_df["screen_decision"] == "include_primary_high").sum()),
            "included_contextual": int((screening_df["screen_decision"] == "include_contextual").sum()),
            "excluded_low_relevance": int((screening_df["screen_decision"] == "exclude_low_relevance").sum()),
            "run_date": RUN_DATE,
            "date_window": {"mindate": MINDATE, "maxdate": MAXDATE},
            "retmax_per_cell": RETMAX_SCREEN,
        }
    else:
        summary = {}
    (SUPP / "screening_summary.json").write_text(json.dumps(summary, indent=2))

    return cargo_df, path_df, screening_df


def build_priority_and_roadmap(path_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mat = (
        path_df.pivot(index="pathway", columns="context", values="count")
        .fillna(0)
        .astype(int)
        .reset_index()
    )

    rows = []
    for _, r in mat.iterrows():
        counts = {
            "neurodegeneration": int(r.get("neurodegeneration", 0)),
            "tumor_metastasis": int(r.get("tumor_metastasis", 0)),
            "cardiac_repair": int(r.get("cardiac_repair", 0)),
        }
        max_ctx = max(counts, key=counts.get)
        min_ctx = min(counts, key=counts.get)
        max_count = counts[max_ctx]
        min_count = counts[min_ctx]
        ratio = (max_count + 1) / (min_count + 1)
        ge100 = sum(1 for v in counts.values() if v >= 100)
        priority_index = ratio * (ge100 / 3)

        rows.append(
            {
                "pathway": r["pathway"],
                "tumor_metastasis_count": counts["tumor_metastasis"],
                "cardiac_repair_count": counts["cardiac_repair"],
                "neurodegeneration_count": counts["neurodegeneration"],
                "tumor_attention": attention_level(counts["tumor_metastasis"]),
                "cardiac_attention": attention_level(counts["cardiac_repair"]),
                "neuro_attention": attention_level(counts["neurodegeneration"]),
                "max_context": max_ctx,
                "min_context": min_ctx,
                "max_to_min_ratio": round(ratio, 2),
                "contexts_ge100": ge100,
                "negative_space_priority_index": round(priority_index, 2),
            }
        )

    priority_df = pd.DataFrame(rows).sort_values(
        "negative_space_priority_index", ascending=False
    )
    priority_df.to_csv(ROOT / "negative_space_priority_index.csv", index=False)

    roadmap_rows = []
    for _, r in priority_df.iterrows():
        pathway = r["pathway"]
        note = ROADMAP_NOTES.get(pathway, {})
        roadmap_rows.append(
            {
                "rank": len(roadmap_rows) + 1,
                "pathway": pathway,
                "target_underexplored_context": r["min_context"],
                "most_developed_context": r["max_context"],
                "negative_space_priority_index": r["negative_space_priority_index"],
                "translational_hypothesis": note.get("translational_hypothesis", ""),
                "safety_gate": note.get("safety_gate", ""),
                "starter_models": note.get("starter_models", ""),
            }
        )

    roadmap_df = pd.DataFrame(roadmap_rows)
    roadmap_df.to_csv(ROOT / "cross_field_phase1_roadmap.csv", index=False)

    return priority_df, roadmap_df


def mine_year_trends():
    rows = []
    for year in range(2010, 2027):
        for ctx_name, ctx_cfg in CONTEXTS.items():
            params = {
                "db": "pubmed",
                "term": f"{EV_CLAUSE} AND {ctx_cfg['query']}",
                "retmode": "json",
                "rettype": "count",
                "mindate": f"{year}/01/01",
                "maxdate": f"{year}/12/31",
                "datetype": "pdat",
            }
            data = _req("esearch.fcgi", params)
            cnt = int(data["esearchresult"].get("count", "0"))
            rows.append({"year": year, "context": ctx_name, "count": cnt})
            time.sleep(SLEEP)
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "pubmed_context_year_counts_2010_2026_strict.csv", index=False)
    return df


def make_figures(cargo_df: pd.DataFrame, path_df: pd.DataFrame, priority_df: pd.DataFrame, trend_df: pd.DataFrame):
    sns.set_theme(style="whitegrid")

    # Figure 1: pathway heatmap
    mat = path_df.pivot(index="pathway", columns="context", values="count").loc[
        :, ["neurodegeneration", "tumor_metastasis", "cardiac_repair"]
    ]
    plt.figure(figsize=(9, 6))
    sns.heatmap(mat, annot=True, fmt="d", cmap="YlGnBu", linewidths=0.5)
    plt.title("EV Pathway Evidence Density by Context (2015-01-01 to 2026-02-15)")
    plt.xlabel("Context")
    plt.ylabel("Pathway")
    plt.tight_layout()
    plt.savefig(FIG / "figure1_pathway_heatmap_strict.png", dpi=320)
    plt.close()

    # Figure 2: cargo bars
    plt.figure(figsize=(9, 5.6))
    tmp = cargo_df.copy()
    sns.barplot(data=tmp, x="cargo_axis", y="count", hue="context")
    plt.title("EV Cargo-Class Literature Density by Context (2015-2026)")
    plt.xlabel("Cargo axis")
    plt.ylabel("PubMed records")
    plt.xticks(rotation=12)
    plt.tight_layout()
    plt.savefig(FIG / "figure2_cargo_context_bars_strict.png", dpi=320)
    plt.close()

    # Figure 3: priority index
    plt.figure(figsize=(9, 5.8))
    sns.barplot(
        data=priority_df,
        y="pathway",
        x="negative_space_priority_index",
        color="#247ba0",
    )
    plt.title("Negative-Space Priority Index by EV Pathway")
    plt.xlabel("Priority index")
    plt.ylabel("Pathway")
    plt.tight_layout()
    plt.savefig(FIG / "figure3_negative_space_priority.png", dpi=320)
    plt.close()

    # Figure 4: trend lines
    plt.figure(figsize=(9, 5.6))
    sns.lineplot(data=trend_df, x="year", y="count", hue="context", marker="o")
    plt.title("Annual EV Literature Growth by Context (2010-2026)")
    plt.xlabel("Year")
    plt.ylabel("PubMed records")
    plt.tight_layout()
    plt.savefig(FIG / "figure4_annual_growth_contexts_strict.png", dpi=320)
    plt.close()

    # Supplementary PRISMA-like flow diagram from screening summary
    summary = json.loads((SUPP / "screening_summary.json").read_text())
    if summary:
        fig, ax = plt.subplots(figsize=(8.2, 6.6))
        fig.patch.set_facecolor("white")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

        def box(text, x, y, w=0.8, h=0.12):
            rect = plt.Rectangle(
                (x - w / 2, y - h / 2),
                w,
                h,
                fill=True,
                facecolor="#f8fafc",
                edgecolor="#1f2937",
                linewidth=2.0,
            )
            ax.add_patch(rect)
            ax.text(x, y, text, ha="center", va="center", fontsize=10, color="#111827")

        box(f"Records identified (top-N across cells)\nn={summary['records_identified_topN']}", 0.5, 0.88)
        box(f"Unique PMIDs screened\nn={summary['unique_pmids_screened']}", 0.5, 0.68)
        box(f"Included: primary high relevance\nn={summary['included_primary_high']}", 0.32, 0.44, w=0.36)
        box(f"Included: contextual\nn={summary['included_contextual']}", 0.68, 0.44, w=0.36)
        box(f"Excluded: low relevance\nn={summary['excluded_low_relevance']}", 0.5, 0.22)

        arrowprops = dict(arrowstyle="->", linewidth=1.8, color="#1f2937")
        ax.annotate("", xy=(0.5, 0.74), xytext=(0.5, 0.82), arrowprops=arrowprops)
        ax.annotate("", xy=(0.32, 0.50), xytext=(0.47, 0.62), arrowprops=arrowprops)
        ax.annotate("", xy=(0.68, 0.50), xytext=(0.53, 0.62), arrowprops=arrowprops)
        ax.annotate("", xy=(0.5, 0.28), xytext=(0.5, 0.38), arrowprops=arrowprops)

        plt.title("Supplementary Figure S1. Screening Workflow Summary", fontsize=12, pad=18, color="#111827")
        plt.tight_layout()
        plt.savefig(SUPP / "figureS1_screening_flow.png", dpi=320)
        plt.close()


def main():
    cargo_df, path_df, screening_df = mine_counts_and_screening()
    priority_df, roadmap_df = build_priority_and_roadmap(path_df)
    trend_df = mine_year_trends()
    make_figures(cargo_df, path_df, priority_df, trend_df)

    method_manifest = {
        "run_date": RUN_DATE,
        "window": {"mindate": MINDATE, "maxdate": MAXDATE},
        "ev_clause": EV_CLAUSE,
        "contexts": CONTEXTS,
        "cargo_axes": CARGO,
        "pathways": PATHWAYS,
        "retmax_screen": RETMAX_SCREEN,
        "outputs": [
            "pubmed_cargo_context_counts_strict_2015_2026.csv",
            "pubmed_pathway_context_counts_strict_2015_2026.csv",
            "supplement/screening_table_pathway_context_top40.csv",
            "negative_space_priority_index.csv",
            "cross_field_phase1_roadmap.csv",
            "pubmed_context_year_counts_2010_2026_strict.csv",
        ],
    }
    (SUPP / "methods_manifest.json").write_text(json.dumps(method_manifest, indent=2))

    print("Pipeline complete.")
    print(f"Pathway-context rows: {len(path_df)}")
    print(f"Screening rows: {len(screening_df)}")
    print(f"Roadmap rows: {len(roadmap_df)}")


if __name__ == "__main__":
    main()
