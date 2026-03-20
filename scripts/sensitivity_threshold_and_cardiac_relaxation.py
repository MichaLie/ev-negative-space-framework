#!/usr/bin/env python3
"""Generate threshold-sweep, continuous-variant, and cardiac-only relaxation
sensitivity analyses (Tables S11 and S12).

B) Maturity-threshold sweep (50, 75, 100, 150) and continuous variant
   -> sensitivity_priority_threshold.csv  (Table S12)
   Uses already-collected strict counts; no API calls required.

D) Cardiac-only relaxation (relax cardiac query while holding neuro/tumor strict)
   -> sensitivity_cardiac_query.csv  (Table S11)
   Requires PubMed E-utilities API calls for relaxed/broadest cardiac queries
   while neuro and tumor contexts remain at strict stringency.

Usage:
    python sensitivity_threshold_and_cardiac_relaxation.py           # threshold sweep only (offline)
    python sensitivity_threshold_and_cardiac_relaxation.py --api     # threshold sweep + cardiac relaxation (API)
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Dict

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
SUPP = ROOT / "supplement"

# PubMed E-utilities configuration
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MINDATE = "2015/01/01"
MAXDATE = "2026/02/15"
RUN_DATE = "2026-02-15"
SLEEP = 0.34

EV_CLAUSE = "(extracellular vesicle OR exosome OR exosomes)"

# Context queries: neuro and tumor are ALWAYS strict for cardiac-only relaxation
NEURO_STRICT = (
    "(neurodegeneration OR neurodegenerative OR Alzheimer OR Parkinson "
    'OR "amyotrophic lateral sclerosis" OR ALS OR frontotemporal dementia OR Huntington)'
)
TUMOR_STRICT = (
    "((cancer OR tumor OR tumour) AND "
    "(metastasis OR metastatic OR premetastatic OR organotropism OR dissemination))"
)

# Cardiac query tiers for relaxation analysis
CARDIAC_TIERS = {
    "strict": (
        "((myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury) "
        "AND (repair OR regeneration OR cardioprotection OR post-infarction))"
    ),
    "relaxed": (
        "(myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury)"
    ),
    "broadest": (
        "(cardiac OR heart OR myocardial OR cardiovascular OR cardiomyopathy)"
    ),
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


# ---------------------------------------------------------------------------
# PubMed helper
# ---------------------------------------------------------------------------

def esearch_count(term: str) -> int:
    """Query PubMed E-utilities and return result count."""
    params = {
        "db": "pubmed",
        "term": term,
        "retmode": "json",
        "mindate": MINDATE,
        "maxdate": MAXDATE,
        "datetype": "pdat",
        "retmax": "0",
    }
    r = requests.get(f"{BASE}/esearch.fcgi", params=params, timeout=45)
    r.raise_for_status()
    return int(r.json()["esearchresult"].get("count", "0"))


# ---------------------------------------------------------------------------
# Shared helper: priority-index formula (matches STAR Methods exactly)
#   priority_index = ((max + s) / (min + s)) * (contexts_ge_threshold / 3)
# ---------------------------------------------------------------------------

def compute_priority(
    counts: Dict[str, int],
    smoothing: float = 1.0,
    maturity_threshold: int = 100,
) -> float:
    """Threshold-based priority index as described in STAR Methods."""
    max_count = max(counts.values())
    min_count = min(counts.values())
    ratio = (max_count + smoothing) / (min_count + smoothing)
    ge = sum(1 for v in counts.values() if v >= maturity_threshold)
    return ratio * (ge / 3.0)


def compute_continuous(
    counts: Dict[str, int],
    smoothing: float = 1.0,
) -> float:
    """Continuous variant: raw smoothed ratio without hard maturity cutoff."""
    max_count = max(counts.values())
    min_count = min(counts.values())
    return (max_count + smoothing) / (min_count + smoothing)


# ---------------------------------------------------------------------------
# B: Threshold sweep + continuous variant  -> Table S12
# ---------------------------------------------------------------------------

def build_threshold_sensitivity() -> pd.DataFrame:
    """Vary maturity threshold from 50 to 150 and compute continuous variant.

    Reads strict pathway-context counts from the main pipeline output
    (pubmed_pathway_context_counts_strict_2015_2026.csv). No API calls needed.
    """
    counts_df = pd.read_csv(ROOT / "pubmed_pathway_context_counts_strict_2015_2026.csv")
    mat = counts_df.pivot(index="pathway", columns="context", values="count")

    rows = []
    for pathway, row in mat.iterrows():
        counts = {
            "neurodegeneration": int(row["neurodegeneration"]),
            "tumor_metastasis": int(row["tumor_metastasis"]),
            "cardiac_repair": int(row["cardiac_repair"]),
        }
        entry = {"pathway": pathway}
        for threshold in [50, 75, 100, 150]:
            pi = compute_priority(counts, smoothing=1.0, maturity_threshold=threshold)
            entry[f"priority_threshold_{threshold}"] = round(pi, 2)
        entry["priority_continuous"] = round(compute_continuous(counts, smoothing=1.0), 2)
        rows.append(entry)

    df = pd.DataFrame(rows).sort_values("priority_threshold_100", ascending=False)
    df.to_csv(ROOT / "sensitivity_priority_threshold.csv", index=False)
    print("Saved: sensitivity_priority_threshold.csv")
    print(df.to_string(index=False))
    return df


# ---------------------------------------------------------------------------
# D: Cardiac-only relaxation  -> Table S11
# ---------------------------------------------------------------------------

def build_cardiac_relaxation() -> pd.DataFrame:
    """Relax only the cardiac query tier while holding neuro and tumor at strict.

    This isolates the effect of cardiac query stringency on evidence counts,
    independent of changes to the other two contexts.
    Requires live PubMed API access (--api flag).
    """
    # Strict counts from main pipeline (already collected)
    strict_df = pd.read_csv(ROOT / "pubmed_pathway_context_counts_strict_2015_2026.csv")
    strict_mat = strict_df.pivot(index="pathway", columns="context", values="count")

    rows = []
    for pathway, pathway_clause in PATHWAYS.items():
        neuro = int(strict_mat.loc[pathway, "neurodegeneration"])
        tumor = int(strict_mat.loc[pathway, "tumor_metastasis"])
        cardiac_strict = int(strict_mat.loc[pathway, "cardiac_repair"])

        # Query cardiac under relaxed and broadest tiers
        cardiac_relaxed = esearch_count(
            f"{EV_CLAUSE} AND {CARDIAC_TIERS['relaxed']} AND {pathway_clause}"
        )
        time.sleep(SLEEP)

        cardiac_broadest = esearch_count(
            f"{EV_CLAUSE} AND {CARDIAC_TIERS['broadest']} AND {pathway_clause}"
        )
        time.sleep(SLEEP)

        rows.append({
            "pathway": pathway,
            "neuro_strict": neuro,
            "tumor_strict": tumor,
            "cardiac_strict": cardiac_strict,
            "cardiac_relaxed": cardiac_relaxed,
            "cardiac_broadest": cardiac_broadest,
        })

    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "sensitivity_cardiac_query.csv", index=False)
    print("\nSaved: sensitivity_cardiac_query.csv")
    print(df.to_string(index=False))
    return df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--api", action="store_true",
        help="Enable PubMed API calls for cardiac-only relaxation (Table S11). "
             "Without this flag, only the threshold sweep (Table S12) is generated."
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Sensitivity Analysis: Threshold Sweep + Cardiac Relaxation")
    print("=" * 60)

    print("\n--- B: Maturity-threshold sweep (Table S12) ---")
    threshold_df = build_threshold_sensitivity()

    if args.api:
        print("\n--- D: Cardiac-only relaxation (Table S11) ---")
        print("(Querying PubMed E-utilities — cardiac tiers only, neuro/tumor held strict)")
        cardiac_df = build_cardiac_relaxation()
    else:
        print("\n--- D: Cardiac-only relaxation (Table S11) ---")
        print("Skipped (requires --api flag for live PubMed queries).")
        print("Existing sensitivity_cardiac_query.csv is retained from original run.")

    print("\n--- Verification ---")
    for t in [50, 75, 100, 150]:
        col = f"priority_threshold_{t}"
        top = threshold_df.sort_values(col, ascending=False).iloc[0]
        status = "PASS" if top["pathway"] == "Integrin/Src" else "CHECK"
        print(f"  Threshold {t}: #1 = {top['pathway']} ({top[col]}) [{status}]")

    top_cont = threshold_df.sort_values("priority_continuous", ascending=False).iloc[0]
    status = "PASS" if top_cont["pathway"] == "Integrin/Src" else "CHECK"
    print(f"  Continuous:   #1 = {top_cont['pathway']} ({top_cont['priority_continuous']}) [{status}]")

    print("\nDone.")


if __name__ == "__main__":
    main()
