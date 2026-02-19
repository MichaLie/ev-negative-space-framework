#!/usr/bin/env python3
"""Build RNA/protein cargo confidence layer for EV negative-space atlas.

Sources:
- ExoCarta (download archives; EV-specific cargo presence)
- miRBase (RNA identifier validity)
- UniProt REST (protein identifier validity)
- PubMed E-utilities (context evidence, biomarker evidence)
- Final full-text adjudication set (mechanistic anchoring)

Outputs are written to synthesis_output/ for manuscript integration.
"""

from __future__ import annotations

import json
import math
import re
import time
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns

ROOT = Path(__file__).resolve().parent.parent  # repo root: Contro_CC/
SUPP = ROOT / 'supplement'
SUPP.mkdir(parents=True, exist_ok=True)

RUN_DATE = '2026-02-15'
MINDATE = '2015/01/01'
MAXDATE = '2026/02/15'
EV_CLAUSE = '(extracellular vesicle OR exosome OR exosomes)'
BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

EXOCARTA_BASE = 'http://www.exocarta.org/Archive'
EXOCARTA_MIRNA = f'{EXOCARTA_BASE}/ExoCarta_miRNA_details_6.txt'
EXOCARTA_PM = f'{EXOCARTA_BASE}/ExoCarta_protein_mRNA_details_6.txt'
MIRBASE_MATURE = 'https://www.mirbase.org/download/mature.fa'
UNIPROT_SEARCH = 'https://rest.uniprot.org/uniprotkb/search'

CONTEXTS = {
    'neurodegeneration': '(neurodegeneration OR neurodegenerative OR Alzheimer OR Parkinson OR "amyotrophic lateral sclerosis" OR ALS OR frontotemporal dementia OR Huntington)',
    'tumor_metastasis': '((cancer OR tumor OR tumour) AND (metastasis OR metastatic OR premetastatic OR organotropism OR dissemination))',
    'cardiac_repair': '((myocardial infarction OR ischemia-reperfusion OR heart failure OR cardiac injury) AND (repair OR regeneration OR cardioprotection OR post-infarction))',
}

# Candidate panel intentionally constrained to pathway-priority and full-text supported terms.
CANDIDATES = [
    # RNA
    {'cargo_type': 'RNA', 'name': 'miR-21', 'query': '("miR-21" OR "microRNA-21")', 'pathway': 'TGF-beta'},
    {'cargo_type': 'RNA', 'name': 'miR-29', 'query': '("miR-29" OR "microRNA-29")', 'pathway': 'TGF-beta'},
    {'cargo_type': 'RNA', 'name': 'miR-30d', 'query': '("miR-30d" OR "microRNA-30d")', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'RNA', 'name': 'miR-105', 'query': '("miR-105" OR "microRNA-105")', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'RNA', 'name': 'miR-122', 'query': '("miR-122" OR "microRNA-122")', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'RNA', 'name': 'miR-124-3p', 'query': '("miR-124" OR "microRNA-124")', 'pathway': 'mTOR'},
    {'cargo_type': 'RNA', 'name': 'miR-126', 'query': '("miR-126" OR "microRNA-126")', 'pathway': 'Wnt'},
    {'cargo_type': 'RNA', 'name': 'miR-133a-3p', 'query': '("miR-133a" OR "microRNA-133a")', 'pathway': 'TGF-beta'},
    {'cargo_type': 'RNA', 'name': 'miR-143-3p', 'query': '("miR-143" OR "microRNA-143")', 'pathway': 'Autophagy'},
    {'cargo_type': 'RNA', 'name': 'miR-181a-5p', 'query': '("miR-181a" OR "microRNA-181a")', 'pathway': 'Autophagy'},
    {'cargo_type': 'RNA', 'name': 'miR-210-3p', 'query': '("miR-210" OR "microRNA-210")', 'pathway': 'Hypoxia/HIF-1'},
    {'cargo_type': 'RNA', 'name': 'miR-934', 'query': '("miR-934" OR "microRNA-934")', 'pathway': 'mTOR'},

    # Protein/markers
    {'cargo_type': 'Protein', 'name': 'ITGA6', 'query': '(ITGA6 OR "integrin alpha 6")', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'Protein', 'name': 'ITGB4', 'query': '(ITGB4 OR "integrin beta 4")', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'Protein', 'name': 'ITGB2', 'query': '(ITGB2 OR "integrin beta 2")', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'Protein', 'name': 'GPR143', 'query': '(GPR143)', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'Protein', 'name': 'CAV1', 'query': '(CAV1 OR caveolin-1)', 'pathway': 'Integrin/Src'},
    {'cargo_type': 'Protein', 'name': 'TGFB1', 'query': '(TGFB1 OR "TGF-beta1")', 'pathway': 'TGF-beta'},
    {'cargo_type': 'Protein', 'name': 'TGFBR2', 'query': '(TGFBR2 OR "TGF-beta receptor")', 'pathway': 'TGF-beta'},
    {'cargo_type': 'Protein', 'name': 'WNT5A', 'query': '(WNT5A OR Wnt5a)', 'pathway': 'Wnt'},
    {'cargo_type': 'Protein', 'name': 'CTNNB1', 'query': '(CTNNB1 OR beta-catenin)', 'pathway': 'Wnt'},
    {'cargo_type': 'Protein', 'name': 'MTOR', 'query': '(MTOR OR mTOR)', 'pathway': 'mTOR'},
    {'cargo_type': 'Protein', 'name': 'AKT1', 'query': '(AKT1 OR AKT)', 'pathway': 'mTOR'},
    {'cargo_type': 'Protein', 'name': 'NFKB1', 'query': '(NFKB1 OR "NF-kappa B")', 'pathway': 'NF-kB'},
    {'cargo_type': 'Protein', 'name': 'HIF1A', 'query': '(HIF1A OR "HIF-1")', 'pathway': 'Hypoxia/HIF-1'},
    {'cargo_type': 'Protein', 'name': 'C3', 'query': '(C3 OR complement C3)', 'pathway': 'Complement'},
    {'cargo_type': 'Protein', 'name': 'C5', 'query': '(C5 OR complement C5)', 'pathway': 'Complement'},
    {'cargo_type': 'Protein', 'name': 'PDCD6IP', 'query': '(PDCD6IP OR ALIX)', 'pathway': 'Autophagy'},
    {'cargo_type': 'Protein', 'name': 'TSG101', 'query': '(TSG101)', 'pathway': 'Autophagy'},
    {'cargo_type': 'Protein', 'name': 'CD63', 'query': '(CD63)', 'pathway': 'Autophagy'},
    {'cargo_type': 'Protein', 'name': 'CD81', 'query': '(CD81)', 'pathway': 'Autophagy'},
]

BIOMARKER_AXIS = '(biomarker OR diagnostic OR prognostic OR plasma OR serum OR blood)'
MECH_AXIS = '(knockdown OR overexpression OR inhibitor OR mediates OR transfer OR delivery OR caus*)'


def req(url: str, params: Dict | None = None, timeout: int = 60, is_json: bool = False):
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json() if is_json else r.text


def pubmed_count(term: str) -> int:
    params = {
        'db': 'pubmed',
        'term': term,
        'retmode': 'json',
        'mindate': MINDATE,
        'maxdate': MAXDATE,
        'datetype': 'pdat',
        'rettype': 'count',
    }
    d = req(f'{BASE}/esearch.fcgi', params=params, is_json=True)
    return int(d['esearchresult'].get('count', '0'))


def download_sources() -> Dict[str, Path]:
    cache = ROOT / 'source_cache'
    cache.mkdir(exist_ok=True)

    files = {
        'exocarta_mirna': cache / 'ExoCarta_miRNA_details_6.txt',
        'exocarta_pm': cache / 'ExoCarta_protein_mRNA_details_6.txt',
        'mirbase_mature': cache / 'mirbase_mature.fa',
    }

    if not files['exocarta_mirna'].exists():
        files['exocarta_mirna'].write_text(req(EXOCARTA_MIRNA), encoding='utf-8', errors='ignore')
    if not files['exocarta_pm'].exists():
        files['exocarta_pm'].write_text(req(EXOCARTA_PM), encoding='utf-8', errors='ignore')
    if not files['mirbase_mature'].exists():
        files['mirbase_mature'].write_text(req(MIRBASE_MATURE), encoding='utf-8', errors='ignore')

    return files


def load_exocarta(files: Dict[str, Path]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    mirna = pd.read_csv(files['exocarta_mirna'], sep='\t', dtype=str)
    pm = pd.read_csv(files['exocarta_pm'], sep='\t', dtype=str)
    mirna.columns = [c.strip() for c in mirna.columns]
    pm.columns = [c.strip() for c in pm.columns]

    mirna['MIRNA ID'] = mirna['MIRNA ID'].fillna('').str.strip()
    mirna['SPECIES'] = mirna['SPECIES'].fillna('').str.strip()
    mirna['EXPERIMENT ID'] = mirna['EXPERIMENT ID'].fillna('').str.strip()

    pm['GENE SYMBOL'] = pm['GENE SYMBOL'].fillna('').str.strip().str.upper()
    pm['SPECIES'] = pm['SPECIES'].fillna('').str.strip()
    pm['EXPERIMENT_ID'] = pm['EXPERIMENT_ID'].fillna('').str.strip()
    pm['CONTENT_TYPE'] = pm['CONTENT_TYPE'].fillna('').str.strip().str.lower()

    return mirna, pm


def load_mirbase_names(files: Dict[str, Path]) -> set:
    names = set()
    for line in files['mirbase_mature'].read_text(encoding='utf-8', errors='ignore').splitlines():
        if line.startswith('>'):
            # observed format: >hsa-miR-21-5p MIMAT0000076 Homo sapiens miR-21-5p
            # historical alternatives can swap token order.
            parts = line[1:].split()
            if not parts:
                continue
            for tok in parts[:2]:
                t = tok.strip().lower()
                if t.startswith('hsa-mir-') or t.startswith('mir-'):
                    names.add(t)
    return names


def normalize_mir_name(name: str) -> List[str]:
    raw = name.lower().strip().replace('_', '-').replace(' ', '')
    raw = raw.replace('microrna-', 'mir-').replace('microrna', 'mir')
    raw = raw.replace('mirna-', 'mir-').replace('mirna', 'mir')
    if raw.startswith('hsa-'):
        raw = raw[4:]
    if raw.startswith('mir-'):
        base = raw[4:]
    elif raw.startswith('mir'):
        base = raw[3:]
    else:
        base = raw

    variants = {
        f'mir-{base}',
        f'hsa-mir-{base}',
    }

    # Add arm-specific variants when arm is unspecified.
    if not re.search(r'-(3p|5p)$', base):
        variants.update({
            f'mir-{base}-3p',
            f'mir-{base}-5p',
            f'hsa-mir-{base}-3p',
            f'hsa-mir-{base}-5p',
        })

    # Add arm-unspecified variant when candidate includes an arm.
    base2 = re.sub(r'-(3p|5p)$', '', base)
    variants.update({
        f'mir-{base2}',
        f'hsa-mir-{base2}',
    })

    return sorted(v.lower() for v in variants)


def mirbase_validate(name: str, mirbase_names: set) -> Tuple[bool, str]:
    candidates = set(normalize_mir_name(name))
    matches = {m for m in candidates if m in mirbase_names}

    # Family-level fallback for names like miR-29 that map to multiple members (e.g., 29a/29b/29c).
    base = name.lower().strip().replace('_', '-').replace(' ', '')
    base = base.replace('microrna-', 'mir-').replace('microrna', 'mir')
    base = base.replace('mirna-', 'mir-').replace('mirna', 'mir')
    if base.startswith('hsa-'):
        base = base[4:]
    if base.startswith('mir-'):
        base = base[4:]
    elif base.startswith('mir'):
        base = base[3:]
    family_prefix = f'hsa-mir-{base}'
    if re.fullmatch(r'\d+', base):
        family_pat = re.compile(rf'^{re.escape(family_prefix)}([a-z]|-|$)')
    else:
        family_pat = re.compile(rf'^{re.escape(family_prefix)}($|-)')
    for mb in mirbase_names:
        if family_pat.match(mb):
            matches.add(mb)

    ordered = sorted(matches)
    return (len(ordered) > 0, ';'.join(ordered[:6])[:180])


def uniprot_validate(symbol: str) -> Tuple[bool, str, str]:
    q = f'gene:{symbol} AND organism_id:9606 AND reviewed:true'
    params = {'query': q, 'format': 'json', 'size': 1}
    try:
        d = req(UNIPROT_SEARCH, params=params, is_json=True)
        res = d.get('results', [])
        if not res:
            return False, '', ''
        r = res[0]
        acc = r.get('primaryAccession', '')
        upid = r.get('uniProtkbId', '')
        return True, acc, upid
    except Exception:
        return False, '', ''


def candidate_in_fulltext(name: str, cargo_type: str, ft_df: pd.DataFrame) -> bool:
    n = name.lower()
    if cargo_type == 'RNA':
        pats = [n, n.replace('mir-', 'miR-').lower(), n.replace('-3p', '').replace('-5p', '')]
    else:
        pats = [n.lower()]
    title_blob = ' '.join(ft_df['title'].fillna('').astype(str).str.lower().tolist())
    return any(p in title_blob for p in pats if p)


def run():
    files = download_sources()
    mirna_exo, pm_exo = load_exocarta(files)
    mirbase_names = load_mirbase_names(files)

    ft = pd.read_csv(ROOT / 'fulltext_mechanistic_adjudication_final.csv')

    rows = []
    query_log = []

    # precompute exocarta normalization denominators
    # will compute on panel values later

    for c in CANDIDATES:
        cargo_type = c['cargo_type']
        name = c['name']
        entry = {
            'run_date': RUN_DATE,
            'cargo_type': cargo_type,
            'candidate': name,
            'pathway': c['pathway'],
            'candidate_query': c['query'],
        }

        # Identifier validity
        if cargo_type == 'RNA':
            ok, detail = mirbase_validate(name, mirbase_names)
            entry['identifier_source'] = 'miRBase mature.fa'
            entry['identifier_valid'] = ok
            entry['identifier_detail'] = detail
            entry['uniprot_accession'] = ''
            entry['uniprot_id'] = ''
        else:
            ok, acc, upid = uniprot_validate(name)
            entry['identifier_source'] = 'UniProt reviewed Homo sapiens'
            entry['identifier_valid'] = ok
            entry['identifier_detail'] = upid
            entry['uniprot_accession'] = acc
            entry['uniprot_id'] = upid

        # ExoCarta evidence
        if cargo_type == 'RNA':
            mir_vars = set(normalize_mir_name(name))
            # Also add family-level matching for ExoCarta (e.g., miR-29 -> miR-29a, miR-29b, miR-29c)
            raw = name.lower().strip().replace('_', '-').replace(' ', '')
            raw = raw.replace('microrna-', 'mir-').replace('microrna', 'mir')
            raw = raw.replace('mirna-', 'mir-').replace('mirna', 'mir')
            if raw.startswith('hsa-'):
                raw = raw[4:]
            if raw.startswith('mir-'):
                base_num = raw[4:]
            elif raw.startswith('mir'):
                base_num = raw[3:]
            else:
                base_num = raw
            # For numeric-only bases (family-level names like miR-29), also match family members
            if re.fullmatch(r'\d+', base_num):
                family_pat = re.compile(rf'^hsa-mir-{re.escape(base_num)}([a-z]|-|$)', re.IGNORECASE)
                exo_ids_lower = mirna_exo['MIRNA ID'].str.lower()
                family_mask = exo_ids_lower.apply(lambda x: bool(family_pat.match(x)))
                mir_vars_mask = exo_ids_lower.isin(mir_vars)
                m = mirna_exo[mir_vars_mask | family_mask].copy()
            else:
                m = mirna_exo[mirna_exo['MIRNA ID'].str.lower().isin(mir_vars)].copy()
            mh = m[m['SPECIES'].str.lower() == 'homo sapiens']
            entry['exocarta_present_any'] = len(m) > 0
            entry['exocarta_present_human'] = len(mh) > 0
            entry['exocarta_experiments_any'] = m['EXPERIMENT ID'].nunique()
            entry['exocarta_experiments_human'] = mh['EXPERIMENT ID'].nunique()
        else:
            p = pm_exo[pm_exo['GENE SYMBOL'].str.upper() == name.upper()].copy()
            ph = p[p['SPECIES'].str.lower() == 'homo sapiens']
            entry['exocarta_present_any'] = len(p) > 0
            entry['exocarta_present_human'] = len(ph) > 0
            entry['exocarta_experiments_any'] = p['EXPERIMENT_ID'].nunique()
            entry['exocarta_experiments_human'] = ph['EXPERIMENT_ID'].nunique()

        # PubMed context evidence + biomarker/mechanistic evidence
        total_hits = 0
        total_biomarker = 0
        total_mech = 0
        contexts_with_hits = 0

        for ctx, ctx_q in CONTEXTS.items():
            base_q = f'{EV_CLAUSE} AND {ctx_q} AND {c["query"]}'
            cnt = pubmed_count(base_q)
            bio = pubmed_count(f'{base_q} AND {BIOMARKER_AXIS}')
            mech = pubmed_count(f'{base_q} AND {MECH_AXIS}')
            query_log.extend([
                {'candidate': name, 'cargo_type': cargo_type, 'context': ctx, 'query_type': 'base', 'query': base_q, 'count': cnt},
                {'candidate': name, 'cargo_type': cargo_type, 'context': ctx, 'query_type': 'biomarker', 'query': f'{base_q} AND {BIOMARKER_AXIS}', 'count': bio},
                {'candidate': name, 'cargo_type': cargo_type, 'context': ctx, 'query_type': 'mechanistic', 'query': f'{base_q} AND {MECH_AXIS}', 'count': mech},
            ])
            entry[f'pubmed_{ctx}'] = cnt
            entry[f'biomarker_{ctx}'] = bio
            entry[f'mechanistic_{ctx}'] = mech
            total_hits += cnt
            total_biomarker += bio
            total_mech += mech
            if cnt > 0:
                contexts_with_hits += 1
            time.sleep(0.2)

        entry['pubmed_total_hits'] = total_hits
        entry['biomarker_total_hits'] = total_biomarker
        entry['mechanistic_total_hits'] = total_mech
        entry['contexts_with_hits'] = contexts_with_hits
        entry['context_breadth_fraction'] = contexts_with_hits / 3.0

        # Full-text support
        entry['fulltext_title_support'] = candidate_in_fulltext(name, cargo_type, ft)

        rows.append(entry)

    df = pd.DataFrame(rows)

    # scoring components
    # identifier score
    df['identifier_score'] = df['identifier_valid'].astype(float)

    # ExoCarta score
    human_log = df['exocarta_experiments_human'].fillna(0).astype(float).map(lambda x: math.log1p(x))
    max_human_log = max(human_log.max(), 1.0)
    df['exocarta_count_norm'] = human_log / max_human_log
    df['exocarta_presence_score'] = df['exocarta_present_human'].astype(float)
    df['exocarta_score'] = 0.5 * df['exocarta_presence_score'] + 0.5 * df['exocarta_count_norm']

    # breadth score
    total_log = df['pubmed_total_hits'].fillna(0).astype(float).map(lambda x: math.log1p(x))
    max_total_log = max(total_log.max(), 1.0)
    df['pubmed_total_norm'] = total_log / max_total_log
    df['breadth_score'] = 0.5 * df['context_breadth_fraction'] + 0.5 * df['pubmed_total_norm']

    # mechanistic score
    mech_log = df['mechanistic_total_hits'].fillna(0).astype(float).map(lambda x: math.log1p(x))
    max_mech_log = max(mech_log.max(), 1.0)
    df['mechanistic_norm'] = mech_log / max_mech_log
    df['fulltext_support_score'] = df['fulltext_title_support'].astype(float)
    df['mechanistic_score'] = 0.6 * df['mechanistic_norm'] + 0.4 * df['fulltext_support_score']

    # biomarker readiness
    bio_log = df['biomarker_total_hits'].fillna(0).astype(float).map(lambda x: math.log1p(x))
    max_bio_log = max(bio_log.max(), 1.0)
    df['biomarker_norm'] = bio_log / max_bio_log

    # Final confidence score (0-100)
    df['cargo_confidence_score'] = 100 * (
        0.25 * df['identifier_score'] +
        0.25 * df['exocarta_score'] +
        0.25 * df['breadth_score'] +
        0.15 * df['mechanistic_score'] +
        0.10 * df['biomarker_norm']
    )

    def tier(x: float) -> str:
        if x >= 70:
            return 'High'
        if x >= 40:
            return 'Moderate'
        return 'Low'

    df['cargo_confidence_tier'] = df['cargo_confidence_score'].map(tier)

    # Save tables
    out_main = ROOT / 'cargo_confidence_rna_protein.csv'
    out_queries = SUPP / 'cargo_confidence_query_log.csv'
    out_pathway = ROOT / 'pathway_cargo_confidence_summary.csv'

    df.sort_values(['cargo_type', 'cargo_confidence_score'], ascending=[True, False]).to_csv(out_main, index=False)
    pd.DataFrame(query_log).to_csv(out_queries, index=False)

    path_summary = df.groupby('pathway', as_index=False).agg(
        n_candidates=('candidate', 'count'),
        mean_confidence=('cargo_confidence_score', 'mean'),
        median_confidence=('cargo_confidence_score', 'median'),
        n_high=('cargo_confidence_tier', lambda s: int((s == 'High').sum())),
        n_moderate=('cargo_confidence_tier', lambda s: int((s == 'Moderate').sum())),
        n_low=('cargo_confidence_tier', lambda s: int((s == 'Low').sum())),
    )
    path_summary.to_csv(out_pathway, index=False)

    # Integrate with negative-space priority
    nsp = pd.read_csv(ROOT / 'negative_space_priority_index.csv')
    merged = nsp.merge(path_summary, on='pathway', how='left')
    merged['mean_confidence'] = merged['mean_confidence'].fillna(0)
    merged['actionability_index'] = merged['negative_space_priority_index'] * (merged['mean_confidence'] / 100.0)
    merged = merged.sort_values('actionability_index', ascending=False)
    merged.to_csv(ROOT / 'negative_space_with_cargo_confidence.csv', index=False)

    # figures
    sns.set_theme(style='whitegrid')

    plt.figure(figsize=(10, 5.8))
    tmp = df.sort_values('cargo_confidence_score', ascending=False)
    sns.barplot(data=tmp, y='candidate', x='cargo_confidence_score', hue='cargo_type', dodge=False)
    # Explicit threshold marker for High vs Moderate tiers
    plt.axvline(70, color='#111827', linestyle='--', linewidth=1.6, alpha=0.9)
    plt.title('Figure 6. RNA/Protein Cargo Confidence Scores (EV-focused multi-source evidence)')
    plt.xlabel('Cargo confidence score (0-100)')
    plt.ylabel('Candidate cargo')
    plt.legend(title='Cargo type', loc='lower right')
    plt.tight_layout()
    plt.savefig(ROOT / 'figure6_cargo_confidence_scores.png', dpi=320)
    plt.close()

    fig, ax = plt.subplots(figsize=(9.4, 6.4))
    sns.scatterplot(
        data=merged,
        x='negative_space_priority_index',
        y='mean_confidence',
        size='n_candidates',
        hue='pathway',
        legend=False,
        sizes=(40, 260),
        ax=ax,
    )
    x_min = float(merged['negative_space_priority_index'].min())
    x_max = float(merged['negative_space_priority_index'].max())
    x_pad = max(0.8, 0.12 * (x_max - x_min))
    ax.set_xlim(x_min - 0.6, x_max + x_pad)

    # Keep labels readable near plot edges by flipping alignment at the right boundary.
    right_edge_threshold = x_max - 0.25
    for _, r in merged.iterrows():
        x = float(r['negative_space_priority_index'])
        y = float(r['mean_confidence'])
        if x >= right_edge_threshold:
            dx, ha = -0.08, 'right'
        else:
            dx, ha = 0.06, 'left'
        ax.text(x + dx, y + 0.2, r['pathway'], fontsize=8, ha=ha, clip_on=False)
    ax.set_title('Figure 7. Pathway Negative-Space vs Cargo Confidence')
    ax.set_xlabel('Negative-space priority index')
    ax.set_ylabel('Mean cargo confidence')
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.12, top=0.92)
    fig.savefig(ROOT / 'figure7_negative_space_vs_cargo_confidence.png', dpi=320)
    plt.close(fig)

    # methods/source manifest
    manifest = {
        'run_date': RUN_DATE,
        'date_window': {'mindate': MINDATE, 'maxdate': MAXDATE},
        'sources': {
            'exocarta_mirna': EXOCARTA_MIRNA,
            'exocarta_protein_mrna': EXOCARTA_PM,
            'mirbase_mature': MIRBASE_MATURE,
            'uniprot_search': UNIPROT_SEARCH,
            'pubmed_eutils': BASE,
            'fulltext_adjudication': 'REPO_ROOT/fulltext_mechanistic_adjudication_final.csv',
            'optional_not_used_for_scoring': [
                'EV-TRACK (no open export endpoint identified in this run)',
                'RNAcentral (endpoint instability/timeouts during this run)',
            ],
        },
        'scoring_weights': {
            'identifier': 0.25,
            'exocarta': 0.25,
            'breadth': 0.25,
            'mechanistic': 0.15,
            'biomarker': 0.10,
        },
        'n_candidates': len(CANDIDATES),
    }
    (SUPP / 'cargo_confidence_manifest.json').write_text(json.dumps(manifest, indent=2))

    print('Saved:', out_main)
    print('Saved:', out_queries)
    print('Saved:', out_pathway)
    print('Saved:', ROOT / 'negative_space_with_cargo_confidence.csv')
    print('Top 10 candidates:')
    print(df.sort_values('cargo_confidence_score', ascending=False)[['candidate', 'cargo_type', 'pathway', 'cargo_confidence_score', 'cargo_confidence_tier']].head(10).to_string(index=False))


if __name__ == '__main__':
    run()
