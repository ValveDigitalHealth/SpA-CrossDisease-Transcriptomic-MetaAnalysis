#!/usr/bin/env python3
"""
Phase 5: Functional Enrichment Analysis and Publication Figures.
Uses g:Profiler for GO/Reactome enrichment and DGIdb for drug interactions.
"""

import os
import json
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
META_DIR = os.path.join(RESULTS_DIR, "meta_analysis")
ENRICHMENT_DIR = os.path.join(RESULTS_DIR, "enrichment")
DRUG_DIR = os.path.join(RESULTS_DIR, "drug_repurposing")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(ENRICHMENT_DIR, exist_ok=True)
os.makedirs(DRUG_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)


def run_gprofiler_enrichment(gene_list, organism="hsapiens"):
    """Run g:Profiler functional enrichment analysis."""
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        "organism": organism,
        "query": gene_list,
        "sources": ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "background": None,
        "combined": False,
        "domain_scope": "annotated",
        "no_evidences": False,
    }
    try:
        resp = requests.post(url, json=payload, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        results = data.get("result", [])
        return results
    except Exception as e:
        print(f"  g:Profiler error: {e}")
        return []


def query_dgidb(gene_list):
    """Query DGIdb for drug-gene interactions."""
    url = "https://dgidb.org/api/v2/interactions.json"
    all_results = []
    for i in range(0, len(gene_list), 20):
        batch = gene_list[i:i+20]
        try:
            resp = requests.get(url, params={"genes": ",".join(batch)}, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                matched = data.get("matchedTerms", [])
                for term in matched:
                    gene = term.get("geneName", "")
                    for interaction in term.get("interactions", []):
                        all_results.append({
                            "gene": gene,
                            "drug": interaction.get("drugName", ""),
                            "interaction_type": interaction.get("interactionTypes", [""])[0] if interaction.get("interactionTypes") else "",
                            "fda_approved": interaction.get("drugApproved", False),
                            "source": interaction.get("sources", [""])[0] if interaction.get("sources") else "",
                        })
        except Exception as e:
            print(f"  DGIdb error: {e}")
    return all_results


def create_enrichment_figures(enrich_results, output_prefix):
    """Create enrichment bar plots and dot plots."""
    if not enrich_results:
        return

    df = pd.DataFrame(enrich_results)
    if df.empty:
        return

    # Select top terms for each source
    sources = df["source"].unique() if "source" in df.columns else []
    top_terms = []

    for source in sources:
        source_df = df[df["source"] == source]
        if len(source_df) > 0:
            top = source_df.nsmallest(min(10, len(source_df)), "p_value")
            top_terms.append(top)

    if not top_terms:
        return

    plot_df = pd.concat(top_terms).reset_index(drop=True)
    plot_df["-log10p"] = -np.log10(plot_df["p_value"].clip(1e-300))

    fig, ax = plt.subplots(figsize=(12, 8))
    colors = plt.cm.Set2(np.linspace(0, 1, len(plot_df["source"].unique())))
    color_map = {s: c for s, c in zip(plot_df["source"].unique(), colors)}
    bar_colors = [color_map[s] for s in plot_df["source"]]

    ax.barh(range(len(plot_df)), plot_df["-log10p"], color=bar_colors)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["name"].str[:60] if "name" in plot_df.columns else plot_df.index, fontsize=8)
    ax.set_xlabel("-log10(p-value)")
    ax.set_title("Functional Enrichment Analysis")
    ax.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.5, label='FDR=0.05')
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_enrichment_bar.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_prefix}_enrichment_bar.png")


def main():
    print("Phase 5: Functional Enrichment & Figures")
    print("=" * 60)

    # Load meta-significant genes
    sig_path = os.path.join(META_DIR, "meta_significant_genes.csv")
    if not os.path.exists(sig_path):
        print("ERROR: No meta-significant genes file found. Run Phase 4 first.")
        return

    sig_df = pd.read_csv(sig_path)
    gene_list = sig_df["gene_symbol"].tolist()
    print(f"Meta-significant genes: {len(gene_list)}")
    print(f"  {gene_list[:10]}...")

    # Functional enrichment
    print("\nRunning g:Profiler enrichment...")
    enrich_results = run_gprofiler_enrichment(gene_list)
    print(f"  Enrichment terms: {len(enrich_results)}")

    if enrich_results:
        enrich_df = pd.DataFrame(enrich_results)
        enrich_path = os.path.join(ENRICHMENT_DIR, "enrichment_results.csv")
        enrich_df.to_csv(enrich_path, index=False)
        print(f"  Saved: {enrich_path}")

        # Create figures
        create_enrichment_figures(enrich_results, os.path.join(FIGURES_DIR, "enrichment"))

    # Drug repurposing
    print("\nQuerying DGIdb for drug interactions...")
    drug_results = query_dgidb(gene_list)
    print(f"  Drug interactions: {len(drug_results)}")

    if drug_results:
        drug_df = pd.DataFrame(drug_results)
        drug_path = os.path.join(DRUG_DIR, "drug_interactions.csv")
        drug_df.to_csv(drug_path, index=False)
        print(f"  Saved: {drug_path}")

        fda_drugs = drug_df[drug_df["fda_approved"] == True] if len(drug_df) > 0 else pd.DataFrame()
        print(f"  FDA-approved drugs: {len(fda_drugs)}")

    print("\nPhase 5 complete.")


if __name__ == "__main__":
    main()
