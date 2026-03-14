#!/usr/bin/env python3
"""
Phase 4: Cross-Disease Transcriptomic Meta-Analysis
Combines DEG p-values across AS datasets using Fisher's combined probability test.
Applies vote counting for directional concordance.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from datetime import datetime

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
DEG_DIR = os.path.join(RESULTS_DIR, "deg")
META_DIR = os.path.join(RESULTS_DIR, "meta_analysis")
os.makedirs(META_DIR, exist_ok=True)

FDR_THRESHOLD = 0.05
MIN_DATASETS = 2
LOG2FC_THRESHOLD = 1.0

AS_DATASETS = ["GSE25101", "GSE73754", "GSE18781", "GSE221786"]
CROSS_DISEASE_DATASETS = ["GSE61281", "GSE59071", "GSE58667"]


def load_deg_results(accession):
    """Load gene-level DEG results."""
    gene_path = os.path.join(DEG_DIR, f"{accession}_deg_genes.csv")
    if os.path.exists(gene_path):
        df = pd.read_csv(gene_path)
        if "gene_symbol" in df.columns:
            df = df.set_index("gene_symbol")
            return df
    return None


def fisher_combined_pvalue(pvalues):
    """Fisher's combined probability test."""
    pvalues = [p for p in pvalues if not np.isnan(p) and 0 < p <= 1]
    if len(pvalues) < 2:
        return np.nan
    chi2_stat = -2 * np.sum(np.log(pvalues))
    df = 2 * len(pvalues)
    combined_p = 1 - stats.chi2.cdf(chi2_stat, df)
    return combined_p


def run_meta_analysis(as_datasets_dfs):
    """Run Fisher's meta-analysis across AS datasets."""
    all_genes = set()
    for df in as_datasets_dfs.values():
        if df is not None:
            all_genes.update(df.index)

    print(f"  Total genes across all AS datasets: {len(all_genes)}")
    results = []

    for gene in all_genes:
        pvalues, log2fcs, directions = [], [], []
        datasets_present = []

        for acc, df in as_datasets_dfs.items():
            if df is not None and gene in df.index:
                row = df.loc[gene]
                pval = row.get("pvalue", np.nan) if hasattr(row, 'get') else row["pvalue"]
                lfc = row.get("log2FC", np.nan) if hasattr(row, 'get') else row["log2FC"]

                if not np.isnan(pval) and 0 < pval <= 1:
                    pvalues.append(pval)
                    log2fcs.append(lfc)
                    directions.append("UP" if lfc > 0 else "DOWN")
                    datasets_present.append(acc)

        if len(pvalues) < 2:
            continue

        combined_p = fisher_combined_pvalue(pvalues)
        n_up = sum(1 for d in directions if d == "UP")
        n_down = sum(1 for d in directions if d == "DOWN")
        dominant = "UP" if n_up >= n_down else "DOWN"
        concordance = max(n_up, n_down) / len(directions)

        results.append({
            "gene_symbol": gene,
            "combined_pvalue": combined_p,
            "n_datasets": len(pvalues),
            "n_up": n_up,
            "n_down": n_down,
            "dominant_direction": dominant,
            "concordance": concordance,
            "mean_log2FC": np.mean(log2fcs),
            "datasets": "|".join(datasets_present),
        })

    results_df = pd.DataFrame(results)
    valid = results_df["combined_pvalue"].notna()
    if valid.sum() > 0:
        _, padj, _, _ = multipletests(results_df.loc[valid, "combined_pvalue"].values, method="fdr_bh")
        results_df.loc[valid, "combined_padj"] = padj
    else:
        results_df["combined_padj"] = np.nan

    return results_df


def identify_meta_significant(results_df):
    """Filter for meta-significant genes with directional concordance."""
    sig_mask = (
        (results_df["combined_padj"] < FDR_THRESHOLD) &
        (results_df["n_datasets"] >= MIN_DATASETS) &
        (results_df["concordance"] >= 0.5) &
        (results_df["mean_log2FC"].abs() >= LOG2FC_THRESHOLD)
    )
    return results_df[sig_mask].copy()


def main():
    print("Phase 4: Meta-Analysis")
    print("=" * 60)
    print(f"Timestamp: {datetime.now().strftime('%Y%m%d_%H%M%S')}")

    # Load AS datasets
    as_dfs = {}
    for acc in AS_DATASETS:
        df = load_deg_results(acc)
        as_dfs[acc] = df
        if df is not None:
            print(f"  Loaded {acc}: {len(df)} genes")
        else:
            print(f"  Missing: {acc}")

    loaded = [k for k, v in as_dfs.items() if v is not None]
    if len(loaded) < 2:
        print("ERROR: Need at least 2 AS datasets")
        return

    print(f"\nRunning Fisher's meta-analysis across {len(loaded)} datasets...")
    results_df = run_meta_analysis(as_dfs)

    results_path = os.path.join(META_DIR, "fisher_meta_results.csv")
    results_df.to_csv(results_path, index=False)
    print(f"Full results: {len(results_df)} genes -> {results_path}")

    sig_df = identify_meta_significant(results_df)
    sig_path = os.path.join(META_DIR, "meta_significant_genes.csv")
    sig_df.to_csv(sig_path, index=False)

    print(f"\nMeta-significant genes: {len(sig_df)}")
    print(f"  Upregulated: {(sig_df['dominant_direction'] == 'UP').sum()}")
    print(f"  Downregulated: {(sig_df['dominant_direction'] == 'DOWN').sum()}")

    up_df = sig_df[sig_df["dominant_direction"] == "UP"]
    down_df = sig_df[sig_df["dominant_direction"] == "DOWN"]
    up_df.to_csv(os.path.join(META_DIR, "meta_upregulated.csv"), index=False)
    down_df.to_csv(os.path.join(META_DIR, "meta_downregulated.csv"), index=False)

    if len(sig_df) > 0:
        print("\nTop meta-significant genes:")
        top = sig_df.nsmallest(20, "combined_padj")
        for _, row in top.iterrows():
            print(f"  {row['gene_symbol']:20s} padj={row['combined_padj']:.2e}  "
                  f"log2FC={row['mean_log2FC']:+.2f}  n={row['n_datasets']}  {row['dominant_direction']}")

    # Cross-disease summary
    print("\nCross-disease DEG summary:")
    for acc in CROSS_DISEASE_DATASETS:
        df = load_deg_results(acc)
        if df is not None and "significant" in df.columns:
            n_sig = df["significant"].sum()
            print(f"  {acc}: {n_sig} significant DEGs")

    print("\nPhase 4 complete.")


if __name__ == "__main__":
    main()
