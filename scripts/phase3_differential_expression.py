#!/usr/bin/env python3
"""
Phase 3: Per-dataset Differential Expression Analysis
Uses limma-equivalent analysis in Python (scipy + statsmodels) for microarray,
and PyDESeq2 or manual Wald test for RNA-seq counts.

For microarray: log2 expression → t-test with BH correction (limma equivalent)
For RNA-seq counts: log2(CPM+1) → t-test with BH correction
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from datetime import datetime

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
PROCESSED_DIR = os.path.join(BASE_DIR, "data/processed")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")
RESULTS_DIR = os.path.join(BASE_DIR, "results/deg")
os.makedirs(RESULTS_DIR, exist_ok=True)

# Thresholds
LOG2FC_THRESHOLD = 1.0   # |log2FC| > 1
PVAL_THRESHOLD = 0.05    # adjusted p-value < 0.05


def load_expression_and_groups(accession, data_type="microarray"):
    """Load expression matrix and sample group assignments."""
    
    if data_type == "microarray":
        expr_path = os.path.join(PROCESSED_DIR, f"{accession}_expression.csv")
    else:
        expr_path = os.path.join(PROCESSED_DIR, f"{accession}_counts.csv")
    
    meta_path = os.path.join(METADATA_DIR, f"{accession}_samples.json")
    
    if not os.path.exists(expr_path) or not os.path.exists(meta_path):
        print(f"  Missing files for {accession}")
        return None, None, None
    
    df = pd.read_csv(expr_path, index_col=0)
    
    with open(meta_path) as f:
        sample_info = json.load(f)
    
    # Handle gene_symbol column for RNA-seq
    gene_symbols = None
    if 'gene_symbol' in df.columns:
        gene_symbols = df['gene_symbol']
        df = df.drop('gene_symbol', axis=1)
    
    # Get disease and control sample IDs
    disease_samples = [s for s, info in sample_info.items() if info["group"] == "disease"]
    control_samples = [s for s, info in sample_info.items() if info["group"] == "control"]
    
    # Match to expression columns
    disease_cols = [c for c in df.columns if c in disease_samples]
    control_cols = [c for c in df.columns if c in control_samples]
    
    # For datasets where sample IDs don't match column names exactly,
    # try matching by position or partial match
    if not disease_cols and not control_cols:
        # Try using column names directly (e.g., GSE221786 has Ctrl_F_0, AS_F_0)
        disease_cols = [c for c in df.columns if any(c.startswith(prefix) for prefix in ['AS_', 'disease'])]
        control_cols = [c for c in df.columns if any(c.startswith(prefix) for prefix in ['Ctrl_', 'HC_', 'control'])]
    
    print(f"  Disease samples matched: {len(disease_cols)}")
    print(f"  Control samples matched: {len(control_cols)}")
    
    return df, disease_cols, control_cols


def run_differential_expression(df, disease_cols, control_cols, accession, data_type="microarray"):
    """Run differential expression analysis."""
    
    if len(disease_cols) < 3 or len(control_cols) < 3:
        print(f"  ERROR: Insufficient samples (disease={len(disease_cols)}, control={len(control_cols)})")
        return None
    
    # Get expression values for each group
    disease_expr = df[disease_cols].values.astype(float)
    control_expr = df[control_cols].values.astype(float)
    
    # For RNA-seq counts, convert to log2(CPM+1)
    if data_type == "rnaseq":
        # CPM normalization
        lib_sizes_d = disease_expr.sum(axis=0)
        lib_sizes_c = control_expr.sum(axis=0)
        
        disease_expr = np.log2(disease_expr / lib_sizes_d * 1e6 + 1)
        control_expr = np.log2(control_expr / lib_sizes_c * 1e6 + 1)
    
    n_genes = df.shape[0]
    results = []
    
    for i in range(n_genes):
        d_vals = disease_expr[i, :]
        c_vals = control_expr[i, :]
        
        # Remove NaN
        d_vals = d_vals[~np.isnan(d_vals)]
        c_vals = c_vals[~np.isnan(c_vals)]
        
        if len(d_vals) < 3 or len(c_vals) < 3:
            results.append({
                "gene_id": df.index[i],
                "log2FC": np.nan,
                "t_stat": np.nan,
                "pvalue": np.nan,
            })
            continue
        
        # Mean expression filter: skip very low expression genes
        mean_expr = (np.mean(d_vals) + np.mean(c_vals)) / 2
        if data_type == "microarray" and mean_expr < 3:  # log2 scale
            results.append({
                "gene_id": df.index[i],
                "log2FC": np.nan,
                "t_stat": np.nan,
                "pvalue": np.nan,
            })
            continue
        if data_type == "rnaseq" and mean_expr < 1:  # log2(CPM+1) scale
            results.append({
                "gene_id": df.index[i],
                "log2FC": np.nan,
                "t_stat": np.nan,
                "pvalue": np.nan,
            })
            continue
        
        # Calculate log2 fold change
        log2fc = np.mean(d_vals) - np.mean(c_vals)
        
        # Welch's t-test (unequal variance)
        try:
            t_stat, pvalue = stats.ttest_ind(d_vals, c_vals, equal_var=False)
        except:
            t_stat, pvalue = np.nan, np.nan
        
        results.append({
            "gene_id": df.index[i],
            "log2FC": log2fc,
            "mean_disease": np.mean(d_vals),
            "mean_control": np.mean(c_vals),
            "t_stat": t_stat,
            "pvalue": pvalue,
        })
    
    results_df = pd.DataFrame(results)
    
    # BH multiple testing correction
    valid_mask = ~results_df["pvalue"].isna()
    if valid_mask.sum() > 0:
        _, padj, _, _ = multipletests(
            results_df.loc[valid_mask, "pvalue"].values,
            method="fdr_bh"
        )
        results_df.loc[valid_mask, "padj"] = padj
    else:
        results_df["padj"] = np.nan
    
    # Classify significance
    results_df["significant"] = (
        (results_df["padj"] < PVAL_THRESHOLD) & 
        (results_df["log2FC"].abs() > LOG2FC_THRESHOLD)
    )
    
    results_df["direction"] = "NS"
    results_df.loc[results_df["significant"] & (results_df["log2FC"] > 0), "direction"] = "UP"
    results_df.loc[results_df["significant"] & (results_df["log2FC"] < 0), "direction"] = "DOWN"
    
    return results_df


def summarize_deg_results(results_df, accession):
    """Print and save summary of DEG results."""
    n_total = len(results_df)
    n_tested = results_df["pvalue"].notna().sum()
    n_sig = results_df["significant"].sum()
    n_up = (results_df["direction"] == "UP").sum()
    n_down = (results_df["direction"] == "DOWN").sum()
    
    print(f"\n  --- {accession} DEG Summary ---")
    print(f"  Total probes/genes: {n_total}")
    print(f"  Tested (passed filters): {n_tested}")
    print(f"  Significant (|log2FC|>{LOG2FC_THRESHOLD}, padj<{PVAL_THRESHOLD}): {n_sig}")
    print(f"    Upregulated: {n_up}")
    print(f"    Downregulated: {n_down}")
    
    # Top DEGs
    if n_sig > 0:
        top_up = results_df[results_df["direction"] == "UP"].nsmallest(10, "padj")
        top_down = results_df[results_df["direction"] == "DOWN"].nsmallest(10, "padj")
        
        print(f"\n  Top 10 Upregulated:")
        for _, row in top_up.iterrows():
            print(f"    {str(row['gene_id']):20s} log2FC={row['log2FC']:+.2f}  padj={row['padj']:.2e}")
        
        print(f"\n  Top 10 Downregulated:")
        for _, row in top_down.iterrows():
            print(f"    {str(row['gene_id']):20s} log2FC={row['log2FC']:+.2f}  padj={row['padj']:.2e}")
    
    return {
        "accession": accession,
        "n_total": int(n_total),
        "n_tested": int(n_tested),
        "n_significant": int(n_sig),
        "n_up": int(n_up),
        "n_down": int(n_down),
    }


def main():
    print("=" * 80)
    print("Phase 3: Per-Dataset Differential Expression Analysis")
    print(f"Timestamp: {datetime.now().strftime('%Y%m%d_%H%M%S')}")
    print(f"Thresholds: |log2FC| > {LOG2FC_THRESHOLD}, padj < {PVAL_THRESHOLD}")
    print("=" * 80)
    
    # Define datasets to analyze
    datasets = [
        # Microarray datasets
        {"accession": "GSE25101", "disease": "AS", "type": "microarray"},
        {"accession": "GSE73754", "disease": "AS", "type": "microarray"},
        {"accession": "GSE18781", "disease": "AS/axSpA", "type": "microarray"},
        {"accession": "GSE61281", "disease": "PsA", "type": "microarray"},
        {"accession": "GSE59071", "disease": "IBD", "type": "microarray"},
        {"accession": "GSE58667", "disease": "jSpA", "type": "microarray"},
        # RNA-seq
        {"accession": "GSE221786", "disease": "AS", "type": "rnaseq"},
    ]
    
    all_summaries = []
    
    for ds in datasets:
        acc = ds["accession"]
        print(f"\n{'='*60}")
        print(f"Analyzing: {acc} ({ds['disease']}, {ds['type']})")
        print(f"{'='*60}")
        
        df, disease_cols, control_cols = load_expression_and_groups(acc, ds["type"])
        
        if df is None:
            print(f"  SKIPPED: Could not load data")
            continue
        
        if len(disease_cols) < 3 or len(control_cols) < 3:
            print(f"  SKIPPED: Insufficient samples")
            continue
        
        results_df = run_differential_expression(df, disease_cols, control_cols, acc, ds["type"])
        
        if results_df is not None:
            # Save full results
            output_path = os.path.join(RESULTS_DIR, f"{acc}_deg_results.csv")
            results_df.to_csv(output_path, index=False)
            
            # Save significant only
            sig_df = results_df[results_df["significant"]]
            sig_path = os.path.join(RESULTS_DIR, f"{acc}_deg_significant.csv")
            sig_df.to_csv(sig_path, index=False)
            
            summary = summarize_deg_results(results_df, acc)
            summary["disease"] = ds["disease"]
            all_summaries.append(summary)
    
    # Overall summary
    print(f"\n{'='*80}")
    print("PHASE 3 OVERALL SUMMARY")
    print(f"{'='*80}")
    
    print(f"\n{'Accession':12s} {'Disease':10s} {'Tested':>8s} {'DEGs':>8s} {'Up':>6s} {'Down':>6s}")
    print("-" * 60)
    for s in all_summaries:
        print(f"{s['accession']:12s} {s['disease']:10s} {s['n_tested']:>8d} {s['n_significant']:>8d} {s['n_up']:>6d} {s['n_down']:>6d}")
    
    total_degs = sum(s["n_significant"] for s in all_summaries)
    print(f"\nTotal DEGs across all datasets: {total_degs}")
    
    # Save summary
    summary_path = os.path.join(RESULTS_DIR, "phase3_summary.json")
    with open(summary_path, "w") as f:
        json.dump(all_summaries, f, indent=2)
    print(f"\nSummary saved to: {summary_path}")


if __name__ == "__main__":
    main()
