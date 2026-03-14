#!/usr/bin/env python3
"""
Phase 3: Per-dataset Differential Expression Analysis
Uses Welch's t-test with BH correction for microarray and RNA-seq.
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

LOG2FC_THRESHOLD = 1.0
PVAL_THRESHOLD = 0.05


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

    disease_samples = [s for s, info in sample_info.items() if info["group"] == "disease"]
    control_samples = [s for s, info in sample_info.items() if info["group"] == "control"]

    disease_cols = [c for c in df.columns if c in disease_samples]
    control_cols = [c for c in df.columns if c in control_samples]

    if not disease_cols and not control_cols:
        disease_cols = [c for c in df.columns if any(c.startswith(p) for p in ['AS_', 'disease'])]
        control_cols = [c for c in df.columns if any(c.startswith(p) for p in ['Ctrl_', 'HC_', 'control'])]

    print(f"  Disease: {len(disease_cols)}, Control: {len(control_cols)}")
    return df, disease_cols, control_cols


def run_differential_expression(df, disease_cols, control_cols, accession, data_type="microarray"):
    """Run Welch's t-test per gene with BH correction."""
    if len(disease_cols) < 3 or len(control_cols) < 3:
        print(f"  ERROR: Insufficient samples")
        return None

    disease_expr = df[disease_cols].values.astype(float)
    control_expr = df[control_cols].values.astype(float)

    if data_type == "rnaseq":
        lib_sizes_d = disease_expr.sum(axis=0)
        lib_sizes_c = control_expr.sum(axis=0)
        disease_expr = np.log2(disease_expr / lib_sizes_d * 1e6 + 1)
        control_expr = np.log2(control_expr / lib_sizes_c * 1e6 + 1)

    n_genes = df.shape[0]
    results = []

    for i in range(n_genes):
        d_vals = disease_expr[i, ~np.isnan(disease_expr[i, :])]
        c_vals = control_expr[i, ~np.isnan(control_expr[i, :])]

        if len(d_vals) < 3 or len(c_vals) < 3:
            results.append({"gene_id": df.index[i], "log2FC": np.nan, "t_stat": np.nan, "pvalue": np.nan})
            continue

        log2fc = np.mean(d_vals) - np.mean(c_vals)

        try:
            t_stat, pvalue = stats.ttest_ind(d_vals, c_vals, equal_var=False)
        except Exception:
            t_stat, pvalue = np.nan, np.nan

        results.append({"gene_id": df.index[i], "log2FC": log2fc,
                        "mean_disease": np.mean(d_vals), "mean_control": np.mean(c_vals),
                        "t_stat": t_stat, "pvalue": pvalue})

    results_df = pd.DataFrame(results)

    valid_mask = ~results_df["pvalue"].isna()
    if valid_mask.sum() > 0:
        _, padj, _, _ = multipletests(results_df.loc[valid_mask, "pvalue"].values, method="fdr_bh")
        results_df.loc[valid_mask, "padj"] = padj
    else:
        results_df["padj"] = np.nan

    results_df["significant"] = (
        (results_df["padj"] < PVAL_THRESHOLD) &
        (results_df["log2FC"].abs() > LOG2FC_THRESHOLD)
    )

    results_df["direction"] = "NS"
    results_df.loc[results_df["significant"] & (results_df["log2FC"] > 0), "direction"] = "UP"
    results_df.loc[results_df["significant"] & (results_df["log2FC"] < 0), "direction"] = "DOWN"

    return results_df


def main():
    print("Phase 3: Differential Expression Analysis")
    print("=" * 60)

    datasets = [
        {"accession": "GSE25101", "disease": "AS", "type": "microarray"},
        {"accession": "GSE73754", "disease": "AS", "type": "microarray"},
        {"accession": "GSE18781", "disease": "AS", "type": "microarray"},
        {"accession": "GSE61281", "disease": "PsA", "type": "microarray"},
        {"accession": "GSE59071", "disease": "IBD", "type": "microarray"},
        {"accession": "GSE58667", "disease": "jSpA", "type": "microarray"},
        {"accession": "GSE221786", "disease": "AS", "type": "rnaseq"},
    ]

    all_summaries = []

    for ds in datasets:
        acc = ds["accession"]
        print(f"\nAnalyzing: {acc} ({ds['disease']}, {ds['type']})")

        df, disease_cols, control_cols = load_expression_and_groups(acc, ds["type"])
        if df is None:
            continue

        results_df = run_differential_expression(df, disease_cols, control_cols, acc, ds["type"])
        if results_df is not None:
            out_path = os.path.join(RESULTS_DIR, f"{acc}_deg_results.csv")
            results_df.to_csv(out_path, index=False)

            n_sig = results_df["significant"].sum()
            n_up = (results_df["direction"] == "UP").sum()
            n_down = (results_df["direction"] == "DOWN").sum()
            print(f"  Significant: {n_sig} (UP: {n_up}, DOWN: {n_down})")
            all_summaries.append({"accession": acc, "disease": ds["disease"],
                                   "n_significant": int(n_sig), "n_up": int(n_up), "n_down": int(n_down)})

    summary_path = os.path.join(RESULTS_DIR, "phase3_summary.json")
    with open(summary_path, "w") as f:
        json.dump(all_summaries, f, indent=2)
    print(f"\nSummary saved: {summary_path}")


if __name__ == "__main__":
    main()
