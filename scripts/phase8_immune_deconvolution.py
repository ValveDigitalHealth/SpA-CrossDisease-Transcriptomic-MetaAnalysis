#!/usr/bin/env python3
"""
Phase 8: Immune Cell Deconvolution
====================================
Cross-Disease Transcriptomic Meta-Analysis of SpA

Applies single-sample Gene Set Enrichment Analysis (ssGSEA) to estimate
relative abundance of 13 immune cell types in each sample.
Compares cell type scores between disease and control groups.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
DECONV_DIR = os.path.join(RESULTS_DIR, "deconvolution")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
METADATA_DIR = os.path.join(BASE_DIR, "data", "metadata")
os.makedirs(DECONV_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# 13 immune cell types with representative gene signatures
IMMUNE_SIGNATURES = {
    "B cells naive": ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "FCER2", "IGHM", "IGHD", "IL4R", "TCL1A"],
    "B cells memory": ["CD19", "MS4A1", "CD27", "CD38", "AIM2", "TNFRSF17", "JCHAIN", "SDC1"],
    "T cells CD4 naive": ["CD3D", "CD3E", "CD4", "CCR7", "SELL", "IL7R", "TCF7", "LEF1", "KLF2"],
    "T cells CD4 memory": ["CD3D", "CD3E", "CD4", "IL7R", "CD44", "FAS", "ICOS", "PDCD1", "LAG3"],
    "T cells CD8": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GZMA", "GNLY"],
    "T cells regulatory (Tregs)": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT", "ENTPD1"],
    "NK cells": ["NCAM1", "KLRD1", "NKG7", "GNLY", "FCER1G", "GZMB", "PRF1", "FCGR3A", "KLRB1"],
    "Monocytes": ["CD14", "LYZ", "CSF1R", "FCGR3A", "MS4A7", "CD68", "S100A9", "S100A8", "VCAN"],
    "Macrophages M0": ["CD68", "MRC1", "CD163", "CSF1R", "ADGRE1", "HLA-DRA"],
    "Macrophages M1": ["CD68", "CD86", "HLA-DRA", "IL1B", "TNF", "CXCL9", "CXCL10", "IDO1"],
    "Macrophages M2": ["CD68", "MRC1", "CD163", "MSR1", "IL10", "CCL18", "TGM2", "TGFB1"],
    "Dendritic cells": ["ITGAX", "HLA-DRA", "CD1C", "FCER1A", "CLEC9A", "XCR1", "SIGLEC6"],
    "Neutrophils": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "SELL", "MMP9", "ELANE"],
}


def ssgsea_score(expr_vector, gene_set, all_genes):
    """Compute ssGSEA enrichment score for one sample."""
    gene_set = [g for g in gene_set if g in all_genes]
    if len(gene_set) < 3:
        return np.nan

    gene_idx = {g: i for i, g in enumerate(all_genes)}
    ranked = np.argsort(expr_vector)[::-1]
    gene_set_idx = set(gene_idx[g] for g in gene_set if g in gene_idx)

    n = len(expr_vector)
    n_set = len(gene_set_idx)
    n_not = n - n_set

    if n_set == 0 or n_not == 0:
        return np.nan

    # Weights: rank-based
    es_up = 0
    es_down = 0
    running_sum = 0
    max_es = 0
    min_es = 0
    in_count = 0
    out_count = 0

    for rank_pos, gene_pos in enumerate(ranked):
        if gene_pos in gene_set_idx:
            in_count += 1
            running_sum += in_count / n_set
        else:
            out_count += 1
            running_sum -= out_count / n_not

        if running_sum > max_es:
            max_es = running_sum
        if running_sum < min_es:
            min_es = running_sum

    es = max_es if abs(max_es) >= abs(min_es) else min_es
    return es


def compute_ssgsea_scores(expr_df):
    """Compute ssGSEA scores for all immune cell types across all samples."""
    all_genes = list(expr_df.index)
    scores = {}

    for cell_type, gene_set in IMMUNE_SIGNATURES.items():
        cell_scores = []
        for sample in expr_df.columns:
            expr_vector = expr_df[sample].values.astype(float)
            score = ssgsea_score(expr_vector, gene_set, all_genes)
            cell_scores.append(score)
        scores[cell_type] = cell_scores

    return pd.DataFrame(scores, index=expr_df.columns)


def compare_groups(scores_df, groups):
    """Compare immune cell type scores between disease and control."""
    disease_samples = [s for s, g in groups.items() if g == 'disease']
    control_samples = [s for s, g in groups.items() if g == 'control']

    disease_samples = [s for s in disease_samples if s in scores_df.index]
    control_samples = [s for s in control_samples if s in scores_df.index]

    results = []
    for cell_type in scores_df.columns:
        d_vals = scores_df.loc[disease_samples, cell_type].dropna()
        c_vals = scores_df.loc[control_samples, cell_type].dropna()

        if len(d_vals) < 3 or len(c_vals) < 3:
            continue

        stat, pval = stats.mannwhitneyu(d_vals, c_vals, alternative='two-sided')

        mean_d = d_vals.mean()
        mean_c = c_vals.mean()
        pooled_std = np.sqrt((d_vals.std()**2 + c_vals.std()**2) / 2)
        cohens_d = (mean_d - mean_c) / pooled_std if pooled_std > 0 else 0

        results.append({
            'cell_type': cell_type,
            'mean_disease': mean_d,
            'mean_control': mean_c,
            'cohens_d': cohens_d,
            'U_statistic': stat,
            'p_value': pval,
            'n_disease': len(d_vals),
            'n_control': len(c_vals),
        })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        _, padj, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
        results_df['padj'] = padj
        results_df['significant'] = results_df['padj'] < 0.05

    return results_df


def create_deconvolution_figures(scores_df, stats_df, accession, disease, groups):
    """Create visualization figures."""
    disease_samples = [s for s, g in groups.items() if g == 'disease' and s in scores_df.index]
    control_samples = [s for s, g in groups.items() if g == 'control' and s in scores_df.index]

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Panel 1: Box plot of significant cell types
    ax = axes[0]
    sig_types = stats_df[stats_df['significant']]['cell_type'].tolist() if 'significant' in stats_df.columns else []
    if not sig_types:
        sig_types = stats_df.nsmallest(5, 'p_value')['cell_type'].tolist() if len(stats_df) > 0 else list(scores_df.columns)[:5]

    plot_data = []
    for cell_type in sig_types[:8]:
        for sample in disease_samples:
            plot_data.append({'cell_type': cell_type, 'score': scores_df.loc[sample, cell_type], 'group': 'Disease'})
        for sample in control_samples:
            plot_data.append({'cell_type': cell_type, 'score': scores_df.loc[sample, cell_type], 'group': 'Control'})

    if plot_data:
        plot_df = pd.DataFrame(plot_data)
        sns.boxplot(data=plot_df, x='cell_type', y='score', hue='group',
                    palette={'Disease': '#e74c3c', 'Control': '#3498db'}, ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
        ax.set_title(f'{accession} ({disease}) - Immune Cell Deconvolution')

    # Panel 2: Effect sizes
    ax2 = axes[1]
    if len(stats_df) > 0 and 'cohens_d' in stats_df.columns:
        stats_sorted = stats_df.sort_values('cohens_d')
        colors = ['#e74c3c' if d > 0 else '#3498db' for d in stats_sorted['cohens_d']]
        ax2.barh(range(len(stats_sorted)), stats_sorted['cohens_d'], color=colors, alpha=0.7)
        ax2.set_yticks(range(len(stats_sorted)))
        ax2.set_yticklabels(stats_sorted['cell_type'], fontsize=9)
        ax2.set_xlabel("Cohen's d")
        ax2.axvline(0, color='black', linewidth=0.5)
        ax2.set_title('Effect sizes (disease vs. control)')

    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, f'deconvolution_{accession}.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()


def load_expression_and_groups(accession):
    """Load expression matrix and sample groups."""
    expr_path = os.path.join(PROCESSED_DIR, f"{accession}_normalized.csv")
    meta_path = os.path.join(METADATA_DIR, f"{accession}_metadata.csv")

    if not os.path.exists(expr_path):
        print(f"  No expression file: {expr_path}")
        return None, None

    expr_df = pd.read_csv(expr_path, index_col=0)

    groups = {}
    if os.path.exists(meta_path):
        meta_df = pd.read_csv(meta_path, index_col=0)
        if 'group' in meta_df.columns:
            for sample in meta_df.index:
                groups[str(sample)] = meta_df.loc[sample, 'group']

    if not groups:
        # Fallback: split samples by position
        n = len(expr_df.columns)
        half = n // 2
        for i, sample in enumerate(expr_df.columns):
            groups[str(sample)] = 'disease' if i < half else 'control'

    return expr_df, groups


def main():
    print("Phase 8: Immune Cell Deconvolution")
    print("=" * 60)

    datasets = [
        ("GSE25101", "AS"),
        ("GSE73754", "AS"),
        ("GSE18781", "AS"),
        ("GSE221786", "AS"),
        ("GSE61281", "PsA"),
        ("GSE59071", "IBD"),
        ("GSE58667", "jSpA"),
    ]

    all_stats = []

    for accession, disease in datasets:
        print(f"\nProcessing {accession} ({disease})...")

        expr_df, groups = load_expression_and_groups(accession)
        if expr_df is None:
            continue

        print(f"  Expression: {expr_df.shape}, Groups: {len(groups)}")

        # Compute ssGSEA scores
        print("  Computing ssGSEA scores...")
        scores_df = compute_ssgsea_scores(expr_df)

        scores_path = os.path.join(DECONV_DIR, f"{accession}_immune_scores.csv")
        scores_df.to_csv(scores_path)

        # Compare groups
        stats_df = compare_groups(scores_df, groups)

        if len(stats_df) > 0:
            stats_path = os.path.join(DECONV_DIR, f"{accession}_deconvolution_stats.csv")
            stats_df.to_csv(stats_path, index=False)

            n_sig = stats_df['significant'].sum() if 'significant' in stats_df.columns else 0
            print(f"  Significant cell types: {n_sig}/13")

            for _, row in stats_df[stats_df.get('significant', pd.Series([False]*len(stats_df)))].iterrows():
                print(f"    {row['cell_type']:30s} d={row['cohens_d']:+.3f}  padj={row['padj']:.3e}")

            all_stats.append({'accession': accession, 'disease': disease,
                              'n_significant': int(n_sig)})

        create_deconvolution_figures(scores_df, stats_df, accession, disease, groups)

    print("\nPhase 8 complete.")


if __name__ == "__main__":
    main()
