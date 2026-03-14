#!/usr/bin/env python3
"""
Phase 8b: Fix PsA (GSE61281) Immune Cell Deconvolution
Re-runs deconvolution with corrected sample group assignments.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
DECONV_DIR = os.path.join(RESULTS_DIR, "deconvolution")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
METADATA_DIR = os.path.join(BASE_DIR, "data", "metadata")

# GSE61281 cell type gene signatures (PsA-relevant)
IMMUNE_SIGNATURES = {
    "CD8 T cells": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7", "GZMA", "GNLY", "EOMES"],
    "CD4 T cells": ["CD3D", "CD3E", "CD4", "IL7R", "TCF7", "CCR7", "SELL", "LEF1", "FOXP3"],
    "Th17 cells": ["IL17A", "IL17F", "RORC", "IL23R", "CCR6", "KLRB1", "CD3D", "RORA", "AHR"],
    "Treg cells": ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18", "TIGIT", "IL10", "ENTPD1"],
    "gamma_delta T": ["TRDC", "TRDV1", "TRDV2", "TRGC1", "TRGV9", "IL17A", "ZBTB16"],
    "B cells": ["CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "IGHM", "IGHD", "TCL1A"],
    "Plasma cells": ["IGHG1", "IGHG2", "JCHAIN", "TNFRSF17", "SDC1", "CD38", "CD27", "MZB1"],
    "NK cells": ["NCAM1", "KLRD1", "NKG7", "GNLY", "FCER1G", "GZMB", "PRF1", "FCGR3A"],
    "Monocytes": ["CD14", "LYZ", "CSF1R", "FCGR3A", "CD68", "S100A9", "S100A8", "VCAN"],
    "Dendritic cells": ["ITGAX", "HLA-DRA", "CD1C", "FCER1A", "CLEC9A", "XCR1", "SIGLEC6"],
    "Neutrophils": ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "SELL", "MMP9", "ELANE"],
    "Mast cells": ["TPSAB1", "CPA3", "GATA2", "KIT", "HDC", "SRGN", "MS4A2"],
    "ILC3": ["RORC", "IL22", "NCR2", "CXCR6", "IL1R1", "KIT", "AHR"],
}


def ssgsea_score(expr_vector, gene_set, all_genes):
    """Compute ssGSEA enrichment score."""
    gene_set = [g for g in gene_set if g in all_genes]
    if len(gene_set) < 3:
        return np.nan

    gene_idx = {g: i for i, g in enumerate(all_genes)}
    n = len(expr_vector)
    n_set = len(gene_set)
    n_not = n - n_set

    if n_set == 0 or n_not == 0:
        return np.nan

    ranked = np.argsort(expr_vector)[::-1]
    gene_set_idx = set(gene_idx[g] for g in gene_set if g in gene_idx)

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
        max_es = max(max_es, running_sum)
        min_es = min(min_es, running_sum)

    return max_es if abs(max_es) >= abs(min_es) else min_es


def load_gse61281_with_correct_groups():
    """Load GSE61281 with corrected PsA vs Control group assignments."""
    expr_path = os.path.join(PROCESSED_DIR, "GSE61281_normalized.csv")
    if not os.path.exists(expr_path):
        print(f"  Expression file not found: {expr_path}")
        return None, None

    expr_df = pd.read_csv(expr_path, index_col=0)
    print(f"  GSE61281 expression: {expr_df.shape}")

    # Load metadata
    meta_path = os.path.join(METADATA_DIR, "GSE61281_metadata.csv")
    groups = {}
    if os.path.exists(meta_path):
        meta_df = pd.read_csv(meta_path, index_col=0)
        if 'group' in meta_df.columns:
            groups = dict(zip(meta_df.index.astype(str), meta_df['group']))
            print(f"  Groups from metadata: {pd.Series(groups).value_counts().to_dict()}")

    if not groups:
        # Manual fallback based on GSE61281 design (20 PsA, 20 controls)
        samples = list(expr_df.columns)
        for i, s in enumerate(samples):
            groups[str(s)] = 'disease' if i < 20 else 'control'
        print(f"  Using positional group assignment")

    return expr_df, groups


def compute_ssgsea_scores(expr_df):
    """Compute ssGSEA scores for all immune cell types."""
    all_genes = list(expr_df.index)
    scores = {}

    for cell_type, gene_set in IMMUNE_SIGNATURES.items():
        cell_scores = []
        for sample in expr_df.columns:
            score = ssgsea_score(expr_df[sample].values.astype(float), gene_set, all_genes)
            cell_scores.append(score)
        scores[cell_type] = cell_scores
        n_found = sum(1 for g in gene_set if g in all_genes)
        print(f"  {cell_type}: {n_found}/{len(gene_set)} signature genes found")

    return pd.DataFrame(scores, index=expr_df.columns)


def statistical_comparison(scores_df, groups):
    """Wilcoxon rank-sum test with BH correction."""
    disease_samples = [s for s, g in groups.items() if g == 'disease' and s in scores_df.index]
    control_samples = [s for s, g in groups.items() if g == 'control' and s in scores_df.index]

    print(f"  Disease: {len(disease_samples)}, Control: {len(control_samples)}")

    results = []
    for cell_type in scores_df.columns:
        d_vals = scores_df.loc[disease_samples, cell_type].dropna()
        c_vals = scores_df.loc[control_samples, cell_type].dropna()

        if len(d_vals) < 3 or len(c_vals) < 3:
            continue

        stat, pval = stats.mannwhitneyu(d_vals, c_vals, alternative='two-sided')
        pooled_std = np.sqrt((d_vals.std()**2 + c_vals.std()**2) / 2)
        cohens_d = (d_vals.mean() - c_vals.mean()) / pooled_std if pooled_std > 0 else 0

        results.append({
            'cell_type': cell_type,
            'mean_disease': d_vals.mean(),
            'mean_control': c_vals.mean(),
            'cohens_d': cohens_d,
            'U_statistic': stat,
            'p_value': pval,
            'n_disease': len(d_vals),
            'n_control': len(c_vals),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, padj, _, _ = multipletests(df['p_value'], method='fdr_bh')
        df['padj'] = padj
        df['significant'] = df['padj'] < 0.05

    return df


def create_figure(scores_df, stats_df, groups):
    """Create PsA deconvolution summary figure."""
    disease_samples = [s for s, g in groups.items() if g == 'disease' and s in scores_df.index]
    control_samples = [s for s, g in groups.items() if g == 'control' and s in scores_df.index]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Effect sizes
    ax = axes[0]
    if len(stats_df) > 0:
        stats_sorted = stats_df.sort_values('cohens_d')
        colors = ['#d32f2f' if d > 0 else '#1565c0' for d in stats_sorted['cohens_d']]
        ax.barh(range(len(stats_sorted)), stats_sorted['cohens_d'], color=colors, alpha=0.8, edgecolor='white')
        ax.set_yticks(range(len(stats_sorted)))
        ax.set_yticklabels(stats_sorted['cell_type'], fontsize=10)
        ax.set_xlabel("Cohen's d (PsA vs Control)", fontsize=12)
        ax.axvline(0, color='black', linewidth=0.5)
        ax.set_title("GSE61281 (PsA)\nImmune Cell Effect Sizes", fontsize=12, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        for i, row in enumerate(stats_sorted.itertuples()):
            note = f"padj={row.padj:.2f}" if hasattr(row, 'padj') else ""
            ax.text(row.cohens_d + 0.01 * np.sign(row.cohens_d), i, f'{row.cohens_d:.2f}',
                    va='center', fontsize=8)

    # Panel B: Box plots
    ax2 = axes[1]
    top_cells = stats_df.reindex(stats_df['p_value'].nsmallest(6).index)['cell_type'].tolist() if len(stats_df) > 0 else []
    if top_cells:
        plot_data = []
        for ct in top_cells:
            for s in disease_samples:
                plot_data.append({'cell_type': ct, 'score': scores_df.loc[s, ct], 'group': 'PsA'})
            for s in control_samples:
                plot_data.append({'cell_type': ct, 'score': scores_df.loc[s, ct], 'group': 'Control'})
        plot_df = pd.DataFrame(plot_data)
        import seaborn as sns
        sns.boxplot(data=plot_df, x='cell_type', y='score', hue='group',
                    palette={'PsA': '#d32f2f', 'Control': '#1565c0'}, ax=ax2, fliersize=3)
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=40, ha='right', fontsize=9)
        ax2.set_title('Top Cell Types by p-value', fontsize=12, fontweight='bold')
        ax2.set_xlabel('')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

    plt.suptitle('GSE61281 (PsA) - Immune Cell Deconvolution Results', fontsize=13)
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, 'deconvolution_GSE61281_PsA_fixed.png')
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {out_path}")


def main():
    print("Phase 8b: Fix PsA Deconvolution (GSE61281)")
    print("=" * 60)

    expr_df, groups = load_gse61281_with_correct_groups()
    if expr_df is None:
        print("Cannot proceed without expression data")
        return

    print("\nComputing ssGSEA scores...")
    scores_df = compute_ssgsea_scores(expr_df)

    scores_path = os.path.join(DECONV_DIR, "GSE61281_immune_scores.csv")
    scores_df.to_csv(scores_path)
    print(f"Scores saved: {scores_path}")

    print("\nStatistical comparison (Wilcoxon + BH)...")
    stats_df = statistical_comparison(scores_df, groups)

    if len(stats_df) > 0:
        stats_path = os.path.join(DECONV_DIR, "PsA_deconvolution_stats.csv")
        stats_df.to_csv(stats_path, index=False)
        print(f"Stats saved: {stats_path}")

        print("\nResults:")
        print(stats_df[['cell_type', 'cohens_d', 'p_value', 'padj', 'significant']].to_string(index=False))

    create_figure(scores_df, stats_df, groups)
    print("\nPhase 8b complete.")


if __name__ == "__main__":
    main()
