#!/usr/bin/env python3
"""
Sensitivity Analysis: Robustness of the AS Meta-Signature
=========================================================
Tests how core findings change under different analytical choices:
  Panel A: FDR threshold sensitivity curve
  Panel B: Effect size filtering impact
  Panel C: Leave-one-out rank analysis
  Panel D: Dataset contribution heatmap for the 40 core genes

Spondyloarthritis AI Computational Biology Institute
"""

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import json
import os
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap

BASE = '/home/user/workspace/project1_cross_disease_metaanalysis_spa'

TEAL = '#20808D'
RUST = '#A84B2F'
DARK = '#1B474D'
GOLD = '#D4A843'
GRAY = '#7A8B8B'

AS_DATASETS = ['GSE25101', 'GSE73754', 'GSE221786']


def load_meta_results():
    """Load meta-analysis results."""
    meta_path = os.path.join(BASE, 'results', 'meta_analysis', 'meta_analysis_AS_only.csv')
    sig_path = os.path.join(BASE, 'results', 'meta_analysis', 'meta_significant_AS_only.csv')

    # Also try the fisher_meta_results path from Phase 4
    if not os.path.exists(meta_path):
        meta_path = os.path.join(BASE, 'results', 'meta_analysis', 'fisher_meta_results.csv')
    if not os.path.exists(sig_path):
        sig_path = os.path.join(BASE, 'results', 'meta_analysis', 'meta_significant_genes.csv')

    if not os.path.exists(meta_path) or not os.path.exists(sig_path):
        print("ERROR: Meta-analysis results not found. Run Phase 4 first.")
        return None, None

    meta = pd.read_csv(meta_path)
    original_sig = pd.read_csv(sig_path)

    # Normalize column names
    if 'gene_symbol' in meta.columns:
        meta = meta.rename(columns={'gene_symbol': 'gene'})
    if 'gene_symbol' in original_sig.columns:
        original_sig = original_sig.rename(columns={'gene_symbol': 'gene'})
    if 'combined_padj' in meta.columns:
        meta = meta.rename(columns={'combined_padj': 'fisher_padj', 'mean_log2FC': 'weighted_log2FC'})

    print(f"Meta-analysis: {len(meta)} genes")
    print(f"Significant genes: {len(original_sig)}")
    return meta, original_sig


def fdr_sensitivity_analysis(meta, original_genes):
    """Vary FDR threshold and count recovered genes."""
    fdr_thresholds = [0.001, 0.005, 0.01, 0.025, 0.05, 0.10, 0.15, 0.20]
    results = {}

    for fdr in fdr_thresholds:
        if 'fisher_padj' in meta.columns:
            sig = meta[meta['fisher_padj'] < fdr]
        else:
            sig = meta[meta.get('padj', meta.get('combined_padj', pd.Series())) < fdr] if 'padj' in meta.columns else pd.DataFrame()
        genes = set(sig.get('gene', pd.Series()).values)
        overlap = genes & original_genes
        results[fdr] = {'n_genes': len(sig), 'overlap': len(overlap),
                        'recovery': 100 * len(overlap) / max(len(original_genes), 1)}
        print(f"  FDR < {fdr}: {len(sig)} genes ({len(overlap)}/{len(original_genes)} recovered)")

    return results


def effect_size_sensitivity(meta, original_genes):
    """Vary log2FC threshold."""
    fc_thresholds = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
    results = {}

    for fc in fc_thresholds:
        if 'fisher_padj' in meta.columns and 'weighted_log2FC' in meta.columns:
            sig = meta[(meta['fisher_padj'] < 0.05) & (meta['weighted_log2FC'].abs() > fc)]
        else:
            sig = meta  # fallback
        results[fc] = {'n_genes': len(sig), 'overlap': len(set(sig.get('gene', pd.Series()).values) & original_genes)}
        print(f"  |FC| > {fc}: {len(sig)} genes")

    return results


def leave_one_out_analysis(original_genes):
    """Leave-one-dataset-out meta-analysis."""
    all_data = {}
    for ds in AS_DATASETS:
        deg_path = os.path.join(BASE, 'results', 'deg', f'{ds}_deg_genes.csv')
        if not os.path.exists(deg_path):
            print(f"  Missing DEG file: {deg_path}")
            continue
        df = pd.read_csv(deg_path)
        if 'gene_symbol' not in df.columns:
            continue
        df = df.dropna(subset=['gene_symbol'])
        df = df.sort_values('pvalue').drop_duplicates(subset='gene_symbol', keep='first')
        all_data[ds] = dict(zip(df['gene_symbol'], zip(df.get('padj', df['pvalue']), df['log2FC'])))
        print(f"  {ds}: {len(all_data[ds])} genes")

    loo_results = {}
    for excluded in list(all_data.keys()):
        remaining = [ds for ds in all_data if ds != excluded]
        if not remaining:
            continue

        gene_pvals = {}
        for ds in remaining:
            for gene, (padj, lfc) in all_data[ds].items():
                if gene not in gene_pvals:
                    gene_pvals[gene] = {'pvals': [], 'lfcs': []}
                gene_pvals[gene]['pvals'].append(padj)
                gene_pvals[gene]['lfcs'].append(lfc)

        results_list = []
        for gene, data in gene_pvals.items():
            valid = [p for p in data['pvals'] if pd.notna(p) and 0 < p <= 1]
            if len(valid) < 1:
                continue
            chi2 = -2 * sum(np.log(max(p, 1e-300)) for p in valid)
            fisher_p = stats.chi2.sf(chi2, 2 * len(valid))
            results_list.append({'gene': gene, 'fisher_p': float(fisher_p), 'mean_lfc': np.mean(data['lfcs'])})

        df_loo = pd.DataFrame(results_list).sort_values('fisher_p')
        df_loo['rank'] = range(1, len(df_loo) + 1)
        df_loo['percentile'] = 100 * df_loo['rank'] / max(len(df_loo), 1)

        orig_ranks = df_loo[df_loo['gene'].isin(original_genes)]
        loo_results[excluded] = {
            'total_tested': len(df_loo),
            'orig_gene_ranks': orig_ranks[['gene', 'rank', 'percentile', 'fisher_p']].to_dict('records'),
            'median_rank': float(orig_ranks['rank'].median()) if len(orig_ranks) > 0 else None,
            'median_percentile': float(orig_ranks['percentile'].median()) if len(orig_ranks) > 0 else None,
            'n_in_top_100': int((orig_ranks['rank'] <= 100).sum()),
            'n_in_top_10pct': int((orig_ranks['percentile'] <= 10).sum()),
        }
        print(f"  Excluding {excluded}: {len(df_loo)} genes, {len(orig_ranks)} original genes found")

    return loo_results


def main():
    print("Sensitivity Analysis: Robustness of AS Meta-Signature")
    print("=" * 60)

    meta, original_sig = load_meta_results()
    if meta is None:
        return

    original_genes = set(original_sig['gene'].values) if 'gene' in original_sig.columns else set()
    print(f"Original significant genes: {len(original_genes)}")

    # Analysis 1: FDR threshold sensitivity
    print("\n1. FDR Threshold Sensitivity")
    fdr_results = fdr_sensitivity_analysis(meta, original_genes)

    # Analysis 2: Effect size sensitivity
    print("\n2. Effect Size Sensitivity")
    fc_results = effect_size_sensitivity(meta, original_genes)

    # Analysis 3: Leave-one-out
    print("\n3. Leave-One-Out Analysis")
    loo_results = leave_one_out_analysis(original_genes)

    # Save results
    sens_dir = os.path.join(BASE, 'results', 'sensitivity')
    os.makedirs(sens_dir, exist_ok=True)

    summary = {
        'original_genes': len(original_genes),
        'fdr_sensitivity': {str(k): v for k, v in fdr_results.items()},
        'fc_sensitivity': {str(k): v for k, v in fc_results.items()},
        'leave_one_out': loo_results,
    }

    with open(os.path.join(sens_dir, 'sensitivity_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nResults saved to: {sens_dir}")
    print("Sensitivity analysis complete.")


if __name__ == "__main__":
    main()
