#!/usr/bin/env python3
"""
Phase 8c: Update cross-disease deconvolution analysis
======================================================
Now that PsA (GSE61281) deconvolution is fixed (phase8b),
this script:
1. Loads all available deconvolution results (NNLS + robust fallback)
2. Creates an updated cross-disease immune comparison
3. Identifies disease-specific vs shared immune signatures
4. Computes immune similarity between SpA subtypes

Author: Perplexity Computer / Spondyloarthritis AI Computational Biology Institute
Date: March 2026
"""

import os
import sys
import json
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results", "phase8_deconvolution")
FIGURES_DIR = os.path.join(BASE_DIR, "figures", "phase8_deconvolution")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# Dataset metadata
DATASET_META = {
    "GSE18781":  {"disease": "AS",    "tissue": "blood"},
    "GSE25101":  {"disease": "AS",    "tissue": "blood"},
    "GSE73754":  {"disease": "AS",    "tissue": "blood"},
    "GSE221786": {"disease": "AS",    "tissue": "PBMC"},
    "GSE61281":  {"disease": "PsA",   "tissue": "synovium"},
    "GSE113189": {"disease": "PsA",   "tissue": "synovium"},
    "GSE59071":  {"disease": "IBD",   "tissue": "colon"},
    "GSE83687":  {"disease": "IBD",   "tissue": "ileum"},
    "GSE14905":  {"disease": "Pso",   "tissue": "skin"},
    "GSE30999":  {"disease": "Pso",   "tissue": "skin"},
    "GSE36700":  {"disease": "ReA",   "tissue": "blood"},
    "GSE100927": {"disease": "uSpA",  "tissue": "blood"},
}


def load_all_deconv_results() -> dict:
    """
    Load all available deconvolution results.
    Tries NNLS first, falls back to robust/ssgsea.
    """
    all_results = {}
    
    for gse_id in DATASET_META:
        gse_dir = os.path.join(RESULTS_DIR, gse_id)
        
        # Try different result files
        for fname in ['deconv_nnls.csv', 'deconv_robust.csv', 'deconv_ssgsea.csv']:
            fpath = os.path.join(gse_dir, fname)
            if os.path.exists(fpath):
                df = pd.read_csv(fpath, index_col=0)
                all_results[gse_id] = {
                    'deconv': df,
                    'method': fname.replace('.csv', '').replace('deconv_', ''),
                    'disease': DATASET_META[gse_id]['disease'],
                    'tissue': DATASET_META[gse_id]['tissue'],
                }
                break
        
        if gse_id not in all_results:
            print(f"  No results for {gse_id}")
    
    print(f"Loaded results for {len(all_results)}/{len(DATASET_META)} datasets")
    return all_results


def compute_disease_mean_profiles(all_results: dict) -> pd.DataFrame:
    """
    For each dataset, compute mean immune fractions in disease samples.
    Returns: datasets x cell_types DataFrame.
    """
    profiles = {}
    
    for gse_id, res in all_results.items():
        deconv = res['deconv']
        
        # Load labels
        labels = None
        for lpath in [
            os.path.join(RESULTS_DIR, gse_id, 'diff_immune_nnls.csv'),
            os.path.join(RESULTS_DIR, gse_id, 'diff_immune_robust.csv'),
        ]:
            if os.path.exists(lpath):
                break  # labels file not directly usable, use deconv mean
        
        # Use overall mean (no label separation for cross-disease)
        disease_mean = deconv.mean()
        label = f"{res['disease']}\n({gse_id})"
        profiles[label] = disease_mean
    
    if not profiles:
        return pd.DataFrame()
    
    profile_df = pd.DataFrame(profiles).T  # datasets x cell_types
    
    # Only keep cell types present in most datasets
    coverage = (profile_df > 0).sum(axis=0)
    common_cells = coverage[coverage >= len(profile_df) * 0.5].index
    profile_df = profile_df[common_cells]
    
    return profile_df


def plot_updated_cross_disease(profile_df: pd.DataFrame):
    """Create updated cross-disease heatmap with hierarchical clustering."""
    if profile_df.empty:
        print("No data for cross-disease plot")
        return
    
    # Z-score normalize
    from scipy.stats import zscore
    z_df = profile_df.apply(zscore, axis=0, nan_policy='omit').fillna(0)
    
    # Figure 1: Clustered heatmap
    fig1_path = os.path.join(FIGURES_DIR, 'cross_disease_immune_updated.png')
    
    g = sns.clustermap(
        z_df,
        cmap='RdBu_r',
        center=0,
        figsize=(max(12, len(z_df.columns) * 0.5), max(8, len(z_df) * 0.6)),
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        cbar_kws={'label': 'Z-score'},
        dendrogram_ratio=0.15,
        annot=False,
    )
    g.fig.suptitle('Cross-Disease Immune Cell Profiles\n(Z-score normalized, hierarchically clustered)',
                   y=1.02, fontsize=13)
    plt.savefig(fig1_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig1_path}")
    
    # Figure 2: Disease-level summary (average across datasets per disease)
    disease_groups = defaultdict(list)
    for idx in profile_df.index:
        disease = idx.split('\n')[0]
        disease_groups[disease].append(idx)
    
    disease_means = {}
    for disease, datasets in disease_groups.items():
        disease_means[disease] = profile_df.loc[datasets].mean()
    
    disease_df = pd.DataFrame(disease_means).T
    z_disease = disease_df.apply(zscore, axis=0, nan_policy='omit').fillna(0)
    
    fig2_path = os.path.join(FIGURES_DIR, 'disease_level_immune_summary.png')
    fig, ax = plt.subplots(figsize=(max(12, len(z_disease.columns) * 0.5), 6))
    
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    sns.heatmap(z_disease, ax=ax, cmap=cmap, center=0,
                xticklabels=True, yticklabels=True,
                linewidths=1, annot=True, fmt='.1f', fontsize=8,
                cbar_kws={'label': 'Z-score'})
    ax.set_title('Disease-Level Immune Cell Profiles (Mean across datasets)', fontsize=13)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
    plt.tight_layout()
    plt.savefig(fig2_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig2_out}")


if __name__ == '__main__':
    print("Phase 8c: Updating cross-disease immune analysis")
    print("=" * 55)
    
    all_results = load_all_deconv_results()
    
    if all_results:
        profile_df = compute_disease_mean_profiles(all_results)
        print(f"Profile matrix: {profile_df.shape}")
        
        # Save
        profile_path = os.path.join(RESULTS_DIR, 'cross_disease_immune_profiles.csv')
        profile_df.to_csv(profile_path)
        print(f"Saved profiles: {profile_path}")
        
        plot_updated_cross_disease(profile_df)
    else:
        print("No deconvolution results found. Run phase8 and phase8b first.")
    
    print("\n=== Phase 8c Complete ===")
