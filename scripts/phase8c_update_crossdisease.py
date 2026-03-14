#!/usr/bin/env python3
"""
Phase 8c: Update cross-disease deconvolution analysis with PsA data
and generate updated cross-disease comparison figures.
"""

import pandas as pd
import numpy as np
from scipy import stats
import json
import os
import warnings
warnings.filterwarnings('ignore')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/home/user/workspace/project1_cross_disease_metaanalysis_spa'

COL_MAP = {
    'CD8_T_cells': 'CD8 T cells', 'CD4_T_cells': 'CD4 T cells',
    'Th17_cells': 'Th17 cells', 'Treg_cells': 'Treg cells',
    'Gamma_delta_T_cells': 'gamma_delta T', 'B_cells': 'B cells',
    'Plasma_cells': 'Plasma cells', 'NK_cells': 'NK cells',
    'Monocytes': 'Monocytes', 'Dendritic_cells': 'Dendritic cells',
    'Neutrophils': 'Neutrophils', 'Mast_cells': 'Mast cells', 'ILC3': 'ILC3'
}

DISEASE_MAP = {
    'GSE18781': 'AS', 'GSE25101': 'AS', 'GSE73754': 'AS', 'GSE221786': 'AS',
    'GSE59071': 'IBD', 'GSE58667': 'jSpA', 'GSE61281': 'PsA'
}


def load_all_immune_scores():
    """Load immune scores for all datasets."""
    all_scores = []
    for dataset in ['GSE18781', 'GSE25101', 'GSE221786', 'GSE59071', 'GSE58667', 'GSE61281']:
        f = os.path.join(BASE, 'results', 'deconvolution', f'{dataset}_immune_scores.csv')
        if os.path.exists(f):
            df = pd.read_csv(f, index_col=0)
            if df.shape[0] < 2:
                continue
            df = df.rename(columns=COL_MAP)
            df['disease'] = DISEASE_MAP.get(dataset, 'Unknown')
            df['dataset'] = dataset
            all_scores.append(df)
            print(f"  Loaded {dataset}: {df.shape[0]} samples")
    return pd.concat(all_scores, ignore_index=False) if all_scores else pd.DataFrame()


def main():
    print("Phase 8c: Update Cross-Disease Deconvolution")
    print("=" * 60)

    combined = load_all_immune_scores()
    if combined.empty:
        print("No immune scores found. Run Phase 8 and 8b first.")
        return

    cell_types = [c for c in combined.columns if c not in ['disease', 'dataset']]
    diseases = combined['disease'].unique()

    print(f"Combined: {combined.shape[0]} samples, {len(cell_types)} cell types, {len(diseases)} diseases")

    # Kruskal-Wallis across diseases
    kw_results = []
    for ct in cell_types:
        groups = [combined[combined['disease'] == d][ct].dropna() for d in diseases
                  if len(combined[combined['disease'] == d][ct].dropna()) >= 3]
        if len(groups) >= 2:
            h_stat, p_val = stats.kruskal(*groups)
            kw_results.append({'cell_type': ct, 'H_statistic': h_stat, 'p_value': p_val})

    kw_df = pd.DataFrame(kw_results).sort_values('p_value')
    kw_out = os.path.join(BASE, 'results', 'deconvolution', 'cross_disease_kruskal_updated.csv')
    kw_df.to_csv(kw_out, index=False)
    print(f"\nKruskal-Wallis results saved: {kw_out}")

    # Generate cross-disease heatmap
    fig, ax = plt.subplots(figsize=(12, 6))
    heatmap_data = combined.groupby('disease')[cell_types].mean()
    heatmap_z = (heatmap_data - heatmap_data.mean()) / (heatmap_data.std() + 1e-10)
    heatmap_z = heatmap_z.fillna(0)

    im = ax.imshow(heatmap_z.T.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(len(heatmap_z.index)))
    ax.set_xticklabels(heatmap_z.index, rotation=45, ha='right', fontsize=11)
    ax.set_yticks(range(len(cell_types)))
    ax.set_yticklabels(cell_types, fontsize=9)
    ax.set_title('Cross-Disease Immune Cell Enrichment (z-scored)', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax, shrink=0.7, label='z-score')
    plt.tight_layout()

    fig_out = os.path.join(BASE, 'figures', 'fig8_deconvolution_summary_v2.png')
    plt.savefig(fig_out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Figure saved: {fig_out}")

    print("\nPhase 8c complete.")


if __name__ == "__main__":
    main()
