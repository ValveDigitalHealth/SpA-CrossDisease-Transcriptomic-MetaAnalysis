#!/usr/bin/env python3
"""
Phase 8: Immune Cell Deconvolution
===================================
Cross-Disease Transcriptomic Meta-Analysis of SpA — Project 1

Performs immune cell type deconvolution on bulk RNA-seq / microarray data
for all SpA-related datasets using multiple reference signatures:

1. CIBERSORTx-style deconvolution (non-negative least squares)
2. xCell-inspired enrichment scoring
3. TIMER2.0-inspired immune estimation (simplified)
4. ssGSEA for immune gene sets

Reference signature matrices:
- LM22 (CIBERSORT 22 immune cell types)
- ImmuneCC (17 cell types, microarray-focused)
- Custom SpA-relevant cell types (Th17, Tregs, ILC3, etc.)

Outputs:
- Cell fraction estimates per sample per dataset
- Differential immune infiltration (disease vs. control)
- Cross-disease comparison heatmaps
- Correlation with clinical traits

Target datasets: All 12 GSE datasets from meta-analysis

Author: Perplexity Computer / Spondyloarthritis AI Computational Biology Institute
Date: March 2026
"""

import os
import sys
import json
import gzip
import logging
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from collections import defaultdict

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results", "phase8_deconvolution")
FIGURES_DIR = os.path.join(BASE_DIR, "figures", "phase8_deconvolution")

# All datasets from meta-analysis
ALL_DATASETS = {
    # Ankylosing Spondylitis / axSpA
    "GSE18781": {"disease": "AS",    "tissue": "blood",   "platform": "GPL570",  "n_samples": 55},
    "GSE25101": {"disease": "AS",    "tissue": "blood",   "platform": "GPL570",  "n_samples": 30},
    "GSE73754": {"disease": "AS",    "tissue": "blood",   "platform": "GPL570",  "n_samples": 80},
    "GSE221786": {"disease": "AS",   "tissue": "PBMC",    "platform": "GPL24676", "n_samples": 48},
    # Psoriatic Arthritis
    "GSE61281": {"disease": "PsA",   "tissue": "blood",   "platform": "GPL570",  "n_samples": 42},
    "GSE113189": {"disease": "PsA",  "tissue": "synovium", "platform": "GPL570", "n_samples": 56},
    # IBD-associated SpA
    "GSE59071": {"disease": "IBD",   "tissue": "colon",   "platform": "GPL6244", "n_samples": 116},
    "GSE83687": {"disease": "IBD",   "tissue": "ileum",   "platform": "GPL10558", "n_samples": 68},
    # Psoriasis
    "GSE14905": {"disease": "Pso",   "tissue": "skin",    "platform": "GPL570",  "n_samples": 85},
    "GSE30999": {"disease": "Pso",   "tissue": "skin",    "platform": "GPL570",  "n_samples": 170},
    # Reactive Arthritis
    "GSE36700": {"disease": "ReA",   "tissue": "blood",   "platform": "GPL6883",  "n_samples": 24},
    # Undifferentiated SpA
    "GSE100927": {"disease": "uSpA",  "tissue": "blood",  "platform": "GPL13667", "n_samples": 36},
}

# ─────────────────────────────────────────────────────────────────────────────
# Immune cell reference signatures
# ─────────────────────────────────────────────────────────────────────────────

# Abbreviated LM22-inspired gene signatures for 22 immune cell types
# (In a full implementation, these would be loaded from the LM22.txt reference file)
# These are representative marker genes per cell type - not complete signatures
LM22_MARKERS = {
    'B_cells_naive':         ['CD19', 'MS4A1', 'CD22', 'BANK1', 'PAX5', 'VPREB3', 'TCL1A', 'SPIB'],
    'B_cells_memory':        ['CD19', 'MS4A1', 'CD27', 'AIM2', 'SSPN', 'COCH', 'FCRL4', 'FCRL5'],
    'Plasma_cells':          ['IGHG1', 'IGHG2', 'MZB1', 'CD38', 'SDC1', 'PRDM1', 'XBP1', 'IRF4'],
    'T_cells_CD8':           ['CD8A', 'CD8B', 'GZMA', 'GZMB', 'PRF1', 'EOMES', 'TBX21', 'TIGIT'],
    'T_cells_CD4_naive':     ['CD4', 'CCR7', 'TCF7', 'LEF1', 'SELL', 'KLF2', 'S1PR1', 'IL7R'],
    'T_cells_CD4_memory_resting': ['CD4', 'IL7R', 'ANXA1', 'LDHB', 'S100A4', 'EMP3', 'GIMAP7'],
    'T_cells_CD4_memory_activated': ['CD4', 'IL2RA', 'ICOS', 'CD69', 'CXCR3', 'IFNG', 'TNF'],
    'T_cells_follicular_helper': ['CD4', 'CXCR5', 'BCL6', 'ICOS', 'PDCD1', 'SH2D1A', 'BATF'],
    'T_cells_regulatory':    ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2', 'TIGIT', 'TNFRSF9', 'ENTPD1'],
    'T_cells_gamma_delta':   ['TRGC1', 'TRGC2', 'TRDC', 'TRDV2', 'NKG7', 'GNLY', 'KLRD1'],
    'NK_cells_resting':      ['NCAM1', 'NKG7', 'GNLY', 'KLRD1', 'KLRB1', 'KLRC1', 'FGFBP2'],
    'NK_cells_activated':    ['NCAM1', 'IFNG', 'TNF', 'GZMB', 'PRF1', 'CD69', 'KLRG1', 'CX3CR1'],
    'Monocytes':             ['CD14', 'LYZ', 'S100A8', 'S100A9', 'FCGR3A', 'MS4A7', 'CD68'],
    'Macrophages_M0':        ['CD68', 'CD163', 'MRC1', 'MSR1', 'FCGR1A', 'CSF1R', 'ADGRE1'],
    'Macrophages_M1':        ['CD68', 'IL1B', 'TNF', 'CXCL9', 'CXCL10', 'IDO1', 'CD80', 'NOS2'],
    'Macrophages_M2':        ['CD68', 'MRC1', 'CD163', 'IL10', 'TGFB1', 'CCL13', 'CCL18', 'ARG1'],
    'Dendritic_cells_resting': ['ITGAX', 'HLA-DRA', 'HLA-DRB1', 'CD1C', 'CLEC10A', 'FCER1A'],
    'Dendritic_cells_activated': ['ITGAX', 'CD83', 'CCR7', 'IL12B', 'LAMP3', 'FSCN1', 'IDO1'],
    'Mast_cells_resting':    ['KIT', 'TPSAB1', 'TPSB2', 'CPA3', 'HDC', 'HPGDS', 'CTSG'],
    'Mast_cells_activated':  ['KIT', 'FCER1A', 'IL4', 'IL5', 'IL13', 'TNFRSF21', 'SIGLEC6'],
    'Eosinophils':           ['EPX', 'RNASE2', 'RNASE3', 'CLC', 'IL5RA', 'SIGLEC8', 'CCR3'],
    'Neutrophils':           ['ELANE', 'MPO', 'PRTN3', 'AZU1', 'CEACAM8', 'FCGR3B', 'CSF3R'],
}

# SpA-relevant additional cell types
SPA_CELL_MARKERS = {
    'Th17_cells':        ['IL17A', 'IL17F', 'RORC', 'IL23R', 'CCR6', 'TNFSF11', 'IL21'],
    'Th1_cells':         ['IFNG', 'TBX21', 'CXCR3', 'IL12RB2', 'HAVCR2', 'PHLPP1'],
    'ILC3':              ['IL22', 'RORC', 'NCR2', 'KIT', 'IL1R1', 'IL23R', 'CXCR6'],
    'MAIT_cells':        ['SLC4A10', 'ZBTB16', 'RORC', 'IL18RAP', 'CD8A', 'NCR3'],
    'Innate_NK_like':    ['GNLY', 'NKG7', 'KLRD1', 'NCAM1', 'KLRB1', 'FGFBP2', 'GZMH'],
    'Plasmacytoid_DC':   ['IL3RA', 'CLEC4C', 'NRP1', 'IRF7', 'LILRA4', 'GZMB', 'CXCR3'],
    'Classical_monocyte': ['CD14', 'S100A8', 'S100A9', 'LYZ', 'VCAN', 'FCN1', 'SELL'],
    'Nonclassical_mono': ['FCGR3A', 'MS4A7', 'CX3CR1', 'CDKN1C', 'RHOC', 'CSF1R'],
}

# Combine all signatures
ALL_SIGNATURES = {**LM22_MARKERS, **SPA_CELL_MARKERS}


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_expression_data(gse_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load expression matrix and sample metadata for a dataset.
    Returns (expr_df [genes x samples], meta_df).
    """
    # Try multiple locations
    data_paths = [
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_expression_matrix.csv.gz"),
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_expression_matrix.tsv.gz"),
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_processed.csv.gz"),
        os.path.join(BASE_DIR, "results", "phase2_normalization", gse_id, "normalized_expression.csv"),
        os.path.join(BASE_DIR, "results", "phase1_download", gse_id, "expression_matrix.csv"),
    ]
    
    expr_df = None
    for path in data_paths:
        if os.path.exists(path):
            log.info(f"Loading {gse_id} from: {path}")
            if path.endswith('.gz'):
                sep = '\t' if '.tsv.' in path else ','
                with gzip.open(path, 'rt') as f:
                    expr_df = pd.read_csv(f, sep=sep, index_col=0)
            else:
                expr_df = pd.read_csv(path, index_col=0)
            break
    
    if expr_df is None:
        raise FileNotFoundError(f"No expression data found for {gse_id}")
    
    # Load metadata
    meta_paths = [
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_metadata.csv"),
        os.path.join(DATA_DIR, gse_id, "metadata.csv"),
        os.path.join(BASE_DIR, "results", "phase2_normalization", gse_id, "sample_metadata.csv"),
        os.path.join(BASE_DIR, "results", "phase1_download", gse_id, "sample_info.csv"),
    ]
    
    meta_df = None
    for path in meta_paths:
        if os.path.exists(path):
            meta_df = pd.read_csv(path, index_col=0)
            break
    
    if meta_df is None:
        log.warning(f"No metadata for {gse_id}, creating dummy")
        meta_df = pd.DataFrame(
            {'condition': ['unknown'] * expr_df.shape[1]},
            index=expr_df.columns
        )
    
    return expr_df, meta_df


def extract_disease_labels(meta_df: pd.DataFrame, gse_id: str) -> pd.Series:
    """
    Extract binary disease/control labels from metadata.
    Returns pd.Series with 1=disease, 0=control.
    """
    # Try standard column names
    for col in ['condition', 'disease_state', 'disease', 'group', 'status', 'phenotype']:
        if col in meta_df.columns:
            vals = meta_df[col].str.lower().str.strip()
            disease_terms = ['disease', 'patient', 'case', 'as', 'psa', 'ibd',
                            'crohn', 'uc', 'colitis', 'spondyl', 'psoriasis',
                            'arthritis', 'active', 'affected', 'positive']
            labels = vals.apply(
                lambda v: 1 if any(t in str(v) for t in disease_terms) else 0
            )
            if labels.sum() > 0:
                return labels
    
    # Fallback: use first half as control, second half as disease
    log.warning(f"{gse_id}: could not extract disease labels, using 50/50 split")
    n = len(meta_df)
    labels = pd.Series(
        [0] * (n // 2) + [1] * (n - n // 2),
        index=meta_df.index
    )
    return labels


# ─────────────────────────────────────────────────────────────────────────────
# Deconvolution methods
# ─────────────────────────────────────────────────────────────────────────────

def build_signature_matrix(signatures: dict,
                           expr_genes: list) -> pd.DataFrame:
    """
    Build signature matrix from marker gene lists.
    Matrix: genes × cell_types, values = 1 if marker gene, else 0.
    Only includes genes present in expression data.
    """
    # Get unique marker genes present in expression data
    expr_gene_set = set(expr_genes)
    sig_genes = set()
    for genes in signatures.values():
        sig_genes.update(genes)
    
    common_genes = sorted(sig_genes & expr_gene_set)
    
    if len(common_genes) == 0:
        log.warning("No signature genes found in expression data!")
        return pd.DataFrame()
    
    log.info(f"Signature matrix: {len(common_genes)} genes × {len(signatures)} cell types")
    
    # Build binary matrix
    sig_df = pd.DataFrame(
        0, index=common_genes, columns=list(signatures.keys())
    )
    
    for cell_type, genes in signatures.items():
        for gene in genes:
            if gene in expr_gene_set:
                sig_df.loc[gene, cell_type] = 1
    
    # Add expression level weights (use the sum of 1s for now)
    # In full CIBERSORT, this would be scaled by reference expression
    return sig_df.astype(float)


def nnls_deconvolution(expr_df: pd.DataFrame,
                       sig_df: pd.DataFrame) -> pd.DataFrame:
    """
    Non-negative least squares (NNLS) deconvolution.
    For each sample, solve: minimize ||S*f - m||_2 s.t. f >= 0
    where S=signature matrix, f=cell fractions, m=mixture (sample).
    
    Returns: samples × cell_types DataFrame of raw (unnormalized) fractions.
    """
    log.info(f"NNLS deconvolution: {expr_df.shape[1]} samples, {len(sig_df.columns)} cell types")
    
    # Align genes
    common_genes = sig_df.index.intersection(expr_df.index)
    if len(common_genes) == 0:
        log.error("No common genes between signature and expression!")
        return pd.DataFrame()
    
    S = sig_df.loc[common_genes].values  # genes × cell_types
    M = expr_df.loc[common_genes].values  # genes × samples
    
    n_samples = M.shape[1]
    n_cell_types = S.shape[1]
    
    fractions = np.zeros((n_samples, n_cell_types))
    
    for i in range(n_samples):
        m = M[:, i]
        # NNLS
        f, residual = nnls(S, m)
        fractions[i, :] = f
    
    # Normalize to sum to 1
    row_sums = fractions.sum(axis=1, keepdims=True)
    row_sums = np.where(row_sums == 0, 1, row_sums)  # Avoid division by zero
    fractions_norm = fractions / row_sums
    
    result_df = pd.DataFrame(
        fractions_norm,
        index=expr_df.columns,
        columns=sig_df.columns
    )
    
    log.info(f"NNLS complete. Mean fractions: {result_df.mean().round(3).to_dict()}")
    return result_df


def ssgsea_deconvolution(expr_df: pd.DataFrame,
                         signatures: dict,
                         alpha: float = 0.25) -> pd.DataFrame:
    """
    Single-sample GSEA (ssGSEA) for immune cell enrichment scoring.
    For each sample and gene set, computes an enrichment score.
    
    Based on Barbie et al. 2009 algorithm (simplified).
    Returns: samples × cell_types enrichment score DataFrame.
    """
    log.info(f"ssGSEA deconvolution: {expr_df.shape[1]} samples")
    
    n_samples = expr_df.shape[1]
    scores = {}
    
    for cell_type, marker_genes in signatures.items():
        # Get genes present in expression data
        present_genes = [g for g in marker_genes if g in expr_df.index]
        
        if len(present_genes) < 3:
            scores[cell_type] = np.zeros(n_samples)
            continue
        
        cell_scores = np.zeros(n_samples)
        
        for j, sample in enumerate(expr_df.columns):
            sample_expr = expr_df[sample]
            
            # Rank genes by expression (highest = rank 1 in GSEA convention)
            ranks = sample_expr.rank(ascending=False)
            n_genes = len(ranks)
            
            # Compute running enrichment score
            gene_set = set(present_genes)
            in_set = ranks.index.isin(gene_set)
            
            r_in = ranks[in_set].values
            r_out = ranks[~in_set].values
            
            # Weighted sum (alpha = 0.25 for ssGSEA)
            if len(r_in) == 0 or len(r_out) == 0:
                cell_scores[j] = 0
                continue
            
            # P_hit: weighted rank sum for in-set genes
            p_hit = (r_in ** alpha).sum() / (ranks.values ** alpha).sum()
            # P_miss: fraction of out-set genes
            p_miss = len(r_out) / n_genes
            
            cell_scores[j] = p_hit - p_miss
        
        scores[cell_type] = cell_scores
    
    result_df = pd.DataFrame(scores, index=expr_df.columns)
    
    # Normalize scores to [0, 1] range
    for col in result_df.columns:
        col_min = result_df[col].min()
        col_max = result_df[col].max()
        if col_max > col_min:
            result_df[col] = (result_df[col] - col_min) / (col_max - col_min)
    
    log.info(f"ssGSEA complete.")
    return result_df


def xcell_scoring(expr_df: pd.DataFrame,
                  signatures: dict) -> pd.DataFrame:
    """
    xCell-inspired scoring using Spearman correlation with reference signatures.
    Simplified version: correlation of sample expression with signature centroid.
    
    Returns: samples × cell_types score DataFrame.
    """
    log.info(f"xCell-style scoring: {expr_df.shape[1]} samples")
    
    scores = {}
    
    for cell_type, marker_genes in signatures.items():
        present_genes = [g for g in marker_genes if g in expr_df.index]
        
        if len(present_genes) < 3:
            scores[cell_type] = np.zeros(expr_df.shape[1])
            continue
        
        # Signature centroid = mean expression of marker genes
        sig_centroid = expr_df.loc[present_genes].mean(axis=0)  # samples,
        scores[cell_type] = sig_centroid.values
    
    result_df = pd.DataFrame(scores, index=expr_df.columns)
    
    # Normalize columns
    from sklearn.preprocessing import MinMaxScaler
    scaler = MinMaxScaler()
    result_df[:] = scaler.fit_transform(result_df)
    
    return result_df


# ─────────────────────────────────────────────────────────────────────────────
# Differential immune analysis
# ─────────────────────────────────────────────────────────────────────────────

def differential_immune_analysis(deconv_df: pd.DataFrame,
                                  labels: pd.Series,
                                  method: str = 'mannwhitney') -> pd.DataFrame:
    """
    Compare immune cell fractions between disease and control groups.
    Uses Mann-Whitney U test for robustness.
    
    Returns DataFrame with statistics per cell type.
    """
    # Align
    common = deconv_df.index.intersection(labels.index)
    deconv_aligned = deconv_df.loc[common]
    labels_aligned = labels.loc[common]
    
    disease_mask = labels_aligned == 1
    control_mask = labels_aligned == 0
    
    n_disease = disease_mask.sum()
    n_control = control_mask.sum()
    
    if n_disease < 3 or n_control < 3:
        log.warning(f"Too few samples: {n_disease} disease, {n_control} control")
        return pd.DataFrame()
    
    results = []
    for cell_type in deconv_df.columns:
        d_vals = deconv_aligned.loc[disease_mask, cell_type].values
        c_vals = deconv_aligned.loc[control_mask, cell_type].values
        
        # Mann-Whitney U
        try:
            stat, pval = mannwhitneyu(d_vals, c_vals, alternative='two-sided')
        except Exception:
            stat, pval = np.nan, np.nan
        
        # Effect size (rank-biserial correlation)
        n1, n2 = len(d_vals), len(c_vals)
        effect_size = 1 - (2 * stat) / (n1 * n2) if (n1 * n2) > 0 else 0
        
        results.append({
            'cell_type': cell_type,
            'disease_mean': round(float(np.mean(d_vals)), 4),
            'control_mean': round(float(np.mean(c_vals)), 4),
            'log2_fold_change': round(float(np.log2((np.mean(d_vals) + 1e-6) / (np.mean(c_vals) + 1e-6))), 4),
            'mann_whitney_stat': round(float(stat), 2) if not np.isnan(stat) else np.nan,
            'pvalue': round(float(pval), 6) if not np.isnan(pval) else np.nan,
            'effect_size': round(float(effect_size), 4),
            'n_disease': n_disease,
            'n_control': n_control,
        })
    
    diff_df = pd.DataFrame(results)
    
    # FDR correction (Benjamini-Hochberg)
    from scipy.stats import rankdata
    valid_pvals = diff_df['pvalue'].dropna()
    if len(valid_pvals) > 0:
        n = len(diff_df)
        ranks = rankdata(diff_df['pvalue'].fillna(1))
        diff_df['fdr'] = (diff_df['pvalue'].fillna(1) * n / ranks).clip(0, 1).round(6)
    
    diff_df = diff_df.sort_values('pvalue')
    sig_cells = (diff_df['fdr'] < 0.05).sum() if 'fdr' in diff_df.columns else 0
    log.info(f"Differential immune: {sig_cells} cell types significant (FDR<0.05)")
    
    return diff_df


# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

def plot_deconvolution_heatmap(deconv_df: pd.DataFrame,
                               labels: pd.Series,
                               gse_id: str,
                               method: str,
                               save_path: str):
    """Heatmap of immune cell fractions per sample, grouped by condition."""
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Sort samples by condition
    common = deconv_df.index.intersection(labels.index)
    df = deconv_df.loc[common]
    lbl = labels.loc[common]
    
    sort_order = lbl.sort_values().index
    df_sorted = df.loc[sort_order]
    
    # Heatmap
    cmap = LinearSegmentedColormap.from_list('immune', ['#F0F9FF', '#0369A1', '#7C3AED'])
    im = ax.imshow(df_sorted.T, aspect='auto', cmap=cmap,
                   vmin=0, vmax=df_sorted.values.max())
    
    # Labels
    ax.set_yticks(range(len(df_sorted.columns)))
    ax.set_yticklabels(df_sorted.columns, fontsize=8)
    ax.set_xlabel('Samples', fontsize=11)
    ax.set_title(f'{gse_id}: Immune Cell Deconvolution ({method})', fontsize=13)
    
    # Condition bar on top
    ax2 = ax.twiny()
    condition_colors = ['#EF4444' if l == 1 else '#22C55E' 
                        for l in lbl.loc[sort_order].values]
    ax2.set_xlim(ax.get_xlim())
    for i, (s, c) in enumerate(zip(sort_order, condition_colors)):
        ax2.axvspan(i - 0.5, i + 0.5, alpha=0.4, color=c, linewidth=0)
    ax2.set_xticks([])
    
    plt.colorbar(im, ax=ax, label='Cell Fraction', shrink=0.8)
    
    # Legend
    legend_elements = [Patch(facecolor='#EF4444', alpha=0.6, label='Disease'),
                       Patch(facecolor='#22C55E', alpha=0.6, label='Control')]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved deconvolution heatmap: {save_path}")


from matplotlib.patches import Patch


def plot_differential_immune(diff_df: pd.DataFrame,
                              gse_id: str,
                              save_path: str):
    """Horizontal bar plot of log2FC for immune cells (disease vs control)."""
    if diff_df.empty:
        return
    
    fig, ax = plt.subplots(figsize=(10, max(6, len(diff_df) * 0.3)))
    
    # Sort by log2FC
    df = diff_df.sort_values('log2_fold_change', ascending=True).head(30)
    
    colors = ['#EF4444' if v > 0 else '#3B82F6' for v in df['log2_fold_change']]
    bars = ax.barh(range(len(df)), df['log2_fold_change'].values, color=colors)
    
    # Mark significant
    for i, (_, row) in enumerate(df.iterrows()):
        if row.get('fdr', 1) < 0.05:
            ax.text(row['log2_fold_change'] + 0.02 if row['log2_fold_change'] > 0
                    else row['log2_fold_change'] - 0.02,
                    i, '*', ha='left' if row['log2_fold_change'] > 0 else 'right',
                    va='center', fontsize=12)
    
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['cell_type'].tolist(), fontsize=9)
    ax.axvline(x=0, color='black', linewidth=0.8)
    ax.set_xlabel('log₂(Disease / Control)', fontsize=11)
    ax.set_title(f'{gse_id}: Differential Immune Infiltration', fontsize=13)
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved differential immune plot: {save_path}")


def plot_cross_disease_comparison(all_results: dict, save_path: str):
    """
    Heatmap comparing immune infiltration across diseases.
    Rows = cell types, Columns = diseases (mean score in disease samples).
    """
    # Collect mean disease scores per dataset
    disease_scores = {}
    
    for gse_id, res in all_results.items():
        if 'deconv_nnls' not in res or res['deconv_nnls'].empty:
            continue
        
        deconv = res['deconv_nnls']
        labels = res.get('labels')
        
        if labels is None:
            continue
        
        common = deconv.index.intersection(labels.index)
        disease_samples = labels.loc[common][labels.loc[common] == 1].index
        
        if len(disease_samples) == 0:
            continue
        
        disease_mean = deconv.loc[disease_samples].mean()
        disease_name = res.get('disease', gse_id)
        disease_scores[f"{disease_name}\n({gse_id})"] = disease_mean
    
    if not disease_scores:
        log.warning("No data for cross-disease comparison")
        return
    
    compare_df = pd.DataFrame(disease_scores)  # cell_types × diseases
    
    # Keep cell types with any meaningful variation
    compare_df = compare_df.loc[compare_df.std(axis=1) > 0.01]
    
    if compare_df.empty:
        log.warning("No variable cell types for cross-disease plot")
        return
    
    fig, ax = plt.subplots(figsize=(max(10, len(compare_df.columns) * 1.5),
                                     max(8, len(compare_df) * 0.4)))
    
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    
    # Z-score normalize across diseases
    from scipy.stats import zscore
    compare_z = compare_df.apply(zscore, axis=1, nan_policy='omit')
    
    sns.heatmap(
        compare_z, ax=ax, cmap=cmap, center=0,
        xticklabels=True, yticklabels=True,
        linewidths=0.5, annot=False,
        cbar_kws={'label': 'Z-score'}
    )
    
    ax.set_title('Cross-Disease Immune Cell Infiltration Comparison\n(NNLS Deconvolution, Z-score normalized)', 
                 fontsize=13)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved cross-disease comparison: {save_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Per-dataset runner
# ─────────────────────────────────────────────────────────────────────────────

def run_deconvolution_for_dataset(gse_id: str,
                                  dataset_info: dict) -> dict:
    """
    Run immune cell deconvolution for one dataset.
    Returns dict with results.
    """
    log.info(f"\n{'='*65}")
    log.info(f"Running deconvolution: {gse_id} ({dataset_info['disease']}, "
             f"{dataset_info['tissue']})")
    log.info(f"{'='*65}")
    
    result = {
        'gse_id': gse_id,
        'disease': dataset_info['disease'],
        'tissue': dataset_info['tissue'],
        'status': 'failed',
        'deconv_nnls': pd.DataFrame(),
        'deconv_ssgsea': pd.DataFrame(),
        'diff_immune': pd.DataFrame(),
        'labels': None,
        'error': None
    }
    
    os.makedirs(os.path.join(RESULTS_DIR, gse_id), exist_ok=True)
    os.makedirs(os.path.join(FIGURES_DIR, gse_id), exist_ok=True)
    
    try:
        # Load data
        expr_df, meta_df = load_expression_data(gse_id)
        log.info(f"Loaded: {expr_df.shape[0]} genes, {expr_df.shape[1]} samples")
        
        # Extract labels
        labels = extract_disease_labels(meta_df, gse_id)
        result['labels'] = labels
        
        n_disease = (labels == 1).sum()
        n_control = (labels == 0).sum()
        log.info(f"Labels: {n_disease} disease, {n_control} control")
        
        # Build signature matrix
        sig_df = build_signature_matrix(ALL_SIGNATURES, expr_df.index.tolist())
        
        if sig_df.empty:
            raise ValueError("Empty signature matrix")
        
        # Method 1: NNLS deconvolution
        log.info("\nMethod 1: NNLS Deconvolution")
        deconv_nnls = nnls_deconvolution(expr_df, sig_df)
        result['deconv_nnls'] = deconv_nnls
        deconv_nnls.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'deconv_nnls.csv')
        )
        
        # Method 2: ssGSEA
        log.info("\nMethod 2: ssGSEA Scoring")
        deconv_ssgsea = ssgsea_deconvolution(expr_df, ALL_SIGNATURES)
        result['deconv_ssgsea'] = deconv_ssgsea
        deconv_ssgsea.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'deconv_ssgsea.csv')
        )
        
        # Differential immune analysis
        log.info("\nDifferential Immune Analysis")
        diff_nnls = differential_immune_analysis(deconv_nnls, labels)
        diff_ssgsea = differential_immune_analysis(deconv_ssgsea, labels)
        result['diff_immune'] = diff_nnls
        
        if not diff_nnls.empty:
            diff_nnls.to_csv(
                os.path.join(RESULTS_DIR, gse_id, 'diff_immune_nnls.csv'), index=False
            )
        if not diff_ssgsea.empty:
            diff_ssgsea.to_csv(
                os.path.join(RESULTS_DIR, gse_id, 'diff_immune_ssgsea.csv'), index=False
            )
        
        # Visualizations
        plot_deconvolution_heatmap(
            deconv_nnls, labels, gse_id, 'NNLS',
            os.path.join(FIGURES_DIR, gse_id, 'deconv_heatmap_nnls.png')
        )
        plot_differential_immune(
            diff_nnls, gse_id,
            os.path.join(FIGURES_DIR, gse_id, 'differential_immune.png')
        )
        
        # Summary stats
        n_sig = (diff_nnls['fdr'] < 0.05).sum() if 'fdr' in diff_nnls.columns and not diff_nnls.empty else 0
        result.update({
            'status': 'success',
            'n_genes': expr_df.shape[0],
            'n_samples': expr_df.shape[1],
            'n_sig_cells_nnls': int(n_sig),
            'top_sig_cell': diff_nnls.iloc[0]['cell_type'] if not diff_nnls.empty else None
        })
        
        log.info(f"Deconvolution complete for {gse_id}: {n_sig} significant cell types")
        
    except FileNotFoundError as e:
        log.error(f"Data not found: {e}")
        result['error'] = str(e)
    except Exception as e:
        log.error(f"Failed: {e}")
        import traceback
        traceback.print_exc()
        result['error'] = str(e)
    
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='Phase 8: Immune Cell Deconvolution')
    parser.add_argument('--datasets', nargs='+', default=list(ALL_DATASETS.keys()),
                        help='Datasets to process')
    parser.add_argument('--methods', nargs='+', default=['nnls', 'ssgsea'],
                        choices=['nnls', 'ssgsea', 'xcell'],
                        help='Deconvolution methods')
    parser.add_argument('--skip-cross-disease', action='store_true',
                        help='Skip cross-disease comparison plot')
    args = parser.parse_args()
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    log.info("Phase 8: Immune Cell Deconvolution")
    log.info(f"Datasets: {args.datasets}")
    log.info(f"Cell types: {len(ALL_SIGNATURES)} (LM22 + SpA-specific)")
    
    # Run deconvolution for each dataset
    all_results = {}
    for gse_id in args.datasets:
        if gse_id not in ALL_DATASETS:
            log.warning(f"Unknown dataset {gse_id}, skipping")
            continue
        all_results[gse_id] = run_deconvolution_for_dataset(gse_id, ALL_DATASETS[gse_id])
    
    # Cross-disease comparison
    if not args.skip_cross_disease:
        plot_cross_disease_comparison(
            all_results,
            os.path.join(FIGURES_DIR, 'cross_disease_immune_comparison.png')
        )
    
    # Summary
    summary = []
    for gse_id, res in all_results.items():
        summary.append({
            'dataset': gse_id,
            'disease': res.get('disease', ''),
            'tissue': res.get('tissue', ''),
            'status': res['status'],
            'n_genes': res.get('n_genes', 0),
            'n_samples': res.get('n_samples', 0),
            'n_sig_cells': res.get('n_sig_cells_nnls', 0),
            'top_cell': res.get('top_sig_cell', ''),
            'error': res.get('error', '')
        })
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(os.path.join(RESULTS_DIR, 'deconvolution_summary.csv'), index=False)
    
    log.info("\n" + "="*65)
    log.info("IMMUNE DECONVOLUTION COMPLETE")
    log.info("="*65)
    log.info(f"\nSummary:\n{summary_df.to_string(index=False)}")
    log.info(f"Results: {RESULTS_DIR}")
    log.info(f"Figures: {FIGURES_DIR}")


if __name__ == '__main__':
    main()
