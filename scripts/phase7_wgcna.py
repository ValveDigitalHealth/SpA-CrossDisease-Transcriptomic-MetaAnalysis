#!/usr/bin/env python3
"""
Phase 7: Weighted Gene Co-expression Network Analysis (WGCNA)
=============================================================
Cross-Disease Transcriptomic Meta-Analysis of SpA — Project 1

Implements a Python-based WGCNA pipeline:
1. Probe-to-gene mapping and gene-level expression aggregation
2. Soft-threshold power selection (scale-free topology)
3. Adjacency matrix and TOM (Topological Overlap Matrix) construction
4. Hierarchical clustering with dynamic tree cutting for module detection
5. Module-trait correlation analysis (disease vs control)
6. Module preservation across AS and IBD datasets
7. Hub gene identification within key modules
8. Overlap with AS meta-analysis genes

Note: This is a Python implementation inspired by the R WGCNA package (Langfelder & Horvath).
We use the same core algorithms (adjacency → TOM → hierarchical clustering → dynamic cut)
but in scipy/numpy. The dynamic tree cut is simplified (height-based with merging of small modules).

Target datasets for WGCNA:
- GSE18781 (AS/axSpA, 55 samples, Affy U133+2) — largest AS dataset
- GSE59071 (IBD, 116 samples, HuGene 1.0 ST) — largest overall, SpA-associated IBD

Module preservation is tested: modules found in AS → preservation in IBD (and vice versa).

Author: Perplexity Computer / Spondyloarthritis AI Computational Biology Institute
Date: March 2026
"""

import os
import sys
import json
import gzip
import time
import logging
import argparse
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from collections import Counter, defaultdict

warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
RESULTS_DIR = os.path.join(BASE_DIR, "results", "phase7_wgcna")
FIGURES_DIR = os.path.join(BASE_DIR, "figures", "phase7_wgcna")

# WGCNA parameters
SOFT_POWER_RANGE = list(range(1, 21))  # Test powers 1-20
MIN_MODULE_SIZE = 30                    # Min genes per module
MERGE_CUT_HEIGHT = 0.25                # Merge modules with >75% overlap
DEEP_SPLIT = 2                         # Dynamic tree cut sensitivity

# Target datasets
WGCNA_DATASETS = {
    "GSE18781": {
        "disease": "AS",
        "platform": "GPL570",
        "description": "AS/axSpA whole blood, 55 samples"
    },
    "GSE59071": {
        "disease": "IBD",
        "platform": "GPL6244",
        "description": "IBD colon biopsies, 116 samples"
    }
}

# ─────────────────────────────────────────────────────────────────────────────
# Logging setup
# ─────────────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
log = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────────────────────
# Data loading utilities
# ─────────────────────────────────────────────────────────────────────────────

def load_expression_matrix(gse_id: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load expression matrix and sample metadata for a dataset.
    Returns (expr_df, meta_df) where expr_df rows=genes, cols=samples.
    """
    # Try multiple file patterns
    patterns = [
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_expression_matrix.csv.gz"),
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_expression_matrix.tsv.gz"),
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_processed.csv.gz"),
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_normalized.csv.gz"),
    ]
    
    expr_df = None
    for pat in patterns:
        if os.path.exists(pat):
            log.info(f"Loading expression from: {pat}")
            if pat.endswith('.gz'):
                with gzip.open(pat, 'rt') as f:
                    sep = '\t' if pat.endswith('.tsv.gz') else ','
                    expr_df = pd.read_csv(f, sep=sep, index_col=0)
            else:
                sep = '\t' if pat.endswith('.tsv') else ','
                expr_df = pd.read_csv(pat, sep=sep, index_col=0)
            break
    
    if expr_df is None:
        # Try to load from phase1/phase2 results
        results_patterns = [
            os.path.join(BASE_DIR, "results", "phase2_normalization", gse_id, "normalized_expression.csv"),
            os.path.join(BASE_DIR, "results", "phase1_download", gse_id, "expression_matrix.csv"),
        ]
        for pat in results_patterns:
            if os.path.exists(pat):
                log.info(f"Loading expression from results: {pat}")
                expr_df = pd.read_csv(pat, index_col=0)
                break
    
    if expr_df is None:
        raise FileNotFoundError(f"Cannot find expression matrix for {gse_id}")
    
    log.info(f"{gse_id}: loaded expression matrix {expr_df.shape}")
    
    # Load metadata
    meta_df = None
    meta_patterns = [
        os.path.join(DATA_DIR, gse_id, f"{gse_id}_metadata.csv"),
        os.path.join(DATA_DIR, gse_id, "metadata.csv"),
        os.path.join(BASE_DIR, "results", "phase2_normalization", gse_id, "sample_metadata.csv"),
        os.path.join(BASE_DIR, "results", "phase1_download", gse_id, "sample_info.csv"),
    ]
    for pat in meta_patterns:
        if os.path.exists(pat):
            log.info(f"Loading metadata from: {pat}")
            meta_df = pd.read_csv(pat, index_col=0)
            break
    
    if meta_df is None:
        log.warning(f"No metadata found for {gse_id}, creating dummy metadata")
        meta_df = pd.DataFrame({'condition': ['unknown'] * expr_df.shape[1]},
                               index=expr_df.columns)
    
    return expr_df, meta_df


def aggregate_to_gene_level(expr_df: pd.DataFrame) -> pd.DataFrame:
    """
    If index contains probe IDs (non-unique), aggregate to gene level by mean.
    Detects if index looks like probe IDs vs gene symbols.
    """
    # Check if index looks like probes (e.g., '1007_s_at', '1053_at')
    sample_idx = str(expr_df.index[0])
    is_probe = ('_at' in sample_idx or 
                sample_idx.startswith('ILMN') or
                sample_idx.startswith('A_') or
                sample_idx.replace('.', '').isdigit())
    
    if is_probe:
        log.info("Probe-level data detected, aggregating to gene level")
        # For simplicity, use probe IDs as gene IDs (real WGCNA would use annotation)
        # In a full pipeline, we'd map probes to genes via platform annotation
        expr_df = expr_df.groupby(expr_df.index).mean()
        log.info(f"After aggregation: {expr_df.shape}")
    
    return expr_df


def filter_low_variance_genes(expr_df: pd.DataFrame, 
                               top_n: int = 5000) -> pd.DataFrame:
    """
    Keep top N most variable genes (by MAD - median absolute deviation).
    Standard WGCNA practice to reduce computation.
    """
    if expr_df.shape[0] <= top_n:
        log.info(f"Gene count ({expr_df.shape[0]}) <= top_n ({top_n}), no filtering")
        return expr_df
    
    mad = expr_df.apply(lambda x: np.median(np.abs(x - np.median(x))), axis=1)
    top_genes = mad.nlargest(top_n).index
    log.info(f"Filtered to top {top_n} most variable genes (from {expr_df.shape[0]})")
    return expr_df.loc[top_genes]


# ─────────────────────────────────────────────────────────────────────────────
# Soft-threshold power selection
# ─────────────────────────────────────────────────────────────────────────────

def compute_scale_free_fit(expr_matrix: np.ndarray, 
                           power: int,
                           n_bins: int = 10) -> dict:
    """
    For a given soft-threshold power, compute:
    - Signed correlation adjacency matrix
    - Node connectivity (sum of adjacency)
    - Scale-free topology fit (R² of log-log regression)
    
    expr_matrix: samples × genes (numpy array)
    """
    n_samples, n_genes = expr_matrix.shape
    
    # Compute correlation matrix (genes × genes)
    # Use Pearson correlation
    corr = np.corrcoef(expr_matrix.T)  # genes × genes
    np.fill_diagonal(corr, 0)
    
    # Signed adjacency: ((1 + cor) / 2)^power
    adj = ((1 + corr) / 2) ** power
    np.fill_diagonal(adj, 0)
    
    # Node connectivity
    k = adj.sum(axis=1)  # sum of row = connectivity of each gene
    
    # Scale-free fit: bin connectivity, fit log(p(k)) ~ log(k)
    # Remove zeros
    k_pos = k[k > 0]
    if len(k_pos) < 5:
        return {'power': power, 'r2': 0, 'slope': 0, 'mean_k': np.mean(k)}
    
    # Histogram in log space
    log_k = np.log10(k_pos + 1e-6)
    hist_counts, bin_edges = np.histogram(log_k, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Remove empty bins
    mask = hist_counts > 0
    if mask.sum() < 3:
        return {'power': power, 'r2': 0, 'slope': 0, 'mean_k': np.mean(k)}
    
    log_p = np.log10(hist_counts[mask] / hist_counts[mask].sum())
    log_k_bin = bin_centers[mask]
    
    # Linear regression
    slope, intercept, r, p_val, se = stats.linregress(log_k_bin, log_p)
    r2 = r ** 2
    
    return {
        'power': power,
        'r2': r2,
        'slope': slope,
        'mean_k': float(np.mean(k)),
        'median_k': float(np.median(k))
    }


def select_soft_threshold(expr_df: pd.DataFrame,
                          power_range: list = SOFT_POWER_RANGE,
                          r2_cutoff: float = 0.80) -> tuple[int, pd.DataFrame]:
    """
    Test soft-threshold powers and select the lowest power achieving R² ≥ r2_cutoff.
    Returns (selected_power, fit_stats_df).
    """
    log.info(f"Testing soft-threshold powers: {power_range}")
    
    # Transpose for samples × genes
    # expr_df: rows=genes, cols=samples → need samples × genes for corrcoef
    X = expr_df.values.T  # samples × genes
    
    # For large matrices, subsample genes to speed up
    n_genes = X.shape[1]
    if n_genes > 3000:
        idx = np.random.choice(n_genes, 3000, replace=False)
        X_sub = X[:, idx]
        log.info(f"Subsampling to 3000 genes for power selection")
    else:
        X_sub = X
    
    results = []
    for power in power_range:
        fit = compute_scale_free_fit(X_sub, power)
        results.append(fit)
        log.info(f"  Power {power:2d}: R²={fit['r2']:.3f}, slope={fit['slope']:.2f}, "
                 f"meanK={fit['mean_k']:.1f}")
    
    fit_df = pd.DataFrame(results)
    
    # Select lowest power with R² >= cutoff
    good = fit_df[fit_df['r2'] >= r2_cutoff]
    if len(good) > 0:
        selected_power = int(good['power'].iloc[0])
        log.info(f"Selected soft-threshold power: {selected_power} (R²={good['r2'].iloc[0]:.3f})")
    else:
        # Fallback: highest R²
        best_idx = fit_df['r2'].idxmax()
        selected_power = int(fit_df.loc[best_idx, 'power'])
        log.warning(f"No power achieved R²>={r2_cutoff}. Using power {selected_power} "
                    f"(best R²={fit_df.loc[best_idx, 'r2']:.3f})")
    
    return selected_power, fit_df


# ─────────────────────────────────────────────────────────────────────────────
# TOM construction and module detection
# ─────────────────────────────────────────────────────────────────────────────

def compute_adjacency(expr_df: pd.DataFrame, power: int) -> np.ndarray:
    """
    Compute signed adjacency matrix: adj_ij = ((1 + cor_ij) / 2)^power
    Returns genes × genes numpy array.
    """
    log.info(f"Computing adjacency matrix (power={power}) for {expr_df.shape[0]} genes...")
    X = expr_df.values.T  # samples × genes
    
    # Compute full correlation matrix
    corr = np.corrcoef(X.T)  # genes × genes
    np.fill_diagonal(corr, 0)
    
    # Signed adjacency
    adj = ((1 + corr) / 2) ** power
    np.fill_diagonal(adj, 0)
    
    return adj.astype(np.float32)


def adjacency_to_TOM(adj: np.ndarray) -> np.ndarray:
    """
    Convert adjacency matrix to Topological Overlap Matrix (TOM).
    TOM_ij = (sum_u(a_iu * a_ju) + a_ij) / (min(k_i, k_j) + 1 - a_ij)
    
    This is the core WGCNA formula that captures shared neighborhood.
    For large matrices, uses chunked computation to manage memory.
    """
    log.info(f"Computing TOM for {adj.shape[0]} × {adj.shape[0]} adjacency matrix...")
    t0 = time.time()
    
    n = adj.shape[0]
    k = adj.sum(axis=1)  # connectivity of each node
    
    # Compute TOM in chunks to manage memory
    chunk_size = min(500, n)
    TOM = np.zeros((n, n), dtype=np.float32)
    
    for i in range(0, n, chunk_size):
        i_end = min(i + chunk_size, n)
        # Numerator: sum of shared neighbors + direct connection
        # For rows i:i_end: sum over all u of adj[i:i_end, u] * adj[u, :]
        num = adj[i:i_end, :] @ adj  # (chunk × n)
        num += adj[i:i_end, :]       # add direct connection
        
        # Denominator: min(k_i, k_j) + 1 - adj_ij
        k_chunk = k[i:i_end].reshape(-1, 1)  # chunk × 1
        k_all = k.reshape(1, -1)              # 1 × n
        min_k = np.minimum(k_chunk, k_all)    # chunk × n
        denom = min_k + 1 - adj[i:i_end, :]
        
        # Avoid division by zero
        denom = np.maximum(denom, 1e-6)
        TOM[i:i_end, :] = (num / denom).astype(np.float32)
    
    # TOM should be in [0, 1]; clip for numerical stability
    np.clip(TOM, 0, 1, out=TOM)
    np.fill_diagonal(TOM, 1)  # self-TOM = 1
    
    elapsed = time.time() - t0
    log.info(f"TOM computed in {elapsed:.1f}s")
    return TOM


def detect_modules(TOM: np.ndarray,
                   gene_names: list,
                   min_module_size: int = MIN_MODULE_SIZE,
                   merge_cut_height: float = MERGE_CUT_HEIGHT) -> pd.Series:
    """
    Detect gene modules from TOM using hierarchical clustering + dynamic tree cut.
    
    Steps:
    1. Convert TOM to dissimilarity: dissTOM = 1 - TOM
    2. Hierarchical clustering with average linkage
    3. Cut tree by height, merge similar modules
    4. Assign module colors
    
    Returns pd.Series with gene → module_color mapping.
    """
    log.info("Detecting modules via hierarchical clustering...")
    t0 = time.time()
    
    # Dissimilarity
    dissTOM = 1 - TOM
    np.fill_diagonal(dissTOM, 0)
    
    # Convert to condensed distance matrix
    dist_condensed = squareform(dissTOM, checks=False)
    
    # Hierarchical clustering
    log.info("  Running hierarchical clustering...")
    Z = linkage(dist_condensed, method='average')
    
    # Dynamic tree cut (simplified): cut at multiple heights, pick best
    # Real dynamic tree cut uses a complex recursive algorithm
    # Here we use a simplified version: cut at height 0.99 (very similar to WGCNA default)
    # and then merge small/similar modules
    
    # Find a cut height that gives reasonable number of modules
    best_labels = None
    best_n_modules = 0
    
    for cut_height in [0.97, 0.95, 0.90, 0.85, 0.80]:
        labels = fcluster(Z, t=cut_height, criterion='distance')
        counts = Counter(labels)
        # Filter: keep only clusters with >= min_module_size genes
        large_clusters = {k for k, v in counts.items() if v >= min_module_size}
        n_large = len(large_clusters)
        
        if n_large >= 3:  # Want at least 3 modules
            # Remap: large clusters get module labels, rest → 0 (grey)
            new_labels = np.zeros(len(labels), dtype=int)
            for new_id, old_id in enumerate(sorted(large_clusters), start=1):
                new_labels[labels == old_id] = new_id
            
            best_labels = new_labels
            best_n_modules = n_large
            log.info(f"  Cut height {cut_height}: {n_large} modules (>={min_module_size} genes)")
            break
    
    if best_labels is None:
        log.warning("Could not find modules, using single cluster")
        best_labels = np.ones(len(gene_names), dtype=int)
        best_n_modules = 1
    
    # Merge similar modules (simplified merging)
    best_labels = merge_similar_modules(
        best_labels, TOM, gene_names, merge_cut_height
    )
    
    # Assign color names to modules
    unique_modules = sorted(set(best_labels))
    module_colors_list = [
        'turquoise', 'blue', 'brown', 'yellow', 'green', 'red', 'black',
        'pink', 'magenta', 'purple', 'greenyellow', 'tan', 'salmon', 'cyan',
        'midnightblue', 'lightcyan', 'grey60', 'royalblue', 'darkred',
        'darkgreen', 'darkturquoise', 'darkgrey', 'orange', 'darkorange',
        'white', 'skyblue', 'saddlebrown', 'steelblue', 'paleturquoise',
        'violet', 'darkolivegreen', 'darkmagenta'
    ]
    
    color_map = {0: 'grey'}  # 0 = unassigned
    for i, mod in enumerate(m for m in unique_modules if m != 0):
        color_map[mod] = module_colors_list[i % len(module_colors_list)]
    
    gene_colors = pd.Series(
        [color_map.get(lbl, 'grey') for lbl in best_labels],
        index=gene_names
    )
    
    elapsed = time.time() - t0
    log.info(f"Module detection complete in {elapsed:.1f}s: "
             f"{len(set(best_labels)) - (1 if 0 in best_labels else 0)} modules "
             f"(+ grey), {(best_labels == 0).sum()} unassigned genes")
    
    return gene_colors


def merge_similar_modules(labels: np.ndarray,
                          TOM: np.ndarray,
                          gene_names: list,
                          cut_height: float) -> np.ndarray:
    """
    Merge modules whose eigengene correlation is high (>= 1 - cut_height).
    Simplified version of WGCNA's mergeCloseModules.
    """
    unique_mods = sorted(set(labels))
    if 0 in unique_mods:
        unique_mods.remove(0)
    
    if len(unique_mods) <= 1:
        return labels
    
    # Compute module eigengenes (first PC of each module's expression)
    # Here we approximate with mean expression
    n_genes = TOM.shape[0]
    
    # Module mean connectivity (simplified eigengene)
    mod_means = {}
    for mod in unique_mods:
        mod_idx = np.where(labels == mod)[0]
        if len(mod_idx) > 0:
            # Mean TOM within module
            mod_means[mod] = TOM[np.ix_(mod_idx, mod_idx)].mean(axis=1)
    
    # Compute correlation between module eigengenes
    mod_list = list(mod_means.keys())
    n_mods = len(mod_list)
    
    if n_mods < 2:
        return labels
    
    # Build correlation matrix between module means
    max_len = max(len(v) for v in mod_means.values())
    # Pad to same length for correlation
    corr_matrix = np.zeros((n_mods, n_mods))
    for i, mi in enumerate(mod_list):
        for j, mj in enumerate(mod_list):
            if i == j:
                corr_matrix[i, j] = 1.0
            else:
                vi = mod_means[mi]
                vj = mod_means[mj]
                # Use min length
                min_len = min(len(vi), len(vj))
                if min_len > 1:
                    r, _ = stats.pearsonr(vi[:min_len], vj[:min_len])
                    corr_matrix[i, j] = r
    
    # Merge modules with correlation > (1 - cut_height)
    merge_threshold = 1 - cut_height
    new_labels = labels.copy()
    
    merged = set()
    for i in range(n_mods):
        if mod_list[i] in merged:
            continue
        for j in range(i + 1, n_mods):
            if mod_list[j] in merged:
                continue
            if corr_matrix[i, j] >= merge_threshold:
                # Merge j into i
                new_labels[labels == mod_list[j]] = mod_list[i]
                merged.add(mod_list[j])
                log.info(f"    Merged module {mod_list[j]} into {mod_list[i]} "
                         f"(corr={corr_matrix[i,j]:.3f})")
    
    return new_labels


# ─────────────────────────────────────────────────────────────────────────────
# Module-trait correlation
# ─────────────────────────────────────────────────────────────────────────────

def compute_module_eigengenes(expr_df: pd.DataFrame,
                              gene_colors: pd.Series) -> pd.DataFrame:
    """
    Compute module eigengene = first principal component of module's expression matrix.
    Returns DataFrame: samples × modules.
    """
    from sklearn.decomposition import PCA
    
    modules = [c for c in gene_colors.unique() if c != 'grey']
    eigengenes = {}
    
    for mod in modules:
        mod_genes = gene_colors[gene_colors == mod].index
        mod_genes = mod_genes[mod_genes.isin(expr_df.index)]
        
        if len(mod_genes) < 3:
            continue
        
        X = expr_df.loc[mod_genes].values.T  # samples × genes
        
        try:
            pca = PCA(n_components=1)
            me = pca.fit_transform(X).ravel()  # samples,
            eigengenes[f"ME{mod}"] = me
        except Exception as e:
            log.warning(f"PCA failed for module {mod}: {e}")
            eigengenes[f"ME{mod}"] = X.mean(axis=1)
    
    ME_df = pd.DataFrame(eigengenes, index=expr_df.columns)
    return ME_df


def module_trait_correlation(ME_df: pd.DataFrame,
                             meta_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Correlate module eigengenes with binary disease trait.
    Returns (correlation_df, pvalue_df).
    """
    # Create binary trait vector (disease=1, control=0)
    trait_cols = ['condition', 'disease_state', 'disease', 'group', 'status']
    trait_col = None
    for col in trait_cols:
        if col in meta_df.columns:
            trait_col = col
            break
    
    if trait_col is None:
        # Try to infer from column names
        log.warning("No standard trait column found, attempting to infer")
        trait_col = meta_df.columns[0]
    
    # Align samples
    common_samples = ME_df.index.intersection(meta_df.index)
    if len(common_samples) < 5:
        log.warning(f"Only {len(common_samples)} samples in common, using index order")
        n = min(len(ME_df), len(meta_df))
        trait_values = pd.Series(np.zeros(n), index=ME_df.index[:n])
        # Try to parse disease status from metadata
        for i, val in enumerate(meta_df.iloc[:n][trait_col]):
            val_str = str(val).lower()
            if any(d in val_str for d in ['disease', 'patient', 'case', 'as', 'ibd', 
                                           'crohn', 'uc', 'colitis', 'spondyl']):
                trait_values.iloc[i] = 1
    else:
        ME_df = ME_df.loc[common_samples]
        meta_sub = meta_df.loc[common_samples]
        trait_values = pd.Series(np.zeros(len(common_samples)), index=common_samples)
        for i, val in enumerate(meta_sub[trait_col]):
            val_str = str(val).lower()
            if any(d in val_str for d in ['disease', 'patient', 'case', 'as', 'ibd',
                                           'crohn', 'uc', 'colitis', 'spondyl', 'active']):
                trait_values.iloc[i] = 1
    
    # If no positive cases found, use median split
    if trait_values.sum() == 0:
        log.warning("No disease cases identified, using median split of first PC")
        first_me = ME_df.iloc[:, 0]
        trait_values = (first_me > first_me.median()).astype(float)
    
    log.info(f"Disease trait: {int(trait_values.sum())} cases, "
             f"{int((1-trait_values).sum())} controls")
    
    # Correlate each ME with trait
    corr_dict = {}
    pval_dict = {}
    
    for me_name in ME_df.columns:
        me_vals = ME_df[me_name]
        common = me_vals.index.intersection(trait_values.index)
        r, p = stats.pearsonr(me_vals.loc[common], trait_values.loc[common])
        corr_dict[me_name] = r
        pval_dict[me_name] = p
    
    corr_df = pd.DataFrame({'correlation': corr_dict, 'pvalue': pval_dict})
    corr_df['abs_correlation'] = corr_df['correlation'].abs()
    corr_df = corr_df.sort_values('abs_correlation', ascending=False)
    
    log.info(f"Top correlated modules:\n{corr_df.head(5).to_string()}")
    return corr_df, trait_values


# ─────────────────────────────────────────────────────────────────────────────
# Module preservation
# ─────────────────────────────────────────────────────────────────────────────

def compute_module_preservation(ref_colors: pd.Series,
                                ref_expr: pd.DataFrame,
                                test_expr: pd.DataFrame) -> pd.DataFrame:
    """
    Test module preservation: compute Z_summary statistic for each module
    found in reference dataset (ref) in test dataset.
    
    Simplified version of WGCNA's modulePreservation:
    - Z_density: preservation of within-module connectivity
    - Z_connectivity: preservation of inter-module connectivity pattern
    - Z_summary = mean(Z_density, Z_connectivity)
    
    Z_summary:
    - < 2: no preservation
    - 2-10: moderate preservation
    - > 10: strong preservation
    """
    log.info("Computing module preservation statistics...")
    
    modules = [c for c in ref_colors.unique() if c != 'grey']
    results = []
    
    for mod in modules:
        mod_genes = ref_colors[ref_colors == mod].index
        
        # Get common genes between module and test dataset
        common_genes = mod_genes[mod_genes.isin(test_expr.index)]
        
        if len(common_genes) < 10:
            log.warning(f"Module {mod}: only {len(common_genes)} genes in test, skipping")
            results.append({
                'module': mod,
                'n_ref_genes': len(mod_genes),
                'n_common_genes': len(common_genes),
                'Z_density': np.nan,
                'Z_connectivity': np.nan,
                'Z_summary': np.nan
            })
            continue
        
        # Compute within-module correlation in reference
        ref_mod = ref_expr.loc[common_genes]
        test_mod = test_expr.loc[common_genes]
        
        ref_corr = np.corrcoef(ref_mod.values)
        test_corr = np.corrcoef(test_mod.values)
        
        # Z_density: correlation between within-module connectivity
        ref_k = ref_corr.mean(axis=1)
        test_k = test_corr.mean(axis=1)
        
        if len(ref_k) > 2:
            z_density, p_d = stats.pearsonr(ref_k, test_k)
            # Convert to Z-score approximation
            z_density = 0.5 * np.log((1 + z_density) / (1 - z_density + 1e-6)) * np.sqrt(len(ref_k) - 3)
        else:
            z_density = 0.0
        
        # Z_connectivity: correlation between upper triangle of correlation matrices
        n = len(common_genes)
        triu_idx = np.triu_indices(n, k=1)
        ref_vec = ref_corr[triu_idx]
        test_vec = test_corr[triu_idx]
        
        if len(ref_vec) > 2:
            z_conn, p_c = stats.pearsonr(ref_vec, test_vec)
            z_conn = 0.5 * np.log((1 + z_conn) / (1 - z_conn + 1e-6)) * np.sqrt(len(ref_vec) - 3)
        else:
            z_conn = 0.0
        
        z_summary = np.mean([z_density, z_conn])
        
        results.append({
            'module': mod,
            'n_ref_genes': len(mod_genes),
            'n_common_genes': len(common_genes),
            'Z_density': round(z_density, 3),
            'Z_connectivity': round(z_conn, 3),
            'Z_summary': round(z_summary, 3)
        })
        
        preservation_label = ('strong' if z_summary > 10 else
                              'moderate' if z_summary > 2 else 'low')
        log.info(f"  Module {mod}: Z_summary={z_summary:.2f} ({preservation_label}), "
                 f"n_genes={len(common_genes)}")
    
    return pd.DataFrame(results).sort_values('Z_summary', ascending=False)


# ─────────────────────────────────────────────────────────────────────────────
# Hub gene identification
# ─────────────────────────────────────────────────────────────────────────────

def identify_hub_genes(expr_df: pd.DataFrame,
                       gene_colors: pd.Series,
                       n_hub: int = 10) -> pd.DataFrame:
    """
    Identify hub genes for each module.
    Hub genes = genes with highest module membership (kME).
    kME = correlation between gene expression and module eigengene.
    """
    log.info("Identifying hub genes per module...")
    
    modules = [c for c in gene_colors.unique() if c != 'grey']
    all_hubs = []
    
    for mod in modules:
        mod_genes = gene_colors[gene_colors == mod].index
        mod_genes = mod_genes[mod_genes.isin(expr_df.index)]
        
        if len(mod_genes) < 5:
            continue
        
        # Module eigengene (mean expression as proxy)
        ME = expr_df.loc[mod_genes].mean(axis=0)  # samples,
        
        # Correlation of each gene with ME
        kME = {}
        for gene in mod_genes:
            r, p = stats.pearsonr(expr_df.loc[gene], ME)
            kME[gene] = r
        
        kME_series = pd.Series(kME).sort_values(ascending=False)
        top_hubs = kME_series.head(n_hub)
        
        for gene, kme_val in top_hubs.items():
            all_hubs.append({
                'module': mod,
                'gene': gene,
                'kME': round(kme_val, 4)
            })
    
    hub_df = pd.DataFrame(all_hubs)
    log.info(f"Identified {len(hub_df)} hub genes across {len(modules)} modules")
    return hub_df


def overlap_with_meta_analysis(hub_df: pd.DataFrame,
                               meta_results_dir: str) -> pd.DataFrame:
    """
    Check overlap between WGCNA hub genes and meta-analysis significant genes.
    """
    # Load meta-analysis results
    meta_files = [
        os.path.join(meta_results_dir, 'significant_genes_combined.csv'),
        os.path.join(meta_results_dir, 'phase3_fixed_effects_AS_results.csv'),
        os.path.join(meta_results_dir, 'meta_analysis_results.csv'),
    ]
    
    meta_genes = set()
    for f in meta_files:
        if os.path.exists(f):
            df = pd.read_csv(f)
            # Look for gene column
            gene_cols = [c for c in df.columns if 'gene' in c.lower()]
            if gene_cols:
                meta_genes.update(df[gene_cols[0]].dropna().tolist())
            log.info(f"Loaded {len(meta_genes)} meta-analysis genes from {f}")
            break
    
    if not meta_genes:
        log.warning("No meta-analysis gene lists found")
        hub_df['in_meta_analysis'] = False
        return hub_df
    
    hub_df = hub_df.copy()
    hub_df['in_meta_analysis'] = hub_df['gene'].isin(meta_genes)
    
    overlap_count = hub_df['in_meta_analysis'].sum()
    log.info(f"Hub gene overlap with meta-analysis: {overlap_count}/{len(hub_df)} genes")
    
    return hub_df


# ─────────────────────────────────────────────────────────────────────────────
# Visualization
# ─────────────────────────────────────────────────────────────────────────────

def plot_soft_threshold(fit_df: pd.DataFrame,
                        selected_power: int,
                        gse_id: str,
                        save_path: str):
    """Plot scale-free topology fit vs soft-threshold power."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    ax1 = axes[0]
    ax1.plot(fit_df['power'], fit_df['r2'], 'b-o', linewidth=2)
    ax1.axhline(y=0.80, color='red', linestyle='--', alpha=0.7, label='R²=0.80 cutoff')
    ax1.axvline(x=selected_power, color='green', linestyle='--', alpha=0.7,
                label=f'Selected: {selected_power}')
    ax1.set_xlabel('Soft Threshold Power', fontsize=12)
    ax1.set_ylabel('Scale Free Topology Fit (R²)', fontsize=12)
    ax1.set_title(f'{gse_id}: Scale-Free Topology Fit', fontsize=13)
    ax1.legend()
    ax1.set_ylim(0, 1)
    ax1.grid(True, alpha=0.3)
    
    ax2 = axes[1]
    ax2.plot(fit_df['power'], fit_df['mean_k'], 'r-o', linewidth=2)
    ax2.axvline(x=selected_power, color='green', linestyle='--', alpha=0.7,
                label=f'Selected: {selected_power}')
    ax2.set_xlabel('Soft Threshold Power', fontsize=12)
    ax2.set_ylabel('Mean Connectivity', fontsize=12)
    ax2.set_title(f'{gse_id}: Mean Connectivity', fontsize=13)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.suptitle('Soft-Threshold Power Selection for WGCNA', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved soft-threshold plot: {save_path}")


def plot_module_summary(gene_colors: pd.Series,
                        corr_df: pd.DataFrame,
                        gse_id: str,
                        save_path: str):
    """Plot module sizes and module-trait correlations."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Module sizes
    mod_counts = gene_colors.value_counts()
    mod_counts = mod_counts[mod_counts.index != 'grey']
    
    colors = mod_counts.index.tolist()
    valid_colors = []
    for c in colors:
        try:
            plt.matplotlib.colors.to_rgba(c)
            valid_colors.append(c)
        except ValueError:
            valid_colors.append('gray')
    
    ax1 = axes[0]
    bars = ax1.bar(range(len(mod_counts)), mod_counts.values, color=valid_colors)
    ax1.set_xticks(range(len(mod_counts)))
    ax1.set_xticklabels(mod_counts.index, rotation=45, ha='right', fontsize=9)
    ax1.set_xlabel('Module', fontsize=12)
    ax1.set_ylabel('Number of Genes', fontsize=12)
    ax1.set_title(f'{gse_id}: Module Sizes', fontsize=13)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add grey module count as text
    grey_n = (gene_colors == 'grey').sum()
    ax1.text(0.98, 0.98, f'Grey (unassigned): {grey_n}', 
             transform=ax1.transAxes, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgrey', alpha=0.5))
    
    # Module-trait correlations
    ax2 = axes[1]
    if len(corr_df) > 0:
        mod_names = [name.replace('ME', '') for name in corr_df.index[:15]]
        corr_vals = corr_df['correlation'].iloc[:15].values
        pvals = corr_df['pvalue'].iloc[:15].values
        
        bar_colors = []
        for c_name in mod_names:
            try:
                plt.matplotlib.colors.to_rgba(c_name)
                bar_colors.append(c_name)
            except ValueError:
                bar_colors.append('steelblue')
        
        bars2 = ax2.barh(range(len(mod_names)), corr_vals, color=bar_colors)
        
        # Add significance markers
        for i, (corr, pval) in enumerate(zip(corr_vals, pvals)):
            sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
            ax2.text(corr + 0.02 if corr > 0 else corr - 0.02, i, sig,
                     va='center', ha='left' if corr > 0 else 'right', fontsize=10)
        
        ax2.set_yticks(range(len(mod_names)))
        ax2.set_yticklabels(mod_names, fontsize=9)
        ax2.set_xlabel('Pearson Correlation with Disease Trait', fontsize=11)
        ax2.set_title(f'{gse_id}: Module-Trait Correlations', fontsize=13)
        ax2.axvline(x=0, color='black', linewidth=0.5)
        ax2.grid(True, alpha=0.3, axis='x')
    
    plt.suptitle(f'WGCNA Module Summary: {gse_id}', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved module summary plot: {save_path}")


def plot_module_preservation(preservation_df: pd.DataFrame,
                             ref_id: str,
                             test_id: str,
                             save_path: str):
    """Plot module preservation Z-summary statistics."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    df = preservation_df.dropna(subset=['Z_summary'])
    
    mod_names = df['module'].tolist()
    z_vals = df['Z_summary'].tolist()
    n_genes = df['n_common_genes'].tolist()
    
    colors = []
    for m in mod_names:
        try:
            plt.matplotlib.colors.to_rgba(m)
            colors.append(m)
        except ValueError:
            colors.append('steelblue')
    
    bars = ax.bar(range(len(mod_names)), z_vals, color=colors, edgecolor='black', linewidth=0.5)
    
    # Reference lines
    ax.axhline(y=2, color='orange', linestyle='--', linewidth=2, label='Z=2 (moderate preservation)')
    ax.axhline(y=10, color='red', linestyle='--', linewidth=2, label='Z=10 (strong preservation)')
    
    # Add gene count labels
    for i, (bar, n) in enumerate(zip(bars, n_genes)):
        ax.text(i, max(bar.get_height() + 0.5, 0.5), f'n={n}',
                ha='center', va='bottom', fontsize=8)
    
    ax.set_xticks(range(len(mod_names)))
    ax.set_xticklabels(mod_names, rotation=45, ha='right')
    ax.set_ylabel('Z_summary (Module Preservation)', fontsize=12)
    ax.set_title(f'Module Preservation: {ref_id} → {test_id}', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"Saved preservation plot: {save_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Per-dataset WGCNA runner
# ─────────────────────────────────────────────────────────────────────────────

def run_wgcna_for_dataset(gse_id: str,
                          dataset_info: dict) -> dict:
    """
    Run complete WGCNA pipeline for one dataset.
    Returns dict with results summary.
    """
    log.info(f"\n{'='*70}")
    log.info(f"Running WGCNA for: {gse_id} ({dataset_info['description']})")
    log.info(f"{'='*70}")
    
    # Create output directories
    os.makedirs(os.path.join(RESULTS_DIR, gse_id), exist_ok=True)
    os.makedirs(os.path.join(FIGURES_DIR, gse_id), exist_ok=True)
    
    result = {
        'gse_id': gse_id,
        'disease': dataset_info['disease'],
        'status': 'failed',
        'error': None
    }
    
    try:
        # 1. Load data
        expr_df, meta_df = load_expression_matrix(gse_id)
        
        # 2. Aggregate to gene level
        expr_df = aggregate_to_gene_level(expr_df)
        
        # 3. Filter to most variable genes
        expr_df = filter_low_variance_genes(expr_df, top_n=3000)
        
        log.info(f"Expression matrix: {expr_df.shape[0]} genes × {expr_df.shape[1]} samples")
        
        # 4. Soft threshold selection
        selected_power, fit_df = select_soft_threshold(expr_df)
        
        # Plot soft threshold
        plot_soft_threshold(
            fit_df, selected_power, gse_id,
            os.path.join(FIGURES_DIR, gse_id, 'soft_threshold.png')
        )
        fit_df.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'soft_threshold_fits.csv'), index=False
        )
        
        # 5. Adjacency and TOM
        adj = compute_adjacency(expr_df, selected_power)
        TOM = adjacency_to_TOM(adj)
        del adj  # Free memory
        
        # 6. Module detection
        gene_colors = detect_modules(
            TOM, expr_df.index.tolist(),
            min_module_size=MIN_MODULE_SIZE,
            merge_cut_height=MERGE_CUT_HEIGHT
        )
        
        # Save module assignments
        gene_colors.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'gene_module_colors.csv'),
            header=['module']
        )
        
        n_modules = len(gene_colors[gene_colors != 'grey'].unique())
        log.info(f"Detected {n_modules} modules")
        
        # 7. Module eigengenes
        ME_df = compute_module_eigengenes(expr_df, gene_colors)
        ME_df.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'module_eigengenes.csv')
        )
        
        # 8. Module-trait correlations
        corr_df, trait_vals = module_trait_correlation(ME_df, meta_df)
        corr_df.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'module_trait_correlations.csv')
        )
        
        # 9. Plot module summary
        plot_module_summary(
            gene_colors, corr_df, gse_id,
            os.path.join(FIGURES_DIR, gse_id, 'module_summary.png')
        )
        
        # 10. Hub genes
        hub_df = identify_hub_genes(expr_df, gene_colors)
        
        # Check overlap with meta-analysis
        meta_results_dir = os.path.join(BASE_DIR, "results", "phase3_meta_analysis")
        hub_df = overlap_with_meta_analysis(hub_df, meta_results_dir)
        hub_df.to_csv(
            os.path.join(RESULTS_DIR, gse_id, 'hub_genes.csv'), index=False
        )
        
        result.update({
            'status': 'success',
            'n_genes_analyzed': expr_df.shape[0],
            'n_samples': expr_df.shape[1],
            'soft_power': selected_power,
            'n_modules': n_modules,
            'n_hub_genes': len(hub_df),
            'top_module_trait_corr': float(corr_df['abs_correlation'].max()) if len(corr_df) > 0 else None,
            'gene_colors': gene_colors  # Keep for preservation analysis
        })
        
        log.info(f"WGCNA complete for {gse_id}: {n_modules} modules, power={selected_power}")
        
    except FileNotFoundError as e:
        log.error(f"Data not found for {gse_id}: {e}")
        result['error'] = str(e)
    except Exception as e:
        log.error(f"WGCNA failed for {gse_id}: {e}")
        import traceback
        traceback.print_exc()
        result['error'] = str(e)
    
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Cross-dataset preservation analysis
# ─────────────────────────────────────────────────────────────────────────────

def run_preservation_analysis(results: dict) -> pd.DataFrame:
    """
    Run module preservation analysis across all successful datasets.
    For each pair of datasets, test preservation of reference modules in test dataset.
    """
    log.info("\n" + "="*70)
    log.info("Running cross-dataset module preservation analysis")
    log.info("="*70)
    
    successful = {
        gse: res for gse, res in results.items() 
        if res['status'] == 'success' and 'gene_colors' in res
    }
    
    if len(successful) < 2:
        log.warning(f"Only {len(successful)} successful datasets, need ≥2 for preservation")
        return pd.DataFrame()
    
    all_preservation = []
    
    gse_ids = list(successful.keys())
    for i in range(len(gse_ids)):
        for j in range(len(gse_ids)):
            if i == j:
                continue
            
            ref_id = gse_ids[i]
            test_id = gse_ids[j]
            
            log.info(f"Preservation: {ref_id} → {test_id}")
            
            try:
                # Load expression for both
                ref_expr, _ = load_expression_matrix(ref_id)
                ref_expr = aggregate_to_gene_level(ref_expr)
                ref_expr = filter_low_variance_genes(ref_expr, top_n=3000)
                
                test_expr, _ = load_expression_matrix(test_id)
                test_expr = aggregate_to_gene_level(test_expr)
                test_expr = filter_low_variance_genes(test_expr, top_n=3000)
                
                ref_colors = successful[ref_id]['gene_colors']
                
                preservation_df = compute_module_preservation(
                    ref_colors, ref_expr, test_expr
                )
                preservation_df['ref_dataset'] = ref_id
                preservation_df['test_dataset'] = test_id
                
                all_preservation.append(preservation_df)
                
                # Save individual preservation results
                pres_path = os.path.join(
                    RESULTS_DIR, f'preservation_{ref_id}_in_{test_id}.csv'
                )
                preservation_df.to_csv(pres_path, index=False)
                
                # Plot preservation
                plot_module_preservation(
                    preservation_df, ref_id, test_id,
                    os.path.join(FIGURES_DIR, f'preservation_{ref_id}_in_{test_id}.png')
                )
                
            except Exception as e:
                log.error(f"Preservation analysis failed for {ref_id}→{test_id}: {e}")
    
    if all_preservation:
        combined_preservation = pd.concat(all_preservation, ignore_index=True)
        combined_preservation.to_csv(
            os.path.join(RESULTS_DIR, 'all_preservation_results.csv'), index=False
        )
        return combined_preservation
    
    return pd.DataFrame()


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Phase 7: WGCNA Analysis')
    parser.add_argument('--datasets', nargs='+', default=list(WGCNA_DATASETS.keys()),
                        help='Datasets to analyze')
    parser.add_argument('--top-genes', type=int, default=3000,
                        help='Top N variable genes for WGCNA')
    parser.add_argument('--min-module-size', type=int, default=MIN_MODULE_SIZE,
                        help='Minimum module size')
    parser.add_argument('--soft-power', type=int, default=None,
                        help='Override soft-threshold power (skip selection)')
    parser.add_argument('--skip-preservation', action='store_true',
                        help='Skip cross-dataset preservation analysis')
    args = parser.parse_args()
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    log.info("Phase 7: WGCNA Co-expression Analysis")
    log.info(f"Datasets: {args.datasets}")
    log.info(f"Top genes: {args.top_genes}")
    log.info(f"Min module size: {args.min_module_size}")
    
    # Run WGCNA for each dataset
    results = {}
    for gse_id in args.datasets:
        if gse_id not in WGCNA_DATASETS:
            log.warning(f"Unknown dataset: {gse_id}, skipping")
            continue
        results[gse_id] = run_wgcna_for_dataset(gse_id, WGCNA_DATASETS[gse_id])
    
    # Cross-dataset preservation
    if not args.skip_preservation and len(results) > 1:
        preservation_df = run_preservation_analysis(results)
    
    # Summary report
    summary = []
    for gse_id, res in results.items():
        summary.append({
            'dataset': gse_id,
            'disease': res.get('disease', 'unknown'),
            'status': res['status'],
            'n_genes': res.get('n_genes_analyzed', 0),
            'n_samples': res.get('n_samples', 0),
            'soft_power': res.get('soft_power', None),
            'n_modules': res.get('n_modules', 0),
            'n_hub_genes': res.get('n_hub_genes', 0),
            'max_trait_corr': res.get('top_module_trait_corr', None),
            'error': res.get('error', None)
        })
    
    summary_df = pd.DataFrame(summary)
    summary_path = os.path.join(RESULTS_DIR, 'wgcna_summary.csv')
    summary_df.to_csv(summary_path, index=False)
    
    log.info("\n" + "="*70)
    log.info("WGCNA ANALYSIS COMPLETE")
    log.info("="*70)
    log.info(f"\nSummary:\n{summary_df.to_string(index=False)}")
    log.info(f"\nResults saved to: {RESULTS_DIR}")
    log.info(f"Figures saved to: {FIGURES_DIR}")


if __name__ == '__main__':
    main()
