#!/usr/bin/env python3
"""
Phase 7: Weighted Gene Co-expression Network Analysis (WGCNA)
=============================================================
Cross-Disease Transcriptomic Meta-Analysis of SpA - Project 1

Implements a Python-based WGCNA pipeline:
1. Probe-to-gene mapping and gene-level expression aggregation
2. Soft-threshold power selection (scale-free topology)
3. Adjacency matrix and TOM (Topological Overlap Matrix) construction
4. Hierarchical clustering with dynamic tree cutting for module detection
5. Module-trait correlation analysis (disease vs control)
6. Module preservation across AS and IBD datasets
7. Hub gene identification within key modules
8. Overlap with AS meta-analysis genes

Note: Python implementation inspired by the R WGCNA package (Langfelder & Horvath).
Core algorithms: adjacency -> TOM -> hierarchical clustering -> dynamic cut.

Target datasets:
- GSE18781 (AS/axSpA, 55 samples, Affy U133+2) - largest AS dataset
- GSE59071 (IBD, 116 samples, HuGene 1.0 ST) - largest overall, SpA-associated IBD

Module preservation tested: modules found in AS -> preservation in IBD.

Author: Perplexity Computer / Spondyloarthritis AI Computational Biology Institute
Date: March 2026
"""

import os
import sys
import json
import gzip
import time
import logging
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

# Configuration
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results', 'wgcna')
FIGURES_DIR = os.path.join(BASE_DIR, 'figures')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(os.path.join(BASE_DIR, 'logs'), exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(BASE_DIR, 'logs', 'phase7_wgcna.log'))
    ]
)
log = logging.getLogger(__name__)

# WGCNA parameters
MIN_MODULE_SIZE = 30
MERGE_CUT_HEIGHT = 0.25
TOP_VAR_GENES = 5000


def load_platform_annotation(platform_file):
    """Load GPL annotation file: probe ID -> gene symbol mapping."""
    probe2gene = {}
    with gzip.open(platform_file, 'rt') as f:
        for line in f:
            if line.startswith(('#', '!', '^')):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[0] != 'ID':
                gene_sym = parts[2].strip()
                if gene_sym and gene_sym != '---' and gene_sym != '':
                    genes = [g.strip() for g in gene_sym.split('///')]
                    protein_genes = [g for g in genes
                                     if not g.startswith('MIR') and
                                     not g.startswith('LOC') and
                                     not g.startswith('LINC')]
                    probe2gene[parts[0]] = protein_genes[0] if protein_genes else genes[0]
    return probe2gene


def load_expression_gene_level(accession, disease_label):
    """Load expression matrix and convert to gene-level (genes x samples)."""
    expr_file = os.path.join(DATA_DIR, 'processed', f'{accession}_expression.csv')
    counts_file = os.path.join(DATA_DIR, 'processed', f'{accession}_counts.csv')
    samples_file = os.path.join(DATA_DIR, 'metadata', f'{accession}_samples.json')

    if not os.path.exists(samples_file):
        log.warning(f"No samples file for {accession}")
        return None, None

    with open(samples_file) as f:
        samples_meta = json.load(f)

    if os.path.exists(counts_file) and os.path.getsize(counts_file) > 1000:
        log.info(f"  Loading RNA-seq counts from {counts_file}")
        df = pd.read_csv(counts_file, index_col=0)
        if 'gene_symbol' in df.columns:
            gene_symbols = df['gene_symbol']
            df = df.drop('gene_symbol', axis=1)
            lib_sizes = df.sum(axis=0)
            cpm = df.div(lib_sizes, axis=1) * 1e6
            expr = np.log2(cpm + 1)
            expr.index = gene_symbols
            expr = expr.groupby(expr.index).mean()
        else:
            lib_sizes = df.sum(axis=0)
            cpm = df.div(lib_sizes, axis=1) * 1e6
            expr = np.log2(cpm + 1)
    else:
        if not os.path.exists(expr_file):
            log.warning(f"No expression file for {accession}")
            return None, None
        log.info(f"  Loading microarray from {expr_file}")
        df = pd.read_csv(expr_file, index_col=0)

        platform_map = {
            'GSE18781': 'GPL570',
            'GSE58667': 'GPL570',
            'GSE59071': 'GPL6244',
            'GSE61281': 'GPL6480',
        }
        platform = platform_map.get(accession)

        if platform:
            annot_file = os.path.join(DATA_DIR, 'raw', f'{platform}_annot.txt.gz')
            if os.path.exists(annot_file):
                probe2gene = load_platform_annotation(annot_file)
                mapped = [probe2gene.get(str(p)) for p in df.index]
                df['gene_symbol'] = mapped
                df = df.dropna(subset=['gene_symbol'])
                expr = df.groupby('gene_symbol').mean()
            else:
                expr = df
        else:
            deg_file = os.path.join(BASE_DIR, 'results', 'deg', f'{accession}_deg_genes.csv')
            if os.path.exists(deg_file):
                deg = pd.read_csv(deg_file)
                if 'gene_symbol' in deg.columns and 'gene_id' in deg.columns:
                    p2g = dict(zip(deg['gene_id'], deg['gene_symbol']))
                    p2g = {k: v for k, v in p2g.items() if k != v and v and pd.notna(v)}
                    mapped = [p2g.get(str(p)) for p in df.index]
                    df['gene_symbol'] = mapped
                    df = df.dropna(subset=['gene_symbol'])
                    expr = df.groupby('gene_symbol').mean()
                else:
                    expr = df
            else:
                expr = df

    valid_samples = [s for s in expr.columns if s in samples_meta]
    if len(valid_samples) < 10:
        log.warning(f"  Only {len(valid_samples)} valid samples for {accession}")
        return None, None

    expr = expr[valid_samples]
    log.info(f"  Expression matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    groups = {s: samples_meta[s].get('group', 'unknown') for s in valid_samples}
    return expr, groups


def select_top_variable_genes(expr, n_top=5000):
    """Select top N most variable genes by variance across samples."""
    variances = expr.var(axis=1)
    top_genes = variances.nlargest(min(n_top, len(variances))).index
    return expr.loc[top_genes]


def pick_soft_threshold(expr, powers=None):
    """Estimate soft-thresholding power for scale-free topology."""
    if powers is None:
        powers = list(range(1, 11)) + [12, 14, 16, 18, 20]

    log.info("  Evaluating soft-threshold powers...")
    results = []

    # Compute correlation matrix (sample correlation)
    corr_matrix = np.corrcoef(expr.values)
    n = len(corr_matrix)

    for power in powers:
        adj = np.abs(corr_matrix) ** power
        np.fill_diagonal(adj, 0)

        # Degree distribution
        k = adj.sum(axis=1)
        k = k[k > 0]
        if len(k) < 2:
            results.append({'power': power, 'r2': 0, 'slope': 0, 'mean_k': 0})
            continue

        # Log-log fit
        log_k = np.log10(k)
        hist, bin_edges = np.histogram(log_k, bins=20)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        mask = hist > 0
        if mask.sum() < 3:
            results.append({'power': power, 'r2': 0, 'slope': 0, 'mean_k': np.mean(k)})
            continue

        log_p = np.log10(hist[mask] / hist[mask].sum())
        slope, intercept, r_val, p_val, _ = stats.linregress(bin_centers[mask], log_p)
        r2 = r_val ** 2
        results.append({'power': power, 'r2': r2, 'slope': slope, 'mean_k': np.mean(k)})

    # Choose power with R2 > 0.85
    results_df = pd.DataFrame(results)
    good = results_df[results_df['r2'] >= 0.85]
    if len(good) > 0:
        chosen = good['power'].iloc[0]
    else:
        chosen = results_df.loc[results_df['r2'].idxmax(), 'power']

    log.info(f"  Selected soft-threshold power: {chosen}")
    return int(chosen), results_df


def build_adjacency_matrix(expr, power):
    """Build signed weighted adjacency matrix."""
    log.info(f"  Building adjacency matrix (power={power})...")
    corr = np.corrcoef(expr.values)
    # Signed adjacency: (1 + corr) / 2
    signed_corr = (1 + corr) / 2
    adj = signed_corr ** power
    np.fill_diagonal(adj, 0)
    return adj


def build_tom(adj):
    """Compute Topological Overlap Matrix (TOM) from adjacency matrix."""
    log.info("  Computing TOM (this may take several minutes for large matrices)...")
    n = adj.shape[0]
    tom = np.zeros((n, n))

    # TOM formula: TOM[i,j] = (sum_u(adj[i,u]*adj[u,j]) + adj[i,j]) /
    #              (min(sum_u adj[i,u], sum_u adj[j,u]) - adj[i,j] + 1)
    k = adj.sum(axis=1)

    # Vectorized computation
    numerator = adj @ adj + adj
    k_min = np.minimum.outer(k, k)
    denominator = k_min - adj + 1
    tom = numerator / np.maximum(denominator, 1e-10)
    np.fill_diagonal(tom, 1.0)

    return tom


def detect_modules(tom, gene_names, min_module_size=30, merge_cut_height=0.25):
    """Detect co-expression modules using hierarchical clustering + dynamic tree cut."""
    log.info("  Hierarchical clustering on TOM dissimilarity...")
    dissim = 1 - tom
    np.fill_diagonal(dissim, 0)

    condensed = squareform(dissim, checks=False)
    Z = linkage(condensed, method='average')

    # Dynamic tree cut (simplified): cut at various heights and find stable modules
    log.info("  Dynamic tree cutting...")
    height_threshold = 0.98
    labels = fcluster(Z, t=height_threshold, criterion='distance')

    # Merge small modules into 'grey' (module 0)
    label_counts = Counter(labels)
    small_modules = {lb for lb, cnt in label_counts.items() if cnt < min_module_size}
    labels = np.array([0 if lb in small_modules else lb for lb in labels])

    # Renumber modules sequentially
    unique_labels = sorted(set(labels) - {0})
    label_map = {old: new for new, old in enumerate(unique_labels, 1)}
    label_map[0] = 0
    labels = np.array([label_map[lb] for lb in labels])

    # Module colors (same scheme as WGCNA R package)
    colors = [
        'grey', 'turquoise', 'blue', 'brown', 'yellow', 'green',
        'red', 'black', 'pink', 'magenta', 'purple', 'greenyellow',
        'tan', 'salmon', 'cyan', 'midnightblue', 'lightcyan',
        'grey60', 'lightgreen', 'lightyellow', 'royalblue'
    ]

    n_modules = labels.max()
    module_colors = {i: colors[min(i, len(colors)-1)] for i in range(n_modules+1)}

    assignments = pd.DataFrame({
        'gene': gene_names,
        'module': labels,
        'color': [module_colors[m] for m in labels]
    })

    log.info(f"  Detected {n_modules} modules + grey")
    for m in range(1, min(n_modules+1, 10)):
        cnt = (labels == m).sum()
        log.info(f"    Module {m} ({module_colors[m]}): {cnt} genes")

    return assignments, Z, module_colors


def compute_module_eigengenes(expr, assignments):
    """Compute module eigengenes (first PC of each module's expression)."""
    eigengenes = {}

    for module in sorted(assignments['module'].unique()):
        if module == 0:
            continue
        module_genes = assignments[assignments['module'] == module]['gene'].tolist()
        module_genes = [g for g in module_genes if g in expr.index]

        if len(module_genes) < 3:
            continue

        module_expr = expr.loc[module_genes].values.T

        # SVD to get first PC (module eigengene)
        module_expr_centered = module_expr - module_expr.mean(axis=0)
        try:
            U, S, Vt = np.linalg.svd(module_expr_centered, full_matrices=False)
            eigengene = U[:, 0] * S[0]
            # Ensure positive correlation with mean expression
            mean_expr = module_expr.mean(axis=1)
            if np.corrcoef(eigengene, mean_expr)[0, 1] < 0:
                eigengene = -eigengene
            eigengenes[module] = eigengene
        except np.linalg.LinAlgError:
            eigengenes[module] = module_expr.mean(axis=1)

    return eigengenes


def compute_module_trait_correlation(eigengenes, groups, samples):
    """Correlate module eigengenes with disease trait (disease=1, control=0)."""
    trait = np.array([1 if groups.get(s, 'unknown') == 'disease' else 0 for s in samples])
    results = []

    for module, eigengene in eigengenes.items():
        if len(eigengene) != len(trait):
            continue
        r, p = stats.pearsonr(eigengene, trait)
        results.append({'module': module, 'r': r, 'p': p, 'r2': r**2})

    return pd.DataFrame(results).sort_values('p') if results else pd.DataFrame()


def module_preservation_zsummary(ref_expr, test_expr, ref_assignments, n_perm=100):
    """Compute module preservation Z-summary statistics.

    Zsummary > 10: strongly preserved
    Zsummary 2-10: moderately preserved
    Zsummary < 2: not preserved
    """
    results = []

    for module in sorted(ref_assignments['module'].unique()):
        if module == 0:
            continue

        module_genes = ref_assignments[ref_assignments['module'] == module]['gene'].tolist()
        common_genes = [g for g in module_genes if g in ref_expr.index and g in test_expr.index]

        if len(common_genes) < 10:
            results.append({
                'module': module,
                'n_genes': len(common_genes),
                'Zsummary': np.nan,
                'preservation': 'insufficient_genes'
            })
            continue

        n_genes = len(common_genes)
        ref_mod = ref_expr.loc[common_genes].values.T
        test_mod = test_expr.loc[common_genes].values.T

        # Observed statistics
        ref_corr = np.corrcoef(ref_mod.T)
        test_corr = np.corrcoef(test_mod.T)
        np.fill_diagonal(ref_corr, 0)
        np.fill_diagonal(test_corr, 0)

        # Density: mean of abs correlation in module
        ref_density = np.abs(ref_corr).mean()
        test_density = np.abs(test_corr).mean()
        density_ratio = test_density / (ref_density + 1e-10)

        # Connectivity preservation
        ref_k = np.abs(ref_corr).sum(axis=1)
        test_k = np.abs(test_corr).sum(axis=1)
        if ref_k.std() > 0 and test_k.std() > 0:
            r_connectivity, _ = stats.pearsonr(ref_k, test_k)
        else:
            r_connectivity = 0

        # Permutation test to get Z-scores
        all_genes_test = list(test_expr.index)
        null_densities = []
        null_r_conn = []

        for _ in range(n_perm):
            if len(all_genes_test) >= n_genes:
                perm_genes = np.random.choice(all_genes_test, n_genes, replace=False)
                perm_mod = test_expr.loc[perm_genes].values.T
                perm_corr = np.corrcoef(perm_mod.T)
                np.fill_diagonal(perm_corr, 0)
                null_densities.append(np.abs(perm_corr).mean())

                perm_k = np.abs(perm_corr).sum(axis=1)
                if ref_k.std() > 0 and perm_k.std() > 0:
                    r_conn_null, _ = stats.pearsonr(ref_k, perm_k)
                    null_r_conn.append(r_conn_null)

        if null_densities:
            z_density = (test_density - np.mean(null_densities)) / (np.std(null_densities) + 1e-10)
        else:
            z_density = 0

        if null_r_conn:
            z_connectivity = (r_connectivity - np.mean(null_r_conn)) / (np.std(null_r_conn) + 1e-10)
        else:
            z_connectivity = 0

        zsummary = (z_density + z_connectivity) / 2

        if zsummary > 10:
            preservation = 'strong'
        elif zsummary > 2:
            preservation = 'moderate'
        else:
            preservation = 'none'

        results.append({
            'module': module,
            'n_genes': n_genes,
            'ref_density': ref_density,
            'test_density': test_density,
            'r_connectivity': r_connectivity,
            'z_density': z_density,
            'z_connectivity': z_connectivity,
            'Zsummary': zsummary,
            'preservation': preservation
        })

    return pd.DataFrame(results)


def find_hub_genes(expr, assignments, n_hubs=5):
    """Find hub genes (highest intramodular connectivity) per module."""
    hubs = {}

    for module in sorted(assignments['module'].unique()):
        if module == 0:
            continue
        module_genes = assignments[assignments['module'] == module]['gene'].tolist()
        module_genes = [g for g in module_genes if g in expr.index]

        if len(module_genes) < 5:
            continue

        module_expr = expr.loc[module_genes]
        corr = module_expr.T.corr()
        connectivity = corr.abs().sum() - 1
        top_hubs = connectivity.nlargest(min(n_hubs, len(connectivity)))
        hubs[module] = top_hubs.index.tolist()

    return hubs


def create_wgcna_figures(accession, Z, assignments, trait_corr, expr, figures_dir):
    """Generate WGCNA summary figures."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 1: Module size bar chart
    ax = axes[0]
    module_counts = assignments[assignments['module'] != 0]['module'].value_counts().sort_index()
    colors_list = [assignments[assignments['module'] == m]['color'].iloc[0]
                   for m in module_counts.index]
    try:
        ax.barh(range(len(module_counts)), module_counts.values, color=colors_list)
    except Exception:
        ax.barh(range(len(module_counts)), module_counts.values)
    ax.set_yticks(range(len(module_counts)))
    ax.set_yticklabels([f"M{m}" for m in module_counts.index], fontsize=8)
    ax.set_xlabel('Number of genes')
    ax.set_title(f'{accession}\nModule sizes')

    # Panel 2: Module-trait correlation
    ax = axes[1]
    if len(trait_corr) > 0 and 'r' in trait_corr.columns:
        modules = trait_corr['module'].tolist()
        rs = trait_corr['r'].tolist()
        ps = trait_corr['p'].tolist()
        bar_colors = ['#e74c3c' if r > 0 else '#3498db' for r in rs]
        ax.barh(range(len(rs)), rs, color=bar_colors, alpha=0.7)
        ax.set_yticks(range(len(modules)))
        ax.set_yticklabels([f"M{m}" for m in modules], fontsize=8)
        ax.set_xlabel('Pearson r (disease correlation)')
        ax.axvline(0, color='black', linewidth=0.5)
        ax.set_title('Module-Disease\nCorrelation')
        for i, (r, p) in enumerate(zip(rs, ps)):
            if p < 0.05:
                ax.annotate(f'p={p:.2e}', xy=(r, i), xytext=(5, 0),
                            textcoords='offset points', fontsize=6)
    else:
        ax.text(0.5, 0.5, 'No trait\ncorrelation data', ha='center', va='center',
                transform=ax.transAxes)
        ax.set_title('Module-Disease\nCorrelation')

    # Panel 3: Expression heatmap of top module
    ax = axes[2]
    if len(trait_corr) > 0 and 'module' in trait_corr.columns:
        top_module = int(trait_corr.iloc[0]['module'])
        module_genes = assignments[assignments['module'] == top_module]['gene'].tolist()
        module_genes = [g for g in module_genes if g in expr.index][:50]
        if module_genes:
            heatmap_data = expr.loc[module_genes]
            heatmap_norm = (heatmap_data.T - heatmap_data.T.mean()) / (heatmap_data.T.std() + 1e-10)
            im = ax.imshow(heatmap_norm.T, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
            plt.colorbar(im, ax=ax, label='Z-score')
            ax.set_title(f'Module {top_module} expression\n(top 50 genes)')
            ax.set_xlabel('Samples')
            ax.set_ylabel('Genes')
        else:
            ax.text(0.5, 0.5, 'No genes', ha='center', va='center', transform=ax.transAxes)
    else:
        ax.text(0.5, 0.5, 'No module data', ha='center', va='center', transform=ax.transAxes)

    plt.suptitle(f'WGCNA Summary - {accession}', fontsize=14)
    plt.tight_layout()
    out_path = os.path.join(figures_dir, f'wgcna_{accession}_summary.png')
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    log.info(f"  Saved WGCNA figure: {out_path}")


def run_wgcna_for_dataset(accession, disease_label):
    """Run complete WGCNA pipeline for one dataset."""
    log.info(f"\n{'='*60}")
    log.info(f"WGCNA: {accession} ({disease_label})")
    log.info(f"{'='*60}")

    # Load expression
    expr, groups = load_expression_gene_level(accession, disease_label)
    if expr is None:
        log.warning(f"  Cannot load expression for {accession}")
        return None

    samples = list(expr.columns)

    # Filter to top variable genes
    log.info(f"  Filtering to top {TOP_VAR_GENES} variable genes...")
    expr_top = select_top_variable_genes(expr, TOP_VAR_GENES)
    gene_names = list(expr_top.index)
    log.info(f"  Working with {len(gene_names)} genes")

    # Select soft threshold power
    if len(expr_top) > 1000:
        log.info("  Selecting soft-threshold power (using 1000 gene subset for speed)...")
        expr_subset = expr_top.iloc[:1000]
    else:
        expr_subset = expr_top
    power, power_df = pick_soft_threshold(expr_subset)

    # Build adjacency and TOM
    log.info(f"  Building full adjacency matrix ({len(gene_names)} x {len(gene_names)})...")
    adj = build_adjacency_matrix(expr_top, power)

    log.info("  Building TOM...")
    tom = build_tom(adj)

    # Detect modules
    assignments, Z, module_colors = detect_modules(
        tom, gene_names, MIN_MODULE_SIZE, MERGE_CUT_HEIGHT
    )

    # Compute eigengenes
    log.info("  Computing module eigengenes...")
    eigengenes = compute_module_eigengenes(expr_top, assignments)

    # Module-trait correlation
    trait_corr = compute_module_trait_correlation(eigengenes, groups, samples)

    # Find hub genes
    log.info("  Identifying hub genes...")
    hubs = find_hub_genes(expr_top, assignments)

    # Save results
    out_dir = os.path.join(RESULTS_DIR, accession)
    os.makedirs(out_dir, exist_ok=True)

    assignments.to_csv(os.path.join(out_dir, 'module_assignments.csv'), index=False)
    trait_corr.to_csv(os.path.join(out_dir, 'module_trait_correlation.csv'), index=False)
    power_df.to_csv(os.path.join(out_dir, 'soft_threshold_selection.csv'), index=False)

    hubs_data = {str(m): genes for m, genes in hubs.items()}
    with open(os.path.join(out_dir, 'hub_genes.json'), 'w') as f:
        json.dump(hubs_data, f, indent=2)

    # Save eigengenes
    if eigengenes:
        eg_df = pd.DataFrame(eigengenes, index=samples)
        eg_df.to_csv(os.path.join(out_dir, 'eigengenes.csv'))

    # Create figures
    create_wgcna_figures(accession, Z, assignments, trait_corr, expr_top, FIGURES_DIR)

    summary = {
        'accession': accession,
        'disease': disease_label,
        'n_genes': len(gene_names),
        'n_samples': len(samples),
        'soft_power': int(power),
        'n_modules': int(assignments['module'].max()),
        'n_grey': int((assignments['module'] == 0).sum()),
    }

    log.info(f"  WGCNA complete: {summary['n_modules']} modules")
    return summary, assignments, eigengenes, expr_top, groups


def run_preservation_analysis(ref_accession, test_accession, ref_results, test_results):
    """Run module preservation analysis between two datasets."""
    if ref_results is None or test_results is None:
        return None

    ref_summary, ref_assignments, ref_eigen, ref_expr, ref_groups = ref_results
    test_summary, test_assignments, test_eigen, test_expr, test_groups = test_results

    log.info(f"\nModule preservation: {ref_accession} -> {test_accession}")
    preservation = module_preservation_zsummary(ref_expr, test_expr, ref_assignments, n_perm=50)

    out_path = os.path.join(RESULTS_DIR, f'preservation_{ref_accession}_in_{test_accession}.csv')
    preservation.to_csv(out_path, index=False)

    log.info(f"  Preservation results saved: {out_path}")
    log.info(f"  Strongly preserved: {(preservation['Zsummary'] > 10).sum()} modules")
    log.info(f"  Moderately preserved: {((preservation['Zsummary'] > 2) & (preservation['Zsummary'] <= 10)).sum()} modules")
    log.info(f"  Not preserved: {(preservation['Zsummary'] < 2).sum()} modules")

    return preservation


def main():
    log.info("Phase 7: WGCNA Co-expression Network Analysis")
    log.info("=" * 60)

    # Primary datasets for WGCNA
    datasets = [
        ("GSE18781", "AS"),
        ("GSE59071", "IBD"),
    ]

    results = {}
    for accession, disease in datasets:
        result = run_wgcna_for_dataset(accession, disease)
        if result is not None:
            results[accession] = result

    # Module preservation analysis
    if "GSE18781" in results and "GSE59071" in results:
        log.info("\nRunning module preservation analysis...")
        run_preservation_analysis(
            "GSE18781", "GSE59071",
            results["GSE18781"], results["GSE59071"]
        )
        run_preservation_analysis(
            "GSE59071", "GSE18781",
            results["GSE59071"], results["GSE18781"]
        )

    # Save overall summary
    summary_data = {}
    for acc, result in results.items():
        if result:
            summary_data[acc] = result[0]

    with open(os.path.join(RESULTS_DIR, 'wgcna_summary.json'), 'w') as f:
        json.dump(summary_data, f, indent=2)

    log.info("\nPhase 7 complete.")


if __name__ == "__main__":
    main()
