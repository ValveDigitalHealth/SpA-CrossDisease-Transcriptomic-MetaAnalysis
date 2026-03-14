#!/usr/bin/env python3
"""
Phase 4: Cross-Disease Transcriptomic Meta-Analysis
Combines DEG results across SpA-spectrum diseases using:
1. Fisher's method for p-value combination
2. Vote counting (consistent direction of change)
3. Effect size integration (weighted mean log2FC)

Uses a more permissive threshold (nominal p < 0.05) for the meta-analysis
to capture subtle but consistent signals across diseases.
"""

import os
import json
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from collections import defaultdict

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results/deg")
META_DIR = os.path.join(BASE_DIR, "results/meta_analysis")
os.makedirs(META_DIR, exist_ok=True)


def load_all_deg_results():
    """Load all gene-level DEG results."""
    datasets = {
        "GSE25101": {"disease": "AS", "n_disease": 16, "n_control": 16},
        "GSE73754": {"disease": "AS", "n_disease": 52, "n_control": 20},
        "GSE18781": {"disease": "AS", "n_disease": 30, "n_control": 25},
        "GSE61281": {"disease": "PsA", "n_disease": 20, "n_control": 12},
        "GSE59071": {"disease": "IBD", "n_disease": 105, "n_control": 11},
        "GSE58667": {"disease": "jSpA", "n_disease": 11, "n_control": 4},
        "GSE221786": {"disease": "AS", "n_disease": 20, "n_control": 8},
    }
    
    all_data = {}
    for acc, info in datasets.items():
        # Try gene-level file first, then raw DEG file
        gene_path = os.path.join(RESULTS_DIR, f"{acc}_deg_genes.csv")
        raw_path = os.path.join(RESULTS_DIR, f"{acc}_deg_results.csv")
        
        filepath = gene_path if os.path.exists(gene_path) else raw_path
        if not os.path.exists(filepath):
            print(f"  WARNING: No DEG results for {acc}")
            continue
        
        df = pd.read_csv(filepath)
        
        # Use gene_symbol if available, otherwise gene_id
        id_col = "gene_symbol" if "gene_symbol" in df.columns else "gene_id"
        
        # Filter to valid entries
        df = df[df["pvalue"].notna() & df[id_col].notna()].copy()
        df["gene"] = df[id_col].astype(str)
        
        # Remove entries that look like probe IDs (not gene symbols)
        # Keep only entries that look like gene symbols (alphabetic start, no underscores with numbers)
        # This is a heuristic — gene symbols typically start with letters
        mask = df["gene"].str.match(r'^[A-Z][A-Z0-9]+$', na=False)
        gene_like = df[mask]
        
        if len(gene_like) < 100:
            # If too few match gene symbol pattern, keep all (probe IDs)
            print(f"  {acc}: Using all {len(df)} entries (probe-level)")
            all_data[acc] = {"df": df, **info}
        else:
            print(f"  {acc}: {len(gene_like)} gene-level entries (from {len(df)} total)")
            all_data[acc] = {"df": gene_like, **info}
    
    return all_data


def fisher_meta_analysis(pvalues):
    """Fisher's method for combining p-values."""
    pvalues = [p for p in pvalues if not np.isnan(p) and p > 0]
    if len(pvalues) < 2:
        return np.nan, len(pvalues)
    
    # Clip very small p-values to avoid log(0)
    pvalues = [max(p, 1e-300) for p in pvalues]
    
    chi2_stat = -2 * sum(np.log(p) for p in pvalues)
    df = 2 * len(pvalues)
    combined_p = 1 - stats.chi2.cdf(chi2_stat, df)
    
    return combined_p, len(pvalues)


def weighted_effect_size(log2fcs, sample_sizes):
    """Weighted mean effect size (inverse variance weighting approximation)."""
    valid = [(fc, n) for fc, n in zip(log2fcs, sample_sizes) if not np.isnan(fc)]
    if not valid:
        return np.nan, 0
    
    fcs, ns = zip(*valid)
    weights = np.array(ns, dtype=float)
    weights = weights / weights.sum()
    
    weighted_fc = sum(w * fc for w, fc in zip(weights, fcs))
    return weighted_fc, len(valid)


def run_meta_analysis(all_data, subset_name="all", subset_diseases=None):
    """Run meta-analysis across specified datasets."""
    
    print(f"\n{'='*60}")
    print(f"Meta-analysis: {subset_name}")
    print(f"{'='*60}")
    
    # Filter to requested diseases
    if subset_diseases:
        data = {k: v for k, v in all_data.items() if v["disease"] in subset_diseases}
    else:
        data = all_data
    
    print(f"  Datasets: {list(data.keys())}")
    print(f"  Diseases: {set(v['disease'] for v in data.values())}")
    
    # Collect all genes across datasets
    gene_results = defaultdict(lambda: {
        "pvalues": [], "log2fcs": [], "sample_sizes": [], 
        "datasets": [], "diseases": [], "directions": []
    })
    
    for acc, info in data.items():
        df = info["df"]
        n_total = info["n_disease"] + info["n_control"]
        
        for _, row in df.iterrows():
            gene = row["gene"]
            p = row.get("pvalue", np.nan)
            fc = row.get("log2FC", np.nan)
            
            gene_results[gene]["pvalues"].append(p)
            gene_results[gene]["log2fcs"].append(fc)
            gene_results[gene]["sample_sizes"].append(n_total)
            gene_results[gene]["datasets"].append(acc)
            gene_results[gene]["diseases"].append(info["disease"])
            
            if not np.isnan(fc):
                gene_results[gene]["directions"].append("UP" if fc > 0 else "DOWN")
    
    # Run meta-analysis for each gene
    meta_results = []
    
    for gene, info in gene_results.items():
        n_datasets = len(info["pvalues"])
        
        # Fisher's combined p-value
        fisher_p, n_fisher = fisher_meta_analysis(info["pvalues"])
        
        # Weighted effect size
        weighted_fc, n_fc = weighted_effect_size(info["log2fcs"], info["sample_sizes"])
        
        # Vote counting: consistency of direction
        n_up = info["directions"].count("UP")
        n_down = info["directions"].count("DOWN")
        n_total_votes = n_up + n_down
        
        if n_total_votes > 0:
            consistency = max(n_up, n_down) / n_total_votes
            consensus_direction = "UP" if n_up > n_down else "DOWN"
        else:
            consistency = 0
            consensus_direction = "NS"
        
        # Number of unique diseases
        n_diseases = len(set(info["diseases"]))
        
        meta_results.append({
            "gene": gene,
            "fisher_p": fisher_p,
            "n_datasets": n_datasets,
            "n_fisher": n_fisher,
            "weighted_log2FC": weighted_fc,
            "n_up": n_up,
            "n_down": n_down,
            "consistency": consistency,
            "consensus_direction": consensus_direction,
            "n_diseases": n_diseases,
            "diseases": ",".join(sorted(set(info["diseases"]))),
            "datasets": ",".join(info["datasets"]),
        })
    
    meta_df = pd.DataFrame(meta_results)
    
    # BH correction on Fisher's combined p-values
    valid_mask = meta_df["fisher_p"].notna()
    if valid_mask.sum() > 0:
        _, padj, _, _ = multipletests(
            meta_df.loc[valid_mask, "fisher_p"].values,
            method="fdr_bh"
        )
        meta_df.loc[valid_mask, "fisher_padj"] = padj
    
    # Define significance: combined padj < 0.05, present in ≥2 datasets, consistency ≥ 0.7
    meta_df["meta_significant"] = (
        (meta_df["fisher_padj"] < 0.05) &
        (meta_df["n_datasets"] >= 2) &
        (meta_df["consistency"] >= 0.7)
    )
    
    # Define "core signature": significant + present across ≥2 disease categories
    meta_df["core_signature"] = (
        meta_df["meta_significant"] &
        (meta_df["n_diseases"] >= 2)
    )
    
    # Sort by significance
    meta_df = meta_df.sort_values("fisher_p")
    
    # Print summary
    n_genes_total = len(meta_df)
    n_sig = meta_df["meta_significant"].sum()
    n_core = meta_df["core_signature"].sum()
    
    print(f"\n  Total genes analyzed: {n_genes_total}")
    print(f"  Meta-significant (padj<0.05, ≥2 datasets, consistency≥0.7): {n_sig}")
    print(f"  Core signature (+ ≥2 diseases): {n_core}")
    
    if n_sig > 0:
        sig_df = meta_df[meta_df["meta_significant"]].head(30)
        print(f"\n  Top 30 meta-significant genes:")
        print(f"  {'Gene':15s} {'Fisher padj':>12s} {'wLog2FC':>8s} {'Datasets':>9s} {'Diseases':>9s} {'Direction':>10s} {'Core':>5s}")
        print("  " + "-" * 75)
        for _, row in sig_df.iterrows():
            core = "✓" if row["core_signature"] else ""
            print(f"  {str(row['gene']):15s} {row['fisher_padj']:.2e}  {row['weighted_log2FC']:+.3f}  "
                  f"{row['n_datasets']:>4d}       {row['n_diseases']:>4d}      {row['consensus_direction']:>5s}  {core:>5s}")
    
    # Save results
    output_path = os.path.join(META_DIR, f"meta_analysis_{subset_name}.csv")
    meta_df.to_csv(output_path, index=False)
    
    # Save significant only
    if n_sig > 0:
        sig_path = os.path.join(META_DIR, f"meta_significant_{subset_name}.csv")
        meta_df[meta_df["meta_significant"]].to_csv(sig_path, index=False)
    
    if n_core > 0:
        core_path = os.path.join(META_DIR, f"core_signature_{subset_name}.csv")
        meta_df[meta_df["core_signature"]].to_csv(core_path, index=False)
    
    return meta_df


def main():
    print("=" * 80)
    print("Phase 4: Cross-Disease Transcriptomic Meta-Analysis")
    print("=" * 80)
    
    # Load all DEG results
    all_data = load_all_deg_results()
    
    if not all_data:
        print("ERROR: No DEG results available")
        return
    
    # 1. Full cross-disease meta-analysis (all datasets)
    full_meta = run_meta_analysis(all_data, "cross_disease")
    
    # 2. AS-only meta-analysis
    as_meta = run_meta_analysis(all_data, "AS_only", {"AS"})
    
    # 3. SpA-specific (excluding IBD)
    spa_meta = run_meta_analysis(all_data, "SpA_spectrum", {"AS", "PsA", "jSpA"})
    
    # Summary comparison
    print(f"\n{'='*80}")
    print("META-ANALYSIS COMPARISON")
    print(f"{'='*80}")
    
    for name, meta_df in [("Cross-disease", full_meta), ("AS-only", as_meta), ("SpA spectrum", spa_meta)]:
        n_sig = meta_df["meta_significant"].sum() if "meta_significant" in meta_df.columns else 0
        n_core = meta_df.get("core_signature", pd.Series(dtype=bool)).sum()
        print(f"  {name}: {n_sig} significant, {n_core} core signature")
    
    # Find genes significant across multiple analyses
    if full_meta["meta_significant"].sum() > 0:
        cross_sig = set(full_meta[full_meta["meta_significant"]]["gene"])
        
        print(f"\n  Cross-disease significant genes: {len(cross_sig)}")
        
        if as_meta["meta_significant"].sum() > 0:
            as_sig = set(as_meta[as_meta["meta_significant"]]["gene"])
            overlap = cross_sig & as_sig
            print(f"  Overlap with AS-only: {len(overlap)}")
    
    # Save comprehensive summary
    summary = {
        "cross_disease": {
            "n_genes": int(len(full_meta)),
            "n_significant": int(full_meta["meta_significant"].sum()),
            "n_core": int(full_meta.get("core_signature", pd.Series(dtype=bool)).sum()),
        },
        "AS_only": {
            "n_genes": int(len(as_meta)),
            "n_significant": int(as_meta["meta_significant"].sum()),
        },
        "SpA_spectrum": {
            "n_genes": int(len(spa_meta)),
            "n_significant": int(spa_meta["meta_significant"].sum()),
        },
    }
    
    with open(os.path.join(META_DIR, "phase4_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
