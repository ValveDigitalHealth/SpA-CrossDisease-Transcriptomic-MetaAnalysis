#!/usr/bin/env python3
"""
Phase 5: Functional Enrichment Analysis and Publication Figures
Uses g:Profiler API for GO/KEGG enrichment
Creates volcano plots, enrichment dot plots, and overlap UpSet-style charts
"""

import os
import json
import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from collections import defaultdict

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
META_DIR = os.path.join(BASE_DIR, "results/meta_analysis")
ENRICH_DIR = os.path.join(BASE_DIR, "results/enrichment")
FIG_DIR = os.path.join(BASE_DIR, "figures")
DEG_DIR = os.path.join(BASE_DIR, "results/deg")
os.makedirs(ENRICH_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# Publication-quality settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


def run_gprofiler_enrichment(gene_list, organism="hsapiens", sources=["GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC"]):
    """Run enrichment analysis using g:Profiler API."""
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    
    payload = {
        "organism": organism,
        "query": gene_list,
        "sources": sources,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
    }
    
    try:
        resp = requests.post(url, json=payload, timeout=60)
        if resp.status_code == 200:
            data = resp.json()
            results = data.get("result", [])
            
            enrichment = []
            for r in results:
                enrichment.append({
                    "source": r.get("source", ""),
                    "term_id": r.get("native", ""),
                    "term_name": r.get("name", ""),
                    "p_value": r.get("p_value", 1),
                    "term_size": r.get("term_size", 0),
                    "intersection_size": r.get("intersection_size", 0),
                    "query_size": r.get("query_size", 0),
                    "precision": r.get("precision", 0),
                    "recall": r.get("recall", 0),
                })
            
            return pd.DataFrame(enrichment)
    except Exception as e:
        print(f"  g:Profiler error: {e}")
    
    return pd.DataFrame()


def create_volcano_plot(deg_path, accession, disease, ax=None):
    """Create a volcano plot for a single dataset."""
    df = pd.read_csv(deg_path)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    
    # Color coding
    df["color"] = "grey"
    if "padj" in df.columns:
        sig_up = (df["padj"] < 0.05) & (df["log2FC"] > 1)
        sig_down = (df["padj"] < 0.05) & (df["log2FC"] < -1)
        df.loc[sig_up, "color"] = "#e74c3c"  # red
        df.loc[sig_down, "color"] = "#3498db"  # blue
    
    # Plot
    for color, group in df.groupby("color"):
        neg_log10p = -np.log10(group["pvalue"].clip(lower=1e-300))
        ax.scatter(group["log2FC"], neg_log10p, c=color, alpha=0.3, s=3, rasterized=True)
    
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{accession} ({disease})")
    ax.axhline(y=-np.log10(0.05), color='grey', linestyle='--', alpha=0.5)
    ax.axvline(x=1, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(x=-1, color='grey', linestyle='--', alpha=0.5)
    
    # Count significant
    n_up = (df.get("color") == "#e74c3c").sum()
    n_down = (df.get("color") == "#3498db").sum()
    ax.text(0.98, 0.98, f"Up: {n_up}\nDown: {n_down}", 
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    return ax


def create_multi_volcano(datasets):
    """Create multi-panel volcano plot."""
    n = len(datasets)
    cols = min(3, n)
    rows = (n + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(5*cols, 4.5*rows))
    if n == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for i, (acc, disease) in enumerate(datasets):
        deg_path = os.path.join(DEG_DIR, f"{acc}_deg_results.csv")
        if os.path.exists(deg_path):
            create_volcano_plot(deg_path, acc, disease, axes[i])
    
    # Hide empty axes
    for i in range(len(datasets), len(axes)):
        axes[i].set_visible(False)
    
    fig.suptitle("Volcano Plots: Differential Expression Across SpA Spectrum", 
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_path = os.path.join(FIG_DIR, "fig1_volcano_plots.png")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def create_enrichment_dotplot(enrich_df, title, output_name, top_n=20):
    """Create publication-quality enrichment dot plot."""
    if enrich_df.empty:
        print(f"  No enrichment results for {title}")
        return
    
    # Top terms by p-value
    plot_df = enrich_df.nsmallest(top_n, "p_value").copy()
    plot_df = plot_df.sort_values("p_value", ascending=False)  # Reverse for plotting
    
    fig, ax = plt.subplots(figsize=(10, max(4, len(plot_df) * 0.35)))
    
    # Color by source
    source_colors = {
        "GO:BP": "#2ecc71",
        "GO:CC": "#e74c3c",
        "GO:MF": "#3498db",
        "KEGG": "#f39c12",
        "REAC": "#9b59b6",
    }
    
    colors = [source_colors.get(s, "grey") for s in plot_df["source"]]
    sizes = plot_df["intersection_size"] * 10
    
    neg_log10p = -np.log10(plot_df["p_value"].clip(lower=1e-300))
    
    # Truncate long term names
    labels = [name[:60] + "..." if len(name) > 60 else name for name in plot_df["term_name"]]
    
    ax.scatter(neg_log10p, range(len(plot_df)), c=colors, s=sizes, alpha=0.7, edgecolors='grey', linewidth=0.5)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(labels)
    ax.set_xlabel("-log10(adjusted p-value)")
    ax.set_title(title)
    
    # Legend for sources
    legend_elements = [mpatches.Patch(facecolor=c, label=s) for s, c in source_colors.items()]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)
    
    plt.tight_layout()
    output_path = os.path.join(FIG_DIR, f"{output_name}.png")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def create_meta_analysis_summary_figure():
    """Create the key meta-analysis summary figure."""
    # Load AS meta-analysis results
    meta_path = os.path.join(META_DIR, "meta_analysis_AS_only.csv")
    if not os.path.exists(meta_path):
        print("  No AS meta-analysis results")
        return
    
    meta_df = pd.read_csv(meta_path)
    sig_df = meta_df[meta_df["meta_significant"]].copy()
    
    if len(sig_df) == 0:
        print("  No significant meta-analysis genes")
        return
    
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig, hspace=0.35, wspace=0.3)
    
    # Panel A: Meta-analysis forest-style plot (top genes)
    ax1 = fig.add_subplot(gs[0, 0])
    top_genes = sig_df.nsmallest(20, "fisher_padj")
    y_pos = range(len(top_genes))
    
    colors = ["#e74c3c" if d == "UP" else "#3498db" for d in top_genes["consensus_direction"]]
    
    ax1.barh(y_pos, top_genes["weighted_log2FC"], color=colors, alpha=0.7, height=0.6)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(top_genes["gene"], fontsize=8)
    ax1.set_xlabel("Weighted log2 Fold Change")
    ax1.set_title("A. Top AS Meta-Analysis Genes", fontweight='bold')
    ax1.axvline(x=0, color='black', linewidth=0.5)
    ax1.invert_yaxis()
    
    # Panel B: Number of datasets supporting each gene
    ax2 = fig.add_subplot(gs[0, 1])
    dataset_counts = sig_df["n_datasets"].value_counts().sort_index()
    ax2.bar(dataset_counts.index, dataset_counts.values, color="#2ecc71", alpha=0.7)
    ax2.set_xlabel("Number of Supporting Datasets")
    ax2.set_ylabel("Number of Genes")
    ax2.set_title("B. Dataset Support", fontweight='bold')
    
    # Panel C: Direction consistency distribution
    ax3 = fig.add_subplot(gs[1, 0])
    consistency_vals = sig_df["consistency"]
    ax3.hist(consistency_vals, bins=20, color="#9b59b6", alpha=0.7, edgecolor='white')
    ax3.set_xlabel("Direction Consistency")
    ax3.set_ylabel("Number of Genes")
    ax3.set_title("C. Direction Consistency", fontweight='bold')
    ax3.axvline(x=0.7, color='red', linestyle='--', label='Threshold (0.7)')
    ax3.legend()
    
    # Panel D: UP vs DOWN breakdown
    ax4 = fig.add_subplot(gs[1, 1])
    n_up = (sig_df["consensus_direction"] == "UP").sum()
    n_down = (sig_df["consensus_direction"] == "DOWN").sum()
    
    ax4.bar(["Upregulated", "Downregulated"], [n_up, n_down], 
            color=["#e74c3c", "#3498db"], alpha=0.7)
    ax4.set_ylabel("Number of Genes")
    ax4.set_title("D. Direction of Change", fontweight='bold')
    
    for i, v in enumerate([n_up, n_down]):
        ax4.text(i, v + 0.5, str(v), ha='center', fontweight='bold')
    
    fig.suptitle("Cross-Disease Transcriptomic Meta-Analysis of Spondyloarthritis", 
                 fontsize=15, fontweight='bold', y=1.02)
    
    output_path = os.path.join(FIG_DIR, "fig2_meta_analysis_summary.png")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def create_deg_comparison_figure():
    """Create a comparison of DEG counts across datasets."""
    # Load all summaries
    summary_path = os.path.join(DEG_DIR, "phase3_summary.json")
    with open(summary_path) as f:
        summaries = json.load(f)
    
    fig, ax = plt.subplots(figsize=(10, 5))
    
    accessions = [s["accession"] for s in summaries]
    diseases = [s["disease"] for s in summaries]
    n_up = [s["n_up"] for s in summaries]
    n_down = [-s["n_down"] for s in summaries]  # Negative for downregulated
    
    disease_colors = {
        "AS": "#e74c3c",
        "AS/axSpA": "#c0392b",
        "PsA": "#3498db",
        "IBD": "#2ecc71",
        "jSpA": "#f39c12",
    }
    
    colors = [disease_colors.get(d, "grey") for d in diseases]
    
    x = range(len(accessions))
    bars_up = ax.bar(x, n_up, color=colors, alpha=0.8, label="Upregulated")
    bars_down = ax.bar(x, n_down, color=colors, alpha=0.5)
    
    ax.set_xticks(x)
    labels = [f"{acc}\n({dis})" for acc, dis in zip(accessions, diseases)]
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Number of DEGs")
    ax.set_title("Differentially Expressed Genes Across SpA Spectrum Datasets\n(|log2FC| > 1, padj < 0.05)", fontweight='bold')
    ax.axhline(y=0, color='black', linewidth=0.5)
    
    # Add count labels
    for bar, count in zip(bars_up, n_up):
        if count > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5, 
                    str(count), ha='center', va='bottom', fontsize=8)
    
    # Legend
    legend_patches = [mpatches.Patch(color=c, label=d) for d, c in disease_colors.items()]
    ax.legend(handles=legend_patches, loc='upper right')
    
    plt.tight_layout()
    output_path = os.path.join(FIG_DIR, "fig3_deg_comparison.png")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def main():
    print("=" * 80)
    print("Phase 5: Functional Enrichment & Publication Figures")
    print("=" * 80)
    
    # 1. Create volcano plots
    print("\n--- Creating Volcano Plots ---")
    datasets = [
        ("GSE25101", "AS"),
        ("GSE73754", "AS"),
        ("GSE18781", "AS/axSpA"),
        ("GSE61281", "PsA"),
        ("GSE59071", "IBD"),
        ("GSE221786", "AS (RNA-seq)"),
    ]
    create_multi_volcano(datasets)
    
    # 2. Create DEG comparison figure
    print("\n--- Creating DEG Comparison Figure ---")
    create_deg_comparison_figure()
    
    # 3. Run enrichment on AS meta-analysis genes
    print("\n--- Running Functional Enrichment (AS meta-genes) ---")
    meta_path = os.path.join(META_DIR, "meta_analysis_AS_only.csv")
    if os.path.exists(meta_path):
        meta_df = pd.read_csv(meta_path)
        sig_genes = meta_df[meta_df["meta_significant"]]["gene"].tolist()
        
        print(f"  Submitting {len(sig_genes)} AS meta-significant genes to g:Profiler...")
        
        if sig_genes:
            enrich_df = run_gprofiler_enrichment(sig_genes)
            
            if not enrich_df.empty:
                print(f"  Got {len(enrich_df)} enrichment terms")
                enrich_df.to_csv(os.path.join(ENRICH_DIR, "AS_meta_enrichment.csv"), index=False)
                
                create_enrichment_dotplot(enrich_df, 
                    "Functional Enrichment: AS Meta-Analysis Genes",
                    "fig4_AS_enrichment")
                
                # Show top terms
                print(f"\n  Top 15 enrichment terms:")
                for _, row in enrich_df.nsmallest(15, "p_value").iterrows():
                    print(f"    {row['source']:8s} {row['term_name'][:50]:50s} p={row['p_value']:.2e} ({row['intersection_size']} genes)")
            else:
                print("  No enrichment results returned")
    
    # 4. Run enrichment on IBD DEGs (GSE59071 has the most DEGs)
    print("\n--- Running Functional Enrichment (IBD top DEGs) ---")
    ibd_deg_path = os.path.join(DEG_DIR, "GSE59071_deg_results.csv")
    if os.path.exists(ibd_deg_path):
        ibd_df = pd.read_csv(ibd_deg_path)
        ibd_sig = ibd_df[ibd_df["significant"]].nsmallest(200, "padj")
        # These are probe IDs — won't work directly with g:Profiler
        # Skip for now — needs platform annotation
        print("  IBD enrichment skipped (probe IDs need annotation)")
    
    # 5. Meta-analysis summary figure
    print("\n--- Creating Meta-Analysis Summary Figure ---")
    create_meta_analysis_summary_figure()
    
    print(f"\n{'='*80}")
    print("Phase 5 Complete")
    print(f"{'='*80}")
    print(f"\nFigures saved to: {FIG_DIR}")
    print(f"Enrichment results saved to: {ENRICH_DIR}")


if __name__ == "__main__":
    main()
