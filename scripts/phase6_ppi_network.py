#!/usr/bin/env python3
"""
Phase 6: PPI Network Construction and Hub Gene Identification
Uses STRING API for protein-protein interactions.
Phase 9: Drug-Gene Interaction via DGIdb API.
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
from collections import defaultdict, Counter
import networkx as nx

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
META_DIR = os.path.join(BASE_DIR, "results/meta_analysis")
NET_DIR = os.path.join(BASE_DIR, "results/networks")
DRUG_DIR = os.path.join(BASE_DIR, "results/drug_repurposing")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(NET_DIR, exist_ok=True)
os.makedirs(DRUG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'figure.dpi': 300,
})


def get_string_interactions(gene_list, species=9606, score_threshold=400):
    """Query STRING API for protein-protein interactions."""
    url = "https://string-db.org/api/json/network"
    
    # STRING API accepts up to 2000 proteins
    params = {
        "identifiers": "%0d".join(gene_list[:200]),
        "species": species,
        "required_score": score_threshold,
        "network_type": "functional",
        "caller_identity": "spacbi_metaanalysis",
    }
    
    try:
        resp = requests.get(url, params=params, timeout=60)
        if resp.status_code == 200:
            interactions = resp.json()
            print(f"  STRING returned {len(interactions)} interactions")
            return interactions
    except Exception as e:
        print(f"  STRING API error: {e}")
    
    return []


def get_string_enrichment(gene_list, species=9606):
    """Get functional enrichment from STRING."""
    url = "https://string-db.org/api/json/enrichment"
    
    params = {
        "identifiers": "%0d".join(gene_list[:200]),
        "species": species,
        "caller_identity": "spacbi_metaanalysis",
    }
    
    try:
        resp = requests.get(url, params=params, timeout=60)
        if resp.status_code == 200:
            enrichment = resp.json()
            return enrichment
    except Exception as e:
        print(f"  STRING enrichment error: {e}")
    
    return []


def build_network(interactions, gene_info=None):
    """Build a NetworkX graph from STRING interactions."""
    G = nx.Graph()
    
    for edge in interactions:
        source = edge.get("preferredName_A", edge.get("stringId_A", ""))
        target = edge.get("preferredName_B", edge.get("stringId_B", ""))
        score = edge.get("score", 0)
        
        G.add_edge(source, target, weight=score)
    
    # Add node attributes if gene_info provided
    if gene_info is not None:
        for gene, info in gene_info.items():
            if gene in G.nodes:
                for key, val in info.items():
                    G.nodes[gene][key] = val
    
    return G


def calculate_hub_genes(G, top_n=20):
    """Calculate hub genes using multiple centrality metrics."""
    if len(G.nodes) == 0:
        return pd.DataFrame()
    
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)
    
    try:
        eigenvector = nx.eigenvector_centrality_numpy(G)
    except:
        eigenvector = {n: 0 for n in G.nodes}
    
    hub_data = []
    for node in G.nodes:
        hub_data.append({
            "gene": node,
            "degree": degree.get(node, 0),
            "betweenness": betweenness.get(node, 0),
            "closeness": closeness.get(node, 0),
            "eigenvector": eigenvector.get(node, 0),
        })
    
    hub_df = pd.DataFrame(hub_data)
    
    # Rank by each metric
    for col in ["degree", "betweenness", "closeness", "eigenvector"]:
        hub_df[f"{col}_rank"] = hub_df[col].rank(ascending=False)
    
    # Composite rank (average of all ranks)
    rank_cols = [c for c in hub_df.columns if c.endswith("_rank")]
    hub_df["composite_rank"] = hub_df[rank_cols].mean(axis=1)
    hub_df = hub_df.sort_values("composite_rank")
    
    return hub_df


def visualize_network(G, hub_df, output_name, title="PPI Network"):
    """Create publication-quality network visualization."""
    if len(G.nodes) == 0:
        print("  Empty network, skipping visualization")
        return
    
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Layout
    if len(G.nodes) > 50:
        pos = nx.spring_layout(G, k=2/np.sqrt(len(G.nodes)), iterations=50, seed=42)
    else:
        pos = nx.spring_layout(G, k=3/np.sqrt(max(len(G.nodes), 1)), iterations=100, seed=42)
    
    # Node sizes based on degree
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees else 1
    node_sizes = [100 + 500 * degrees.get(n, 0) / max_degree for n in G.nodes]
    
    # Node colors based on hub rank
    top_hubs = set(hub_df.nsmallest(10, "composite_rank")["gene"]) if len(hub_df) > 0 else set()
    node_colors = ["#e74c3c" if n in top_hubs else "#95a5a6" for n in G.nodes]
    
    # Draw edges
    nx.draw_networkx_edges(G, pos, alpha=0.15, width=0.5, edge_color='grey', ax=ax)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, 
                          alpha=0.8, edgecolors='white', linewidths=0.5, ax=ax)
    
    # Labels for hub genes only
    label_genes = top_hubs
    labels = {n: n for n in G.nodes if n in label_genes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold', ax=ax)
    
    ax.set_title(f"{title}\n({len(G.nodes)} nodes, {len(G.edges)} edges)", fontweight='bold')
    ax.axis('off')
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#e74c3c', label='Top 10 Hub Genes'),
        mpatches.Patch(facecolor='#95a5a6', label='Other Genes'),
    ]
    ax.legend(handles=legend_elements, loc='lower left')
    
    output_path = os.path.join(FIG_DIR, f"{output_name}.png")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def query_dgidb(gene_list):
    """Query DGIdb for drug-gene interactions."""
    url = "https://dgidb.org/api/v2/interactions.json"
    
    params = {
        "genes": ",".join(gene_list[:100]),
    }
    
    try:
        resp = requests.get(url, params=params, timeout=60)
        if resp.status_code == 200:
            data = resp.json()
            interactions = []
            
            for match in data.get("matchedTerms", []):
                gene = match.get("geneName", "")
                for interaction in match.get("interactions", []):
                    drug = interaction.get("drugName", "")
                    interaction_type = interaction.get("interactionTypes", [])
                    types_str = "; ".join([t.get("type", "") for t in interaction_type]) if interaction_type else ""
                    sources = interaction.get("sources", [])
                    pmids = interaction.get("pmids", [])
                    
                    interactions.append({
                        "gene": gene,
                        "drug": drug,
                        "interaction_type": types_str,
                        "n_sources": len(sources),
                        "sources": "; ".join(sources[:5]),
                        "pmids": "; ".join(pmids[:3]),
                    })
            
            return pd.DataFrame(interactions)
    except Exception as e:
        print(f"  DGIdb error: {e}")
    
    return pd.DataFrame()


def create_hub_gene_figure(hub_df, output_name):
    """Create hub gene ranking figure."""
    if hub_df.empty:
        return
    
    top20 = hub_df.head(20)
    
    fig, axes = plt.subplots(1, 4, figsize=(16, 8), sharey=True)
    
    metrics = ["degree", "betweenness", "closeness", "eigenvector"]
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12"]
    
    for ax, metric, color in zip(axes, metrics, colors):
        ax.barh(range(len(top20)), top20[metric], color=color, alpha=0.7)
        ax.set_xlabel(metric.capitalize())
        ax.set_title(metric.capitalize(), fontweight='bold')
    
    axes[0].set_yticks(range(len(top20)))
    axes[0].set_yticklabels(top20["gene"])
    axes[0].invert_yaxis()
    
    fig.suptitle("Top 20 Hub Genes by Centrality Metrics", fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    output_path = os.path.join(FIG_DIR, f"{output_name}.png")
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close()


def main():
    print("=" * 80)
    print("Phase 6: PPI Network & Hub Genes + Phase 9: Drug Repurposing")
    print("=" * 80)
    
    # Load AS meta-significant genes
    meta_path = os.path.join(META_DIR, "meta_analysis_AS_only.csv")
    meta_df = pd.read_csv(meta_path)
    sig_genes = meta_df[meta_df["meta_significant"]]["gene"].tolist()
    
    print(f"\nAS meta-significant genes: {len(sig_genes)}")
    print(f"  Genes: {sig_genes}")
    
    # Also get the broader set (nominal p < 0.01, consistent direction)
    broader_genes = meta_df[
        (meta_df["fisher_padj"] < 0.1) & 
        (meta_df["n_datasets"] >= 2) &
        (meta_df["consistency"] >= 0.7)
    ]["gene"].tolist()
    
    print(f"\nBroader gene set (padj < 0.1): {len(broader_genes)}")
    
    # ==== PPI Network ====
    print(f"\n--- STRING PPI Network ---")
    
    # Use the broader gene set for network
    query_genes = broader_genes if len(broader_genes) > 10 else sig_genes
    
    interactions = get_string_interactions(query_genes)
    
    if interactions:
        # Build network
        gene_info = {}
        for _, row in meta_df[meta_df["gene"].isin(query_genes)].iterrows():
            gene_info[row["gene"]] = {
                "log2FC": row.get("weighted_log2FC", 0),
                "fisher_padj": row.get("fisher_padj", 1),
                "n_datasets": row.get("n_datasets", 0),
            }
        
        G = build_network(interactions, gene_info)
        print(f"  Network: {len(G.nodes)} nodes, {len(G.edges)} edges")
        
        # Hub genes
        hub_df = calculate_hub_genes(G)
        hub_df.to_csv(os.path.join(NET_DIR, "hub_genes.csv"), index=False)
        
        print(f"\n  Top 10 Hub Genes (composite rank):")
        for _, row in hub_df.head(10).iterrows():
            print(f"    {row['gene']:15s} degree={row['degree']:<4} betweenness={row['betweenness']:.3f} "
                  f"closeness={row['closeness']:.3f}")
        
        # Network statistics
        stats = {
            "n_nodes": len(G.nodes),
            "n_edges": len(G.edges),
            "density": nx.density(G),
            "avg_clustering": nx.average_clustering(G),
            "n_components": nx.number_connected_components(G),
        }
        
        with open(os.path.join(NET_DIR, "network_stats.json"), "w") as f:
            json.dump(stats, f, indent=2)
        
        print(f"\n  Network statistics:")
        for k, v in stats.items():
            print(f"    {k}: {v}")
        
        # Visualize
        visualize_network(G, hub_df, "fig5_ppi_network", "SpA Meta-Analysis PPI Network")
        create_hub_gene_figure(hub_df, "fig6_hub_genes")
        
        # Save edge list
        edges_df = pd.DataFrame([
            {"source": u, "target": v, "weight": d.get("weight", 0)}
            for u, v, d in G.edges(data=True)
        ])
        edges_df.to_csv(os.path.join(NET_DIR, "ppi_edges.csv"), index=False)
        
    else:
        print("  No STRING interactions found")
        hub_df = pd.DataFrame()
    
    # ==== Drug-Gene Interactions ====
    print(f"\n--- Drug-Gene Interactions (DGIdb) ---")
    
    # Query hub genes + significant genes
    drug_query = sig_genes[:50]
    if not hub_df.empty:
        top_hubs = hub_df.head(20)["gene"].tolist()
        drug_query = list(set(drug_query + top_hubs))
    
    print(f"  Querying DGIdb for {len(drug_query)} genes...")
    drug_df = query_dgidb(drug_query)
    
    if not drug_df.empty:
        drug_df.to_csv(os.path.join(DRUG_DIR, "drug_gene_interactions.csv"), index=False)
        
        n_genes_with_drugs = drug_df["gene"].nunique()
        n_drugs = drug_df["drug"].nunique()
        
        print(f"  Found {len(drug_df)} interactions")
        print(f"  {n_genes_with_drugs} genes with drug interactions")
        print(f"  {n_drugs} unique drugs")
        
        # Top drugs by number of target genes
        drug_counts = drug_df.groupby("drug")["gene"].nunique().sort_values(ascending=False)
        print(f"\n  Top drugs by target count:")
        for drug, count in drug_counts.head(15).items():
            targets = drug_df[drug_df["drug"] == drug]["gene"].unique()
            interaction = drug_df[drug_df["drug"] == drug]["interaction_type"].iloc[0]
            print(f"    {drug:30s} targets {count} gene(s): {', '.join(targets[:5])} [{interaction}]")
        
        # Create drug-target figure
        fig, ax = plt.subplots(figsize=(10, 8))
        top_drugs = drug_counts.head(15)
        colors = plt.cm.RdYlBu(np.linspace(0.2, 0.8, len(top_drugs)))
        
        ax.barh(range(len(top_drugs)), top_drugs.values, color=colors)
        ax.set_yticks(range(len(top_drugs)))
        ax.set_yticklabels(top_drugs.index, fontsize=9)
        ax.set_xlabel("Number of Target Genes")
        ax.set_title("Drug-Gene Interactions for SpA Hub Genes (DGIdb)", fontweight='bold')
        ax.invert_yaxis()
        
        plt.tight_layout()
        output_path = os.path.join(FIG_DIR, "fig7_drug_targets.png")
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_path}")
        plt.close()
    else:
        print("  No drug interactions found")
    
    # ==== STRING Enrichment ====
    print(f"\n--- STRING Functional Enrichment ---")
    enrichment = get_string_enrichment(query_genes)
    
    if enrichment:
        enrich_df = pd.DataFrame(enrichment)
        enrich_df.to_csv(os.path.join(NET_DIR, "string_enrichment.csv"), index=False)
        
        # Show top enrichment
        if "category" in enrich_df.columns:
            for cat in enrich_df["category"].unique()[:3]:
                cat_df = enrich_df[enrich_df["category"] == cat].head(5)
                print(f"\n  {cat}:")
                for _, row in cat_df.iterrows():
                    print(f"    {row.get('description', 'N/A')[:60]:60s} p={row.get('p_value', row.get('fdr', 'N/A'))}")
    
    print(f"\n{'='*80}")
    print("Phase 6 + 9 Complete")
    print(f"{'='*80}")


if __name__ == "__main__":
    main()
