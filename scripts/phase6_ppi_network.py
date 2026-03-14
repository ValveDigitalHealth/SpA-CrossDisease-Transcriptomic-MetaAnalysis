#!/usr/bin/env python3
"""
Phase 6: PPI Network Construction and Hub Gene Identification.
Uses STRING API (v12.0) for protein-protein interactions.
Builds NetworkX graph and identifies hub genes by centrality.
"""

import os
import json
import requests
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results")
META_DIR = os.path.join(RESULTS_DIR, "meta_analysis")
NETWORK_DIR = os.path.join(RESULTS_DIR, "networks")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(NETWORK_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

STRING_API = "https://string-db.org/api"
SPECIES = 9606
CONFIDENCE_THRESHOLD = 0.7


def get_string_interactions(gene_list, confidence=0.7):
    """Query STRING API for protein-protein interactions."""
    url = f"{STRING_API}/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": SPECIES,
        "required_score": int(confidence * 1000),
        "caller_identity": "SpA_MetaAnalysis_v1",
    }
    try:
        resp = requests.post(url, data=params, timeout=60)
        resp.raise_for_status()
        interactions = resp.json()
        print(f"  STRING returned {len(interactions)} interactions")
        return interactions
    except Exception as e:
        print(f"  STRING API error: {e}")
        return []


def build_network(interactions, gene_list):
    """Build NetworkX graph from STRING interactions."""
    G = nx.Graph()
    G.add_nodes_from(gene_list)

    for interaction in interactions:
        gene_a = interaction.get("preferredName_A", "")
        gene_b = interaction.get("preferredName_B", "")
        score = interaction.get("score", 0)

        if gene_a and gene_b and score >= CONFIDENCE_THRESHOLD:
            G.add_edge(gene_a, gene_b, weight=score)

    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    return G


def compute_centrality(G):
    """Compute network centrality metrics."""
    centrality = {}
    centrality["degree"] = dict(nx.degree_centrality(G))
    centrality["betweenness"] = dict(nx.betweenness_centrality(G))
    centrality["clustering"] = dict(nx.clustering(G))

    # Combine into DataFrame
    nodes = list(G.nodes())
    df = pd.DataFrame({
        "gene": nodes,
        "degree": [G.degree(n) for n in nodes],
        "degree_centrality": [centrality["degree"].get(n, 0) for n in nodes],
        "betweenness_centrality": [centrality["betweenness"].get(n, 0) for n in nodes],
        "clustering_coefficient": [centrality["clustering"].get(n, 0) for n in nodes],
    })
    df = df.sort_values("degree_centrality", ascending=False)
    return df


def visualize_network(G, hub_genes, output_path):
    """Create network visualization."""
    if len(G.nodes()) == 0:
        return

    fig, ax = plt.subplots(1, 1, figsize=(14, 12))

    pos = nx.spring_layout(G, k=2.0, seed=42)

    # Color nodes by hub status
    hub_set = set(hub_genes)
    node_colors = ["#e74c3c" if n in hub_set else "#3498db" for n in G.nodes()]
    node_sizes = [500 if n in hub_set else 200 for n in G.nodes()]

    nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color="gray", ax=ax)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, ax=ax, alpha=0.8)
    nx.draw_networkx_labels(G, pos, font_size=7, ax=ax)

    ax.set_title("PPI Network - Meta-Significant Genes (STRING)", fontsize=14)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    print("Phase 6: PPI Network Analysis")
    print("=" * 60)

    # Load meta-significant genes
    sig_path = os.path.join(META_DIR, "meta_significant_genes.csv")
    if not os.path.exists(sig_path):
        print("ERROR: No meta-significant genes. Run Phase 4 first.")
        return

    sig_df = pd.read_csv(sig_path)
    gene_list = sig_df["gene_symbol"].tolist()
    print(f"Genes for PPI network: {len(gene_list)}")

    # Query STRING
    print("\nQuerying STRING API...")
    interactions = get_string_interactions(gene_list, CONFIDENCE_THRESHOLD)

    if not interactions:
        print("No interactions found. Creating empty network.")
        G = nx.Graph()
        G.add_nodes_from(gene_list)
    else:
        G = build_network(interactions, gene_list)

    # Compute centrality
    print("\nComputing network centrality...")
    centrality_df = compute_centrality(G)

    hub_path = os.path.join(NETWORK_DIR, "hub_genes.csv")
    centrality_df.to_csv(hub_path, index=False)
    print(f"  Hub genes saved: {hub_path}")

    # Top hub genes (top 10% by degree)
    n_hubs = max(5, len(centrality_df) // 10)
    hub_genes = centrality_df.head(n_hubs)["gene"].tolist()
    print(f"\nTop hub genes:")
    for _, row in centrality_df.head(10).iterrows():
        print(f"  {row['gene']:15s} degree={row['degree']:3d}  centrality={row['degree_centrality']:.3f}")

    # Save network
    graphml_path = os.path.join(NETWORK_DIR, "ppi_network.graphml")
    nx.write_graphml(G, graphml_path)
    print(f"  Network saved: {graphml_path}")

    # Network stats
    stats = {
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
        "density": nx.density(G),
        "n_components": nx.number_connected_components(G),
        "top_hubs": hub_genes[:10],
    }
    stats_path = os.path.join(NETWORK_DIR, "network_stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)

    # Visualize
    print("\nCreating network visualization...")
    viz_path = os.path.join(FIGURES_DIR, "ppi_network.png")
    visualize_network(G, hub_genes, viz_path)

    print("\nPhase 6 complete.")


if __name__ == "__main__":
    main()
