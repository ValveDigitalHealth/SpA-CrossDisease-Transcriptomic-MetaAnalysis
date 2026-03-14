#!/usr/bin/env python3
"""
Phase 3b: Map probe IDs to gene symbols for all datasets.
Downloads platform annotation tables from GEO.
"""

import os
import gzip
import re
import json
import requests
import pandas as pd
import time

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
RESULTS_DIR = os.path.join(BASE_DIR, "results/deg")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")

# Platform to dataset mapping
PLATFORM_MAP = {
    "GPL6947": {"datasets": ["GSE25101"], "name": "Illumina HumanHT-12 V3.0"},
    "GPL6883": {"datasets": ["GSE73754"], "name": "Illumina HumanRef-8 v3.0"},
    "GPL570": {"datasets": ["GSE18781"], "name": "Affymetrix HG-U133 Plus 2.0"},
    "GPL6244": {"datasets": ["GSE61281"], "name": "Affymetrix HuGene 1.0 ST"},
    "GPL13667": {"datasets": ["GSE59071"], "name": "Illumina HumanHT-12 V4.0"},
    "GPL96": {"datasets": ["GSE58667"], "name": "Affymetrix HG-U133A"},
}


def download_platform_annotation(gpl_id):
    """Download platform annotation from GEO."""
    output_path = os.path.join(RAW_DIR, f"{gpl_id}_annot.txt.gz")
    
    if os.path.exists(output_path):
        print(f"  Already exists: {output_path}")
        return output_path
    
    # Try SOFT format
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:-3]}nnn/{gpl_id}/annot/{gpl_id}.annot.gz"
    print(f"  Downloading: {url}")
    
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(resp.content)
            print(f"  Downloaded: {os.path.getsize(output_path)/1024:.1f} KB")
            return output_path
    except Exception as e:
        print(f"  Error: {e}")
    
    # Try table format
    url2 = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:-3]}nnn/{gpl_id}/soft/{gpl_id}_family.soft.gz"
    output_path2 = os.path.join(RAW_DIR, f"{gpl_id}_soft.gz")
    
    print(f"  Trying SOFT format: {url2}")
    try:
        resp = requests.get(url2, timeout=120)
        if resp.status_code == 200:
            with open(output_path2, "wb") as f:
                f.write(resp.content)
            print(f"  Downloaded SOFT: {os.path.getsize(output_path2)/1024:.1f} KB")
            return output_path2
    except Exception as e:
        print(f"  Error: {e}")
    
    return None


def parse_platform_annotation(filepath, gpl_id):
    """Parse platform annotation to extract probe → gene symbol mapping."""
    if filepath is None:
        return {}
    
    probe_to_gene = {}
    
    try:
        opener = gzip.open if filepath.endswith(".gz") else open
        
        with opener(filepath, "rt", errors="replace") as f:
            # Skip header lines starting with ! or #
            in_table = False
            header = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith("!platform_table_begin") or line.startswith("#"):
                    in_table = True
                    continue
                elif line.startswith("!platform_table_end"):
                    break
                
                if not in_table and not line.startswith("!") and not line.startswith("#"):
                    # This might be the annotation file format (tab-separated with header)
                    if header is None and ("ID" in line or "Gene Symbol" in line or "Symbol" in line):
                        header = line.split("\t")
                        continue
                    elif header:
                        fields = line.split("\t")
                        if len(fields) >= len(header):
                            # Find ID and gene symbol columns
                            id_idx = None
                            gene_idx = None
                            for i, h in enumerate(header):
                                h_lower = h.lower().strip()
                                if h_lower == "id":
                                    id_idx = i
                                elif h_lower in ["gene symbol", "symbol", "gene_symbol", "gene_assignment"]:
                                    gene_idx = i
                            
                            if id_idx is not None and gene_idx is not None:
                                probe_id = fields[id_idx].strip()
                                gene = fields[gene_idx].strip()
                                if gene and gene != "---" and gene != "":
                                    # Handle multiple gene symbols (take first)
                                    gene = gene.split("///")[0].split("//")[0].strip()
                                    probe_to_gene[probe_id] = gene
                        continue
                
                if in_table:
                    if header is None:
                        header = line.split("\t")
                    else:
                        fields = line.split("\t")
                        if len(fields) >= 2:
                            id_idx = None
                            gene_idx = None
                            for i, h in enumerate(header):
                                h_lower = h.lower().strip()
                                if h_lower == "id":
                                    id_idx = i
                                elif h_lower in ["gene symbol", "symbol", "gene_symbol", "gene_assignment", "gene"]:
                                    gene_idx = i
                            
                            if id_idx is not None and gene_idx is not None and gene_idx < len(fields):
                                probe_id = fields[id_idx].strip()
                                gene = fields[gene_idx].strip()
                                if gene and gene != "---":
                                    gene = gene.split("///")[0].split("//")[0].strip()
                                    probe_to_gene[probe_id] = gene
    
    except Exception as e:
        print(f"  Error parsing {filepath}: {e}")
    
    print(f"  Mapped {len(probe_to_gene)} probes to gene symbols for {gpl_id}")
    return probe_to_gene


def map_deg_to_genes(accession, probe_to_gene):
    """Map DEG results from probe IDs to gene symbols."""
    deg_path = os.path.join(RESULTS_DIR, f"{accession}_deg_results.csv")
    
    if not os.path.exists(deg_path):
        print(f"  No DEG results for {accession}")
        return None
    
    df = pd.read_csv(deg_path)
    
    # Map probe IDs to gene symbols
    df["gene_symbol"] = df["gene_id"].astype(str).map(probe_to_gene)
    
    # Count mappings
    n_mapped = df["gene_symbol"].notna().sum()
    n_total = len(df)
    print(f"  {accession}: {n_mapped}/{n_total} probes mapped to gene symbols ({n_mapped/n_total*100:.1f}%)")
    
    # For genes with multiple probes, keep the one with smallest p-value
    mapped_df = df[df["gene_symbol"].notna()].copy()
    if len(mapped_df) > 0:
        mapped_df = mapped_df.sort_values("pvalue").drop_duplicates(subset="gene_symbol", keep="first")
        print(f"  After deduplication: {len(mapped_df)} unique genes")
    
    # Save gene-level results
    output_path = os.path.join(RESULTS_DIR, f"{accession}_deg_genes.csv")
    mapped_df.to_csv(output_path, index=False)
    
    return mapped_df


def main():
    print("=" * 80)
    print("Phase 3b: Probe-to-Gene Symbol Mapping")
    print("=" * 80)
    
    all_mappings = {}
    
    for gpl_id, info in PLATFORM_MAP.items():
        print(f"\n--- {gpl_id}: {info['name']} ---")
        
        filepath = download_platform_annotation(gpl_id)
        mapping = parse_platform_annotation(filepath, gpl_id)
        
        if not mapping:
            print(f"  WARNING: No mappings extracted for {gpl_id}")
            # Try alternative: use the expression file's probe IDs as gene symbols
            # (some platforms use gene symbols directly)
        
        all_mappings[gpl_id] = mapping
        
        # Map DEGs for associated datasets
        for acc in info["datasets"]:
            if mapping:
                map_deg_to_genes(acc, mapping)
        
        time.sleep(0.5)
    
    # Handle GSE221786 (RNA-seq with ENSEMBL IDs) - needs ENSEMBL to gene symbol mapping
    print(f"\n--- GSE221786 (RNA-seq, ENSEMBL IDs) ---")
    deg_path = os.path.join(RESULTS_DIR, "GSE221786_deg_results.csv")
    if os.path.exists(deg_path):
        df = pd.read_csv(deg_path)
        # The gene_id has version numbers (ENSG00000006075.11)
        # Strip version for mapping
        df["ensembl_id"] = df["gene_id"].astype(str).str.split(".").str[0]
        
        # We need to map ENSEMBL to gene symbols
        # Use the gene_symbol column from the counts file
        counts_path = os.path.join(BASE_DIR, "data/processed/GSE221786_counts.csv")
        if os.path.exists(counts_path):
            counts_df = pd.read_csv(counts_path, index_col=0, usecols=[0, 1], nrows=30000)
            if "gene_symbol" in counts_df.columns:
                ensembl_to_gene = dict(zip(counts_df.index.astype(str), counts_df["gene_symbol"]))
                df["gene_symbol"] = df["gene_id"].astype(str).map(ensembl_to_gene)
                n_mapped = df["gene_symbol"].notna().sum()
                print(f"  GSE221786: {n_mapped}/{len(df)} mapped via counts file")
                
                mapped_df = df[df["gene_symbol"].notna()].copy()
                mapped_df = mapped_df.sort_values("pvalue").drop_duplicates(subset="gene_symbol", keep="first")
                output_path = os.path.join(RESULTS_DIR, "GSE221786_deg_genes.csv")
                mapped_df.to_csv(output_path, index=False)
                print(f"  Saved: {len(mapped_df)} unique genes")
    
    # Summary
    print(f"\n{'='*80}")
    print("PROBE MAPPING SUMMARY")
    print(f"{'='*80}")
    
    gene_files = [f for f in os.listdir(RESULTS_DIR) if f.endswith("_deg_genes.csv")]
    for f in sorted(gene_files):
        filepath = os.path.join(RESULTS_DIR, f)
        df = pd.read_csv(filepath)
        n_sig = df.get("significant", pd.Series()).sum() if "significant" in df.columns else 0
        print(f"  {f}: {len(df)} genes, {n_sig} significant")
    
    # Save mapping stats
    stats = {}
    for gpl_id, mapping in all_mappings.items():
        stats[gpl_id] = {"n_probes_mapped": len(mapping)}
    
    stats_path = os.path.join(METADATA_DIR, "probe_mapping_stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    main()
