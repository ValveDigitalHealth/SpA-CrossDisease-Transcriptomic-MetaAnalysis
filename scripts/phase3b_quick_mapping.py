#!/usr/bin/env python3
"""
Phase 3b (quick): Map probe IDs to gene symbols using g:Profiler API.
More reliable than downloading GEO platform annotations.
"""

import os
import json
import requests
import pandas as pd
import numpy as np
from collections import defaultdict

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results/deg")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")


def map_illumina_probes_from_expression(accession):
    """For Illumina platforms, probe IDs in series matrix are often ILMN_ ids.
    We can use the expression file and the GEO annotation embedded in series matrix."""
    deg_path = os.path.join(RESULTS_DIR, f"{accession}_deg_results.csv")
    df = pd.read_csv(deg_path)
    
    # Check if gene_id looks like gene symbols already
    sample_ids = df["gene_id"].astype(str).head(20).tolist()
    
    # For Illumina, probe IDs start with ILMN_
    has_ilmn = any("ILMN_" in str(x) for x in sample_ids)
    # For Affymetrix, IDs are like 1007_s_at or numeric
    has_affy = any("_at" in str(x) or "_s_at" in str(x) for x in sample_ids)
    has_agilent = any("A_23_P" in str(x) or "A_24_P" in str(x) for x in sample_ids)
    # For HumanHT-12 V4, IDs are numeric
    all_numeric = all(str(x).isdigit() for x in sample_ids)
    
    print(f"  {accession}: ILMN={has_ilmn}, Affy={has_affy}, Agilent={has_agilent}, Numeric={all_numeric}")
    return df


def use_gprofiler_convert(probe_ids, namespace, organism="hsapiens"):
    """Use g:Profiler gconvert API to map IDs to gene symbols."""
    url = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
    
    # Process in batches of 500
    all_mappings = {}
    
    for i in range(0, len(probe_ids), 500):
        batch = probe_ids[i:i+500]
        
        payload = {
            "organism": organism,
            "target": "ENSG",
            "query": batch,
        }
        
        try:
            resp = requests.post(url, json=payload, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                results = data.get("result", [])
                for r in results:
                    incoming = r.get("incoming", "")
                    name = r.get("name", "")
                    if name and name != "N/A":
                        all_mappings[incoming] = name
        except Exception as e:
            print(f"  g:Profiler batch error: {e}")
    
    return all_mappings


def use_biomart_mapping(probe_type="illumina"):
    """Use Ensembl BioMart REST API for probe mapping."""
    # This is more reliable for large-scale mapping
    base_url = "https://rest.ensembl.org"
    
    # We'll use a different approach: extract gene symbols from GEO series matrix
    # metadata which sometimes includes gene annotations
    pass


def extract_gene_symbols_from_series_matrix(accession):
    """Try to extract gene symbol annotations from series matrix metadata."""
    import gzip
    
    matrix_path = os.path.join(BASE_DIR, f"data/raw/{accession}_series_matrix.txt.gz")
    if not os.path.exists(matrix_path):
        return {}
    
    probe_to_gene = {}
    
    with gzip.open(matrix_path, "rt", errors="replace") as f:
        for line in f:
            # Look for gene annotation table (some series matrices include this)
            if line.startswith("!platform_table_begin"):
                # Read platform table
                header = next(f).strip().split("\t")
                
                id_idx = None
                gene_idx = None
                for j, h in enumerate(header):
                    h_lower = h.lower().strip()
                    if h_lower == "id":
                        id_idx = j
                    elif "symbol" in h_lower or "gene" in h_lower:
                        gene_idx = j
                
                if id_idx is not None and gene_idx is not None:
                    for pline in f:
                        if pline.startswith("!platform_table_end"):
                            break
                        fields = pline.strip().split("\t")
                        if len(fields) > max(id_idx, gene_idx):
                            probe = fields[id_idx]
                            gene = fields[gene_idx].split("///")[0].strip()
                            if gene and gene != "---":
                                probe_to_gene[probe] = gene
    
    return probe_to_gene


def create_gene_level_results_with_relaxed_approach(accession, disease):
    """Create gene-level DEG results, handling probe mapping pragmatically."""
    deg_path = os.path.join(RESULTS_DIR, f"{accession}_deg_results.csv")
    if not os.path.exists(deg_path):
        return None
    
    df = pd.read_csv(deg_path)
    
    # Check what kind of IDs we have
    sample_ids = df["gene_id"].astype(str).head(5).tolist()
    print(f"\n  {accession} sample IDs: {sample_ids}")
    
    # Strategy 1: IDs might already be gene symbols
    # Strategy 2: Use ENSEMBL IDs from RNA-seq  
    # Strategy 3: Use the probe IDs as-is for now (we'll aggregate at meta-analysis)
    
    # For RNA-seq with ENSEMBL IDs (GSE221786)
    if any("ENSG" in str(x) for x in sample_ids):
        counts_path = os.path.join(BASE_DIR, "data/processed/GSE221786_counts.csv")
        if os.path.exists(counts_path):
            counts_df = pd.read_csv(counts_path, index_col=0, usecols=[0, 1])
            if "gene_symbol" in counts_df.columns:
                mapping = dict(zip(counts_df.index.astype(str), counts_df["gene_symbol"]))
                df["gene_symbol"] = df["gene_id"].astype(str).map(mapping)
                n_mapped = df["gene_symbol"].notna().sum()
                print(f"    ENSEMBL mapping: {n_mapped}/{len(df)}")
    
    # For Illumina ILMN_ probes, try g:Profiler
    elif any("ILMN_" in str(x) for x in sample_ids):
        probe_ids = df["gene_id"].astype(str).tolist()
        print(f"    Attempting g:Profiler for {len(probe_ids)} ILMN probes...")
        mapping = use_gprofiler_convert(probe_ids[:2000], "ILLUMINA_HUMANHT_12_V3")
        if mapping:
            df["gene_symbol"] = df["gene_id"].astype(str).map(mapping)
            n_mapped = df["gene_symbol"].notna().sum()
            print(f"    g:Profiler mapping: {n_mapped}/{len(df)}")
        else:
            print(f"    g:Profiler failed, using probe IDs as-is")
            df["gene_symbol"] = df["gene_id"].astype(str)
    
    # For numeric IDs (Illumina HT-12 V4) or Affy IDs
    else:
        # Use probe IDs as-is for now — will be resolved in meta-analysis
        df["gene_symbol"] = df["gene_id"].astype(str)
    
    # Add disease category
    df["disease"] = disease
    
    # Deduplicate: keep probe with smallest p-value per gene
    if "gene_symbol" in df.columns and df["gene_symbol"].notna().any():
        valid = df[df["gene_symbol"].notna() & (df["gene_symbol"] != "")].copy()
        valid = valid.sort_values("pvalue").drop_duplicates(subset="gene_symbol", keep="first")
        
        output_path = os.path.join(RESULTS_DIR, f"{accession}_deg_genes.csv")
        valid.to_csv(output_path, index=False)
        
        n_sig = valid["significant"].sum() if "significant" in valid.columns else 0
        print(f"    Saved: {len(valid)} unique genes, {n_sig} significant")
        return valid
    
    return None


def main():
    print("=" * 80)
    print("Phase 3b: Probe-to-Gene Mapping (Quick Method)")
    print("=" * 80)
    
    datasets = [
        ("GSE25101", "AS"),
        ("GSE73754", "AS"),
        ("GSE18781", "AS"),
        ("GSE61281", "PsA"),
        ("GSE59071", "IBD"),
        ("GSE58667", "jSpA"),
        ("GSE221786", "AS"),
    ]
    
    all_results = {}
    for acc, disease in datasets:
        result = create_gene_level_results_with_relaxed_approach(acc, disease)
        if result is not None:
            all_results[acc] = result
    
    # Summary
    print(f"\n{'='*80}")
    print("GENE MAPPING SUMMARY")
    print(f"{'='*80}")
    
    for acc, df in all_results.items():
        n_with_symbol = df["gene_symbol"].notna().sum()
        n_sig = df["significant"].sum() if "significant" in df.columns else 0
        print(f"  {acc}: {len(df)} genes, {n_with_symbol} with symbols, {n_sig} significant")


if __name__ == "__main__":
    main()
