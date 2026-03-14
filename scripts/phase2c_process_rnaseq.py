#!/usr/bin/env python3
"""
Phase 2c: Process RNA-seq count files and prepare final dataset matrix.
"""

import os
import pandas as pd
import numpy as np
import gzip
import json

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
PROCESSED_DIR = os.path.join(BASE_DIR, "data/processed")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")


def process_gse141646():
    """Process GSE141646 RNA-seq counts (AS pre/post TNFi)."""
    print("\n=== Processing GSE141646 (AS RNA-seq) ===")
    
    filepath = os.path.join(RAW_DIR, "GSE141646_AS-TNF_counts.txt.gz")
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)
    
    # Remove QC rows
    qc_rows = [r for r in df.index if r.startswith('N_')]
    df = df.drop(qc_rows, errors='ignore')
    
    print(f"  Raw counts: {df.shape[0]} genes × {df.shape[1]} samples")
    print(f"  Columns: {list(df.columns[:6])}")
    
    # Classify samples: PRE vs POST
    sample_info = {}
    pre_samples = []
    post_samples = []
    
    for col in df.columns:
        if '_PRE' in col:
            sample_info[col] = {"group": "disease", "detail": f"AS pre-treatment: {col}"}
            pre_samples.append(col)
        elif '_POST' in col:
            sample_info[col] = {"group": "disease_treated", "detail": f"AS post-treatment: {col}"}
            post_samples.append(col)
    
    print(f"  PRE-treatment (disease): {len(pre_samples)}")
    print(f"  POST-treatment (treated): {len(post_samples)}")
    
    # For meta-analysis, use PRE samples as "AS disease" — no healthy controls in this dataset
    # We'll note this limitation — it can still contribute via effect size relative to post-treatment
    
    # Save full counts
    df.to_csv(os.path.join(PROCESSED_DIR, "GSE141646_counts.csv"))
    
    # Save sample info
    with open(os.path.join(METADATA_DIR, "GSE141646_samples.json"), "w") as f:
        json.dump(sample_info, f, indent=2)
    
    print(f"  Saved: {df.shape[0]} genes × {df.shape[1]} samples")
    return df, sample_info


def process_gse221786():
    """Process GSE221786 RNA-seq (AS sex-stratified)."""
    print("\n=== Processing GSE221786 (AS sex-stratified RNA-seq) ===")
    
    filepath = os.path.join(RAW_DIR, "GSE221786_matrix.txt.gz")
    df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0)
    
    # Has a 'gene' column with gene names
    gene_col = None
    if 'gene' in df.columns:
        gene_col = df['gene']
        df = df.drop('gene', axis=1)
    
    print(f"  Counts: {df.shape[0]} genes × {df.shape[1]} samples")
    print(f"  Columns: {list(df.columns[:10])}")
    
    # Classify: Ctrl_F, Ctrl_M = control; AS_F, AS_M = disease
    sample_info = {}
    for col in df.columns:
        if col.startswith('Ctrl'):
            sample_info[col] = {"group": "control", "detail": f"Healthy control: {col}"}
        elif col.startswith('AS'):
            sample_info[col] = {"group": "disease", "detail": f"AS patient: {col}"}
        else:
            sample_info[col] = {"group": "unknown", "detail": col}
    
    n_disease = sum(1 for s in sample_info.values() if s["group"] == "disease")
    n_control = sum(1 for s in sample_info.values() if s["group"] == "control")
    print(f"  AS (disease): {n_disease}, Ctrl (control): {n_control}")
    
    # Save with gene names as additional column
    if gene_col is not None:
        df.insert(0, 'gene_symbol', gene_col.values)
    df.to_csv(os.path.join(PROCESSED_DIR, "GSE221786_counts.csv"))
    
    with open(os.path.join(METADATA_DIR, "GSE221786_samples.json"), "w") as f:
        json.dump(sample_info, f, indent=2)
    
    print(f"  Saved: {df.shape}")
    return df, sample_info


def summarize_all_datasets():
    """Create a final summary of all usable datasets."""
    print("\n" + "=" * 80)
    print("FINAL DATASET SUMMARY FOR META-ANALYSIS")
    print("=" * 80)
    
    datasets = []
    
    # Microarray datasets (from series matrix)
    microarray_info = {
        "GSE25101": {"disease": "AS", "tissue": "Whole blood", "platform": "Illumina HumanHT-12 V3.0"},
        "GSE73754": {"disease": "AS", "tissue": "Whole blood", "platform": "Illumina HumanRef-8 v3.0"},
        "GSE18781": {"disease": "AS/axSpA", "tissue": "Whole blood", "platform": "Affymetrix HG-U133 Plus 2.0"},
        "GSE61281": {"disease": "PsA", "tissue": "Whole blood", "platform": "Affymetrix HuGene 1.0 ST"},
        "GSE59071": {"disease": "IBD", "tissue": "Whole blood", "platform": "Illumina HumanHT-12 V4.0"},
        "GSE58667": {"disease": "jSpA", "tissue": "Blood", "platform": "Affymetrix HG-U133 Plus 2.0"},
    }
    
    for acc, info in microarray_info.items():
        expr_path = os.path.join(PROCESSED_DIR, f"{acc}_expression.csv")
        meta_path = os.path.join(METADATA_DIR, f"{acc}_samples.json")
        
        if os.path.exists(expr_path):
            df = pd.read_csv(expr_path, index_col=0, nrows=2)
            n_genes = sum(1 for _ in open(expr_path)) - 1
            n_samples = df.shape[1]
            
            n_disease, n_control = 0, 0
            if os.path.exists(meta_path):
                with open(meta_path) as f:
                    samples = json.load(f)
                n_disease = sum(1 for s in samples.values() if s["group"] == "disease")
                n_control = sum(1 for s in samples.values() if s["group"] == "control")
            
            entry = {
                "accession": acc,
                "disease": info["disease"],
                "data_type": "Microarray",
                "tissue": info["tissue"],
                "platform": info["platform"],
                "n_genes": n_genes,
                "n_samples": n_samples,
                "n_disease": n_disease,
                "n_control": n_control,
                "usable": n_disease >= 5 and n_control >= 3,
            }
            datasets.append(entry)
    
    # RNA-seq datasets
    rnaseq_info = {
        "GSE221786": {"disease": "AS", "tissue": "Blood", "platform": "Illumina"},
        "GSE141646": {"disease": "AS", "tissue": "PBMC", "platform": "Illumina HiSeq 4000"},
    }
    
    for acc, info in rnaseq_info.items():
        counts_path = os.path.join(PROCESSED_DIR, f"{acc}_counts.csv")
        meta_path = os.path.join(METADATA_DIR, f"{acc}_samples.json")
        
        if os.path.exists(counts_path):
            df = pd.read_csv(counts_path, index_col=0, nrows=2)
            n_genes = sum(1 for _ in open(counts_path)) - 1
            
            # Exclude gene_symbol column if present
            data_cols = [c for c in df.columns if c != 'gene_symbol']
            n_samples = len(data_cols)
            
            n_disease, n_control = 0, 0
            if os.path.exists(meta_path):
                with open(meta_path) as f:
                    samples = json.load(f)
                n_disease = sum(1 for s in samples.values() if s["group"] == "disease")
                n_control = sum(1 for s in samples.values() if s["group"] == "control")
            
            entry = {
                "accession": acc,
                "disease": info["disease"],
                "data_type": "RNA-seq",
                "tissue": info["tissue"],
                "platform": info["platform"],
                "n_genes": n_genes,
                "n_samples": n_samples,
                "n_disease": n_disease,
                "n_control": n_control,
                "usable": n_disease >= 5 and n_control >= 3,
            }
            datasets.append(entry)
    
    # Print summary table
    print(f"\n{'Accession':12s} {'Disease':8s} {'Type':10s} {'Tissue':15s} {'Genes':>7s} {'Samples':>8s} {'Disease':>8s} {'Control':>8s} {'Usable':>7s}")
    print("-" * 100)
    
    usable_count = 0
    for d in datasets:
        usable = "✓" if d["usable"] else "✗"
        if d["usable"]:
            usable_count += 1
        print(f"{d['accession']:12s} {d['disease']:8s} {d['data_type']:10s} {d['tissue']:15s} "
              f"{d['n_genes']:>7d} {d['n_samples']:>8d} {d['n_disease']:>8d} {d['n_control']:>8d} {usable:>7s}")
    
    print(f"\nUsable datasets: {usable_count}/{len(datasets)}")
    print(f"Total samples across usable datasets: {sum(d['n_samples'] for d in datasets if d['usable'])}")
    
    # Save final inventory
    inv_path = os.path.join(METADATA_DIR, "final_dataset_inventory.json")
    with open(inv_path, "w") as f:
        json.dump(datasets, f, indent=2)
    print(f"\nFinal inventory saved to: {inv_path}")
    
    return datasets


if __name__ == "__main__":
    process_gse141646()
    process_gse221786()
    datasets = summarize_all_datasets()
