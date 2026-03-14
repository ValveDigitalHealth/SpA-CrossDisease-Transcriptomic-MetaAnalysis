#!/usr/bin/env python3
"""
Phase 3b (quick): Map probe IDs to gene symbols using g:Profiler API.
More reliable than downloading GEO platform annotations.
"""

import os
import requests
import pandas as pd
import numpy as np

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RESULTS_DIR = os.path.join(BASE_DIR, "results/deg")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")


def use_gprofiler_convert(probe_ids, organism="hsapiens"):
    """Use g:Profiler gconvert API to map IDs to gene symbols."""
    url = "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
    all_mappings = {}

    for i in range(0, len(probe_ids), 500):
        batch = probe_ids[i:i+500]
        payload = {"organism": organism, "target": "ENSG", "query": batch}
        try:
            resp = requests.post(url, json=payload, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                for r in data.get("result", []):
                    incoming = r.get("incoming", "")
                    name = r.get("name", "")
                    if name and name != "N/A":
                        all_mappings[incoming] = name
        except Exception as e:
            print(f"  g:Profiler error: {e}")

    return all_mappings


def create_gene_level_results(accession, disease):
    """Create gene-level DEG results with probe mapping."""
    deg_path = os.path.join(RESULTS_DIR, f"{accession}_deg_results.csv")
    if not os.path.exists(deg_path):
        return None

    df = pd.read_csv(deg_path)
    sample_ids = df["gene_id"].astype(str).head(5).tolist()
    print(f"  {accession} sample IDs: {sample_ids}")

    if any("ENSG" in str(x) for x in sample_ids):
        counts_path = os.path.join(BASE_DIR, "data/processed/GSE221786_counts.csv")
        if os.path.exists(counts_path):
            counts_df = pd.read_csv(counts_path, index_col=0, usecols=[0, 1])
            if "gene_symbol" in counts_df.columns:
                mapping = dict(zip(counts_df.index.astype(str), counts_df["gene_symbol"]))
                df["gene_symbol"] = df["gene_id"].astype(str).map(mapping)
    elif any("ILMN_" in str(x) for x in sample_ids):
        probe_ids = df["gene_id"].astype(str).tolist()
        mapping = use_gprofiler_convert(probe_ids[:2000])
        if mapping:
            df["gene_symbol"] = df["gene_id"].astype(str).map(mapping)
        else:
            df["gene_symbol"] = df["gene_id"].astype(str)
    else:
        df["gene_symbol"] = df["gene_id"].astype(str)

    df["disease"] = disease

    if "gene_symbol" in df.columns:
        valid = df[df["gene_symbol"].notna() & (df["gene_symbol"] != "")].copy()
        valid = valid.sort_values("pvalue").drop_duplicates(subset="gene_symbol", keep="first")
        out_path = os.path.join(RESULTS_DIR, f"{accession}_deg_genes.csv")
        valid.to_csv(out_path, index=False)
        print(f"  Saved: {len(valid)} unique genes")
        return valid
    return None


def main():
    print("Phase 3b (quick): Probe-to-Gene Mapping")
    datasets = [
        ("GSE25101", "AS"), ("GSE73754", "AS"), ("GSE18781", "AS"),
        ("GSE61281", "PsA"), ("GSE59071", "IBD"), ("GSE58667", "jSpA"),
        ("GSE221786", "AS"),
    ]
    for acc, disease in datasets:
        print(f"\n{acc}:")
        create_gene_level_results(acc, disease)


if __name__ == "__main__":
    main()
