#!/usr/bin/env python3
"""
Phase 3b: Map probe IDs to gene symbols for all datasets.
Downloads platform annotation tables from GEO FTP.
"""

import os
import gzip
import requests
import pandas as pd
import time

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
RESULTS_DIR = os.path.join(BASE_DIR, "results/deg")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")

PLATFORM_MAP = {
    "GPL6947": {"datasets": ["GSE25101"], "name": "Illumina HumanHT-12 V3.0"},
    "GPL6883": {"datasets": ["GSE73754"], "name": "Illumina HumanRef-8 v3.0"},
    "GPL570":  {"datasets": ["GSE18781"], "name": "Affymetrix HG-U133 Plus 2.0"},
    "GPL6244": {"datasets": ["GSE61281"], "name": "Affymetrix HuGene 1.0 ST"},
    "GPL13667":{"datasets": ["GSE59071"], "name": "Illumina HumanHT-12 V4.0"},
    "GPL96":   {"datasets": ["GSE58667"], "name": "Affymetrix HG-U133A"},
}


def download_platform_annotation(gpl_id):
    """Download platform annotation from GEO FTP."""
    output_path = os.path.join(RAW_DIR, f"{gpl_id}_annot.txt.gz")
    if os.path.exists(output_path):
        return output_path

    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{gpl_id[:-3]}nnn/{gpl_id}/annot/{gpl_id}.annot.gz"
    print(f"  Downloading: {url}")
    try:
        resp = requests.get(url, timeout=60)
        if resp.status_code == 200:
            with open(output_path, "wb") as f:
                f.write(resp.content)
            return output_path
    except Exception as e:
        print(f"  Error: {e}")
    return None


def parse_platform_annotation(filepath, gpl_id):
    """Parse platform annotation to get probe -> gene symbol mapping."""
    if not filepath:
        return {}

    probe_to_gene = {}
    try:
        opener = gzip.open if filepath.endswith(".gz") else open
        with opener(filepath, "rt", errors="replace") as f:
            in_table = False
            header = None
            id_idx = gene_idx = None

            for line in f:
                line = line.strip()
                if line.startswith("!platform_table_begin"):
                    in_table = True
                    continue
                elif line.startswith("!platform_table_end"):
                    break

                if not in_table:
                    continue

                fields = line.split("\t")
                if header is None:
                    header = [h.lower().strip() for h in fields]
                    for i, h in enumerate(header):
                        if h == "id":
                            id_idx = i
                        elif h in ["gene symbol", "symbol", "gene_symbol"]:
                            gene_idx = i
                    continue

                if id_idx is not None and gene_idx is not None and gene_idx < len(fields):
                    probe = fields[id_idx].strip()
                    gene = fields[gene_idx].strip().split("///")[0].split("//")[0].strip()
                    if gene and gene != "---":
                        probe_to_gene[probe] = gene
    except Exception as e:
        print(f"  Parse error: {e}")

    print(f"  {gpl_id}: {len(probe_to_gene)} probe-gene mappings")
    return probe_to_gene


def map_deg_to_genes(accession, probe_to_gene):
    """Map DEG probe IDs to gene symbols, deduplicate by best p-value."""
    deg_path = os.path.join(RESULTS_DIR, f"{accession}_deg_results.csv")
    if not os.path.exists(deg_path):
        return None

    df = pd.read_csv(deg_path)
    df["gene_symbol"] = df["gene_id"].astype(str).map(probe_to_gene)

    n_mapped = df["gene_symbol"].notna().sum()
    print(f"  {accession}: {n_mapped}/{len(df)} mapped")

    mapped_df = df[df["gene_symbol"].notna()].copy()
    if len(mapped_df) > 0:
        mapped_df = mapped_df.sort_values("pvalue").drop_duplicates(subset="gene_symbol", keep="first")

    out_path = os.path.join(RESULTS_DIR, f"{accession}_deg_genes.csv")
    mapped_df.to_csv(out_path, index=False)
    return mapped_df


def main():
    print("Phase 3b: Probe-to-Gene Mapping")
    print("=" * 60)

    for gpl_id, info in PLATFORM_MAP.items():
        print(f"\n{gpl_id}: {info['name']}")
        filepath = download_platform_annotation(gpl_id)
        mapping = parse_platform_annotation(filepath, gpl_id)
        for acc in info["datasets"]:
            map_deg_to_genes(acc, mapping)
        time.sleep(0.5)

    print("\nPhase 3b complete.")


if __name__ == "__main__":
    main()
