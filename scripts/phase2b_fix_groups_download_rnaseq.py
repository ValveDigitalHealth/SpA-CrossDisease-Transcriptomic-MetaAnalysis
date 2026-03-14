#!/usr/bin/env python3
"""
Phase 2b: Fix sample group assignments and download RNA-seq supplementary data.
"""

import os
import json
import requests
import gzip
import time
import re
import pandas as pd

BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
PROCESSED_DIR = os.path.join(BASE_DIR, "data/processed")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")


def fix_gse73754():
    """Fix GSE73754 groups — AS vs HC is in the title field."""
    print("\n=== Fixing GSE73754 sample groups ===")
    with open(os.path.join(METADATA_DIR, "GSE73754_metadata.json")) as f:
        meta = json.load(f)
    
    sample_info = {}
    for i, (acc, title) in enumerate(zip(meta["geo_accession"], meta["title"])):
        if title.startswith("AS"):
            sample_info[acc] = {"group": "disease", "detail": title}
        elif title.startswith("HC"):
            sample_info[acc] = {"group": "control", "detail": title}
        else:
            sample_info[acc] = {"group": "unknown", "detail": title}
    
    n_disease = sum(1 for s in sample_info.values() if s["group"] == "disease")
    n_control = sum(1 for s in sample_info.values() if s["group"] == "control")
    print(f"  AS (disease): {n_disease}, HC (control): {n_control}")
    
    with open(os.path.join(METADATA_DIR, "GSE73754_samples.json"), "w") as f:
        json.dump(sample_info, f, indent=2)
    print("  Saved corrected sample groups")


def fix_gse141646():
    """Fix GSE141646 groups — pre-treatment AS vs post-treatment."""
    print("\n=== Fixing GSE141646 sample groups ===")
    with open(os.path.join(METADATA_DIR, "GSE141646_metadata.json")) as f:
        meta = json.load(f)
    
    sample_info = {}
    for key in ["characteristics_ch1", "source_name_ch1", "title", "description"]:
        if key in meta:
            print(f"  {key}: {list(set(meta[key]))[:5]}")
    
    # Classify based on title/characteristics
    for i, acc in enumerate(meta.get("geo_accession", [])):
        title = meta["title"][i] if "title" in meta else ""
        char = meta["characteristics_ch1"][i] if "characteristics_ch1" in meta else ""
        
        detail = f"{title} | {char}"
        
        if "pre" in title.lower() or "before" in title.lower() or "baseline" in title.lower():
            sample_info[acc] = {"group": "disease", "detail": detail, "subgroup": "AS_pre_treatment"}
        elif "post" in title.lower() or "after" in title.lower():
            sample_info[acc] = {"group": "disease_treated", "detail": detail, "subgroup": "AS_post_treatment"}
        else:
            sample_info[acc] = {"group": "unknown", "detail": detail}
    
    groups = {}
    for s in sample_info.values():
        g = s["group"]
        groups[g] = groups.get(g, 0) + 1
    print(f"  Groups: {groups}")
    
    with open(os.path.join(METADATA_DIR, "GSE141646_samples.json"), "w") as f:
        json.dump(sample_info, f, indent=2)


def download_supplementary_files(accession):
    """Download supplementary files for RNA-seq datasets."""
    print(f"\n=== Downloading supplementary files for {accession} ===")
    
    # Try to list supplementary files
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/suppl/"
    
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            # Find file links
            files = re.findall(r'href="([^"]+)"', resp.text)
            data_files = [f for f in files if any(ext in f.lower() for ext in 
                          ['.txt.gz', '.csv.gz', '.tsv.gz', '.xlsx', '_counts', '_expression', '_fpkm', '_tpm'])]
            
            print(f"  Found {len(data_files)} potential data files:")
            for f in data_files:
                print(f"    {f}")
            
            # Download the most likely expression file
            for fname in data_files:
                if any(kw in fname.lower() for kw in ['count', 'expression', 'fpkm', 'tpm', 'gene']):
                    full_url = url + fname
                    out_path = os.path.join(RAW_DIR, fname)
                    
                    if os.path.exists(out_path):
                        print(f"  Already exists: {out_path}")
                        return out_path
                    
                    print(f"  Downloading: {full_url}")
                    resp2 = requests.get(full_url, timeout=120, stream=True)
                    if resp2.status_code == 200:
                        with open(out_path, "wb") as f:
                            for chunk in resp2.iter_content(chunk_size=8192):
                                f.write(chunk)
                        size_mb = os.path.getsize(out_path) / 1024 / 1024
                        print(f"  Downloaded: {size_mb:.1f} MB")
                        return out_path
            
            # If no keyword match, download the largest file
            if data_files:
                fname = data_files[0]
                full_url = url + fname
                out_path = os.path.join(RAW_DIR, fname)
                
                if not os.path.exists(out_path):
                    print(f"  Downloading first file: {full_url}")
                    resp2 = requests.get(full_url, timeout=120, stream=True)
                    if resp2.status_code == 200:
                        with open(out_path, "wb") as f:
                            for chunk in resp2.iter_content(chunk_size=8192):
                                f.write(chunk)
                        size_mb = os.path.getsize(out_path) / 1024 / 1024
                        print(f"  Downloaded: {size_mb:.1f} MB")
                        return out_path
        else:
            print(f"  HTTP {resp.status_code} for supplementary listing")
    except Exception as e:
        print(f"  Error: {e}")
    
    return None


def try_parse_supplementary(filepath, accession):
    """Try to parse a supplementary expression file."""
    if filepath is None:
        return None
    
    print(f"  Trying to parse: {os.path.basename(filepath)}")
    
    try:
        if filepath.endswith('.gz'):
            df = pd.read_csv(filepath, sep='\t', compression='gzip', index_col=0, nrows=5)
        elif filepath.endswith('.xlsx'):
            df = pd.read_excel(filepath, index_col=0, nrows=5)
        else:
            df = pd.read_csv(filepath, sep='\t', index_col=0, nrows=5)
        
        print(f"  Preview: {df.shape[0]} rows × {df.shape[1]} columns")
        print(f"  Columns: {list(df.columns[:10])}")
        print(f"  Index: {list(df.index[:5])}")
        return True
    except Exception as e:
        print(f"  Parse error: {e}")
        # Try comma-separated
        try:
            if filepath.endswith('.gz'):
                df = pd.read_csv(filepath, sep=',', compression='gzip', index_col=0, nrows=5)
            else:
                df = pd.read_csv(filepath, sep=',', index_col=0, nrows=5)
            print(f"  CSV preview: {df.shape[0]} rows × {df.shape[1]} columns")
            return True
        except Exception as e2:
            print(f"  CSV parse also failed: {e2}")
            return False


def main():
    print("=" * 80)
    print("Phase 2b: Fix Groups & Download RNA-seq Data")
    print("=" * 80)
    
    # Fix group assignments
    fix_gse73754()
    fix_gse141646()
    
    # RNA-seq datasets that need supplementary files
    rnaseq_datasets = ["GSE141646", "GSE221786", "GSE137510", "GSE194060", "GSE195501"]
    
    for acc in rnaseq_datasets:
        filepath = download_supplementary_files(acc)
        if filepath:
            try_parse_supplementary(filepath, acc)
        time.sleep(1)
    
    print("\n" + "=" * 80)
    print("Phase 2b Complete")
    print("=" * 80)
    
    # Summary of what we have
    print("\nDatasets with expression data ready for analysis:")
    expr_files = [f for f in os.listdir(PROCESSED_DIR) if f.endswith("_expression.csv")]
    for f in sorted(expr_files):
        filepath = os.path.join(PROCESSED_DIR, f)
        try:
            # Read just header to count columns
            df = pd.read_csv(filepath, index_col=0, nrows=2)
            n_genes = sum(1 for _ in open(filepath)) - 1
            print(f"  {f}: {n_genes} genes × {df.shape[1]} samples")
        except:
            print(f"  {f}: (error reading)")
    
    print("\nSupplementary files downloaded:")
    for f in sorted(os.listdir(RAW_DIR)):
        size = os.path.getsize(os.path.join(RAW_DIR, f)) / 1024 / 1024
        print(f"  {f}: {size:.1f} MB")


if __name__ == "__main__":
    main()
