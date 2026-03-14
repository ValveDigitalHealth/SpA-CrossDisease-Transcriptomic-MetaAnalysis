#!/usr/bin/env python3
"""
Phase 2: Download and preprocess GEO expression datasets for SpA meta-analysis.
Uses NCBI GEO SOFT format to get processed expression matrices.
Focuses on bulk datasets with disease vs control design.
"""

import os
import gzip
import io
import re
import csv
import json
import requests
import time
import numpy as np
import pandas as pd
from datetime import datetime

# Configuration
BASE_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa"
RAW_DIR = os.path.join(BASE_DIR, "data/raw")
PROCESSED_DIR = os.path.join(BASE_DIR, "data/processed")
METADATA_DIR = os.path.join(BASE_DIR, "data/metadata")
LOG_DIR = os.path.join(BASE_DIR, "logs")

for d in [RAW_DIR, PROCESSED_DIR, METADATA_DIR, LOG_DIR]:
    os.makedirs(d, exist_ok=True)

# Priority datasets for meta-analysis (bulk, blood/PBMC, disease vs control)
PRIORITY_DATASETS = [
    # AS
    "GSE25101",   # 16 AS + 16 HC, microarray, whole blood
    "GSE73754",   # 52 AS + 20 HC, microarray, whole blood
    "GSE141646",  # 22 AS pre-treatment + 22 post, RNA-seq, PBMC
    "GSE221786",  # ~18 AS + 10 HC, RNA-seq, blood
    "GSE18781",   # axial SpA blood, microarray
    # PsA
    "GSE61281",   # PsA blood, microarray
    "GSE137510",  # PsA PBMC, RNA-seq
    # IBD
    "GSE59071",   # IBD blood, microarray
    # Uveitis
    "GSE194060",  # HLA-B27+ AAU, DCs from blood, RNA-seq
    "GSE195501",  # Uveitis subtypes, blood, RNA-seq
    # SpA general
    "GSE58667",   # Juvenile SpA, blood, microarray
]


def download_geo_series_matrix(accession):
    """Download series matrix file from GEO."""
    # Try series matrix first (contains processed expression data)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/matrix/{accession}_series_matrix.txt.gz"
    
    output_path = os.path.join(RAW_DIR, f"{accession}_series_matrix.txt.gz")
    
    if os.path.exists(output_path):
        print(f"  Already downloaded: {output_path}")
        return output_path
    
    print(f"  Downloading series matrix from: {url}")
    try:
        resp = requests.get(url, timeout=120, stream=True)
        if resp.status_code == 200:
            with open(output_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
            size_mb = os.path.getsize(output_path) / 1024 / 1024
            print(f"  Downloaded: {size_mb:.1f} MB")
            return output_path
        else:
            print(f"  HTTP {resp.status_code} — trying alternate URL...")
    except Exception as e:
        print(f"  Download error: {e}")
    
    # Try alternate URL pattern (some series have multiple matrix files)
    url2 = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/matrix/"
    try:
        resp = requests.get(url2, timeout=30)
        if resp.status_code == 200:
            # Find matrix file names
            import re
            files = re.findall(r'href="([^"]+_series_matrix\.txt\.gz)"', resp.text)
            if files:
                for fname in files:
                    full_url = url2 + fname
                    print(f"  Trying: {full_url}")
                    resp2 = requests.get(full_url, timeout=120, stream=True)
                    if resp2.status_code == 200:
                        out = os.path.join(RAW_DIR, fname)
                        with open(out, "wb") as f:
                            for chunk in resp2.iter_content(chunk_size=8192):
                                f.write(chunk)
                        size_mb = os.path.getsize(out) / 1024 / 1024
                        print(f"  Downloaded: {size_mb:.1f} MB")
                        return out
    except Exception as e:
        print(f"  Alternate URL error: {e}")
    
    print(f"  FAILED to download {accession}")
    return None


def parse_series_matrix(filepath):
    """Parse a GEO series matrix file into expression matrix + sample metadata."""
    if filepath is None:
        return None, None
    
    print(f"  Parsing: {os.path.basename(filepath)}")
    
    metadata_lines = []
    data_lines = []
    in_data = False
    
    opener = gzip.open if filepath.endswith(".gz") else open
    
    with opener(filepath, "rt", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line.startswith("!"):
                metadata_lines.append(line)
            elif line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                in_data = True
                data_lines.append(line)
            elif in_data and line and not line.startswith("!"):
                data_lines.append(line)
    
    # Parse metadata
    metadata = {}
    for line in metadata_lines:
        if line.startswith("!Sample_"):
            key = line.split("\t")[0].replace("!Sample_", "")
            values = line.split("\t")[1:]
            values = [v.strip('"') for v in values]
            if key not in metadata:
                metadata[key] = values
    
    # Parse expression data
    if not data_lines:
        print(f"  WARNING: No expression data found in {filepath}")
        return None, metadata
    
    try:
        data_str = "\n".join(data_lines)
        df = pd.read_csv(io.StringIO(data_str), sep="\t", index_col=0)
        # Clean column names
        df.columns = [c.strip('"') for c in df.columns]
        # Remove completely empty rows
        df = df.dropna(how="all")
        print(f"  Expression matrix: {df.shape[0]} probes × {df.shape[1]} samples")
        return df, metadata
    except Exception as e:
        print(f"  ERROR parsing expression data: {e}")
        return None, metadata


def extract_sample_groups(metadata, accession):
    """Try to identify disease vs control groups from metadata."""
    # Check various metadata fields for group info
    group_fields = ["characteristics_ch1", "source_name_ch1", "title", "description"]
    
    sample_info = {}
    n_samples = 0
    
    for field in group_fields:
        if field in metadata:
            values = metadata[field]
            n_samples = len(values)
            for i, v in enumerate(values):
                v_lower = v.lower()
                sample_id = metadata.get("geo_accession", [f"sample_{i}"])[i] if "geo_accession" in metadata else f"sample_{i}"
                
                if any(kw in v_lower for kw in ["control", "healthy", "normal", "hc", "non-disease"]):
                    sample_info[sample_id] = {"group": "control", "detail": v}
                elif any(kw in v_lower for kw in ["ankylosing", "spondylitis", " as ", "as_", "axial spa", 
                                                     "psoriatic arthritis", "psa", "ibd", "crohn", "colitis",
                                                     "uveitis", "spondyloarth", "disease", "patient", "case"]):
                    sample_info[sample_id] = {"group": "disease", "detail": v}
    
    if not sample_info and n_samples > 0:
        # If we couldn't auto-classify, store all info for manual review
        for field in group_fields:
            if field in metadata:
                for i, v in enumerate(metadata[field]):
                    sample_id = metadata.get("geo_accession", [f"sample_{i}"])[i] if "geo_accession" in metadata else f"sample_{i}"
                    if sample_id not in sample_info:
                        sample_info[sample_id] = {"group": "unknown", "detail": v}
    
    return sample_info


def log_dataset_status(accession, status, details=""):
    """Log download/processing status."""
    log_path = os.path.join(LOG_DIR, "phase2_download_log.csv")
    exists = os.path.exists(log_path)
    with open(log_path, "a", newline="") as f:
        writer = csv.writer(f)
        if not exists:
            writer.writerow(["timestamp", "accession", "status", "details"])
        writer.writerow([datetime.now().isoformat(), accession, status, details])


def main():
    print("=" * 80)
    print("Phase 2: Download and Preprocess GEO Datasets")
    print(f"Timestamp: {datetime.now().strftime('%Y%m%d_%H%M%S')}")
    print(f"Priority datasets: {len(PRIORITY_DATASETS)}")
    print("=" * 80)
    
    results = {}
    
    for acc in PRIORITY_DATASETS:
        print(f"\n{'='*60}")
        print(f"Processing: {acc}")
        print(f"{'='*60}")
        
        # Download
        filepath = download_geo_series_matrix(acc)
        
        if filepath is None:
            log_dataset_status(acc, "DOWNLOAD_FAILED")
            results[acc] = {"status": "failed", "reason": "download"}
            continue
        
        # Parse
        expr_df, metadata = parse_series_matrix(filepath)
        
        if expr_df is None:
            log_dataset_status(acc, "PARSE_FAILED", "No expression data in matrix file")
            results[acc] = {"status": "failed", "reason": "parse"}
            continue
        
        # Extract sample groups
        sample_info = extract_sample_groups(metadata, acc)
        
        # Save processed expression matrix
        expr_path = os.path.join(PROCESSED_DIR, f"{acc}_expression.csv")
        expr_df.to_csv(expr_path)
        
        # Save sample metadata
        if sample_info:
            meta_path = os.path.join(METADATA_DIR, f"{acc}_samples.json")
            with open(meta_path, "w") as f:
                json.dump(sample_info, f, indent=2)
        
        # Save full metadata
        meta_full_path = os.path.join(METADATA_DIR, f"{acc}_metadata.json")
        with open(meta_full_path, "w") as f:
            json.dump(metadata, f, indent=2, default=str)
        
        n_disease = sum(1 for s in sample_info.values() if s["group"] == "disease")
        n_control = sum(1 for s in sample_info.values() if s["group"] == "control")
        n_unknown = sum(1 for s in sample_info.values() if s["group"] == "unknown")
        
        results[acc] = {
            "status": "success",
            "shape": list(expr_df.shape),
            "n_disease": n_disease,
            "n_control": n_control,
            "n_unknown": n_unknown,
            "expr_path": expr_path,
        }
        
        log_dataset_status(acc, "SUCCESS", f"shape={expr_df.shape}, disease={n_disease}, control={n_control}")
        
        print(f"  Saved expression matrix: {expr_df.shape}")
        print(f"  Groups: disease={n_disease}, control={n_control}, unknown={n_unknown}")
        
        time.sleep(1)  # Rate limiting
    
    # Summary
    print(f"\n{'='*80}")
    print("PHASE 2 SUMMARY")
    print(f"{'='*80}")
    
    success = sum(1 for r in results.values() if r["status"] == "success")
    failed = sum(1 for r in results.values() if r["status"] == "failed")
    
    print(f"Successfully processed: {success}/{len(PRIORITY_DATASETS)}")
    print(f"Failed: {failed}")
    
    for acc, r in results.items():
        status_icon = "✓" if r["status"] == "success" else "✗"
        if r["status"] == "success":
            print(f"  {status_icon} {acc}: {r['shape'][0]} genes × {r['shape'][1]} samples "
                  f"(disease={r['n_disease']}, control={r['n_control']}, unknown={r['n_unknown']})")
        else:
            print(f"  {status_icon} {acc}: FAILED ({r.get('reason', 'unknown')})")
    
    # Save results summary
    summary_path = os.path.join(METADATA_DIR, "phase2_results.json")
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {summary_path}")


if __name__ == "__main__":
    main()
