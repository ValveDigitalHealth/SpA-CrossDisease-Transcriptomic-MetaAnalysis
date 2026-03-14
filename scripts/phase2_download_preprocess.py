#!/usr/bin/env python3
"""
Phase 2: Download and preprocess GEO expression datasets.
Downloads series matrix files and applies normalization.
"""

import os
import gzip
import urllib.request
import numpy as np
import pandas as pd
from pathlib import Path

BASE_DIR = Path("/home/user/workspace/project1_cross_disease_metaanalysis_spa")
RAW_DIR = BASE_DIR / "data" / "raw"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
METADATA_DIR = BASE_DIR / "data" / "metadata"

RAW_DIR.mkdir(parents=True, exist_ok=True)
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

DATASET_INFO = {
    "GSE25101": {"disease": "AS", "platform": "microarray"},
    "GSE73754": {"disease": "AS", "platform": "microarray"},
    "GSE18781": {"disease": "AS", "platform": "microarray"},
    "GSE61281": {"disease": "PsA", "platform": "microarray"},
    "GSE59071": {"disease": "IBD", "platform": "microarray"},
    "GSE58667": {"disease": "jSpA", "platform": "microarray"},
}


def download_geo_matrix(accession, output_dir):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:6]}nnn/{accession}/matrix/{accession}_series_matrix.txt.gz"
    out_path = output_dir / f"{accession}_series_matrix.txt.gz"
    if out_path.exists():
        print(f"  {accession}: Already downloaded")
        return out_path
    print(f"  Downloading {accession}...")
    try:
        urllib.request.urlretrieve(url, out_path)
        print(f"  {accession}: Downloaded")
        return out_path
    except Exception as e:
        print(f"  ERROR downloading {accession}: {e}")
        return None


def parse_series_matrix(filepath):
    """Parse a GEO series matrix file into expression matrix and metadata."""
    samples = []
    sample_chars = {}
    expression_data = []
    probe_ids = []
    in_table = False

    opener = gzip.open if str(filepath).endswith('.gz') else open
    with opener(filepath, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('!Series_sample_id'):
                samples = line.split('\t')[1:]
                samples = [s.strip('"') for s in samples]
            elif line.startswith('!Sample_geo_accession'):
                gsms = [s.strip('"') for s in line.split('\t')[1:]]
                for gsm in gsms:
                    if gsm not in sample_chars:
                        sample_chars[gsm] = {}
                    sample_chars[gsm]['gsm'] = gsm
            elif line.startswith('!Sample_characteristics_ch1'):
                chars = [s.strip('"') for s in line.split('\t')[1:]]
                for i, gsm in enumerate(gsms if 'gsms' in dir() else []):
                    if i < len(chars):
                        char = chars[i]
                        if ':' in char:
                            key, val = char.split(':', 1)
                            key = key.strip().lower().replace(' ', '_')
                            val = val.strip()
                            if gsm not in sample_chars:
                                sample_chars[gsm] = {}
                            sample_chars[gsm][key] = val
            elif line.startswith('!dataset_table_begin') or line == '"ID_REF"\t' or line.startswith('ID_REF'):
                in_table = True
                continue
            elif line.startswith('!dataset_table_end') or line.startswith('!series_matrix_table_end'):
                in_table = False
            elif in_table and line.strip():
                parts = line.split('\t')
                if parts[0].strip('"') == 'ID_REF':
                    continue
                probe_id = parts[0].strip('"')
                try:
                    values = [float(v) if v.strip() and v.strip() != 'null' else np.nan
                              for v in parts[1:]]
                    if len(values) > 0:
                        probe_ids.append(probe_id)
                        expression_data.append(values)
                except (ValueError, IndexError):
                    pass

    if not expression_data:
        return None, None

    expr_df = pd.DataFrame(expression_data, index=probe_ids)
    if samples:
        expr_df.columns = samples[:expr_df.shape[1]]

    meta_df = pd.DataFrame.from_dict(sample_chars, orient='index')

    return expr_df, meta_df


def normalize_expression(expr_df):
    """Apply log2 transformation if needed and quantile normalization."""
    # Check if already log-transformed
    max_val = expr_df.max().max()
    if max_val > 100:
        expr_df = np.log2(expr_df + 1)
        print("    Applied log2 transformation")

    # Quantile normalization
    rank_mean = expr_df.stack().groupby(expr_df.stack().apply(lambda x: pd.Series(expr_df[expr_df.columns[0]].rank(method='first')).iloc[0] if False else x)).mean()
    # Simple approach: sort each column, average, reorder
    sorted_df = pd.DataFrame(np.sort(expr_df.values, axis=0), columns=expr_df.columns)
    row_means = sorted_df.mean(axis=1)
    quantile_norm = expr_df.rank(method='first').apply(lambda col: col.map(lambda r: row_means.iloc[int(r)-1] if not np.isnan(r) else np.nan))
    print("    Applied quantile normalization")
    return quantile_norm


def main():
    print("Phase 2: Download and Preprocess GEO Datasets")
    print("=" * 60)

    for accession, info in DATASET_INFO.items():
        print(f"\nProcessing {accession} ({info['disease']})...")

        # Download
        raw_path = download_geo_matrix(accession, RAW_DIR)
        if raw_path is None:
            print(f"  FAILED: Could not download {accession}")
            continue

        # Parse
        print(f"  Parsing matrix file...")
        expr_df, meta_df = parse_series_matrix(raw_path)
        if expr_df is None:
            print(f"  FAILED: Could not parse {accession}")
            continue

        print(f"  Matrix shape: {expr_df.shape}")

        # Normalize
        print(f"  Normalizing...")
        norm_df = normalize_expression(expr_df)

        # Save
        out_path = PROCESSED_DIR / f"{accession}_normalized.csv"
        norm_df.to_csv(out_path)
        print(f"  Saved: {out_path}")

        if meta_df is not None and not meta_df.empty:
            meta_path = METADATA_DIR / f"{accession}_metadata.csv"
            meta_df.to_csv(meta_path)
            print(f"  Metadata: {meta_path}")

    print("\nPhase 2 complete.")


if __name__ == "__main__":
    main()
