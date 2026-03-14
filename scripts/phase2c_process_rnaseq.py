#!/usr/bin/env python3
"""
Phase 2c: Process RNA-seq count files and prepare final expression matrix.
Applies TMM normalization and log2-CPM transformation for GSE221786.
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path

BASE_DIR = Path("/home/user/workspace/project1_cross_disease_metaanalysis_spa")
RAW_DIR = BASE_DIR / "data" / "raw" / "GSE221786"
PROCESSED_DIR = BASE_DIR / "data" / "processed"
METADATA_DIR = BASE_DIR / "data" / "metadata"


def compute_tmm_factors(count_matrix):
    """Compute TMM normalization factors."""
    # Use geometric mean as reference sample
    log_counts = np.log2(count_matrix + 0.5)
    ref_sample = log_counts.mean(axis=1)

    tmm_factors = []
    for col in count_matrix.columns:
        sample = count_matrix[col]
        # M values (log ratio)
        m = np.log2((sample + 0.5) / (count_matrix.sum(axis=1) / sample.sum() + 0.5))
        # A values (mean average)
        a = 0.5 * (np.log2(sample + 0.5) + np.log2(count_matrix.sum(axis=1) / sample.sum() + 0.5))
        # Filter by 30th/70th percentile of M
        m_30 = np.percentile(m.dropna(), 30)
        m_70 = np.percentile(m.dropna(), 70)
        a_95 = np.percentile(a.dropna(), 95)
        mask = (m >= m_30) & (m <= m_70) & (a <= a_95)
        tmm_factor = np.exp(np.log(2) * m[mask].mean()) if mask.sum() > 0 else 1.0
        tmm_factors.append(tmm_factor)

    return np.array(tmm_factors)


def log2_cpm(count_matrix, tmm_factors=None):
    """Convert counts to log2-CPM."""
    lib_sizes = count_matrix.sum(axis=0)
    if tmm_factors is not None:
        lib_sizes = lib_sizes * tmm_factors
    # CPM
    cpm = count_matrix.divide(lib_sizes, axis=1) * 1e6
    # log2 + 1
    log2cpm = np.log2(cpm + 1)
    return log2cpm


def main():
    print("Phase 2c: Process RNA-seq Data (GSE221786)")
    print("=" * 60)

    count_files = list(RAW_DIR.glob("*.txt")) + list(RAW_DIR.glob("*.tsv")) + list(RAW_DIR.glob("*.csv"))
    if not count_files:
        print(f"No count files found in {RAW_DIR}")
        print("Please download GSE221786 count matrix from GEO first")
        return

    print(f"Found {len(count_files)} count file(s)")
    count_file = count_files[0]
    print(f"Using: {count_file}")

    # Load count matrix
    count_df = pd.read_csv(count_file, sep='\t', index_col=0)
    print(f"Count matrix shape: {count_df.shape}")

    # Remove version suffix from Ensembl IDs if present
    count_df.index = count_df.index.str.split('.').str[0]

    # Filter low-count genes (keep genes with CPM > 1 in at least half the samples)
    cpm_raw = count_df.divide(count_df.sum(axis=0), axis=1) * 1e6
    keep = (cpm_raw > 1).sum(axis=1) >= (count_df.shape[1] / 2)
    count_df = count_df[keep]
    print(f"After filtering low-count genes: {count_df.shape}")

    # Compute TMM factors
    print("Computing TMM normalization factors...")
    tmm_factors = compute_tmm_factors(count_df)

    # Compute log2-CPM
    print("Computing log2-CPM...")
    log2cpm_df = log2_cpm(count_df, tmm_factors)

    # Save
    out_path = PROCESSED_DIR / "GSE221786_normalized.csv"
    log2cpm_df.to_csv(out_path)
    print(f"Saved: {out_path}")
    print(f"Final expression matrix: {log2cpm_df.shape}")

    print("\nPhase 2c complete.")


if __name__ == "__main__":
    main()
