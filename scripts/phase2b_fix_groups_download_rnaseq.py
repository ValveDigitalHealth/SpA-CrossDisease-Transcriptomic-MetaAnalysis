#!/usr/bin/env python3
"""
Phase 2b: Fix sample group assignments and download RNA-seq data (GSE221786).
Corrects ambiguous GEO metadata sample labels.
"""

import os
import json
import pandas as pd
from pathlib import Path

BASE_DIR = Path("/home/user/workspace/project1_cross_disease_metaanalysis_spa")
METADATA_DIR = BASE_DIR / "data" / "metadata"
RAW_DIR = BASE_DIR / "data" / "raw"

# Manual group assignments based on GEO metadata review
GROUP_ASSIGNMENTS = {
    "GSE25101": {
        "disease_keywords": ["ankylosing spondylitis", "AS patient", "case"],
        "control_keywords": ["healthy control", "normal", "control"],
    },
    "GSE73754": {
        "disease_keywords": ["AS", "ankylosing spondylitis"],
        "control_keywords": ["HC", "healthy", "control"],
    },
    "GSE18781": {
        "disease_keywords": ["AS", "patient"],
        "control_keywords": ["control", "healthy"],
    },
    "GSE61281": {
        "disease_keywords": ["psoriatic arthritis", "PsA"],
        "control_keywords": ["healthy", "control"],
    },
    "GSE59071": {
        "disease_keywords": ["IBD", "Crohn", "ulcerative colitis", "patient"],
        "control_keywords": ["control", "healthy"],
    },
    "GSE58667": {
        "disease_keywords": ["SpA", "juvenile", "patient"],
        "control_keywords": ["control", "healthy"],
    },
}


def fix_sample_groups(accession, metadata_df):
    """Fix sample group assignments using keyword matching."""
    if accession not in GROUP_ASSIGNMENTS:
        print(f"  No manual assignments for {accession}")
        return metadata_df

    assignments = GROUP_ASSIGNMENTS[accession]
    disease_kws = [k.lower() for k in assignments["disease_keywords"]]
    control_kws = [k.lower() for k in assignments["control_keywords"]]

    def assign_group(row):
        row_str = ' '.join(str(v).lower() for v in row.values)
        if any(kw in row_str for kw in disease_kws):
            return "disease"
        elif any(kw in row_str for kw in control_kws):
            return "control"
        return "unknown"

    metadata_df['group'] = metadata_df.apply(assign_group, axis=1)
    print(f"  Groups: {metadata_df['group'].value_counts().to_dict()}")
    return metadata_df


def main():
    print("Phase 2b: Fix Sample Groups and Download RNA-seq")
    print("=" * 60)

    # Fix metadata for all datasets
    for accession in GROUP_ASSIGNMENTS:
        meta_path = METADATA_DIR / f"{accession}_metadata.csv"
        if not meta_path.exists():
            print(f"  {accession}: No metadata file found, skipping")
            continue

        print(f"\nFixing groups for {accession}...")
        meta_df = pd.read_csv(meta_path, index_col=0)
        meta_df = fix_sample_groups(accession, meta_df)
        meta_df.to_csv(meta_path)
        print(f"  Saved updated metadata: {meta_path}")

    # RNA-seq download note
    print("\nNote: GSE221786 (RNA-seq) requires separate download from GEO.")
    print("Download count matrix from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221786")
    print("Place in data/raw/GSE221786/")

    print("\nPhase 2b complete.")


if __name__ == "__main__":
    main()
