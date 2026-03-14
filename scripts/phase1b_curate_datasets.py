#!/usr/bin/env python3
"""
Phase 1b: Curate the final dataset list for SpA meta-analysis.
Filters for GSE-level series with expression profiling from blood/PBMC.
"""

import csv
import json
import os

OUTPUT_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa/data/metadata"

CURATED_DATASETS = [
    {"accession": "GSE25101", "disease": "AS", "tissue": "Whole blood", "platform": "Illumina HumanRef-8 v2.0", "data_type": "Microarray", "n_disease": 16, "n_control": 16, "year": 2011, "pmid": "21854295"},
    {"accession": "GSE73754", "disease": "AS", "tissue": "Whole blood", "platform": "Illumina HumanRef-8 v2.0", "data_type": "Microarray", "n_disease": 52, "n_control": 20, "year": 2016, "pmid": ""},
    {"accession": "GSE18781", "disease": "AS", "tissue": "Peripheral blood", "platform": "Affymetrix HG-U133 Plus 2.0", "data_type": "Microarray", "n_disease": 18, "n_control": 18, "year": 2010, "pmid": "20920073"},
    {"accession": "GSE221786", "disease": "AS", "tissue": "Blood", "platform": "Illumina NovaSeq 6000", "data_type": "RNA-seq", "n_disease": 30, "n_control": 30, "year": 2023, "pmid": ""},
    {"accession": "GSE61281", "disease": "PsA", "tissue": "Peripheral blood", "platform": "Agilent 4x44K", "data_type": "Microarray", "n_disease": 20, "n_control": 20, "year": 2014, "pmid": "25303302"},
    {"accession": "GSE59071", "disease": "IBD", "tissue": "Peripheral blood", "platform": "Affymetrix HG-U133 Plus 2.0", "data_type": "Microarray", "n_disease": 60, "n_control": 40, "year": 2015, "pmid": ""},
    {"accession": "GSE58667", "disease": "jSpA", "tissue": "Peripheral blood", "platform": "Illumina HumanHT-12 V4.0", "data_type": "Microarray", "n_disease": 30, "n_control": 20, "year": 2015, "pmid": ""},
]


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    csv_path = os.path.join(OUTPUT_DIR, "dataset_inventory.csv")
    json_path = os.path.join(OUTPUT_DIR, "dataset_inventory.json")

    with open(csv_path, "w", newline="") as f:
        fieldnames = ["accession", "disease", "tissue", "platform", "data_type", "n_disease", "n_control", "year", "pmid"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(CURATED_DATASETS)

    with open(json_path, "w") as f:
        json.dump(CURATED_DATASETS, f, indent=2)

    print(f"Dataset inventory saved: {len(CURATED_DATASETS)} datasets")
    print(f"  CSV: {csv_path}")
    print(f"  JSON: {json_path}")

    for ds in CURATED_DATASETS:
        print(f"  {ds['accession']}: {ds['disease']} | {ds['n_disease']}D/{ds['n_control']}C | {ds['platform']}")


if __name__ == "__main__":
    main()
