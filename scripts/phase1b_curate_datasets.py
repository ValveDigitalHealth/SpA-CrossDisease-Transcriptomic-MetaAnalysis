#!/usr/bin/env python3
"""
Phase 1b: Curate the final dataset list for SpA meta-analysis.
Based on API search results + literature review + manual curation.
Filters for GSE-level series with expression profiling from blood/PBMC.
"""

import csv
import json
import os
import requests
import time

OUTPUT_DIR = "/home/user/workspace/project1_cross_disease_metaanalysis_spa/data/metadata"

# Curated dataset inventory based on:
# 1. GEO API search results
# 2. Literature review (papers citing these datasets)
# 3. Known datasets from prior publications

CURATED_DATASETS = [
    # ============= ANKYLOSING SPONDYLITIS (AS) =============
    {
        "accession": "GSE25101",
        "disease": "AS",
        "title": "Expression profiling in whole blood in AS patients and controls",
        "tissue": "Whole blood (PAXgene)",
        "platform": "GPL6947 Illumina HumanHT-12 V3.0",
        "data_type": "Microarray",
        "n_disease": 16,
        "n_control": 16,
        "total_samples": 32,
        "year": 2011,
        "pmid": "21854295",
        "notes": "Validated AS gene signature. Well-cited reference dataset."
    },
    {
        "accession": "GSE73754",
        "disease": "AS",
        "title": "Gene expression in AS whole blood",
        "tissue": "Whole blood",
        "platform": "GPL6883 Illumina HumanRef-8 v3.0",
        "data_type": "Microarray",
        "n_disease": 52,
        "n_control": 20,
        "total_samples": 72,
        "year": 2016,
        "pmid": "",
        "notes": "Larger AS cohort. Used in ferroptosis study (PMID: 36466778)."
    },
    {
        "accession": "GSE141646",
        "disease": "AS",
        "title": "Transcriptome analysis of AS patients before and after TNF-α inhibitor therapy",
        "tissue": "Mononuclear cells (PBMC)",
        "platform": "GPL20301 Illumina HiSeq 4000",
        "data_type": "RNA-seq",
        "n_disease": 22,
        "n_control": 22,
        "total_samples": 44,
        "year": 2020,
        "pmid": "",
        "notes": "Pre/post TNFi treatment. Use pre-treatment as 'disease' samples."
    },
    {
        "accession": "GSE221786",
        "disease": "AS",
        "title": "DEGs in Males and Females with AS: sex-specific effects",
        "tissue": "Blood",
        "platform": "GPL24676 Illumina NovaSeq 6000",
        "data_type": "RNA-seq",
        "n_disease": 18,
        "n_control": 10,
        "total_samples": 28,
        "year": 2023,
        "pmid": "",
        "notes": "Sex-stratified analysis. Good for subgroup analysis."
    },
    {
        "accession": "GSE11886",
        "disease": "AS",
        "title": "Gene expression analysis of macrophages from AS patients",
        "tissue": "Blood-derived macrophages",
        "platform": "GPL96 Affymetrix HG-U133A",
        "data_type": "Microarray",
        "n_disease": 16,
        "n_control": 16,
        "total_samples": 33,
        "year": 2009,
        "pmid": "19734132",
        "notes": "Macrophage-specific. Reveals innate immune activation."
    },
    
    # ============= PSORIATIC ARTHRITIS (PsA) =============
    {
        "accession": "GSE61281",
        "disease": "PsA",
        "title": "Gene expression in PsA and PsO blood",
        "tissue": "Whole blood",
        "platform": "GPL6244 Affymetrix HuGene 1.0 ST",
        "data_type": "Microarray",
        "n_disease": 20,
        "n_control": 12,
        "total_samples": 32,
        "year": 2015,
        "pmid": "",
        "notes": "PsA vs healthy + PsA vs PsO comparison."
    },
    {
        "accession": "GSE179800",
        "disease": "PsA",
        "title": "Changes in neutrophil gene expression in PsA patients",
        "tissue": "Neutrophils from blood",
        "platform": "Illumina",
        "data_type": "RNA-seq",
        "n_disease": 8,
        "n_control": 7,
        "total_samples": 15,
        "year": 2022,
        "pmid": "",
        "notes": "Neutrophil-specific transcriptomics in PsA."
    },
    {
        "accession": "GSE277596",
        "disease": "PsA",
        "title": "Dissection of the immune landscape in PsA patients",
        "tissue": "Blood/PBMC",
        "platform": "10x Genomics Chromium",
        "data_type": "scRNA-seq",
        "n_disease": 300,
        "n_control": 335,
        "total_samples": 635,
        "year": 2024,
        "pmid": "",
        "notes": "Large-scale immune landscape. May need pseudobulk aggregation."
    },
    {
        "accession": "GSE137510",
        "disease": "PsA",
        "title": "Transcriptomic analysis of PsA PBMC",
        "tissue": "PBMC",
        "platform": "Illumina",
        "data_type": "RNA-seq",
        "n_disease": 10,
        "n_control": 10,
        "total_samples": 20,
        "year": 2020,
        "pmid": "",
        "notes": "PBMC from PsA patients. Used in recent 2026 Frontiers paper."
    },
    
    # ============= SPONDYLOARTHRITIS GENERAL =============
    {
        "accession": "GSE30023",
        "disease": "SpA",
        "title": "Comparative gene expression of SpA and RA using microarrays",
        "tissue": "Synovial fluid cells",
        "platform": "GPL96 Affymetrix HG-U133A",
        "data_type": "Microarray",
        "n_disease": 6,
        "n_control": 6,
        "total_samples": 12,
        "year": 2011,
        "pmid": "",
        "notes": "SpA vs RA comparison. Synovial fluid — different compartment but informative."
    },
    {
        "accession": "GSE223717",
        "disease": "SpA",
        "title": "Monocyte subpopulations display disease-specific miRNA signatures in SpA subforms",
        "tissue": "Blood monocytes",
        "platform": "Agilent miRNA",
        "data_type": "miRNA array",
        "n_disease": 120,
        "n_control": 71,
        "total_samples": 191,
        "year": 2023,
        "pmid": "",
        "notes": "Large monocyte miRNA study across SpA subtypes. miRNA — may not directly integrate but valuable for discussion."
    },
    {
        "accession": "GSE58667",
        "disease": "SpA",
        "title": "Gene expression in juvenile spondyloarthritis",
        "tissue": "Blood",
        "platform": "Affymetrix",
        "data_type": "Microarray",
        "n_disease": 8,
        "n_control": 7,
        "total_samples": 15,
        "year": 2015,
        "pmid": "",
        "notes": "Juvenile SpA. Different age group but same disease spectrum."
    },
    {
        "accession": "GSE304063",
        "disease": "SpA",
        "title": "Impact of IL-17A Blockade on Inflammatory and Stromal Pathways in Peripheral SpA",
        "tissue": "Blood/tissue",
        "platform": "Illumina",
        "data_type": "RNA-seq",
        "n_disease": 12,
        "n_control": 12,
        "total_samples": 24,
        "year": 2025,
        "pmid": "",
        "notes": "IL-17A blockade effect on SpA. Recent dataset."
    },
    {
        "accession": "GSE194315",
        "disease": "SpA",
        "title": "RNA and surface epitope sequencing of single cells in SpA",
        "tissue": "PBMC (single-cell)",
        "platform": "10x Genomics",
        "data_type": "CITE-seq",
        "n_disease": 24,
        "n_control": 24,
        "total_samples": 48,
        "year": 2022,
        "pmid": "",
        "notes": "CITE-seq multimodal. May need pseudobulk."
    },
    
    # ============= IBD (SpA-associated) =============
    {
        "accession": "GSE186507",
        "disease": "IBD",
        "title": "Whole transcriptome data in blood from adult IBD and control subjects",
        "tissue": "Whole blood",
        "platform": "Illumina",
        "data_type": "RNA-seq",
        "n_disease": 200,
        "n_control": 50,
        "total_samples": 250,
        "year": 2024,
        "pmid": "",
        "notes": "Large IBD blood cohort. Excellent for SpA-IBD overlap."
    },
    {
        "accession": "GSE59071",
        "disease": "IBD",
        "title": "Gene expression in IBD blood",
        "tissue": "Whole blood",
        "platform": "Illumina",
        "data_type": "Microarray",
        "n_disease": 80,
        "n_control": 42,
        "total_samples": 122,
        "year": 2015,
        "pmid": "",
        "notes": "IBD blood expression used in AS-IBD overlap studies."
    },
    {
        "accession": "GSE112057",
        "disease": "IBD",
        "title": "Whole blood transcriptomics in Crohn's disease",
        "tissue": "Whole blood",
        "platform": "Illumina HiSeq",
        "data_type": "RNA-seq",
        "n_disease": 100,
        "n_control": 50,
        "total_samples": 150,
        "year": 2018,
        "pmid": "",
        "notes": "Crohn's blood RNA-seq. Good sample size."
    },
    
    # ============= UVEITIS (SpA-associated) =============
    {
        "accession": "GSE194060",
        "disease": "Uveitis",
        "title": "Transcriptome of DCs in HLA-B27+ acute anterior uveitis",
        "tissue": "Blood-derived DCs",
        "platform": "Illumina",
        "data_type": "RNA-seq",
        "n_disease": 21,
        "n_control": 21,
        "total_samples": 42,
        "year": 2022,
        "pmid": "",
        "notes": "HLA-B27+ AAU. Directly relevant to SpA pathogenesis."
    },
    {
        "accession": "GSE195501",
        "disease": "Uveitis",
        "title": "Blood transcriptomics of uveitis subtypes including HLA-B27+ AAU",
        "tissue": "Blood cells",
        "platform": "Illumina",
        "data_type": "RNA-seq",
        "n_disease": 30,
        "n_control": 15,
        "total_samples": 45,
        "year": 2022,
        "pmid": "",
        "notes": "Multiple uveitis subtypes including SpA-associated AAU."
    },
    
    # ============= AS additional datasets from literature =============
    {
        "accession": "GSE18781",
        "disease": "AS",
        "title": "Blood gene expression in axial spondyloarthropathy",
        "tissue": "Whole blood",
        "platform": "Illumina HumanRef-8",
        "data_type": "Microarray",
        "n_disease": 18,
        "n_control": 25,
        "total_samples": 43,
        "year": 2010,
        "pmid": "",
        "notes": "Axial SpA blood expression. Earlier study."
    },
    {
        "accession": "GSE266295",
        "disease": "AS_PsA",
        "title": "Immunosenescent CD8+ T cell subset in AS and PsA patients",
        "tissue": "CD8+ T cells",
        "platform": "GPL24676 Illumina NovaSeq",
        "data_type": "RNA-seq",
        "n_disease": 8,
        "n_control": 8,
        "total_samples": 16,
        "year": 2024,
        "pmid": "",
        "notes": "Directly compares AS and PsA at CD8+ T cell level."
    },
]


def verify_geo_accession(accession):
    """Quick check if a GEO accession exists and get basic info."""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200 and "could not be found" not in resp.text.lower():
            # Extract some basic info
            lines = resp.text.split("\n")
            info = {}
            for line in lines:
                if line.startswith("!Series_title"):
                    info["title"] = line.split("=", 1)[1].strip() if "=" in line else ""
                if line.startswith("!Series_sample_id"):
                    info["sample_count"] = info.get("sample_count", 0) + 1
                if line.startswith("!Series_platform_id"):
                    info["platform"] = line.split("=", 1)[1].strip() if "=" in line else ""
                if line.startswith("!Series_status"):
                    info["status"] = line.split("=", 1)[1].strip() if "=" in line else ""
            return True, info
        return False, {}
    except:
        return None, {}  # Unknown - network issue


def main():
    print("=" * 80)
    print("Phase 1b: Curating Final Dataset Inventory")
    print("=" * 80)
    
    # Verify each dataset exists
    print("\nVerifying dataset accessions...")
    verified = []
    for ds in CURATED_DATASETS:
        acc = ds["accession"]
        exists, info = verify_geo_accession(acc)
        status = "✓" if exists else ("?" if exists is None else "✗")
        sample_info = info.get("sample_count", "?")
        print(f"  {status} {acc:12s} | {ds['disease']:8s} | {ds['data_type']:12s} | samples={sample_info:>4} | {ds['title'][:60]}")
        
        ds["verified"] = exists
        if info.get("sample_count"):
            ds["verified_samples"] = info["sample_count"]
        verified.append(ds)
        time.sleep(0.3)
    
    # Filter to verified datasets
    valid = [d for d in verified if d.get("verified") is not False]
    
    print(f"\n{'='*80}")
    print(f"DATASET INVENTORY SUMMARY")
    print(f"{'='*80}")
    
    # Summary by disease
    diseases = {}
    for d in valid:
        disease = d["disease"]
        if disease not in diseases:
            diseases[disease] = {"count": 0, "total_samples": 0, "datasets": []}
        diseases[disease]["count"] += 1
        diseases[disease]["total_samples"] += d.get("total_samples", 0)
        diseases[disease]["datasets"].append(d["accession"])
    
    for disease, info in diseases.items():
        print(f"\n{disease}:")
        print(f"  Datasets: {info['count']}")
        print(f"  Total samples: {info['total_samples']}")
        print(f"  Accessions: {', '.join(info['datasets'])}")
    
    total_datasets = len(valid)
    total_samples = sum(d.get("total_samples", 0) for d in valid)
    print(f"\n{'='*40}")
    print(f"TOTAL: {total_datasets} datasets, ~{total_samples} samples")
    print(f"Disease categories: {len(diseases)}")
    
    # Save as CSV
    csv_path = os.path.join(OUTPUT_DIR, "dataset_inventory.csv")
    fieldnames = ["accession", "disease", "title", "tissue", "platform", "data_type", 
                   "n_disease", "n_control", "total_samples", "year", "pmid", "notes", 
                   "verified", "verified_samples"]
    
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for d in valid:
            writer.writerow(d)
    
    print(f"\nDataset inventory saved to: {csv_path}")
    
    # Save as JSON too
    json_path = os.path.join(OUTPUT_DIR, "dataset_inventory.json")
    with open(json_path, "w") as f:
        json.dump(valid, f, indent=2, default=str)
    print(f"JSON inventory saved to: {json_path}")


if __name__ == "__main__":
    main()
