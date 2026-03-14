"""
Configuration file for Cross-Disease Transcriptomic Meta-Analysis of Spondyloarthritis.

All hardcoded parameters, thresholds, dataset definitions, and API endpoints
are centralized here. Individual phase scripts should import from this module
rather than defining their own constants.

Usage:
    from config import CONFIG, DATASETS, DEG_PARAMS, META_PARAMS
"""

import os
from pathlib import Path

# ==============================================================================
# Project Directory Structure
# ==============================================================================

BASE_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = BASE_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
METADATA_DIR = DATA_DIR / "metadata"
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = BASE_DIR / "figures"
DOCS_DIR = BASE_DIR / "docs"
LOGS_DIR = BASE_DIR / "logs"

# Subdirectories for results
DEG_DIR = RESULTS_DIR / "deg"
META_DIR = RESULTS_DIR / "meta_analysis"
ENRICHMENT_DIR = RESULTS_DIR / "enrichment"
NETWORK_DIR = RESULTS_DIR / "networks"
WGCNA_DIR = RESULTS_DIR / "wgcna"
DECONV_DIR = RESULTS_DIR / "deconvolution"
DRUG_DIR = RESULTS_DIR / "drug_repurposing"
SENSITIVITY_DIR = RESULTS_DIR / "sensitivity"

# ==============================================================================
# GEO Dataset Definitions
# ==============================================================================

DATASETS = {
    "GSE25101": {
        "disease": "AS",
        "tissue": "Peripheral blood",
        "platform": "GPL6104",
        "platform_name": "Illumina HumanRef-8 v2.0",
        "technology": "Microarray",
        "description": "Ankylosing spondylitis whole blood expression profiling",
    },
    "GSE73754": {
        "disease": "AS",
        "tissue": "Peripheral blood",
        "platform": "GPL6104",
        "platform_name": "Illumina HumanRef-8 v2.0",
        "technology": "Microarray",
        "description": "Ankylosing spondylitis blood gene expression",
    },
    "GSE18781": {
        "disease": "AS",
        "tissue": "Peripheral blood",
        "platform": "GPL570",
        "platform_name": "Affymetrix HG-U133 Plus 2.0",
        "technology": "Microarray",
        "description": "Ankylosing spondylitis Affymetrix expression profiling",
    },
    "GSE221786": {
        "disease": "AS",
        "tissue": "Peripheral blood",
        "platform": "RNA-seq",
        "platform_name": "Illumina RNA-seq",
        "technology": "RNA-seq",
        "description": "Ankylosing spondylitis RNA-seq expression profiling",
    },
    "GSE61281": {
        "disease": "PsA",
        "tissue": "Peripheral blood",
        "platform": "GPL6480",
        "platform_name": "Agilent Whole Human Genome Microarray 4x44K",
        "technology": "Microarray",
        "description": "Psoriatic arthritis blood gene expression",
    },
    "GSE59071": {
        "disease": "IBD",
        "tissue": "Peripheral blood",
        "platform": "GPL570",
        "platform_name": "Affymetrix HG-U133 Plus 2.0",
        "technology": "Microarray",
        "description": "Inflammatory bowel disease blood expression profiling",
    },
    "GSE58667": {
        "disease": "jSpA",
        "tissue": "Peripheral blood",
        "platform": "GPL10558",
        "platform_name": "Illumina HumanHT-12 V4.0",
        "technology": "Microarray",
        "description": "Juvenile spondyloarthritis blood gene expression",
    },
}

# Ordered list of dataset accessions for consistent iteration
DATASET_ORDER = [
    "GSE25101", "GSE73754", "GSE18781", "GSE221786",
    "GSE61281", "GSE59071", "GSE58667",
]

# Disease categories
DISEASES = ["AS", "PsA", "IBD", "jSpA"]
AS_DATASETS = ["GSE25101", "GSE73754", "GSE18781", "GSE221786"]
CROSS_DISEASE_DATASETS = ["GSE61281", "GSE59071", "GSE58667"]

# ==============================================================================
# Differential Expression Parameters
# ==============================================================================

DEG_PARAMS = {
    "log2fc_threshold": 1.0,
    "padj_threshold": 0.05,
    "test_method": "welch",
    "correction_method": "fdr_bh",
    "min_samples_per_group": 3,
}

# ==============================================================================
# Meta-Analysis Parameters
# ==============================================================================

META_PARAMS = {
    "fdr_threshold": 0.05,
    "min_datasets": 2,
    "combination_method": "fisher",
    "vote_threshold": 0.5,
}

# ==============================================================================
# PPI Network Parameters (STRING)
# ==============================================================================

PPI_PARAMS = {
    "confidence_threshold": 0.7,
    "network_type": "physical",
    "species_taxid": 9606,
    "max_nodes": 500,
    "hub_percentile": 90,
}

# ==============================================================================
# WGCNA Parameters
# ==============================================================================

WGCNA_PARAMS = {
    "soft_threshold_power": 6,
    "min_module_size": 30,
    "merge_cut_height": 0.25,
    "deep_split": 2,
    "network_type": "signed",
    "tom_type": "signed",
    "reassign_threshold": 0.25,
    "n_top_genes": 5000,
}

# ==============================================================================
# Immune Deconvolution Parameters (ssGSEA)
# ==============================================================================

DECONV_PARAMS = {
    "method": "ssgsea",
    "n_cell_types": 13,
    "cell_types": [
        "B cells naive",
        "B cells memory",
        "T cells CD4 naive",
        "T cells CD4 memory",
        "T cells CD8",
        "T cells regulatory (Tregs)",
        "NK cells",
        "Monocytes",
        "Macrophages M0",
        "Macrophages M1",
        "Macrophages M2",
        "Dendritic cells",
        "Neutrophils",
    ],
    "permutations": 1000,
    "min_gene_set_size": 10,
}

# ==============================================================================
# Drug Repurposing Parameters (DGIdb)
# ==============================================================================

DRUG_PARAMS = {
    "min_interaction_score": 0.5,
    "interaction_types": [
        "inhibitor", "antagonist", "agonist", "modulator",
        "antibody", "activator", "blocker",
    ],
    "fda_approved_only": False,
    "max_genes_query": 50,
}

# ==============================================================================
# Sensitivity Analysis Parameters
# ==============================================================================

SENSITIVITY_PARAMS = {
    "fdr_thresholds": [0.01, 0.05, 0.10, 0.20],
    "log2fc_thresholds": [0.5, 1.0, 1.5, 2.0],
    "leave_one_out": True,
    "bootstrap_iterations": 1000,
    "random_seed": 42,
}

# ==============================================================================
# Figure and Visualization Parameters
# ==============================================================================

FIGURE_PARAMS = {
    "dpi": 300,
    "format": "png",
    "figsize_single": (8, 6),
    "figsize_wide": (14, 6),
    "figsize_tall": (8, 12),
    "figsize_large": (14, 12),
    "font_family": "Arial",
    "font_size": 12,
    "title_size": 14,
    "label_size": 11,
    "colormap": "RdBu_r",
    "palette": "Set2",
}

# ==============================================================================
# External API Endpoints
# ==============================================================================

API_ENDPOINTS = {
    "ncbi_esearch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "ncbi_esummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    "ncbi_efetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "geo_query": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
    "string_api": "https://string-db.org/api",
    "string_network": "https://string-db.org/api/json/network",
    "string_enrichment": "https://string-db.org/api/json/enrichment",
    "string_interaction_partners": "https://string-db.org/api/json/interaction_partners",
    "gprofiler": "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
    "gprofiler_convert": "https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
    "gprofiler_orth": "https://biit.cs.ut.ee/gprofiler/api/orth/orth/",
    "dgidb": "https://dgidb.org/api/v2/interactions.json",
    "dgidb_genes": "https://dgidb.org/api/v2/genes.json",
    "dgidb_drugs": "https://dgidb.org/api/v2/drugs.json",
}

API_SETTINGS = {
    "timeout": 30,
    "max_retries": 3,
    "retry_delay": 5,
    "rate_limit_delay": 0.34,
    "ncbi_api_key": None,
    "user_agent": "SpA_MetaAnalysis/1.0 (Cross-Disease Transcriptomic Study)",
}


def ensure_dirs():
    """Create all necessary output directories if they don't exist."""
    dirs = [
        DATA_DIR, RAW_DIR, PROCESSED_DIR, METADATA_DIR,
        RESULTS_DIR, DEG_DIR, META_DIR, ENRICHMENT_DIR,
        NETWORK_DIR, WGCNA_DIR, DECONV_DIR, DRUG_DIR, SENSITIVITY_DIR,
        FIGURES_DIR, DOCS_DIR, LOGS_DIR,
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def get_dataset_info(accession):
    """Return dataset metadata dictionary for a given GEO accession."""
    if accession not in DATASETS:
        raise ValueError(f"Unknown dataset: {accession}. Valid: {list(DATASETS.keys())}")
    return DATASETS[accession]


def get_datasets_by_disease(disease):
    """Return list of GEO accessions for a given disease category."""
    return [acc for acc, info in DATASETS.items() if info["disease"] == disease]


if __name__ == "__main__":
    print("SpA Cross-Disease Meta-Analysis Configuration")
    print("=" * 55)
    print(f"Base directory:    {BASE_DIR}")
    print(f"Data directory:    {DATA_DIR}")
    print(f"Results directory: {RESULTS_DIR}")
    print(f"Figures directory: {FIGURES_DIR}")
