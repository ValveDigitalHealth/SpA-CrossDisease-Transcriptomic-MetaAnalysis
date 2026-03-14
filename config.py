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
    "log2fc_threshold": 1.0,          # |log2FC| cutoff for significance
    "padj_threshold": 0.05,           # BH-adjusted p-value cutoff
    "test_method": "welch",           # Welch's t-test (unequal variance)
    "correction_method": "fdr_bh",    # Benjamini-Hochberg FDR correction
    "min_samples_per_group": 3,       # Minimum samples in each group
}

# ==============================================================================
# Meta-Analysis Parameters
# ==============================================================================

META_PARAMS = {
    "fdr_threshold": 0.05,            # FDR threshold for meta-significance
    "min_datasets": 2,                # Minimum datasets with concordant direction
    "combination_method": "fisher",   # Fisher's combined probability test
    "vote_threshold": 0.5,            # Fraction of datasets for vote counting
}

# ==============================================================================
# PPI Network Parameters (STRING)
# ==============================================================================

PPI_PARAMS = {
    "confidence_threshold": 0.7,      # STRING combined score >= 700 (high confidence)
    "network_type": "physical",       # Type of interactions
    "species_taxid": 9606,            # Homo sapiens NCBI taxonomy ID
    "max_nodes": 500,                 # Maximum network nodes for visualization
    "hub_percentile": 90,             # Percentile for hub gene identification
}

# ==============================================================================
# WGCNA Parameters
# ==============================================================================

WGCNA_PARAMS = {
    "soft_threshold_power": 6,        # Soft-thresholding power (scale-free topology)
    "min_module_size": 30,            # Minimum genes per module
    "merge_cut_height": 0.25,         # Module merge distance threshold
    "deep_split": 2,                  # Tree cutting sensitivity (0-4)
    "network_type": "signed",         # Signed network for direction preservation
    "tom_type": "signed",             # Topological overlap matrix type
    "reassign_threshold": 0.25,       # Module reassignment threshold
    "n_top_genes": 5000,              # Top variable genes for WGCNA input
}

# ==============================================================================
# Immune Deconvolution Parameters (ssGSEA)
# ==============================================================================

DECONV_PARAMS = {
    "method": "ssgsea",               # Single-sample GSEA
    "n_cell_types": 13,               # Number of immune cell types
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
    "permutations": 1000,             # Permutations for enrichment scoring
    "min_gene_set_size": 10,          # Minimum genes per cell type signature
}

# ==============================================================================
# Drug Repurposing Parameters (DGIdb)
# ==============================================================================

DRUG_PARAMS = {
    "min_interaction_score": 0.5,     # Minimum drug-gene interaction score
    "interaction_types": [
        "inhibitor", "antagonist", "agonist", "modulator",
        "antibody", "activator", "blocker",
    ],
    "fda_approved_only": False,       # Include non-approved candidates
    "max_genes_query": 50,            # Maximum genes per DGIdb query
}

# ==============================================================================
# Sensitivity Analysis Parameters
# ==============================================================================

SENSITIVITY_PARAMS = {
    "fdr_thresholds": [0.01, 0.05, 0.10, 0.20],
    "log2fc_thresholds": [0.5, 1.0, 1.5, 2.0],
    "leave_one_out": True,            # Run leave-one-dataset-out analysis
    "bootstrap_iterations": 1000,     # Bootstrap resampling iterations
    "random_seed": 42,
}

# ==============================================================================
# Figure and Visualization Parameters
# ==============================================================================

FIGURE_PARAMS = {
    "dpi": 300,                       # Publication quality
    "format": "png",                  # Default output format
    "figsize_single": (8, 6),         # Single panel figure
    "figsize_wide": (14, 6),          # Wide figure
    "figsize_tall": (8, 12),          # Tall figure
    "figsize_large": (14, 12),        # Large multi-panel figure
    "font_family": "Arial",
    "font_size": 12,
    "title_size": 14,
    "label_size": 11,
    "colormap": "RdBu_r",            # Diverging colormap for heatmaps
    "palette": "Set2",                # Categorical palette
}

# ==============================================================================
# External API Endpoints
# ==============================================================================

API_ENDPOINTS = {
    # NCBI E-utilities
    "ncbi_esearch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "ncbi_esummary": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    "ncbi_efetch": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "geo_query": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi",

    # STRING API (protein-protein interactions)
    "string_api": "https://string-db.org/api",
    "string_network": "https://string-db.org/api/json/network",
    "string_enrichment": "https://string-db.org/api/json/enrichment",
    "string_interaction_partners": "https://string-db.org/api/json/interaction_partners",

    # g:Profiler (functional enrichment)
    "gprofiler": "https://biit.cs.ut.ee/gprofiler/api/gost/profile/",
    "gprofiler_convert": "https://biit.cs.ut.ee/gprofiler/api/convert/convert/",
    "gprofiler_orth": "https://biit.cs.ut.ee/gprofiler/api/orth/orth/",

    # DGIdb (drug-gene interactions)
    "dgidb": "https://dgidb.org/api/v2/interactions.json",
    "dgidb_genes": "https://dgidb.org/api/v2/genes.json",
    "dgidb_drugs": "https://dgidb.org/api/v2/drugs.json",
}

# API request settings
API_SETTINGS = {
    "timeout": 30,                    # Request timeout in seconds
    "max_retries": 3,                 # Maximum retry attempts
    "retry_delay": 5,                 # Seconds between retries
    "rate_limit_delay": 0.34,         # NCBI: max 3 requests/second
    "ncbi_api_key": None,             # Set via NCBI_API_KEY environment variable
    "user_agent": "SpA_MetaAnalysis/1.0 (Cross-Disease Transcriptomic Study)",
}

# ==============================================================================
# Pipeline Phase Definitions
# ==============================================================================

PIPELINE_PHASES = [
    {
        "phase": "1",
        "name": "Dataset Discovery",
        "script": "phase1_dataset_discovery.py",
        "description": "Query NCBI GEO for SpA-related transcriptomic datasets",
    },
    {
        "phase": "1b",
        "name": "Dataset Curation",
        "script": "phase1b_curate_datasets.py",
        "description": "Curate candidates to final dataset inventory",
    },
    {
        "phase": "2",
        "name": "Download & Preprocessing",
        "script": "phase2_download_preprocess.py",
        "description": "Download GEO data and normalize expression matrices",
    },
    {
        "phase": "2b",
        "name": "Fix Groups & RNA-seq Download",
        "script": "phase2b_fix_groups_download_rnaseq.py",
        "description": "Fix sample group assignments and download RNA-seq data",
    },
    {
        "phase": "2c",
        "name": "Process RNA-seq",
        "script": "phase2c_process_rnaseq.py",
        "description": "Process and normalize RNA-seq expression data",
    },
    {
        "phase": "3",
        "name": "Differential Expression",
        "script": "phase3_differential_expression.py",
        "description": "Run Welch's t-test with BH correction per dataset",
    },
    {
        "phase": "3b",
        "name": "Probe Mapping",
        "script": "phase3b_probe_mapping.py",
        "description": "Map microarray probes to gene symbols",
    },
    {
        "phase": "3b_quick",
        "name": "Quick Probe Mapping",
        "script": "phase3b_quick_mapping.py",
        "description": "Rapid probe-to-gene mapping for validation",
    },
    {
        "phase": "4",
        "name": "Meta-Analysis",
        "script": "phase4_meta_analysis.py",
        "description": "Fisher's combined probability test across datasets",
    },
    {
        "phase": "5",
        "name": "Enrichment & Figures",
        "script": "phase5_enrichment_and_figures.py",
        "description": "GO/Reactome enrichment via g:Profiler and figure generation",
    },
    {
        "phase": "6",
        "name": "PPI Network",
        "script": "phase6_ppi_network.py",
        "description": "STRING protein-protein interaction network analysis",
    },
    {
        "phase": "7",
        "name": "WGCNA",
        "script": "phase7_wgcna.py",
        "description": "Weighted gene co-expression network analysis",
    },
    {
        "phase": "8",
        "name": "Immune Deconvolution",
        "script": "phase8_immune_deconvolution.py",
        "description": "ssGSEA-based immune cell type deconvolution",
    },
    {
        "phase": "8b",
        "name": "Fix PsA Deconvolution",
        "script": "phase8b_fix_psa_deconvolution.py",
        "description": "Correct PsA deconvolution results",
    },
    {
        "phase": "8c",
        "name": "Update Cross-Disease",
        "script": "phase8c_update_crossdisease.py",
        "description": "Update cross-disease comparison tables and figures",
    },
    {
        "phase": "sensitivity",
        "name": "Sensitivity Analysis",
        "script": "sensitivity_analysis.py",
        "description": "Robustness analysis across thresholds and configurations",
    },
    {
        "phase": "supplementary",
        "name": "Supplementary Tables",
        "script": "generate_supplementary_tables.py",
        "description": "Generate supplementary Excel tables for manuscript",
    },
    {
        "phase": "abstract",
        "name": "Graphical Abstract",
        "script": "graphical_abstract.py",
        "description": "Generate graphical abstract for manuscript",
    },
    {
        "phase": "figures",
        "name": "Inject Figures",
        "script": "inject_figures.py",
        "description": "Inject figures into manuscript document",
    },
    {
        "phase": "figures_v2",
        "name": "Inject Figures v2",
        "script": "inject_figures_v2.py",
        "description": "Updated figure injection for revised manuscript",
    },
]

# ==============================================================================
# Utility Functions
# ==============================================================================

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
    print()
    print(f"Datasets ({len(DATASETS)}):")
    for acc in DATASET_ORDER:
        info = DATASETS[acc]
        print(f"  {acc}: {info['disease']:4s} | {info['platform_name']}")
    print()
    print(f"Pipeline phases: {len(PIPELINE_PHASES)}")
    print(f"DEG thresholds:  |log2FC| >= {DEG_PARAMS['log2fc_threshold']}, "
          f"padj < {DEG_PARAMS['padj_threshold']}")
    print(f"Meta-analysis:   FDR < {META_PARAMS['fdr_threshold']}, "
          f"min datasets = {META_PARAMS['min_datasets']}")
    print(f"PPI confidence:  >= {PPI_PARAMS['confidence_threshold']}")
