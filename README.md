# Cross-Disease Transcriptomic Meta-Analysis of Spondyloarthritis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Datasets: 7 GEO](https://img.shields.io/badge/datasets-7%20GEO-green.svg)](https://www.ncbi.nlm.nih.gov/geo/)
[![Samples: 370](https://img.shields.io/badge/samples-370-orange.svg)]()

> First cross-disease transcriptomic meta-analysis spanning the full spondyloarthritis spectrum — from ankylosing spondylitis to psoriatic arthritis, inflammatory bowel disease, and juvenile spondyloarthritis.

---

## Citation

Kuklik P. Cross-Disease Transcriptomic Meta-Analysis Identifies RNA Processing Dysregulation and Drug Repurposing Candidates for Ankylosing Spondylitis. *Spondyloarthritis AI Computational Biology Institute, Lisbon, Portugal.* Submitted to *Annals of the Rheumatic Diseases*, 2026.

---

## Overview

This repository contains a fully reproducible computational pipeline for cross-disease transcriptomic meta-analysis of spondyloarthritis (SpA). The study integrates **7 Gene Expression Omnibus (GEO) datasets** encompassing **370 peripheral blood samples** across four SpA-related diseases: ankylosing spondylitis (AS), psoriatic arthritis (PsA), inflammatory bowel disease (IBD), and juvenile spondyloarthritis (jSpA).

The pipeline implements a 10-phase workflow spanning dataset discovery, differential expression analysis, Fisher's meta-analysis, functional enrichment, protein-protein interaction network construction, weighted gene co-expression network analysis (WGCNA), immune cell deconvolution, drug repurposing, and comprehensive sensitivity analyses.

All computational analyses were performed by Perplexity Computer under the direction of Pawel Kuklik at the Spondyloarthritis AI Computational Biology Institute (Lisbon, Portugal).

---

## Key Findings

- **40 consensus AS meta-significant genes** identified (35 downregulated, 5 upregulated) through Fisher's combined probability test across 4 AS datasets
- **RNA processing and splicing pathway enrichment** — Gene Ontology and Reactome analyses converge on mRNA splicing, RNA processing, and spliceosome-related pathways as the central dysregulated programs in AS
- **Hub genes**: PTBP1, SRSF10, EXOSC10, SNIP1, GFM1 — splicing regulators and RNA processing factors dominate the PPI network
- **Drug repurposing via DGIdb**: DCK (deoxycytidine kinase) yielded 23 drug-gene interactions including 7 FDA-approved drugs; cladribine emerges as a top repurposing candidate given its established use in multiple sclerosis
- **IBD shows massive immune activation**: 11 of 13 immune cell types reach FDR significance in IBD, compared to selective perturbation in AS, PsA, and jSpA
- **Limited WGCNA module preservation** across diseases, indicating distinct molecular programs despite shared clinical SpA classification
- **Sensitivity analyses** confirm robustness of core findings across FDR thresholds, effect size cutoffs, and leave-one-out dataset removal

---

## Datasets

| GEO Accession | Disease | Tissue | Platform | Technology | Samples (Disease/Control) |
|:---|:---|:---|:---|:---|:---|
| [GSE25101](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25101) | AS | Peripheral blood | Illumina HumanRef-8 v2.0 | Microarray | 16 / 16 |
| [GSE73754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73754) | AS | Peripheral blood | Illumina HumanRef-8 v2.0 | Microarray | 52 / 20 |
| [GSE18781](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18781) | AS | Peripheral blood | Affymetrix HG-U133 Plus 2.0 | Microarray | 18 / 18 |
| [GSE221786](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221786) | AS | Peripheral blood | Illumina RNA-seq | RNA-seq | 30 / 30 |
| [GSE61281](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61281) | PsA | Peripheral blood | Agilent Whole Genome 4x44K | Microarray | 20 / 20 |
| [GSE59071](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59071) | IBD | Peripheral blood | Affymetrix HG-U133 Plus 2.0 | Microarray | 60 / 40 |
| [GSE58667](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58667) | jSpA | Peripheral blood | Illumina HumanHT-12 V4.0 | Microarray | 30 / 20 |

**Total**: 7 datasets, 4 diseases, ~370 samples, 4 microarray platforms + 1 RNA-seq platform

---

## Pipeline Architecture

The analysis pipeline consists of 11 logical phases executed across 20 Python scripts:

```
Phase 0   Environment Setup & Dependency Installation
    │
Phase 1   Dataset Discovery & Curation
    │      ├── phase1_dataset_discovery.py    (NCBI E-utilities → GEO search)
    │      └── phase1b_curate_datasets.py     (21 candidates → 7 final datasets)
    │
Phase 2   Data Download & Preprocessing
    │      ├── phase2_download_preprocess.py          (GEO matrix download + normalization)
    │      ├── phase2b_fix_groups_download_rnaseq.py  (Sample group fixes + RNA-seq)
    │      └── phase2c_process_rnaseq.py              (RNA-seq processing pipeline)
    │
Phase 3   Differential Expression Analysis
    │      ├── phase3_differential_expression.py      (Welch's t-test + BH correction)
    │      ├── phase3b_probe_mapping.py               (Probe → gene symbol mapping)
    │      └── phase3b_quick_mapping.py               (Rapid mapping for validation)
    │
Phase 4   Fisher's Meta-Analysis
    │      └── phase4_meta_analysis.py                (Combined p-values + vote counting)
    │
Phase 5   Functional Enrichment & Visualization
    │      └── phase5_enrichment_and_figures.py       (g:Profiler → GO, Reactome)
    │
Phase 6   PPI Network Analysis
    │      └── phase6_ppi_network.py                  (STRING → NetworkX hub detection)
    │
Phase 7   WGCNA Co-expression Network
    │      └── phase7_wgcna.py                        (Module detection + preservation)
    │
Phase 8   Immune Cell Deconvolution
    │      ├── phase8_immune_deconvolution.py         (ssGSEA, 13 cell types)
    │      ├── phase8b_fix_psa_deconvolution.py       (PsA result correction)
    │      └── phase8c_update_crossdisease.py         (Cross-disease comparisons)
    │
Phase 9   Drug-Gene Interaction Analysis
    │      (integrated within enrichment and network scripts via DGIdb API)
    │
Phase 10  Sensitivity Analysis
    │      └── sensitivity_analysis.py                (FDR thresholds, leave-one-out)
    │
Phase 11  Figures, Manuscript & Submission
           ├── graphical_abstract.py
           ├── generate_supplementary_tables.py
           ├── inject_figures.py
           └── inject_figures_v2.py
```

---

## Quick Start

### Prerequisites

- Python 3.9 or higher
- Internet connection (for GEO data download and API queries)
- ~5 GB free disk space (for raw data and results)

### Installation

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/SpA_CrossDisease_MetaAnalysis.git
cd SpA_CrossDisease_MetaAnalysis

# Create virtual environment (recommended)
python -m venv .venv
source .venv/bin/activate   # Linux/macOS
# .venv\Scripts\activate    # Windows

# Install dependencies
pip install -r requirements.txt
```

### Run the Full Pipeline

```bash
python run_pipeline.py --all
```

### Run Individual Phases

```bash
# Run a specific phase
python run_pipeline.py --phase 3

# Run multiple specific phases
python run_pipeline.py --phase 1 1b 2

# Run from a specific phase onwards
python run_pipeline.py --from-phase 4

# List all available phases
python run_pipeline.py --list

# Run with logging enabled
python run_pipeline.py --all --log

# Stop on first error
python run_pipeline.py --all --stop-on-error
```

---

## Step-by-Step Reproduction Guide

Below is a detailed walkthrough of each pipeline phase. Someone with basic Python knowledge should be able to follow these instructions to fully reproduce the analysis.

### Phase 0: Environment Setup

Before running any analysis, set up the Python environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Verify the configuration:

```bash
python config.py
```

This prints the project directory structure, dataset inventory, and key parameters. Confirm that 7 datasets are listed and all paths resolve correctly.

---

### Phase 1: Dataset Discovery

```bash
python scripts/phase1_dataset_discovery.py
python scripts/phase1b_curate_datasets.py
```

**What it does**: Queries the NCBI GEO database via E-utilities for spondyloarthritis-related transcriptomic datasets. The discovery script identifies ~21 candidate datasets. The curation script applies inclusion/exclusion criteria (peripheral blood samples, case-control design, minimum sample size) to select the final 7 datasets.

**Output files**:
- `data/metadata/dataset_inventory.csv` — Tabular inventory of selected datasets
- `data/metadata/dataset_inventory.json` — Structured metadata for programmatic access

**How to verify**: Open `dataset_inventory.csv` and confirm that exactly 7 datasets are listed with correct GEO accessions (GSE25101, GSE73754, GSE18781, GSE221786, GSE61281, GSE59071, GSE58667), disease labels, and platform annotations.

---

### Phase 2: Data Download & Preprocessing

```bash
python scripts/phase2_download_preprocess.py
python scripts/phase2b_fix_groups_download_rnaseq.py
python scripts/phase2c_process_rnaseq.py
```

**What it does**:
1. **phase2**: Downloads GEO series matrix files for all 6 microarray datasets. Applies log2 transformation (if needed) and quantile normalization. Extracts sample metadata and assigns disease/control group labels.
2. **phase2b**: Corrects sample group assignments where GEO metadata is ambiguous. Downloads the RNA-seq dataset (GSE221786) which requires separate handling.
3. **phase2c**: Processes RNA-seq count data — applies TMM normalization and log2-CPM transformation to make RNA-seq data comparable to microarray datasets.

**Output files**:
- `data/raw/` — Downloaded series matrix files (one per dataset)
- `data/processed/` — Normalized expression matrices (CSV, one per dataset)
- `data/metadata/` — Updated sample metadata with group assignments

**How to verify**:
- Check that `data/raw/` contains downloaded files for each GEO accession
- Open any processed expression matrix and confirm values are log2-scale (typically range 2–16 for microarray)
- Verify sample counts match the expected numbers in the Datasets table above

**Notes**: Downloads require internet access. GEO servers may be slow; allow up to 30 minutes. If a download fails, re-running the script will retry failed datasets.

---

### Phase 3: Differential Expression Analysis

```bash
python scripts/phase3_differential_expression.py
python scripts/phase3b_probe_mapping.py
python scripts/phase3b_quick_mapping.py
```

**What it does**:
1. **phase3**: Performs per-dataset differential expression analysis using Welch's t-test (unequal variance assumption) for each gene/probe. Applies Benjamini-Hochberg FDR correction. Uses thresholds of |log2FC| ≥ 1.0 and adjusted p-value < 0.05.
2. **phase3b**: Maps microarray probe IDs to HGNC gene symbols using platform annotation files. Resolves many-to-one mappings by selecting the probe with lowest p-value per gene.
3. **phase3b_quick**: A rapid mapping variant for iterative validation and QC.

**Output files**:
- `results/deg/` — Per-dataset DEG result tables (CSV) with columns: gene, log2FC, p-value, adjusted p-value, t-statistic
- `results/deg/deg_summary.csv` — Summary counts of up/down-regulated genes per dataset

**How to verify**:
- Each dataset should produce a DEG results file in `results/deg/`
- Check that the summary shows AS datasets have moderate DEG counts and IBD shows a higher count
- Spot-check that significant genes have |log2FC| ≥ 1.0 and padj < 0.05

---

### Phase 4: Fisher's Meta-Analysis

```bash
python scripts/phase4_meta_analysis.py
```

**What it does**: Combines per-dataset p-values using Fisher's combined probability test across the 4 AS datasets (GSE25101, GSE73754, GSE18781, GSE221786). Also performs vote counting to ensure directional concordance — a gene must show consistent up- or down-regulation in at least 2 AS datasets. Applies FDR correction (threshold 0.05) to the combined p-values.

**Output files**:
- `results/meta_analysis/fisher_meta_results.csv` — Full meta-analysis results for all genes
- `results/meta_analysis/meta_significant_genes.csv` — The 40 consensus meta-significant genes
- `results/meta_analysis/meta_upregulated.csv` — 5 upregulated meta-significant genes
- `results/meta_analysis/meta_downregulated.csv` — 35 downregulated meta-significant genes

**How to verify**:
- Confirm `meta_significant_genes.csv` contains 40 genes
- Check that 35 are downregulated and 5 are upregulated
- Verify that PTBP1, SRSF10, EXOSC10, SNIP1, GFM1 appear in the results

---

### Phase 5: Functional Enrichment & Figures

```bash
python scripts/phase5_enrichment_and_figures.py
```

**What it does**: Submits the 40 meta-significant genes to g:Profiler for functional enrichment analysis against Gene Ontology (Biological Process, Molecular Function, Cellular Component), Reactome pathways, and KEGG pathways. Generates publication-quality figures including volcano plots, enrichment dot plots, and heatmaps.

**Output files**:
- `results/enrichment/` — Enrichment result tables (CSV) with term ID, name, p-value, gene list
- `figures/` — Volcano plots, enrichment bar/dot plots, heatmaps (PNG, 300 dpi)

**How to verify**:
- Enrichment results should show RNA processing/splicing terms at the top
- Check for GO terms like "mRNA processing", "RNA splicing", "spliceosome"
- Figures should render correctly with readable labels

---

### Phase 6: PPI Network Analysis

```bash
python scripts/phase6_ppi_network.py
```

**What it does**: Queries the STRING database (v12.0) for protein-protein interactions among the 40 meta-significant genes at high confidence (combined score ≥ 0.7). Constructs a NetworkX graph object and computes network topology metrics (degree centrality, betweenness centrality, clustering coefficient). Identifies hub genes based on degree centrality.

**Output files**:
- `results/networks/ppi_network.graphml` — Network in GraphML format
- `results/networks/hub_genes.csv` — Ranked hub genes with centrality scores
- `results/networks/network_stats.json` — Network-level statistics
- `figures/ppi_network.png` — Network visualization

**How to verify**:
- Hub genes should include PTBP1, SRSF10, EXOSC10, SNIP1, GFM1
- Network should show a connected component around RNA processing genes
- STRING confidence scores should all be ≥ 0.7

---

### Phase 7: WGCNA Co-expression Network

```bash
python scripts/phase7_wgcna.py
```

**What it does**: Performs Weighted Gene Co-expression Network Analysis on the top 5,000 most variable genes within each dataset. Identifies co-expression modules using dynamic tree cutting. Tests module preservation across diseases using permutation-based statistics. Correlates modules with disease status.

**Output files**:
- `results/wgcna/` — Module assignments, eigengene values, preservation statistics
- `figures/wgcna_*.png` — Module dendrograms, eigengene heatmaps, preservation plots

**How to verify**:
- Each dataset should produce identifiable co-expression modules (colored dendrogram)
- Module preservation across diseases should be limited (Zsummary < 10 for most cross-disease comparisons)
- At least some modules should correlate significantly with disease status

---

### Phase 8: Immune Cell Deconvolution

```bash
python scripts/phase8_immune_deconvolution.py
python scripts/phase8b_fix_psa_deconvolution.py
python scripts/phase8c_update_crossdisease.py
```

**What it does**:
1. **phase8**: Applies single-sample Gene Set Enrichment Analysis (ssGSEA) to estimate the relative abundance of 13 immune cell types in each sample. Compares cell type scores between disease and control groups using Wilcoxon rank-sum tests with FDR correction.
2. **phase8b**: Corrects specific issues in PsA deconvolution results.
3. **phase8c**: Generates cross-disease comparison tables and figures showing immune cell perturbation patterns across AS, PsA, IBD, and jSpA.

**Cell types analyzed** (13 total):
B cells naive, B cells memory, T cells CD4 naive, T cells CD4 memory, T cells CD8, T cells regulatory (Tregs), NK cells, Monocytes, Macrophages M0, Macrophages M1, Macrophages M2, Dendritic cells, Neutrophils

**Output files**:
- `results/deconvolution/` — Per-dataset ssGSEA scores, statistical comparisons, cross-disease summary
- `figures/deconvolution_*.png` — Box plots, heatmaps of immune cell proportions

**How to verify**:
- IBD should show the broadest immune activation (11/13 cell types FDR-significant)
- AS datasets should show more selective perturbation
- Cross-disease comparison should reveal distinct immune signatures per disease

---

### Phase 9: Drug Repurposing (via DGIdb)

Drug-gene interaction analysis is integrated within the enrichment and network analysis scripts. It queries the DGIdb (Drug-Gene Interaction Database) API for druggable targets among the meta-significant genes.

**Key result**: DCK (deoxycytidine kinase) yields 23 drug-gene interactions including 7 FDA-approved drugs. Cladribine, an FDA-approved drug for multiple sclerosis, emerges as a top repurposing candidate.

**Output files**:
- `results/drug_repurposing/` — Drug-gene interaction tables, druggability scores

---

### Phase 10: Sensitivity Analysis

```bash
python scripts/sensitivity_analysis.py
```

**What it does**: Tests the robustness of the meta-analysis results across multiple dimensions:
1. **FDR threshold sensitivity**: Varies the FDR cutoff (0.01, 0.05, 0.10, 0.20) and counts recovered meta-significant genes
2. **Effect size sensitivity**: Varies the log2FC threshold (0.5, 1.0, 1.5, 2.0) and measures impact on gene lists
3. **Leave-one-out analysis**: Removes each AS dataset in turn and re-runs Fisher's meta-analysis to identify genes robust to single-dataset removal
4. **Bootstrap resampling**: 1,000 iterations for confidence interval estimation

**Output files**:
- `results/sensitivity/` — Robustness tables, leave-one-out results, threshold sweep data
- `figures/sensitivity_*.png` — Threshold sensitivity curves, overlap heatmaps

**How to verify**:
- Core findings (RNA processing pathway, hub genes) should persist across reasonable threshold ranges
- Leave-one-out should show >70% gene retention when removing any single dataset

---

### Phase 11: Figures & Manuscript Support

```bash
python scripts/graphical_abstract.py
python scripts/generate_supplementary_tables.py
python scripts/inject_figures.py
python scripts/inject_figures_v2.py
```

**What it does**:
1. **graphical_abstract.py**: Generates the graphical abstract summarizing the study design and key findings
2. **generate_supplementary_tables.py**: Compiles all results into formatted Excel supplementary tables for journal submission
3. **inject_figures.py / inject_figures_v2.py**: Embeds generated figures into the manuscript document at designated placeholders

**Output files**:
- `figures/graphical_abstract.png` — Publication-ready graphical abstract
- `docs/` — Supplementary tables (XLSX), updated manuscript with figures

---

## Output Structure

After running the full pipeline, the repository will contain:

```
data/
├── raw/                    # Downloaded GEO series matrix files
│   ├── GSE25101/           # Raw data per dataset
│   ├── GSE73754/
│   ├── GSE18781/
│   ├── GSE221786/
│   ├── GSE61281/
│   ├── GSE59071/
│   └── GSE58667/
├── processed/              # Normalized expression matrices (CSV)
│   ├── GSE25101_normalized.csv
│   ├── ...
│   └── GSE58667_normalized.csv
└── metadata/               # Dataset inventories and sample metadata
    ├── dataset_inventory.csv
    └── dataset_inventory.json

results/
├── deg/                    # Per-dataset DEG results
│   ├── GSE25101_deg.csv
│   ├── ...
│   └── deg_summary.csv
├── meta_analysis/          # Fisher's combined meta-analysis
│   ├── fisher_meta_results.csv
│   ├── meta_significant_genes.csv
│   ├── meta_upregulated.csv
│   └── meta_downregulated.csv
├── enrichment/             # GO/Reactome enrichment results
├── networks/               # PPI network, hub genes
│   ├── ppi_network.graphml
│   ├── hub_genes.csv
│   └── network_stats.json
├── wgcna/                  # Co-expression modules and preservation
├── deconvolution/          # ssGSEA immune cell scores
├── drug_repurposing/       # DGIdb interaction tables
└── sensitivity/            # Robustness analysis results

figures/                    # All publication-quality figures (PNG, 300 dpi)

docs/                       # Manuscript, cover letter, supplementary tables
```

---

## Dependencies

| Package | Version | Purpose |
|:---|:---|:---|
| pandas | ≥ 1.5.0 | Data manipulation, CSV/Excel I/O |
| numpy | ≥ 1.24.0 | Numerical computing, array operations |
| scipy | ≥ 1.10.0 | Statistical tests (t-test, Fisher's method, Wilcoxon) |
| statsmodels | ≥ 0.14.0 | Multiple testing correction (BH FDR) |
| matplotlib | ≥ 3.7.0 | Publication-quality figure generation |
| seaborn | ≥ 0.12.0 | Statistical visualization (heatmaps, box plots) |
| networkx | ≥ 3.0 | PPI network construction and analysis |
| requests | ≥ 2.28.0 | HTTP requests to external APIs |
| openpyxl | ≥ 3.1.0 | Excel file generation for supplementary tables |

All dependencies are specified in `requirements.txt` and can be installed via:

```bash
pip install -r requirements.txt
```

---

## External APIs Used

The pipeline queries four external APIs during execution. All APIs are free for academic use and do not require authentication keys (though an NCBI API key is recommended for faster GEO queries).

| API | URL | Purpose | Phase |
|:---|:---|:---|:---|
| NCBI E-utilities | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ | GEO dataset search and series matrix download | 1, 2 |
| STRING | https://string-db.org/api | Protein-protein interaction network retrieval | 6 |
| g:Profiler | https://biit.cs.ut.ee/gprofiler/ | GO and Reactome functional enrichment | 5 |
| DGIdb | https://dgidb.org/api/ | Drug-gene interaction queries | 5, 6 |

**Rate limits**: NCBI limits to 3 requests/second without an API key (10/second with a key). The pipeline includes built-in delays to respect these limits. STRING and g:Profiler have generous rate limits for typical academic queries.

**Optional**: Set an NCBI API key to increase rate limits:

```bash
export NCBI_API_KEY="your_key_here"
```

---

## Methodological Notes

### Differential Expression
- **Test**: Welch's t-test (does not assume equal variances between groups)
- **Correction**: Benjamini-Hochberg FDR (controls false discovery rate at 5%)
- **Thresholds**: |log2FC| ≥ 1.0 and adjusted p-value < 0.05
- **Rationale**: Welch's t-test is robust for small, unequal sample sizes typical of GEO datasets. The log2FC threshold of 1.0 (2-fold change) ensures biological meaningfulness.

### Meta-Analysis
- **Method**: Fisher's combined probability test — combines independent p-values from 4 AS datasets into a single test statistic following a chi-squared distribution
- **Vote counting**: Requires concordant direction (up or down) in ≥ 2 datasets to prevent conflicting signals
- **FDR**: Applied to combined p-values to control the meta-level false discovery rate

### Immune Deconvolution
- **Method**: Single-sample Gene Set Enrichment Analysis (ssGSEA) — a rank-based method that scores each sample against reference gene signatures for 13 immune cell types
- **Comparison**: Wilcoxon rank-sum test (non-parametric) comparing disease vs. control scores per cell type per dataset
- **Advantage**: ssGSEA works on bulk expression data without requiring a reference mixture, making it applicable across different microarray platforms

### WGCNA
- **Network construction**: Signed weighted correlation network using soft-thresholding to achieve approximate scale-free topology
- **Module detection**: Dynamic tree cutting on the topological overlap matrix (TOM) dendrogram
- **Preservation**: Zsummary statistic from permutation-based module preservation analysis — values > 10 indicate strong preservation, < 2 indicate no preservation

### Drug Repurposing
- **Source**: Drug-Gene Interaction Database (DGIdb) — curated from over 30 source databases
- **Approach**: Query meta-significant genes, especially hub genes, for known drug interactions
- **Validation**: Cross-reference with FDA approval status and existing clinical indications

---

## Configuration

All analysis parameters are centralized in `config.py`. Key configuration sections:

```python
from config import (
    DATASETS,           # GEO dataset definitions
    DEG_PARAMS,         # Differential expression thresholds
    META_PARAMS,        # Meta-analysis settings
    PPI_PARAMS,         # STRING network parameters
    WGCNA_PARAMS,       # Co-expression analysis settings
    DECONV_PARAMS,      # Immune deconvolution cell types
    DRUG_PARAMS,        # Drug repurposing filters
    SENSITIVITY_PARAMS, # Robustness analysis ranges
    API_ENDPOINTS,      # External API URLs
)
```

To modify thresholds or parameters, edit `config.py` before running the pipeline. Changes propagate to all downstream scripts automatically.

---

## Troubleshooting

### Common Issues

**GEO download failures**
- GEO servers can be intermittent. Re-run phase 2 scripts — they will retry failed downloads.
- If persistent, check https://www.ncbi.nlm.nih.gov/geo/ for maintenance notices.

**API rate limiting**
- The pipeline includes automatic delays, but high-traffic periods may cause 429 errors.
- Set `NCBI_API_KEY` environment variable for higher rate limits.
- For STRING/g:Profiler timeouts, increase `API_SETTINGS["timeout"]` in `config.py`.

**Memory issues with large datasets**
- GSE59071 (IBD) and GSE221786 (RNA-seq) are the largest datasets.
- Ensure at least 8 GB RAM available. Close other applications if needed.
- If issues persist, process datasets sequentially rather than in parallel.

**Missing probe annotations**
- Phase 3b requires platform annotation files. These are downloaded automatically.
- If annotation download fails, manually download `.annot` files from GEO and place in `data/`.

---

## Project Structure

```
SpA_CrossDisease_MetaAnalysis/
├── README.md                 # This file
├── LICENSE                   # MIT License
├── requirements.txt          # Python dependencies
├── config.py                 # Centralized configuration and parameters
├── run_pipeline.py           # Master pipeline runner (argparse CLI)
├── .gitignore                # Git ignore rules
├── scripts/                  # Analysis scripts (20 total)
│   ├── phase1_dataset_discovery.py
│   ├── phase1b_curate_datasets.py
│   ├── phase2_download_preprocess.py
│   ├── phase2b_fix_groups_download_rnaseq.py
│   ├── phase2c_process_rnaseq.py
│   ├── phase3_differential_expression.py
│   ├── phase3b_probe_mapping.py
│   ├── phase3b_quick_mapping.py
│   ├── phase4_meta_analysis.py
│   ├── phase5_enrichment_and_figures.py
│   ├── phase6_ppi_network.py
│   ├── phase7_wgcna.py
│   ├── phase8_immune_deconvolution.py
│   ├── phase8b_fix_psa_deconvolution.py
│   ├── phase8c_update_crossdisease.py
│   ├── sensitivity_analysis.py
│   ├── graphical_abstract.py
│   ├── generate_supplementary_tables.py
│   ├── inject_figures.py
│   └── inject_figures_v2.py
├── data/
│   └── metadata/
│       └── .gitkeep
├── results/
│   └── .gitkeep
└── figures/
    └── .gitkeep
```

---

## About

This research was conducted by the **Spondyloarthritis AI Computational Biology Institute** (Lisbon, Portugal). Computational analyses were performed by Perplexity Computer under the direction of **Pawel Kuklik**.

The study represents the first cross-disease transcriptomic meta-analysis spanning the full spondyloarthritis spectrum. By integrating data across AS, PsA, IBD, and jSpA, it reveals shared and disease-specific molecular programs in peripheral blood, identifies RNA processing dysregulation as a central feature of AS, and proposes actionable drug repurposing candidates.

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

Copyright (c) 2026 Pawel Kuklik, Spondyloarthritis AI Computational Biology Institute
