# CCIcompass 

<!-- badges -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **CCIcompass** is an R package for prioritizing and ranking cell–cell interactions (CCIs) from single-cell RNA-seq data using the **Composite Interaction Score (CIS)** — an ensemble metric based on Rank-Biased Precision (RBP).

---

## Table of Contents

1. [Background](#background)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Full Tutorial](#full-tutorial)
5. [Repository Structure](#repository-structure)
6. [Reproducing Paper Figures](#reproducing-paper-figures)
7. [Session Info](#session-info)
8. [Citation](#citation)

---

## Background

Single-cell RNA-seq enables the inference of cell–cell communication, but different tools (CellPhoneDB, CellChat, NATMI, etc.) often disagree on which interactions are most important. **CCIcompass** resolves this by computing a **Composite Interaction Score (CIS)** that:

- Runs six CCI inference methods through the [LIANA](https://github.com/saezlab/liana) framework
- Ranks each interaction per method using **Rank-Biased Precision (RBP)** — a weighted rank aggregation that emphasises top-ranked interactions
- Sums RBP scores across all methods to produce a single, robust CIS per interaction
- Enables cross-tissue and cross-compartment comparison of conserved and tissue-specific CCIs

This package accompanies the paper:

> *"A Composite Interaction Score: prioritizing cell-cell interactions from single-cell RNAseq with application to pre-menopausal epithelial barriers"*
> Kholod O. et al., under review (2026)

---

## Installation

### 1. Install dependencies

```r
# CRAN
install.packages(c("devtools", "dplyr", "ggplot2", "readr", "tibble",
                   "patchwork", "ggrepel", "igraph"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

# LIANA (required for CCI inference)
devtools::install_github("saezlab/liana")
```

### 2. Install CCIcompass

```r
devtools::install_github("olhakholod/CCIcompass")
```

---

## Quick Start

```r
library(CCIcompass)

# Load the built-in example dataset (pre-computed LIANA outputs)
data("example_cci_list")

# Compute CIS across all methods (patience parameter p = 0.8)
cis_results <- compute_cis(example_cci_list, p = 0.8)

# View top 10 interactions
head(cis_results, 10)

# Visualise as a dot plot
plot_cis_dotplot(cis_results, top_n = 20)
```

**Expected output:**

| ligand | receptor | source | target | CIS |
|--------|----------|--------|--------|-----|
| MIF    | CD74     | Epi-Intestine | DC | 4.82 |
| APP    | CD74     | Stromal | Mac | 4.61 |
| ...    | ...      | ...    | ...    | ... |

---

## Full Tutorial

A step-by-step tutorial covering the complete workflow — from raw Seurat object to ranked, cross-tissue CCI comparison — is available in:

📄 **[`vignettes/CIS_tutorial.Rmd`](vignettes/CIS_tutorial.Rmd)**

The tutorial covers:

1. Preparing a Seurat object for LIANA
2. Running all six CCI inference methods
3. Filtering CCIs by donor (same-patient source–target pairs)
4. Retaining CCIs common across donors
5. Computing CIS with `compute_cis()`
6. Filtering to interactions present in ≥2 donors
7. Visualising results: dot plots, Venn diagrams, PPI hub networks
8. Cross-tissue comparison of conserved and tissue-specific CCIs

---

## Repository Structure

```
CCIcompass/
├── R/                        # Package source functions
│   ├── compute_cis.R         #   compute_cis()   – core RBP scoring
│   ├── filter_ccis.R         #   filter_ccis()   – same-patient & donor filters
│   ├── plot_dotplot.R        #   plot_cis_dotplot()
│   ├── plot_network.R        #   plot_ppi_network()
│   └── utils.R               #   internal helpers
│
├── vignettes/
│   └── CIS_tutorial.Rmd      # Full worked tutorial with example data
│
├── data/
│   └── example_cci_list.rda  # Small example dataset (LIANA outputs)
│
├── scripts/
│   └── figures/              # One script per main figure
│       ├── Figure_3_pipeline.R
│       ├── Figure_4_pipeline.R
│       └── Figure_5_pipeline.R
│
├── inst/extdata/             # Raw example input files
├── man/                      # Auto-generated documentation
├── tests/                    # Unit tests (testthat)
├── DESCRIPTION
├── NAMESPACE
├── session_info.txt          # Exact R + package versions used in paper
└── README.md
```

---

## Reproducing Paper Figures

All main figures can be reproduced using the scripts in `scripts/figures/`.
Each script is self-contained: it loads processed data, runs the relevant analysis, and outputs the figure as a PDF.

| Figure | Script | Description |
|--------|--------|-------------|
| Figure 3 | [`scripts/figures/Figure_3_pipeline.R`](scripts/figures/Figure_3_pipeline.R) | Conserved CCIs across epithelial barriers (Venn + dot plots + PPI networks) |
| Figure 4 | [`scripts/figures/Figure_4_pipeline.R`](scripts/figures/Figure_4_pipeline.R) | MIF & CD74 expression validation (GTEx violin plots + HPA bar charts) |
| Figure 5 | [`scripts/figures/Figure_5_pipeline.R`](scripts/figures/Figure_5_pipeline.R) | Tissue-specific CCIs (highlighted Venn + dot plots + PPI networks) |

> **Data availability:** Processed Seurat objects and LIANA outputs are deposited under data. The interactive Shiny app is available at https://cellinteractionsdb.shinyapps.io/shinyapp_ro3_scrna-seq/<img width="468" height="14" alt="image" src="https://github.com/user-attachments/assets/9ddd6433-fdc2-4ad1-9eee-4accd8261e40" />


---

## Session Info

The exact R and package versions used to generate all paper results are recorded in [`session_info.txt`](session_info.txt).

Key versions:
- R 4.3.1
- liana 1.0.0
- Seurat 5.0.0
- ggplot2 3.4.4

---

## Citation

If you use CCIcompass, please cite:

```bibtex
@article{kholod2025ccicompass,
  title   = {A Composite Interaction Score: prioritizing cell-cell interactions
             from single-cell RNAseq with application to pre-menopausal epithelial barriers},
  author  = {Kholod, Olha and others},
  journal = {under review},
  year    = {2026},
  doi     = {10.XXXX/XXXXXXX}
}
```
