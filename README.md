# CCIcompass

> **CCIcompass** is an R package for prioritising cell-cell interactions (CCIs)
> from single-cell RNA-seq data using the **Composite Interaction Score (CIS)**
> — an ensemble metric that aggregates six CCI inference methods into one
> robust, ranked output.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## How it works

CCIcompass wraps the [LIANA](https://github.com/saezlab/liana) framework to run
six CCI methods (connectome, logfc, natmi, sca, cellphonedb, cellchat), filters
interactions by donor, intersects across methods, and ranks everything using
**Rank-Biased Precision (RBP)**:

$$\text{CIS} = \sum_{m=1}^{6} p^{r_m - 1}$$

where $r_m$ is the rank within method $m$ and $p = 0.9$ (default).

The whole pipeline runs in **one function call** from a Seurat object.

---

## Installation

### Requirements
- R ≥ 4.2.0
- macOS (Apple Silicon): run `xcode-select --install` in Terminal first

```r
# Step 1 — CRAN packages
install.packages(c("devtools", "BiocManager", "dplyr", "ggplot2",
                   "readr", "tibble"))

# Step 2 — Bioconductor dependencies
BiocManager::install(c("BiocGenerics", "S4Vectors",
                       "SingleCellExperiment", "SummarizedExperiment"))

# Step 3 — Seurat
install.packages("Seurat")

# Step 4 — LIANA
devtools::install_github("saezlab/liana")

# Step 5 — CCIcompass
devtools::install_github("olhakholod/CCIcompass")
```

---

## Quick start

```r
library(CCIcompass)

obj <- readRDS("my_seurat_object.rds")

# Option A: donor ID already in cell-type label
# e.g. meta.data$celltype.combo == "Epi-Intestine_donor1"
results <- run_pipeline(
  seurat_obj = obj,
  idents_col = "celltype.combo",
  min_donors = 2
)

# Option B: donor ID in a separate column
# e.g. meta.data$celltype == "Epi-Intestine"
#      meta.data$patient_id == "donor1"
results <- run_pipeline(
  seurat_obj = obj,
  idents_col = "celltype",
  donor_col  = "patient_id",
  min_donors = 2
)

# View top results
head(results[, c("source", "target", "ligand", "receptor", "total_rbp_score")])

# Plot
plot_cis_dotplot(results, top_n = 20)
```

### Expected output

| source | target | ligand | receptor | total_rbp_score |
|--------|--------|--------|----------|-----------------|
| Epi-Intestine | DC | MIF | CD74 | 4.82 |
| Stromal | Mac | APP | CD74 | 4.61 |
| ... | ... | ... | ... | ... |

---

## Try it without your own data

```r
library(CCIcompass)

# Load the built-in example dataset (pre-computed LIANA outputs)
data("example_consensus")

# Compute CIS directly
cis <- compute_cis(example_consensus, p = 0.8)
head(cis[, c("ligand", "receptor", "source", "target", "total_rbp_score")])

# Plot
plot_cis_dotplot(cis, top_n = 15)
```

---

## Repository structure

```
CCIcompass/
├── R/
│   ├── run_pipeline.R      # Main function: Seurat → ranked CIS table
│   ├── compute_cis.R       # CIS computation (RBP aggregation)
│   ├── plot_dotplot.R      # Dot plot visualisation
│   └── data.R              # Dataset documentation
├── inst/extdata/           # Example CSV files (used to build .rda datasets)
├── vignettes/
│   └── CIS_tutorial.Rmd   # Full walkthrough with example data
├── scripts/figures/        # Figure reproduction scripts for the paper
│   ├── Figure_3_pipeline.R
│   ├── Figure_4_pipeline.R
│   └── Figure_5_pipeline.R
├── data/                   # Built .rda example datasets
├── DESCRIPTION
├── NAMESPACE
└── session_info.txt
```

---

## Reproducing paper figures

Scripts for all main figures are in `scripts/figures/`. Each script loads
processed data and reproduces the figure panels.

| Figure | Script | Description |
|--------|--------|-------------|
| Figure 3 | `scripts/figures/Figure_3_pipeline.R` | Conserved CCIs — Venn, dot plots, PPI networks |
| Figure 4 | `scripts/figures/Figure_4_pipeline.R` | MIF & CD74 expression validation |
| Figure 5 | `scripts/figures/Figure_5_pipeline.R` | Tissue-specific CCIs |

---

## Citation

```bibtex
@article{kholod2025ccicompass,
  title   = {A Composite Interaction Score: prioritizing cell-cell interactions
             from single-cell RNAseq with application to pre-menopausal
             epithelial barriers},
  author  = {Kholod, Olha and others},
  journal = {Journal of Advanced Research},
  year    = {2026},
  doi     = {10.1016/j.jare.2026.03.046}
}
```

---

## Full tutorial

See [`vignettes/CIS_tutorial.Rmd`](vignettes/CIS_tutorial.Rmd) for a
step-by-step walkthrough with visualised outputs.
