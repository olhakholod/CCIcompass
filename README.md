# CCIcompass

## Description
**CCIcompass** is an R package designed to **prioritize and rank cell-cell interactions** based on composite interaction score (CIS). CIS based on ranked-biased precision (RBP), which prioritizes top-ranked interactions from methods implemented in LIANA package and calculates a weighted composite score that reflects agreement across tools, while accounting for rank position. It streamlines analysis and comparison across approaches, helping researchers identify biologically meaningful interactions in single-cell datasets.

## Features
- Harmonizes outputs from different CCI inference tools (e.g., LIANA, CellPhoneDB, CellChat)
- Computes composite scores (CIS) to rank ligand-receptor interactions
- Compatible with Seurat and SingleCellExperiment objects

## Dependencies
- dplyr
- ggplot2
- Seurat
- readr
- tibble

## Installation
install.packages("devtools")
devtools::install_github("olhakholod/CCIcompass")

## Usage
```r
library(ccicompass)
library(ggplot2)

# Load your preprocessed single-cell dataset (e.g., Seurat object)
seu <- readRDS("your_seurat_object.rds")

# Import precomputed CCI results from another tool (e.g., LIANA, CellPhoneDB)
cci_results <- import_results("path/to/your_liana_output.csv")

# Rank interactions using composite scoring
ranked <- rank_ccis(cci_results)

# Select top interactions for plotting
top_ranked <- ranked %>%
  dplyr::arrange(desc(score)) %>%
  dplyr::slice_head(n = 30)

# Create dot plot
ggplot(top_ranked, aes(x = source_celltype, y = target_celltype, size = score, color = score)) +
  geom_point() +
  scale_color_viridis_c() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top Ranked Cell-Cell Interactions",
    x = "Ligand Source",
    y = "Receptor Target",
    size = "CCI Score",
    color = "CCI Score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```

## Citation
#### ...
