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
library(CCIcompass)

# Import precomputed CCI results from LIANA package
cci_results <- import_results("path/to/cci_output.csv")

# Rank interactions using composite scoring
ranked <- rank_ccis(cci_results)

# Select top interactions for plotting
top_ranked <- ranked %>%
  dplyr::arrange(desc(score)) %>%
  dplyr::slice_head(n = 10)

# Create dot plot
ggplot(top_ranked, aes(x = source_celltype, y = target_celltype, size = score)) +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Top ranked CCIs by CIS metric",
    x = "source-target",
    y = "ligand-receptor",
    size = "CIS score",
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
```

## Citation
#### ...
