
# =============================================================================
# make_data.R  –  Run this ONCE on your local machine to create the
#                  package example datasets (.rda files in data/)
#
# Usage:
#   Rscript data-raw/make_data.R
#   (or source("data-raw/make_data.R") from within the package project)
# =============================================================================

library(readr)
library(usethis)

DATA_RAW <- "inst/extdata"   # where the CSVs live in the package

# ── 1. example_cci_list  (named list, one df per method) ─────────────────────
methods <- c("cellchat","cellphonedb","connectome","logfc","natmi","sca")

example_cci_list <- lapply(methods, function(m) {
  read_csv(file.path(DATA_RAW, paste0("example_cci_", m, ".csv")),
           show_col_types = FALSE)
})
names(example_cci_list) <- methods

# ── 2. example_consensus  (joined table, input to compute_cis()) ─────────────
example_consensus <- read_csv(
  file.path(DATA_RAW, "example_consensus.csv"),
  show_col_types = FALSE
)

# ── 3. HPA bar chart data ─────────────────────────────────────────────────────
example_hpa_mif  <- read_csv(file.path(DATA_RAW, "HPA_MIF_positive_cells.csv"),
                              show_col_types = FALSE)
example_hpa_cd74 <- read_csv(file.path(DATA_RAW, "HPA_CD74_positive_cells.csv"),
                              show_col_types = FALSE)

# ── Save as .rda in data/ ────────────────────────────────────────────────────
usethis::use_data(example_cci_list,  overwrite = TRUE)
usethis::use_data(example_consensus, overwrite = TRUE)
usethis::use_data(example_hpa_mif,   overwrite = TRUE)
usethis::use_data(example_hpa_cd74,  overwrite = TRUE)

message("Done. Four datasets saved to data/")
