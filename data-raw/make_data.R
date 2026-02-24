# data-raw/make_data.R
# Run this ONCE on your local machine after uploading the CSV files to
# inst/extdata/. It converts them into .rda datasets in data/.
#
# Usage (from within the CCIcompass project in RStudio):
#   source("data-raw/make_data.R")
#
# Or from the terminal:
#   Rscript data-raw/make_data.R

library(readr)

# Path to CSV files inside the package
ext <- system.file("extdata", package = "CCIcompass")
if (ext == "") ext <- "inst/extdata"   # fallback when running from source

# ── example_consensus ────────────────────────────────────────────────────────
example_consensus <- read_csv(
  file.path(ext, "example_consensus.csv"),
  show_col_types = FALSE
)

# ── example_cci_list ─────────────────────────────────────────────────────────
methods <- c("cellchat", "cellphonedb", "connectome", "logfc", "natmi", "sca")
example_cci_list <- lapply(methods, function(m) {
  read_csv(file.path(ext, paste0("example_cci_", m, ".csv")),
           show_col_types = FALSE)
})
names(example_cci_list) <- methods

# ── HPA bar chart data ───────────────────────────────────────────────────────
example_hpa_mif <- read_csv(
  file.path(ext, "HPA_MIF_positive_cells.csv"),
  show_col_types = FALSE
)
example_hpa_cd74 <- read_csv(
  file.path(ext, "HPA_CD74_positive_cells.csv"),
  show_col_types = FALSE
)

# ── Save as .rda files ───────────────────────────────────────────────────────
dir.create("data", showWarnings = FALSE)
save(example_consensus, file = "data/example_consensus.rda", compress = "xz")
save(example_cci_list,  file = "data/example_cci_list.rda",  compress = "xz")
save(example_hpa_mif,   file = "data/example_hpa_mif.rda",   compress = "xz")
save(example_hpa_cd74,  file = "data/example_hpa_cd74.rda",  compress = "xz")

message("Done. Four datasets saved to data/")
message("Now commit and push the data/ folder to GitHub.")
