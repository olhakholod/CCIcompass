# ============================================================
# CIS Sensitivity Analysis: RBP parameter p
# Generates modified Figure 1 (skin) and Supplementary Figure 1 (PBMC)
# Each shows CIS curves for p ∈ {0.5, 0.6, 0.7, 0.8, 0.9} + average rank
#
# Input CSVs must contain per-method rank columns:
#   cellchat_rank, cellphonedb_rank, connectome_rank,
#   logfc_rank, natmi_rank, sca_rank
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# ─── USER SETTINGS ────────────────────────────────────────────────────────────

SKIN_RANKED_PATH <- "humanSkin_ranked.csv"   # <-- set your path
PBMC_RANKED_PATH <- "pbmc3k_ranked.csv"      # <-- set your path

# Path to your positives (overexpressed interactions ground-truth) files
# These must have columns: source, target, ligand, receptor
SKIN_POS_PATH <- "humanSkin_overexpressed_LR_pairs.csv"   # <-- set your path
PBMC_POS_PATH <- "pbmc3k_overexpressed_LR_pairs_with_bins.csv"  # <-- set your path

# p values to evaluate for sensitivity analysis
P_VALUES <- c(0.5, 0.6, 0.7, 0.8, 0.9)

# k values for precision/recall/fisher/nDCG
K_VALUES <- c(10, 20, 50, 100, 200)

# Output file names
OUT_FIGURE1    <- "Figure_1_sensitivity.pdf"
OUT_SUPP_FIG1  <- "Supplementary_Figure_1_sensitivity.pdf"

# ─── COLOUR / LINE TYPE PALETTE ───────────────────────────────────────────────
# Gold for average rank, shades of purple for each p value
# Solid line = p=0.9 (the default); others use distinct line types

PURPLE_SHADES <- c(
  "0.5" = "#d8b4fe",   # lightest purple
  "0.6" = "#a855f7",
  "0.7" = "#7c3aed",
  "0.8" = "#5b21b6",
  "0.9" = "#3b0764"    # darkest purple (default)
)

LINE_TYPES <- c(
  "0.5" = "dotted",
  "0.6" = "dotdash",
  "0.7" = "dashed",
  "0.8" = "longdash",
  "0.9" = "solid"
)

AVG_RANK_COLOR <- "#D4A017"   # gold

# ─── HELPER FUNCTIONS ─────────────────────────────────────────────────────────

rank_cols <- c("cellchat_rank", "cellphonedb_rank", "connectome_rank",
               "logfc_rank", "natmi_rank", "sca_rank")

#' Compute CIS for a given p value from individual method rank columns
#' Formula used: sum of p^(r-1) across all 6 methods (no (1-p) normalisation factor,
#' matching the implementation in the stored CIS column of the input CSVs)
compute_cis <- function(df, p) {
  rbp_scores <- map_dfc(rank_cols, function(col) {
    r <- df[[col]]
    p^(r - 1)
  })
  rowSums(rbp_scores)
}

#' Build a ranked evaluation table for a single method
build_ranked_table <- function(df, method_label, score_vec, positives_key, higher_better = TRUE) {
  df <- df %>%
    mutate(.score = score_vec,
           .key   = paste(source, target, ligand, receptor, sep = "||"))
  if (higher_better) {
    df <- df %>% arrange(desc(.score))
  } else {
    df <- df %>% arrange(.score)
  }
  df %>%
    mutate(
      method     = method_label,
      position   = row_number(),
      is_positive = as.integer(.key %in% positives_key),
      cum_tp     = cumsum(is_positive),
      precision  = cum_tp / position,
      recall     = cum_tp / length(positives_key)
    )
}

dcg_at_k <- function(rel, k) {
  idx <- seq_len(min(k, length(rel)))
  sum(rel[idx] / log2(idx + 1))
}
idcg_at_k <- function(n_pos, k) {
  rel <- rep(1, min(n_pos, k))
  sum(rel / log2(seq_along(rel) + 1))
}

fisher_neglog10p_at_k <- function(rel, k, N, P) {
  tp <- sum(rel[seq_len(k)])
  fp <- k - tp
  fn <- P - tp
  tn <- N - k - fn
  mat <- matrix(c(tp, fp, fn, tn), nrow = 2)
  ft  <- suppressWarnings(fisher.test(mat))
  -log10(ft$p.value)
}

#' Compute all metrics at each k for one ranked table
compute_metrics_at_k <- function(ranked_df, ks, P, N) {
  rel <- ranked_df$is_positive
  map_dfr(ks[ks <= N], function(k) {
    ndcg <- if (idcg_at_k(P, k) > 0) dcg_at_k(rel, k) / idcg_at_k(P, k) else NA_real_
    tibble(
      k                = k,
      precision_at_k   = sum(rel[seq_len(k)]) / k,
      recall_at_k      = sum(rel[seq_len(k)]) / P,
      fisher_neglog10p = fisher_neglog10p_at_k(rel, k, N, P),
      nDCG_at_k        = ndcg
    )
  })
}

# ─── MAIN PIPELINE ────────────────────────────────────────────────────────────

run_sensitivity <- function(ranked_path, pos_path, p_values, k_values,
                             dataset_label) {

  ranked <- read_csv(ranked_path, show_col_types = FALSE)
  pos_raw <- read_csv(pos_path, show_col_types = FALSE)

  key_cols <- c("source", "target", "ligand", "receptor")
  positives_key <- pos_raw %>%
    distinct(across(all_of(key_cols))) %>%
    mutate(.key = paste(source, target, ligand, receptor, sep = "||")) %>%
    pull(.key)

  P <- length(positives_key)
  N <- nrow(ranked)
  ks <- k_values[k_values <= N]

  # Build all ranked tables: one per p value + average rank
  all_ranked <- list()

  # CIS for each p
  for (p in p_values) {
    label <- paste0("CIS (p=", p, ")")
    cis_scores <- compute_cis(ranked, p)
    all_ranked[[label]] <- build_ranked_table(ranked, label, cis_scores,
                                               positives_key, higher_better = TRUE)
  }

  # Average rank baseline (lower average rank = better, so invert score)
  avg_rank_score <- -rowMeans(ranked[, rank_cols], na.rm = TRUE)
  all_ranked[["average rank"]] <- build_ranked_table(ranked, "average rank",
                                                       avg_rank_score, positives_key,
                                                       higher_better = TRUE)

  # Compute metrics for all methods
  metrics <- imap_dfr(all_ranked, function(df, method) {
    compute_metrics_at_k(df, ks, P, N) %>%
      mutate(method = method)
  })

  list(metrics = metrics, ks = ks, dataset = dataset_label)
}

# ─── PLOTTING FUNCTION ────────────────────────────────────────────────────────

make_sensitivity_plots <- function(result, p_values, panel_labels = c("B","C","D","E"),
                                    title_suffix = "") {

  metrics <- result$metrics
  ks      <- result$ks

  # Define aesthetics for each method
  # average rank = gold solid; CIS lines = purple shades + line types
  p_labels  <- paste0("CIS (p=", p_values, ")")
  all_labels <- c(p_labels, "average rank")

  color_map <- c(
    setNames(PURPLE_SHADES[as.character(p_values)], p_labels),
    "average rank" = AVG_RANK_COLOR
  )
  ltype_map <- c(
    setNames(LINE_TYPES[as.character(p_values)], p_labels),
    "average rank" = "solid"
  )
  size_map <- c(
    setNames(rep(0.8, length(p_labels)), p_labels),
    "average rank" = 1.2
  )

  # Order factor for legend: average rank first, then p from low to high
  factor_levels <- c("average rank", rev(p_labels))
  metrics <- metrics %>%
    mutate(method = factor(method, levels = factor_levels))

  base_theme <- theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(colour = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position  = "right",
      legend.key.width = unit(1.5, "cm"),
      plot.title       = element_text(size = 11, hjust = 0.5)
    )

  # Precision at k
  pB <- metrics %>%
    ggplot(aes(x = k, y = precision_at_k,
               color = method, linetype = method, linewidth = method)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = ks) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = ltype_map) +
    scale_linewidth_manual(values = size_map) +
    labs(x = "k", y = "Precision at k",
         title = paste0(panel_labels[1], ". Precision at k"),
         color = "Method", linetype = "Method", linewidth = "Method") +
    base_theme

  # Recall at k
  pC <- metrics %>%
    ggplot(aes(x = k, y = recall_at_k,
               color = method, linetype = method, linewidth = method)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = ks) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = ltype_map) +
    scale_linewidth_manual(values = size_map) +
    labs(x = "k", y = "Recall at k",
         title = paste0(panel_labels[2], ". Recall at k"),
         color = "Method", linetype = "Method", linewidth = "Method") +
    base_theme

  # Fisher's exact test enrichment
  pD <- metrics %>%
    ggplot(aes(x = k, y = fisher_neglog10p,
               color = method, linetype = method, linewidth = method)) +
    geom_hline(yintercept = -log10(0.05),  linetype = "dotted", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.001), linetype = "dashed",  color = "black", linewidth = 0.4) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = ks) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = ltype_map) +
    scale_linewidth_manual(values = size_map) +
    labs(x = "k", y = expression(-log[10](p-value)),
         title = paste0(panel_labels[3], ". Fisher's exact test enrichment at top k"),
         color = "Method", linetype = "Method", linewidth = "Method") +
    base_theme

  # nDCG at k
  pE <- metrics %>%
    ggplot(aes(x = k, y = nDCG_at_k,
               color = method, linetype = method, linewidth = method)) +
    geom_line() +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = ks) +
    scale_color_manual(values = color_map) +
    scale_linetype_manual(values = ltype_map) +
    scale_linewidth_manual(values = size_map) +
    labs(x = "k", y = "nDCG at k",
         title = paste0(panel_labels[4], ". nDCG at k"),
         color = "Method", linetype = "Method", linewidth = "Method") +
    base_theme

  list(precision = pB, recall = pC, fisher = pD, ndcg = pE)
}

# ─── RUN ──────────────────────────────────────────────────────────────────────

message("Running sensitivity analysis on skin dataset...")
skin_result <- run_sensitivity(
  ranked_path  = SKIN_RANKED_PATH,
  pos_path     = SKIN_POS_PATH,
  p_values     = P_VALUES,
  k_values     = K_VALUES,
  dataset_label = "Human Skin"
)

message("Running sensitivity analysis on PBMC dataset...")
pbmc_result <- run_sensitivity(
  ranked_path  = PBMC_RANKED_PATH,
  pos_path     = PBMC_POS_PATH,
  p_values     = P_VALUES,
  k_values     = K_VALUES,
  dataset_label = "PBMC 3K"
)

# ─── FIGURE 1 (Skin — panels B–E, replaces current panels) ───────────────────
skin_plots <- make_sensitivity_plots(
  skin_result,
  p_values      = P_VALUES,
  panel_labels  = c("B", "C", "D", "E")
)

fig1_bottom <- (skin_plots$precision | skin_plots$recall) /
               (skin_plots$fisher    | skin_plots$ndcg) +
  plot_annotation(
    title    = "Figure 1 (panels B–E): CIS sensitivity to RBP parameter p — Human Skin dataset",
    subtitle = "Gold line = average rank baseline. Purple lines = CIS at p ∈ {0.5, 0.6, 0.7, 0.8, 0.9}.",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

ggsave(OUT_FIGURE1, fig1_bottom, width = 12, height = 8, device = "pdf")
message("Saved: ", OUT_FIGURE1)

# ─── SUPPLEMENTARY FIGURE 1 (PBMC — panels A–D) ──────────────────────────────
pbmc_plots <- make_sensitivity_plots(
  pbmc_result,
  p_values      = P_VALUES,
  panel_labels  = c("A", "B", "C", "D")
)

supp1 <- (pbmc_plots$precision | pbmc_plots$recall) /
          (pbmc_plots$fisher    | pbmc_plots$ndcg) +
  plot_annotation(
    title    = "Supplementary Figure 1: CIS sensitivity to RBP parameter p — PBMC 3K dataset",
    subtitle = "Gold line = average rank baseline. Purple lines = CIS at p ∈ {0.5, 0.6, 0.7, 0.8, 0.9}.",
    theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                     plot.subtitle = element_text(size = 10))
  )

ggsave(OUT_SUPP_FIG1, supp1, width = 10, height = 8, device = "pdf")
message("Saved: ", OUT_SUPP_FIG1)

message("Done. Note: set SKIN_POS_PATH and PBMC_POS_PATH to your ground-truth positives files.")
