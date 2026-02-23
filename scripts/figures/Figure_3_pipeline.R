# =============================================================================
# Figure 3 – Complete Pipeline
# Cell-Cell Interaction (CCI) analysis across epithelial, immune, and stromal
# compartments in intestine, skin, and uterus tissues.
#
# Pipeline overview:
#   Step 1:  Run LIANA CCI inference on Seurat objects (per tissue)
#   Step 2:  Filter CCIs to same-patient source–target pairs
#   Step 3:  Retain CCIs common across all donors (≥3 or ≥2, see steps 3/6)
#   Step 4:  Intersect CCIs common across all 6 methods
#   Step 5:  Rank CCIs by Rank-Biased Precision (RBP) score
#   Step 6:  Filter to interactions present in ≥2 donors (final filter)
#   Step 7:  Venn diagrams – shared CCIs across tissues (panels A–C)
#   Step 8:  Dot plots – top CCIs per compartment (panels D–F)
#   Step 9:  PPI hub-gene network plots (panels G–I)
#   Step 10: Assemble and export Figure 3 PDF
#
# Input files required (update paths in the CONFIG section below):
#   - Seurat RDS objects for intestine, skin, and uterus
#   - Pre-ranked CCI CSVs per compartment (epi/immune/stromal × 3 tissues)
#   - STRING interaction TSV files for hub-gene networks (3 compartments)
#
# Output:
#   - Figure_3.pdf  (3 × 3 panel layout)
# =============================================================================


# =============================================================================
# 0. LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
  library(liana)
  library(dplyr)
  library(ggplot2)
  library(ggvenn)       # or VennDiagram; used for panels A-C
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(readr)
  library(scales)
  library(ggrepel)
  library(packcircles)
  library(patchwork)    # final figure assembly
})


# =============================================================================
# 1. USER CONFIG  ← update all paths here
# =============================================================================

# --- Seurat objects (Step 1) -------------------------------------------------
SEURAT_INTESTINE <- "path/to/female_intestine_final.rds"
SEURAT_SKIN      <- "path/to/female_skin_final.rds"
SEURAT_UTERUS    <- "path/to/female_uterus_final.rds"

# --- Donor IDs present in each dataset (Step 3) -----------------------------
INTESTINE_DONORS <- c("A26 (386C)", "A38 (432C)", "A30 (398B)")
SKIN_DONORS      <- c("donor1", "donor2", "donor3")   # update as needed
UTERUS_DONORS    <- c("donor1", "donor2", "donor3")   # update as needed

# --- Pre-ranked CCI CSVs for dot plots (Steps 5/6 outputs, Step 8 input) ----
# Each CSV must have columns: source, target, ligand, receptor,
#                              total_rbp_score, patient_id
DOTPLOT_FILES <- list(
  epithelial = list(
    intestine = "path/to/01_CIS_ranked_100_CCIs/01_epi_intestine.csv",
    skin      = "path/to/01_CIS_ranked_100_CCIs/02_keratinocytes_skin.csv",
    uterus    = "path/to/01_CIS_ranked_100_CCIs/03_epi_ciliated_uterus.csv"
  ),
  immune = list(
    intestine = "path/to/01_CIS_ranked_100_CCIs/01_mac_intestine.csv",
    skin      = "path/to/01_CIS_ranked_100_CCIs/02_macro_skin.csv",
    uterus    = "path/to/01_CIS_ranked_100_CCIs/03_myeloid_uterus.csv"
  ),
  stromal = list(
    intestine = "path/to/01_CIS_ranked_100_CCIs/01_stromal_intestine.csv",
    skin      = "path/to/01_CIS_ranked_100_CCIs/02_fibroblasts_skin.csv",
    uterus    = "path/to/01_CIS_ranked_100_CCIs/03_fibroblasts_uterus.csv"
  )
)

# --- STRING TSV files for hub-gene networks (Step 9) ------------------------
STRING_FILES <- list(
  epithelial = "path/to/hub_genes/01_epi_common_string_interactions.tsv",
  immune     = "path/to/hub_genes/02_immune_common_string_interactions.tsv",
  stromal    = "path/to/hub_genes/03_stromal_common_string_interactions.tsv"
)

# --- Output ------------------------------------------------------------------
OUT_PDF <- "Figure_3.pdf"

# --- Venn diagram counts (panels A-C) ----------------------------------------
# If you prefer to re-derive these counts programmatically from your data,
# replace the hard-coded values below with the actual set sizes from your
# filtered CCI lists.  The values shown reproduce the published figure.
VENN_EPI    <- list(intestine = 39, skin = 37, uterus = 49,
                    int_skin = 18, int_uter = 16, skin_uter = 12, all = 7)
VENN_IMMUNE <- list(intestine = 21, skin = 27, uterus = 36,
                    int_skin = 16, int_uter = 14, skin_uter = 16, all = 10)
VENN_STROM  <- list(intestine = 27, skin = 38, uterus = 41,
                    int_skin = 19, int_uter = 19, skin_uter = 13, all = 8)

# RBP patience parameter (Step 5)
RBP_P <- 0.8


# =============================================================================
# STEP 1 – Run LIANA on each tissue
# =============================================================================
# Produces per-method CSV files used in subsequent steps.
# Skip / comment out this block if LIANA outputs already exist on disk.

run_liana <- function(rds_path, tissue_label) {
  message("[ Step 1 ] Running LIANA for: ", tissue_label)
  obj <- readRDS(rds_path)

  # Four rank-based methods
  res_multi <- liana_wrap(obj,
                          idents_col = "celltype.combo",
                          method     = c("connectome", "logfc", "natmi", "sca"),
                          resource   = "CellPhoneDB")
  for (m in names(res_multi)) {
    write.csv(res_multi[[m]],
              paste0(tissue_label, "_", m, ".csv"),
              row.names = FALSE)
  }

  # CellPhoneDB permutation-based method
  options(future.globals.maxSize = 10000 * 1024^2)
  res_cpdb <- liana_wrap(obj,
                         idents_col         = "celltype.combo",
                         method             = "cellphonedb",
                         resource           = "CellPhoneDB",
                         permutation.params = list(nperms    = 100,
                                                   parallelize = FALSE,
                                                   workers   = 4),
                         expr_prop          = 0.05)
  write.csv(res_cpdb,
            paste0(tissue_label, "_cellphonedb.csv"),
            row.names = FALSE)

  # CellChat
  res_cc <- liana_wrap(obj,
                       idents_col = "celltype.combo",
                       method     = "call_cellchat",
                       resource   = "CellPhoneDB")
  write.csv(res_cc,
            paste0(tissue_label, "_cellchat.csv"),
            row.names = FALSE)

  invisible(NULL)
}

# Uncomment to run:
# run_liana(SEURAT_INTESTINE, "female_intestine")
# run_liana(SEURAT_SKIN,      "female_skin")
# run_liana(SEURAT_UTERUS,    "female_uterus")


# =============================================================================
# STEP 2 – Filter CCIs to same-patient source–target pairs
# =============================================================================

filter_same_patient <- function(csv_path, out_path) {
  message("[ Step 2 ] Filtering by patient: ", basename(csv_path))
  data <- read.csv(csv_path)

  extract_patient_id <- function(x) sub(".*_(.*?)\\s.*", "\\1", x)

  filtered <- data %>%
    mutate(src_pid = extract_patient_id(source),
           tgt_pid = extract_patient_id(target)) %>%
    filter(src_pid == tgt_pid) %>%
    select(-src_pid, -tgt_pid)

  write.csv(filtered, out_path, row.names = FALSE)
  invisible(filtered)
}

# Example usage (repeat for every tissue × method combination):
# filter_same_patient("female_intestine_sca.csv",
#                     "filtered_female_intestine_sca.csv")


# =============================================================================
# STEP 3 – Retain CCIs common across all donors (present in ≥N donors)
# =============================================================================

filter_common_across_donors <- function(csv_path, donors, out_path,
                                        n_donors = NULL) {
  message("[ Step 3 ] Common across donors: ", basename(csv_path))
  data    <- read.csv(csv_path)
  n_req   <- if (is.null(n_donors)) length(donors) else n_donors

  filtered <- data %>%
    filter(patient_id %in% donors) %>%
    group_by(source, target, ligand, receptor) %>%
    filter(n() == n_req) %>%
    ungroup()

  write.csv(filtered, out_path, row.names = FALSE)
  invisible(filtered)
}

# Example usage:
# filter_common_across_donors("filtered_female_intestine_sca.csv",
#                              INTESTINE_DONORS,
#                              "common_filtered_female_intestine_sca.csv")


# =============================================================================
# STEP 4 – Intersect CCIs across all 6 methods
# =============================================================================

intersect_across_methods <- function(tissue_label, out_path) {
  message("[ Step 4 ] Intersecting methods for: ", tissue_label)
  methods <- c("cellchat", "cellphonedb", "connectome", "logfc", "natmi", "sca")
  dfs     <- lapply(methods, function(m) {
    path <- paste0("common_filtered_", tissue_label, "_", m, ".csv")
    df   <- read.csv(path)
    df   %>% rename_with(~paste0(m, "_", .), -c(patient_id, source, target,
                                                  ligand, receptor))
  })

  common <- dfs[[1]]
  for (i in 2:length(dfs)) {
    common <- inner_join(common, dfs[[i]],
                         by = c("patient_id", "source", "target",
                                "ligand", "receptor"))
  }

  write.csv(common, out_path, row.names = FALSE)
  invisible(common)
}

# Example usage:
# intersect_across_methods("female_intestine",
#                           "intestine_common_interactions_across_methods.csv")


# =============================================================================
# STEP 5 – Rank CCIs by Rank-Biased Precision (RBP) score
# =============================================================================

rank_by_rbp <- function(csv_path, out_path, p = RBP_P) {
  message("[ Step 5 ] RBP ranking: ", basename(csv_path))
  data <- read.csv(csv_path)

  rbp_score <- function(rank_vec, patience) patience^(rank_vec - 1)

  ranked <- data %>%
    mutate(
      cellchat_rank    = rank(-cellchat_prob,           ties.method = "min"),
      cellphonedb_rank = rank(-cellphonedb_lr.mean,     ties.method = "min"),
      connectome_rank  = rank(-connectome_weight_sc,    ties.method = "min"),
      logfc_rank       = rank(-logfc_logfc_comb,        ties.method = "min"),
      natmi_rank       = rank(-natmi_prod_weight,       ties.method = "min"),
      sca_rank         = rank(-sca_LRscore,             ties.method = "min")
    ) %>%
    mutate(
      rbp_cellchat    = rbp_score(cellchat_rank,    p),
      rbp_cellphonedb = rbp_score(cellphonedb_rank, p),
      rbp_connectome  = rbp_score(connectome_rank,  p),
      rbp_logfc       = rbp_score(logfc_rank,        p),
      rbp_natmi       = rbp_score(natmi_rank,        p),
      rbp_sca         = rbp_score(sca_rank,          p)
    ) %>%
    mutate(
      total_rbp_score = rbp_cellchat + rbp_cellphonedb + rbp_connectome +
                        rbp_logfc    + rbp_natmi        + rbp_sca
    ) %>%
    arrange(desc(total_rbp_score))

  write.csv(ranked, out_path, row.names = FALSE)
  invisible(ranked)
}

# Example usage:
# rank_by_rbp("intestine_common_interactions_across_methods.csv",
#             "ranked_interactions_by_rbp_intestine.csv")


# =============================================================================
# STEP 6 – Filter to interactions present in ≥2 donors (final filter)
# =============================================================================

filter_min_donors <- function(csv_path, out_path, min_donors = 2) {
  message("[ Step 6 ] Filter ≥", min_donors, " donors: ", basename(csv_path))
  data <- read.csv(csv_path)

  filtered <- data %>%
    group_by(source, target, ligand, receptor) %>%
    filter(n_distinct(patient_id) >= min_donors) %>%
    slice(1) %>%
    ungroup()

  write.csv(filtered, out_path, row.names = FALSE)
  invisible(filtered)
}

# Example usage:
# filter_min_donors("ranked_interactions_by_rbp_intestine.csv",
#                   "final_filtered_intestine.csv")


# =============================================================================
# STEP 7 – Venn diagrams (panels A, B, C)
# =============================================================================
# Helper to draw a 3-set Euler/Venn diagram from pre-computed counts.
# The counts are defined in the CONFIG section (VENN_EPI, etc.).

draw_venn <- function(counts, labels = c("intestine", "skin", "uterus"),
                      fill_cols = c("#E5C4A1", "#C8DFA4", "#A4C8DF")) {
  # Build a named list suitable for ggvenn
  n      <- counts
  region <- list(
    intestine          = seq_len(n$intestine),
    skin               = seq_len(n$skin)     + n$intestine,
    uterus             = seq_len(n$uterus)   + n$intestine + n$skin
  )
  # ggvenn expects sets as named lists of element IDs
  # We synthesise artificial element IDs that reproduce the overlap counts.
  total  <- n$intestine + n$skin + n$uterus +
            n$int_skin + n$int_uter + n$skin_uter + n$all
  i <- 0L
  mk <- function(k) { r <- seq_len(k) + i; i <<- i + k; r }

  sets <- list(
    intestine = c(mk(n$all), mk(n$int_skin), mk(n$int_uter), mk(n$intestine)),
    skin      = c(sets$intestine[seq_len(n$all + n$int_skin)],
                  mk(n$skin_uter), mk(n$skin)),
    uterus    = c(sets$intestine[seq_len(n$all)],
                  sets$skin[seq_len(n$all + n$int_skin + n$skin_uter)][
                    tail(seq_len(n$all + n$int_skin + n$skin_uter),
                         n$all + n$skin_uter)],
                  mk(n$int_uter), mk(n$uterus))
  )
  # Simpler approach: use ggvenn with synthetic integer vectors
  all_ids    <- mk2 <- local({
    ctr <- 0L
    function(k) { r <- ctr + seq_len(k); ctr <<- ctr + k; r }
  })
  shared_all      <- mk2(n$all)
  shared_is       <- c(shared_all, mk2(n$int_skin  - n$all))   # int ∩ skin (excl. triple)
  shared_iu       <- c(shared_all, mk2(n$int_uter  - n$all))
  shared_su       <- c(shared_all, mk2(n$skin_uter - n$all))
  only_int        <- mk2(n$intestine)
  only_skin       <- mk2(n$skin)
  only_uter       <- mk2(n$uterus)

  sets_final <- list(
    intestine = c(shared_all, shared_is, shared_iu, only_int),
    skin      = c(shared_all, shared_is, shared_su, only_skin),
    uterus    = c(shared_all, shared_iu, shared_su, only_uter)
  )
  names(sets_final) <- labels

  ggvenn(sets_final,
         fill_color    = fill_cols,
         stroke_size   = 0.5,
         set_name_size = 4,
         text_size     = 3.5)
}

panel_A <- draw_venn(VENN_EPI,
                     labels = c("epithelial intestine",
                                "epithelial skin",
                                "epithelial uterus"))
panel_B <- draw_venn(VENN_IMMUNE,
                     labels = c("immune intestine",
                                "immune skin",
                                "immune uterus"),
                     fill_cols = c("#D4A4E0", "#A4D4E0", "#E0D4A4"))
panel_C <- draw_venn(VENN_STROM,
                     labels = c("stromal intestine",
                                "stromal skin",
                                "stromal uterus"),
                     fill_cols = c("#A4E0B4", "#E0A4A4", "#A4A4E0"))


# =============================================================================
# STEP 8 – Dot plots (panels D, E, F)
# =============================================================================

# Ligand-receptor pairs of interest per compartment
LR_EPI <- c("CD99-PILRA", "APP-CD74", "HLA-E-KLRC1",
             "RPS19-C5AR1", "MIF-CD74", "FAM3C-HLA-C", "LTB-LTBR")

LR_IMMUNE <- c("COPA-CD74", "APP-CD74", "CCL3-CCR1", "LILRB4-LAIR1",
                "CCL5-CCR1", "C3-C3AR1", "RPS19-C5AR1", "MIF-CD74",
                "CSF1-CSF1R", "TNFSF13B-HLA-DPB1")

LR_STROM <- c("COPA-CD74", "APP-CD74", "HLA-E-KLRC1", "FAM3C-HLA-C",
               "C3-C3AR1", "RPS19-C5AR1", "MIF-CD74", "CXCL12-CXCR4")

# Tissue colour palette used in legend
TISSUE_COLS <- c(intestine = "#E5835A", skin = "#82B0D9", uterus = "#7DC99E")

build_dotplot <- function(file_list, lr_of_interest, title_str) {
  data <- map_dfr(names(file_list), function(tiss) {
    read_csv(file_list[[tiss]], show_col_types = FALSE) %>%
      mutate(tissue = tiss,
             LR     = paste(ligand, receptor, sep = "-"))
  })

  plot_data <- data %>%
    filter(LR %in% lr_of_interest) %>%
    mutate(source_target = paste(source, target, sep = "-"),
           LR = factor(LR, levels = rev(sort(unique(LR)))),
           tissue = factor(tissue, levels = c("intestine", "skin", "uterus")))

  ggplot(plot_data, aes(x = LR, y = source_target)) +
    geom_point(aes(size = total_rbp_score, colour = tissue), alpha = 0.85) +
    scale_size_continuous(name   = "CIS",
                          range  = c(1, 7),
                          breaks = seq(0.5, 4.5, 0.5)) +
    scale_colour_manual(values = TISSUE_COLS) +
    coord_flip() +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
      axis.text.y  = element_text(size = 7),
      strip.text   = element_text(size = 9, face = "bold"),
      panel.grid   = element_line(colour = "grey92"),
      legend.key.size = unit(0.4, "cm")
    ) +
    labs(x = "CCIs", y = "source–target cells",
         colour = NULL, title = title_str)
}

panel_D <- build_dotplot(DOTPLOT_FILES$epithelial, LR_EPI,    "D")
panel_E <- build_dotplot(DOTPLOT_FILES$immune,     LR_IMMUNE, "E.")
panel_F <- build_dotplot(DOTPLOT_FILES$stromal,    LR_STROM,  "F.")


# =============================================================================
# STEP 9 – PPI hub-gene network plots (panels G, H, I)
# =============================================================================

# ---------- 9a. Helpers (from hub_gene_analysis_drl_V2.R) --------------------

BET_LIMITS <- c(0, 0.5)
BET_BREAKS <- seq(0, 0.5, by = 0.1)

DRL_OPTIONS <- list(
  init.iterations = 300L,
  maxiter         = 2000L,
  cool.factor     = 0.95
)

read_edges_string <- function(path) {
  tb     <- read_tsv(path, show_col_types = FALSE)
  cn     <- names(tb); lc <- tolower(cn); map <- setNames(cn, lc)
  p1     <- map[["protein1"]]; p2 <- map[["protein2"]]
  if (is.na(p1) || is.na(p2)) {
    p1 <- map[["preferredname_a"]]; p2 <- map[["preferredname_b"]]
  }
  scorecol <- map[["combined_score"]]
  if (is.na(scorecol)) scorecol <- map[["combinedscore"]]
  stopifnot(!is.na(p1), !is.na(p2))

  ed <- tb %>%
    transmute(
      protein1      = .data[[p1]],
      protein2      = .data[[p2]],
      combined_score = if (!is.na(scorecol))
                         suppressWarnings(as.numeric(.data[[scorecol]]))
                       else NA_real_
    ) %>%
    filter(!is.na(protein1), !is.na(protein2), protein1 != protein2) %>%
    distinct()

  if (any(!is.na(ed$combined_score))) {
    mx <- suppressWarnings(max(ed$combined_score, na.rm = TRUE))
    if (is.finite(mx) && mx > 10) ed$combined_score <- ed$combined_score / 1000
  }
  ed
}

graph_from_edges <- function(ed) {
  if (!all(is.na(ed$combined_score))) {
    ed <- ed %>% mutate(distance = ifelse(combined_score > 0,
                                          1 / combined_score, Inf))
    graph_from_data_frame(
      ed %>% select(protein1, protein2,
                    weight = combined_score, distance),
      directed = FALSE)
  } else {
    graph_from_data_frame(ed %>% select(protein1, protein2),
                          directed = FALSE)
  }
}

prune_to_lcc <- function(g, min_comp_size = 2) {
  comp   <- components(g)
  keep   <- which(comp$csize >= min_comp_size)
  if (!length(keep)) return(delete_vertices(g, V(g)))
  lcc_id <- keep[which.max(comp$csize[keep])]
  induced_subgraph(g, vids = V(g)[comp$membership == lcc_id])
}

node_metrics <- function(g, dataset_name) {
  use_dist <- "distance" %in% edge_attr_names(g)
  deg      <- degree(g)
  bet      <- betweenness(g, directed = FALSE, normalized = TRUE,
                          weights = if (use_dist) E(g)$distance else NULL)
  eig      <- eigen_centrality(
                g, directed = FALSE,
                weights = if ("weight" %in% edge_attr_names(g))
                            E(g)$weight else NULL)$vector
  tibble(gene = names(deg),
         degree      = as.numeric(deg),
         betweenness = as.numeric(bet),
         eigenvector = as.numeric(eig)) %>%
    mutate(r_deg  = rank(-degree,      ties.method = "min"),
           r_bet  = rank(-betweenness, ties.method = "min"),
           r_eig  = rank(-eigenvector, ties.method = "min"),
           mean_rank     = rowMeans(cbind(r_deg, r_bet, r_eig), na.rm = TRUE),
           consensus_rank = rank(mean_rank, ties.method = "min"),
           dataset = dataset_name) %>%
    arrange(consensus_rank, desc(degree))
}

theme_pub_void <- function(base_size = 11) {
  theme_void(base_size = base_size) +
    theme(plot.title      = element_text(hjust = 0.5, face = "bold"),
          legend.title    = element_text(face = "bold"),
          plot.margin     = margin(6, 6, 6, 6))
}

# ---------- 9b. Build a single PPI panel -------------------------------------

build_ppi_panel <- function(path, dataset_name,
                             top_n_label = 10, seed = 7L) {
  ed <- read_edges_string(path)
  stopifnot(nrow(ed) > 0)

  g  <- graph_from_edges(ed)
  if ("weight" %in% edge_attr_names(g)) {
    w       <- E(g)$weight
    w       <- scales::rescale(sqrt(pmax(w, 1e-6)), to = c(0.15, 0.75))
    E(g)$weight <- w
  }

  g <- prune_to_lcc(g, min_comp_size = 2)
  if (gorder(g) == 0 || gsize(g) == 0)
    stop("Empty graph after pruning for ", dataset_name)
  if (!"weight" %in% edge_attr_names(g)) E(g)$weight <- 0.5

  nm          <- node_metrics(g, dataset_name)
  label_genes <- head(nm$gene, min(top_n_label, nrow(nm)))

  TG <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    left_join(nm, by = c("name" = "gene")) %>%
    mutate(col_bet = pmin(pmax(betweenness,
                               BET_LIMITS[1]), BET_LIMITS[2]))

  set.seed(seed)
  lay <- create_layout(TG, layout = "igraph", algorithm = "drl",
                       weights = if ("weight" %in% edge_attr_names(g))
                                   E(g)$weight else NULL,
                       options = DRL_OPTIONS)

  # Post-layout de-overlap
  node_size_range <- c(3, 16)
  node_sizes <- scales::rescale(lay$degree, to = node_size_range,
                                from = range(lay$degree, na.rm = TRUE))
  node_sizes[!is.finite(node_sizes)] <- mean(node_size_range)
  radii <- node_sizes * 0.04

  cr     <- packcircles::circleRepelLayout(
              data.frame(x = lay$x, y = lay$y, r = radii),
              xlim    = range(lay$x, na.rm = TRUE),
              ylim    = range(lay$y, na.rm = TRUE),
              maxiter = 4000, wrap = FALSE)
  lay$x  <- cr$layout$x
  lay$y  <- cr$layout$y

  edge_w_range <- if ("weight" %in% edge_attr_names(g)) c(0.4, 2.2) else c(0.6, 0.6)

  ggraph(lay) +
    geom_edge_link(aes(width = weight), alpha = 0.25,
                   colour = "grey40", show.legend = TRUE) +
    geom_node_point(aes(size = degree, color = col_bet),
                    alpha = 0.95, show.legend = TRUE) +
    geom_text_repel(
      data        = subset(lay, name %in% label_genes),
      aes(x = x, y = y, label = name),
      size        = 3.5, fontface = "bold", max.overlaps = Inf,
      box.padding = 0.5, point.padding = 0.4,
      force = 2, max.time = 2, min.segment.length = 0
    ) +
    scale_color_gradient(
      name   = "Betweenness\ncentrality",
      low    = "#f7f5bc", high  = "#025043",
      limits = BET_LIMITS, breaks = BET_BREAKS, oob = squish
    ) +
    scale_size_continuous(name  = "Degree",
                          range = node_size_range) +
    ggraph::scale_edge_width_continuous(
      name   = "STRING\nscore",
      range  = edge_w_range,
      limits = c(0, 1),
      breaks = pretty_breaks(4)
    ) +
    coord_equal(expand = TRUE, clip = "off") +
    theme_pub_void()
}

panel_G <- build_ppi_panel(STRING_FILES$epithelial, "Epithelial")
panel_H <- build_ppi_panel(STRING_FILES$immune,     "Immune")
panel_I <- build_ppi_panel(STRING_FILES$stromal,    "Stromal")


# =============================================================================
# STEP 10 – Assemble and export Figure 3 PDF
# =============================================================================

figure_3 <- (panel_A | panel_B | panel_C) /
            (panel_D | panel_E | panel_F) /
            (panel_G | panel_H | panel_I) +
  plot_annotation(tag_levels = list(c("A.", "B", "C.", "D", "E.", "F.",
                                      "G.", "H.", "I."))) &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(OUT_PDF,
       plot   = figure_3,
       width  = 22,
       height = 24,
       units  = "in",
       device = "pdf")

message("Done. Figure 3 saved to: ", OUT_PDF)
