# =============================================================================
# Figure 5 – Complete Pipeline
# Tissue-specific (unique) epithelial CCI analysis across intestine, skin,
# and uterus, plus a summary schematic panel (panel J).
#
# Figure layout:
#   Row 1: Venn diagrams highlighting tissue-unique CCIs (panels A–C)
#            - A: intestine-unique highlighted
#            - B: skin-unique highlighted
#            - C: uterus-unique highlighted
#   Row 2: Dot plots – top tissue-unique CCIs per epithelial compartment
#            - D: intestine-unique epithelial CCIs
#            - E: skin-unique (keratinocyte) CCIs
#            - F: uterus-unique (Epi-Ciliated) CCIs
#   Row 3: PPI hub-gene networks for tissue-unique CCIs (panels G–I)
#            - G: intestine-unique
#            - H: skin-unique
#            - I: uterus-unique
#   Row 4: Summary schematic – tissue-conserved vs tissue-specific CCIs (J)
#            (rendered as a ggplot annotation panel)
#
# Pipeline overview (Steps 1-6 are shared with Figure 3; run them first):
#   Step 1–6: same CCI inference, filtering, ranking as Figure 3 pipeline
#   Step 7:   Extract tissue-UNIQUE CCIs (present in only one tissue)
#   Step 8:   Venn diagrams with per-tissue highlighting (panels A–C)
#   Step 9:   Dot plots for unique CCIs per tissue (panels D–F)
#   Step 10:  PPI hub-gene network plots (panels G–I)
#   Step 11:  Summary schematic panel (panel J)
#   Step 12:  Assemble and export Figure_5.pdf
#
# Input files required (update paths in CONFIG below):
#   - Final ranked + filtered CCI CSVs per tissue × epithelial compartment
#     (outputs of Step 6 from Figure_3_pipeline.R)
#   - STRING interaction TSV files for each tissue-unique hub network
#
# Output:
#   - Figure_5.pdf
# =============================================================================


# =============================================================================
# 0. LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggvenn)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(readr)
  library(scales)
  library(ggrepel)
  library(packcircles)
  library(patchwork)
  library(grid)
  library(ggtext)        # for panel J rich-text labels
})


# =============================================================================
# 1. USER CONFIG  ← update all paths here
# =============================================================================

# --- Final ranked/filtered CCI CSVs (output of Figure 3 pipeline Step 6) ----
# Epithelial compartment, one file per tissue.
# Required columns: source, target, ligand, receptor, total_rbp_score, patient_id
EPI_CCI_FILES <- list(
  intestine = "path/to/final_filtered_epi_intestine.csv",
  skin      = "path/to/final_filtered_epi_skin.csv",
  uterus    = "path/to/final_filtered_epi_uterus.csv"
)

# --- STRING TSV files for tissue-unique hub networks (panels G–I) -----------
STRING_FILES_UNIQUE <- list(
  intestine = "path/to/hub_genes/intestine_unique_string_interactions.tsv",
  skin      = "path/to/hub_genes/skin_unique_string_interactions.tsv",
  uterus    = "path/to/hub_genes/uterus_unique_string_interactions.tsv"
)

# --- Venn diagram counts (same values as Figure 3, panels A–C) --------------
# These reproduce the published numbers; replace with your actual counts.
VENN_EPI <- list(
  intestine = 23, skin = 23, uterus = 28,
  int_skin  = 9,  int_uter = 5, skin_uter = 2, all = 7
)

# --- Output ------------------------------------------------------------------
OUT_PDF <- "Figure_5.pdf"

# RBP patience (should match Figure 3 pipeline)
RBP_P <- 0.8


# =============================================================================
# STEP 7 – Extract tissue-UNIQUE CCIs
# =============================================================================
# A CCI (ligand-receptor + source-target) is "tissue-unique" if it appears
# in exactly one tissue after the final filtering step.

extract_unique_ccis <- function(cci_file_list) {
  # Load all tissues
  combined <- map_dfr(names(cci_file_list), function(tiss) {
    read_csv(cci_file_list[[tiss]], show_col_types = FALSE) %>%
      mutate(tissue = tiss,
             LR     = paste(ligand, receptor, sep = "-"))
  })

  # Identify LR pairs that appear in only one tissue
  lr_tissue_counts <- combined %>%
    distinct(LR, tissue) %>%
    count(LR)

  unique_lrs <- lr_tissue_counts %>% filter(n == 1) %>% pull(LR)

  unique_ccis <- combined %>%
    filter(LR %in% unique_lrs)

  message("Tissue-unique LR pairs: ", length(unique_lrs))
  unique_ccis
}

unique_epi <- extract_unique_ccis(EPI_CCI_FILES)

# Split back by tissue for downstream use
unique_intestine <- unique_epi %>% filter(tissue == "intestine")
unique_skin      <- unique_epi %>% filter(tissue == "skin")
unique_uterus    <- unique_epi %>% filter(tissue == "uterus")


# =============================================================================
# STEP 8 – Venn diagrams with per-tissue highlighting (panels A, B, C)
# =============================================================================
# Each panel shows the same 3-set Venn but with a different tissue's circle
# highlighted (bold label + darker fill) to indicate the tissue-unique set.

build_venn_sets <- function(counts) {
  # Synthesise integer-ID element sets that reproduce the overlap counts
  mk <- local({ i <- 0L; function(k) { r <- i + seq_len(k); i <<- i + k; r } })

  shared_all <- mk(counts$all)
  shared_is  <- c(shared_all, mk(counts$int_skin  - counts$all))
  shared_iu  <- c(shared_all, mk(counts$int_uter  - counts$all))
  shared_su  <- c(shared_all, mk(counts$skin_uter - counts$all))
  only_int   <- mk(counts$intestine)
  only_skin  <- mk(counts$skin)
  only_uter  <- mk(counts$uterus)

  list(
    intestine = c(shared_all, shared_is, shared_iu, only_int),
    skin      = c(shared_all, shared_is, shared_su, only_skin),
    uterus    = c(shared_all, shared_iu, shared_su, only_uter)
  )
}

VENN_SETS <- build_venn_sets(VENN_EPI)

# Base fill colours (muted)
BASE_FILLS <- c(intestine = "#E5C4A1", skin = "#C8DFA4", uterus = "#A4C8DF")

draw_venn_highlight <- function(sets, highlight_tissue,
                                tissue_labels = c("epithelial intestine",
                                                  "epithelial skin",
                                                  "epithelial uterus")) {
  fills <- BASE_FILLS
  fills[highlight_tissue] <- colorspace::darken(fills[highlight_tissue], 0.25)

  named_sets         <- sets
  names(named_sets)  <- tissue_labels

  ggvenn(named_sets,
         fill_color    = unname(fills),
         stroke_size   = 0.6,
         set_name_size = 3.5,
         text_size     = 3.2) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

# Fall back to base ggvenn if colorspace is unavailable
if (!requireNamespace("colorspace", quietly = TRUE)) {
  draw_venn_highlight <- function(sets, highlight_tissue,
                                  tissue_labels = c("epithelial intestine",
                                                    "epithelial skin",
                                                    "epithelial uterus")) {
    named_sets        <- sets
    names(named_sets) <- tissue_labels
    fills             <- unname(BASE_FILLS)
    ggvenn(named_sets,
           fill_color    = fills,
           stroke_size   = 0.6,
           set_name_size = 3.5,
           text_size     = 3.2)
  }
}

panel_A <- draw_venn_highlight(VENN_SETS, "intestine")
panel_B <- draw_venn_highlight(VENN_SETS, "skin")
panel_C <- draw_venn_highlight(VENN_SETS, "uterus")


# =============================================================================
# STEP 9 – Dot plots for tissue-unique CCIs (panels D, E, F)
# =============================================================================

# Top LR pairs to display (from published figure; update as needed)
LR_UNIQUE_INTESTINE <- c(
  "CCL3-CCR1", "CCL5-CCR1", "CD31-PECAM1", "CEACAM1-CEACAM1",
  "DDR1-TTR", "DLL1-NOTCH1", "EFNB2-EPHB2", "FAM3C-PDC.D4",
  "FGFR1-PTPRR", "GRN-SORT1", "GRN-TNFRSF1B", "GUCA2A-GUCY2C",
  "HLA-E-KLRC1", "LGALS9-LRP1", "LGALS9-MRC1", "MDK-LRP1",
  "SCT-VIPR1", "TNFSF9-TNFRSFA", "VEGFa-EPHB2"
)

LR_UNIQUE_SKIN <- c(
  "ANXA1-FPR1", "AREG-EGFR", "AREG-ICAM1", "CCL2-CCR5",
  "CLECB4-KLRF1", "DSC1-DSG1", "EFNA1-EPHA1", "ENTPD1-ADORA2B",
  "GRN-EGFR", "HEBP1-FPR1", "HGF-CD44", "HGF-MET",
  "HLA-A-KIR3DL1", "HLA-F-KIR3DL1", "IL1B-ADRB2",
  "KLRF1-CLECB4", "LGALS9-MET", "SELE-CD44", "TGFB1-EGFR",
  "TIGIT-NECTIN1", "TNFSF14-LTBR", "VEGFa-NRP1"
)

LR_UNIQUE_UTERUS <- c(
  "ALOX5-ALOX5R", "ANXA1-FPR1", "BAG-NCR1", "CCL2-SLC7A1",
  "CD2-CD58", "CLECBD-KLRB1", "COPa-PRY2", "CSF1-CSF1R",
  "CXCL1-DPP4", "EDN1-EDNRA", "EREG-ERBB2", "FAM3C-LAMP1",
  "FASLG-TNFRSF6B", "GAL1-HLA-DPa1", "IGF1-IGF1R", "IGF2-IGF1R",
  "JAG1-NOTCH1", "JAG2-NOTCH1", "JAG3-NOTCH1", "JAG4-NOTCH2",
  "LGALS9-CD44", "SEMA4F-NRP1", "SPP1-CD44", "SPP1-PTGER4",
  "TNF-CELSR1", "TNF-NOTCH1", "TNF-PTPRS", "TNF-RIPK"
)

TISSUE_COLS <- c(intestine = "#E5835A", skin = "#82B0D9", uterus = "#7DC99E")

build_unique_dotplot <- function(data_tiss, lr_of_interest, tissue_col) {
  plot_data <- data_tiss %>%
    filter(LR %in% lr_of_interest) %>%
    mutate(
      source_target = paste(source, target, sep = "-"),
      LR            = factor(LR, levels = rev(sort(unique(LR[LR %in% lr_of_interest]))))
    )

  ggplot(plot_data, aes(x = LR, y = source_target)) +
    geom_point(aes(size = total_rbp_score), colour = tissue_col, alpha = 0.85) +
    scale_size_continuous(
      name   = "CIS",
      range  = c(1, 6),
      breaks = seq(0.5, 2.5, 0.5)
    ) +
    coord_flip() +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
      axis.text.y  = element_text(size = 6),
      panel.grid   = element_line(colour = "grey92"),
      legend.key.size = unit(0.35, "cm"),
      legend.text  = element_text(size = 7),
      legend.title = element_text(size = 7)
    ) +
    labs(x = "CCIs", y = "source–target cells")
}

panel_D <- build_unique_dotplot(unique_intestine, LR_UNIQUE_INTESTINE,
                                 TISSUE_COLS["intestine"])
panel_E <- build_unique_dotplot(unique_skin,      LR_UNIQUE_SKIN,
                                 TISSUE_COLS["skin"])
panel_F <- build_unique_dotplot(unique_uterus,    LR_UNIQUE_UTERUS,
                                 TISSUE_COLS["uterus"])


# =============================================================================
# STEP 10 – PPI hub-gene network plots (panels G, H, I)
# =============================================================================
# Re-uses the same DRL-layout helper functions from Figure 3 pipeline.

BET_LIMITS  <- c(0, 0.5)
BET_BREAKS  <- seq(0, 0.5, by = 0.1)
DRL_OPTIONS <- list(init.iterations = 300L, maxiter = 2000L, cool.factor = 0.95)

read_edges_string <- function(path) {
  tb       <- read_tsv(path, show_col_types = FALSE)
  cn       <- names(tb); lc <- tolower(cn); map <- setNames(cn, lc)
  p1       <- map[["protein1"]]; p2 <- map[["protein2"]]
  if (is.na(p1) || is.na(p2)) {
    p1 <- map[["preferredname_a"]]; p2 <- map[["preferredname_b"]]
  }
  scorecol <- map[["combined_score"]]
  if (is.na(scorecol)) scorecol <- map[["combinedscore"]]
  stopifnot(!is.na(p1), !is.na(p2))

  ed <- tb %>%
    transmute(
      protein1       = .data[[p1]],
      protein2       = .data[[p2]],
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
    ed <- ed %>% mutate(distance = ifelse(combined_score > 0, 1/combined_score, Inf))
    graph_from_data_frame(
      ed %>% select(protein1, protein2, weight = combined_score, distance),
      directed = FALSE)
  } else {
    graph_from_data_frame(ed %>% select(protein1, protein2), directed = FALSE)
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
    weights = if ("weight" %in% edge_attr_names(g)) E(g)$weight else NULL)$vector
  tibble(gene = names(deg),
         degree      = as.numeric(deg),
         betweenness = as.numeric(bet),
         eigenvector = as.numeric(eig)) %>%
    mutate(r_deg  = rank(-degree,      ties.method = "min"),
           r_bet  = rank(-betweenness, ties.method = "min"),
           r_eig  = rank(-eigenvector, ties.method = "min"),
           mean_rank      = rowMeans(cbind(r_deg, r_bet, r_eig), na.rm = TRUE),
           consensus_rank = rank(mean_rank, ties.method = "min"),
           dataset        = dataset_name) %>%
    arrange(consensus_rank, desc(degree))
}

theme_pub_void <- function(base_size = 11) {
  theme_void(base_size = base_size) +
    theme(plot.title   = element_text(hjust = 0.5, face = "bold"),
          legend.title = element_text(face = "bold"),
          plot.margin  = margin(6, 6, 6, 6))
}

build_ppi_panel <- function(path, dataset_name, top_n_label = 10, seed = 7L) {
  ed <- read_edges_string(path)
  stopifnot(nrow(ed) > 0)

  g <- graph_from_edges(ed)
  if ("weight" %in% edge_attr_names(g)) {
    w <- scales::rescale(sqrt(pmax(E(g)$weight, 1e-6)), to = c(0.15, 0.75))
    E(g)$weight <- w
  }

  g <- prune_to_lcc(g, min_comp_size = 2)
  if (gorder(g) == 0 || gsize(g) == 0)
    stop("Empty graph after pruning for: ", dataset_name)
  if (!"weight" %in% edge_attr_names(g)) E(g)$weight <- 0.5

  nm          <- node_metrics(g, dataset_name)
  label_genes <- head(nm$gene, min(top_n_label, nrow(nm)))

  TG <- as_tbl_graph(g) %>%
    activate(nodes) %>%
    left_join(nm, by = c("name" = "gene")) %>%
    mutate(col_bet = pmin(pmax(betweenness, BET_LIMITS[1]), BET_LIMITS[2]))

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
    xlim = range(lay$x, na.rm = TRUE), ylim = range(lay$y, na.rm = TRUE),
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
      size = 3.5, fontface = "bold", max.overlaps = Inf,
      box.padding = 0.5, point.padding = 0.4,
      force = 2, max.time = 2, min.segment.length = 0
    ) +
    scale_color_gradient(
      name   = "Betweenness\ncentrality",
      low    = "#f7f5bc", high = "#025043",
      limits = BET_LIMITS, breaks = BET_BREAKS, oob = squish
    ) +
    scale_size_continuous(name = "Degree", range = node_size_range) +
    ggraph::scale_edge_width_continuous(
      name   = "STRING\nscore",
      range  = edge_w_range, limits = c(0, 1), breaks = pretty_breaks(4)
    ) +
    coord_equal(expand = TRUE, clip = "off") +
    theme_pub_void()
}

panel_G <- build_ppi_panel(STRING_FILES_UNIQUE$intestine, "Intestine-unique")
panel_H <- build_ppi_panel(STRING_FILES_UNIQUE$skin,      "Skin-unique")
panel_I <- build_ppi_panel(STRING_FILES_UNIQUE$uterus,    "Uterus-unique")


# =============================================================================
# STEP 11 – Summary schematic panel (panel J)
# =============================================================================
# Panel J is a conceptual diagram. It is reproduced here as a ggplot
# annotation layout. Replace with your own illustration or import an image
# if you have a pre-made graphic (e.g. a PNG exported from BioRender).
#
# To import a pre-made image instead:
#   library(png); library(grid)
#   img <- readPNG("path/to/panel_J.png")
#   panel_J <- ggplot() +
#     annotation_custom(rasterGrob(img, interpolate = TRUE)) +
#     theme_void()

make_schematic_panel <- function() {

  # --- Layout data -----------------------------------------------------------
  conserved <- tibble(
    xmin = c(0.02, 0.35, 0.68),
    xmax = c(0.32, 0.65, 0.98),
    ymin = 0.08, ymax = 0.72,
    fill  = c("#E5C4A1", "#C8DFA4", "#A4C8DF"),
    cell1 = c("Epi-Intestine", "Keratinocytes", "Epi-Ciliated"),
    cell2 = c("DC",            "Langerhans",    "Myeloid"),
    func  = c("Epithelial\nrestitution", "Wound\nhealing", "Tissue\nremodeling")
  )

  specific <- tibble(
    xmin = c(0.02, 0.35, 0.68),
    xmax = c(0.32, 0.65, 0.98),
    ymin = 0.08, ymax = 0.72,
    fill  = c("#E5C4A1", "#C8DFA4", "#A4C8DF"),
    lr    = c("GUCA2A/GUCA2B-\nGUCY2C", "HLA-A/HLA-F-\nKIR3DL1", "SPP1-PTGER4"),
    cell1 = c("Epi-Intestine -\nEnteroendocrine\nL cells",
              "Keratinocytes -\nNK cells",
              "Epi-Ciliated -\nMyeloid"),
    func  = c("Mediates epithelial-\nendocrine crosstalk",
              "Keratinocytes limit\nNK cell cytotoxicity",
              "Epithelial-myeloid\ncrosstalk in resolving\ninflammation")
  )

  # Build the plot with two sections side by side
  p <- ggplot() +
    # ---- outer box ----
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "white", colour = "grey30", linewidth = 0.6) +

    # ---- section headers ----
    annotate("text", x = 0.25, y = 0.93,
             label = "Tissue-conserved CCIs",
             fontface = "bold", size = 4.5, hjust = 0.5) +
    annotate("text", x = 0.75, y = 0.93,
             label = "Tissue-specific CCIs",
             fontface = "bold", size = 4.5, hjust = 0.5) +
    annotate("segment", x = 0.5, xend = 0.5, y = 0.02, yend = 0.98,
             colour = "grey50", linewidth = 0.5, linetype = "dashed") +

    # ---- conserved: MIF-CD74 banner ----
    annotate("rect", xmin = 0.01, xmax = 0.49, ymin = 0.78, ymax = 0.90,
             fill = "#F0EAD6", colour = "grey60", linewidth = 0.4) +
    annotate("text", x = 0.25, y = 0.84,
             label = "MIF - CD74", fontface = "bold.italic", size = 4) +

    # ---- conserved: three tissue boxes ----
    {
      map(seq_len(nrow(conserved)), function(i) {
        r <- conserved[i, ]
        list(
          annotate("rect",
                   xmin = r$xmin * 0.48 + 0.01,
                   xmax = r$xmax * 0.48 + 0.01,
                   ymin = r$ymin, ymax = r$ymax,
                   fill = r$fill, colour = "grey70", linewidth = 0.3),
          annotate("text",
                   x = ((r$xmin + r$xmax) / 2) * 0.48 + 0.01,
                   y = 0.55,
                   label = paste0(r$cell1, " –\n", r$cell2),
                   size = 2.8, hjust = 0.5),
          annotate("text",
                   x = ((r$xmin + r$xmax) / 2) * 0.48 + 0.01,
                   y = 0.18,
                   label = r$func,
                   size = 2.6, hjust = 0.5, colour = "grey20",
                   fontface = "italic")
        )
      })
    } +

    # ---- specific: three tissue boxes ----
    {
      map(seq_len(nrow(specific)), function(i) {
        r <- specific[i, ]
        list(
          annotate("rect",
                   xmin = r$xmin * 0.48 + 0.51,
                   xmax = r$xmax * 0.48 + 0.51,
                   ymin = r$ymin, ymax = r$ymax,
                   fill = r$fill, colour = "grey70", linewidth = 0.3),
          annotate("text",
                   x = ((r$xmin + r$xmax) / 2) * 0.48 + 0.51,
                   y = 0.75,
                   label = r$lr,
                   size  = 2.5, hjust = 0.5, fontface = "bold.italic"),
          annotate("text",
                   x = ((r$xmin + r$xmax) / 2) * 0.48 + 0.51,
                   y = 0.45,
                   label = r$cell1,
                   size  = 2.6, hjust = 0.5),
          annotate("text",
                   x = ((r$xmin + r$xmax) / 2) * 0.48 + 0.51,
                   y = 0.18,
                   label = r$func,
                   size  = 2.5, hjust = 0.5, colour = "grey20",
                   fontface = "italic")
        )
      })
    } +

    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(4, 4, 4, 4))

  p
}

panel_J <- make_schematic_panel()


# =============================================================================
# STEP 12 – Assemble and export Figure 5 PDF
# =============================================================================

# Rows 1–3: 3-column panels (A–I)
top_row    <- panel_A | panel_B | panel_C
middle_row <- panel_D | panel_E | panel_F
network_row <- panel_G | panel_H | panel_I

# Row 4: panel J spans full width
figure_5 <- (top_row / middle_row / network_row / panel_J) +
  plot_layout(heights = c(1.1, 2.2, 2.2, 1.5)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C",
                        "D", "E.", "F.",
                        "G", "H", "I.",
                        "J."))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(OUT_PDF,
       plot   = figure_5,
       width  = 22,
       height = 28,
       units  = "in",
       device = "pdf")

message("Done. Figure 5 saved to: ", OUT_PDF)
