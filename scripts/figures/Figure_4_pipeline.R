# =============================================================================
# Figure 4 – MIF & CD74 expression validation
#
# Panels:
#   A. MIF expression violin plot (GTEx nTPM, pre vs post-menopause)
#   B. CD74 expression violin plot (GTEx nTPM, pre vs post-menopause)
#   C. % MIF-positive cells per donor – bar chart (HPA data)
#   D. MIF IHC images (3 donors × 3 ROIs) – imported from file
#   E. % CD74-positive cells per donor – bar chart (HPA data)
#   F. CD74 IHC images (2 donors × 3 ROIs) – imported from file
#
# INPUT FILES (update paths in the CONFIG section):
#   GTEx CSVs  – one per gene × tissue × menopause status; must contain nTPM
#   HPA CSVs   – one per gene; must contain: donor, cell_type, pct_positive, sd
#   IHC images – PNG/JPG files for panels D and F
#
# OUTPUT: Figure_4.pdf
# =============================================================================


# =============================================================================
# 0. LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(ggpattern)   # install.packages("ggpattern") – for hatched violins
  library(png)         # for importing IHC images
  library(jpeg)
  library(grid)
})


# =============================================================================
# 1. USER CONFIG  ← update all paths here
# =============================================================================

# --- GTEx CSVs (panels A & B) ------------------------------------------------
# Each file must have at minimum an `nTPM` column.
GTEX_MIF <- list(
  uterus_pre   = "path/to/MIF_endometrium_pre.csv",
  uterus_post  = "path/to/MIF_endometrium_post.csv",
  skin_pre     = "path/to/MIF_skin_pre.csv",
  skin_post    = "path/to/MIF_skin_post.csv",
  intestine_pre  = "path/to/MIF_small_intestine_pre.csv",
  intestine_post = "path/to/MIF_small_intestine_post.csv"
)

GTEX_CD74 <- list(
  uterus_pre   = "path/to/CD74_endometrium_pre.csv",
  uterus_post  = "path/to/CD74_endometrium_post.csv",
  skin_pre     = "path/to/CD74_skin_pre.csv",
  skin_post    = "path/to/CD74_skin_post.csv",
  intestine_pre  = "path/to/CD74_small_intestine_pre.csv",
  intestine_post = "path/to/CD74_small_intestine_post.csv"
)

# --- HPA bar chart CSVs (panels C & E) ---------------------------------------
# Required columns: donor (character), cell_type (character),
#                   pct_positive (numeric 0-100), sd (numeric)
HPA_MIF_CSV  <- "path/to/HPA_MIF_positive_cells.csv"
HPA_CD74_CSV <- "path/to/HPA_CD74_positive_cells.csv"

# Cell-type display labels and fill colours
MIF_CELLTYPES <- c(
  "Glandular cells of intestine" = "#E5835A",
  "Keratinocytes of skin"        = "#82B0D9",
  "Glandular cells of uterus"    = "#7DC99E"
)

CD74_CELLTYPES <- c(
  "Glandular cells of intestine" = "#E5835A",
  "Glandular cells of uterus"    = "#7DC99E"
)

# --- IHC image files (panels D & F) -----------------------------------------
# Supply paths to PNG or JPG images.
# Panel D: MIF IHC – 3 donors × 3 ROIs (9 images, row-major order)
IHC_MIF_IMAGES <- list(
  "donor 3 (intestine)" = c("path/to/MIF_donor3_ROI1.png",
                             "path/to/MIF_donor3_ROI2.png",
                             "path/to/MIF_donor3_ROI3.png"),
  "donor 4 (skin)"      = c("path/to/MIF_donor4_ROI1.png",
                             "path/to/MIF_donor4_ROI2.png",
                             "path/to/MIF_donor4_ROI3.png"),
  "donor 7 (endometrium)" = c("path/to/MIF_donor7_ROI1.png",
                               "path/to/MIF_donor7_ROI2.png",
                               "path/to/MIF_donor7_ROI3.png")
)

# Panel F: CD74 IHC – 2 donors × 3 ROIs (6 images)
IHC_CD74_IMAGES <- list(
  "donor 2 (intestine)"   = c("path/to/CD74_donor2_ROI1.png",
                               "path/to/CD74_donor2_ROI2.png",
                               "path/to/CD74_donor2_ROI3.png"),
  "donor 6 (endometrium)" = c("path/to/CD74_donor6_ROI1.png",
                               "path/to/CD74_donor6_ROI2.png",
                               "path/to/CD74_donor6_ROI3.png")
)

# --- Output ------------------------------------------------------------------
OUT_PDF <- "Figure_4.pdf"


# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

# Load a GTEx CSV list and bind into one data frame
load_gtex <- function(file_list) {
  map_dfr(names(file_list), function(key) {
    parts  <- str_split(key, "_", n = 2)[[1]]
    tissue <- str_to_title(parts[1])
    status <- ifelse(parts[2] == "pre", "Premenopause", "Postmenopause")
    read_csv(file_list[[key]], show_col_types = FALSE) %>%
      mutate(tissue = factor(tissue,
                             levels = c("Intestine", "Skin", "Uterus")),
             status = factor(status,
                             levels = c("Premenopause", "Postmenopause")))
  })
}

# Read a PNG or JPEG and return a rasterGrob
load_image_grob <- function(path) {
  ext <- tolower(tools::file_ext(path))
  img <- if (ext == "png") readPNG(path) else readJPEG(path)
  rasterGrob(img, interpolate = TRUE)
}

# Assemble a grid of IHC images into a single patchwork-compatible ggplot
build_ihc_panel <- function(image_list,
                             col_labels = c("ROI 1", "ROI 2", "ROI 3")) {
  donor_names <- names(image_list)
  n_rows <- length(donor_names)
  n_cols <- length(col_labels)

  # One ggplot per cell
  plots <- vector("list", n_rows * n_cols)
  k <- 1L
  for (i in seq_len(n_rows)) {
    for (j in seq_len(n_cols)) {
      path <- image_list[[i]][j]
      grob <- load_image_grob(path)
      p <- ggplot() +
        annotation_custom(grob,
                          xmin = -Inf, xmax = Inf,
                          ymin = -Inf, ymax = Inf) +
        theme_void() +
        theme(
          plot.margin = margin(1, 1, 1, 1),
          # Row label on leftmost cell only
          axis.title.y = if (j == 1)
            element_text(size = 8, angle = 90,
                         hjust = 0.5, face = "bold",
                         margin = margin(r = 4))
          else element_blank()
        ) +
        # Column header on top row only
        { if (i == 1) labs(title = col_labels[j]) else labs(title = NULL) } +
        { if (j == 1) labs(y = donor_names[i]) else labs(y = NULL) }

      if (i == 1) {
        p <- p + theme(plot.title = element_text(size = 8, hjust = 0.5,
                                                  face = "bold"))
      }
      plots[[k]] <- p
      k <- k + 1L
    }
  }

  wrap_plots(plots, nrow = n_rows, ncol = n_cols) &
    theme(plot.margin = margin(1, 1, 1, 1))
}


# =============================================================================
# 3. PANEL A – MIF violin plot
# =============================================================================

mif_data <- load_gtex(GTEX_MIF)

# Tissue fill colours (solid for Premenopause; hatched pattern for Postmenopause)
tissue_fills <- c(intestine = "#E5835A", skin = "#82B0D9", uterus = "#7DC99E")

panel_A <- ggplot(mif_data,
                  aes(x = tissue, y = nTPM,
                      fill    = tissue,
                      pattern = status)) +
  geom_violin_pattern(
    aes(pattern_density = status),
    position     = position_dodge(width = 0.75),
    width        = 0.75,
    trim         = FALSE,
    colour       = "black",
    linewidth    = 0.35,
    pattern_fill = "black",
    pattern_colour = NA,
    pattern_spacing = 0.03
  ) +
  geom_boxplot(
    aes(group = interaction(tissue, status)),
    position     = position_dodge(width = 0.75),
    width        = 0.12,
    fill         = "white",
    colour       = "black",
    linewidth    = 0.35,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    name   = "tissue",
    values = c(Intestine = "#E5835A", Skin = "#82B0D9", Uterus = "#7DC99E")
  ) +
  scale_pattern_manual(
    name   = "status",
    values = c(Premenopause = "none", Postmenopause = "crosshatch")
  ) +
  scale_pattern_density_manual(
    name   = "status",
    values = c(Premenopause = 0, Postmenopause = 0.15)
  ) +
  scale_x_discrete(labels = tolower) +
  labs(title = "MIF expression", y = "nTPM", x = NULL) +
  guides(
    fill    = guide_legend(override.aes = list(pattern = "none")),
    pattern = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 8, face = "bold"),
    axis.text.x     = element_text(size = 9)
  )


# =============================================================================
# 4. PANEL B – CD74 violin plot
# =============================================================================

cd74_data <- load_gtex(GTEX_CD74)

panel_B <- ggplot(cd74_data,
                  aes(x = tissue, y = nTPM,
                      fill    = tissue,
                      pattern = status)) +
  geom_violin_pattern(
    aes(pattern_density = status),
    position       = position_dodge(width = 0.75),
    width          = 0.75,
    trim           = FALSE,
    colour         = "black",
    linewidth      = 0.35,
    pattern_fill   = "black",
    pattern_colour = NA,
    pattern_spacing = 0.03
  ) +
  geom_boxplot(
    aes(group = interaction(tissue, status)),
    position      = position_dodge(width = 0.75),
    width         = 0.12,
    fill          = "white",
    colour        = "black",
    linewidth     = 0.35,
    outlier.shape = NA
  ) +
  scale_fill_manual(
    name   = "tissue",
    values = c(Intestine = "#E5835A", Skin = "#82B0D9", Uterus = "#7DC99E")
  ) +
  scale_pattern_manual(
    name   = "status",
    values = c(Premenopause = "none", Postmenopause = "crosshatch")
  ) +
  scale_pattern_density_manual(
    name   = "status",
    values = c(Premenopause = 0, Postmenopause = 0.15)
  ) +
  scale_x_discrete(labels = tolower) +
  labs(title = "CD74 expression", y = "nTPM", x = NULL) +
  guides(
    fill    = guide_legend(override.aes = list(pattern = "none")),
    pattern = guide_legend(override.aes = list(fill = "white"))
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.key.size = unit(0.5, "cm"),
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 8, face = "bold"),
    axis.text.x     = element_text(size = 9)
  )


# =============================================================================
# 5. PANEL C – MIF % positive cells bar chart (HPA)
# =============================================================================

hpa_mif <- read_csv(HPA_MIF_CSV, show_col_types = FALSE) %>%
  mutate(
    donor     = factor(donor,
                       levels = unique(donor)),   # preserves CSV order
    cell_type = factor(cell_type, levels = names(MIF_CELLTYPES))
  )

panel_C <- ggplot(hpa_mif,
                  aes(x = donor, y = pct_positive, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.7),
           width    = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = pct_positive - sd,
        ymax = pct_positive + sd),
    position = position_dodge(width = 0.7),
    width    = 0.25, linewidth = 0.4
  ) +
  scale_fill_manual(values = MIF_CELLTYPES, name = NULL) +
  scale_y_continuous(limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0)) +
  labs(
    title = "% of MIF positive cells per donor (HPA)",
    x     = "Donor",
    y     = "% of positive cells (mean ± SD)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  )


# =============================================================================
# 6. PANEL D – MIF IHC images
# =============================================================================

panel_D <- build_ihc_panel(IHC_MIF_IMAGES)


# =============================================================================
# 7. PANEL E – CD74 % positive cells bar chart (HPA)
# =============================================================================

hpa_cd74 <- read_csv(HPA_CD74_CSV, show_col_types = FALSE) %>%
  mutate(
    donor     = factor(donor, levels = unique(donor)),
    cell_type = factor(cell_type, levels = names(CD74_CELLTYPES))
  )

panel_E <- ggplot(hpa_cd74,
                  aes(x = donor, y = pct_positive, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.7),
           width    = 0.65, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = pct_positive - sd,
        ymax = pct_positive + sd),
    position = position_dodge(width = 0.7),
    width    = 0.25, linewidth = 0.4
  ) +
  scale_fill_manual(values = CD74_CELLTYPES, name = NULL) +
  scale_y_continuous(limits = c(0, 110),
                     breaks = c(0, 25, 50, 75, 100),
                     expand = c(0, 0)) +
  labs(
    title = "% of CD74 positive cells per donor (HPA)",
    x     = "Donor",
    y     = "% of positive cells (mean ± SD)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.45, "cm")
  )


# =============================================================================
# 8. PANEL F – CD74 IHC images
# =============================================================================

panel_F <- build_ihc_panel(IHC_CD74_IMAGES)


# =============================================================================
# 9. ASSEMBLE & EXPORT Figure 4
# =============================================================================

# Row 1: A (left) | B (right)
row1 <- panel_A | panel_B

# Row 2: C (left, narrower) | D (right, wider – 3×3 image grid)
row2 <- panel_C | panel_D

# Row 3: E (left) | F (right)
row3 <- panel_E | panel_F

figure_4 <- row1 / row2 / row3 +
  plot_layout(heights = c(1.2, 1.5, 1.3)) +
  plot_annotation(
    tag_levels = list(c("A.", "B.", "C.", "D.", "E.", "F."))
  ) &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave(OUT_PDF,
       plot   = figure_4,
       width  = 16,
       height = 20,
       units  = "in",
       device = "pdf")

message("Done. Figure 4 saved to: ", OUT_PDF)


# =============================================================================
# APPENDIX – Expected CSV formats
# =============================================================================
#
# GTEx files (panels A & B):
#   nTPM
#   123.4
#   56.7
#   ...
#
# HPA MIF file (panel C):
#   donor,cell_type,pct_positive,sd
#   donor 1,Glandular cells of intestine,26,3.1
#   donor 2,Glandular cells of intestine,93,2.0
#   donor 3,Glandular cells of intestine,5,0.8
#   donor 4,Keratinocytes of skin,16,4.2
#   ...
#
# HPA CD74 file (panel E):
#   donor,cell_type,pct_positive,sd
#   donor 1,Glandular cells of intestine,91,2.3
#   ...
