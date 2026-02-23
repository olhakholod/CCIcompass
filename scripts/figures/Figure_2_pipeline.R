###############################################################################
##########            Integration for epithelial barriers            ##########
###############################################################################
library(Seurat)
library(ggplot2)
library(patchwork)

# Load annotated Seurat objects
intestine <- readRDS("/Users/f006qpk/Desktop/PLOS\ Comp\ Biol\ paper/data/intestine/02_female_intestine_annotated_processed.rds")
skin <- readRDS("/Users/f006qpk/Desktop/PLOS\ Comp\ Biol\ paper/data/skin/02_female_skin_annotated_processed.rds")
uterus <- readRDS("/Users/f006qpk/Desktop/PLOS\ Comp\ Biol\ paper/data/uterus/02_female_uterus_annotated_processed.rds")

# Example to merge more than two Seurat objects
merged_obj <- merge(x = intestine, y = list(skin, uterus))

# this object resides on Terra in Olha_Playground directory
saveRDS(merged_obj, "merged_obj_RO3.rds")

###############################################################################
# Use Terra for this part
# Set environment variables in a Jupyter Notebook
project <- Sys.getenv('WORKSPACE_NAMESPACE')
workspace <- Sys.getenv('WORKSPACE_NAME')
bucket <- Sys.getenv('WORKSPACE_BUCKET')

# Loads some required packages
library(Seurat)
library(dplyr)
library(glmGamPoi)
library(harmony)

system("mkdir $(pwd)/RO3")
system("gsutil cp -r gs://fc-bc13c4b2-a50a-47ff-aeb6-1261785eda08/RO3/merged_obj_RO3.rds ./RO3")
list.files(path = "./RO3")

# Loads the saved .rds file into a Seurat object
merged_obj <- readRDS("RO3/merged_obj_RO3.rds")

# Run Seurat pipeline
merged_obj <- SCTransform(merged_obj, vst.flavor = "v2")
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

# Run integration with Harmony
merged_obj <- RunHarmony(merged_obj, c("donor_id"))
merged_obj <- RunUMAP(object = merged_obj, dims = 1:30)

# Dimensional reduction plot
DimPlot(object = merged_obj, reduction = "umap", group.by = "cell_type")

# Save .rds file
saveRDS(merged_obj,file="integrated_RO3_dataset.rds")

# Copy all files generated in the notebook into the bucket
system(paste0("gsutil cp integrated_RO3_dataset.rds gs://fc-bc13c4b2-a50a-47ff-aeb6-1261785eda08/RO3/"),intern=TRUE)

# Run list command to see if file is in the bucket
system(paste0("gsutil ls ",bucket),intern=TRUE)

###############################################################################
# Start with Terra output to harmonize cell_type
merged_obj <- readRDS("/Users/f006qpk/Desktop/PLOS\ Comp\ Biol\ paper/data/integrated/integrated_annotated_RO3_dataset.rds")

# Make some levels in metadata as a factor
# cell_type
merged_obj@meta.data$cell_type <- factor(merged_obj@meta.data$cell_type)
levels(merged_obj@meta.data$cell_type)

# Harmonize levels in cell_type
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="activated CD4-positive, alpha-beta T cell"] <- "T-CD4"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="B cell"] <- "B"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="CD4-positive, alpha-beta T cell"] <- "T-CD4"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="CD8-positive, alpha-beta T cell"] <- "T-CD8"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="colon epithelial cell"] <- "Epi-Colon"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="cytotoxic T cell"] <- "CTL"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="endothelial cell of artery"] <- "Endo"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="endothelial cell of uterus"] <- "Endo"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="enterocyte of colon"] <- "Enterocyte"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="fibroblast"] <- "Fibro"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="glandular epithelial cell"] <- "Epi-Glandular"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="group 3 innate lymphoid cell"] <- "ILC3"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="IgA plasma cell"] <- "Plasma"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="inflammatory macrophage"] <- "M1 Mac"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="intestinal tuft cell"] <- "Tuft"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="keratinocyte"] <- "Keratinocyte"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="leukocyte"] <- "Leukocyte"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="macrophage"] <- "Mac"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="melanocyte"] <- "Melanocyte"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="monocyte"] <- "Mono"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="mucosal invariant T cell"] <- "T-Mucosal"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="myofibroblast cell"] <- "Myofibro"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="natural killer cell"] <- "NK"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="plasma cell"] <- "Plasma"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="regulatory T cell"] <- "T-reg"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="Schwann cell"] <- "Schwann"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="stromal cell of ovary"] <- "Stromal-Ovary"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="T-helper 1 cell"] <- "Th1"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="transit amplifying cell"] <- "TAC"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="type EC enteroendocrine cell"] <- "EC cells"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="type N enteroendocrine cell"] <- "N cells"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="uterine smooth muscle cell"] <- "SMC-Uterus"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="activated CD8-positive, alpha-beta T cell"] <- "T-CD8"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="capillary endothelial cell"] <- "Endo"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="CD8-positive, alpha-beta memory T cell"] <- "T-CD8"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="ciliated cell"] <- "Epi-Ciliated"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="conventional dendritic cell"] <- "DC"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="dendritic cell"] <- "DC"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="epithelial cell"] <- "Epi"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="endothelial cell of vascular tree"] <- "Endo"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="gamma-delta T cell"] <- "T-gamma/delta"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="glial cell"] <- "Glial"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="helper T cell"] <- "Th"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="IgG plasma cell"] <- "Plasma"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="innate lymphoid cell"] <- "ILC3"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="intestine goblet cell"] <- "Goblet"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="Langerhans cell"] <- "Langerhans"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="M cell of gut"] <- "M cells"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="mast cell"] <- "Mast"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="memory B cell"] <- "B-mem"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="monocyte-derived dendritic cell"] <- "mo-DC"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="naive B cell"] <- "B-naïve"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="pericyte"] <- "Pericyte"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="plasmacytoid dendritic cell"] <- "DC"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="reticular cell"] <- "Reticular cell"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="skin fibroblast"] <- "Fibro-skin"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="stromal cell"] <- "Stromal"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="supporting cell"] <- "Supporting cell"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="T-helper 17 cell"] <- "Th17"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="type D enteroendocrine cell"] <- "D cells"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="type L enteroendocrine cell"] <- "L cells"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="unknown"] <- "Unknown"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="vein endothelial cell"] <- "Endo"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="stem cell"] <- "Stem cell"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="endothelial cell of lymphatic vessel"] <- "Endo"
levels(merged_obj@meta.data$cell_type)[levels(merged_obj@meta.data$cell_type)=="myeloid cell"] <- "Myeloid"

DimPlot(object = merged_obj, reduction = "umap", group.by = "cell_type",
        raster = TRUE, raster.dpi = c(800, 800)) + theme_void() +
  theme(legend.position = "right",
        plot.title = element_blank())

# Save annotated object
#saveRDS(merged_obj, "integrated_annotated_RO3_dataset.rds")

###############################################################################
# Figure 2A: Pie charts for Donor, Sample, and # of cells by organ
library(ggplot2)

# Define organ colors (intestine = orange, skin = blue/purple, uterus = green)
organ_colors_pie <- c("intestine" = "#fc8d62", "skin" = "#8da0cb", "uterus" = "#66c2a5")

# Add organ column based on donor_id if not already present
donor_organ_map <- list(
  "uterus"    = c("A13", "A30", "E1", "E2", "E3", "SAMN15049042", "SAMN15049043",
                  "SAMN15049044", "SAMN15049045", "SAMN15049046", "SAMN15049047",
                  "SAMN15049048", "SAMN15049049", "SAMN15049050", "SAMN15049051"),
  "intestine" = c("A26 (386C)", "A30 (398B)", "A38 (432C)"),
  "skin"      = c("S1", "S2", "S3", "S4", "S5")
)

merged_obj@meta.data$organ <- sapply(merged_obj@meta.data$donor_id, function(x) {
  org <- names(donor_organ_map)[sapply(donor_organ_map, function(y) x %in% y)]
  if (length(org) > 0) org else NA
})

meta <- merged_obj@meta.data

# Helper to make a pie chart
make_pie <- function(df_col, title) {
  counts <- as.data.frame(table(df_col))
  colnames(counts) <- c("organ", "n")
  counts$organ <- as.character(counts$organ)
  ggplot(counts, aes(x = "", y = n, fill = organ)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y") +
    scale_fill_manual(values = organ_colors_pie,
                      labels = c("intestine", "skin", "uterus"),
                      name = title) +
    theme_void() +
    theme(legend.title = element_text(size = 12),
          legend.text  = element_text(size = 11))
}

# Pie 1: Number of unique donors per organ
donor_organ_df <- unique(meta[, c("donor_id", "organ")])
pie_donor <- make_pie(donor_organ_df$organ, "Donor")

# Pie 2: Number of unique samples per organ (using tissue/sample_id column if available)
# Using tissue as proxy for sample
pie_sample <- make_pie(meta$organ, "Sample")

# Pie 3: Number of cells per organ
pie_cells <- make_pie(meta$organ, "# of cells")

# Combine the three pie charts (Panel A)
library(patchwork)
panel_A <- pie_donor + pie_sample + pie_cells + plot_layout(nrow = 1)
print(panel_A)

###############################################################################
# Figure 2D: UMAP colored by cell_type
# UMAP based on cell_type variable
library(ggplot2)

# Define the cell type categories and levels
cell_type_categories <- list(
  "B cells" = c("B", "B-mem", "B-naïve"),
  "Endothelial" = c("Endo"),
  "Enteroendocrine" = c("EC cells", "N cells", "M cells", "D cells", "L cells"),
  "Epithelial" = c("Epi-Colon", "Enterocyte", "Epi-Glandular", "Tuft", "Keratinocyte", 
                   "Epi-Ciliated", "Epi", "Goblet"),
  "Fibroblasts" = c("Fibro", "Myofibro", "Fibro-skin"),
  "Lymphoid cells" = c("ILC3", "NK", "Plasma"),
  "Myeloid cells" = c("M1 Mac", "Leukocyte", "Mac", "Mono", "DC", "Langerhans", 
                      "Mast", "mo-DC", "Myeloid"),
  "Neural" = c("Schwann", "Glial"),
  "Other" = c("Stem cell", "TAC", "Pericyte", "Reticular cell", "Supporting cell", 
              "Unknown", "Melanocyte", "SMC-Uterus"),
  "Stromal" = c("Stromal-Ovary", "Stromal"),
  "T cells" = c("T-CD4", "T-CD8", "CTL", "T-Mucosal", "T-reg", "Th1", "T-gamma/delta", 
                "Th", "Th17")
)

# Create a new column in metadata with the general categories
merged_obj@meta.data$cell_category <- sapply(merged_obj@meta.data$cell_type, function(x) {
  cat <- names(cell_type_categories)[sapply(cell_type_categories, function(y) x %in% y)]
  if (length(cat) > 0) {
    return(cat)
  } else {
    return(NA)  # Assign NA if no category is found
  }
})

# Ensure that the cell_category column is a factor (optional)
merged_obj@meta.data$cell_category <- factor(merged_obj@meta.data$cell_category)

# Generate colors for categories and levels
category_colors <- c(
  "B cells" = "#66c2a5",
  "Endothelial" = "#fc8d62",
  "Enteroendocrine" = "#8da0cb",
  "Epithelial" = "#e78ac3",
  "Fibroblasts" = "#a6d854",
  "Lymphoid cells" = "#ffd92f",
  "Myeloid cells" = "#e5c494",
  "Neural" = "#b3b3b3",
  "Other" = "#8c6bb1",
  "Stromal" = "#6baed6",
  "T cells" = "#1f78b4"
)

# Assign colors to each level
level_colors <- unlist(lapply(names(cell_type_categories), function(cat) {
  levels <- cell_type_categories[[cat]]
  rep(category_colors[cat], length(levels))
}))

names(level_colors) <- unlist(cell_type_categories)

# Reorder factor levels in cell_type to ensure correct order in the legend
merged_obj@meta.data$cell_type <- factor(merged_obj@meta.data$cell_type, levels = unlist(cell_type_categories))

# Custom function to create headers in the legend
legend_labels <- unlist(lapply(names(cell_type_categories), function(cat) {
  c(cat, cell_type_categories[[cat]])
}))

# Create a combined vector with category colors for easier grouping
combined_colors <- c(category_colors, level_colors)
combined_colors <- combined_colors[legend_labels]

# Plot UMAP with custom legend and no axes or title
umap_plot <- DimPlot(merged_obj, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 6, repel = TRUE, 
                     raster = TRUE, raster.dpi = c(1200, 1200)) + 
  scale_color_manual(values = combined_colors, breaks = legend_labels, labels = legend_labels) +
  theme(legend.position = "right", legend.title = element_blank(),
        axis.line = element_blank(), # Remove axis lines
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank(), # Remove axis titles
        plot.title = element_blank()) + # Remove plot title
  guides(col = guide_legend(
    override.aes = list(size = 7), 
    ncol = 2, 
    byrow = FALSE,
    title = NULL,
    label.position = "right"
  ))

# Adjust the legend to group by category and create space for headers
umap_plot <- umap_plot + theme(
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  legend.box = "vertical",
  legend.spacing.y = unit(0.1, 'cm')
)

# Print the UMAP plot with the grouped custom legend
print(umap_plot)

###############################################################################
# UMAP based on donor_id variable
# donor_id
merged_obj@meta.data$donor_id <- factor(merged_obj@meta.data$donor_id)
levels(merged_obj@meta.data$donor_id)
DimPlot(object = merged_obj, reduction = "umap", group.by = "donor_id") + theme_void() +
  theme(legend.position = "right",
        plot.title = element_blank())

# Define the organ categories and their corresponding colors
organ_colors <- c(
  "uterus" = "#66c2a5",
  "intestine" = "#fc8d62",
  "skin" = "#8da0cb"
)

# Define the donor_organ list
donor_organ <- list(
  "uterus" = c("A13", "A30", "E1", "E2", "E3", "SAMN15049042", "SAMN15049043", "SAMN15049044", "SAMN15049045", "SAMN15049046", "SAMN15049047", "SAMN15049048", "SAMN15049049", "SAMN15049050", "SAMN15049051"),
  "intestine" = c("A26 (386C)", "A30 (398B)", "A38 (432C)"),
  "skin" = c("S1", "S2", "S3", "S4", "S5")
)

# Assign colors to each donor_id based on the organ category
donor_colors <- unlist(lapply(names(donor_organ), function(org) {
  donors <- donor_organ[[org]]
  setNames(rep(organ_colors[org], length(donors)), donors)
}))

# Ensure that the organ column is added to the Seurat object's metadata
merged_obj@meta.data$organ <- sapply(merged_obj@meta.data$donor_id, function(x) {
  cat <- names(donor_organ)[sapply(donor_organ, function(y) x %in% y)]
  if (length(cat) > 0) {
    return(cat)
  } else {
    return(NA)  # Assign NA if no category is found
  }
})

# Ensure that the organ column is a factor (optional)
merged_obj@meta.data$organ <- factor(merged_obj@meta.data$organ)

# Combine organ and donor labels for a grouped legend
legend_labels <- unlist(lapply(names(donor_organ), function(org) {
  c(donor_organ[[org]])
}))

# Combine colors including organ headers for proper grouping in the legend
combined_colors <- donor_colors[legend_labels]

# Plot UMAP with custom legend and no axes or title
umap_plot <- DimPlot(merged_obj, reduction = "umap", group.by = "donor_id", label = TRUE, label.size = 6, repel = TRUE, 
                     raster = TRUE, raster.dpi = c(1200, 1200)) + 
  scale_color_manual(values = combined_colors, breaks = legend_labels, labels = legend_labels) +
  theme(legend.position = "right", legend.title = element_blank(),
        axis.line = element_blank(), # Remove axis lines
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank(), # Remove axis titles
        plot.title = element_blank()) + # Remove plot title
  guides(col = guide_legend(
    override.aes = list(size = 7), 
    ncol = 1, 
    byrow = FALSE, # Set to FALSE to group by column
    title = NULL,
    label.position = "right"
  ))

# Adjust the legend to group by category and create space for headers
umap_plot <- umap_plot + theme(
  legend.text = element_text(size = 8, face = "bold"),
  legend.title = element_text(size = 10, face = "bold"),
  legend.box = "vertical",
  legend.spacing.y = unit(0.1, 'cm')
)

# Print the UMAP plot with the grouped custom legend
print(umap_plot)

###############################################################################
# Figure 2B: UMAP colored by tissue
# tissue
merged_obj@meta.data$tissue <- factor(merged_obj@meta.data$tissue)
levels(merged_obj@meta.data$tissue)
DimPlot(object = merged_obj, reduction = "umap", group.by = "tissue") + theme_void() +
  theme(legend.position = "right",
        plot.title = element_blank())

# Define the organ categories and their corresponding colors
organ_colors <- c(
  "uterus" = "#66c2a5",
  "intestine" = "#fc8d62",
  "skin" = "#8da0cb"
)

# Define the tissue_organ list
tissue_organ <- list(
  "skin" = c("skin epidermis", "dermis"),
  "intestine" = c("duodenum", "ascending colon", "caecum", "descending colon", "jejunum", "rectum", "transverse colon", 
                  "sigmoid colon", "vermiform appendix", "large intestine", "mesenteric lymph node"),
  "uterus" = c("endometrium")
)

# Assign colors to each tissue based on the organ category
tissue_colors <- unlist(lapply(names(tissue_organ), function(org) {
  tissues <- tissue_organ[[org]]
  setNames(rep(organ_colors[org], length(tissues)), tissues)
}))

# Ensure that the organ column is added to the Seurat object's metadata
merged_obj@meta.data$organ <- sapply(merged_obj@meta.data$tissue, function(x) {
  cat <- names(tissue_organ)[sapply(tissue_organ, function(y) x %in% y)]
  if (length(cat) > 0) {
    return(cat)
  } else {
    return(NA)  # Assign NA if no category is found
  }
})

# Ensure that the organ column is a factor (optional)
merged_obj@meta.data$organ <- factor(merged_obj@meta.data$organ)

# Combine tissue labels for a grouped legend
legend_labels <- unlist(lapply(names(tissue_organ), function(org) {
  c(tissue_organ[[org]])
}))

# Combine colors including organ headers for proper grouping in the legend
combined_colors <- tissue_colors[legend_labels]

# Plot UMAP with custom legend and no axes or title
umap_plot <- DimPlot(merged_obj, reduction = "umap", group.by = "tissue", label = TRUE, label.size = 6, repel = TRUE, 
                     raster = TRUE, raster.dpi = c(1200, 1200)) + 
  scale_color_manual(values = combined_colors, breaks = legend_labels, labels = legend_labels) +
  theme(legend.position = "right", legend.title = element_blank(),
        axis.line = element_blank(), # Remove axis lines
        axis.text = element_blank(), # Remove axis text
        axis.ticks = element_blank(), # Remove axis ticks
        axis.title = element_blank(), # Remove axis titles
        plot.title = element_blank()) + # Remove plot title
  guides(col = guide_legend(
    override.aes = list(size = 7), 
    ncol = 1, 
    byrow = FALSE, # Set to FALSE to group by column
    title = NULL,
    label.position = "right"
  ))

# Adjust the legend to group by category and create space for headers
umap_plot <- umap_plot + theme(
  legend.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  legend.box = "vertical",
  legend.spacing.y = unit(0.1, 'cm')
)

# Print the UMAP plot with the grouped custom legend
print(umap_plot)

###############################################################################
# Figure 2E: Heatmap for top 3 markers per cell_type
Idents(object=merged_obj) <- "cell_type"
DefaultAssay(merged_obj) <- "RNA"

# Find all markers
merged_obj.markers <- FindAllMarkers(merged_obj, test.use = "wilcox", 
                                     only.pos = TRUE, min.pct =0.2, logfc.threshold = 0.5)

write.table(merged_obj.markers,file="markers_cell_type_epi_barriers.txt",sep="\t",quote=FALSE,col.names=NA)

# Top 3 markers
markers <- read.csv("markers_cell_type_epi_barriers.csv")
head(markers)

top3.sub <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
top3.sub$gene

# Heatmaps
merged_obj.avg <- AverageExpression(merged_obj, return.seurat = TRUE)
DoHeatmap(merged_obj.avg,features=top3.sub$gene, raster = FALSE) + 
  theme(text = element_text(size = 5)) + scale_fill_gradientn(colors = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(256))

###############################################################################
# StackedBar plot 
# for donors (not shown in Figure 2; retained for reference)
# Define the colors for each category
category_colors <- c(
  "B cells" = "#66c2a5",
  "Endothelial" = "#fc8d62",
  "Enteroendocrine" = "#8da0cb",
  "Epithelial" = "#e78ac3",
  "Fibroblasts" = "#a6d854",
  "Lymphoid cells" = "#ffd92f",
  "Myeloid cells" = "#e5c494",
  "Neural" = "#b3b3b3",
  "Other" = "#8c6bb1",
  "Stromal" = "#6baed6",
  "T cells" = "#1f78b4"
)

# Define the correspondence between cell types and categories
cell_type_categories <- list(
  "B cells" = c("B", "B-mem", "B-naïve"),
  "Endothelial" = c("Endo"),
  "Enteroendocrine" = c("EC cells", "N cells", "M cells", "D cells", "L cells"),
  "Epithelial" = c("Epi-Colon", "Enterocyte", "Epi-Glandular", "Tuft", "Keratinocyte", 
                   "Epi-Ciliated", "Epi", "Goblet"),
  "Fibroblasts" = c("Fibro", "Myofibro", "Fibro-skin"),
  "Lymphoid cells" = c("ILC3", "NK", "Plasma"),
  "Myeloid cells" = c("M1 Mac", "Leukocyte", "Mac", "Mono", "DC", "Langerhans", 
                      "Mast", "mo-DC", "Myeloid"),
  "Neural" = c("Schwann", "Glial"),
  "Other" = c("Stem cell", "TAC", "Pericyte", "Reticular cell", "Supporting cell", 
              "Unknown", "Melanocyte", "SMC-Uterus"),
  "Stromal" = c("Stromal-Ovary", "Stromal"),
  "T cells" = c("T-CD4", "T-CD8", "CTL", "T-Mucosal", "T-reg", "Th1", "T-gamma/delta", 
                "Th", "Th17")
)

# Create a named vector for colors directly associated with each cell type
cell_type_colors <- unlist(lapply(names(cell_type_categories), function(cat) {
  setNames(rep(category_colors[cat], length(cell_type_categories[[cat]])), cell_type_categories[[cat]])
}))

# Map these colors in your data frame
merged_obj@meta.data$fill_color <- cell_type_colors[merged_obj@meta.data$cell_type]

# Plot using ggplot with manual scale for fill and unique legend labels
ggplot(merged_obj@meta.data, aes(x = donor_id, fill = cell_type)) + 
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(
    values = cell_type_colors, 
    name = "", 
    breaks = unique(unlist(cell_type_categories)), 
    labels = names(cell_type_colors)[match(unique(unlist(cell_type_categories)), names(cell_type_colors))]
  ) +
  labs(y = "Proportion", x = "Donors") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),   # x-axis labels
    axis.text.y = element_text(size = 14),                         # y-axis tick numbers
    axis.title.x = element_text(size = 16),                        # x-axis title "Donors"
    axis.title.y = element_text(size = 16),                        # y-axis title "Proportion"
    legend.text = element_text(size = 14),                         # legend items
    legend.title = element_text(size = 14)                         # legend title (if present)
  )

# Figure 2C: Stacked bar plot by tissue
# for tissue
ggplot(merged_obj@meta.data, aes(x = tissue, fill = cell_type)) + 
  geom_bar(stat = "count", position = "fill") +
  scale_fill_manual(
    values = cell_type_colors, 
    name = "", 
    breaks = unique(unlist(cell_type_categories)), 
    labels = names(cell_type_colors)[match(unique(unlist(cell_type_categories)), names(cell_type_colors))]
  ) +
  labs(y = "Proportion", x = "Tissues") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # x-axis labels
    axis.text.y = element_text(size = 14),                         # y-axis labels
    axis.title.x = element_text(size = 16),                        # x-axis title
    axis.title.y = element_text(size = 16),                        # y-axis title
    legend.text = element_text(size = 14),                         # legend items
    legend.title = element_text(size = 14)                         # legend title (if not blank)
  )

###############################################################################
# Violin plot

# Load the dplyr library for data manipulation
library(dplyr)

# Create a mapping vector for new general categories
new_categories <- c(
  "B" = "B cells", "B-mem" = "B cells", "B-naïve" = "B cells",
  "Endo" = "Endothelial",
  "EC cells" = "Enteroendocrine", "N cells" = "Enteroendocrine", "M cells" = "Enteroendocrine", 
  "D cells" = "Enteroendocrine", "L cells" = "Enteroendocrine",
  "Epi-Colon" = "Epithelial", "Enterocyte" = "Epithelial", "Epi-Glandular" = "Epithelial", 
  "Tuft" = "Epithelial", "Keratinocyte" = "Epithelial", "Epi-Ciliated" = "Epithelial", 
  "Epi" = "Epithelial", "Goblet" = "Epithelial",
  "Fibro" = "Fibroblasts", "Myofibro" = "Fibroblasts", "Fibro-skin" = "Fibroblasts",
  "ILC3" = "Lymphoid cells", "NK" = "Lymphoid cells", "Plasma" = "Lymphoid cells",
  "M1 Mac" = "Myeloid cells", "Leukocyte" = "Myeloid cells", "Mac" = "Myeloid cells", 
  "Mono" = "Myeloid cells", "DC" = "Myeloid cells", "Langerhans" = "Myeloid cells", 
  "Mast" = "Myeloid cells", "mo-DC" = "Myeloid cells", "Myeloid" = "Myeloid cells",
  "Schwann" = "Neural", "Glial" = "Neural",
  "Stem cell" = "Other", "TAC" = "Other", "Pericyte" = "Other", "Reticular cell" = "Other", 
  "Supporting cell" = "Other", "Unknown" = "Other", "Melanocyte" = "Other", "SMC-Uterus" = "Other",
  "Stromal-Ovary" = "Stromal", "Stromal" = "Stromal",
  "T-CD4" = "T cells", "T-CD8" = "T cells", "CTL" = "T cells", "T-Mucosal" = "T cells", 
  "T-reg" = "T cells", "Th1" = "T cells", "T-gamma/delta" = "T cells", "Th" = "T cells", "Th17" = "T cells"
)

# Extract the metadata from your Seurat object
metadata <- merged_obj@meta.data

# Map the new categories safely
metadata$cell_type_general <- sapply(metadata$cell_type, function(x) {
  if (x %in% names(new_categories)) {
    new_categories[x]
  } else {
    NA  # This will help identify if there are any cell types not mapped
  }
})

# Check for any NA values which indicate unmapped cell types
sum(is.na(metadata$cell_type_general))

# Add the updated metadata back to the Seurat object
merged_obj <- AddMetaData(merged_obj, metadata = metadata)

###############################################################################
# NOTE: The violin plot and feature plot sections below are supplementary
# analyses not shown in Figure 2. They are retained here for reference.
# The 'plots' object below should be generated via VlnPlot() before use.
# Example:
# plots <- VlnPlot(merged_obj, features = c("EPCAM", "VIM", "PTPRC"),
#                  group.by = "cell_type_general", pt.size = 0, combine = FALSE)

plots <- lapply(
  X = plots,
  FUN = function(p) p + 
    ggplot2::scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3", "#8c6bb1", "#6baed6", "#1f78b4")) +
    ggplot2::theme(
      legend.text = element_text(size = 16),
      legend.spacing.y = unit(0.6, 'cm')  # <-- increase to add more space between items
    )
)


library(patchwork)

# Combine with wrap_plots and set legend position
combined_plot <- wrap_plots(plots, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

combined_plot

###############################################################################
# Feature plots
# Define genes of interest
genes <- c("EPCAM", "VIM", "PTPRC")

# Generate feature plots
feature_plots <- FeaturePlot(
  object = merged_obj,
  features = genes,
  reduction = "umap",
  pt.size = 0.3,
  raster = TRUE, raster.dpi = c(500, 500),
  order = TRUE,
  max.cutoff = "q95"  # Clip at 95th percentile
)

# Show the plots
feature_plots
