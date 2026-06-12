# Clear the R environment
rm(list = ls())


######### Plot Cell Type Heatmap ############

# Load required packages
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Create output directory
output <- paste(outdir, "celltype", sep = "/")
dir.create(output, showWarnings = FALSE)

# Read the RDS file, assumed to be a Seurat object
file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

# Set the default assay to RNA if it has not been specified
DefaultAssay(scedata) <- "RNA"

# Define genes of interest
cellmarker <- c(
  "ELL2", "MALAT1", "BRAF", "PDE4B", "TNFAIP3", "LINC-PINT", "NFKBIA",
  "ARID4B", "ATP1B3", "HSPD1", "HSP90AB1", "HSP90AA1", "HSPA1B", "HSPA1A"
)

# Filter out genes that are not present in the dataset
cellmarker_present <- intersect(cellmarker, rownames(scedata))

# Calculate average expression of selected genes in each cell type
avg_expr <- AverageExpression(
  scedata,
  features = cellmarker_present,
  group.by = "celltype"
)$RNA

# Convert the expression data into matrix format
expr_mat <- as.matrix(avg_expr)

# Standardize the entire expression matrix using Z-score
# mean_all <- mean(expr_mat, na.rm = TRUE)
# sd_all <- sd(expr_mat, na.rm = TRUE)
# expr_mat_scaled <- (expr_mat - mean_all) / sd_all

# Perform row-wise Z-score normalization for each gene
expr_mat_scaled <- t(scale(t(expr_mat)))

# Extract column names as cell type names
celltypes <- colnames(expr_mat_scaled)

# Assign colors to cell types
library(RColorBrewer)

celltype_colors <- setNames(
  colorRampPalette(brewer.pal(8, "Set2"))(length(celltypes)),
  celltypes
)

# Create column annotation
col_anno <- HeatmapAnnotation(
  Celltype = celltypes,
  col = list(Celltype = celltype_colors),
  annotation_name_gp = gpar(fontsize = 10)
)

# Set output file path
heatmap_output <- file.path(output, "cellmarker_expression_heatmap.pdf")

# Plot heatmap and save as PDF
pdf(heatmap_output, width = 8, height = 6)

Heatmap(
  expr_mat_scaled,
  name = "Z-score",
  top_annotation = col_anno,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 14),
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_columns = TRUE,
  cluster_rows = FALSE,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 14),
    labels_gp = gpar(fontsize = 13)
  ),
  col = colorRamp2(
    c(-2, 0, 2),
    c("#3399CC", "white", "#FF3366")
  )
)

dev.off()


