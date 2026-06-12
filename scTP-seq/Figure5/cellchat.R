# Clear environment
rm(list = ls())

############## Merge multiple subclusters back to major clusters ##############

library(Seurat)
library(ggplot2)
library(grid)

# Color palette
col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1", "#FF9999", "#66CCCC", "#4F6272", "#FF3366",
         "#CC0066", "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
         "#E5D2DD", "#58A4C3", "#F9BB72", "#F3B1A0", "#57C3F3",
         "#E59CC4", "#437eb8", "#66CCCC", "#99CCFF", "#3399CC",
         "#FF3366", "#CC0066", "#FF9933", "#CCFFCC", "#00CC66",
         "#99FFFF", "#FF3300", "#6699CC", "#9999FF", "#CCCCFF",
         "#FF6699", "#6699CC", "#FFFFCC")

# Output directory
outdir <- "./out"
merge_output <- file.path(outdir, "merge_subclusters")
dir.create(merge_output, recursive = TRUE, showWarnings = FALSE)

# Input files
# If your folders still use Chinese names, change these two paths back to:
# file.path(outdir, "DCs细胞", "celltype.rds")
# file.path(outdir, "T细胞", "celltype.rds")
dc_subcluster_rds <- file.path(outdir, "DCs_cells", "celltype.rds")
tcell_subcluster_rds <- file.path(outdir, "T_cells", "celltype.rds")
major_cluster_rds <- file.path(outdir, "celltype.rds")

# Read data
DCs <- readRDS(dc_subcluster_rds)          # DC subcluster annotation
Tcells <- readRDS(tcell_subcluster_rds)   # T cell subcluster annotation
ScRNA <- readRDS(major_cluster_rds)       # Major cell type annotation for all cells

# Set identities
Idents(DCs) <- "celltype"
Idents(Tcells) <- "celltype"
Idents(ScRNA) <- "celltype"

# Extract selected major cell types
selected_celltypes <- c("T cells", "DCs")
all <- subset(ScRNA, idents = selected_celltypes)

table(all@meta.data$celltype)

# Merge subcluster annotations into the major-cluster object
Idents(all, cells = colnames(DCs)) <- Idents(DCs)
Idents(all, cells = colnames(Tcells)) <- Idents(Tcells)

# Create a new metadata column for merged cell type annotation
all$celltype_merged <- Idents(all)

# Update identities
Idents(all) <- "celltype_merged"

# Check merged cell type statistics
table(all@meta.data$celltype_merged)

# Save merged Seurat object
merged_rds <- file.path(outdir, "celltype_merged_subclusters.rds")
saveRDS(all, file = merged_rds)

# UMAP plot by merged cell type: PDF
pdf(file.path(merge_output, "merged_celltype_umap.pdf"), width = 10, height = 6)
DimPlot(
  object = all,
  group.by = "celltype_merged",
  reduction = "umap",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_blank()
  )
dev.off()

# UMAP plot by merged cell type: SVG
svg(file.path(merge_output, "merged_celltype_umap.svg"), width = 10, height = 6)
DimPlot(
  object = all,
  group.by = "celltype_merged",
  reduction = "umap",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_blank()
  )
dev.off()

# UMAP plot by treatment: PDF
pdf(file.path(merge_output, "treatment_umap.pdf"), width = 6, height = 4)
DimPlot(
  object = all,
  group.by = "treatment",
  reduction = "umap",
  pt.size = 0.1,
  label = FALSE,
  label.size = 6,
  repel = TRUE,
  cols = col
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.title = element_blank()
  )
dev.off()

# UMAP plot by treatment: SVG
svg(file.path(merge_output, "treatment_umap.svg"), width = 6, height = 4)
DimPlot(
  object = all,
  group.by = "treatment",
  reduction = "umap",
  pt.size = 0.1,
  label = FALSE,
  label.size = 6,
  repel = TRUE,
  cols = col
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.title = element_blank()
  )
dev.off()

# UMAP split by treatment: PDF
pdf(file.path(merge_output, "split_treatment_umap.pdf"), width = 13, height = 5)
DimPlot(
  all,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()

# UMAP split by treatment: SVG
svg(file.path(merge_output, "split_treatment_umap.svg"), width = 13, height = 5)
DimPlot(
  all,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
) +
  theme(
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()



############################# Cell-cell communication comparison #############################
############################# Normal vs Tumor groups #############################

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)
library(ComplexHeatmap)
library(cowplot)
library(gridExtra)
library(grid)
library(circlize)

options(stringsAsFactors = FALSE)

# Define colors
col <- c("#A4CDE1", "#FF9999", "#66CCCC", "#FFCCCC", "#CCFFCC",
         "#E5D2DD", "#4F6272", "#CC99CC", "#F9BB72", "#57C3F3",
         "#E59CC4", "#437eb8", "#99CCFF", "#3399CC", "#FF3366",
         "#FF9933", "#9999FF", "#00CC66", "#99FFFF", "#FF3300",
         "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC")

# CellChat output directory
cellchat_output <- file.path(outdir, "cellchat")
dir.create(cellchat_output, recursive = TRUE, showWarnings = FALSE)

# Read merged Seurat object
seuratdata <- readRDS(merged_rds)

head(seuratdata)
View(seuratdata@meta.data)

# Extract Normal and Tumor samples
normal_samples <- subset(seuratdata, subset = treatment %in% c("Normal"))
tumor_samples <- subset(seuratdata, subset = treatment %in% c("Tumor"))

# Create CellChat objects
cellchat_normal <- createCellChat(
  object = normal_samples@assays$RNA@data,
  meta = normal_samples@meta.data,
  group.by = "celltype_merged"
)

cellchat_tumor <- createCellChat(
  object = tumor_samples@assays$RNA@data,
  meta = tumor_samples@meta.data,
  group.by = "celltype_merged"
)

# Remove unused factor levels
cellchat_normal@idents <- droplevels(cellchat_normal@idents)
cellchat_tumor@idents <- droplevels(cellchat_tumor@idents)

# Analyze Normal group
cellchat_normal@DB <- CellChatDB.human
cellchat_normal <- subsetData(cellchat_normal)
cellchat_normal <- identifyOverExpressedGenes(cellchat_normal)
cellchat_normal <- identifyOverExpressedInteractions(cellchat_normal)
cellchat_normal <- computeCommunProb(cellchat_normal, raw.use = TRUE, population.size = TRUE)
cellchat_normal <- computeCommunProbPathway(cellchat_normal)
cellchat_normal <- aggregateNet(cellchat_normal)
cellchat_normal <- netAnalysis_computeCentrality(cellchat_normal, slot.name = "netP")
saveRDS(cellchat_normal, file.path(cellchat_output, "cellchat_normal.rds"))

# Analyze Tumor group
cellchat_tumor@DB <- CellChatDB.human
cellchat_tumor <- subsetData(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedGenes(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedInteractions(cellchat_tumor)
cellchat_tumor <- computeCommunProb(cellchat_tumor, raw.use = TRUE, population.size = TRUE)
cellchat_tumor <- computeCommunProbPathway(cellchat_tumor)
cellchat_tumor <- aggregateNet(cellchat_tumor)
cellchat_tumor <- netAnalysis_computeCentrality(cellchat_tumor, slot.name = "netP")
saveRDS(cellchat_tumor, file.path(cellchat_output, "cellchat_tumor.rds"))

# Merge CellChat objects
cellchat_list <- list(normal = cellchat_normal, tumor = cellchat_tumor)
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list), cell.prefix = TRUE)
saveRDS(cellchat, file.path(cellchat_output, "cellchat_merged.rds"))

# Reload CellChat objects
cellchat_tumor <- readRDS(file.path(cellchat_output, "cellchat_tumor.rds"))
cellchat_normal <- readRDS(file.path(cellchat_output, "cellchat_normal.rds"))
cellchat <- readRDS(file.path(cellchat_output, "cellchat_merged.rds"))

# Overview of interaction number and strength
gg1 <- compareInteractions(
  cellchat,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "count",
  color.use = col
) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.caption = element_text(size = 14)
  )

gg2 <- compareInteractions(
  cellchat,
  show.legend = FALSE,
  group = c(1, 2),
  measure = "weight",
  color.use = col
) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.caption = element_text(size = 14)
  )

p <- gg1 + gg2
ggsave(file.path(cellchat_output, "overview_interaction_number_strength.pdf"), p, width = 6, height = 4)



################### Ligand-receptor comparison analysis ###################

# Check cell type identity levels
levels(cellchat@idents$joint)
levels(cellchat_tumor@idents)
levels(cellchat_normal@idents)

##### DC-related signaling in Tumor group #####
p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(17, 10, 12, 13, 15),
  targets.use = c(2, 3, 4, 5, 6, 7, 8),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Tumor_DCs.pdf"), p, width = 16, height = 20)


##### DC-related signaling in Normal group #####
p <- netVisual_bubble(
  cellchat_normal,
  sources.use = c(8, 10, 12, 13, 14),
  targets.use = c(1, 2, 3, 4, 5, 7),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Normal_DCs.pdf"), p, width = 16, height = 20)


# Select ligand-receptor pairs of interest
pairLR.use <- as.data.frame(c(
  "BTLA_TNFRSF14", "CD80_CTLA4", "CD86_CTLA4",
  "ICOSL_CTLA4", "ICOSL_CD28", "PD-L1_PD-1", "HLA-E_CD8A",
  "HLA-E_KLRC1", "MIF_CD74", "MIF_CD74_CXCR4",
  "MIF_(CD74+CXCR4)", "MIF_(CD74+CD44)", "TIGIT_NECTIN2",
  "LGALS9_HAVCR2", "SEMA4A_PLXNB1", "SEMA4D_PLXNB1",
  "LGALS9_HAVCR2", "IL1A-IL1R2", "IL1B_IL1R2", "TGFA_EGFR",
  "CLEC2B_KLRB1", "CLEC2C_KLRB1", "IL1B_IL1R2",
  "NECTIN2_TIGIT", "SPP1_ITGA4_ITGB1", "SPP1_CD44"
))

colnames(pairLR.use) <- "interaction_name"

# Bubble plot for selected ligand-receptor pairs
p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(17, 10, 12, 13, 15),
  targets.use = c(2, 3, 4, 5, 6, 7, 8),
  pairLR.use = pairLR.use,
  angle.x = 45
) +
  theme(
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )

ggsave(file.path(cellchat_output, "LR_bubble_selected_pairs_Tumor_DCs.pdf"), p, width = 16, height = 6)


##### T cell-related signaling in Tumor group #####
levels(cellchat_tumor@idents)

p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(3, 5, 6, 8),
  targets.use = c(1, 2, 4, 7, 9, 10),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Tumor_T_cells.pdf"), p, width = 10, height = 10)


##### Exhausted CD8 T cell-related signaling in Tumor group #####
levels(cellchat_tumor@idents)

p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(5),
  targets.use = c(1, 2, 3, 4, 6, 7, 8, 9, 10),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Tumor_Exhausted_CD8T.pdf"), p, width = 8, height = 10)


##### Exhausted CD4 T cell-related signaling in Tumor group #####
levels(cellchat_tumor@idents)

p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(3),
  targets.use = c(1, 2, 4, 5, 6, 7, 8, 9, 10),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Tumor_Exhausted_CD4T.pdf"), p, width = 10, height = 10)


##### Treg 1-related signaling in Tumor group #####
levels(cellchat_tumor@idents)

p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(8),
  targets.use = c(1, 2, 3, 4, 5, 6, 7, 9, 10),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Tumor_Treg1.pdf"), p, width = 10, height = 10)


##### Treg 2-related signaling in Tumor group #####
levels(cellchat_tumor@idents)

p <- netVisual_bubble(
  cellchat_tumor,
  sources.use = c(6),
  targets.use = c(1, 2, 3, 4, 5, 7, 8, 9, 10),
  angle.x = 45
)

p <- p + theme(
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 12)
)

ggsave(file.path(cellchat_output, "LR_bubble_Tumor_Treg2.pdf"), p, width = 10, height = 10)


# Show upregulated and downregulated ligand-receptor pairs
p1 <- netVisual_bubble(
  cellchat,
  sources.use = c(8),
  targets.use = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11),
  comparison = c(1, 2),
  max.dataset = 2,
  title.name = "Increased signaling",
  angle.x = 45,
  remove.isolate = TRUE
)

p2 <- netVisual_bubble(
  cellchat,
  sources.use = c(8),
  targets.use = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11),
  comparison = c(1, 2),
  max.dataset = 1,
  title.name = "Decreased signaling",
  angle.x = 45,
  remove.isolate = TRUE
)

pc <- p1 + p2
ggsave(file.path(cellchat_output, "LR_regulated_DCs_T_cells.pdf"), pc, width = 12, height = 8)



################### Network and heatmap analysis ###################

# Differential interaction network: interaction number
png(file.path(cellchat_output, "diff_interaction_count.png"), width = 800, height = 800)
par(cex = 2)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
dev.off()

# Differential interaction network: interaction strength
png(file.path(cellchat_output, "diff_interaction_weight.png"), width = 800, height = 800)
par(cex = 2)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

p1_img <- ggdraw() + draw_image(file.path(cellchat_output, "diff_interaction_count.png"))
p2_img <- ggdraw() + draw_image(file.path(cellchat_output, "diff_interaction_weight.png"))
combined_plot <- plot_grid(p1_img, p2_img, ncol = 2)

ggsave(
  file.path(cellchat_output, "diff_combined_interaction.pdf"),
  combined_plot,
  width = 8,
  height = 4
)


# Differential heatmap: interaction number
png(file.path(cellchat_output, "diff_heatmap_count.png"), width = 900, height = 800, res = 150)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex = 1.5)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

# Differential heatmap: interaction strength
png(file.path(cellchat_output, "diff_heatmap_weight.png"), width = 900, height = 800, res = 150)
par(cex.axis = 1.5, cex.lab = 1.5, cex.main = 2, cex = 1.5)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

p1_img <- ggdraw() + draw_image(file.path(cellchat_output, "diff_heatmap_count.png"))
p2_img <- ggdraw() + draw_image(file.path(cellchat_output, "diff_heatmap_weight.png"))
combined_plot <- plot_grid(p1_img, p2_img, ncol = 2)

ggsave(
  file.path(cellchat_output, "diff_combined_heatmap.pdf"),
  combined_plot,
  width = 12,
  height = 5
)


# Circle plots comparing interaction number between groups
weight.max <- getMaxWeight(cellchat_list, attribute = c("org.idents", "count"))

for (i in seq_along(cellchat_list)) {
  pdf_filename <- file.path(
    cellchat_output,
    paste0("circle_compare_interaction_number_", names(cellchat_list)[i], ".pdf")
  )
  
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Number of interactions - ", names(cellchat_list)[i])
  
  par(oma = c(0, 0, 0, 0))
  par(cex = 2)
  
  netVisual_circle(
    cellchat_list[[i]]@net$count,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max[2],
    edge.width.max = 10
  )
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off()
}


# Circle plots comparing interaction strength between groups
weight.max.weight <- getMaxWeight(cellchat_list, attribute = c("org.idents", "weight"))

for (i in seq_along(cellchat_list)) {
  pdf_filename <- file.path(
    cellchat_output,
    paste0("circle_compare_interaction_strength_", names(cellchat_list)[i], ".pdf")
  )
  
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Interaction strength - ", names(cellchat_list)[i])
  
  par(oma = c(0, 0, 0, 0))
  par(cex = 2)
  
  netVisual_circle(
    cellchat_list[[i]]@net$weight,
    weight.scale = TRUE,
    label.edge = FALSE,
    edge.weight.max = weight.max.weight[2],
    edge.width.max = 30
  )
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off()
}



################### Chord diagram analysis ###################

# Chord diagrams comparing interaction number between groups
weight.max <- getMaxWeight(cellchat_list, attribute = c("org.idents", "count"))

for (i in seq_along(cellchat_list)) {
  pdf_filename <- file.path(
    cellchat_output,
    paste0("chord_compare_interaction_number_", names(cellchat_list)[i], ".pdf")
  )
  
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Number of interactions - ", names(cellchat_list)[i])
  
  interaction_data <- cellchat_list[[i]]@net$count
  
  chordDiagram(
    interaction_data,
    preAllocate = 1,
    direction.type = c("diffHeight"),
    link.arr.type = "big.arrow",
    annotationTrack = c("name", "grid"),
    link.border = NA,
    transparency = 0.5
  )
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off()
}


# Chord diagrams comparing interaction strength between groups
weight.max.weight <- getMaxWeight(cellchat_list, attribute = c("org.idents", "weight"))

for (i in seq_along(cellchat_list)) {
  pdf_filename <- file.path(
    cellchat_output,
    paste0("chord_compare_interaction_strength_", names(cellchat_list)[i], ".pdf")
  )
  
  pdf(pdf_filename, width = 10, height = 10)
  plot_title <- paste0("Interaction strength - ", names(cellchat_list)[i])
  
  interaction_weight_data <- cellchat_list[[i]]@net$weight
  
  chordDiagram(
    interaction_weight_data,
    preAllocate = 1,
    direction.type = c("diffHeight"),
    link.arr.type = "big.arrow",
    annotationTrack = c("name", "grid"),
    link.border = NA,
    transparency = 0.5
  )
  
  title(main = plot_title, line = 1, cex.main = 1.5)
  dev.off()
}



################### Identification and visualization of conserved and specific signaling pathways ###################

gg1 <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)

p <- gg1 + gg2

ggsave(
  file.path(cellchat_output, "compare_pathway_strength.pdf"),
  p,
  width = 10,
  height = 10
)

saveRDS(cellchat, file.path(cellchat_output, "cellchat_compare.rds"))



################### Comparison of cell signaling patterns ###################

pathway.union <- Reduce(
  union,
  list(
    cellchat_list[[1]]@netP$pathways,
    cellchat_list[[2]]@netP$pathways
  )
)

# All signaling roles
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat_list[[1]],
  pattern = "all",
  signaling = pathway.union,
  title = names(cellchat_list)[1],
  width = 6,
  height = 22
)

ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat_list[[2]],
  pattern = "all",
  signaling = pathway.union,
  title = names(cellchat_list)[2],
  width = 6,
  height = 22
)

ht1_grob <- grid.grabExpr(draw(ht1))
ht2_grob <- grid.grabExpr(draw(ht2))

combined_plot <- grid.arrange(
  ht1_grob,
  ht2_grob,
  ncol = 2,
  widths = unit.c(unit(12, "cm"), unit(11, "cm"))
)

pdf(file.path(cellchat_output, "combined_signaling_role_heatmap_all.pdf"), width = 10, height = 12)
grid.draw(combined_plot)
dev.off()


# Outgoing signaling roles
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat_list[[1]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(cellchat_list)[1],
  width = 6,
  height = 22
)

ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat_list[[2]],
  pattern = "outgoing",
  signaling = pathway.union,
  title = names(cellchat_list)[2],
  width = 6,
  height = 22
)

ht1_grob <- grid.grabExpr(draw(ht1))
ht2_grob <- grid.grabExpr(draw(ht2))

combined_plot <- grid.arrange(
  ht1_grob,
  ht2_grob,
  ncol = 2,
  widths = unit.c(unit(12, "cm"), unit(11, "cm"))
)

pdf(file.path(cellchat_output, "combined_signaling_role_heatmap_outgoing.pdf"), width = 10, height = 12)
grid.draw(combined_plot)
dev.off()


# Incoming signaling roles
ht1 <- netAnalysis_signalingRole_heatmap(
  cellchat_list[[1]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(cellchat_list)[1],
  width = 6,
  height = 22
)

ht2 <- netAnalysis_signalingRole_heatmap(
  cellchat_list[[2]],
  pattern = "incoming",
  signaling = pathway.union,
  title = names(cellchat_list)[2],
  width = 6,
  height = 22
)

ht1_grob <- grid.grabExpr(draw(ht1))
ht2_grob <- grid.grabExpr(draw(ht2))

combined_plot <- grid.arrange(
  ht1_grob,
  ht2_grob,
  ncol = 2,
  widths = unit.c(unit(12, "cm"), unit(11, "cm"))
)

pdf(file.path(cellchat_output, "combined_signaling_role_heatmap_incoming.pdf"), width = 10, height = 12)
grid.draw(combined_plot)
dev.off()



################### Comparison of selected signaling pathways ###################

# Check available signaling pathways in Tumor group
df.net <- subsetCommunication(cellchat_tumor)
table(df.net$pathway_name)

# Check available signaling pathways in Normal group
df.net <- subsetCommunication(cellchat_normal)
table(df.net$pathway_name)

# Signaling pathways of interest
pathways_list <- c(
  "BTLA", "CD226", "CD137", "PD-L1", "TIGIT", "TNF", "IFN-II", "TGFb",
  "CD137", "CD80", "CD86", "TIGIT", "CDH1", "MIF", "NOTCH"
)

# Generate plots for each selected signaling pathway
for (pathways.show in pathways_list) {
  
  # Check whether the pathway is available in any CellChat object
  available_pathways <- unique(unlist(lapply(cellchat_list, function(x) x@netP$pathways)))
  
  if (!(pathways.show %in% available_pathways)) {
    message(paste("Skipping pathway:", pathways.show, "- not found in any CellChat object"))
    next
  }
  
  # Calculate maximum pathway weight
  weight.max <- tryCatch(
    {
      getMaxWeight(cellchat_list, slot.name = "netP", attribute = pathways.show)
    },
    error = function(e) {
      message(paste("Error in getMaxWeight for pathway:", pathways.show))
      return(NULL)
    }
  )
  
  if (is.null(weight.max)) {
    next
  }
  
  # Generate gene expression heatmaps, contribution plots, and interaction heatmaps
  for (i in seq_along(cellchat_list)) {
    
    sample_output <- paste0("compare_", pathways.show, "_", names(cellchat_list)[i])
    
    # Gene expression heatmap
    png(
      file = file.path(cellchat_output, paste0(sample_output, "_gene_expression_heatmap.png")),
      width = 1200,
      height = 1400,
      res = 300
    )
    p3 <- plotGeneExpression(cellchat_list[[i]], signaling = pathways.show)
    print(p3)
    dev.off()
    
    # Contribution analysis
    png(
      file = file.path(cellchat_output, paste0(sample_output, "_contribution.png")),
      width = 1200,
      height = 800,
      res = 300
    )
    p4 <- netAnalysis_contribution(cellchat_list[[i]], signaling = pathways.show)
    print(p4)
    dev.off()
    
    # Interaction heatmap
    png(
      file = file.path(cellchat_output, paste0(sample_output, "_netVisual_heatmap.png")),
      width = 1200,
      height = 1200,
      res = 300
    )
    p5 <- netVisual_heatmap(
      cellchat_list[[i]],
      signaling = pathways.show,
      color.heatmap = "Reds"
    )
    print(p5)
    dev.off()
  }
  
  # Combine gene expression heatmaps
  gene_expr_filenames <- lapply(
    seq_along(cellchat_list),
    function(i) file.path(
      cellchat_output,
      paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_gene_expression_heatmap.png")
    )
  )
  
  gene_expr_imgs <- lapply(gene_expr_filenames, function(x) ggdraw() + draw_image(x))
  combined_gene_expr_plot <- plot_grid(plotlist = gene_expr_imgs, ncol = 2)
  
  ggsave(
    file.path(cellchat_output, paste0("combined_compare_", pathways.show, "_gene_expression_heatmap.pdf")),
    combined_gene_expr_plot,
    width = 10,
    height = 8
  )
  
  # Combine contribution plots
  contribution_filenames <- lapply(
    seq_along(cellchat_list),
    function(i) file.path(
      cellchat_output,
      paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_contribution.png")
    )
  )
  
  contribution_imgs <- lapply(contribution_filenames, function(x) ggdraw() + draw_image(x))
  combined_contribution_plot <- plot_grid(plotlist = contribution_imgs, ncol = 2)
  
  ggsave(
    file.path(cellchat_output, paste0("combined_compare_", pathways.show, "_contribution.pdf")),
    combined_contribution_plot,
    width = 8,
    height = 8
  )
  
  # Combine interaction heatmaps
  heatmap_filenames <- lapply(
    seq_along(cellchat_list),
    function(i) file.path(
      cellchat_output,
      paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_netVisual_heatmap.png")
    )
  )
  
  heatmap_imgs <- lapply(heatmap_filenames, function(x) ggdraw() + draw_image(x))
  combined_heatmap_plot <- plot_grid(plotlist = heatmap_imgs, ncol = 2)
  
  ggsave(
    file.path(cellchat_output, paste0("combined_compare_", pathways.show, "_netVisual_heatmap.pdf")),
    combined_heatmap_plot,
    width = 8,
    height = 8
  )
  
  # Circle network plots for selected signaling pathway
  for (i in seq_along(cellchat_list)) {
    
    pdf_filename <- file.path(
      cellchat_output,
      paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_network.pdf")
    )
    
    pdf(pdf_filename, width = 10, height = 10)
    
    plot_title <- paste0(pathways.show, " - ", names(cellchat_list)[i])
    
    par(oma = c(0, 0, 0, 0))
    par(cex = 2)
    
    netVisual_aggregate(
      cellchat_list[[i]],
      signaling = pathways.show,
      layout = "circle",
      edge.weight.max = weight.max[1],
      edge.width.max = 30,
      vertex.label.cex = 1.2
    )
    
    title(main = plot_title, line = 1, cex.main = 1.5)
    dev.off()
  }
  
  # Chord plots for selected signaling pathway
  for (i in seq_along(cellchat_list)) {
    
    png_filename <- file.path(
      cellchat_output,
      paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_chord.png")
    )
    
    png(png_filename, width = 800, height = 800)
    
    par(oma = c(0, 0, 0, 0))
    par(cex = 2)
    
    netVisual_aggregate(
      cellchat_list[[i]],
      signaling = pathways.show,
      layout = "chord",
      pt.title = 3,
      title.space = 0.05,
      signaling.name = paste(pathways.show, names(cellchat_list)[i]),
      vertex.label.cex = 1,
      font.main = 2
    )
    
    dev.off()
  }
  
  # Combine chord plots
  chord_filenames <- lapply(
    seq_along(cellchat_list),
    function(i) file.path(
      cellchat_output,
      paste0("compare_", pathways.show, "_", names(cellchat_list)[i], "_chord.png")
    )
  )
  
  chord_imgs <- lapply(chord_filenames, function(x) ggdraw() + draw_image(x))
  combined_chord_plot <- plot_grid(plotlist = chord_imgs, ncol = 2)
  
  ggsave(
    file.path(cellchat_output, paste0("combined_compare_", pathways.show, "_chord.pdf")),
    combined_chord_plot,
    width = 8,
    height = 8
  )
}


