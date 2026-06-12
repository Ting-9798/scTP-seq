# Clear the R environment
rm(list = ls())

# Load required packages
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)

####################### Seurat Analysis #####################

# Define color palette
col <- c("#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
         '#E5D2DD', '#4F6272', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699",
         "#6699CC", "#FFFFCC")

setwd("./out(rna)/AT2/")
outdir <- "./out(rna)/AT2/"

MergeOUT <- paste(outdir, "Merge", sep = "/")
dir.create(MergeOUT, recursive = TRUE, showWarnings = FALSE)
OUTPUT <- paste(MergeOUT, "Multiple_", sep = "/")

data <- readRDS("./out(rna)/celltype.rds")

# Select AT2 cells
Cells.sub <- subset(data@meta.data, celltype == "AT2")
summary(Cells.sub$celltype)

ScRNA <- subset(data, cells = row.names(Cells.sub))
View(ScRNA@meta.data)


#### 6. Scaling and PCA Dimensionality Reduction ####

# Scale data
ScRNA <- ScaleData(ScRNA)

# Run PCA
ScRNA <- RunPCA(ScRNA, npcs = 30)

pdf(paste(OUTPUT, "Dimplot.pdf"), width = 9, height = 6)
p1 <- DimPlot(
  object = ScRNA,
  reduction = "pca",
  pt.size = .1,
  group.by = "treatment",
  cols = col
)
CombinePlots(plots = list(p1))
dev.off()

pdf(paste(OUTPUT, "vlnplot.pdf"), width = 9, height = 6)
p2 <- VlnPlot(
  object = ScRNA,
  features = "PC_1",
  group.by = "treatment",
  pt.size = 0,
  cols = col
)
CombinePlots(plots = list(p2))
dev.off()

# PCA heatmap visualization
pdf(paste(OUTPUT, "DimHeatmap.pdf"), width = 9, height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Evaluate the number of principal components
pdf(paste0(OUTPUT, "PCA-ElbowPlot.pdf"), width = 6, height = 5)
ElbowPlot(ScRNA)
dev.off()

save(ScRNA, file = "ScRNA_before_clustering.RData")


#### 7. Cell Clustering and Annotation ####

col <- c('#E5D2DD', '#FF6666', '#6A4C93', "#BC8F8F", '#FFCC99', '#FF9999',
         "#FFCCCC", '#4F6272', '#58A4C3',
         "#66CCCC", '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#CCFFCC", "#00CC66", "#99FFFF",
         "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#FF3300", "#6699CC", "#9999FF", "#CCCCFF",
         "#CC99CC", "#FF6699", "#6699CC", "#FFFFCC")

load("ScRNA_before_clustering.RData")

# Cell clustering
ScRNA <- ScRNA %>%
  RunUMAP(dims = 1:20) %>%
  RunTSNE(dims = 1:20) %>%
  FindNeighbors(dims = 1:20)

ScRNA <- FindClusters(
  ScRNA,
  resolution = seq(from = 0.1, to = 1.0, by = 0.1)
)

# pdf(paste(OUTPUT, "clustree.pdf"), width = 10, height = 9)
# library(clustree)
# clustree(ScRNA)
# dev.off()

# Idents(ScRNA) <- "integrated_snn_res.1"
Idents(ScRNA) <- "RNA_snn_res.1"

# Select resolution based on the clustering tree
ScRNA$seurat_clusters <- ScRNA@active.ident
table(Idents(ScRNA))

# Set treatment factor order if needed
# ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Visualize clusters split by treatment
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 8, height = 4)
DimPlot(
  ScRNA,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "treatment",
  label.size = 5,
  cols = col
) +
  theme(
    strip.text = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )
dev.off()

# Generate UMAP plot colored by cluster
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 5.5, height = 5)
DimPlot(
  ScRNA,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 6,
  cols = col
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22)
  )
dev.off()

pdf(paste(OUTPUT, "cluster_umap1.pdf"), width = 6, height = 4)
DimPlot(
  ScRNA,
  reduction = "umap",
  label = FALSE,
  repel = TRUE,
  cols = col,
  group.by = "treatment"
) +
  ggtitle(NULL) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22)
  )
dev.off()

# Visualize treatment groups together
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"), width = 6, height = 6)
DimPlot(
  ScRNA,
  repel = TRUE,
  reduction = "umap",
  group.by = "treatment"
) +
  scale_color_manual(values = col) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"),
    legend.position = c(.01, .1)
  ) +
  labs(title = "Sample Origin")
dev.off()

saveRDS(ScRNA, "ScRNA_clustered.rds")


# Visualize marker gene expression patterns across clusters
output <- paste(outdir, "cell_localization", sep = "/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

file_path <- file.path(outdir, "ScRNA_clustered.rds")
ScRNA <- readRDS(file_path)

# Identify marker genes
ScRNA.markers <- FindAllMarkers(
  ScRNA,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(ScRNA.markers, paste0(output, "./ScRNA.all.markers.csv"))

dim.use <- 1:30

top5 <- ScRNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

write.csv(
  top5,
  file = paste0(output, "/top20_marker_genes_tsne_", max(dim.use), "PC.csv")
)

pdf(
  paste0(output, "/Heatmap_all_cluster_tsne_", max(dim.use), "PC.pdf"),
  width = 25,
  height = 20
)

DoHeatmap(ScRNA, features = top5$gene, size = 5) +
  scale_fill_gradientn(colors = c("#437eb8", "white", "#FF3366")) +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text.x = element_text(size = 30)
  )

dev.off()

pdf(
  paste0(output, "/DotPlot_all_cluster_tsne_", max(dim.use), "PC.pdf"),
  width = 70,
  height = 10
)

DotPlot(ScRNA, features = unique(top5$gene)) +
  RotatedAxis() +
  scale_color_gradientn(colors = c('#FF9999', "white", "#FF3366")) +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

dev.off()

dpi <- 300

png(
  paste0(output, "/DotPlot_all_cluster_tsne_", max(dim.use), "PC.png"),
  w = 70 * dpi,
  h = 10 * dpi,
  units = "px",
  res = dpi,
  type = "cairo"
)

DotPlot(ScRNA, features = unique(top5$gene)) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("#FFCCCC", "white", "#FF3366")) +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

dev.off()


########### Manual Cell Annotation ###########

col <- c('#FF6666', '#E5D2DD', '#6A4C93', "#BC8F8F", '#FFCC99', '#FF9999',
         "#FFCCCC", '#4F6272', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699",
         "#6699CC", "#FFFFCC")

# Set output directory
setwd("./out(rna)/AT2/")
outdir <- "./out(rna)/AT2/"

output <- paste(outdir, "celltype", sep = "/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

file_path <- file.path(outdir, "ScRNA_clustered.rds")
scedata <- readRDS(file_path)

library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)

col <- c('#CCCCCC', "#FF3366", "#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC",
         "#CCFFCC", "#FFFFCC", '#E5D2DD', '#4F6272', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699",
         "#6699CC", "#FFFFCC")

file_path <- file.path(outdir, "ScRNA_clustered.rds")
scedata <- readRDS(file_path)

# Annotate clusters with cell type labels
scedata <- RenameIdents(scedata, c(
  "0" = "Cyfra 21−1 low",
  "1" = "Cyfra 21−1 low",
  "2" = "Cyfra 21−1 low",
  "3" = "Cyfra 21−1 low",
  "4" = "Cyfra 21−1 low",
  "5" = "Cyfra 21−1 low",
  "6" = "Cyfra 21−1 low",
  "7" = "Cyfra 21−1 low",
  "8" = "Cyfra 21−1 high",
  "9" = "Cyfra 21−1 low",
  "10" = "Cyfra 21−1 low",
  "11" = "Cyfra 21−1 low",
  "12" = "Cyfra 21−1 low",
  "13" = "Cyfra 21−1 high",
  "14" = "Cyfra 21−1 low",
  "15" = "Cyfra 21−1 low"
))

# Add cell type annotation to metadata
scedata$celltype <- scedata@active.ident
head(scedata@meta.data)

saveRDS(scedata, "celltype.rds")


library(ggsci)

# Plot UMAP colored by cell type
pdf(paste(output, "ann_umap.pdf", sep = "/"), width = 5.5, height = 5)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  cols = col,
  label.box = TRUE
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.position = "none",
    # legend.title = element_text(size = 18),
    # legend.text = element_text(size = 18),
    plot.title = element_blank()
  )

dev.off()

# legend.position = c(0.99, 0.12)
# legend.justification = c("right", "bottom")

# Plot UMAP colored by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = "/"), width = 5.5, height = 5)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
  label.size = 5,
  repel = TRUE,
  cols = col,
  label.box = TRUE
) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.position = "none",
    # legend.title = element_text(size = 18),
    # legend.text = element_text(size = 18),
    plot.title = element_blank()
  )

# legend.position = c(0.99, 0.12)
# legend.justification = c("right", "bottom")
dev.off()


pdf(
  paste(output, "ann-diff-umap.pdf", sep = "/"),
  width = 6 * length(unique(scedata$treatment)),
  height = 5
)

DimPlot(
  scedata,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

dev.off()

svg(
  paste(output, "ann-diff-umap.svg", sep = "/"),
  width = 6 * length(unique(scedata$treatment)),
  height = 5
)

DimPlot(
  scedata,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
) +
  theme(
    strip.text = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

dev.off()


####### Calculate Cell Proportions ###########

col <- c('#E5D2DD', "#FF3366", "#DC050C", "#FB8072", "#1965B0", "#7BAFDE",
         "#882E72", "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066",
         "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC")

output <- paste(outdir, "celltype", sep = "/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scedata <- readRDS(file_path)

table(scedata$seurat_clusters)

######## Calculate Cell Counts for Each Cell Type Across All Samples ########

cell_counts <- as.data.frame(table(Idents(scedata)))
colnames(cell_counts) <- c("CellType", "Counts")

# Sort cell types by cell count in descending order
cell_counts <- cell_counts[order(-cell_counts$Counts), ]

# Save cell counts for all cell types
write.csv(cell_counts, paste(output, "cell_counts.csv", sep = "/"), row.names = FALSE)

# Save the top 11 cell types
cell_counts_top11 <- head(cell_counts, 11)
write.csv(cell_counts_top11, paste(output, "cell_counts_top11.csv", sep = "/"), row.names = FALSE)

# Load required package
library(ggplot2)

# Plot cell type distribution bar plot
p <- ggplot(cell_counts, aes(x = reorder(CellType, -Counts), y = Counts, fill = CellType)) +
  geom_bar(stat = "identity") +
  labs(x = "Cell Type", y = "Counts", title = "") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  scale_fill_manual(values = col)

# Save plot
ggsave(paste(output, "cell_type_distribution.pdf", sep = "/"), plot = p, width = 7, height = 6, dpi = 800)
ggsave(paste(output, "cell_type_distribution.svg", sep = "/"), plot = p, width = 7, height = 6, dpi = 800)


# Calculate cell counts of each cell type in each group
# Count each cell type by sample
cell_counts_group <- as.data.frame(table(scedata$orig.ident, Idents(scedata)))
colnames(cell_counts_group) <- c("Sample", "CellType", "Counts")

# Add group information using treatment metadata
meta_data <- scedata@meta.data
group_info <- unique(meta_data[, c("orig.ident", "treatment")])
cell_counts_group <- merge(cell_counts_group, group_info, by.x = "Sample", by.y = "orig.ident")

# Calculate the proportion of each cell type within each sample
cell_counts_group <- cell_counts_group %>%
  group_by(Sample) %>%
  mutate(Ratio = Counts / sum(Counts))

p <- ggplot(cell_counts_group, aes(x = Sample, y = Counts, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = "", y = "Counts") +
  scale_fill_manual(values = col) +
  # scale_x_discrete(labels = c("WF-1", "WF-2")) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 10, height = 8, dpi = 800)

file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 10, height = 8, dpi = 800)


# Plot cell type proportions
p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = "", y = "Ratio") +
  scale_fill_manual(values = col) +
  # scale_x_discrete(labels = c("WF-1", "WF-2")) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 18, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 10, height = 8, dpi = 800)

file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 10, height = 8, dpi = 800)


############ Group-Level Cell Type Composition ############

cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate the proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))

# cell_counts_treatment$Treatment <- factor(cell_counts_treatment$Treatment, levels = c("0d", "3d", "7d", "14d"))

########## Plot Stacked Bar Plot of Cell Counts ##########

p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = "", y = "Counts") +
  scale_fill_manual(values = col) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 3.5 * length(unique(scedata$treatment)), height = 6, dpi = 800)

file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 3.5 * length(unique(scedata$treatment)), height = 6, dpi = 800)


########## Plot Stacked Bar Plot of Cell Proportions ##########

p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = '#222222') +
  theme_classic() +
  labs(x = "", y = "Ratio") +
  scale_fill_manual(values = col) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 3.5 * length(unique(scedata$treatment)), height = 6, dpi = 800)

file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 3.5 * length(unique(scedata$treatment)), height = 6, dpi = 800)


##### 19. Differential Expression Analysis Between Two Groups #####

library(scRNAtoolVis)
library(ggsci)
library(patchwork)
library(tidyverse)
library(ggrepel)

col <- c('#437eb8', '#FF6666', "#FFFFCC", '#FFCC99', '#FF9999',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300", "#FFCCCC",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699", "#6699CC")

# Create output directory for differential expression analysis
output <- paste(outdir, "differential_analysis", sep = "/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

file_path <- file.path(outdir, "celltype.rds")
scRNAsub <- readRDS(file_path)

# Identify differentially expressed genes between Cyfra 21−1 high and Cyfra 21−1 low groups
logFCfilter <- 0.25
adjPvalFilter <- 0.05

# Perform differential expression analysis
scRNAsub.cluster.markers <- FindMarkers(
  object = scRNAsub,
  ident.1 = "Cyfra 21−1 high",
  ident.2 = "Cyfra 21−1 low",
  group.by = "celltype",
  logfc.threshold = 0,
  min.pct = 0.25,
  test.use = "wilcox"
)

scRNAsub.cluster.markers$gene <- rownames(scRNAsub.cluster.markers)

# Add significance annotation
scRNAsub.cluster.markers <- scRNAsub.cluster.markers %>%
  mutate(
    Significance = ifelse(
      p_val_adj < adjPvalFilter & abs(avg_log2FC) > logFCfilter,
      ifelse(avg_log2FC > 0, "Up", "Down"),
      "Normal"
    )
  )

write.table(
  scRNAsub.cluster.markers,
  file = file.path(output, "sig.markers_Cyfra21_high_vs_low.txt"),
  sep = "\t",
  row.names = TRUE,
  quote = FALSE
)

saveRDS(scRNAsub.cluster.markers, file = file.path(output, "ScRNA.sig.markers.rds"))

# Save upregulated and downregulated genes separately
upregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")

downregulated_genes <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")

write.csv(
  upregulated_genes,
  file = file.path(output, "upregulated_genes_Cyfra21_high_vs_low.csv"),
  row.names = TRUE,
  quote = FALSE
)

write.csv(
  downregulated_genes,
  file = file.path(output, "downregulated_genes_Cyfra21_high_vs_low.csv"),
  row.names = TRUE,
  quote = FALSE
)

# Count upregulated and downregulated genes
upregulated_genes <- sum(scRNAsub.cluster.markers$Significance == "Up")
downregulated_genes <- sum(scRNAsub.cluster.markers$Significance == "Down")
total_diff_genes <- upregulated_genes + downregulated_genes

# Save upregulated and downregulated gene data frames
upregulated_genes_df <- scRNAsub.cluster.markers %>%
  filter(Significance == "Up")

downregulated_genes_df <- scRNAsub.cluster.markers %>%
  filter(Significance == "Down")

# Select top upregulated and downregulated genes for labeling
top_genes_upregulated <- upregulated_genes_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  head(15)

top_genes_downregulated <- downregulated_genes_df %>%
  filter(p_val_adj < 0.05 & avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  head(15)


# Define genes of interest
genes <- c(
  # Upregulated tumor-promoting genes
  "TP63", "HLA-DRB5", "PITPNC1", "PVT1", "FAT1",
  "MMP7", "TCF4", "CEACAM6", "ITGA2", "TNC", "EZR",
  "LGALS1", "FAT1", "ETS2",
  "FLNA", "FOSL2", "BCL2L1", "LMNA", "AHNAK",
  "S100A10", "TNS1", "OSMR", "PIGR", "ANXA2", "HSPG2",
  
  # Downregulated tumor-suppressive genes
  "DLC1", "CYLD", "PTEN", "CAVIN1", "CTNNA1",
  "NFKBIA", "WIF1", "HHIP", "RBBP6", "ATF3"
)

# Filter genes of interest
interested_genes <- scRNAsub.cluster.markers %>%
  filter(gene %in% genes)

# Plot volcano plot
p <- ggplot(scRNAsub.cluster.markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = Significance), size = 2, shape = 18) +
  scale_color_manual(values = c("#339999", "#FFCCCC", "#FF0066")) +
  geom_hline(yintercept = -log10(adjPvalFilter), linetype = "dashed") +
  geom_vline(xintercept = c(-logFCfilter, logFCfilter), linetype = "dashed") +
  # geom_text_repel(data = top_genes_upregulated, aes(label = gene),
  #                 size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  # geom_text_repel(data = top_genes_downregulated, aes(label = gene),
  #                 size = 4, fontface = "bold", max.overlaps = 50, box.padding = 0.6) +
  geom_text_repel(
    data = interested_genes,
    aes(label = gene),
    size = 5,
    fontface = "bold",
    max.overlaps = 50,
    box.padding = 0.6
  ) +
  theme_classic() +
  labs(
    title = "Cyfra 21−1 high vs Cyfra 21−1 low",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Significance"
  ) +
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = seq(-1.5, 2, by = 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 18)
  )

# Save volcano plot
ggsave(file.path(output, "Cyfra21_high_vs_low_volcano_plot.svg"), p, width = 8, height = 7, dpi = 300)
ggsave(file.path(output, "Cyfra21_high_vs_low_volcano_plot.pdf"), p, width = 8, height = 7, dpi = 300)


############### Expression of Differential Genes Across Cell Types ##################

library(Seurat)
library(tidyverse)
library(ggsci)

# Keep only genes present in the expression matrix
genes <- genes[genes %in% rownames(ScRNA)]

# Extract cell type information
ScRNA$celltype <- as.factor(ScRNA@meta.data$celltype)
treatment_groups <- levels(ScRNA$celltype)

# Extract expression matrix
expr_matrix <- GetAssayData(ScRNA, slot = "data")[genes, ]

avg_expr <- AverageExpression(ScRNA, features = genes, group.by = "celltype")$RNA
avg_expr_selected <- avg_expr[, treatment_groups]

# Save average expression table
avg_expr_df <- avg_expr_selected %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene")

write.table(
  avg_expr_df,
  file = paste0(output, "/differential_gene_expression.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Plot gene expression DotPlot grouped by cell type
plot <- DotPlot(ScRNA, features = unique(genes), group.by = "celltype") +
  RotatedAxis() +
  coord_flip() +
  scale_color_gradientn(colors = c('#CCCCCC', "white", "#FF3366")) +
  theme(
    axis.text = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20)
  )

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_by_celltype.pdf", sep = "/"), plot = plot, width = 7, height = 12)
ggsave(filename = paste(output, "marker_DotPlot_by_celltype.svg", sep = "/"), plot = plot, width = 7, height = 12)


##### GSEA and Enrichment Analysis After Differential Expression Analysis #####

# library(org.Mm.eg.db)  # Mouse annotation database
library(org.Hs.eg.db)    # Human annotation database
library(clusterProfiler)
library(enrichplot)
library(DOSE)

scRNAsub.cluster.markers <- readRDS(file.path(output, "ScRNA.sig.markers.rds"))


####### GO and KEGG Enrichment Analysis #######

# Extract upregulated and downregulated genes
gene_up <- scRNAsub.cluster.markers$gene[scRNAsub.cluster.markers$Significance == "Up"]
gene_down <- scRNAsub.cluster.markers$gene[scRNAsub.cluster.markers$Significance == "Down"]

# Convert SYMBOL to ENTREZID
gene_up_entrez <- as.character(na.omit(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_up,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  )[, 2]
))

gene_down_entrez <- as.character(na.omit(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gene_down,
    columns = "ENTREZID",
    keytype = "SYMBOL"
  )[, 2]
))

# Perform GO enrichment analysis
go_up <- enrichGO(
  gene = gene_up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

go_down <- enrichGO(
  gene = gene_down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

# Convert gene IDs from ENTREZID to SYMBOL
go_up@result$geneID <- sapply(strsplit(go_up@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )$SYMBOL
  paste(symbols, collapse = "/")
})

go_down@result$geneID <- sapply(strsplit(go_down@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )$SYMBOL
  paste(symbols, collapse = "/")
})

write.csv(as.data.frame(go_up), file = file.path(output, "go_up_results.csv"))
write.csv(as.data.frame(go_down), file = file.path(output, "go_down_results.csv"))

# Plot GO enrichment results
go_plot_up <- dotplot(go_up) + ggtitle("Upregulated Genes GO Enrichment")
go_plot_down <- dotplot(go_down) + ggtitle("Downregulated Genes GO Enrichment")

# Save GO enrichment plots
ggsave(file.path(output, "go_enrich_up_dot.pdf"), plot = go_plot_up, width = 6, height = 6)
ggsave(file.path(output, "go_enrich_up_dot.svg"), plot = go_plot_up, width = 6, height = 6)

ggsave(file.path(output, "go_enrich_down_dot.pdf"), plot = go_plot_down, width = 6, height = 5)
ggsave(file.path(output, "go_enrich_down_dot.svg"), plot = go_plot_down, width = 6, height = 5)

# Perform KEGG enrichment analysis
kegg_up <- enrichKEGG(
  gene = gene_up_entrez,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

kegg_down <- enrichKEGG(
  gene = gene_down_entrez,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

# Convert gene IDs from ENTREZID to SYMBOL
kegg_up@result$geneID <- sapply(strsplit(kegg_up@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )$SYMBOL
  paste(symbols, collapse = "/")
})

kegg_down@result$geneID <- sapply(strsplit(kegg_down@result$geneID, "/"), function(ids) {
  symbols <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ids,
    columns = "SYMBOL",
    keytype = "ENTREZID"
  )$SYMBOL
  paste(symbols, collapse = "/")
})

write.csv(as.data.frame(kegg_up), file = file.path(output, "kegg_up_results.csv"))
write.csv(as.data.frame(kegg_down), file = file.path(output, "kegg_down_results.csv"))

# Plot KEGG enrichment results
kegg_plot_up <- dotplot(kegg_up, showCategory = 20) +
  ggtitle("Upregulated Genes KEGG Enrichment")

kegg_plot_down <- dotplot(kegg_down, showCategory = 20) +
  ggtitle("Downregulated Genes KEGG Enrichment")

# Save KEGG enrichment plots
ggsave(file.path(output, "kegg_enrich_up_dot.pdf"), plot = kegg_plot_up, width = 7, height = 8)
ggsave(file.path(output, "kegg_enrich_up_dot.svg"), plot = kegg_plot_up, width = 7, height = 8)

ggsave(file.path(output, "kegg_enrich_down_dot.pdf"), plot = kegg_plot_down, width = 8, height = 8)
ggsave(file.path(output, "kegg_enrich_down_dot.svg"), plot = kegg_plot_down, width = 8, height = 8)

# Combine GO and KEGG enrichment plots
combined_plot_GO <- go_plot_up + go_plot_down + plot_layout(guides = "collect")
combined_plot_KEGG <- kegg_plot_up + kegg_plot_down + plot_layout(guides = "collect")

# Save combined plots
ggsave(file.path(output, "combined_GO_dot.pdf"), plot = combined_plot_GO, width = 13, height = 10)
ggsave(file.path(output, "combined_KEGG_dot.pdf"), plot = combined_plot_KEGG, width = 13, height = 12)


# Prepare visualization data
# Extract GO enrichment results
go_up_dt <- as.data.frame(go_up)
go_down_dt <- as.data.frame(go_down)

# Define GO terms of interest
interested_terms <- c(
  "small GTPase-mediated signal transduction",
  "protein localization to plasma membrane",
  "Wnt signaling pathway",
  "Rho protein signal transduction",
  "cellular response to epidermal growth factor stimulus",
  "regulation of autophagy",
  "response to epidermal growth factor",
  "ERK1 and ERK2 cascade",
  "positive regulation of MAPK cascade",
  "positive regulation of cell projection organization",
  "substrate-dependent cell migration",
  "receptor-mediated endocytosis",
  "epidermal growth factor receptor signaling pathway",
  "positive regulation of autophagy",
  "Notch signaling pathway",
  "intracellular receptor signaling pathway",
  "ERBB signaling pathway",
  "cellular response to fibroblast growth factor stimulus"
)

# Filter GO enrichment results by terms of interest
go_up_dt <- go_up_dt[go_up_dt$Description %in% interested_terms, ]

# Extract KEGG enrichment results
kegg_up_dt <- as.data.frame(kegg_up)
kegg_down_dt <- as.data.frame(kegg_down)

# Define KEGG pathways of interest
interested_terms <- c(
  "Focal adhesion",
  "Ubiquitin mediated proteolysis",
  "Regulation of actin cytoskeleton",
  "Hedgehog signaling pathway",
  "Wnt signaling pathway",
  "Phospholipase D signaling pathway",
  "Autophagy - animal",
  "ECM-receptor interaction",
  "Sphingolipid signaling pathway",
  "Proteoglycans in cancer",
  "Focal adhesion",
  "Cell adhesion molecules",
  "Adherens junction",
  "Rap1 signaling pathway",
  "Regulation of actin cytoskeleton",
  "Ubiquitin mediated proteolysis",
  "Leukocyte transendothelial migration",
  "Pathways in cancer",
  "Non-small cell lung cancer",
  "Small cell lung cancer",
  "PD-L1 expression and PD-1 checkpoint pathway in cancer"
)

# Filter KEGG enrichment results by pathways of interest
kegg_up_dt <- kegg_up_dt[kegg_up_dt$Description %in% interested_terms, ]

# Define colors for pathway classification
classification_colors <- c('#437eb8', '#FF6666', '#FFCC99', '#FF9999', '#80c5d8',
                           "#9999FF", "#FFCCCC", "#99CCFF", "#FF3366", "#CCCCFF",
                           "#CC0066", "#FFFFCC", "#66CCCC", "#FF9933", "#CCFFCC",
                           "#00CC66", "#99FFFF", "#FF3300", "#6699CC", "#CC99CC",
                           "#FF6699", "#FF0000", "#6666CC", "#FF9966", "#669999",
                           "#CC99FF", "#FFCCFF")

# Wrap long text labels
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

# Define bar plot function for GO or KEGG enrichment results
plot_GO_bar <- function(dt, title) {
  
  # Sort by p-value
  dt <- dt[order(-dt$pvalue, decreasing = TRUE), ]
  # dt <- dt[order(dt$Count, decreasing = TRUE), ]
  
  # Select top 20 pathways
  dt <- head(dt, 20)
  dt$Description <- factor(wrap_text(dt$Description), levels = wrap_text(dt$Description))
  
  # Left panel: enrichment p-value
  p1 <- ggplot(dt, aes(x = Description, y = log10(p.adjust), fill = Description)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = classification_colors) +
    coord_flip() +
    ylab("-log10(P-value)") +
    xlab("") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 16),
      axis.title.x = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    ggtitle(title)
  
  # Right panel: gene count
  p2 <- ggplot(dt, aes(x = Description, y = Count)) +
    geom_bar(stat = "identity", fill = "#66CCCC") +
    coord_flip() +
    ylab("Gene Count") +
    xlab("") +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 16),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      axis.title.x = element_text(size = 16, face = "bold"),
      plot.margin = margin(10, 10, 10, 10),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  # Combine two panels
  p_combined <- p1 + p2 + plot_layout(widths = c(2, 1.5))
  return(p_combined)
}

# Plot GO and KEGG enrichment bar plots
go_up_plot <- plot_GO_bar(go_up_dt, "Upregulated Genes GO Enrichment")
go_down_plot <- plot_GO_bar(go_down_dt, "Downregulated Genes GO Enrichment")

kegg_up_plot <- plot_GO_bar(kegg_up_dt, "Upregulated Genes KEGG Enrichment")
kegg_down_plot <- plot_GO_bar(kegg_down_dt, "Downregulated Genes KEGG Enrichment")

# Save enrichment bar plots
ggsave(file.path(output, "go_enrich_up_bar.pdf"), plot = go_up_plot, width = 8, height = 6)
ggsave(file.path(output, "go_enrich_down_bar.pdf"), plot = go_down_plot, width = 8, height = 6)

ggsave(file.path(output, "kegg_enrich_up_bar.pdf"), plot = kegg_up_plot, width = 8, height = 5)
ggsave(file.path(output, "kegg_enrich_down_bar.pdf"), plot = kegg_down_plot, width = 8, height = 5)


