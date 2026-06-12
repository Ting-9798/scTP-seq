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
# install.packages("tidydr")

# Set working directory
setwd("./data/")

# Define data directory
data_dir <- "./data/"

folders <- c("A")

# Initialize Seurat object list
seurat.list <- list()

# Loop through each subfolder, read expression matrix files, and perform Seurat preprocessing
for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # Read data in 10X format
  sample_data <- Read10X(data.dir = folder_path, gene.column = 2)
  sample_name <- basename(folder)
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = sample_data,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # Add mitochondrial and ribosomal gene percentages
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rp"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPL-")
  
  # Quality control filtering
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15 & percent.rp < 15
  )
  # seurat_obj <- subset(seurat_obj, subset = percent.mt < 25)
  
  # Data normalization and highly variable gene identification
  # seurat_obj <- NormalizeData(seurat_obj)
  # seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  
  # Add treatment column
  seurat_obj$treatment <- sample_name
  
  # Assign treatment group based on folder name
  if (str_detect(sample_name, "P-CK-1")) {
    seurat_obj$treatment <- "TIP-seq"
    seurat_obj$orig.ident <- "TIP-seq-1"
  } else if (str_detect(sample_name, "P-CK-2")) {
    seurat_obj$treatment <- "TIP-seq"
    seurat_obj$orig.ident <- "TIP-seq-2"
  } else if (str_detect(sample_name, "QJY-10X")) {
    seurat_obj$treatment <- "10X"
    seurat_obj$orig.ident <- "10X"
  }
  
  # Add processed Seurat object to the list
  seurat.list[[sample_name]] <- seurat_obj
}

# Integrate data
# anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:20)
# combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:20)

combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)
combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
View(combined_seurat@meta.data)

# Save merged Seurat object
saveRDS(combined_seurat, file = "./out/combined_seurat.rds")


####################### Seurat Analysis #####################

# Define color palette
col <- c(
  "#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD',
  '#4F6272', '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
  "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066", "#FF9933", "#CCFFCC",
  "#00CC66", "#99FFFF", "#FF3300", "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC",
  "#FF6699", "#6699CC", "#FFFFCC"
)

col <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
  "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066", "#00CC66",
  "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD', '#58A4C3',
  '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
  "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
  "#FF9933", "#CCFFCC", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC"
)

setwd("./out")
outdir <- "./out"

MergeOUT <- paste(outdir, "Merge", sep = "/")
dir.create(MergeOUT)
OUTPUT <- paste(MergeOUT, "Multiple_", sep = "/")

# Construct full file path
file_path <- file.path(outdir, "combined_seurat.rds")
ScRNA <- readRDS(file_path)


# Set output directory for scatter plots
plot_dir <- paste(outdir, "scatter_plots", sep = "/")
dir.create(plot_dir, showWarnings = FALSE)

# Extract all sample names from orig.ident
metadata <- ScRNA@meta.data
samples <- unique(metadata$treatment)

# Get count matrix as a sparse matrix
expr_data <- ScRNA@assays$RNA@counts

# Calculate total counts for each cell
total_counts <- Matrix::colSums(expr_data)

# Construct TP10K matrix in sparse format
TP10K <- t(t(expr_data) / total_counts) * 10000

# Convert to log2(TP10K + 1) while preserving sparsity
log2_TP10K <- log1p(TP10K) / log(2)


library(cowplot)

# Generate pairwise correlation scatter plots among groups
generate_combined_plots <- function(samples, plot_dir) {
  plot_list <- list()
  plot_count <- 1
  
  # Loop through all pairwise sample combinations
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      sample1 <- samples[i]
      sample2 <- samples[j]
      
      # Extract gene expression values from two samples
      sample1_data <- log2_TP10K[, metadata$treatment == sample1]
      sample2_data <- log2_TP10K[, metadata$treatment == sample2]
      
      # Calculate mean expression values
      sample1_mean <- rowMeans(sample1_data)
      sample2_mean <- rowMeans(sample2_data)
      
      # Create plotting data frame
      plot_data <- data.frame(Sample1 = sample1_mean, Sample2 = sample2_mean)
      
      # Calculate correlation coefficient
      correlation <- cor(sample1_mean, sample2_mean)
      
      # Fit linear regression model
      fit <- lm(Sample2 ~ Sample1, data = plot_data)
      slope <- coef(fit)[2]
      intercept <- coef(fit)[1]
      r2 <- summary(fit)$r.squared
      
      # Generate scatter plot
      plot <- ggplot(data = plot_data, aes(x = Sample1, y = Sample2)) +
        geom_hex(bins = 100) +
        scale_fill_gradientn(
          colors = c("#FF3366", "#1E90FF"),
          name = "Count"
        ) +
        # geom_abline(slope = 1, size = 1, intercept = 0, linetype = "dashed", color = "black") +
        geom_smooth(method = "lm", se = FALSE, color = "#DC050C", size = 1.2) +
        theme_minimal(base_size = 14) +
        labs(
          title = paste(sample1, "vs", sample2),
          subtitle = paste0(
            "R² = ", round(r2, 4),
            " | y = ", round(slope, 4), "x ",
            ifelse(intercept >= 0, "+ ", "- "),
            abs(round(intercept, 4))
          ),
          x = paste(sample1, "log2(TP10K)"),
          y = paste(sample2, "log2(TP10K)")
        ) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
          plot.subtitle = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          # legend.title = element_text(size = 20),
          # legend.text = element_text(size = 18),
          legend.position = "none",
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 2)
        )
      
      # Add current plot to the plot list
      plot_list[[plot_count]] <- plot
      plot_count <- plot_count + 1
    }
  }
  
  combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = "v")
  
  # Save output
  ggsave(
    filename = paste0(plot_dir, "/All_Samples_combined_scatter_plots.pdf"),
    limitsize = FALSE,
    plot = combined_plot,
    width = 14,
    height = 5.5 * ceiling(length(plot_list) / 3)
  )
  
  ggsave(
    filename = paste0(plot_dir, "/All_Samples_combined_scatter_plots.svg"),
    limitsize = FALSE,
    plot = combined_plot,
    width = 14,
    height = 5.5 * ceiling(length(plot_list) / 3)
  )
}

# Generate pairwise scatter plots for all samples
generate_combined_plots(samples, plot_dir)


# Generate violin plot for quality control metrics
pdf(paste(OUTPUT, "QC-VlnPlot.pdf"), width = 9, height = 6)
VlnPlot(
  ScRNA,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  group.by = "treatment",
  pt.size = 0,
  cols = col
)
dev.off()

# Generate box plots for quality control metrics
svg(paste(OUTPUT, "QC-BoxPlot.svg"), width = 8, height = 6)

p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, color = treatment)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, color = treatment)) +
  geom_boxplot(size = 1.2) +
  scale_color_manual(values = col) +
  labs(title = "nCount_RNA", x = "", y = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

CombinePlots(plots = list(p1, p2))
dev.off()


pdf(paste(OUTPUT, "QC-ViolinPlot.pdf"), width = 12, height = 6)

# Violin plot 1: nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nFeature_RNA", x = "", y = "") +
  # scale_fill_manual(values = c("#A4CDE1", "#67A4CC", "#277FB8")) +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Violin plot 2: nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "nCount_RNA", x = "", y = "") +
  # scale_fill_manual(values = c("#A4CDE1", "#67A4CC", "#277FB8")) +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Violin plot 3: percent.mt
p3 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = percent.mt, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  labs(title = "percent.mt", x = "", y = "") +
  # scale_fill_manual(values = c("#A4CDE1", "#67A4CC", "#277FB8")) +
  scale_fill_manual(values = col) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 18, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Combine plots
library(patchwork)
(p1 | p2 | p3) + plot_layout(ncol = 3)
dev.off()


# QC correlation plots: gene number, mitochondrial percentage, and RNA count
pdf(paste(OUTPUT, "cor-plot.pdf"), width = 15, height = 6)
plot1 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = col)
plot2 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = col)
CombinePlots(plots = list(plot1, plot2), legend = "right")
dev.off()

# Calculate cell cycle scores
pdf(paste(OUTPUT, "cellcycle.pdf"), width = 9, height = 6)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ScRNA <- CellCycleScoring(ScRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
VlnPlot(ScRNA, features = c("S.Score", "G2M.Score"), group.by = "treatment", pt.size = 1, cols = col)
dev.off()


#### 5. Expression Normalization ####

ScRNA <- NormalizeData(
  ScRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identify highly variable genes
ScRNA <- FindVariableFeatures(
  ScRNA,
  selection.method = "vst",
  nfeatures = 2000
)

# Visualize highly variable genes
pdf(paste(OUTPUT, "variable_gene.pdf"), width = 9, height = 6)
top10 <- head(VariableFeatures(ScRNA), 10)
plot1 <- VariableFeaturePlot(ScRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size = 3)
CombinePlots(plots = list(plot1, plot2), legend = "bottom")
dev.off()


#### 6. Scaling and PCA Dimensionality Reduction ####

# Scale data
ScRNA <- ScaleData(ScRNA)

# Run PCA
ScRNA <- RunPCA(ScRNA, npcs = 30)

pdf(paste(OUTPUT, "Dimplot.pdf"), width = 9, height = 6)
DimPlot(object = ScRNA, reduction = "pca", pt.size = .1, group.by = "treatment", cols = col)
dev.off()

pdf(paste(OUTPUT, "vlnplot.pdf"), width = 9, height = 6)
VlnPlot(object = ScRNA, features = "PC_1", group.by = "treatment", pt.size = 0, cols = col)
dev.off()

# PCA heatmap visualization
pdf(paste(OUTPUT, "DimHeatmap.pdf"), width = 9, height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Evaluate the number of principal components
pdf(paste0(OUTPUT, "PCA-ElbowPlot.pdf"), width = 6, height = 5)
ElbowPlot(ScRNA)
dev.off()

# Batch correction
# install.packages("harmony")
library(harmony)

ScRNA <- RunHarmony(
  ScRNA,
  group.by.vars = c("orig.ident"),
  plot_convergence = TRUE,
  verbose = TRUE
)

pdf(paste(OUTPUT, "Dimplot-correct.pdf"), width = 12, height = 6)
DimPlot(
  object = ScRNA,
  reduction = "harmony",
  pt.size = 0.1,
  group.by = "orig.ident"
)
dev.off()

pdf(paste(OUTPUT, "vlnplot-correct.pdf"), width = 12, height = 6)
VlnPlot(
  object = ScRNA,
  features = "harmony_1",
  group.by = "orig.ident",
  pt.size = 0
)
dev.off()

save(ScRNA, file = "ScRNA_batch_corrected_before_clustering.RData")


#### 7. Cell Clustering and Annotation ####

col <- c(
  '#FF6666', '#E5D2DD', '#6A4C93', "#BC8F8F", '#FFCC99', '#FF9999',
  "#FFCCCC", '#4F6272', '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3',
  '#E59CC4', '#437eb8', "#66CCCC", "#99CCFF", '#3399CC', "#FF3366",
  "#CC0066", "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699", "#6699CC",
  "#FFFFCC"
)

load("ScRNA_batch_corrected_before_clustering.RData")

# Cell clustering
# ScRNA <- ScRNA %>%
#   RunUMAP(dims = 1:20) %>%
#   RunTSNE(dims = 1:20) %>%
#   FindNeighbors(dims = 1:20)

ScRNA <- ScRNA %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  RunTSNE(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30)

ScRNA <- FindClusters(
  ScRNA,
  resolution = seq(from = 0.1, to = 1, by = 0.1)
)

# pdf(paste(OUTPUT, "clustree.pdf"), width = 10, height = 9)
# library(clustree)
# clustree(ScRNA)
# dev.off()

# Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.0.7"

# Select the resolution based on the clustering tree
ScRNA$seurat_clusters <- ScRNA@active.ident
table(Idents(ScRNA))

# Ensure treatment factor levels are ordered if needed
# ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Visualize clusters split by treatment
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(
  ScRNA,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "treatment",
  cols = col
)
dev.off()

# Generate individual UMAP plots
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)

DimPlot(ScRNA, reduction = "umap", label = TRUE, repel = TRUE, cols = col) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

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
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

dev.off()

# Visualize tumor and normal groups together
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


# Select top marker genes to visualize expression patterns across clusters
output <- paste(outdir, "cell_localization", sep = "/")
dir.create(output)

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

DoHeatmap(ScRNA, features = top5$gene, size = 2) +
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


########### Manual Cell Annotation ###########

col <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
  "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066", "#00CC66",
  "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD', '#58A4C3',
  '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
  "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
  "#FF9933", "#CCFFCC", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC"
)

setwd("./out/")
outdir <- "./out/"

output <- paste(outdir, "celltype", sep = "/")
dir.create(output)

file_path <- file.path(outdir, "ScRNA_clustered.rds")
scedata <- readRDS(file_path)

# Define marker genes for different cell types
cellmarker <- c(
  "IL7R", "CCR7",                                # Naive CD4+ T cells
  "CD14", "LYZ",                                # CD14+ monocytes
  "IL7R", "CD44", "S100A4", "NEAT1",           # Memory CD4+ T cells
  "MS4A1",                                      # B cells
  "IGLC2", "IGHA1", "MZB1", "JCHAIN",          # Plasma cells
  "CD8A", "CD8B",                              # CD8+ T cells
  "FCGR3A", "MS4A7",                            # FCGR3A+ monocytes
  "GNLY", "NKG7",                               # NK cells
  "FCER1A", "CST3", "LGALS2",                   # Dendritic cells
  "MYB", "STMN1", "CLEC4C", "CUX2", "IL3RA", "LILRA4", "TCF4",  # pDCs
  "PPBP",                                       # Platelets
  "KLRB1",                                      # Th17 cells
  "FOXO1", "FOXP3", "RTKN2", "IL2RA",          # Treg cells
  "CD25", "CTLA4", "TNFRSF18", "IL10", "TGFB", "GITR", "IKZF2",
  "CD3D", "CD3E", "CD3G",                      # T cells
  "FOSB", "LCP1", "ZFP36", "LRP",              # Monocytes
  "PCNA", "TOP2A", "CCNA2", "CDK1"             # Proliferating cells
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize immune cell marker gene expression using DotPlot
library(ggplot2)

plot <- DotPlot(scedata, features = unique(cellmarker)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 18),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c("#330066", "#336699", "#66CC66", "#FFCC33")
  )

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep = "/"), plot = plot, width = 14, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep = "/"), plot = plot, width = 14, height = 8)


library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)

col <- c(
  "#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
  '#E5D2DD', '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4',
  '#437eb8', "#99CCFF", '#3399CC', "#FF3366", "#CC0066", "#FF9933",
  "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300", "#6699CC", "#9999FF",
  "#CCCCFF", "#CC99CC", "#FF6699", "#6699CC", "#FFFFCC"
)

col <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
  "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066", "#00CC66",
  "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD', '#58A4C3',
  '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
  "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
  "#FF9933", "#CCFFCC", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC"
)

file_path <- file.path(outdir, "ScRNA_clustered.rds")
scedata <- readRDS(file_path)

# Annotate clusters with cell type labels
scedata <- RenameIdents(scedata, c(
  "0" = "CD14+ Mono",
  "1" = "CD8+ T",
  "2" = "Naive CD4+ T",
  "3" = "Naive CD4+ T",
  "4" = "NK",
  "5" = "NK",
  "6" = "CD14+ Mono",
  "7" = "B",
  "8" = "CD14+ Mono",
  "9" = "Memory CD4+ T",
  "10" = "Th17",
  "11" = "Naive CD4+ T",
  "12" = "B",
  "13" = "FCGR3A+ Mono",
  "14" = "CD14+ Mono",
  "15" = "NK",
  "16" = "DCs",
  "17" = "Treg",
  "18" = "NK",
  "19" = "Platelets",
  "20" = "pDCs",
  "21" = "CD14+ Mono",
  "22" = "B",
  "23" = "Plasma Cells"
))

# Add cell type annotation to metadata
scedata$celltype <- scedata@active.ident
head(scedata@meta.data)

# Extract UMAP coordinates and cell type annotations
umap_coords <- as.data.frame(Embeddings(scedata, reduction = "umap"))
umap_coords$celltype <- scedata$celltype
# umap_coords$CB <- scedata$CB

# Save as CSV file
write.csv(umap_coords, "celltype_umap.csv", row.names = TRUE)

saveRDS(scedata, "celltype.rds")


library(ggsci)

# Plot UMAP colored by cell type
pdf(paste(output, "ann_umap.pdf", sep = "/"), width = 7, height = 4)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.01,
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_blank()
  )

dev.off()

# legend.position = c(0.99, 0.12)
# legend.justification = c("right", "bottom")

# Plot UMAP colored by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = "/"), width = 7, height = 4)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.01,
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_blank()
  )

# legend.position = c(0.99, 0.12)
# legend.justification = c("right", "bottom")
dev.off()


pdf(
  paste(output, "ann-diff-umap.pdf", sep = "/"),
  width = 5.5 * length(unique(scedata$treatment)),
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
  width = 5.5 * length(unique(scedata$treatment)),
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


##### Add Total Cell Number for Each Group #####

# Count cells in each treatment group
cell_counts <- scedata@meta.data %>%
  as_tibble() %>%
  dplyr::count(treatment, name = "n") %>%
  mutate(label = sprintf("%s (n = %s cells)", treatment, format(n, big.mark = ",")))

# Construct named vector to replace facet labels
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# Save PDF output
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
  cols = col
) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
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

# Save SVG output
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
  cols = col
) +
  facet_wrap(~treatment, labeller = labeller(treatment = label_map)) +
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


############### Plot TIP-seq PBMC Data ###############

library(ggsci)

# Filter PBMC data from the TIP-seq group
pbmc_data <- subset(scedata, treatment == "TIP-seq")

# Plot UMAP for the TIP-seq PBMC group
pdf(paste(output, "pbmc_TIP_umap.pdf", sep = "/"), width = 7, height = 4)

DimPlot(
  object = pbmc_data,
  group.by = "celltype",
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_blank()
  )

dev.off()

# Plot UMAP for the TIP-seq PBMC group and save as SVG
svg(paste(output, "pbmc_TIP_umap.svg", sep = "/"), width = 7, height = 4)

DimPlot(
  object = pbmc_data,
  group.by = "celltype",
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_blank()
  )

dev.off()


############### Plot 10X PBMC Data ###############

library(ggsci)

# Filter PBMC data from the 10X group
pbmc_data <- subset(scedata, treatment == "10X")

# Plot UMAP for the 10X PBMC group
pdf(paste(output, "pbmc_10X_umap.pdf", sep = "/"), width = 7, height = 4)

DimPlot(
  object = pbmc_data,
  group.by = "celltype",
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_blank()
  )

dev.off()

# Plot UMAP for the 10X PBMC group and save as SVG
svg(paste(output, "pbmc_10X_umap.svg", sep = "/"), width = 7, height = 4)

DimPlot(
  object = pbmc_data,
  group.by = "celltype",
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
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    plot.title = element_blank()
  )

dev.off()


####### Calculate Cell Proportions ###########

col <- c(
  "#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
  '#E5D2DD', '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4',
  '#437eb8', "#99CCFF", '#3399CC', "#FF3366", "#CC0066", "#FF9933",
  "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300", "#6699CC", "#9999FF",
  "#CCCCFF", "#CC99CC", "#FF6699", "#6699CC", "#FFFFCC"
)

col <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
  "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066", "#00CC66",
  "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD', '#58A4C3',
  '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
  "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
  "#FF9933", "#CCFFCC", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC"
)

output <- paste(outdir, "celltype", sep = "/")
dir.create(output)

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
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = "#222222") +
  theme_classic() +
  labs(x = "", y = "Counts") +
  scale_fill_manual(values = col) +
  # scale_x_discrete(labels = c("WF-1", "WF-2")) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

# Add cell count labels
p <- p + geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7)

file_path <- paste0(output, "/genecount.pdf")
ggsave(file_path, plot = p, width = 4 * length(unique(scedata$orig.ident)), height = 8, dpi = 800)

file_path <- paste0(output, "/genecount.svg")
ggsave(file_path, plot = p, width = 4 * length(unique(scedata$orig.ident)), height = 8, dpi = 800)


p <- ggplot(cell_counts_group, aes(x = Sample, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = "#222222") +
  theme_classic() +
  labs(x = "", y = "Ratio") +
  scale_fill_manual(values = col) +
  # scale_x_discrete(labels = c("WF-1", "WF-2")) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 22, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22),
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

# Add cell proportion labels
p <- p + geom_text(
  aes(label = scales::percent(Ratio, accuracy = 0.1)),
  position = position_stack(vjust = 0.5),
  size = 7
)

file_path <- paste0(output, "/geneRatio.pdf")
ggsave(file_path, plot = p, width = 4 * length(unique(scedata$orig.ident)), height = 8, dpi = 800)

file_path <- paste0(output, "/geneRatio.svg")
ggsave(file_path, plot = p, width = 4 * length(unique(scedata$orig.ident)), height = 8, dpi = 800)


############ Group-Level Cell Type Composition ############

cell_counts_treatment <- as.data.frame(table(scedata$treatment, Idents(scedata)))
colnames(cell_counts_treatment) <- c("Treatment", "CellType", "Counts")

# Calculate the proportion of each cell type within each treatment group
cell_counts_treatment <- cell_counts_treatment %>%
  group_by(Treatment) %>%
  mutate(Ratio = Counts / sum(Counts))


########## Plot Stacked Bar Plot of Cell Counts ##########

p1 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Counts, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = "#222222") +
  theme_classic() +
  labs(x = "", y = "Counts") +
  scale_fill_manual(values = col) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 24, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

file_path <- paste0(output, "/genecount_treatment.pdf")
ggsave(file_path, plot = p1, width = 4 * length(unique(scedata$treatment)), height = 8, dpi = 800)

file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 4 * length(unique(scedata$treatment)), height = 8, dpi = 800)


########## Plot Stacked Bar Plot of Cell Proportions ##########

p2 <- ggplot(cell_counts_treatment, aes(x = Treatment, y = Ratio, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = "#222222") +
  theme_classic() +
  labs(x = "", y = "Ratio") +
  scale_fill_manual(values = col) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
    axis.text.x = element_text(size = 24, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 24),
    axis.title.y = element_text(size = 26),
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

file_path <- paste0(output, "/geneRatio_treatment.pdf")
ggsave(file_path, plot = p2, width = 4 * length(unique(scedata$treatment)), height = 8, dpi = 800)

file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 4 * length(unique(scedata$treatment)), height = 8, dpi = 800)


