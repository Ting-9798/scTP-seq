####################### Seurat Analysis #####################

# Set output color palette
col <- c("#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
         '#E5D2DD', '#4F6272', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699",
         "#6699CC", "#FFFFCC")

setwd("./out/")
outdir <- "./out/"

MergeOUT <- paste(outdir, "Merge", sep = "/")
dir.create(MergeOUT, recursive = TRUE, showWarnings = FALSE)
OUTPUT <- paste(MergeOUT, "Multiple_", sep = "/")

# Read input Seurat objects
This_study <- readRDS("./out/combined_seurat.rds")
Ref <- readRDS("./out/combined_seurat.rds")

# Add treatment information to metadata
This_study$treatment <- "This_study"
Ref$treatment <- "Ref"

# Merge This_study and Ref Seurat objects
combined_ScRNA <- merge(x = This_study, y = Ref)

# Check merged metadata
View(combined_ScRNA@meta.data)

# Save merged Seurat object
saveRDS(combined_ScRNA, file = paste(outdir, "combined_ScRNA.rds", sep = "/"))


# Read merged Seurat object
file_path <- file.path(outdir, "combined_ScRNA.rds")
ScRNA <- readRDS(file_path)

# Set output directory for scatter plots
plot_dir <- paste(outdir, "scatter_plots", sep = "/")
dir.create(plot_dir, showWarnings = FALSE)

# Extract all treatment groups
metadata <- ScRNA@meta.data
samples <- unique(metadata$treatment)

# Get raw count matrix as a sparse matrix
expr_data <- ScRNA@assays$RNA@counts

# Calculate total counts for each cell
total_counts <- Matrix::colSums(expr_data)

# Construct TP10K matrix in sparse format
TP10K <- t(t(expr_data) / total_counts) * 10000

# Convert to log2(TP10K + 1) while preserving sparsity
log2_TP10K <- log1p(TP10K) / log(2)


library(cowplot)

# Generate pairwise correlation scatter plots among treatment groups
generate_combined_plots <- function(samples, plot_dir) {
  plot_list <- list()
  plot_count <- 1
  
  # Loop through all pairwise sample combinations
  for (i in 1:(length(samples) - 1)) {
    for (j in (i + 1):length(samples)) {
      sample1 <- samples[i]
      sample2 <- samples[j]
      
      # Extract gene expression values for the two groups
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
          colors = c("#1E90FF", "#FF3366"),
          name = "Count"
        ) +
        # geom_abline(slope = 1, size = 1, intercept = 0,
        #             linetype = "dashed", color = "black")
        geom_smooth(method = "lm", se = FALSE, color = "#3366CC", size = 1.2) +
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
      
      # Add current plot to plot list
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


########## Add Mean Value Labels ##########

library(ggplot2)
library(dplyr)
library(patchwork)
# install.packages("shadowtext")
# library(shadowtext)

# Define color mapping to ensure consistent colors across treatment groups
treatment_levels <- unique(ScRNA@meta.data$treatment)
treatment_colors <- c()

treatment_colors <- sapply(treatment_levels, function(lvl) {
  if (grepl("Ref$", lvl)) {
    "#1E90FF"
  } else if (grepl("This_study$", lvl)) {
    "#FF3366"
  } else {
    "#9E9E9E"
  }
})

names(treatment_colors) <- treatment_levels

# Calculate mean gene number and mean UMI count for each treatment group
means_df <- ScRNA@meta.data %>%
  group_by(treatment) %>%
  summarise(
    mean_features = mean(nFeature_RNA, na.rm = TRUE),
    mean_counts = mean(nCount_RNA, na.rm = TRUE),
    .groups = "drop"
  )

# Format numbers as integers or one decimal place
fmt_num <- function(x) {
  ifelse(x >= 100, sprintf("%.0f", x), sprintf("%.1f", x))
}

# Save QC violin plot with mean labels
pdf(paste(OUTPUT, "QC-ViolinPlot.pdf"), width = 8, height = 5)

# Violin plot 1: nFeature_RNA
p1 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nFeature_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  geom_label(
    data = means_df,
    aes(
      x = treatment,
      y = mean_features,
      label = paste0(fmt_num(mean_features)),
      fill = NULL,
      color = treatment
    ),
    vjust = -0.6,
    size = 5,
    inherit.aes = FALSE,
    fontface = "bold",
    label.size = 0,
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines"),
    fill = "grey90"
  ) +
  labs(title = "nFeature_RNA", x = "", y = "") +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 5500)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Violin plot 2: nCount_RNA
p2 <- ggplot(data = ScRNA@meta.data, aes(x = treatment, y = nCount_RNA, fill = treatment)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 6, color = "black") +
  geom_label(
    data = means_df,
    aes(
      x = treatment,
      y = mean_counts,
      label = paste0(fmt_num(mean_counts)),
      fill = NULL,
      color = treatment
    ),
    vjust = -0.6,
    size = 5,
    inherit.aes = FALSE,
    fontface = "bold",
    label.size = 0,
    label.r = unit(0.1, "lines"),
    label.padding = unit(0.15, "lines"),
    fill = "grey90"
  ) +
  labs(title = "nCount_RNA", x = "", y = "") +
  scale_fill_manual(values = treatment_colors) +
  scale_color_manual(values = treatment_colors, guide = "none") +
  coord_cartesian(ylim = c(0, 40000)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 20, angle = 30, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 22, face = "bold"),
    plot.title = element_text(size = 24, hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# Combine plots
(p1 | p2) + plot_layout(ncol = 2)
dev.off()


# QC correlation plots between gene number, mitochondrial percentage, and RNA count
pdf(paste(OUTPUT, "cor-plot.pdf"), width = 15, height = 6)
plot1 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = col)
plot2 <- FeatureScatter(ScRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = col)
CombinePlots(plots = list(plot1, plot2), legend = "right")
dev.off()


# Calculate cell cycle scores
pdf(paste(OUTPUT, "cellcycle.pdf"), width = 9, height = 6)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ScRNA <- CellCycleScoring(
  ScRNA,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)
VlnPlot(
  ScRNA,
  features = c("S.Score", "G2M.Score"),
  group.by = "treatment",
  pt.size = 1,
  cols = col
)
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
DimPlot(
  object = ScRNA,
  reduction = "pca",
  pt.size = .1,
  group.by = "treatment",
  cols = col
)
dev.off()

pdf(paste(OUTPUT, "vlnplot.pdf"), width = 9, height = 6)
VlnPlot(
  object = ScRNA,
  features = "PC_1",
  group.by = "treatment",
  pt.size = 0,
  cols = col
)
dev.off()

# PCA heatmap visualization
pdf(paste(OUTPUT, "DimHeatmap.pdf"), width = 9, height = 6)
DimHeatmap(ScRNA, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

# Evaluate principal components
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

col <- c("#1E90FF", "#FF3366", "#DC050C", "#FB8072", "#1965B0", "#7BAFDE",
         "#882E72", "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
         '#E5D2DD', '#4F6272', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699",
         "#6699CC", "#FFFFCC")

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
  resolution = seq(from = 0.1, to = 1.0, by = 0.1)
)

# pdf(paste(OUTPUT, "clustree.pdf"), width = 10, height = 9)
# library(clustree)
# clustree(ScRNA)
# dev.off()

# Idents(ScRNA) <- "integrated_snn_res.0.7"
Idents(ScRNA) <- "RNA_snn_res.0.7"

# Select the resolution based on clustering structure
ScRNA$seurat_clusters <- ScRNA@active.ident
table(Idents(ScRNA))

# Set treatment factor order if needed
# ScRNA$treatment <- factor(ScRNA$treatment, levels = c("Non-infected", "Infected"))

# Visualize clusters split by treatment
pdf(paste(OUTPUT, "split.by_cluster_umap.pdf"), width = 10, height = 5)
DimPlot(
  ScRNA,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "treatment",
  cols = col,
  raster = FALSE
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

# Visualize clusters split by original sample
pdf(paste(OUTPUT, "split.by_cluster_umap_sample.pdf"), width = 30, height = 5)
DimPlot(
  ScRNA,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "orig.ident",
  cols = col,
  raster = FALSE
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

# Generate individual UMAP plots
pdf(paste(OUTPUT, "cluster_umap.pdf"), width = 6, height = 5)

DimPlot(
  ScRNA,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  cols = col,
  raster = FALSE
) +
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
  group.by = "treatment",
  raster = FALSE
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

# Visualize different treatment groups together
pdf(paste(OUTPUT, "cluster-diff_umap.pdf"), width = 6, height = 6)
DimPlot(
  ScRNA,
  repel = TRUE,
  reduction = "umap",
  group.by = "treatment",
  raster = FALSE
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

pdf(
  paste0(output, "/DotPlot_all_cluster_tsne_", max(dim.use), "PC.pdf"),
  width = 100,
  height = 10
)

DotPlot(ScRNA, features = unique(top5$gene)) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("#FF9999", "white", "#FF3366")) +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

dev.off()

dpi <- 300

png(
  paste0(output, "/DotPlot_all_cluster_tsne_", max(dim.use), "PC.png"),
  w = 100 * dpi,
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

col <- c("#1E90FF", "#FF3366", "#DC050C", "#FB8072", "#1965B0", "#7BAFDE",
         "#882E72", "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1", '#FF9999', "#66CCCC", "#FFCCCC", "#CCFFCC", "#FFFFCC",
         '#E5D2DD', '#4F6272', '#58A4C3',
         '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699",
         "#6699CC", "#FFFFCC")

setwd("./out/")
outdir <- "./out/"

output <- paste(outdir, "celltype", sep = "/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

file_path <- file.path(outdir, "ScRNA_clustered.rds")
scedata <- readRDS(file_path)

# Define marker genes for different cell types
cellmarker <- c(
  "PECAM1", "VWF", "CDH5",                         # Endothelial cells
  "PROX1", "LYVE1", "CCL21",                      # Lymphatic endothelial cells
  # "PDGFRB", "RGS5",                              # Pericytes
  # "SULF2", "ATP1B3", "KRT19", "WFDC2",          # Cancer cells
  # "FMNL2", "PCSK2", "CACNA2D1",                 # Neuroendocrine cells
  # "PSD3", "FOXI1", "CFTR", "ATP6V0D2",
  # "ATP6V1B1", "ASCL3", "TP63", "EZR",           # Pulmonary ionocytes
  "TPSAB1", "TPSB2",                               # Mast cells
  # "PCNA", "TOP2A", "CDK1",                      # Proliferating cells
  "COL1A1", "COL1A2", "DCN",                      # Fibroblasts / stromal cells
  # "MYH11", "CNN1",                               # Smooth muscle cells
  # "KRT5", "TP63", "KRT14", "ITGA6", "ITGB4", "NGFR",  # Basal cells
  "CD79A", "MS4A1",                                # B cells
  "MZB1", "JCHAIN",                                # Plasma cells
  "GZMB", "LILRA4", "SPIB",                       # pDCs
  "HLA-DPB1", "HLA-DRA", "HLA-DRB1",              # Dendritic cells
  "CD68", "LGALS3", "ITGAM", "PPARG",             # Macrophages
  "FCN1", "CD300E", "NLRP3", "TBC1D8",            # Monocytes
  "S100A8", "S100A9",                              # Neutrophils
  "NKG7", "CCL5", "KLRB1", "GZMA", "KLRF1", "PRF1",  # NK cells
  "BRCA1", "MCM6", "HELLS", "CDT1", "DTL",        # NK T cells
  # "IL2RA",                                       # Treg cells
  # "CD8A", "CD8B",                                # CD8+ T cells
  "GZMA", "GZMB", "IFNG", "CCL4", "CCL5", "IL2", "TBX21",  # Cytotoxic T cells
  "TRBC1", "TRBC2", "CCL5", "CD2", "CD3E", "CD3G",          # T cells
  # "CCR7", "CD3D", "CD3E", "CD4", "CD8A",
  # "SELL", "TCF7", "LEF1",                       # Naive T cells
  # "S100A8", "S100A9", "IL1R2",                  # Neutrophils
  "SCGB1A1", "SCGB3A2",                            # Club cells
  "KCNE1", "FOXJ1", "TPPP3", "TUBB4B", "TUBB", "TP73", "CCDC7",  # Ciliated cells
  "DNAH5", "DNAH9", "TEKT1",                       # Ciliated epithelial cells
  "EPCAM", "KRT18", "KRT19",                       # Epithelial cells
  "AGER", "CAV1", "CLIC5", "HOPX", "SEMA3E", "COL4A3",  # AT1 cells
  "LAMP3", "ABCA3", "SLC34A2"                      # AT2 cells
)

cellmarker <- cellmarker[cellmarker %in% rownames(scedata)]

# Visualize marker gene expression using DotPlot
library(ggplot2)

plot <- DotPlot(scedata, features = unique(cellmarker)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 1, vjust = 0.5, angle = 90, size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(x = NULL, y = NULL) +
  guides(size = guide_legend(order = 3)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colours = c("#330066", "#336699", "#66CC66", "#FFCC33")
  )

# Save DotPlot
ggsave(filename = paste(output, "marker_DotPlot_1.pdf", sep = "/"), plot = plot, width = 15, height = 8)
ggsave(filename = paste(output, "marker_DotPlot_1.svg", sep = "/"), plot = plot, width = 15, height = 8)


library("Seurat")
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(ggsci)

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066",
         "#00CC66", "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD',
         '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#99FFFF", "#FF3300",
         "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC")

# file_path <- file.path(outdir, "ScRNA_clustered.rds")
# scedata <- readRDS(file_path)

# Annotate clusters with cell type labels
scedata <- RenameIdents(scedata, c(
  "0" = "T",
  "1" = "NK",
  "2" = "Club",
  "3" = "Neutrophils",
  "4" = "T",
  "5" = "Macrophages",
  "6" = "Endothelial",
  "7" = "Mast",
  "8" = "Fibroblasts",
  "9" = "B",
  "10" = "DCs",
  "11" = "DCs",
  "12" = "DCs",
  "13" = "Fibroblasts",
  "14" = "Club",
  "15" = "AT2",
  "16" = "Plasma",
  "17" = "Macrophages",
  "18" = "Endothelial",
  "19" = "Monocytes",
  "20" = "AT1",
  "21" = "T",
  "22" = "Ciliated",
  "23" = "Fibroblasts",
  "24" = "AT2",
  "25" = "Lymphatic endothelial",
  "26" = "AT2",
  "27" = "NK T",
  "28" = "pDCs",
  "29" = "AT1",
  "30" = "Mast",
  "31" = "Fibroblasts",
  "32" = "Mast",
  "33" = "Fibroblasts",
  "34" = "Lymphatic endothelial"
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
pdf(paste(output, "ann_umap.pdf", sep = "/"), width = 7, height = 6)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
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
    legend.position = "none",
    # legend.title = element_text(size = 18),
    # legend.text = element_text(size = 18),
    plot.title = element_blank()
  )

dev.off()

# legend.position = c(0.99, 0.12)
# legend.justification = c("right", "bottom")


# Plot UMAP colored by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = "/"), width = 10, height = 6)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
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
    legend.position = "none",
    # legend.title = element_text(size = 18),
    # legend.text = element_text(size = 18),
    plot.title = element_blank()
  )

# legend.position = c(0.99, 0.12)
# legend.justification = c("right", "bottom")
dev.off()


pdf(paste(output, "ann-diff-umap.pdf", sep = "/"), width = 12, height = 5)
DimPlot(
  scedata,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
)
dev.off()

svg(paste(output, "ann-diff-umap.svg", sep = "/"), width = 12, height = 5)
DimPlot(
  scedata,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
)
dev.off()


########## Add Total Cell Number ##########

library(ggsci)

# Calculate total cell number
total_cells <- ncol(scedata)

# Construct title label
title_label <- paste("( n =", total_cells, "cells )")

# Plot UMAP colored by cell type and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = "/"), width = 7, height = 6)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
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
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  ggtitle(title_label)

dev.off()

# Plot UMAP colored by cell type and save as SVG
svg(paste(output, "ann_umap.svg", sep = "/"), width = 7, height = 6)

DimPlot(
  object = scedata,
  group.by = "celltype",
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
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
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  ) +
  ggtitle(title_label)

dev.off()


###### Add Cell Type Labels ##########

library(ggsci)
library(Seurat)
library(dplyr)
library(ggrepel)

# Calculate total cell number
total_cells <- ncol(scedata)

# Construct title label
title_label <- paste("( n =", total_cells, "cells )")

# Define cell type colors
celltype_colors <- col

# Extract UMAP coordinates
umap_coords <- Embeddings(scedata, "umap")
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype

# Create a data frame containing the center coordinates of each cell type
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Calculate center coordinates for each cell type
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Calculate the boundary range of each cell type to guide label placement
celltype_ranges <- umap_df %>%
  group_by(celltype) %>%
  summarise(
    min_x = min(UMAP1),
    max_x = max(UMAP1),
    min_y = min(UMAP2),
    max_y = max(UMAP2),
    width = max_x - min_x,
    height = max_y - min_y
  )

# Merge center coordinates with boundary information
celltype_labels <- left_join(celltype_centers, celltype_ranges, by = "celltype")

# Determine label placement direction based on cluster shape
celltype_labels <- celltype_labels %>%
  mutate(
    direction_x = ifelse(width > height, 1, 0),
    direction_y = ifelse(height > width, 1, 0),
    nudge_x = ifelse(UMAP1 > mean(UMAP1), 2, -2),
    nudge_y = ifelse(UMAP2 > mean(UMAP2), 2, -2)
  )

# Plot UMAP with cell type labels and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = "/"), width = 7, height = 6)

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels,
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7,
    fontface = "bold",
    box.padding = 1.5,
    point.padding = 0.8,
    nudge_x = celltype_labels$nudge_x,
    nudge_y = celltype_labels$nudge_y,
    min.segment.length = 1,
    segment.size = 0.5,
    segment.color = "grey40",
    segment.alpha = 0.7,
    force = 2,
    max.iter = 10000,
    direction = "both",
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  ggtitle(title_label)

dev.off()

# Plot UMAP with cell type labels and save as SVG
svg(paste(output, "ann_umap.svg", sep = "/"), width = 7, height = 6)

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
  geom_point(size = 0.1) +
  ggrepel::geom_text_repel(
    data = celltype_labels,
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 7,
    fontface = "bold",
    box.padding = 1.5,
    point.padding = 0.8,
    nudge_x = celltype_labels$nudge_x,
    nudge_y = celltype_labels$nudge_y,
    min.segment.length = 1,
    segment.size = 0.5,
    segment.color = "grey40",
    segment.alpha = 0.7,
    force = 2,
    max.iter = 10000,
    direction = "both",
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  ggtitle(title_label)

dev.off()


##### Add Borders to Text Labels #####

library(ggsci)
library(Seurat)
library(dplyr)
library(ggrepel)

# Calculate total cell number
total_cells <- ncol(scedata)

# Construct title label
title_label <- paste("( n =", total_cells, "cells )")

# Define cell type colors
celltype_colors <- col

# Extract UMAP coordinates
umap_coords <- Embeddings(scedata, "umap")
umap_data <- as.data.frame(umap_coords)
umap_data$celltype <- scedata$celltype
umap_data$label <- rownames(umap_data)

# Create a data frame containing the center coordinates of each cell type
umap_df <- as.data.frame(Embeddings(scedata, "umap"))
umap_df$celltype <- scedata$celltype
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# Calculate center coordinates for each cell type
celltype_centers <- umap_df %>%
  group_by(celltype) %>%
  summarise(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

# Plot UMAP with bordered cell type labels and save as PDF
pdf(paste(output, "ann_umap.pdf", sep = "/"), width = 7, height = 6)

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype, label = label)) +
  geom_point(size = 0.1) +
  ggrepel::geom_label_repel(
    data = celltype_centers,
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey30",
    label.size = 0.5,
    label.r = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  ggtitle(title_label)

dev.off()

# Plot UMAP with bordered cell type labels and save as SVG
svg(paste(output, "ann_umap.svg", sep = "/"), width = 7, height = 6)

ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = celltype, label = label)) +
  geom_point(size = 0.1) +
  ggrepel::geom_label_repel(
    data = celltype_centers,
    aes(x = UMAP1, y = UMAP2, label = celltype, color = celltype),
    size = 5,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey30",
    label.size = 0.5,
    label.r = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = celltype_colors) +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = 2, hjust = 0.03),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20)
  ) +
  ggtitle(title_label)

dev.off()


##### Add Total Cell Number for Each Group #####

# Count cells in each treatment group
cell_counts <- scedata@meta.data %>%
  group_by(treatment) %>%
  summarise(n = n()) %>%
  mutate(label = paste0(treatment, " ( n = ", n, " cells)"))

# Construct named vector to replace facet labels
label_map <- setNames(cell_counts$label, cell_counts$treatment)

# Save PDF output
pdf(
  paste(output, "ann-diff-umap.pdf", sep = "/"),
  width = 6.5 * length(unique(scedata$treatment)),
  height = 5
)

DimPlot(
  scedata,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  cols = col,
  raster = FALSE
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
  width = 6.5 * length(unique(scedata$treatment)),
  height = 5
)

DimPlot(
  scedata,
  reduction = "umap",
  split.by = "treatment",
  pt.size = 0.1,
  label = FALSE,
  cols = col,
  raster = FALSE
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


# Select target samples
selected_samples <- c("T1", "T2", "T3", "T4")

# Filter cells with orig.ident belonging to T1, T2, T3, and T4
filtered_data <- subset(scedata, subset = orig.ident %in% selected_samples)

# Plot selected samples and save as PDF
pdf(paste(output, "ann-diff-umap-sample-selected.pdf", sep = "/"), width = 22, height = 5)

DimPlot(
  filtered_data,
  reduction = "umap",
  split.by = "orig.ident",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
)

dev.off()

# Plot selected samples and save as SVG
svg(paste(output, "ann-diff-umap-sample-selected.svg", sep = "/"), width = 22, height = 5)

DimPlot(
  filtered_data,
  reduction = "umap",
  split.by = "orig.ident",
  pt.size = 0.1,
  label = FALSE,
  label.size = 5,
  repel = TRUE,
  cols = col
)

dev.off()


####### Calculate Cell Proportions ###########

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066",
         "#00CC66", "#CC99CC", "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD',
         '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
         "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
         "#FF9933", "#CCFFCC", "#99FFFF", "#FF3300",
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
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = "#222222") +
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
  geom_bar(stat = "identity", width = 0.9, size = 0.5, colour = "#222222") +
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
ggsave(file_path, plot = p1, width = 4 * length(unique(scedata$treatment)), height = 7, dpi = 800)

file_path <- paste0(output, "/genecount_treatment.svg")
ggsave(file_path, plot = p1, width = 4 * length(unique(scedata$treatment)), height = 7, dpi = 800)


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
ggsave(file_path, plot = p2, width = 4 * length(unique(scedata$treatment)), height = 7, dpi = 800)

file_path <- paste0(output, "/geneRatio_treatment.svg")
ggsave(file_path, plot = p2, width = 4 * length(unique(scedata$treatment)), height = 7, dpi = 800)



