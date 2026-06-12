## Visualization
rm(list = ls())

library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize)
library(reshape2)
library(cowplot)

setwd("/pyscenic/")

## Input files
input_clustered_rds <- "./ScRNA_clustered.rds"
input_celltype_rds <- "./celltype.rds"
input_loom <- "out_SCENIC.loom"

## Output files
scenic_input_csv <- "for_scenic_data.csv"
rss_visual_rdata <- "for_rss_and_visual.Rdata"



############################
## Prepare files for pySCENIC
############################

scenic_input <- readRDS(input_clustered_rds)

# The count matrix must be transposed before running pySCENIC.
# Otherwise, pySCENIC may report an input-format error.
write.csv(
  t(as.matrix(scenic_input@assays$RNA@counts)),
  file = scenic_input_csv
)



############################
## 1. Extract information from out_SCENIC.loom
############################

loom <- open_loom(input_loom)

regulons_incid_mat <- get_regulons(loom, column.attr.name = "Regulons")
regulons_incid_mat[1:4, 1:4]

regulons <- regulonsToGeneLists(regulons_incid_mat)

regulonAUC <- get_regulons_AUC(
  loom,
  column.attr.name = "RegulonsAUC"
)

regulonAucThresholds <- get_regulon_thresholds(loom)

tail(
  regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))]
)

embeddings <- get_embeddings(loom)

close_loom(loom)

rownames(regulonAUC)
names(regulons)



############################
## 2. Read Seurat object and match cells with SCENIC results
############################

seurat.data <- readRDS(input_celltype_rds)
seurat.data

seurat.data@meta.data$Barcode <- colnames(seurat.data)

# Keep cells that are present in the regulon AUC matrix.
seurat.data <- subset(
  seurat.data,
  Barcode %in% colnames(regulonAUC)
)

seurat.data

table(seurat.data@meta.data$Seu_Clusters)

DimPlot(
  seurat.data,
  reduction = "umap",
  label = TRUE
)



############################
## 3. AUC visualization preparation
############################

sub_regulonAUC <- regulonAUC[, match(colnames(seurat.data), colnames(regulonAUC))]

dim(sub_regulonAUC)

# Confirm that the cell order is consistent.
identical(colnames(sub_regulonAUC), colnames(seurat.data))

cellClusters <- data.frame(
  row.names = colnames(seurat.data),
  seurat_clusters = as.character(seurat.data$T_cell_activation_3class)
)

# Add treatment information to each cell type.
seurat.data@meta.data$Celltype_Group <- paste0(
  seurat.data@meta.data$T_cell_activation_3class,
  "_",
  seurat.data@meta.data$treatment
)

table(seurat.data@meta.data$Celltype_Group)

cellTypes <- data.frame(
  row.names = colnames(seurat.data),
  celltype = seurat.data$T_cell_activation_3class
)

Celltype_Group <- data.frame(
  row.names = colnames(seurat.data),
  celltype = seurat.data$Celltype_Group
)

head(cellTypes)
head(Celltype_Group)

sub_regulonAUC[1:4, 1:4]

# Save intermediate objects for RSS analysis and visualization.
save(
  sub_regulonAUC,
  cellTypes,
  Celltype_Group,
  cellClusters,
  seurat.data,
  file = rss_visual_rdata
)



############################
## 4.1. Mean transcription factor activity
############################

# Calculate the average regulon activity across different single-cell subgroups.
selectedResolution <- "celltype"

cellsPerGroup <- split(
  rownames(cellTypes),
  cellTypes[, selectedResolution]
)

# Remove extended regulons.
sub_regulonAUC <- sub_regulonAUC[
  onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),
]

dim(sub_regulonAUC)

# Calculate mean regulon activity by group.
regulonActivity_byGroup <- sapply(
  cellsPerGroup,
  function(cells) {
    rowMeans(getAUC(sub_regulonAUC)[, cells])
  }
)

# Scale regulon activity.
# The scale function normalizes columns, so the matrix is transposed first.
regulonActivity_byGroup_Scaled <- t(
  scale(
    t(regulonActivity_byGroup),
    center = TRUE,
    scale = TRUE
  )
)

dim(regulonActivity_byGroup_Scaled)

regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

Heatmap(
  regulonActivity_byGroup_Scaled,
  name = "z-score",
  col = colorRamp2(
    seq(from = -2, to = 2, length = 11),
    rev(brewer.pal(11, "Spectral"))
  ),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE
)



############################
## 4.2. RSS analysis for identifying cell-type-specific TFs
############################

rss <- calcRSS(
  AUC = getAUC(sub_regulonAUC),
  cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution]
)

rss <- na.omit(rss)

# Convert RSS results to matrix.
rss <- as.matrix(rss)

# Plot RSS results.
# Setting order_rows = FALSE keeps the original row order.
rssPlot <- plotRSS(
  rss,
  order_rows = FALSE
)

rss_treatment <- calcRSS(
  AUC = getAUC(sub_regulonAUC),
  cellAnnotation = Celltype_Group[colnames(sub_regulonAUC), "celltype"]
)

# Optional interactive visualization.
# rssPlot <- plotRSS(rss_treatment)
# plotly::ggplotly(rssPlot$plot)



############################
## 5. Select top cell-type-specific TFs based on RSS Z-score
############################

rss_data <- rssPlot$plot$data

head(rss_data)
dim(rss_data)

# Convert RSS plot data to wide format.
rss_data <- dcast(
  rss_data,
  Topic ~ cellType,
  value.var = "Z"
)

rownames(rss_data) <- rss_data[, 1]
rss_data <- rss_data[, -1]

head(rss_data)
dim(rss_data)
colnames(rss_data)

# Assign each regulon to the cluster where it has the highest Z-score.
rss_summary <- data.frame(
  TF = rownames(rss_data),
  Max_Zscore = apply(rss_data, 1, max),
  MaxCluster = colnames(rss_data)[apply(rss_data, 1, which.max)]
)

head(rss_summary)
table(rss_summary$MaxCluster)

rss_summary$gene <- gsub(
  "\\(\\+\\)|\\(\\-\\)",
  "",
  rss_summary$TF
)

rss_summary$regulation <- ifelse(
  str_detect(rss_summary$TF, "\\(\\+\\)"),
  "Pos",
  "Neg"
)

head(rss_summary)

write.csv(
  rss_summary,
  file = "rss_summary.csv",
  row.names = FALSE
)

top_positive_TFs_df <- rss_summary[rss_summary$regulation == "Pos", ]

top_positive_TFs_df <- top_positive_TFs_df %>%
  group_by(MaxCluster) %>%
  slice_max(order_by = Max_Zscore, n = 5, with_ties = FALSE)

head(top_positive_TFs_df)
table(top_positive_TFs_df$MaxCluster)

plot_celltypes <- unique(top_positive_TFs_df$MaxCluster)

top_positive_TFs_df <- top_positive_TFs_df %>%
  mutate(
    MaxCluster = factor(
      MaxCluster,
      levels = plot_celltypes,
      ordered = TRUE
    )
  ) %>%
  arrange(MaxCluster)

top_positive_TFs <- top_positive_TFs_df[
  top_positive_TFs_df$MaxCluster %in% plot_celltypes,
]$TF

top_positive_TFs

# Visualize RSS values for selected TFs.
sub_rssPlot <- plotRSS(
  rss[top_positive_TFs, plot_celltypes, drop = FALSE],
  varName = "Cell_type",
  thr = 0.01,
  col.low = "#330066",
  col.mid = "#66CC66",
  col.high = "#FFCC33"
)

sub_rssPlot <- sub_rssPlot$plot

ggsave(
  "rss_dotplot.pdf",
  plot = sub_rssPlot,
  width = 3,
  height = 4
)



############################
## 6. Heatmap of top regulon activity
############################

rss_scaled <- regulonActivity_byGroup_Scaled

head(rss_scaled)

df <- do.call(
  rbind,
  lapply(seq_len(ncol(rss_scaled)), function(i) {
    data.frame(
      regulon = rownames(rss_scaled),
      cluster = colnames(rss_scaled)[i],
      activity_in_cluster = rss_scaled[, i],
      activity_in_other_clusters = apply(rss_scaled[, -i, drop = FALSE], 1, median)
    )
  })
)

df$activity_difference <- df$activity_in_cluster - df$activity_in_other_clusters

top5 <- df %>%
  group_by(cluster) %>%
  slice_max(order_by = activity_difference, n = 5, with_ties = FALSE) %>%
  ungroup()

# Remove duplicated regulons to avoid duplicated row names in the annotation.
top5 <- top5[!duplicated(top5$regulon), ]

row_annotation <- data.frame(
  cluster = top5$cluster
)

rownames(row_annotation) <- top5$regulon

heatmap_matrix <- rss_scaled[top5$regulon, , drop = FALSE]

# Heatmap color scale.
heatmap_colors <- colorRampPalette(
  c("#4DBBD5FF", "white", "#C71000FF")
)(100)

# Row annotation colors.
col40 <- c(
  "#B22222", "#F4A460", "#FED439FF", "#91D1C2FF", "#79AF97FF",
  "#00A087FF", "#4DBBD5FF", "#3B4992CC", "#3C5488FF", "#8968CD",
  "#CD96CD", "#B09C85FF", "#CD8C95", "#FD8CC1FF", "#eacb85",
  "#FFF68F", "#A2CD5A", "#6E8B3D", "#20B2AA", "#6CA6CD",
  "#3A5FCD", "#925E9FFF", "#7D26CD", "#8B475D", "#008B45CC",
  "#631879CC", "#008280CC", "#BB0021CC", "#5F559BCC", "#A20056CC",
  "#808180CC", "#1B1919CC", "#FF6F00FF", "#C71000FF", "#008EA0FF",
  "#8A4198FF", "#5A9599FF", "#FF6348FF", "#84D7E1FF", "#F7B6D2CC"
)

unique_clusters <- unique(row_annotation$cluster)

annotation_colors <- list(
  cluster = setNames(
    col40[seq_along(unique_clusters)],
    unique_clusters
  )
)

pdf("regulon_activity_heatmap.pdf", width = 4, height = 3)

ComplexHeatmap::pheatmap(
  heatmap_matrix,
  annotation_row = row_annotation,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 10,
  fontsize_col = 12,
  heatmap_legend_param = list(
    title = "Relative value",
    title_gp = gpar(
      fontsize = 11,
      fontface = "bold",
      col = "black"
    )
  ),
  color = heatmap_colors,
  annotation_colors = annotation_colors,
  border_color = NA
)

dev.off()



############################
## 7. Dot plot of TF activity
############################

# Add regulon AUC matrix to Seurat metadata.
add_regulon_auc_to_metadata <- function(seurat_obj, auc_obj) {
  auc_mat <- getAUC(auc_obj)
  auc_mat <- auc_mat[, colnames(seurat_obj), drop = FALSE]
  
  new_cols <- setdiff(rownames(auc_mat), colnames(seurat_obj@meta.data))
  
  if (length(new_cols) > 0) {
    seurat_obj@meta.data <- cbind(
      seurat_obj@meta.data,
      t(auc_mat[new_cols, , drop = FALSE])
    )
  }
  
  return(seurat_obj)
}

seurat.data <- add_regulon_auc_to_metadata(
  seurat.data,
  sub_regulonAUC
)

top_positive_TFs <- top_positive_TFs_df$TF

# Calculate average TF activity from metadata.
calculate_tf_activity_from_metadata <- function(
    seurat_obj,
    tf_list,
    group_var = "treatment"
) {
  metadata <- seurat_obj@meta.data
  
  tf_columns <- tf_list[tf_list %in% colnames(metadata)]
  
  if (length(tf_columns) == 0) {
    stop("None of the selected TFs were found in the metadata.")
  }
  
  groups <- unique(metadata[[group_var]])
  
  result <- matrix(
    NA,
    nrow = length(tf_columns),
    ncol = length(groups)
  )
  
  rownames(result) <- tf_columns
  colnames(result) <- groups
  
  for (tf in tf_columns) {
    for (group in groups) {
      group_cells <- which(metadata[[group_var]] == group)
      result[tf, group] <- mean(metadata[group_cells, tf], na.rm = TRUE)
    }
  }
  
  return(result)
}

# Calculate TF activity by treatment.
tf_activity <- calculate_tf_activity_from_metadata(
  seurat.data,
  top_positive_TFs,
  group_var = "treatment"
)

print("TF activity matrix:")
print(tf_activity)

# Sort TFs by Tumor - Normal difference.
if ("Tumor" %in% colnames(tf_activity) & "Normal" %in% colnames(tf_activity)) {
  
  tf_differences <- tf_activity[, "Tumor"] - tf_activity[, "Normal"]
  sorted_tfs <- names(sort(tf_differences, decreasing = TRUE))
  
  # Set treatment group order.
  treatment_order <- c("Tumor", "Normal")
  seurat.data$treatment <- factor(
    seurat.data$treatment,
    levels = treatment_order
  )
  
  cat("TF ranking based on Tumor - Normal activity difference:\n")
  
  for (i in seq_along(sorted_tfs)) {
    tf <- sorted_tfs[i]
    tumor_act <- tf_activity[tf, "Tumor"]
    normal_act <- tf_activity[tf, "Normal"]
    diff <- tumor_act - normal_act
    
    cat(
      sprintf(
        "%2d. %s: Tumor = %.3f, Normal = %.3f, Difference = %.3f\n",
        i,
        tf,
        tumor_act,
        normal_act,
        diff
      )
    )
  }
  
} else {
  
  cat("Detected treatment groups:\n")
  print(colnames(tf_activity))
  
  sorted_tfs <- top_positive_TFs
}

# Dot plot of selected TFs.
p_sorted <- DotPlot(
  seurat.data,
  features = sorted_tfs,
  group.by = "T_cell_activation_3class"
) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      hjust = 1,
      vjust = 1,
      angle = 45,
      size = 12
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 12
    ),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  coord_flip() +
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0
  )

print(p_sorted)

ggsave(
  "TF_activity_by_treatment_sorted_dotplot.pdf",
  plot = p_sorted,
  width = 5,
  height = 5
)



############################
## 8. Select the top five positive TFs for each treatment group and draw dot plots
############################

# Calculate RSS values for each cell-type-treatment group.
treatment_rss <- calcRSS(
  AUC = getAUC(sub_regulonAUC),
  cellAnnotation = Celltype_Group[colnames(sub_regulonAUC), "celltype"]
)

treatment_rss <- na.omit(treatment_rss)

# Extract positive regulons.
positive_tfs <- rownames(treatment_rss)[
  grepl("\\(\\+\\)", rownames(treatment_rss))
]

treatment_rss_pos <- treatment_rss[positive_tfs, , drop = FALSE]

# Get treatment groups.
treatment_groups <- unique(seurat.data@meta.data$treatment)

# Select the top five positive TFs for each treatment group.
top_tfs_per_treatment <- list()

for (treatment_group in treatment_groups) {
  
  treatment_columns <- grep(
    paste0("_", treatment_group, "$"),
    colnames(treatment_rss_pos),
    value = TRUE
  )
  
  treatment_data <- treatment_rss_pos[, treatment_columns, drop = FALSE]
  
  if (ncol(treatment_data) > 0) {
    
    if (ncol(treatment_data) > 1) {
      treatment_means <- rowMeans(treatment_data, na.rm = TRUE)
    } else {
      treatment_means <- treatment_data[, 1]
    }
    
    sorted_tfs <- names(sort(treatment_means, decreasing = TRUE))
    top_tfs <- head(sorted_tfs, 5)
    
    top_tfs_per_treatment[[treatment_group]] <- top_tfs
    
    cat("Treatment:", treatment_group, "\n")
    cat("Top 5 positive TFs:", paste(top_tfs, collapse = ", "), "\n")
    cat(
      "RSS values:",
      paste(round(treatment_means[top_tfs], 3), collapse = ", "),
      "\n\n"
    )
  }
}

# Merge top TFs across treatment groups and remove duplicates.
all_top_tfs <- unique(unlist(top_tfs_per_treatment))

cat(
  "Total unique top TFs across all treatments:",
  length(all_top_tfs),
  "\n"
)

print(all_top_tfs)

# If too few TFs are selected, supplement with globally high-RSS positive TFs.
if (length(all_top_tfs) < 5) {
  
  overall_means <- rowMeans(treatment_rss_pos, na.rm = TRUE)
  
  additional_tfs <- names(
    sort(overall_means, decreasing = TRUE)
  )[seq_len(min(10, length(overall_means)))]
  
  additional_tfs <- setdiff(additional_tfs, all_top_tfs)
  
  all_top_tfs <- c(
    all_top_tfs,
    head(additional_tfs, 5 - length(all_top_tfs))
  )
}

# Make sure regulon AUC values are available in metadata.
seurat.data <- add_regulon_auc_to_metadata(
  seurat.data,
  sub_regulonAUC
)

# Set treatment order.
treatment_order <- sort(unique(seurat.data@meta.data$treatment))

seurat.data$treatment <- factor(
  seurat.data$treatment,
  levels = treatment_order
)

# Dot plot without axis flipping.
p_treatment <- DotPlot(
  seurat.data,
  features = all_top_tfs,
  group.by = "T_cell_activation_3class"
) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      hjust = 1,
      angle = 45,
      vjust = 1,
      size = 12
    ),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    ),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0,
    name = "Average\nExpression"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "Percent\nExpressed"
  )

print(p_treatment)

ggsave(
  "TF_activity_by_treatment_top5_dotplot_horizontal.pdf",
  plot = p_treatment,
  width = 8,
  height = 2.8
)

# Dot plot with flipped axes.
p_treatment_flipped <- DotPlot(
  seurat.data,
  features = all_top_tfs,
  group.by = "T_cell_activation_3class"
) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      hjust = 1,
      angle = 45,
      vjust = 1,
      size = 14
    ),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold"
    ),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    color = "Average\nExpression",
    size = "Percent\nExpressed"
  ) +
  scale_color_gradient2(
    low = "#4DBBD5FF",
    mid = "white",
    high = "#C71000FF",
    midpoint = 0,
    name = "Average\nExpression"
  ) +
  scale_size_continuous(
    range = c(2, 6),
    name = "Percent\nExpressed"
  ) +
  coord_flip()

print(p_treatment_flipped)

ggsave(
  "TF_activity_by_treatment_top5_dotplot.pdf",
  plot = p_treatment_flipped,
  width = 5.5,
  height = 5
)



############################
## 9. Optional: draw separate dot plots for each treatment group
############################

for (treatment_group in treatment_groups) {
  
  if (!is.null(top_tfs_per_treatment[[treatment_group]])) {
    
    treatment_tfs <- top_tfs_per_treatment[[treatment_group]]
    
    treatment_subset <- subset(
      seurat.data,
      treatment == treatment_group
    )
    
    if (length(unique(treatment_subset$celltype)) > 1) {
      
      p_single <- DotPlot(
        treatment_subset,
        features = treatment_tfs,
        group.by = "celltype"
      ) +
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(
            hjust = 1,
            vjust = 1,
            angle = 45,
            size = 10
          ),
          axis.text.y = element_text(
            face = "italic",
            size = 12
          ),
          plot.title = element_text(
            hjust = 0.5,
            size = 12,
            face = "bold"
          )
        ) +
        labs(
          x = NULL,
          y = NULL,
          title = paste("Top TFs in", treatment_group)
        ) +
        scale_color_gradient2(
          low = "#4DBBD5FF",
          mid = "white",
          high = "#C71000FF",
          midpoint = 0
        )
      
      ggsave(
        paste0("TF_activity_", treatment_group, "_by_celltype_dotplot.pdf"),
        plot = p_single,
        width = 6,
        height = 5
      )
    }
  }
}

