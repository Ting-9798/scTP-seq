# Clear the R environment
rm(list = ls())

########### Protein Data Processing #############

# Load required packages
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)
library(grDevices)

# install.packages("tidydr")

setwd("./pro/")

# Define the data directory
data_dir <- "./pro/"

folders <- c("A")

# Initialize the Seurat object list
seurat.list <- list()

# Loop through each subfolder, read the expression matrix, and perform Seurat analysis
suffix_counter <- 1  # Used to generate a unique suffix

for (folder in folders) {
  folder_path <- file.path(data_dir, folder)
  
  # Read 10X-format data
  sample_data <- Read10X(data.dir = folder_path)
  sample_name <- basename(folder)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = sample_data,
    project = sample_name,
    assay = "Protein"
  )
  
  # Add the percentage of mitochondrial features
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Quality control filtering
  # seurat_obj <- subset(seurat_obj, subset = percent.mt < 25)
  seurat_obj <- subset(seurat_obj, subset = nCount_Protein < 600 & percent.mt < 25)
  
  # Normalize data and identify highly variable features
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = 2000
  )
  
  # Assign molecule group
  seurat_obj$mole <- "Pro"
  
  # Add the processed Seurat object to the list
  seurat.list[[sample_name]] <- seurat_obj
}

# Merge all Seurat objects
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat.list)

combined_seurat@meta.data$CB <- rownames(combined_seurat@meta.data)
# View(combined_seurat@meta.data)

# Save the merged Seurat object
saveRDS(combined_seurat, file = "./out/combined_seurat.rds")





################# Visualize Each Protein Marker #################

# Clear the R environment
rm(list = ls())

####################### Merge Transcriptomic and Protein Data for UMAP Visualization ############

col <- c(
  '#FF6666', '#E5D2DD', '#6A4C93', "#BC8F8F", '#FFCC99', '#FF9999',
  "#FFCCCC", '#4F6272', '#58A4C3', '#F9BB72', '#F3B1A0', '#57C3F3',
  '#E59CC4', '#437eb8', "#66CCCC", "#99CCFF", '#3399CC', "#FF3366",
  "#CC0066", "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699", "#6699CC",
  "#FFFFCC"
)

# Load required packages
library(Seurat)
library(dplyr)
library(data.table)
library(stringr)
library(Matrix)
library(ggplot2)
library(tidydr)

# Create output directory
output <- "./out/Merge/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Load transcriptomic and protein data
rna_data <- readRDS("./out/ScRNA.rds")
pro_data <- readRDS("./out/combined_seurat.rds")

# Ensure that the CB column is consistent between RNA and protein data
rna_data@meta.data$CB <- rownames(rna_data@meta.data)
pro_data@meta.data$CB <- rownames(pro_data@meta.data)

# Extract UMAP coordinates and nCount information
rna_umap <- data.frame(
  CB = rownames(Embeddings(rna_data, "umap")),
  umap_1 = Embeddings(rna_data, "umap")[, 1],
  umap_2 = Embeddings(rna_data, "umap")[, 2],
  nCount_RNA = rna_data@meta.data$nCount_RNA
)

# Define protein marker names
protein_dict <- list(
  "barcode-1" = "AKR1B10",
  "barcode-2" = "CD276",
  "barcode-3" = "FGL1",
  "barcode-4" = "IL4I1"
)

# Extract protein count data
protein_data <- as.data.frame(GetAssayData(pro_data, assay = "Protein", slot = "counts"))
protein_data$PB <- rownames(protein_data)

library(tidyr)

# Convert protein data to long format
protein_data <- protein_data %>%
  gather(key = "CB", value = "counts", -PB)

# Annotate protein names according to the PB column
protein_data$PB <- sapply(protein_data$PB, function(pb) {
  if (pb %in% names(protein_dict)) {
    return(protein_dict[[pb]])
  } else {
    return("Others")
  }
})

# Summarize protein expression values for each CB
pro_umap <- protein_data %>%
  group_by(CB, PB) %>%
  summarize(nCount_Protein = sum(counts), .groups = "drop")

# Merge transcriptomic and protein data
merged_data <- merge(rna_umap, pro_umap, by = "CB", all.x = TRUE)

# Replace NA values with 0
merged_data[is.na(merged_data)] <- 0

# merged_data <- merged_data[complete.cases(merged_data), ]

# Replace 0 values in the PB column with "Others"
merged_data$PB <- ifelse(merged_data$PB == 0, "Others", merged_data$PB)


# Plot transcriptomic UMAP using RNA feature counts
# Sort cells by expression level from low to high
merged_data <- merged_data[order(merged_data$nCount_RNA), ]

pdf(file.path(output, "UMAP_RNA_FeaturePlot(single)11.pdf"), width = 6.5, height = 5)

ggplot(merged_data, aes(x = umap_1, y = umap_2, color = nCount_RNA)) +
  geom_point(size = 0.8) +
  scale_color_gradientn(
    colors = c("#E5D2DD", "#0099CC", "#003399"),
    limits = c(0, max(merged_data$nCount_RNA, na.rm = TRUE))
  ) +
  # scale_color_gradientn(colors = c("#E5D2DD", "#FF9999", "#0099CC")) +
  theme_minimal() +
  theme_dr(
    xlength = 0.15,
    ylength = 0.15,
    arrow = arrow(length = unit(0.1, "inches"), type = "closed")
  ) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", hjust = 0.03),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "", x = "UMAP_1", y = "UMAP_2", color = "nCount RNA")

dev.off()


# Plot protein UMAP using protein feature counts
# Sort cells by expression level from low to high
merged_data <- merged_data[order(merged_data$nCount_Protein), ]

pdf(file.path(output, "UMAP_Protein_FeaturePlot(single)11.pdf"), width = 6, height = 5)

ggplot(merged_data, aes(x = umap_1, y = umap_2, color = nCount_Protein)) +
  geom_point(size = 0.8) +
  scale_color_gradientn(colors = c("#CCCCCC", "#33CCCC", "#3333FF", "#FF6666")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),       # Remove axis text
    axis.ticks = element_blank(),      # Remove axis ticks
    axis.title = element_blank(),
    # axis.title = element_text(face = "bold", hjust = 0.03),
    legend.position = "none"
    # legend.title = element_text(size = 14),
    # legend.text = element_text(size = 12)
  ) +
  labs(title = "", x = "UMAP_1", y = "UMAP_2", color = "nCount Protein")

dev.off()


# Merge treatment information
treatment_info <- data.frame(
  CB = rownames(pro_data@meta.data),
  treatment = pro_data@meta.data$treatment
)

merged_data <- merge(merged_data, treatment_info, by = "CB", all.x = TRUE)

# Replace NA values with "Unknown"
merged_data$treatment[is.na(merged_data$treatment)] <- "Unknown"

# Loop through each treatment group and plot UMAP separately
group_list <- unique(merged_data$treatment)

for (group in group_list) {
  group_data <- merged_data[merged_data$treatment == group, ]
  
  # Sort cells by protein expression level
  group_data <- group_data[order(group_data$nCount_Protein), ]
  
  pdf(
    file.path(output, paste0("UMAP_Protein_FeaturePlot_", group, ".pdf")),
    width = 6,
    height = 5
  )
  
  print(
    ggplot(group_data, aes(x = umap_1, y = umap_2, color = nCount_Protein)) +
      geom_point(size = 0.8) +
      scale_color_gradientn(colors = c("#CCCCCC", "#33CCCC", "#3333FF", "#FF6666")) +
      theme_minimal() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none"
      ) +
      labs(title = "", color = "nCount Protein")
  )
  
  dev.off()
}





################# Plot Protein Expression UMAP for Each Protein Separately ###############

# Get all protein markers from the PB column and remove "Others"
protein_markers <- unique(merged_data$PB)
protein_markers <- protein_markers[!protein_markers %in% c("Others")]

# Define protein markers to plot
protein_markers <- c("AKR1B10", "CD276", "FGL1", "IL4I1")

library(gridExtra)

# Create an empty list to store plots
plots <- list()

# Loop through each protein marker and generate UMAP plots
for (protein in protein_markers) {
  
  # Subset data for the current protein
  protein_data <- merged_data[merged_data$PB == protein, ]
  
  # Sort cells by protein expression level from low to high
  protein_data <- protein_data[order(protein_data$nCount_Protein), ]
  
  p <- ggplot(protein_data, aes(x = umap_1, y = umap_2, color = nCount_Protein)) +
    geom_point(size = 0.8) +
    scale_color_gradientn(colors = c("#CCCCCC", "#33CCCC", "#3333FF", "#FF6666")) +
    theme_minimal() +
    theme_dr(
      xlength = 0.15,
      ylength = 0.15,
      arrow = arrow(length = unit(0.1, "inches"), type = "closed")
    ) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", hjust = 0.03),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    ) +
    labs(
      title = paste(protein),
      x = "UMAP_1",
      y = "UMAP_2",
      color = "nCount Protein"
    )
  
  # Store the plot object in the list
  plots[[protein]] <- p
}

# Combine all protein marker UMAP plots into a four-column layout
pdf(file.path(output, "UMAP_Protein_FeaturePlot_combined.pdf"), width = 22, height = 5)
do.call(grid.arrange, c(plots, ncol = 4))
dev.off()


