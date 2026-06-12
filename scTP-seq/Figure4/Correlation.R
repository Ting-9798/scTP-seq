######################## Extract RNA and Protein Data for Correlation Analysis ########################

########### Calculate Correlations by Treatment Group ##############

# Clear the R environment
rm(list = ls())

# Load required R packages
library(Seurat)
library(dplyr)

# Set output directory
output <- "./out(T)/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Load RNA and protein data
rna_data <- readRDS("./out(T)/celltype.rds")
# View(rna_data@meta.data)
pro_data <- readRDS("./out(T)/combined_seurat.rds")

# Subset data by treatment group
rna_data_Tumor <- subset(rna_data, subset = treatment == "Tumor")
rna_data_Normal <- subset(rna_data, subset = treatment == "Normal")
pro_data_Tumor <- subset(pro_data, subset = treatment == "Tumor")
pro_data_Normal <- subset(pro_data, subset = treatment == "Normal")

# Extract expression matrices
rna_Tumor_expr <- as.data.frame(
  t(as.matrix(GetAssayData(rna_data_Tumor, assay = "RNA", layer = "data")))
)
rna_Normal_expr <- as.data.frame(
  t(as.matrix(GetAssayData(rna_data_Normal, assay = "RNA", layer = "data")))
)

pro_Tumor_expr <- as.data.frame(
  t(as.matrix(GetAssayData(pro_data_Tumor, assay = "Protein", slot = "data")))
)
pro_Normal_expr <- as.data.frame(
  t(as.matrix(GetAssayData(pro_data_Normal, assay = "Protein", slot = "data")))
)

# Define target protein list
target_proteins <- c("AKR1B10", "CD276", "FGL1", "IL4I1")

# Define protein barcode dictionary
protein_dict <- list(
  "barcode-1" = "AKR1B10",
  "barcode-2" = "CD276",
  "barcode-3" = "FGL1",
  "barcode-4" = "IL4I1"
)

# Function to rename protein columns
rename_proteins <- function(df, dict) {
  colnames(df) <- sapply(
    colnames(df),
    function(x) ifelse(x %in% names(dict), dict[[x]], x)
  )
  return(df)
}

pro_Tumor_expr <- rename_proteins(pro_Tumor_expr, protein_dict)
pro_Normal_expr <- rename_proteins(pro_Normal_expr, protein_dict)

# Extract target protein expression values if present
common_proteins <- intersect(target_proteins, colnames(pro_Tumor_expr))
pro_Tumor_expr <- pro_Tumor_expr[, common_proteins, drop = FALSE]
pro_Normal_expr <- pro_Normal_expr[, common_proteins, drop = FALSE]

# Ensure cell barcode consistency between RNA and protein data
common_cells_Tumor <- intersect(rownames(rna_Tumor_expr), rownames(pro_Tumor_expr))
common_cells_Normal <- intersect(rownames(rna_Normal_expr), rownames(pro_Normal_expr))

rna_Tumor_expr <- rna_Tumor_expr[common_cells_Tumor, , drop = FALSE]
pro_Tumor_expr <- pro_Tumor_expr[common_cells_Tumor, , drop = FALSE]
rna_Normal_expr <- rna_Normal_expr[common_cells_Normal, , drop = FALSE]
pro_Normal_expr <- pro_Normal_expr[common_cells_Normal, , drop = FALSE]

# Function to calculate the gene-protein correlation matrix
calc_correlation_matrix <- function(rna_mat, pro_mat, protein_list) {
  cor_matrix <- matrix(NA, nrow = length(protein_list), ncol = ncol(rna_mat))
  rownames(cor_matrix) <- protein_list
  colnames(cor_matrix) <- colnames(rna_mat)
  
  for (protein in protein_list) {
    for (gene in colnames(rna_mat)) {
      cor_value <- tryCatch({
        cor(rna_mat[, gene], pro_mat[, protein], method = "pearson")
      }, error = function(e) NA)
      
      cor_matrix[protein, gene] <- cor_value
    }
  }
  
  return(as.data.frame(cor_matrix))
}

# Calculate correlation matrices
cor_Tumor <- calc_correlation_matrix(rna_Tumor_expr, pro_Tumor_expr, common_proteins)
cor_Normal <- calc_correlation_matrix(rna_Normal_expr, pro_Normal_expr, common_proteins)

# Transpose results: rows represent genes and columns represent proteins
cor_Tumor_t <- t(cor_Tumor)
cor_Normal_t <- t(cor_Normal)

# Save correlation results
write.csv(cor_Tumor_t, file = file.path(output, "Tumor_Correlation.csv"))
write.csv(cor_Normal_t, file = file.path(output, "Normal_Correlation.csv"))


########### Calculate Gene-Protein Correlations for Each Protein After Merging Samples ############

# Clear the R environment
rm(list = ls())

# Load required R packages
library(Seurat)
library(dplyr)

# Set output directory
output <- "./out(T)/by_group/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Load RNA and protein data
rna_data <- readRDS("./out(T)/celltype.rds")
# View(rna_data@meta.data)
pro_data <- readRDS("./out(T)/combined_seurat.rds")

# Define sample groups
sample_groups <- list(
  Tumor = "Tumor",
  Normal = "Normal"
)

# Define target proteins
target_proteins <- c("AKR1B10", "CD276", "FGL1", "IL4I1")

# Define protein barcode dictionary
protein_dict <- list(
  "barcode-1" = "AKR1B10",
  "barcode-2" = "CD276",
  "barcode-3" = "FGL1",
  "barcode-4" = "IL4I1"
)

# Function to rename protein columns
rename_proteins <- function(df, dict) {
  colnames(df) <- sapply(
    colnames(df),
    function(x) ifelse(x %in% names(dict), dict[[x]], x)
  )
  return(df)
}

# Function to calculate correlation vector
calc_correlation_vector <- function(rna_mat, protein_vec) {
  sapply(rna_mat, function(gene_expr) {
    tryCatch({
      cor(gene_expr, protein_vec, method = "pearson")
    }, error = function(e) NA)
  })
}

# Initialize result list
protein_results <- list()
for (p in target_proteins) {
  protein_results[[p]] <- list()
}

# Loop through merged sample groups
for (sample_name in names(sample_groups)) {
  message("Processing ", sample_name)
  
  group_ids <- sample_groups[[sample_name]]
  
  # Subset RNA and protein data by treatment group
  rna_subset <- subset(rna_data, subset = treatment %in% group_ids)
  pro_subset <- subset(pro_data, subset = treatment %in% group_ids)
  
  # Extract expression matrices
  rna_expr <- as.data.frame(
    t(as.matrix(GetAssayData(rna_subset, assay = "RNA", layer = "data")))
  )
  pro_expr <- as.data.frame(
    t(as.matrix(GetAssayData(pro_subset, assay = "Protein", slot = "data")))
  )
  
  # Transform, rename, and filter expression matrices
  rna_expr <- log10(rna_expr + 1)
  pro_expr <- log10(pro_expr + 1)
  pro_expr <- rename_proteins(pro_expr, protein_dict)
  
  common_proteins <- intersect(target_proteins, colnames(pro_expr))
  pro_expr <- pro_expr[, common_proteins, drop = FALSE]
  
  common_cells <- intersect(rownames(rna_expr), rownames(pro_expr))
  rna_expr <- rna_expr[common_cells, , drop = FALSE]
  pro_expr <- pro_expr[common_cells, , drop = FALSE]
  
  if (nrow(rna_expr) == 0 || ncol(pro_expr) == 0) {
    warning("Skipping ", sample_name, ": no overlapping cells or proteins.")
    next
  }
  
  # Calculate correlation vector for each protein
  for (protein in common_proteins) {
    cor_vec <- calc_correlation_vector(rna_expr, pro_expr[[protein]])
    protein_results[[protein]][[sample_name]] <- cor_vec
  }
}

# Summarize and export results
for (protein in names(protein_results)) {
  result_list <- protein_results[[protein]]
  
  if (length(result_list) == 0) next
  
  cor_df <- do.call(cbind, result_list)
  cor_df <- na.omit(cor_df)
  
  write.csv(
    cor_df,
    file = file.path(output, paste0(protein, "_Correlation_ByGroup.csv")),
    row.names = TRUE
  )
}


###################### Scatter Plot with Protein-Specific Colors #######################

library(ggplot2)
library(ggrepel)

input_dir <- "./out(T)/by_group/"
plot_dir <- file.path(input_dir, "plots_40_patient")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

protein_files <- list.files(input_dir, pattern = "_Correlation_ByGroup.csv$", full.names = TRUE)
sample_pairs <- combn(c("Normal", "Tumor"), 2, simplify = FALSE)

# Define fixed color order
fixed_colors <- c("#58A4C3", "#A4CDE1", "#FF9999", "#66CCCC", "#FFCCCC", "#CCFFCC")

# Extract protein names and assign colors
protein_names <- gsub("_Correlation_ByGroup.csv", "", basename(protein_files))
n_proteins <- length(protein_names)
protein_colors <- setNames(fixed_colors[1:n_proteins], protein_names)

for (file in protein_files) {
  protein_name <- sub("_Correlation_ByGroup.csv", "", basename(file))
  cor_df <- read.csv(file, row.names = 1)
  
  all_top_genes <- data.frame()
  protein_dir <- file.path(plot_dir, protein_name)
  dir.create(protein_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (pair in sample_pairs) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    if (!(s1 %in% colnames(cor_df)) || !(s2 %in% colnames(cor_df))) next
    
    df_plot <- data.frame(
      Gene = rownames(cor_df),
      X = cor_df[[s1]],
      Y = cor_df[[s2]]
    )
    
    # Select the top 20 genes with the highest Y-axis correlation
    top_y <- df_plot[order(-df_plot$Y), ][1:min(20, nrow(df_plot)), ]
    top_y$Sample1 <- s1
    top_y$Sample2 <- s2
    top_y$Protein <- protein_name
    all_top_genes <- rbind(all_top_genes, top_y)
    
    # Add protein column for color mapping
    df_plot$Protein <- protein_name
    
    # Plot scatter plot
    p <- ggplot(df_plot, aes(x = X, y = Y, color = Protein)) +
      geom_point(size = 1.2, alpha = 0.8) +
      geom_text_repel(
        data = top_y,
        aes(label = Gene),
        color = "black",
        size = 5.5,
        max.overlaps = 100,
        show.legend = FALSE
      ) +
      scale_color_manual(values = protein_colors) +
      ggtitle(paste0(protein_name)) +
      xlab(paste0(s1, " correlation")) +
      ylab(paste0(s2, " correlation")) +
      xlim(c(-0.06, NA)) +
      ylim(c(-0.02, NA)) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 18),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
      )
    
    ggsave(
      filename = file.path(protein_dir, paste0(protein_name, "_", s1, "_vs_", s2, "_TopY20.pdf")),
      plot = p,
      width = 6,
      height = 6.5
    )
  }
  
  # Save integrated top gene table
  write.csv(
    all_top_genes,
    file.path(protein_dir, paste0(protein_name, "_TopPos20_AllPairs.csv")),
    row.names = FALSE
  )
}


######## Part 2: Read Top20 Files, Extract Expression Data, and Generate Functional Analysis ########

library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# Set base directory and load RNA data
base_dir <- "./rna"
ScRNA <- readRDS(file.path(base_dir, "celltype.rds"))
ScRNA$orig.ident <- as.factor(ScRNA@meta.data$orig.ident)
treatment_groups <- levels(ScRNA$orig.ident)

# Set output directory for functional analysis results
analysis_dir <- file.path(base_dir, "out(T)", "by_group", "functional_analysis")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

# Define target protein names
target_proteins <- c("AKR1B10", "CD276", "FGL1", "IL4I1")

# Get target protein directories
all_protein_dirs <- list.files(
  file.path(base_dir, "out(T)", "by_group", "plots_40_patient"),
  full.names = TRUE
)

target_dirs <- all_protein_dirs[basename(all_protein_dirs) %in% target_proteins]

# Merge all Top20 genes
all_genes <- c()
gene_source <- list()

for (protein_dir in target_dirs) {
  protein_name <- basename(protein_dir)
  csv_file <- file.path(protein_dir, paste0(protein_name, "_TopPos20_AllPairs.csv"))
  
  if (!file.exists(csv_file)) next
  
  gene_info <- read.csv(csv_file)
  gene_list <- unique(gene_info$Gene)
  gene_list <- gene_list[gene_list %in% rownames(ScRNA)]
  
  if (length(gene_list) == 0) next
  
  all_genes <- c(all_genes, gene_list)
  gene_source[[protein_name]] <- gene_list
}

# Merge and deduplicate genes
all_genes <- unique(all_genes)

if (length(all_genes) == 0) {
  stop("No valid genes were found for analysis.")
}

# Extract expression matrix
expr_matrix <- GetAssayData(ScRNA, slot = "data")[all_genes, ]
avg_expr <- AverageExpression(ScRNA, features = all_genes, group.by = "orig.ident")$RNA
avg_expr_selected <- avg_expr[, treatment_groups]

# Save expression data
write.table(
  avg_expr_selected,
  file = file.path(analysis_dir, "Top5Proteins_Merged_AvgExpr.txt"),
  sep = "\t",
  quote = FALSE
)

# Generate DotPlot
genes_to_plot <- all_genes[all_genes %in% rownames(ScRNA)]

# Calculate average expression across treatment groups
avg_exp <- Seurat::AverageExpression(
  ScRNA,
  features = genes_to_plot,
  group.by = "treatment",
  assays = "RNA"
)$RNA

# Calculate normalized expression of the Tumor group relative to other groups
other_groups <- setdiff(colnames(avg_exp), "Tumor")
tumor_normalized <- avg_exp[, "Tumor"] - rowMeans(avg_exp[, other_groups, drop = FALSE])

# Sort genes by normalized expression
genes_positive <- names(sort(tumor_normalized[tumor_normalized > 0], decreasing = TRUE))
genes_negative <- names(sort(tumor_normalized[tumor_normalized <= 0], decreasing = FALSE))

genes_ordered <- c(genes_positive, genes_negative)
genes_to_plot <- intersect(genes_ordered, genes_to_plot)

# Create PDF DotPlot
pdf(file.path(analysis_dir, "Top5Proteins_Grouped_DotPlot.pdf"), width = 6, height = 10)

DotPlot(
  ScRNA,
  features = genes_to_plot,
  group.by = "treatment"
) +
  scale_color_gradientn(colors = c("#47CFD1", "white", "#FF3366")) +
  RotatedAxis() +
  coord_flip() +
  labs(
    x = NULL,
    y = NULL,
    color = "Avg. Expression",
    size = "Pct. Expressed"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 18),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    panel.background = element_blank()
  )

dev.off()


######## Three-Set Venn Analysis: Upregulated Top20 vs Downregulated Top20 vs All-Patient TopPos20 ##########

suppressPackageStartupMessages({
  library(ggvenn)
  library(data.table)
  library(tidyverse)
})

# Set paths
analysis_dir <- file.path(base_dir, "out(T)", "by_group", "functional_analysis")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)

# Differential analysis output directory
output_dir <- file.path(analysis_dir, "differential_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Upregulated and downregulated gene files
up_file <- file.path(output_dir, "upregulated_genes_Tumor_vs_Normal.csv")
down_file <- file.path(output_dir, "downregulated_genes_Tumor_vs_Normal.csv")

# Root directory containing protein-specific TopPos20 files
patients_root <- file.path(base_dir, "out(T)", "by_group", "plots_40_patient")

# Helper function for safely reading the first column
safe_read_first_col <- function(filepath, n_top = NULL) {
  if (!file.exists(filepath)) {
    warning(sprintf("[WARN] File does not exist: %s", filepath))
    return(character(0))
  }
  
  dt <- tryCatch(
    suppressWarnings(
      fread(filepath, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
    ),
    error = function(e) {
      warning(sprintf("[WARN] Failed to read file: %s; %s", filepath, e$message))
      return(NULL)
    }
  )
  
  if (is.null(dt) || ncol(dt) < 1) return(character(0))
  
  vec <- dt[[1]]
  vec <- vec[!is.na(vec) & nzchar(trimws(vec))]
  
  if (!is.null(n_top)) {
    vec <- head(vec, n_top)
  }
  
  # Remove duplicated values while preserving the original order
  vec[!duplicated(vec)]
}

# Read upregulated and downregulated genes and extract the top 20 genes from the first column
up_top20 <- safe_read_first_col(up_file, n_top = 20)
down_top20 <- safe_read_first_col(down_file, n_top = 20)

# Export upregulated and downregulated Top20 gene lists
write.table(
  up_top20,
  file = file.path(output_dir, "Upregulated_Top20_first_col.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  down_top20,
  file = file.path(output_dir, "Downregulated_Top20_first_col.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Summarize first-column genes from all *_TopPos20_AllPairs.csv files
if (!dir.exists(patients_root)) {
  stop(sprintf("Directory does not exist: %s", patients_root))
}

subdirs <- list.dirs(patients_root, full.names = TRUE, recursive = FALSE)

all_patient_files <- unlist(
  lapply(subdirs, function(sd) {
    list.files(
      sd,
      pattern = "_TopPos20_AllPairs\\.csv$",
      full.names = TRUE,
      recursive = FALSE
    )
  }),
  use.names = FALSE
)

if (length(all_patient_files) == 0) {
  warning("No *_TopPos20_AllPairs.csv files were found in any subdirectory.")
}

patient_genes_list <- lapply(all_patient_files, safe_read_first_col)
patient_genes_vec <- unique(unlist(patient_genes_list, use.names = FALSE))
patient_genes_vec <- patient_genes_vec[!is.na(patient_genes_vec) & nzchar(trimws(patient_genes_vec))]

# Export merged all-patient TopPos20 gene list
write.table(
  patient_genes_vec,
  file = file.path(output_dir, "AllPatients_TopPos20_first_col_merged_unique.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Draw three-set Venn diagram
venn_sets <- list(
  "Top20 highly consensus correlated" = patient_genes_vec,
  "Top20 Upregulated" = up_top20,
  "Top20 Downregulated" = down_top20
)

venn_plot <- ggvenn(
  venn_sets,
  fill_color = c("#FFCCCC", "#A4CDE1", "#C8E6C9"),
  fill_alpha = 0.7,
  stroke_linetype = "longdash",
  set_name_size = 7,
  text_size = 7,
  show_percentage = FALSE
) +
  ggplot2::ggtitle("") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    text = element_text(size = 14)
  )

# Save Venn diagram
ggplot2::ggsave(
  file.path(output_dir, "Venn_3sets_TopPos20_vs_UpTop20_vs_DownTop20.pdf"),
  plot = venn_plot,
  width = 8,
  height = 6
)

ggplot2::ggsave(
  file.path(output_dir, "Venn_3sets_TopPos20_vs_UpTop20_vs_DownTop20.svg"),
  plot = venn_plot,
  width = 8,
  height = 6
)


######## Define Venn Regions ########

A_name <- "Top20 highly consensus correlated"
B_name <- "Top20 Upregulated"
C_name <- "Top20 Downregulated"

A <- unique(patient_genes_vec)
B <- unique(up_top20)
C <- unique(down_top20)

# Unique regions
A_only <- setdiff(A, union(B, C))
B_only <- setdiff(B, union(A, C))
C_only <- setdiff(C, union(A, B))

# Pairwise intersection excluding the third set
AB_only <- setdiff(intersect(A, B), C)

# Create gene sets
gene_sets <- list(
  "A_only" = A_only,
  "B_only" = B_only,
  "C_only" = C_only,
  "AB_only" = AB_only
)

# Create display name mapping
set_names <- c(
  "A_only" = paste0(A_name, " only"),
  "B_only" = paste0(B_name, " only"),
  "C_only" = paste0(C_name, " only"),
  "AB_only" = paste0(A_name, " ∩ ", B_name, " only")
)


######## Draw Gene Expression DotPlots for Selected Venn Regions ########

selected_regions <- c("A_only", "B_only", "C_only", "AB_only")

for (region_key in selected_regions) {
  genes_in_region <- gene_sets[[region_key]]
  region_display_name <- set_names[region_key]
  
  # Keep genes present in the ScRNA object
  genes_to_plot <- genes_in_region[genes_in_region %in% rownames(ScRNA)]
  
  if (length(genes_to_plot) > 0) {
    cat(sprintf(
      "Drawing DotPlot for %s region with %d genes\n",
      region_display_name,
      length(genes_to_plot)
    ))
    
    # Calculate average gene expression
    avg_expression <- AverageExpression(
      ScRNA,
      features = genes_to_plot,
      group.by = "treatment",
      assays = "RNA"
    )
    
    expr_matrix <- avg_expression$RNA
    
    # Apply different sorting strategies according to region type
    if (region_key == "A_only") {
      if ("Tumor" %in% colnames(expr_matrix) && sum(colnames(expr_matrix) != "Tumor") > 0) {
        other_groups <- setdiff(colnames(expr_matrix), "Tumor")
        tumor_normalized <- expr_matrix[, "Tumor"] -
          rowMeans(expr_matrix[, other_groups, drop = FALSE])
        
        genes_positive <- names(sort(tumor_normalized[tumor_normalized > 0], decreasing = TRUE))
        genes_negative <- names(sort(tumor_normalized[tumor_normalized <= 0], decreasing = FALSE))
        
        genes_sorted <- c(genes_positive, genes_negative)
        genes_sorted <- intersect(genes_sorted, genes_to_plot)
        
        cat("A_only region is sorted by Tumor-vs-other differential expression\n")
      } else {
        gene_means <- rowMeans(expr_matrix, na.rm = TRUE)
        genes_sorted <- names(sort(gene_means, decreasing = TRUE))
        cat("A_only region is sorted by mean expression because Tumor group is missing\n")
      }
    } else {
      gene_means <- rowMeans(expr_matrix, na.rm = TRUE)
      genes_sorted <- names(sort(gene_means, decreasing = TRUE))
      cat(sprintf("%s region is sorted by mean expression\n", region_key))
    }
    
    pdf_file <- file.path(output_dir, sprintf("DotPlot_%s.pdf", region_key))
    
    pdf(pdf_file, width = 7, height = max(6, length(genes_to_plot) * 0.3))
    
    p <- DotPlot(
      ScRNA,
      features = genes_sorted,
      group.by = "treatment"
    ) +
      scale_color_gradientn(colors = c("#47CFD1", "white", "#FF3366")) +
      RotatedAxis() +
      coord_flip() +
      labs(
        title = region_display_name,
        x = NULL,
        y = NULL,
        color = "Avg. Expression",
        size = "Pct. Expressed"
      ) +
      theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 20),
        axis.text.y = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16)
      )
    
    print(p)
    dev.off()
    
    # Save sorted gene list and sorting information
    if (region_key == "A_only" && exists("tumor_normalized")) {
      sorted_genes_df <- data.frame(
        gene = genes_sorted,
        tumor_normalized_expression = tumor_normalized[genes_sorted],
        avg_expression = rowMeans(expr_matrix[genes_sorted, , drop = FALSE]),
        stringsAsFactors = FALSE
      )
    } else {
      gene_means <- rowMeans(expr_matrix, na.rm = TRUE)
      sorted_genes_df <- data.frame(
        gene = genes_sorted,
        avg_expression = gene_means[genes_sorted],
        stringsAsFactors = FALSE
      )
    }
    
    sorted_file <- file.path(output_dir, sprintf("Sorted_Genes_%s.tsv", region_key))
    fwrite(sorted_genes_df, sorted_file, sep = "\t", quote = FALSE)
    
    cat(sprintf("Sorted gene list has been saved to: %s\n", sorted_file))
  } else {
    cat(sprintf(
      "Warning: no genes in %s region are present in the ScRNA object. Skipping plot.\n",
      region_display_name
    ))
  }
}


######## Combined DotPlot for Selected Venn Regions ########

selected_regions <- c("B_only", "AB_only", "A_only", "C_only")

all_selected_genes <- unique(unlist(gene_sets[selected_regions]))
genes_to_plot_combined <- all_selected_genes[all_selected_genes %in% rownames(ScRNA)]

if (length(genes_to_plot_combined) > 0 && length(genes_to_plot_combined) <= 100) {
  
  avg_expression_combined <- AverageExpression(
    ScRNA,
    features = genes_to_plot_combined,
    group.by = "treatment",
    assays = "RNA"
  )
  
  expr_matrix_combined <- avg_expression_combined$RNA
  
  can_tumor_diff <- "Tumor" %in% colnames(expr_matrix_combined) &&
    sum(colnames(expr_matrix_combined) != "Tumor") > 0
  
  if (can_tumor_diff) {
    other_groups_combined <- setdiff(colnames(expr_matrix_combined), "Tumor")
    tumor_normalized_all <- expr_matrix_combined[, "Tumor"] -
      rowMeans(expr_matrix_combined[, other_groups_combined, drop = FALSE])
  }
  
  region_sorted_lists <- list()
  
  for (rk in selected_regions) {
    genes_r <- intersect(gene_sets[[rk]], rownames(expr_matrix_combined))
    
    if (length(genes_r) == 0) {
      region_sorted_lists[[rk]] <- character(0)
      next
    }
    
    sub_mat <- expr_matrix_combined[genes_r, , drop = FALSE]
    
    if (rk == "A_only") {
      if (isTRUE(can_tumor_diff)) {
        tn <- tumor_normalized_all[genes_r]
        pos <- names(sort(tn[tn > 0], decreasing = TRUE))
        neg <- names(sort(tn[tn <= 0], decreasing = FALSE))
        genes_sorted_r <- c(pos, neg)
      } else {
        gene_means_r <- rowMeans(sub_mat, na.rm = TRUE)
        genes_sorted_r <- names(sort(gene_means_r, decreasing = TRUE))
      }
    } else {
      gene_means_r <- rowMeans(sub_mat, na.rm = TRUE)
      genes_sorted_r <- names(sort(gene_means_r, decreasing = TRUE))
    }
    
    region_sorted_lists[[rk]] <- genes_sorted_r
  }
  
  genes_sorted_combined <- unique(unlist(region_sorted_lists))
  genes_sorted_combined <- intersect(genes_sorted_combined, genes_to_plot_combined)
  
  pdf(
    file.path(output_dir, "DotPlot_All_Selected_Regions.pdf"),
    width = 7,
    height = max(8, length(genes_sorted_combined) * 0.3)
  )
  
  p_combined <- DotPlot(
    ScRNA,
    features = genes_sorted_combined,
    group.by = "treatment"
  ) +
    scale_color_gradientn(colors = c("#47CFD1", "white", "#FF3366")) +
    RotatedAxis() +
    coord_flip() +
    labs(
      title = "All Selected Regions",
      x = NULL,
      y = NULL,
      color = "Avg. Expression",
      size = "Pct. Expressed"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, size = 20),
      axis.text.y = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 18),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18)
    )
  
  print(p_combined)
  dev.off()
  
  region_of_gene <- vapply(
    genes_sorted_combined,
    function(g) {
      regs <- names(which(sapply(gene_sets[selected_regions], function(x) g %in% x)))
      paste(regs, collapse = ", ")
    },
    character(1)
  )
  
  avg_expr_vec <- rowMeans(expr_matrix_combined[genes_sorted_combined, , drop = FALSE], na.rm = TRUE)
  
  if (isTRUE(can_tumor_diff)) {
    tn_vec <- tumor_normalized_all[genes_sorted_combined]
    
    sorted_combined_df <- data.frame(
      gene = genes_sorted_combined,
      tumor_normalized_expression = tn_vec,
      avg_expression = avg_expr_vec,
      region = region_of_gene,
      stringsAsFactors = FALSE
    )
  } else {
    sorted_combined_df <- data.frame(
      gene = genes_sorted_combined,
      avg_expression = avg_expr_vec,
      region = region_of_gene,
      stringsAsFactors = FALSE
    )
  }
  
  fwrite(
    sorted_combined_df,
    file.path(output_dir, "Sorted_Genes_All_Regions.tsv"),
    sep = "\t",
    quote = FALSE
  )
  
} else if (length(genes_to_plot_combined) > 100) {
  cat(sprintf(
    "Note: the combined plot contains %d genes, which is too many. The combined plot was skipped.\n",
    length(genes_to_plot_combined)
  ))
}


######## Enrichment Analysis for Each Venn Region ########

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(data.table)
library(dplyr)

# Function for enrichment analysis
perform_enrichment_analysis <- function(gene_symbols, set_name, output_dir) {
  
  # Create subdirectory
  analysis_dir <- file.path(output_dir, gsub("[^[:alnum:]]", "_", set_name))
  
  if (!dir.exists(analysis_dir)) {
    dir.create(analysis_dir, recursive = TRUE)
  }
  
  # Check whether the gene list is valid
  if (length(gene_symbols) == 0) {
    message(sprintf("Gene set %s is empty. Skipping analysis.", set_name))
    return(NULL)
  }
  
  # Gene annotation
  gene_anno <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  if (nrow(gene_anno) == 0) {
    message(sprintf("Failed to annotate SYMBOL to ENTREZID for gene set %s", set_name))
    return(NULL)
  }
  
  results <- list()
  
  # GSEA KEGG analysis
  if (exists("avg_expr_selected") && length(gene_symbols) > 10) {
    tryCatch({
      common_genes <- intersect(gene_symbols, rownames(avg_expr_selected))
      
      if (length(common_genes) > 5) {
        logfc_vec <- rowMeans(avg_expr_selected[common_genes, , drop = FALSE])
        geneList <- logfc_vec
        names(geneList) <- gene_anno$ENTREZID[match(common_genes, gene_anno$SYMBOL)]
        geneList <- na.omit(geneList)
        geneList <- sort(geneList, decreasing = TRUE)
        
        if (length(geneList) > 10) {
          kk_gse <- gseKEGG(
            geneList = geneList,
            organism = "hsa",
            nPerm = 1000,
            minGSSize = 10,
            pvalueCutoff = 0.25,
            verbose = FALSE
          )
          
          if (nrow(kk_gse) > 0) {
            kk_gse <- setReadable(kk_gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
            
            write.csv(
              as.data.frame(kk_gse),
              file = file.path(analysis_dir, "GSEA_KEGG_results.csv")
            )
            
            pdf(file.path(analysis_dir, "GSEA_ridgeplot.pdf"), width = 10, height = 8)
            print(
              ridgeplot(kk_gse, showCategory = 15, fill = "pvalue") +
                theme(axis.text = element_text(size = 12)) +
                ggtitle(paste("GSEA KEGG -", set_name))
            )
            dev.off()
            
            if (nrow(kk_gse) >= 3) {
              pdf(file.path(analysis_dir, "GSEA_dotplot.pdf"), width = 10, height = 8)
              print(
                dotplot(kk_gse, showCategory = 15) +
                  ggtitle(paste("GSEA KEGG Dotplot -", set_name))
              )
              dev.off()
            }
            
            results[["GSEA_KEGG"]] <- kk_gse
          }
        }
      }
    }, error = function(e) {
      message(sprintf("GSEA KEGG analysis failed (%s): %s", set_name, e$message))
    })
  }
  
  # GO enrichment analysis
  tryCatch({
    go_bp <- enrichGO(
      gene = gene_anno$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
    
    if (nrow(go_bp) > 0) {
      go_bp <- setReadable(go_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      write.csv(
        as.data.frame(go_bp),
        file = file.path(analysis_dir, "GO_BP_results.csv")
      )
      
      pdf(file.path(analysis_dir, "GO_dotplot.pdf"), width = 10, height = 8)
      print(
        dotplot(go_bp, showCategory = 20) +
          ggtitle(paste("GO Biological Process -", set_name))
      )
      dev.off()
      
      pdf(file.path(analysis_dir, "GO_barplot.pdf"), width = 10, height = 8)
      print(
        barplot(go_bp, showCategory = 20) +
          ggtitle(paste("GO Biological Process -", set_name))
      )
      dev.off()
      
      results[["GO_BP"]] <- go_bp
    }
  }, error = function(e) {
    message(sprintf("GO enrichment analysis failed (%s): %s", set_name, e$message))
  })
  
  # KEGG pathway enrichment analysis
  tryCatch({
    kegg_enrich <- enrichKEGG(
      gene = gene_anno$ENTREZID,
      organism = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1,
      qvalueCutoff = 0.2
    )
    
    if (nrow(kegg_enrich) > 0) {
      kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      write.csv(
        as.data.frame(kegg_enrich),
        file = file.path(analysis_dir, "KEGG_pathway_results.csv")
      )
      
      pdf(file.path(analysis_dir, "KEGG_dotplot.pdf"), width = 10, height = 8)
      print(
        dotplot(kegg_enrich, showCategory = 20) +
          ggtitle(paste("KEGG Pathway -", set_name))
      )
      dev.off()
      
      pdf(file.path(analysis_dir, "KEGG_barplot.pdf"), width = 10, height = 8)
      print(
        barplot(kegg_enrich, showCategory = 20) +
          ggtitle(paste("KEGG Pathway -", set_name))
      )
      dev.off()
      
      if (nrow(kegg_enrich) >= 5 && nrow(kegg_enrich) <= 30) {
        pdf(file.path(analysis_dir, "KEGG_emapplot.pdf"), width = 12, height = 10)
        print(
          emapplot(pairwise_termsim(kegg_enrich)) +
            ggtitle(paste("KEGG Pathway Network -", set_name))
        )
        dev.off()
      }
      
      results[["KEGG"]] <- kegg_enrich
    }
  }, error = function(e) {
    message(sprintf("KEGG enrichment analysis failed (%s): %s", set_name, e$message))
  })
  
  # Save gene list
  write.csv(
    data.frame(
      Gene_Symbol = gene_symbols,
      ENTREZID = gene_anno$ENTREZID[match(gene_symbols, gene_anno$SYMBOL)]
    ),
    file = file.path(analysis_dir, "gene_list.csv")
  )
  
  return(results)
}

# Run enrichment analysis for all gene sets
all_enrichment_results <- list()

for (set_id in names(gene_sets)) {
  set_name <- set_names[set_id]
  genes <- gene_sets[[set_id]]
  
  message(sprintf("Analyzing gene set: %s (%d genes)", set_name, length(genes)))
  
  results <- perform_enrichment_analysis(genes, set_name, output_dir)
  
  if (!is.null(results)) {
    all_enrichment_results[[set_id]] <- results
    message(sprintf("Analysis completed: %s", set_name))
  } else {
    message(sprintf("Analysis skipped: %s", set_name))
  }
}

# Save analysis summary
summary_data <- data.frame(
  GeneSet = names(set_names),
  SetName = set_names,
  NumGenes = sapply(gene_sets, length),
  AnalysisCompleted = names(gene_sets) %in% names(all_enrichment_results)
)

write.csv(summary_data, file.path(output_dir, "enrichment_analysis_summary.csv"))

# Generate comparison plots if multiple gene sets were successfully analyzed
successful_sets <- names(all_enrichment_results)

if (length(successful_sets) > 1) {
  tryCatch({
    go_list <- list()
    
    for (set_id in successful_sets) {
      if ("GO_BP" %in% names(all_enrichment_results[[set_id]])) {
        go_list[[set_names[set_id]]] <- all_enrichment_results[[set_id]][["GO_BP"]]
      }
    }
    
    if (length(go_list) > 1) {
      pdf(file.path(output_dir, "GO_comparison_dotplot.pdf"), width = 12, height = 10)
      print(
        dotplot(go_list, showCategory = 10) +
          ggtitle("GO Biological Process Comparison")
      )
      dev.off()
    }
  }, error = function(e) {
    message("An error occurred while generating comparison plots: ", e$message)
  })
}

message("Enrichment analysis completed. Results were saved to: ", output_dir)


######## Overall Enrichment Analysis for All Candidate Genes ########

# Gene annotation
gene_anno <- bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

if (nrow(gene_anno) == 0) {
  stop("Failed to annotate gene SYMBOL to ENTREZID.")
}

# Construct geneList using average expression values
logfc_vec <- rowMeans(avg_expr_selected[gene_anno$SYMBOL, , drop = FALSE])
geneList <- logfc_vec
names(geneList) <- gene_anno$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

# GSEA KEGG analysis
kk_gse <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  nPerm = 1000,
  minGSSize = 10,
  pvalueCutoff = 0.25,
  verbose = FALSE
)

if (nrow(kk_gse) > 0) {
  kk_gse <- setReadable(kk_gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  write.csv(
    as.data.frame(kk_gse),
    file = file.path(analysis_dir, "GSEA_KEGG_results.csv")
  )
  
  pdf(file.path(analysis_dir, "GSEA_ridgeplot.pdf"), width = 10, height = 8)
  print(
    ridgeplot(kk_gse, showCategory = 15, fill = "pvalue") +
      theme(axis.text = element_text(size = 12))
  )
  dev.off()
}

# GO enrichment analysis
go_res <- enrichGO(
  gene = gene_anno$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

if (nrow(go_res) > 0) {
  go_res <- setReadable(go_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  write.csv(
    as.data.frame(go_res),
    file = file.path(analysis_dir, "GO_BP_results.csv")
  )
  
  pdf(file.path(analysis_dir, "GO_dotplot.pdf"), width = 8, height = 6)
  print(
    dotplot(go_res, showCategory = 20) +
      ggtitle("GO Biological Process Enrichment")
  )
  dev.off()
}

# KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(
  gene = gene_anno$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.1
)

if (nrow(kegg_enrich) > 0) {
  kegg_enrich <- setReadable(kegg_enrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  write.csv(
    as.data.frame(kegg_enrich),
    file = file.path(analysis_dir, "KEGG_pathway_results.csv")
  )
  
  pdf(file.path(analysis_dir, "KEGG_dotplot.pdf"), width = 8, height = 6)
  print(
    dotplot(kegg_enrich, showCategory = 20) +
      ggtitle("KEGG Pathway Enrichment")
  )
  dev.off()
}


