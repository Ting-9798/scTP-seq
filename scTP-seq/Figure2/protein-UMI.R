####################### Proteomics - UMI #####################
############### Merge Replicates and Plot ####################
######## Functional Proteins ########

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
library(tidyr)
# install.packages("tidydr")

# Set working directory
setwd("./data/pro/")

# Define output directory
outdir <- "./out(pro)/umi/"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
OUTPUT <- paste(outdir)

# Define protein barcode dictionary
protein_dict <- list(
  "barcode-1" = "SCCA",
  "barcode-2" = "CEA",
  "barcode-3" = "NSE",
  "barcode-4" = "Cyfra 21-1",
  "barcode-5" = "ProGRP"
)

# Define input folders
folders <- c("A")

# Initialize merged data frame
merged_data <- data.frame()

# Iterate through each folder
for (folder in folders) {
  
  # Construct file paths
  barcode_file <- paste0(folder, "/", "barcodes.tsv")
  feature_file <- paste0(folder, "/", "features.tsv")
  if (!file.exists(feature_file)) feature_file <- file.path(folder, "/", "genes.tsv")
  matrix_file <- paste0(folder, "/", "matrix.mtx")
  
  # Read single-cell expression matrix data
  barcodes <- read.delim(barcode_file, header = FALSE, stringsAsFactors = FALSE)
  features <- read.delim(feature_file, header = FALSE, stringsAsFactors = FALSE)
  matrix <- readMM(matrix_file)
  
  # Assign row names and column names to the matrix
  rownames(matrix) <- features$V1
  colnames(matrix) <- barcodes$V1
  
  # Extract PB and UMI data
  data <- as.data.frame(as.matrix(matrix))
  data$barcode <- rownames(data)
  
  # Convert data into long format
  data <- data %>%
    gather(key = "CB", value = "UMI", -barcode) %>%
    rename(PB = barcode)
  
  # Annotate PB based on protein barcode dictionary
  data$Pro <- sapply(data$PB, function(pb) {
    if (pb %in% names(protein_dict)) {
      return(protein_dict[[pb]])
    } else {
      return("Others")
    }
  })
  
  # Remove rows annotated as "Others"
  data <- data %>% filter(Pro != "Others")
  
  # Merge into the final data frame
  merged_data <- rbind(merged_data, data)
}

# Set protein order as an ordered factor
# merged_data$Pro <- factor(
#   merged_data$Pro,
#   levels = c("SCCA", "CEA", "NSE", "Cyfra 21-1", "ProGRP")
# )

merged_data_matrix <- merged_data

merged_data_matrix$Pro <- factor(
  merged_data_matrix$Pro,
  levels = c("SCCA", "CEA", "NSE", "Cyfra 21-1", "ProGRP")
)

# Convert UMI counts to numeric values
merged_data_matrix$UMI <- as.numeric(merged_data_matrix$UMI)

merged_data <- merged_data_matrix

# Define color palette
col <- c(
  '#57C3F3', "#FF3366", "#66CCCC", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
  "#A4CDE1", '#FF9999', "#66CCCC", '#4F6272', "#FF3366", "#CC0066", "#CC99CC",
  "#FFCCCC", "#CCFFCC", "#FFFFCC", '#E5D2DD', '#58A4C3',
  '#F9BB72', '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8',
  "#66CCCC", "#99CCFF", '#3399CC', "#FF3366", "#CC0066",
  "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300",
  "#6699CC", "#9999FF", "#CCCCFF", "#FF6699", "#6699CC", "#FFFFCC"
)

##################### Plot Total Protein UMI for Each Sample in Each Group ###################

# Calculate total UMI counts for each cell barcode within each group and sample
group_total <- merged_data %>%
  group_by(group, CB, sample) %>%
  summarise(TotalUMI = sum(UMI), .groups = "drop")

# Define group order
group_order <- c("NTC", "3T3", "A549")

# Set quantile filtering thresholds
quantiles <- c(0.95, 0.999)
quantile_labels <- c("5", "0.1")

for (i in seq_along(quantiles)) {
  
  q <- quantiles[i]
  label <- quantile_labels[i]
  
  # Filter cell barcodes within each group based on the selected quantile
  filtered_total <- group_total %>%
    group_by(group) %>%
    filter(TotalUMI <= quantile(TotalUMI, q)) %>%
    ungroup()
  
  # Set sample order
  filtered_total$sample <- factor(
    filtered_total$sample,
    levels = unique(filtered_total$sample[filtered_total$group %in% group_order])
  )
  
  # Plot violin plot with boxplot overlay for each sample
  output_pdf <- paste0(outdir, "UMI_Total_ViolinPlot_group_", label, ".pdf")
  pdf(output_pdf, width = 6, height = 6)
  
  # Create base plot
  p <- ggplot(filtered_total, aes(x = sample, y = TotalUMI, fill = group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = col) +
    labs(
      title = "",
      x = "Sample",
      y = "Total UMI"
    ) +
    ylim(-100, 1000) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  
  # Calculate statistical significance among samples within each group
  groups <- unique(filtered_total$group)
  
  # Initialize data frame for significant comparisons
  sig_data <- data.frame()
  
  for (grp in groups) {
    
    # Extract data for the current group
    group_data <- filtered_total %>% filter(group == grp)
    
    # Get all samples in the current group
    samples <- unique(group_data$sample)
    
    # Perform statistical analysis only when there are at least three samples within the group
    if (length(samples) >= 3) {
      
      # Perform Kruskal-Wallis test
      kw_test <- kruskal.test(TotalUMI ~ sample, data = group_data)
      
      # If the overall test is significant, perform pairwise comparisons
      if (kw_test$p.value < 0.05) {
        
        # Perform Dunn's test for multiple comparisons
        if (requireNamespace("dunn.test", quietly = TRUE)) {
          
          dunn_result <- dunn.test::dunn.test(
            group_data$TotalUMI,
            group_data$sample,
            method = "bonferroni"
          )
          
          # Extract significant pairwise comparisons
          for (j in 1:length(dunn_result$comparisons)) {
            
            if (dunn_result$P.adjusted[j] < 0.05) {
              
              # Parse sample pairs
              comp <- unlist(strsplit(dunn_result$comparisons[j], " - "))
              
              # Add significant comparison result
              sig_data <- rbind(
                sig_data,
                data.frame(
                  group = grp,
                  sample1 = comp[1],
                  sample2 = comp[2],
                  p_value = dunn_result$P.adjusted[j]
                )
              )
            }
          }
          
        } else {
          
          # If dunn.test is not installed, use pairwise Wilcoxon test instead
          warning("Package 'dunn.test' is not installed. Pairwise Wilcoxon test will be used instead.")
          
          # Generate all pairwise sample comparisons
          sample_pairs <- combn(samples, 2, simplify = FALSE)
          
          for (pair in sample_pairs) {
            
            sample1_data <- group_data %>% filter(sample == pair[1])
            sample2_data <- group_data %>% filter(sample == pair[2])
            
            # Perform Wilcoxon test
            wilcox_test <- wilcox.test(sample1_data$TotalUMI, sample2_data$TotalUMI)
            
            # Apply Bonferroni correction
            adj_p <- p.adjust(
              wilcox_test$p.value,
              method = "bonferroni",
              n = length(sample_pairs)
            )
            
            if (adj_p < 0.05) {
              
              sig_data <- rbind(
                sig_data,
                data.frame(
                  group = grp,
                  sample1 = pair[1],
                  sample2 = pair[2],
                  p_value = adj_p
                )
              )
            }
          }
        }
      }
    }
  }
  
  # Add significance annotations to the plot if significant comparisons are detected
  if (nrow(sig_data) > 0) {
    
    # Install and load ggsignif for significance annotations
    if (!requireNamespace("ggsignif", quietly = TRUE)) {
      install.packages("ggsignif")
    }
    library(ggsignif)
    
    # Add significance annotations for each significant comparison
    for (k in 1:nrow(sig_data)) {
      
      # Determine significance level
      p_val <- sig_data$p_value[k]
      significance <- ""
      
      if (p_val < 0.001) {
        significance <- "***"
      } else if (p_val < 0.01) {
        significance <- "**"
      } else if (p_val < 0.05) {
        significance <- "*"
      }
      
      if (significance != "") {
        
        # Get x-axis positions of the compared samples
        x1 <- which(levels(filtered_total$sample) == sig_data$sample1[k])
        x2 <- which(levels(filtered_total$sample) == sig_data$sample2[k])
        
        # Calculate y-axis position for significance annotation
        base_height <- 1000
        
        # Start from 800 and decrease by 100 for each comparison to reduce overlap
        y_pos <- base_height - (k * 100)
        
        # Prevent the annotation from being placed too low
        if (y_pos < 850) {
          y_pos <- 850 + (k %% 3) * 50
        }
        
        # Add significance annotation
        p <- p + geom_signif(
          comparisons = list(c(sig_data$sample1[k], sig_data$sample2[k])),
          annotations = significance,
          y_position = y_pos,
          tip_length = 0.01,
          vjust = 0.5,
          textsize = 4
        )
      }
    }
  }
  
  print(p)
  dev.off()
  
  # Save significance results to a text file
  if (nrow(sig_data) > 0) {
    
    sig_output <- paste0(outdir, "UMI_Total_Significance_", label, ".txt")
    
    write.table(
      sig_data,
      sig_output,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    cat(paste("Significance results saved to:", sig_output, "\n"))
  }
}


