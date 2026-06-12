########################## Merged Plotting #######################
######### Expression Patterns of Marker Genes Across Cells #########

col <- c('#FF6666','#E5D2DD',"#BC8F8F",'#FFCC99','#4F6272','#58A4C3',"#CC0066",
         '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  "#FF3366",'#FF9999',
         "#66CCCC","#99CCFF", '#3399CC',
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300", '#F9BB72', 
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699","#6699CC","#FFFFCC")

setwd("./out")
outdir <- "./out/"


####################### Tumor Markers ##########################
######### Expression Patterns of Marker Genes Across Cells #########

col <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
         "#A4CDE1",'#FF9999',"#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",
         '#E5D2DD','#58A4C3',
         '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',
         "#99CCFF", '#3399CC',"#FF3366","#CC0066",
         "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
         "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699",
         "#6699CC","#FFFFCC")

# Select differentially expressed cell populations for visualization
output <- paste(outdir, "marker", sep = "/")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

# Read the annotated Seurat object
file_path <- file.path(outdir, "celltype.rds")
ScRNA <- readRDS(file_path)


# Define tumor marker genes
cellmarker <- c(
  "SERPINB3", "CEACAM5", "ENO2", "KRT19", "GRP"
)

# Retain only marker genes present in the dataset
cellmarker <- cellmarker[cellmarker %in% rownames(ScRNA)]


# Define sample colors
col_sample <- c('#FF9999',"#A4CDE1","#66CCCC","#FFCCCC","#CCFFCC","#FFFFCC",
                '#E5D2DD','#4F6272','#58A4C3',
                '#F9BB72', '#F3B1A0','#57C3F3', '#E59CC4','#437eb8',  
                "#66CCCC","#99CCFF", '#3399CC',"#FF3366","#CC0066",
                "#FF9933","#CCFFCC","#00CC66","#99FFFF","#FF3300",
                "#6699CC","#9999FF","#CCCCFF","#CC99CC","#FF6699",
                "#6699CC","#FFFFCC")


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


# Load required packages
library(ggpubr)
library(ggplot2)
library(dplyr)

# Create a list to store box-violin plots
boxviolin_plots_sample <- list()

for (gene in cellmarker) {
  
  # Extract expression data for the current gene
  gene_exp <- ScRNA@assays$RNA@data[gene, ]
  
  # Calculate the 1st and 99th percentiles
  q01 <- quantile(gene_exp, 0.01, na.rm = TRUE)
  q99 <- quantile(gene_exp, 0.99, na.rm = TRUE)
  
  # Remove the lowest 1% of cells based on expression
  filtered_cells <- which(gene_exp > q01)
  
  # Create a new Seurat object containing the filtered cells
  ScRNA_filtered <- subset(ScRNA, cells = names(filtered_cells))
  
  # Extract expression data and metadata for statistical testing
  exp_data <- as.numeric(ScRNA_filtered@assays$RNA@data[gene, ])
  treatment_groups <- as.character(ScRNA_filtered@meta.data$treatment)
  
  # Create a data frame for plotting and statistical testing
  plot_data <- data.frame(
    expression = exp_data,
    treatment = factor(treatment_groups)
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$expression), ]
  
  # Check whether there are enough groups and data points for statistical testing
  unique_treatments <- levels(plot_data$treatment)
  
  # Calculate y-axis range for adjusting p-value label positions
  y_max <- max(plot_data$expression, na.rm = TRUE)
  y_range <- y_max - min(plot_data$expression, na.rm = TRUE)
  
  # Dynamically adjust label height parameters
  step_increase_factor <- 0.15
  bracket_adjustment <- 0.08 * y_range
  label_y_adjustment <- 0.1 * y_range
  
  # Create box-violin plot
  p <- ggplot(plot_data, aes(x = treatment, y = expression, fill = treatment)) +
    
    # Violin plot with transparency
    geom_violin(alpha = 0.7, scale = "width", trim = TRUE) +
    
    # Boxplot overlay
    geom_boxplot(
      width = 0.2,
      alpha = 0.8,
      outlier.shape = NA,
      fill = "white",
      color = "black"
    ) +
    
    # Set fill colors
    scale_fill_manual(values = treatment_colors) +
    
    # Theme settings
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 20, margin = margin(r = 10)),
      axis.title.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12)
    ) +
    labs(
      title = gene,
      y = "Expression Level"
    )
  
  # Add statistical tests according to the number of treatment groups
  if (length(unique_treatments) >= 2 && nrow(plot_data) > 0) {
    
    # Check whether each group has enough data points
    group_counts <- table(plot_data$treatment)
    valid_groups <- names(group_counts)[group_counts >= 3]
    
    if (length(valid_groups) >= 2) {
      
      # Keep only valid groups
      plot_data_filtered <- plot_data %>%
        filter(treatment %in% valid_groups)
      
      # Calculate y-axis range after filtering
      y_max_filtered <- max(plot_data_filtered$expression, na.rm = TRUE)
      y_min_filtered <- min(plot_data_filtered$expression, na.rm = TRUE)
      y_range_filtered <- y_max_filtered - y_min_filtered
      
      # Choose statistical method according to the number of groups
      if (length(valid_groups) == 2) {
        
        # Two-group comparison: Wilcoxon rank-sum test
        test_result <- wilcox.test(expression ~ treatment, data = plot_data_filtered)
        p_val <- test_result$p.value
        
        # Calculate y-axis position for the significance label
        label_y_base <- y_max_filtered + 0.1 * y_range_filtered
        
        # Add significance annotation
        p <- p +
          stat_compare_means(
            method = "wilcox.test",
            comparisons = list(valid_groups),
            label = "p.signif",
            size = 6,
            label.y = label_y_base,
            tip.length = 0.01,
            bracket.size = 0.6,
            vjust = 0.5
          )
        
        # Generate exact p-value label
        p_label <- ifelse(
          p_val < 0.001, "p < 0.001",
          ifelse(
            p_val < 0.01, "p < 0.01",
            sprintf("p = %.3f", p_val)
          )
        )
        
      } else {
        
        # Multi-group comparison: Kruskal-Wallis test
        kruskal_result <- kruskal.test(expression ~ treatment, data = plot_data_filtered)
        p_val_kruskal <- kruskal_result$p.value
        
        # Add global test result as subtitle
        p_label_kruskal <- ifelse(
          p_val_kruskal < 0.001, "p < 0.001",
          ifelse(
            p_val_kruskal < 0.01, "p < 0.01",
            sprintf("p = %.3f", p_val_kruskal)
          )
        )
        
        p <- p + labs(subtitle = paste("Kruskal-Wallis test:", p_label_kruskal))
        
        # Perform pairwise comparisons if the global test is significant
        if (p_val_kruskal < 0.05) {
          
          # Generate all pairwise comparisons
          comparisons_list <- combn(valid_groups, 2, simplify = FALSE)
          
          # Calculate the number of comparisons
          n_comparisons <- length(comparisons_list)
          
          # Dynamically adjust step.increase according to the number of comparisons
          step_increase <- min(0.15, 0.6 / n_comparisons)
          
          # Calculate the base y-axis position for the first comparison
          base_y <- y_max_filtered + 0.15 * y_range_filtered
          
          # Add pairwise significance annotations
          p <- p +
            stat_compare_means(
              method = "wilcox.test",
              comparisons = comparisons_list,
              label = "p.signif",
              size = 6,
              tip.length = 0.01,
              bracket.size = 0.5,
              step.increase = step_increase,
              label.y.npc = "top",
              label.y = base_y,
              hide.ns = TRUE
            )
        }
      }
    }
  }
  
  # Adjust y-axis upper limit to leave space for significance labels
  current_ymax <- layer_scales(p)$y$range$range[2]
  new_ymax <- current_ymax * 1.1
  
  # Set y-axis limits
  p <- p + coord_cartesian(
    ylim = c(min(plot_data$expression, na.rm = TRUE), new_ymax)
  )
  
  # Add plot to the list
  boxviolin_plots_sample[[gene]] <- p
}

# Save plots as PDF
pdf(
  paste0(output, "/cellmarker_BoxViolinPlot_bySample_filtered_with_stats.pdf"),
  width = 16,
  height = 4
)
print(cowplot::plot_grid(plotlist = boxviolin_plots_sample, ncol = 5))
dev.off()

# Save plots as SVG
svg(
  paste0(output, "/cellmarker_BoxViolinPlot_bySample_filtered_with_stats.svg"),
  width = 16,
  height = 4
)
print(cowplot::plot_grid(plotlist = boxviolin_plots_sample, ncol = 5))
dev.off()




