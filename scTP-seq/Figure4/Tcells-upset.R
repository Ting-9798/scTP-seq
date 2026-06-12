########################## Combined Plotting #######################
######### Expression Trend of Marker Genes Across Cells ###########

col <- c(
  '#FF6666', '#E5D2DD', "#BC8F8F", '#FFCC99', '#4F6272', '#58A4C3', "#CC0066",
  '#F3B1A0', '#57C3F3', '#E59CC4', '#437eb8', "#FF3366", '#FF9999',
  "#66CCCC", "#99CCFF", '#3399CC',
  "#FF9933", "#CCFFCC", "#00CC66", "#99FFFF", "#FF3300", '#F9BB72',
  "#6699CC", "#9999FF", "#CCCCFF", "#CC99CC", "#FF6699", "#6699CC", "#FFFFCC"
)

# Set working directory and output directory
setwd("./out/T_cells/")
outdir <- "./out/T_cells/"

# install.packages("UpSetR")
library(UpSetR)
library(RColorBrewer)
library(ggplot2)

# Import the gene grouping matrix file
asvs <- read.csv("network_venn_gene_group_matrix.csv")

# Check plotting data
head(asvs)
colnames(asvs)

# Convert plotting data to binary 0-1 format if needed
# asvs[asvs > 0] <- 1

# Basic UpSetR plot
upset(asvs)

# Save the UpSet plot as PDF
pdf(file = paste0(outdir, "upset_plot.pdf"), width = 8, height = 4.5)

upset(
  asvs,
  nset = 6,
  nintersects = 20,
  order.by = c("degree", "freq"),
  decreasing = c(TRUE, TRUE),
  sets.bar.color = c("#9BD4D3", "#FBDABE", "#C9E6D7"),
  set_size.show = FALSE,
  show.numbers = "yes",
  sets.x.label = "Gene Count",
  mainbar.y.label = "Intersection Size",
  number.angles = 0,
  set_size.angles = 0,
  main.bar.color = "#57C3F3",
  matrix.color = "#9933FF",
  point.size = 2.8,
  line.size = 1.2,
  mb.ratio = c(0.7, 0.3),
  shade.color = "gray86",
  shade.alpha = 0.5,
  matrix.dot.alpha = 0.5,
  text.scale = c(2, 2, 2, 2, 2, 2),
  empty.intersections = "on",  # Show empty intersections
  queries = list(
    list(
      query = intersects,
      params = list("Highly.consensus.correlated", "Top20.Upregulated"),
      active = TRUE,
      color = "#DC050C"
    )
  )
)

dev.off()


