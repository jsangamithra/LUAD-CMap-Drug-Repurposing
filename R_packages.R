# Run this once to install all required R packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install(c(
  "DESeq2", "TCGAbiolinks", "EnhancedVolcano",
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA",
  "ComplexHeatmap", "limma", "edgeR"
), ask = FALSE)

install.packages(c(
  "tidyverse", "ggplot2", "ggrepel", "pheatmap",
  "RColorBrewer", "patchwork", "viridis", "ggpubr",
  "writexl", "data.table"
))

cat("All R packages installed successfully!\n")
