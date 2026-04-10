#########################
# Version Numbers:
# R 4.5.1
# Seurat 5.4.0
# SeuratObject 5.3.0
# dplyr 1.2.1
# ggplot2 4.0.2
# patchwork 1.3.2
# pheatmap 1.0.13
# SingleR 2.12.0
# celldex 1.20.0
# SingleCellExperiment 1.32.0
# DESeq2 1.50.2
# clusterProfiler 4.18.4
# org.Mm.eg.db 3.22.0
# BiocParallel 1.44.0
# Bioconductor 3.22
#########################

# Install helpers
install.packages("remotes")
install.packages("BiocManager")

# Set Bioconductor version
BiocManager::install(version = "3.22", ask = FALSE)

# CRAN packages (exact versions) 
remotes::install_version("Seurat", version = "5.4.0")
remotes::install_version("SeuratObject", version = "5.3.0")
remotes::install_version("dplyr", version = "1.2.1")
remotes::install_version("ggplot2", version = "4.0.2")
remotes::install_version("patchwork", version = "1.3.2")
remotes::install_version("pheatmap", version = "1.0.13")

# Bioconductor packages 
BiocManager::install(c(
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "DESeq2",
  "clusterProfiler",
  "org.Mm.eg.db",
  "BiocParallel"
), ask = FALSE)

