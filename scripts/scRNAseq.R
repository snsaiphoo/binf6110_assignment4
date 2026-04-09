library(Seurat) 
library(dplyr)
library(patchwork)
library(ggplot2)

# load in Seurat object 
nasal <- readRDS("../seurat_ass4.rds")

nasal

set.seed(123)

table(nasal$orig.ident)

# Mitochondrial, Hemoglobin, and Ribosomal protein Percentage 
nasal[["percent.mt"]] <- PercentageFeatureSet(nasal, pattern = "^mt-")

nasal[["percent.hb"]] <- PercentageFeatureSet(nasal, pattern = "^Hb[ab]")

nasal[["percent.ribo"]] <- PercentageFeatureSet(nasal, pattern = "^Rp[sl]")


# Quality Control Visualizations

# Violin Plots
vln_plot_1 <- VlnPlot(
  nasal,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2
)

vln_plot_2 <- VlnPlot(
  nasal,
  features = c("percent.mt", "percent.ribo", "percent.hb"),
  ncol = 3
)

# ggsave(
#   filename = "../figures/QC_violin_plot_1.png",
#   plot = vln_plot_1,
#   width = 10,
#   height = 4,
#   dpi = 300
# )
# 
# ggsave(
#   filename = "../figures/QC_violin_plot_2.png",
#   plot = vln_plot_2,
#   width = 10,
#   height = 4,
#   dpi = 300
# )

# Knee Plot
nasal <- CalculateBarcodeInflections(nasal)

knee_plot <- BarcodeInflectionsPlot(nasal)

# ggsave(
#   filename = "../figures/QC_knee_plot.png",
#   plot = knee_plot,
#   width = 6,
#   height = 5,
#   dpi = 300
# )

# Based on QC plots:
# nFeature_RNA: remove low-quality cells (<500) and likely doublets (>4000)
# nCount_RNA: remove near-empty droplets (<1000) and outlier high-count cells (>10000)
# percent.mt: <15% threshold, slightly lenient for infection context
# percent.hb: <5% to remove rare RBC-contaminated cells
# percent.ribo: <50% to exclude extreme ribosomal outliers 

nasal <- subset(
  nasal,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 4000 &
    nCount_RNA > 1000 &
    nCount_RNA < 10000 &
    percent.mt < 15 &
    percent.hb < 5 &
    percent.ribo < 50
)

table(nasal$orig.ident)

vln_plot_3 <- VlnPlot(nasal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb", "percent.ribo"), ncol = 5)

# ggsave(
#   filename = "../figures/QC_violin_plot_3.png",
#   plot = vln_plot_3,
#   width = 10,
#   height = 4,
#   dpi = 300
# )

# QC Results.
# Filtering removed 603 cells (0.39%) from 156,572 total, low removal rate
# confirming data was already high quality prior to filtering.
# Pre vs Post filtering cell counts:
# D02: 40,148 -> 40,003
# D05: 26,344 -> 26,235
# D08: 29,523 -> 29,416
# D14: 30,203 -> 30,054 
# Naive: 30,354 -> 30,261

# Normalization - log normalize 
nasal <- NormalizeData(nasal, normalization.method = "LogNormalize")

# Selecting Highly Variable Genes 
nasal <- FindVariableFeatures(nasal, selection.method = "vst", nfeatures = 2000)

# Scaling 
all.genes <- rownames(nasal)

# Scaling HVGs
nasal <- ScaleData(nasal)

# PCA 
nasal <- RunPCA(nasal)
elbow <- ElbowPlot(nasal, ndims = 50) 

ggsave(
  filename = "../figures/elbow_plot.png",
  plot = elbow,
  width = 10,
  height = 4,
  dpi = 300
)

# Get variance explained
var_explained <- nasal[["pca"]]@stdev^2
var_percent <- round(100 * var_explained / sum(var_explained), 1)

pca <- DimPlot(nasal, reduction = "pca") + ggtitle("PC1-PC2") +   
  xlab(paste0("PC1 (", var_percent[1], "%)")) +
  ylab(paste0("PC2 (", var_percent[2], "%)")) +
  theme(
    plot.title = element_text(hjust = 0.7)  # center the title
  )

# PCA shows clear structure with cells forming distinct branches, suggesting
# underlying biological heterogeneity. Samples are well mixed across the space,
# indicating minimal batch effects and that variation is driven by biology.

ggsave(
  filename = "../figures/pca_plot.png",
  plot = pca,
  width = 8,
  height = 6,
  dpi = 300
)

# Through evaluating the Elbow plot the first 40 PCs will be used 
# Find Neighbors using the first 40 PCs
set.seed(123)

# Through evaluating the Elbow plot, the first 40 PCs will be used
nasal <- FindNeighbors(nasal, dims = 1:40)

# FindClusters at multiple resolutions to determine optimal clustering
nasal <- FindClusters(nasal, resolution = seq(0.1, 0.8, 0.1))

# Run UMAP using the first 40 PCs
nasal <- RunUMAP(nasal, dims = 1:40)

# Plot UMAPs for each resolution for comparison
resolutions <- seq(0.1, 0.8, 0.1)

umap_list <- lapply(resolutions, function(res) {
  DimPlot(
    nasal,
    group.by = paste0("RNA_snn_res.", res),
    label = TRUE
  ) + ggtitle(paste("Resolution", res))
})

# Combine plots
combined_umap <- wrap_plots(umap_list, ncol = 3)

combined_umap

ggsave(
  "../figures/umap_resolutions.png",
  combined_umap,
  width = 18,
  height = 12,
  dpi = 300
)

# Clustree to visualize cluster stability across resolutions
library(clustree)
clust <- clustree(nasal, prefix = "RNA_snn_res.")

ggsave(
  "../figures/cluster_tree.png",
  clust,
  width = 18,
  height = 12,
  dpi = 300
)

# Resolution of 0.4 was picked
Idents(nasal) <- "RNA_snn_res.0.4"

# Visual check for batch effects
umap_batch <- DimPlot(nasal, group.by = "orig.ident")

ggsave("../figures/umap_orig.png", umap_batch, width = 8, height = 6, dpi = 300)

umap_full <- DimPlot(nasal, reduction = "umap", label = TRUE)
ggsave("../figures/umap_full.png", umap_full, width = 8, height = 6, dpi = 300)
# UMAP colored by sample shows cells from all conditions are well mixed within 
# clusters, indicating minimal batch effects. No batch correction was applied.

# Cell-type annotations
# Find all marker genes

nasal.markers <- FindAllMarkers(
  nasal,
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = 0.2,
  max.cells.per.ident = 1000  # chosen a bit lower to reduce runtime and match paper's methods 
)

saveRDS(nasal.markers, file = "../results/nasal_markers.rds")
saveRDS(nasal, file = "../results/nasal_processed.rds")

#---------------------------------------------------------------
nasal <- readRDS("../results/nasal_processed.rds")
nasal.markers <- readRDS("../results/nasal_markers.rds")

set.seed(123)

top5 <- nasal.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

# Feature Plots
# For Cluster 0
fp_0 <- FeaturePlot(
  nasal,
  features = c("Ncam2", "Cidea", "Gldc", "Mdga2"),
  reduction = "umap"
) 

ggsave(
  "../figures/fp_0.png",
  fp_0,
  width = 18,
  height = 12,
  dpi = 300
)

# For Cluster 1  
fp_1 <- FeaturePlot(
  nasal,
  features = c("Nqo1", "Cartpt", "Acsm4", "Hs3st5"),
  reduction = "umap"
) 

ggsave(
  "../figures/fp_1.png",
  fp_1,
  width = 18,
  height = 12,
  dpi = 300
)

# For Cluster 2 - 
fp_2 <- FeaturePlot(
  nasal,
  features = c("Pf4", "Ms4a7", "Mrc1", "Trem2"),
  reduction = "umap"
) 

ggsave(
  "../figures/fp_2.png",
  fp_2,
  width = 18,
  height = 12,
  dpi = 300
)
# > top5 %>% filter(cluster==17) the top 5 genes were highlighted 
# per chosen cluster and manually annotated 
# the group of genes were researched and then compared to the 
# annotated labels found here: https://singlecell.broadinstitute.org/single_cell/study/SCP2216/primary-nasal-viral-infection-rewires-the-tissue-scale-memory-response-primary-infection-dataset
# Manual annotation of main clusters

levels(nasal)

# create a label vector
new.cluster.ids <- levels(nasal)

# manually assign new cluster labels based on manual annotation comparison
new.cluster.ids[as.numeric(new.cluster.ids) == 0] <- "Epithelial cells"
new.cluster.ids[as.numeric(new.cluster.ids) == 2] <- "Macrophages"
new.cluster.ids[as.numeric(new.cluster.ids) == 3] <- "Basal epithelial cells"
new.cluster.ids[as.numeric(new.cluster.ids) == 4] <- "B cells"
new.cluster.ids[as.numeric(new.cluster.ids) == 5] <- "Dendritic cells"
new.cluster.ids[as.numeric(new.cluster.ids) == 6] <- "Endothelial cells"
new.cluster.ids[as.numeric(new.cluster.ids) == 9] <- "Fibroblasts"
new.cluster.ids[as.numeric(new.cluster.ids) == 10] <- "T cells"
new.cluster.ids[as.numeric(new.cluster.ids) == 25] <- "Natural killer cells"

# apply new labels
names(new.cluster.ids) <- levels(nasal)
nasal <- RenameIdents(nasal, new.cluster.ids)

# save to the metadata
nasal$cell_type_manual <- Idents(nasal)

# Set identities to manual labels
Idents(nasal) <- "cell_type_manual"


# Replot
p_manual <- DimPlot(nasal, label = TRUE, repel = TRUE) + NoLegend()

ggsave("../figures/umap_manual_annotations.png",
       p_manual, width = 10, height = 8, dpi = 300)


nasal$cell_type_manual <- as.character(nasal$cell_type_manual)
nasal$cell_type_manual[grepl("^[0-9]+$", nasal$cell_type_manual)] <- "Other"
nasal$cell_type_manual <- as.factor(nasal$cell_type_manual)
Idents(nasal) <- "cell_type_manual"
table(nasal$cell_type_manual)

# Automated Annotation 
library(SingleR)
library(celldex)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(nasal)
ref <- MouseRNAseqData()

pred <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main
)

nasal$cell_type_singleR <- pred$labels

table(nasal$cell_type_manual, nasal$cell_type_singleR)

p_singleR <- DimPlot(
  nasal,
  group.by = "cell_type_singleR",
  label = TRUE,
  repel = TRUE
) +
  theme(
    text = element_text(size = 14),
    legend.text = element_text(size = 10)
  )

ggsave(
  "../figures/umap_singleR.png",
  p_singleR,
  width = 12,
  height = 9,
  dpi = 300
)

p_manual <- DimPlot(
  nasal,
  group.by = "cell_type_manual",
  label = TRUE,
  repel = TRUE
) + ggtitle("Manual Annotation")

p_singleR <- DimPlot(
  nasal,
  group.by = "cell_type_singleR",
  label = TRUE,
  repel = TRUE
) + ggtitle("SingleR Annotation")

combined <- p_manual + p_singleR

ggsave(
  "../figures/umap_comparison.png",
  combined,
  width = 16,
  height = 8,
  dpi = 300
)

# Differential Analysis 
# comparing conditions 5dpi vs 14 dpi 
# looking at macrophage, dendritic, NK, and T cells 

# remove missing mouse_ids
nasal$mouse_id[is.na(nasal$mouse_id)] <- "unknown"
nasal <- subset(nasal, subset = mouse_id != "unknown")

# subset comparison 
nasal_sub <- subset(nasal, subset = time %in% c("D05", "D14"))

cells_of_interest <- c("T cells", "Natural killer cells", "Dendritic cells", "Macrophages")

nasal_sub <- subset(nasal_sub, idents = cells_of_interest)

# pseudobulk aggregation
pseudo_nasal <- AggregateExpression(
  nasal_sub,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("time", "mouse_id", "cell_type_manual")
)

# create identities
pseudo_nasal$celltype.time <- paste(
  pseudo_nasal$cell_type_manual,
  pseudo_nasal$time,
  sep = "_"
)

Idents(pseudo_nasal) <- "celltype.time"

levels(Idents(pseudo_nasal))

library(DESeq2)

# Running DESeq2
# DE for the Macrophages
de_mac <- FindMarkers(
  pseudo_nasal,
  ident.1 = "Macrophages_D14",
  ident.2 = "Macrophages_D05",
  test.use = "DESeq2"
)

head(de_mac)

# DE for the T cells
de_t <- FindMarkers(
  pseudo_nasal,
  ident.1 = "T cells_D14",
  ident.2 = "T cells_D05",
  test.use = "DESeq2"
)

# DE for the NK cells
de_nk <- FindMarkers(
  pseudo_nasal,
  ident.1 = "Natural killer cells_D14",
  ident.2 = "Natural killer cells_D05",
  test.use = "DESeq2"
)

# DE for the Dendritic cells
de_dc <- FindMarkers(
  pseudo_nasal,
  ident.1 = "Dendritic cells_D14",
  ident.2 = "Dendritic cells_D05",
  test.use = "DESeq2"
)

set.seed(123)
# heatmap for macrophage cluster
plot_heatmap_for_celltype <- function(de_results, celltype_name, pseudo_obj) {
  
  library(dplyr)
  library(pheatmap)
  
  # Get top genes
  top_up_high <- de_results %>%
    filter(p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    head(10)
  
  top_up_low <- de_results %>%
    filter(p_val_adj < 0.05) %>%
    arrange(avg_log2FC) %>%
    head(10)
  
  top_genes <- c(rownames(top_up_high), rownames(top_up_low))
  
  # subset
  pseudo_sub <- subset(pseudo_obj, subset = cell_type_manual == celltype_name)
  
  # rrder time
  pseudo_sub$time <- factor(pseudo_sub$time, levels = c("D05", "D14"))
  pseudo_sub <- pseudo_sub[, order(pseudo_sub$time)]
  
  # matrix information
  mat <- GetAssayData(pseudo_sub, layer = "data")[top_genes, ]
  mat_scaled <- t(scale(t(mat)))
  
  # annotation
  annotation_col <- data.frame(time = pseudo_sub$time)
  rownames(annotation_col) <- colnames(mat_scaled)
  
  # clean filename
  file_name <- paste0(
    "../figures/",
    gsub(" ", "_", tolower(celltype_name)),
    "_heatmap.png"
  )
  
  # 7. Save directly
  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = paste(celltype_name, ": D05 vs D14"),
    filename = file_name,
    width = 8,
    height = 6
  )
}

plot_heatmap_for_celltype(de_mac, "Macrophages", pseudo_nasal)
plot_heatmap_for_celltype(de_t, "T cells", pseudo_nasal)
plot_heatmap_for_celltype(de_nk, "Natural killer cells", pseudo_nasal)
plot_heatmap_for_celltype(de_dc, "Dendritic cells", pseudo_nasal)

# ORA - Macrophages
run_ORA <- function(de_results, celltype_name) {
  
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  
  sig_genes <- de_results %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
  
  gene_symbols <- rownames(sig_genes)
  
  gene_df <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )
  
  ego <- enrichGO(
    gene = gene_df$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  p <- dotplot(ego, showCategory = 10) +
    ggtitle(paste("GO Enrichment in", celltype_name))
  
  print(p)
  
  ggsave(
    paste0("../figures/", gsub(" ", "_", tolower(celltype_name)), "_GO.png"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

run_ORA(de_mac, "Macrophages")
run_ORA(de_t, "T cells")
run_ORA(de_nk, "Natural killer cells")
run_ORA(de_dc, "Dendritic cells")

# save final object
saveRDS(nasal, file = "../results/nasal_final.rds")
