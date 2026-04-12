library(Seurat) 
library(dplyr)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(org.Mm.eg.db)
library(BiocParallel)

# load in Seurat object 
nasal <- readRDS("../seurat_ass4.rds")

##### Nasal Object #####
nasal

set.seed(123)

table(nasal$orig.ident)

# Mitochondrial, Hemoglobin, and Ribosomal protein Percentage 
nasal[["percent.mt"]] <- PercentageFeatureSet(nasal, pattern = "^mt-")

nasal[["percent.hb"]] <- PercentageFeatureSet(nasal, pattern = "^Hb[ab]")

nasal[["percent.ribo"]] <- PercentageFeatureSet(nasal, pattern = "^Rp[sl]")


##### Quality Control Visualizations #####

# Violin Plots
vln_plot_1 <- VlnPlot(
  nasal,
  features = c("nFeature_RNA", "nCount_RNA"),
  ncol = 2,
  pt.size = 0
)

p1 <- VlnPlot(nasal, features = "percent.mt", pt.size = 0) + NoLegend()
p2 <- VlnPlot(nasal, features = "percent.ribo", pt.size = 0) + NoLegend()
p3 <- VlnPlot(nasal, features = "percent.hb", pt.size = 0.5,) + 
  ylim(0, 6) 

vln_plot_2 <- p1 | p2 | p3
vln_plot_2

ggsave(
  filename = "../figures/QC_violin_plot_1.png",
  plot = vln_plot_1,
  width = 10,
  height = 4,
  dpi = 300
)

ggsave(
  filename = "../figures/QC_violin_plot_2.png",
  plot = vln_plot_2,
  width = 10,
  height = 4,
  dpi = 300
)

# Knee Plot
nasal <- CalculateBarcodeInflections(nasal)

knee_plot <- BarcodeInflectionsPlot(nasal)

ggsave(
  filename = "../figures/QC_knee_plot.png",
  plot = knee_plot,
  width = 6,
  height = 5,
  dpi = 300
)

# Based on QC plots:
# nFeature_RNA: remove low-quality cells (<500) and likely doublets (>4000)
# nCount_RNA: remove near-empty droplets (<1000) and outlier high-count cells (>10000)
# percent.mt: <15% threshold, to account for infection and based on violin plots
# percent.hb: <5% to remove rare RBC-contaminated cells
# percent.ribo: <50% to exclude extreme ribosomal outliers 

nasal <- subset(
  nasal,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 3800 &
    nCount_RNA > 1000 &
    nCount_RNA < 10000 &
    percent.mt < 15 &
    percent.hb < 5 &
    percent.ribo < 50
)

table(nasal$orig.ident)

p1 <- VlnPlot(nasal, features = "nFeature_RNA", pt.size = 0) + NoLegend()
p2 <- VlnPlot(nasal, features = "nCount_RNA", pt.size = 0) + NoLegend()
p3 <- VlnPlot(nasal, features = "percent.mt", pt.size = 0) + NoLegend()
p4 <- VlnPlot(nasal, features = "percent.hb", pt.size = 0.1) + ylim(0, 6) + NoLegend()
p5 <- VlnPlot(nasal, features = "percent.ribo", pt.size = 0)

vln_plot_3 <- p1 | p2 | p3 | p4 | p5
vln_plot_3

ggsave(
  filename = "../figures/QC_violin_plot_3.png",
  plot = vln_plot_3,
  width = 18,
  height = 4,
  dpi = 300
)

# QC Results
# Filtering removed 606 cells (0.39%) from 156,572 total, low removal rate
# confirming data was already high quality prior to filtering.
# Pre vs Post filtering cell counts
# D02: 40,148 -> 40,002
# D05: 26,344 -> 26,235
# D08: 29,523 -> 29,412
# D14: 30,203 -> 30,054 
# Naive: 30,354 -> 30,261

##### Filtering  #####

# Normalization - log normalize 
nasal <- NormalizeData(nasal, normalization.method = "LogNormalize")

# Selecting Highly Variable Genes 
nasal <- FindVariableFeatures(nasal, selection.method = "vst", nfeatures = 2000)

# Scaling 
all.genes <- rownames(nasal)

# Scaling HVGs
nasal <- ScaleData(nasal)

##### PCA #####
nasal <- RunPCA(nasal)
elbow <- ElbowPlot(nasal, ndims = 50) 

ggsave(
  filename = "../figures/elbow_plot.png",
  plot = elbow,
  width = 8,
  height = 5,
  dpi = 300
)

# Get variance explained
var_explained <- nasal[["pca"]]@stdev^2
var_percent <- round(100 * var_explained / sum(var_explained), 1)

pca <- DimPlot(nasal, reduction = "pca") + ggtitle("PC1-PC2") +   
  xlab(paste0("PC1 (", var_percent[1], "%)")) +
  ylab(paste0("PC2 (", var_percent[2], "%)")) +
  theme(
    plot.title = element_text(hjust = 0.5) 
  )

# PCA shows clear structure with cells forming distinct branches, 
# suggesting meaningful patterns of variation were captured 
# also indicating minimal batch effects and that variation is driven by biology.

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

#### Determining UMAP Resolution ####
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

# Adjust mouse_id values so the unassigned has a label 
nasal$mouse_id <- as.character(nasal$mouse_id)

nasal$mouse_id[nasal$mouse_id == "" | is.na(nasal$mouse_id)] <- "Unassigned"

nasal$mouse_id <- as.factor(nasal$mouse_id)

table(nasal$mouse_id)

##### Batch Effects #####
# Visual check for batch effects
umap_batch <- DimPlot(nasal, group.by = "mouse_id")
ggsave("../figures/umap_batch.png", umap_batch, width = 18, height = 12, dpi = 300)

##### UMAP #####
umap_full <- DimPlot(nasal, reduction = "umap", label = TRUE)
ggsave("../figures/umap_full.png", umap_full, width = 8, height = 6, dpi = 300)
# UMAP colored by sample shows cells from all conditions are well mixed within 
# clusters, indicating minimal batch effects. No batch correction was applied.

##### UMAP across Tissues ####
p1 <- DimPlot(nasal, label = TRUE) + 
  ggtitle("Base UMAP") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) + NoLegend()

p2 <- DimPlot(nasal, group.by = "organ_custom") + 
  ggtitle("Tissue Type")  + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) 

ggsave("../figures/umap_tissue.png", p2, width = 8, height = 6, dpi = 300) 

p1 + p2

ggsave("../figures/umap_cluster_tissue.png", p1+p2, width = 18, height = 10, dpi = 300)


# saveRDS(nasal, file = "../results/nasal_premarkers.rds")
##### Cell-type annotations #####
# Find all marker genes
nasal <- readRDS("../results/nasal_premarkers.rds")


nasal.markers <- FindAllMarkers(
  nasal,
  only.pos = TRUE,
  logfc.threshold = 0.5,
  min.pct = 0.2,
  max.cells.per.ident = 1000  # chosen a bit lower to reduce runtime and match paper's methods 
)

#saveRDS(nasal.markers, file = "../results/nasal_markers.rds")
#saveRDS(nasal, file = "../results/nasal_processed.rds")

# nasal <- readRDS("../results/nasal_processed.rds")
# nasal.markers <- readRDS("../results/nasal_markers.rds")

# saveRDS(nasal.markers, file = "../results/nasal_markers2.rds")
nasal <- readRDS("../results/nasal_premarkers.rds")
nasal.markers <- readRDS("../results/nasal_markers2.rds")

set.seed(123)

top5 <- nasal.markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

# Cluster 0 - Olfactory Neurons
top5 %>% filter(cluster == 0) 

# show a true biological signal with a high pct.1 and low pct.2
genes_0 <- c("Ncam2", "Gldc", "Cidea", "Mdga2")

# Feature plots 
fp_0 <- FeaturePlot(nasal, features = genes_0, reduction = "umap", ncol = 2) 

ggsave("../figures/fp_0.png", fp_0, width = 10, height = 8, dpi = 300)

# Cluster 1 - Sustenacular cells
top5 %>% filter(cluster == 1) 

fp_1 <- FeaturePlot(
  nasal,
  features = c("Nqo1", "Cartpt", "Acsm4", "Kirrel3"),
  reduction = "umap"
) 

ggsave("../figures/fp_1.png", fp_1, width = 10, height = 8, dpi = 300)

# For Cluster 2 - Macrophages
top5 %>% filter(cluster == 2) 

fp_2 <- FeaturePlot(
  nasal,
  features = c("Fcrls", "Ms4a7", "Mrc1", "Trem2"),
  reduction = "umap"
) 

ggsave(
  "../figures/fp_2.png",
  fp_2,
  width = 18,
  height = 12,
  dpi = 300
)

# Inspection was done for more clusters but not shown in FeaturePlots
top5 %>% filter(cluster == 3)
top5 %>% filter(cluster == 4)
top5 %>% filter(cluster == 5)
top5 %>% filter(cluster == 6)
top5 %>% filter(cluster == 7)
top5 %>% filter(cluster == 8)
top5 %>% filter(cluster == 9)
top5 %>% filter(cluster == 10)
top5 %>% filter(cluster == 11)
top5 %>% filter(cluster == 12)
top5 %>% filter(cluster == 13)
top5 %>% filter(cluster == 14)
top5 %>% filter(cluster == 15)
top5 %>% filter(cluster == 19)
top5 %>% filter(cluster == 25)


# the group of genes were researched and then compared to the 
# annotated labels found here: https://singlecell.broadinstitute.org/single_cell/study/SCP2216/primary-nasal-viral-infection-rewires-the-tissue-scale-memory-response-primary-infection-dataset
# Manual annotation of main clusters

levels(nasal)

# create a label vector
new.cluster.ids <- levels(nasal)

#### Manual Annotation #####
# manually assign new cluster labels based on manual annotation comparison
# some clusters were merged together for simplicity as they contain the same cell type at a high level
new.cluster.ids[new.cluster.ids == "0"]  <- "Olfactory Neurons"
new.cluster.ids[new.cluster.ids %in% c("1", "7", "8", "14")] <- "Sustenacular"
new.cluster.ids[new.cluster.ids == "2"]  <- "Macrophages"
new.cluster.ids[new.cluster.ids == "3"]  <- "Basal cells"
new.cluster.ids[new.cluster.ids == "4"]  <- "B cells"
new.cluster.ids[new.cluster.ids == "5"]  <- "Monocytes"
new.cluster.ids[new.cluster.ids == "6"]  <- "Endothelial"
new.cluster.ids[new.cluster.ids == "9"]  <- "Fibroblasts"
new.cluster.ids[new.cluster.ids == "10"] <- "T cells"
new.cluster.ids[new.cluster.ids %in% c("12", "18", "23")] <- "Neutrophils"
new.cluster.ids[new.cluster.ids == "13"] <- "Dendritic"
new.cluster.ids[new.cluster.ids == "19"] <- "Pericytes"
new.cluster.ids[new.cluster.ids == "25"] <- "NK Cells"

# apply new labels
names(new.cluster.ids) <- levels(nasal)
nasal <- RenameIdents(nasal, new.cluster.ids)

# save to the metadata
nasal$cell_type_manual <- Idents(nasal)

# ensure the unidentified clusters are "Other:
nasal$cell_type_manual <- as.character(nasal$cell_type_manual)
nasal$cell_type_manual[grepl("^[0-9]+$", nasal$cell_type_manual)] <- "Other"
nasal$cell_type_manual <- as.factor(nasal$cell_type_manual)

# Set identities to manual labels
Idents(nasal) <- "cell_type_manual"

# Replot with new annotations
p_manual <- DimPlot(nasal, label = TRUE, repel = TRUE, label.size = 4.5, label.color = "black") + ggtitle("Manual Annotation") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) 
p_manual

ggsave("../figures/umap_manual_annotations.png", p_manual, width = 11, height = 8, dpi = 300)

# to double check 
table(nasal$cell_type_manual)

#### Automated Annotation #####
# convert nasal to the right format
sce <- as.SingleCellExperiment(nasal)

ref_mouse <- MouseRNAseqData()

ref_immune  <- ImmGenData()

param <- SnowParam(workers = 4, type = "SOCK") 

pred <- SingleR(
  test = sce,
  ref = list(Gen = ref_mouse, Imm = ref_immune),
  labels = list(ref_mouse$label.main, ref_immune$label.main),
  BPPARAM = param
)

nasal$cell_type_singleR <- pred$labels

# add the automated annotations to the object 
annotation_table <- table(nasal$cell_type_manual, nasal$cell_type_singleR)
annotation_table

table(nasal$cell_type_manual, nasal$cell_type_singleR)

p_singleR <- DimPlot(
  nasal,
  group.by = "cell_type_singleR",
  label = TRUE,
  repel = TRUE,
  label.size = 4.5
) + ggtitle("SingleR Annotation") + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)) 


ggsave(
  "../figures/umap_singleR.png",
  p_singleR,
  width = 12,
  height = 9,
  dpi = 300
)

combined <- p_manual + p_singleR

ggsave(
  "../figures/umap_comparison.png",
  combined,
  width = 16,
  height = 8,
  dpi = 300
)

##### Differential Analysis ####
# comparing conditions 5dpi vs 14 dpi, chosen as this is right in the immune response phase and then recovery 
# looking at macrophage, dendritic, monocytes, neutrophils, NK, and T cells 

# remove missing mouse_ids
nasal <- subset(nasal, subset = mouse_id != "Unassigned")

# to check they have been removed
DimPlot(nasal, group.by = "mouse_id")

# subset comparison 
nasal_sub <- subset(nasal, subset = time %in% c("D05", "D14"))

cells_of_interest <- c("T cells", "NK Cells", "Dendritic", "Macrophages", "Monocytes", "Neutrophils")

nasal_sub <- subset(nasal_sub, idents = cells_of_interest)

# Comparison over time
dpi5_14 <- DimPlot(
  nasal_sub, 
  split.by = "time", 
  group.by = "cell_type_manual", 
  label = TRUE, 
  repel = TRUE
) + 
  ggtitle("DPI 5 vs DPI 14") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  "../figures/dpi_5_14_comparison.png",
  dpi5_14,
  width = 10,
  height = 8,
  dpi = 300
)

# pseudobulk aggregation
pseudo_nasal <- AggregateExpression(
  nasal_sub,
  assays = "RNA",
  return.seurat = TRUE,
  group.by = c("time", "mouse_id", "cell_type_manual")
)

# create identities
pseudo_nasal$celltype.time <- paste(pseudo_nasal$cell_type_manual, pseudo_nasal$time, sep = "_")

Idents(pseudo_nasal) <- "celltype.time"

levels(Idents(pseudo_nasal))

# Plotting heatmap function
plot_heatmap <- function(de_results, celltype_name, pseudo_obj, top_n = 10) {

  # Filter for significant genes (padj < 0.05)
  sig_res <- de_results %>% 
    filter(p_val_adj < 0.05)
  
  if (nrow(sig_res) < 2) {
    print(paste(celltype_name, "- no significant genes found."))
    return(NULL)
  }
  
  # get top 10 Up and top 10 Down of the DPI states 5 and 14
  # arranging by log2FC score
  top_up <- sig_res %>% arrange(desc(avg_log2FC)) %>% head(top_n)
  top_down <- sig_res %>% arrange(avg_log2FC) %>% head(top_n)
  top_genes <- c(rownames(top_up), rownames(top_down))
  
  # Subset the pseudobulk for this specific cell type
  pseudo_sub <- subset(pseudo_obj, subset = cell_type_manual == celltype_name)
  
  # order 5 on the left and 14 on the right
  pseudo_sub$time <- factor(pseudo_sub$time, levels = c("D05", "D14"))
  pseudo_sub <- pseudo_sub[, order(pseudo_sub$time)]
  
  # extract and scale for visualizations
  mat <- GetAssayData(pseudo_sub, layer = "data")[top_genes, , drop = FALSE]
  mat_scaled <- t(scale(t(as.matrix(mat))))
  
  # plotting
  file_path <- paste0("../figures/", gsub(" ", "_", tolower(celltype_name)), "_DE_heatmap.png")
  
  pheatmap(
    mat_scaled,
    annotation_col = data.frame(time = pseudo_sub$time, row.names = colnames(mat_scaled)),
    cluster_cols = FALSE,
    show_colnames = FALSE,
    main = paste(celltype_name, ": DPI 5 vs DPI 14"),
    filename = file_path,
    width = 7,
    height = 6
  )
}

# Run the differential analysis 

# macrophages
de_mac <- FindMarkers(pseudo_nasal, ident.1 = "Macrophages_D14", ident.2 = "Macrophages_D05", test.use = "DESeq2")
plot_heatmap(de_mac, "Macrophages", pseudo_nasal)

# monocytes
de_mono <- FindMarkers(pseudo_nasal, ident.1 = "Monocytes_D14", ident.2 = "Monocytes_D05", test.use = "DESeq2")
plot_heatmap(de_mono, "Monocytes", pseudo_nasal)

# neutrophils
de_neutro <- FindMarkers(pseudo_nasal, ident.1 = "Neutrophils_D14", ident.2 = "Neutrophils_D05", test.use = "DESeq2")
plot_heatmap(de_neutro, "Neutrophils", pseudo_nasal)

# NK cells
de_nk <- FindMarkers(pseudo_nasal, ident.1 = "NK Cells_D14", ident.2 = "NK Cells_D05", test.use = "DESeq2")
plot_heatmap(de_nk, "NK Cells", pseudo_nasal)

# T cells
de_t <- FindMarkers(pseudo_nasal, ident.1 = "T cells_D14", ident.2 = "T cells_D05", test.use = "DESeq2")
plot_heatmap(de_t, "T cells", pseudo_nasal)

# another feature plot showing the differences in time for the macrophage gene Isg15 and Hspa1a
# Isg15 - High D5 and Hes1 - High D14
fp_I_H <- FeaturePlot(
  nasal, 
  features = c("Isg15", "Mrc1"), 
  split.by = "time", 
  order = TRUE, 
  cols = c("lightgrey", "mediumvioletred")
) + theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  strip.text = element_text(face = "bold", size = 12) )

ggsave(
  "../figures/fp_I_H_comparison.png",
  fp_I_H,
  width = 15,
  height = 8,
  dpi = 300
)

fp_neutro <- FeaturePlot(
  nasal, 
  features = c("Phf11d", "Hk2"), 
  split.by = "time", 
  order = TRUE, 
  cols = c("lightgrey", "mediumvioletred")
) + theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  strip.text = element_text(face = "bold", size = 12) )

ggsave(
  "../figures/fp_neutro_comp.png",
  fp_neutro,
  width = 15,
  height = 8,
  dpi = 300
)

vln_tissue <- VlnPlot(
  nasal, 
  features = c("Phf11d", "Hk2"), 
  group.by = "organ_custom"
)

ggsave("../figures/vln_neutro_tissue.png", vln_tissue, width = 10, height = 6)

##### Enrichment #####
# ORA - Macrophages
run_ora <- function(de_results, celltype_name) {

  # filter for D14 Up-regulated genes 
  sig_genes <- de_results %>% 
    filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>%
    rownames()
  
  # Entrez Conversion
  gene_df <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  # Enrichment
  ego <- enrichGO(gene = gene_df$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", readable = TRUE)
  
  # Plot
  p <- dotplot(ego, showCategory = 10) + 
    ggtitle(paste(celltype_name, "ORA - D14 vs D05)")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  file_name <- paste0("../figures/", gsub(" ", "_", tolower(celltype_name)), "_final_GO.png")
  ggsave(file_name, plot = p, width = 8, height = 6, dpi = 300)
  
  return(p)
}

run_ora(de_mac, "Macrophages")
run_ora(de_t, "T cells")

saveRDS(nasal, file = "../results/nasal_final.rds")
nasal <- readRDS("../results/nasal_final.rds")
