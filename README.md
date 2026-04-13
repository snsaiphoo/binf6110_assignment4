# Comparative scRNA-seq of Murine Nasal Mucosa

## Introduction 
A primary function of the nasal airway is to continuously filter pathogens, as this is an entry point for viral infections [1]. This is considered an immune defense site, which makes the nasal mucosa an important area of research [1]. More specifically, the epithelial tissue lining the nasal pathway actively signals and recruits immune cells during infection [2], along with undergoing structural and cellular changes to better protect the mucosa [3]. However, many viruses, including influenza, are capable of infiltrating epithelial cells and using them as a primary site of infection [2]. Understanding the cell states and transcriptional changes of these cells can provide insight into host–virus interactions and support the development of vaccines and therapeutic strategies [3].

Many acute viral infections follow a characteristic pattern, starting with a rapid increase to peak viral levels, a slower decline as the immune response begins to control the infection, and finally a rapid decrease associated with viral clearance [4, 5].

In this investigation, a single-cell RNA sequencing (scRNA-seq) analysis will be conducted to visualize changes in cellular populations over the course of influenza A infection and recovery. The dataset consists of respiratory mucosa (RM), olfactory mucosa (OM), and lateral nasal gland (LNG) tissues collected at five stages: uninfected (naive), peak viral load (2 days post-infection, dpi), myeloid recruitment (5 dpi), T-cell infiltration (8 dpi), and viral clearance (14 dpi) [3]. The aim of this study is to examine how cell populations and their gene expression change across these stages of infection, to better understand how epithelial and immune cells respond to and recover from influenza A infection.

### Design Rationale
A typical scRNA-seq dataset is high-dimensional and requires several quality control steps, including filtering, normalization, and dimensionality reduction. Once preprocessing is complete, the data will be ready for clustering, cell type annotation, differential expression analysis, and functional enrichment. In this study, these steps will be carried out sequentially to characterize how cell populations change across infection stages.

Selecting an appropriate scRNA-seq analysis framework involves balancing computational efficiency, accessibility, and biological reliability. Commonly used frameworks include Seurat (R), Scanpy (Python), and Bioconductor-based workflows such as OSCA. One benefit of Scanpy is that it handles large datasets efficiently and integrates well within the Python ecosystem [6]. While OSCA provides a more statistically rigorous framework through Bioconductor [7]. Benchmarking studies indicate that differences in performance between these tools are generally small for typical dataset sizes, even when accounting for the speed advantages of GPU-based methods [8]. Seurat, however, offers a more integrated and accessible workflow that covers preprocessing, clustering, and visualization within a single framework, and recent updates have improved its scalability and support for multimodal data [9]. Its widespread adoption and accessible documentation also support reproducibility. For these reasons, Seurat will be used for this analysis.

Cell type annotation allows clusters to be given biological meaning, but manual annotation is often subjective and inconsistent. Several automated tools have been developed to address this, including Seurat, SingleR, scmap, and CHETAH [9, 10, 11, 12]. Benchmarking studies show that most methods perform reasonably well overall, but accuracy tends to drop when distinguishing rare or closely related cell types [10]. Seurat identifies major cell types reliably but is less suited to resolving subtle or rare populations [10]. SingleR, by contrast, assigns cell identities by comparing expression profiles to reference transcriptomic datasets, which produces more consistent and objective annotations [13]. For these reasons, SingleR will be used to annotate cell types in this study.

Differential gene expression (DGE) analysis identifies genes that change between conditions, but results can be affected by variability at both the cell and subject levels. Standard cell-level approaches, such as those used in Seurat, treat each cell as an independent observation, which can inflate false discovery rates [9,14]. Pseudobulk methods avoid this problem by summing counts across cells within each sample before performing statistical testing, making the analysis more reliable [14]. DGE analysis will therefore be carried out using a pseudobulk approach, in which DESeq2 is called through the Seurat `FindMarkers` function [9], combining the accessibility of Seurat’s workflow with a pseudobulk-based differential expression strategy.

To determine the biological significance of the differentially expressed genes identified between infection time points, over-representation analysis (ORA) will be performed using clusterProfiler [15]. For each comparison between infection stages, ORA will identify which biological pathways and processes are significantly enriched, helping to clarify how the cellular response changes as infection progresses. This will make it possible to pinpoint stage-specific changes, such as shifts in immune activation or inflammatory signaling, that would not be obvious from looking at individual gene results alone.

## Methods
The scRNA-seq analysis consists of a single R script located in the [`scripts`](scripts/) folder, which can be run to execute the steps detailed below. The dataset can be accessed through the article "Primary nasal influenza infection rewires tissue-scale memory response dynamics" by Kazer et al. [3], available through the [original publication](https://doi.org/10.1016/j.immuni.2024.06.005).

### 1.0 - Data & Tools

#### 1.1 - Data Acquisition
The dataset used in this analysis was provided as a pre-made Seurat object, consisting of 156,572 cells and 25,129 genes across three tissue types collected at five infection stages.

#### 1.2 - R Environment
Analysis was performed using R version 4.5.1. All required CRAN and Bioconductor dependencies are specified in the [`environment.R`](scripts/environment.R) script, which can be used to reproduce the R environment used for this analysis, including package version numbers.

### 2.0 - Data Preprocessing

#### 2.1 - Removing Low Quality Cells
Multiple vignettes were consulted to determine appropriate preprocessing methods for scRNA-seq analysis, including the Seurat vignette [16]. Quality control (QC) filtering was performed to remove low-quality cells prior to downstream analysis. Commonly used QC metrics include the number of detected genes and the proportion of mitochondrial gene expression, where high mitochondrial proportions indicate dying or stressed cells and low gene counts suggest damaged cells [16]. A knee plot was used to identify a natural inflection point to guide filtering thresholds [17]. Cells with high hemoglobin gene expression were filtered to remove red blood cell contamination [18], and cells with extreme ribosomal gene expression were excluded as an additional quality criterion [16].

Metrics were visualized using violin plots and a knee plot to inform threshold selection, and cutoffs were chosen based on these distributions. The following criteria were applied:

- `nFeature_RNA` between 500 and 3,800
- `nCount_RNA` between 1,000 and 10,000
- `percent.mt` below 15%
- `percent.hb` below 5%
- `percent.ribo` below 50%

#### 2.2 - Normalization and Selecting Highly Variable Genes
Following filtering, gene expression counts were normalized using log normalization via `NormalizeData()`, which scales each cell to a total count of 10,000 before applying a log transformation. Highly variable genes (HVGs) were then identified using `FindVariableFeatures()` with the variance-stabilizing transformation (VST) method, selecting the top 2,000 most variable genes for downstream analysis. The data was subsequently scaled using `ScaleData()` to ensure equal contribution of each gene during dimensionality reduction [16].

#### 2.3 - Dimensionality Reduction with PCA
Dimensionality reduction is a necessary step in scRNA-seq analysis, as this data is extremely high-dimensional. The reduction makes the data more manageable to work with and reveals the underlying biological structure [16]. Principal component analysis (PCA) is a dimensionality reduction method that captures linear variation in the data. PCA was performed on the scaled HVGs using RunPCA(), reducing the data to a lower-dimensional representation. After, an elbow plot was generated to determine how many principal components (PCs) captured meaningful variation. The goal is to use the fewest PCs that capture the most variation for downstream analysis. 

### 3.0 - Clustering and UMAP

#### 3.1 - FindNeighbors()
A shared nearest neighbour (SNN) graph was constructed using `FindNeighbors()` based on the first 40 principal components, as determined by visual inspection of the elbow plot. This graph forms the basis for graph-based clustering in the following step [16].

#### 3.2 - FindClusters() and UMAP
Clustering was performed using `FindClusters()` across a range of resolutions (0.1 to 0.8) to identify the optimal granularity of clusters [19]. To determine the most appropriate resolution, the clustree package was used to visualize cluster stability across resolutions, as outlined in a scRNA-seq clustering tutorial [19]. The clustree plot was evaluated at each resolution to identify the point at which clusters remained stable without producing an excessive number of random branches, indicating over-clustering. A resolution of 0.4 was selected based on this assessment, as it represented the point where cluster structure was stable and biologically interpretable. UMAP was then run using `RunUMAP()` on the first 40 PCs to produce the final visualization. 

#### 3.3 - Checking for Batch Effects
Batch effects were assessed by visualizing the UMAP coloured by `mouse_id`. UMAP plots were additionally generated, coloured by tissue type to examine whether clustering reflected tissue of origin or broader biological variation.

### 4.0 - Cell Type Annotation

#### 4.1 - FindAllMarkers()
To identify marker genes for each cluster, `FindAllMarkers()` was run using only positive markers, with a log2 fold-change threshold of 0.5, a minimum cell expression proportion of 0.2, and a maximum of 1,000 cells per identity, consistent with the approach used in the original study [3, 20]. The log2 fold-change threshold ensures that only genes meaningfully enriched in a cluster relative to all others are retained, as broadly expressed genes are not useful for cluster identification. The top 5 marker genes per cluster were extracted and inspected, prioritizing genes with high expression within the cluster of interest (pct.1) and low expression elsewhere (pct.2), as together these metrics capture how localized a gene's expression is to a specific cluster. For selected clusters, `FeaturePlot()` was used to visually confirm marker gene expression on the UMAP.

#### 4.2 - Manual Annotation
The top 4 genes per cluster were selected based on the pct.1 and pct.2 criteria described above and researched to determine their known biological functions and infer cluster identity. Clusters were then manually assigned cell type labels, with several sharing the same high-level identity merged for simplicity. Not all clusters were annotated leaving the remaining clusters to be grouped as "Other" on the UMAP. Once complete, annotations were cross-referenced against the Broad Institute Single Cell Portal, where the data from the original study has been publicly deposited [3], to confirm labelling decisions. The final annotations were saved to the metadata and visualized on the UMAP.

#### 4.3 - Automated Annotation
Automated cell type annotation was performed using SingleR to compare with the manual annotations [13]. The Seurat object was converted to a SingleCellExperiment object and annotated against two reference datasets: `MouseRNAseqData` for general mouse cell types [22] and `ImmGenData` for mouse immune cell types [23], with the use of both informed by the immune and non-immune populations identified during manual annotation. Parallelization was applied using BiocParallel with `SnowParam()` across 4 workers to manage the increased computational load, as all processes were run on a local machine [24]. The resulting labels were added to the metadata and a comparison between manual and SingleR annotations was produced to assess agreement between the two approaches.

### 5.0 - Differential Gene Expression (DGE) Analysis

#### 5.1 - Time Point and Cell Type Selection
DGE analysis was performed comparing cells at 5 days post-infection (D05) and 14 days post-infection (D14), representing the peak immune response and viral clearance stages. The analysis focused on six immune cell types: T cells, NK cells, Dendritic cells, Macrophages, Monocytes, and Neutrophils. These were chosen from analyzing the annotated UMAPs and identifying cell types that clustered together. Cells with unassigned mouse IDs were removed beforehand to ensure all observations could be attributed to a specific biological replicate.

#### 5.2 - Pseudobulk Aggregation
Gene expression counts were aggregated across cells within each combination of time point, mouse ID, and cell type using `AggregateExpression()`. Without this step, each individual cell would be treated as an independent observation, artificially inflating the sample size and leading to an elevated false positive rate. By reducing the data to one aggregated value per mouse per cell type per condition, the analysis is performed at the level of biological replicates rather than individual cells. Combined cell type and time point identities were then created to enable pairwise comparisons between D05 and D14 for each cell type.

#### 5.3 - DESeq2 Differential Expression
DGE between D05 and D14 was tested for each cell type using Seurat’s `FindMarkers()` function with the DESeq2 statistical test applied to pseudobulk-aggregated data. For each comparison, D14 was set as `ident.1` and D05 as `ident.2`, meaning positive log2 fold-change values reflect genes upregulated at D14 and negative values reflect genes upregulated at D05. Results were filtered using an adjusted p-value threshold of 0.05. Heatmaps were generated for each cell type showing the top 10 most upregulated and top 10 most downregulated genes, ranked by average log2 fold-change, with gene expression values scaled prior to plotting.

#### 5.4 - Feature Plot Validation
To further support the differential expression findings, `FeaturePlot()` was used to visualize select genes of interest across time all the time points. For Macrophages, _Isg15_ and _Mrc1_ were plotted as representative markers of D05 and D14 activity. For Neutrophils, _Phf11d_ and _Hk2_ were visualized across time points.

### 6.0 - Functional Enrichment Analysis

#### 6.1 - Over-Representation Analysis (ORA)
To identify the biological processes associated with the differentially expressed genes, over-representation analysis (ORA) was performed using the `enrichGO()` function from the clusterProfiler package [15]. For each cell type, genes were filtered to retain only those that were significantly upregulated at D14 relative to D05, using an adjusted p-value threshold of 0.05 and a log2 fold-change threshold of 0.5. This focused the analysis on genes associated with the transition from peak immune response to viral clearance.

Gene symbols were first converted to Entrez IDs using `bitr()` with the mouse annotation database `org.Mm.eg.db`, as this is required for GO enrichment testing. ORA was then performed against the Biological Process (BP) ontology, which captures the higher-level biological processes that genes are involved in. The top 10 enriched GO terms were visualized as dot plots for each cell type.

## Results
### Quality Control 
<div align="center">
<img src="figures/QC_violin_plot_1.png" width="600"/>
<br>
<b>Figure 1. Pre-filtered Violin plots for nFeature_RNA and nCount_RNA. </b> Violin plots showing the number of detected genes (nFeature_RNA, left) and total counts (nCount_RNA, right) per cell across five samples prior to filtering. Distributions are consistent across all samples, reflecting high baseline data quality.
</div>

<br>
Violin plots were used to evaluate the quality of the data in the provided Seurat object. Figure 1 demonstrates the distributions of nFeature_RNA and nCount_RNA, which are the number of detected genes and the total RNA counts per cell across the five time points. These violin plots show a consistent distributions across all samples. The filtering thresholds based on this plot were set to nFeature_RNA between 500 and 3,800 and nCount_RNA between 1,000 and 10,000. The nFeature_RNA threshold was chosen to remove low-quality cells while retaining higher-quality cells with greater gene complexity. The nCount_RNA threshold was chosen to preserve the majority of cells while still capturing the variation in the upper tails across the different time points.
</br>
<br>
<div align="center">
<img src="figures/QC_violin_plot_2.png" width="600"/>
<br>
<b>Figure 2. Pre-filtered Violin plots for mitochondrial, ribosomal, and hemoglobin gene percentages. </b> Violin plots showing the percentage of mitochondrial (percent.mt, left), ribosomal (percent.ribo, middle), and hemoglobin (percent.hb, right) gene counts per cell across five samples prior to filtering. Individual cells are shown for percent.hb to highlight the sparse distribution of hemoglobin-expressing cells.
</div>

<br>
Figure 2 shows the distributions of mitochondrial, ribosomal, and hemoglobin gene percentages across all five time points. The majority of cells show low mitochondrial and hemoglobin percentages, with only rare extreme outliers. To reduce potential noise, cells were filtered at thresholds of less than 15% mitochondrial, less than 50% ribosomal, and less than 5% hemoglobin reads [18]. Overall, these distributions indicate relatively clean data, suggesting that only minor filtering was required.
<br>
<br>
<div align="center">
<img src="figures/QC_knee_plot.png" width="400" height="400"/>
<br>
<b>Figure 3. Pre-filtering knee plot of cells ranked by total RNA counts. </b> Cells from all five time points (Naive, D02, D05, D08, D14) are ranked by total count (log10 nCount_RNA). The gradual and consistent decline across all samples indicates that the vast majority of cells are of high quality, supporting the use of minimal filtering thresholds
</div>

<br>
Figure 3 shows the ranking of cells by total RNA count across all time points. The gradual decline in the curve, rather than a sharp drop-off, indicates that most detected cells are high quality, with no clear inflection point separating low-quality droplets. This supports the filtering thresholds chosen from the violin plots and suggests that only minimal filtering was necessary.
</br>
<br>

<div align="center">
<img src="figures/QC_violin_plot_3.png" width=900" height="450"/>
<br>
<b>Figure 4. Post-filtered Violin plots for all QC conditions. </b> Violin plots showing the distribution of detected genes, total counts, mitochondrial, hemoglobin, and ribosomal gene percentages per cell across five sampling times following quality control filtering. Individual cells are shown for percent.hb to highlight the sparse distribution of hemoglobin-expressing cells.
</div>

</br>
Figure 4 shows the distributions of key quality control metrics across all five time points following filtering. Compared to the pre-filtered data, the distributions are tighter and more consistent, with reduced extreme values. Low mitochondrial and hemoglobin percentages are maintained across samples, indicating that the applied filtering thresholds effectively improved overall data quality.
<br>
<br>
<div align="center">
<b>Table 1.</b> Pre- and post-filtering cell counts per condition.

| Condition | Pre-filtering | Post-filtering | Cells Removed | % Removed |
|-----------|--------------|----------------|---------------|-----------|
| D02 | 40,148 | 40,002 | 146 | 0.36% |
| D05 | 26,344 | 26,235 | 109 | 0.41% |
| D08 | 29,523 | 29,414 | 109 | 0.37% |
| D14 | 30,203 | 30,054 | 149 | 0.49% |
| Naive | 30,354 | 30,261 | 93 | 0.31% |
| **Total** | **156,572** | **155,966** | **606** | **0.39%** |
</div>

<br>
Table 1 summarizes the filtering results. Only a small number of cells were removed during filtering (0.39% overall), with similar percentages across all conditions. This shows that the dataset was already of high quality, and filtering mainly removed a few low-quality cells without changing the overall data structure.

### Single-Cell Analysis 
<br>

<div align="center">
<img src="figures/elbow_plot.png" width=600" height="400"/>
<br>
<b>Figure 5. Elbow plot of first 50 PCs. </b> Standard deviation explained by each of the top 50 principal components (PCs). The curve gradually flattens beyond PC 40, suggesting that the majority of meaningful biological variation is captured within the first 40 PCs, which were selected for downstream clustering and dimensionality reduction.
</div>

<br>
Following normalization, HVG selection, and scaling, PCA was performed on the scaled HVGs to reduce the dimensionality of the data. Figure 5 displays the elbow plot of the first 50 principal components. The curve shows a gradual decline in standard deviation, beginning to flatten beyond PC 40, suggesting that the majority of meaningful biological variation is captured within the first 40 PCs. Based on this, the first 40 PCs were selected for downstream clustering and dimensionality reduction, which is consistent with the selection made by the original authors [3, 20].
<br>

<div align="center">
<img src="figures/pca_plot.png" width=600" height="500"/>
<br>
<b>Figure 6. PC1-PC2 plot of cells across all samples. </b> Scatter plot of cells projected onto the first two principal components (PC1: 11.2%, PC2: 9.8%), colored by sample identity. Cells from all five samples intermix across the principal component space.
</div>
</br>
Figure 6 shows the PCA projection of cells across all time points. Cells from different samples are well mixed across the PCA space, suggesting minimal batch effects. The visible structure likely reflects underlying biological variation rather than technical differences between samples.
<br>
</br>

<div align="center">
<img src="figures/umap_resolutions.png" width=600" height="500"/>
<br>
<b>Figure 7. Clustering resolution comparison for UMAPs. </b> UMAP visualizations of cells clustered at resolutions 0.1 through 0.8, yielding 22 to 53 clusters respectively. Increasing resolution results in progressive subdivision of cell populations. These plots were used to identify an optimal clustering resolution that balances biological interpretability with sufficient granularity for downstream cell type annotation.
</div>

<br>

<div align="center">
<img src="figures/cluster_tree.png" width=600" height="500"/>
<br>
<b>Figure 8. Clustree diagram showing cluster stability across resolutions. </b> Clustree plot displaying how clusters evolve across resolutions 0.1 to 0.8. The optimal resolution was selected at the point where clusters remain stable with few spurious branching events, balancing cluster granularity with interpretability.
</div>

<br>
UMAP is a dimensionality reduction method that captures non-linear patterns of variation, allowing for the visualization of complex cell population structure. Following the construction of the SNN graph using the first 40 PCs, clustering was performed and visualized using UMAP. Figure 7 displays UMAP plots across resolutions 0.1 to 0.8, showing an increasing number of clusters with higher resolutions. Figure 8 displays the clustree plot, where increasing random branching at higher resolutions indicated over-clustering, while resolution 0.4 showed stable clusters without any random branches. Based on this, a resolution of 0.4 was selected as the point where clusters were stable and biologically interpretable. The original authors used a higher resolution of 0.6, however this was supported by multiple rounds of filtering that produced a more refined dataset [3, 20]. Given that this analysis used a single filtering round, a more conservative resolution of 0.4 was chosen.
<br>
</br>

<div align="center">
<img src="figures/umap_batch.png" width=600" height="500"/>
<br>
<b>Figure 9. UMAP colored by mouse and tissue identity. </b> Clustree diagram showing cluster stability across resolutions. </b> UMAP visualization of all cells colored by individual mouse, tissue type (ET, RT, Sinus), and timepoint (Naive, D02, D05, D08, D14).
</div>

<br>
Figure 9 supports minimal batch effects, as cells from all three mice and conditions are well mixed throughout the UMAP rather than clustering by sample identity, justifying the decision to proceed without batch correction.
<br>
</br>

<div align="center">
<img src="figures/umap_cluster_tissue.png" width=700" height="500"/>
<br>
<b>Figure 10. UMAP visualization of cell clusters and tissue type distribution. </b> UMAP plots showing all cells colored by cluster identity (left, 38 clusters) and tissue type (right; LNG, OM, RM). 
</div>

<br>
Figure 10 displays the base UMAP with 38 clusters alongside the UMAP colored by tissue type (LNG, OM, RM). While some clusters show tissue-specific enrichment, the majority of clusters contain cells from multiple tissue types, suggesting that clustering is primarily driven by cell type rather than tissue of origin. This further supports that the clustering reflects meaningful biological variation rather than technical differences between tissues.
<br>
</br>

<div align="center">
<img src="figures/fp_0.png" width=700" height="500"/>
<br>
<b>Figure 11. Feature plots of top marker genes used to annotate cluster 0.</b> UMAP feature plots showing normalized expression of the four most significant marker genes for cluster 0 (<i>Ncam2</i>, <i>Gldc</i>, <i>Cidea</i>, and <i>Mdga2</i>). Expression is specifically enriched in a distinct region of the UMAP, supporting the annotation of this cluster as olfactory neurons. This approach was repeated for all major clusters to guide manual cell type annotation.
</div>

<br>
Once the optimal resolution of 0.4 was identified, cell type annotation was performed. Both manual and automated annotation were carried out to allow for a comparison between the two approaches. Manual annotation was performed by identifying the top marker genes for each cluster and researching their known biological functions to infer cell type identity. Figure 11 displays an example of this process for cluster 0, showing feature plots of the top four marker genes (<i>Ncam2</i>, <i>Gldc</i>, <i>Cidea</i>, and <i>Mdga2</i>), where expression is specifically enriched in a distinct region of the UMAP, supporting the annotation of this cluster as olfactory neurons. This approach was repeated across all major clusters.
<br>

<div align="center">
<img src="figures/umap_manual_annotations.png" width=700" height="500"/>
<br>
<b>Figure 12. UMAP of manually annotated cell types. </b> UMAP visualization of all cells following manual annotation, identifying 14 distinct cell populations including immune cells (neutrophils, macrophages, monocytes, dendritic cells, T cells, NK cells, B cells), structural cells (endothelial, fibroblasts, pericytes, basal cells), and nasal-specific cell types (olfactory neurons, sustenacular cells). Cell type identity was assigned based on known marker gene expression. Unresolved clusters are labeled as "Other."
</div>

<br>

Figure 12 displays the resulting manually annotated UMAP, where distinct cell populations form well-separated clusters, reflecting the biological diversity of the nasal mucosa. Both immune and structural cell types were successfully identified, along with nasal-specific populations that are characteristic of this tissue. The presence of these cell types is consistent with the three tissue types present in this dataset, being the RM, OM, and LNG.

<br>

<div align="center">
<img src="figures/umap_singleR.png" width=700" height="500"/>
<br>
<b>Figure 13.  UMAP of SingleR automated cell type annotation. </b> UMAP visualization of all cells annotated using SingleR with the mouse immune and general mouse reference datasets, identifying 27 predicted cell types.
</div>

<br>
Figure 13 displays the SingleR automated annotation, which predicted 27 cell types using the `MouseRNAseqData` and `ImmGenData` reference datasets. While several major populations were consistent with the manual annotation, SingleR predicted cell types not expected in nasal tissue such as hepatocytes and adipocytes, likely as a result of the broad MouseRNAseqData reference not accurately reflecting the tissue being studied. As a result, manual annotation was used to guide the final cell type assignments.
<br>

### DGE analysis

<br>

<div align="center">
<img src="figures/macrophages_DE_heatmap.png" width=600" height="600"/>
<br>
<b>Figure 14. Heatmap of differentially expressed genes in macrophages between DPI 5 and DPI 14.</b> Heatmap showing scaled expression of the top differentially expressed genes in macrophages at day 5 (teal) and day 14 (pink) viral clearance. Interferon-stimulated genes (<i>Ifit2</i>, <i>Rsad2</i>, <i>Isg15</i>, <i>Irf7</i>) are enriched at DPI 5, while heat shock and regulatory genes (<i>Hspa1a</i>, <i>Hspa1b</i>, <i>Hes1</i>) are more prominent at DPI 14.
</div>
<br>

<div align="center">
<img src="figures/neutrophils_DE_heatmap.png" width=600" height="600"/>
<br>
<b>Figure 15. Heatmap of differentially expressed genes in neutrophils between DPI 5 and DPI 14.</b> Heatmap showing scaled expression of the top differentially expressed genes in neutrophils at day 5 (teal) and day 14 (pink) viral clearance. Some interferon-related genes (<i>Ifi211</i>, <i>Bst2</i>) show higher expression at DPI 5, while metabolic genes (<i>Ndufv3</i>, <i>Ndufa3</i>) are modestly elevated at DPI 14, though the overall transcriptional differences between timepoints are less pronounced than those observed in macrophages.
</div>
<br>

<div align="center">
<img src="figures/fp_I_H_comparison.png" width=800" height="400"/>
<br>
<b>Figure 16. Feature plots of <i>Isg15</i> and <i>Hes1</i> expression across timepoints in macrophages.</b> UMAP feature plots showing expression of <i>Isg15</i> (top) and <i>Hes1</i> (bottom) across all five timepoints (D02, D05, D08, D14, Naive). <i>Isg15</i> expression is notably elevated at D05, consistent with the myeloid recruitment stage of infection, and declines by D14 following viral clearance. <i>Hes1</i> expression appears slightly lighter at D05 and becomes more broadly and strongly expressed by D14, consistent with its role as a regulatory factor during the resolution of infection.
</div>

<br>

<div align="center">
<img src="figures/fp_neutro_comp.png" width=800" height="400"/>
<br>
<b>Figure 17. Feature plots of <i>Phf11d</i> and <i>Hk2</i> expression across timepoints in neutrophils.</b> UMAP feature plots showing expression of <i>Phf11d</i> (top) and <i>Hk2</i> (bottom) across all five timepoints (D02, D05, D08, D14, Naive). <i>Phf11d</i> expression appears elevated at D05 during the myeloid recruitment stage and declines by D14 following viral clearance. <i>Hk2</i> expression appears highest at D02 during peak viral load and becomes more broadly distributed across timepoints, reflecting the more modest transcriptional differences observed in neutrophils compared to macrophages.
</div>

<br>

Following cell type annotation, differential gene expression analysis was performed comparing cells at D05 and D14. Gene expression counts were first aggregated across cells within each combination of time point, mouse ID, and cell type using `AggregateExpression()`, allowing the analysis to be performed at the level of biological replicates rather than individual cells. DGE was then tested using Seurat’s `FindMarkers()` function with the DESeq2 statistical test applied to pseudobulk-aggregated data, where positive log2 fold-change values reflect genes upregulated at D14 and negative values reflect genes upregulated at D05, using an adjusted p-value threshold of 0.05.

Figures 14 and 15 display the heatmaps of the top differentially expressed genes in macrophages and neutrophils between D05 and D14, respectively. In both cell types, interferon-related genes show higher expression at D05, while distinct gene sets are elevated at D14. Transcriptional differences between time points appear more pronounced in macrophages than in neutrophils.

Figure 16 displays the feature plots of Isg15 and Hes1 in macrophages across all five time points. _Isg15_ expression is notably elevated at D05 and declines by D14, while _Hes1_ expression appears slightly higher at D14 compared to D05. Figure 17 displays the feature plots of _Phf11d_ and _Hk2_ in neutrophils across all five time points. _Phf11d_ expression appears elevated at D05 and declines by D14, while _Hk2_ expression is higher at D14, particularly within the neutrophil cluster region.

### Over-representation analysis (ORA)
<br>

<div align="center">
<img src="figures/macrophages_final_GO.png" width=600" height="700"/>
<br>
<b>Figure 18.  Gene Ontology (GO) ORA of macrophage pathways at DPI 14 versus DPI 5. </b> Dot plot showing GO biological pathways significantly enriched in macrophages at DPI 14 compared to DPI 5. Dot size represents the number of genes in each pathway and color indicates statistical significance.
</div>

<br>

Figure 18 shows the GO ORA dot plot for macrophages at D14 compared to D05. The most enriched pathways are related to tumor necrosis factor (TNF) response, which have the highest gene ratios and lowest adjusted p-values. Other enriched pathways include those involved in handling protein stress, such as the unfolded protein response, as well as pathways related to reducing apoptotic (cell death) signaling.

## Discussion 
### Macrophages across viral infection stages
As shown in the infection timeline, D05 represents the myeloid recruitment stage and D14 represents viral clearance, providing a window to examine how macrophage transcriptional programs shift across these two distinct phases of infection [3]. Macrophages showed the most pronounced transcriptional differences between D05 and D14 and were highlighted due to the pronounced differences observed between the two time points. Additional heatmaps were generated for other cell type clusters and can be found in the [`figures`](figures/) directory. At D05, interferon-stimulated genes including _Isg15_, _Ifit2_, and _Ifit3_ were highly expressed, reflecting an active antiviral response consistent with the myeloid recruitment stage of infection. Myeloid cells are crucial for detecting pathogens and initiating the production of cytokines and interferons to mitigate infections [25]. Interferon-stimulated genes are a group of gene products that coordinately combat viral infections [26]. When a pathogen has been detected, interferons are rapidly produced and bind to receptors on nearby cells, activating the Janus kinase/signal transducer and activator of transcription (JAK-STAT) pathway, which then drives the rapid transcription of ISGs [26]. These ISGs then carry out a variety of antiviral functions such as degrading viral RNA, blocking viral replication, and preventing further viral spread [27]. This coordinated interferon response at D05 is consistent with the known role of innate immune cells in the early antiviral defense of the nasal mucosa [3].

This response is reflected in the heatmap of Figure 14 and the feature plot of Figure 16. These plots show _Isg15_ expression is visibly elevated and widespread at D05, consistent with the myeloid recruitment stage, and is clearly diminished by D14, consistent with viral clearance. _Isg15_ is recognized as one of the most highly upregulated proteins during active viral infection, functioning both as an intracellular antiviral effector and an extracellular signaling molecule that activates immune cells to promote viral clearance [28]. Its decline by D14 is therefore consistent with the resolution of the acute antiviral response, as the signaling that drives its expression is likely reduced [28].

By D14, representing viral clearance, the transcriptional profile of macrophages suggests a shift away from active antiviral signaling towards a resolution state, consistent with the known regulatory and immunosuppressive roles that myeloid cells can play following the elimination of viruses [25]. Heat shock proteins _Hspa1a_ and _Hspa1b_ were more prominently expressed at D14 as seen in the heatmap of Figure 14, which are known to be involved in regulating inflammation [29]. Their upregulation at D14 is consistent with macrophages transitioning towards a resolution phase following viral clearance [29]. The upregulation of _Hes1_ at D14 further supports this shift. _Hes1_ is a transcription repressor that suppresses inflammatory gene expression in macrophages, suggesting that the immune response is being actively dampened as the infection resolves by D14 [30]. This is also evident in the feature plot of Figure 16, where _Hes1_ expression appears slightly lighter at D05 than at D14, becoming slightly more expressed across the UMAP, consistent with its role as a sustained regulatory factor that increases during the resolution of infection [30].

### Neutrophils across viral infection stages
In contrast to macrophages, neutrophils showed less pronounced transcriptional differences between D05 and D14. Neutrophils are the most abundant innate immune cell and are rapidly recruited to the site of infection [25]. Upon recruitment, they can directly engulf and destroy influenza A virus particles, and release toxic substances from their granules that help neutralize the virus [31]. This is consistent with the elevated expression of interferon-related genes, like _Ifi211_, at D05 observed in Figure 15, reflecting an ongoing antiviral response during the myeloid recruitment stage of infection. However, the overall transcriptional differences between D05 and D14 in neutrophils are more modest compared to macrophages, which may reflect the fact that neutrophils are short-lived cells whose peak activity likely happens before D05 [25, 32]. As seen in the feature plot of Figure 17, _Hk2_ expression appears highest at D02 during peak viral load, which is consistent with previous findings that influenza A virus infection is associated with increased _Hk2_ expression and glycolytic activity [33, 34]. By D05 and into D14, _Hk2_ expression becomes more broadly distributed across time points, suggesting a shift in metabolic activity as infection progresses and the virus is cleared [33]. Overall, the neutrophil transcriptional response reflects their importance as an early responder during influenza infection [25, 31, 32].

### Macrophage GO BP enrichment
Figure 18 displays the GO ORA dot plot for macrophages at D14 compared to D05, showing which biological pathways are most enriched during the resolution phase of infection.

The two most significantly enriched pathways are the cellular response to tumor necrosis factor (TNF) and response to tumor necrosis factor, which show the largest dot sizes and lowest adjusted p-values in the plot. TNF is a signaling molecule that is produced by macrophages in response to viral infection and plays an important role in regulating cell survival, proliferation, and immune responses [35, 36]. TNF has also been shown to inhibit the replication of viruses such as influenza, suggesting that its enrichment at D14 may reflect an ongoing antiviral and regulatory role in macrophages during the resolution of infection [35]. Given that macrophages are major producers and responders to TNFα [36], the enrichment of TNF response pathways at D14 suggests continued involvement of TNF-related signaling pathways during this stage [35, 36].

Additional enriched pathways are related to the unfolded protein response, which is a cellular mechanism for managing protein misfolding and stress during and after infection [37]. This response has been associated with macrophage survival and sustained immune function, and its enrichment at D14 suggests that macrophages may be adapting to cellular stress as they transition out of an active infection state [37].

Finally, the enrichment of pathways involved in the negative regulation of apoptosis suggests that macrophages may be promoting their own survival during the resolution phase. The enrichment of TNF-related pathways also indicates that TNF signaling is still active, and since TNF can influence both cell death and survival, it may help macrophages persist at this stage. Overall, this suggests that TNF signaling may contribute to macrophage survival and their role in tissue repair and restoring normal function after viral clearance [35, 38].

## References
[1] N. Zhang, K. Van Crombruggen, E. Gevaert, and C. Bachert, “Barrier function of the nasal mucosa in health and type-2 biased airway diseases,” Allergy, vol. 71, no. 3, pp. 295–307, Jan. 2016, doi: https://doi.org/10.1111/all.12809. <br/>
[2] F. Zhu et al., “H1N1 Influenza Virus-Infected Nasal Mucosal Epithelial Progenitor Cells Promote Dendritic Cell Recruitment and Maturation,” Frontiers in Immunology, vol. 13, no. 879575, Apr. 2022, doi: https://doi.org/10.3389/fimmu.2022.879575.  <br/>
[3] S. W. Kazer et al., “Primary nasal influenza infection rewires tissue-scale memory response dynamics,” Immunity, vol. 57, no. 8, pp. 1955-1974.e8, Aug. 2024, doi: https://doi.org/10.1016/j.immuni.2024.06.005.  <br/>
[4] C. Contreras, J. M. Newby, and T. Hillen, “Personalized Virus Load Curves for Acute Viral Infections,” Viruses, vol. 13, no. 9, p. 1815, Sep. 2021, doi: https://doi.org/10.3390/v13091815.  <br/>
[5] A. Iwasaki and P. S. Pillai, “Innate immunity to influenza virus infection,” Nature Reviews Immunology, vol. 14, no. 5, pp. 315–328, Apr. 2014, doi: https://doi.org/10.1038/nri3665.  <br/>
[6] F. A. Wolf, P. Angerer, and F. J. Theis, “SCANPY: large-scale single-cell gene expression data analysis,” Genome Biology, vol. 19, no. 1, Feb. 2018, doi: https://doi.org/10.1186/s13059-017-1382-0.  <br/>
[7] R. A. Amezquita et al., “Orchestrating single-cell analysis with Bioconductor,” Nature Methods, vol. 17, no. 2, Dec. 2019, doi: https://doi.org/10.1038/s41592-019-0654-x.  <br/>
[8] I. Billato et al., “Benchmarking large-scale single-cell RNA-seq analysis,” Oct. 2025, doi: https://doi.org/10.1101/2025.10.28.681564.  <br/>
[9] R. Satija, “Seurat,” Satijalab.org, 2019. https://satijalab.org/seurat/  <br/>
[10] D. Aran et al., “Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage,” Nature Immunology, vol. 20, no. 2, pp. 163–172, Feb. 2019, doi: https://doi.org/10.1038/s41590-018-0276-y.   <br/>
[11] V. Y. Kiselev, A. Yiu, and M. Hemberg, “scmap: projection of single-cell RNA-seq data across data sets,” Nature Methods, vol. 15, no. 5, pp. 359–362, Apr. 2018, doi: https://doi.org/10.1038/nmeth.4644.  <br/>
[12] J. K. de Kanter, P. Lijnzaad, T. Candelli, T. Margaritis, and F. C. P. Holstege, “CHETAH: a selective, hierarchical cell type identification method for single-cell RNA sequencing,” Nucleic Acids Research, vol. 47, no. 16, Jun. 2019, doi: https://doi.org/10.1093/nar/gkz543.  <br/> 
[13] Q. Huang, Y. Liu, Y. Du, and L. X. Garmire, “Evaluation of Cell Type Annotation R Packages on Single-cell RNA-seq Data,” Genomics, Proteomics & Bioinformatics, vol. 19, no. 2, Dec. 2020, doi: https://doi.org/10.1016/j.gpb.2020.07.004.  <br/>
[14] A. L. Thurman, J. A. Ratcliff, M. S. Chimenti, and A. A. Pezzulo, “Differential gene expression analysis for multi-subject single-cell RNA-sequencing studies with aggregateBioVar,” Bioinformatics, vol. 37, no. 19, pp. 3243–3251, Apr. 2021, doi: https://doi.org/10.1093/bioinformatics/btab337.  <br/>
[15] T. Wu et al., “clusterProfiler 4.0: A universal enrichment tool for interpreting omics data,” The Innovation, vol. 2, no. 3, p. 100141, Aug. 2021. <br/>
[16] “Seurat - Guided Clustering Tutorial,” satijalab.org, Oct. 31, 2023. https://satijalab.org/seurat/articles/pbmc3k_tutorial.html <br/>
[17] “Plot the Barcode Distribution and Calculated Inflection Points — BarcodeInflectionsPlot,” Satijalab.org, 2026. https://satijalab.org/seurat/reference/barcodeinflectionsplot (accessed Apr. 8, 2026). <br/>
[18] M. Su et al., “Data analysis guidelines for single-cell RNA-seq in biomedical studies and clinical applications,” Military Medical Research, vol. 9, no. 1, Dec. 2022, doi: https://doi.org/10.1186/s40779-022-00434-8. <br/>
[19] “Chapter 4 Clustering | scRNAseq Analysis in R with Seurat,” Github.io, 2025. https://swbioinf.github.io/scRNAseqInR_Doco/clustering.html (accessed Apr. 8, 2026).  <br/>
[20] jo-m-lab, “GitHub - jo-m-lab/IAV-nasal-sc-atlas at v1.0.0,” GitHub, May 02, 2024. https://github.com/jo-m-lab/IAV-nasal-sc-atlas/tree/v1.0.0 (accessed Apr. 8, 2026).  <br/>
[21] “Chapter 5 Using multiple references | Assigning cell types with SingleR,” Bioconductor.org, 2025. https://bioconductor.org/books/release/SingleRBook/using-multiple-references.html (accessed Apr. 8, 2026). <br/>
[22] LTLA, “Obtain mouse bulk expression data of sorted cell populations...,” Rdrr.io, Jun. 03, 2024. https://rdrr.io/github/LTLA/celldex/man/MouseRNAseqData.html (accessed Apr. 8, 2026).  <br/>
[23] LTLA, “Obtain mouse bulk expression data from the Immunologic Genome...,” Rdrr.io, Jun. 03, 2024. https://rdrr.io/github/LTLA/celldex/man/ImmGenData.html (accessed Apr. 8, 2026).  <br/>
[24] “1. Introduction to BiocParallel,” Bioconductor.org, 2022. https://www.bioconductor.org/packages//release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.html#clusters-of-independent-processes-with-snowparam (accessed Apr. 8, 2026). <br/>
[25] W. Wang, L. Xu, J. Su, M. P. Peppelenbosch, and Q. Pan, “Transcriptional Regulation of Antiviral Interferon-Stimulated Genes,” Trends in Microbiology, vol. 25, no. 7, pp. 573–584, Jul. 2017, doi: https://doi.org/10.1016/j.tim.2017.01.001.  <br/>
[26] J. W. Schoggins and C. M. Rice, “Interferon-stimulated genes and their antiviral effector functions,” Current Opinion in Virology, vol. 1, no. 6, pp. 519–525, Dec. 2011, doi: https://doi.org/10.1016/j.coviro.2011.10.008.  <br/>
[27] B. T. Freitas, F. E. M. Scholte, É. Bergeron, and S. D. Pegan, “How ISG15 combats viral infection,” Virus Research, vol. 286, p. 198036, Sep. 2020, doi: https://doi.org/10.1016/j.virusres.2020.198036.  <br/>
[28] Helena Trevisan Schroeder, C. Henrique, Thiago Gomes Heck, M. Krause, and Paulo, “Heat shock response during the resolution of inflammation and its progressive suppression in chronic-degenerative inflammatory diseases,” Cell Stress & Chaperones, vol. 29, no. 1, Jan. 2024, doi: https://doi.org/10.1016/j.cstres.2024.01.002. <br/>
[29] Y. Shang et al., “The transcriptional repressor Hes1 attenuates inflammation by regulating transcription elongation,” Nature Immunology, vol. 17, no. 8, pp. 930–937, Jun. 2016, doi: https://doi.org/10.1038/ni.3486. <br/>
[30] A. A. Stegelmeier et al., “Myeloid Cells during Viral Infections and Inflammation,” Viruses, vol. 11, no. 2, Feb. 2019, doi: https://doi.org/10.3390/v11020168. <br/>
[31] C. Johansson and F. C. M. Kirsebom, “Neutrophils in respiratory viral infections,” Mucosal Immunology, vol. 14, no. 4, pp. 815–827, Jul. 2021, doi: https://doi.org/10.1038/s41385-021-00397-4.  <br/>
[32] D. J. Eddins et al., “Transcriptional reprogramming of infiltrating neutrophils drives lung pathology in severe COVID-19 despite low viral load,” Blood Advances, vol. 7, no. 5, pp. 778–799, Nov. 2022, doi: https://doi.org/10.1182/bloodadvances.2022008834.  <br/>
[33] Z. Bie, Y. Tan, C. Ye, and K. Wei, “Metabolic hijacking styles: a review of how viral life cycles dictate glucose metabolism reprogramming,” Frontiers in Microbiology, vol. 16, Dec. 2025, doi: https://doi.org/10.3389/fmicb.2025.1690133.  <br/>
[34] J. Li, Y. Wang, H. Deng, L. Su, and H. Qiu, “Cellular metabolism hijacked by viruses for immunoevasion: potential antiviral targets,” Frontiers in Immunology, vol. 14, Jul. 2023, doi: https://doi.org/10.3389/fimmu.2023.1228811. <br/>
[35] K. Pant, A. Chandrasekaran, C. J. Chang, A. Vageesh, A. J. Popkov, and J. B. Weinberg, “Effects of tumor necrosis factor on viral replication and pulmonary inflammation during acute mouse adenovirus type 1 respiratory infection,” Virology, vol. 547, pp. 12–19, Aug. 2020, doi: https://doi.org/10.1016/j.virol.2020.05.004. <br/>
[36] N. Parameswaran and S. Patial, “Tumor Necrosis Factor-α Signaling in Macrophages,” Critical ReviewsTM in Eukaryotic Gene Expression, vol. 20, no. 2, pp. 87–103, 2010, doi: https://doi.org/10.1615/critreveukargeneexpr.v20.i2.10. <br/>
[37] J. Grootjans, A. Kaser, R. J. Kaufman, and R. S. Blumberg, “The unfolded protein response in immunity and inflammation,” Nature Reviews Immunology, vol. 16, no. 8, pp. 469–484, Jun. 2016, doi: https://doi.org/10.1038/nri.2016.62.  <br/>
[38] K. Högner et al., “Macrophage-expressed IFN-β Contributes to Apoptotic Alveolar Epithelial Cell Injury in Severe Influenza Virus Pneumonia,” PLoS Pathogens, vol. 9, no. 2, p. e1003188, Feb. 2013, doi: https://doi.org/10.1371/journal.ppat.1003188. <br/>
