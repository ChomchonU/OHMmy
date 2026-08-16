# Advanced Visualization and Marker Analysis

## Introduction

Once clustering is complete, interpreting the biological identity of
those clusters requires flexible, high-quality visualizations. `OHMmy`
provides robust wrappers around popular plotting libraries (like
`ggplot2` and `Nebulosa`) to automate the generation of UMAP overlays,
Dot Plots, Violin Plots, and Co-expression overlays.

This guide walks through the major visualization functions using a
processed Seurat object.

``` r


library(Seurat)
library(OHMmy)

# Load your clustered Seurat object
seurat_obj <- readRDS("data/seurat_clustered.rds")
out_dir <- "results/visualizations/"
```

## Part 1: Dimensionality Reduction & Metadata Overlays

To quickly visualize how different metadata categories (like predicted
cell types, patient phases, or clusters) map onto your UMAP, use
[`PlotDimByFactors()`](https://chomchonu.github.io/OHMmy/reference/PlotDimByFactors.md).
This function automatically handles formatting and saving.

``` r


# Visualize clusters and cell type predictions on the UMAP
PlotDimByFactors(
  seurat_obj = seurat_obj,
  factors = c("seurat_clusters", "predicted.celltype.l2"),
  reduction = "umap",
  sample_name = "Cohort_A",
  output_dir = paste0(out_dir, "umap_overlays/")
)
```

## Part 2: Feature Density Mapping (`plot_combined`)

The
[`plot_combined()`](https://chomchonu.github.io/OHMmy/reference/plot_combined.md)
function automatically loops through custom lists—whether they are gene
expression signatures or continuous metadata features (like sequencing
depth or clinical metrics)—and maps their gradients across your UMAP
space.

``` r


# 1. Define lists containing marker genes or numerical QC/clinical metadata
feature_subsets <- list(
  Markers = c("CD8A", "GZMB", "NCAM1", "FCGR3A"),
  QC_Clinical = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt", "log10_viremia")
)

# 2. Generate combined feature density plots on UMAP
plot_combined(
  seurat_obj = seurat_obj,
  cluster_col = "seurat_clusters",
  sample_name = "Cohort_A",
  reduction = "umap",
  genes_list = feature_subsets,
  output_dir = paste0(out_dir, "feature_density/"),
  save_format = "jpg",
  width = 20,
  dpi = 300
)
```

## Part 3: Profiling QC and Marker Distributions

### 3.1 Metadata & QC Violin Plots

To inspect how numerical QC features, cell cycle scores, or clinical
metadata variables (e.g., `age`, `log10_viremia`) are distributed across
your clusters, use
[`plot_violin_qc_single()`](https://chomchonu.github.io/OHMmy/reference/plot_violin_qc_single.md)
by passing them into `qc_meta_features`.

``` r


qc_metrics <- c("nCount_RNA", "nFeature_RNA", "pct_counts_mt", "age", "log10_viremia", "S.Score", "G2M.Score")

plot_violin_qc_single(
  seurat_obj = seurat_obj,
  sample_name = "Cohort_A",
  qc_meta_features = qc_metrics,
  res_col = "seurat_clusters",
  assay = "RNA",
  output_dir = paste0(out_dir, "violin_qc/"),
  save_format = "jpg",
  width = 20,
  dpi = 300
)
```

### 3.2 Marker Gene Violin Plots

You can similarly pass a list of specific gene features to evaluate
marker expression across clusters.

``` r


marker_genes <- c("CD8A", "GZMB", "NCAM1", "CCR7", "IL7R", "MKI67")

plot_violin_qc_single(
  seurat_obj = seurat_obj,
  sample_name = "Cohort_A",
  gene_features = marker_genes,
  res_col = "seurat_clusters",
  assay = "RNA",
  output_dir = paste0(out_dir, "violin_genes/")
)
```

### 3.3 Split Dot Plots with Dendrograms

For a macroscopic view of gene expression across subsets, the
[`plot_dot_dendro_split()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_split.md)
function creates beautiful Dot Plots grouped by feature categories.

``` r


# 1. Create a data frame assigning genes to functional groups
feature_df <- data.frame(
  feature = c("CD8A", "CD8B", "NCAM1", "GZMB", "GZMK", "CCR7", "IL7R"),
  group = c("Lineage", "Lineage", "Lineage", "Effector", "Effector", "Naive", "Naive")
)

# 2. Generate the grouped Dot Plot
plot_dot_dendro_split(
  seurat_obj = seurat_obj, 
  meta_col = "seurat_clusters", 
  feature_df = feature_df, 
  scale = TRUE,
  prefix = "Marker_Panel",
  output_dir = paste0(out_dir, "dot_plots/")
)
```

## Part 4: Gene Co-Expression and Blend Plots

Understanding how two genes interact (or are mutually exclusive) in the
same cell requires specialized density and correlation plots.

### 4.1 Statistical Correlation

Calculate the Pearson or Spearman correlation between specific gene
pairs within your clusters.

``` r


# Calculate expression correlation between key functional pairs
plot_gene_pair_correlations(
  seurat_obj = seurat_obj,
  cluster_col = "seurat_clusters",
  gene_pairs = list(
    c("CD8A", "GZMB"),
    c("IL7R", "TCF7"),
    c("IFNG", "TNF")
  ),
  cor_method = "spearman",
  sample_name = "Marker_Correlations",
  output_dir = paste0(out_dir, "correlations/")
)
```

### 4.2 Nebulosa Blend Overlays

To visually identify dual-expressing cells (e.g., cells expressing both
`CD8A` and a specific receptor), `OHMmy` uses `Nebulosa` to generate
smoothed density blend plots on your UMAP.

``` r


# Generate density blend plots for specific receptor/effector pairs
plot_blend_nebulosa(
  seurat_obj = seurat_obj,
  cluster_col = "seurat_clusters",
  gene_pairs = list(
    c("CX3CR1", "GZMB"),
    c("PDCD1", "TIGIT")
  ),
  reduction = "umap",
  sample_name = "Dual_Expression",
  blend_threshold = 0.1,
  output_dir = paste0(out_dir, "blend_plots/")
)
```

## Part 5: Clinical and Technical Demographics

Finally, after defining your clusters, you often need to see how
clinical metadata (like disease severity or batch) is distributed across
them.

``` r


# Plot the distribution of patient metadata across your defined clusters
metadata_to_check <- c("batch", "severity", "Phase")

for (meta in metadata_to_check) {
  plot_cluster_distributions(
    seurat_obj = seurat_obj, 
    cluster_col = "seurat_clusters",
    batch_col = meta,
    file_prefix = paste0("Dist_", meta),
    output_dir = paste0(out_dir, "distributions/")
  )
}
```
