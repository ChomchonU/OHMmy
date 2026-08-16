# Comprehensive single-cell workflow with OHMmy

## Introduction

Welcome to `OHMmy`. This package is designed to streamline highly
complex single-cell RNA sequencing workflows. From ambient RNA
decontamination and technical noise regression to differential
expression and pathway enrichment, `OHMmy` provides robust, automated
wrappers around industry-standard tools (like Seurat, SoupX, Nebulosa,
and DESeq2).

This comprehensive guide walks through a complete analysis pipeline.

``` r


library(Seurat)
library(OHMmy)

# Define a global output directory for this pipeline
out_dir <- "results/comprehensive_workflow/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
```

## 1. Ambient RNA Decontamination (SoupX)

Cell-free ambient RNA can severely confound downstream clustering.
`OHMmy` provides seamless wrappers to load standard outputs (e.g., from
CellRanger) and accurately remove this contamination before building our
Seurat object.

``` r


# 1. Process SoupX across all samples in your directories
soup_results <- process_soupx_samples(
  filtered_dir = "data/filtered/", 
  raw_dir = "data/raw/", 
  multiFac = 1 
)

# 2. Transfer the cleaned counts and metadata back to your raw Seurat object
seurat_raw <- readRDS("data/seurat_raw.rds")
seurat_clean <- addSoupXMetaToSeurat(
  original_seurat = seurat_raw,
  soupx_seurat_list = soup_results$final_seurat
)
```

## 2. Processing, Integration, and Clustering

Next, we normalize the data, regress out technical noise (like cell
cycle or mitochondrial percentages), integrate across batches, and
systematically determine the best clustering resolution.

``` r


# 1. Log-Normalization, technical regression, and Batch Integration
seurat_integrated <- ProcessSeuratLOG(
  seurat_clean,
  batch_col = "batch",
  vars_to_regress = c("pct_counts_mt", "S.Score", "G2M.Score"),
  tcr_bcr_patterns = "^TR[ABDG]|^IG[HKL]", # Exclude V(D)J genes from clustering
  integration_method = "CCAIntegration",
  reduction_name = "pca.log",
  integration_reduction = "integrated.cca.log",
  dims = 1:40
) 

# 2. Run UMAP and test clustering across multiple resolutions
cluster_results <- ClusterAndUMAP(
  seurat_obj = seurat_integrated,
  reduction = "integrated.cca.log",
  dims = 1:40,
  cluster_resolutions = seq(0.1, 2, by = 0.1),
  final_resolution = 1.0, 
  sample_name = "Cohort_A",
  plot_dir = paste0(out_dir, "clustree/")
)

# Extract the fully processed object
seurat_obj <- cluster_results$seurat
```

## 3. Quality Control & Principal Component Diagnostics

Before annotating, we verify that our clusters are biologically driven
and not artifacts of technical noise (e.g., sequencing depth or
ribosomal genes).

``` r


# 1. Check QC distributions across your new clusters
plot_violin_qc_single(
  seurat_obj = seurat_obj,
  res_col = "seurat_clusters",
  qc_meta_features = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt"),
  sample_name = "Cohort_A",
  output_dir = paste0(out_dir, "QC/")
)

# 2. Quantify exactly how much of each PC is driven by technical gene families
plot_technical_contribution(
  seurat_obj = seurat_obj,
  sample_name = "Cohort_A",
  technical_keywords = c("^MT-", "^RPL", "^RPS"),
  max_pcs = 40,
  cutoff = 15,
  output_dir = paste0(out_dir, "Diagnostics/")
)
```

## 4. Advanced Visualization and Annotation

Interpret the biological identity of your clusters using flexible UMAP
overlays, density mapping, and gene co-expression tools.

``` r


# 1. Visualize clusters and metadata overlays on the UMAP
PlotDimByFactors(
  seurat_obj = seurat_obj, 
  factors = c("seurat_clusters", "predicted.celltype.l2"),
  sample_name = "Cohort_A",
  output_dir = paste0(out_dir, "UMAP/")
)

# 2. Advanced Dotplots (Grouped with Dendrograms)
feature_df <- data.frame(
  feature = c("CD8A", "NCAM1", "GZMB", "PRF1", "CCR7", "IL7R"),
  group = c("Lineage", "Lineage", "Effector", "Effector", "Naive", "Naive")
)

plot_dot_dendro_split(
  seurat_obj = seurat_obj, 
  meta_col = "seurat_clusters",
  feature_df = feature_df,
  prefix = "Marker_Panel",
  output_dir = paste0(out_dir, "Dotplots/")
)

# 3. Density Blend Overlays (Nebulosa) for Co-Expression
plot_blend_nebulosa(
  seurat_obj = seurat_obj, 
  cluster_col = "seurat_clusters",
  gene_pairs = list(c("CX3CR1", "GZMB")),
  reduction = "umap",
  sample_name = "Cohort_A",
  output_dir = paste0(out_dir, "Blends/")
)
```

## 5. Statistical Abundance & Confounder Analysis

Ensure your cohorts are clinically balanced, then test for significant
shifts in cell type proportions between your conditions (e.g., Healthy
vs. Infected).

``` r


# 1. Confounder Check: Ensure patient metadata is balanced
plot_metadata_stats(
  seurat_obj = seurat_obj,
  sample_col = "donor_id",       
  condition_col = "disease_state", 
  metadata_vars = c("age", "viral_load", "sex"), 
  output_dir = paste0(out_dir, "Stats/")
)

# 2. Cell Type Abundance: Test for expansion/contraction of clusters
plot_cell_abundance(
  seurat_obj = seurat_obj,
  sample_col = "donor_id",       
  condition_col = "disease_state", 
  celltype_col = "seurat_clusters", 
  facet_by_cluster = TRUE,
  pairwise_label = "p.adj",
  output_dir = paste0(out_dir, "Abundance/")
)
```

## 6. Differential Expression & Functional Enrichment

Extract biological meaning by discovering defining cluster markers,
running Gene Set Enrichment Analysis (GSEA), and visualizing
condition-specific DESeq2 results.

``` r


library(msigdbr)

# 1. Find all cluster markers and automatically generate a heatmap
marker_results <- FindTopMarkersAndHeatmap(
  seurat_obj = seurat_obj, 
  sample_name = "Cohort_A",
  output_dir_base = paste0(out_dir, "Markers/")
)

# Extract the marker dataframe for downstream GSEA/ORA
my_markers_df <- marker_results$markers

# 2. PATHWAY ENRICHMENT (GSEA)
m_t2g_hallmark <- msigdbr(species = "Homo sapiens", category = "H")[, c("gs_name", "gene_symbol")]

gsea_results <- run_global_gsea(
  GSEA_df = my_markers_df, 
  m_t2g = m_t2g_hallmark, 
  title_prefix = "Hallmark",
  output_dir = paste0(out_dir, "GSEA/")
)

# 3. PSEUDO-BULK DIFFERENTIAL EXPRESSION VISUALIZATIONS (DESeq2)
# Note: These functions visualize pre-computed DESeq2 results objects (res_obj)

# Volcano Plots
generate_volcano_trio(
  res_obj = deseq2_disease_vs_healthy, 
  plot_title = "Severe Disease vs Healthy Controls", 
  target_genes = c("IFNG", "TNF", "GZMB")
)

# DESeq2 Heatmaps across samples
generate_and_save_heatmap(
  res_obj = deseq2_disease_vs_healthy,
  comp_name = "Disease_vs_Healthy",
  comp_title = "Severe Disease vs Healthy Controls",
  vsd_data = vsd_matrix,               
  anno_col = annotation_dataframe, 
  ordered_samps = sample_order_list,
  out_dir = paste0(out_dir, "DESeq2/")
)
```

**End of Workflow.** By following these steps, you have successfully
decontaminated, processed, integrated, evaluated, annotated, and
statistically analyzed your scRNA-seq dataset!
