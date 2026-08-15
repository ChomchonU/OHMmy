## ----include = FALSE----------------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE,
  eval = FALSE # Set to FALSE for fast package building
)

## -----------------------------------------------------------------------------
# 
# library(Seurat)
# library(OHMmy)
# 
# # 1. Load data and extract metadata
# counts_list <- load_counts(data_dir = "path/to/data")
# sample_names <- get_sample_names(counts_list)
# 
# # 2. Visualize initial Quality Control metrics
# plot_violin_qc_single(counts_list[[1]], feature = "nFeature_RNA")

## -----------------------------------------------------------------------------
# 
# # Step-by-step SoupX preparation
# soupx_inputs <- prepare_soupx_inputs(raw_counts = counts_list$raw, filtered_counts = counts_list$filtered)
# soup_channels <- create_soup_channels(soupx_inputs$raw, soupx_inputs$filtered)
# 
# # Pre-clustering for contamination estimation
# seurat_pre <- create_seurat_for_clustering(soupx_inputs$filtered)
# soup_channels <- estimate_contamination(soup_channels, clusters = seurat_pre$seurat_clusters)
# 
# # Or, use the overarching wrapper to process and build the final clean object
# clean_matrices <- process_soupx_samples(soup_channels)
# seurat_obj <- create_final_seurat(clean_matrices)
# 
# # Post-clustering ambient RNA refinement and metadata tracking
# seurat_obj <- run_soupx_post_clustering(seurat_obj)
# seurat_obj <- addSoupXMetaToSeurat(seurat_obj, metadata_df)

## -----------------------------------------------------------------------------
# 
# # Choose your processing pipeline:
# # seurat_obj <- ProcessSeuratLOG(seurat_obj) # Standard log-normalization
# seurat_obj <- ProcessSeuratSCT(seurat_obj)   # SCTransform (Recommended)
# 
# # Run PCA, UMAP, and clustering
# seurat_obj <- ClusterAndUMAP(seurat_obj, resolution = 0.5)
# 
# # Clean up any stale or unnecessary reductions to save memory
# seurat_obj <- CleanSeuratReductions(seurat_obj)
# 
# # Create a custom cohesive color palette for your new clusters
# cluster_colors <- make_cluster_palette(seurat_obj)

## -----------------------------------------------------------------------------
# 
# # Check PCA loadings and metadata correlations
# plot_vizdimloadings(seurat_obj, dims = 1:2)
# plot_pc_metadata_correlation(seurat_obj, metadata_cols = c("nCount_RNA", "percent.mt"))
# 
# # Assess technical contribution per cluster
# plot_technical_contribution(seurat_obj, feature = "percent.mt")
# plot_stacked_technical_contribution(seurat_obj, features = c("percent.mt", "percent.ribo"))

## -----------------------------------------------------------------------------
# 
# # Visualize UMAP split by experimental factors
# PlotDimByFactors(seurat_obj, group.by = "Condition")
# plot_combined(seurat_obj, reduction = "umap", group.by = "Condition")
# 
# # Statistical evaluation of cluster distributions and metadata
# plot_cluster_distributions(seurat_obj, group.by = "SampleID")
# plot_cell_abundance(seurat_obj, group.by = "Condition", test = "wilcox")
# plot_metadata_stats(seurat_obj, feature = "score", group.by = "Cluster")

## -----------------------------------------------------------------------------
# 
# # Visualize gene pair correlations and overlap
# plot_gene_pair_correlations(seurat_obj, feature1 = "CD8A", feature2 = "XCL1")
# plot_blend_nebulosa(seurat_obj, features = c("CD8A", "XCL1"))
# 
# # Internal V5 blend wrapper for strict layer management
# # .blend_feature_plot_v5(seurat_obj, features = c("CD4", "IL2RA"))

## -----------------------------------------------------------------------------
# 
# # Generate standard dotplots with hierarchical clustering trees
# plot_dot_dendro(seurat_obj, features = top_markers)
# 
# # Split dotplots by condition or metadata
# plot_dot_dendro_split(seurat_obj, features = top_markers, split.by = "Condition")
# plot_dot_dendro_multi(seurat_obj, feature_list = list(Tcells = c("CD3D"), Bcells = c("CD79A")))
# plot_split_dotplots_by_gene_cluster(seurat_obj, features = top_markers)
# plot_gene_markers_with_dendro(seurat_obj, features = top_markers)

## -----------------------------------------------------------------------------
# 
# # Extract binned expression and plot top mixed genes
# binned_expr <- extract_binned_expression(seurat_obj, bins = 10)
# top_mixed <- get_top_mixed_genes(seurat_obj)
# 
# # Comprehensive heatmaps
# FindTopMarkersAndHeatmap(seurat_obj, group.by = "seurat_clusters")
# generate_and_save_heatmap(seurat_obj, features = top_markers, filename = "heatmap.png")
# generate_dimheatmaps(seurat_obj, dims = 1:5)

## -----------------------------------------------------------------------------
# 
# # Generate a trio of volcano plots comparing multiple conditions
# generate_volcano_trio(seurat_obj, ident.1 = "Infected", ident.2 = "Control")
# 
# # Run global enrichment profiling on your markers
# # ora_results <- run_global_ora(seurat_obj, database = "GO")
# # gsea_results <- run_global_gsea(seurat_obj, database = "KEGG")

