# Embedding overlays, distributions, dot plots and co-expression -------------------

test_that("make_cluster_palette returns n distinct colours for every palette", {
  for (p in c("auto", "kelly", "alphabet", "polychrome")) {
    cols <- make_cluster_palette(12, palette = p)
    expect_length(cols, 12)
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}", cols)))
  }
  # interpolates beyond the palette size
  expect_length(make_cluster_palette(60, palette = "kelly"), 60)
})

test_that("PlotDimByFactors returns and saves one plot per factor", {
  obj <- small_seurat()
  out <- test_out_dir("dimplots")

  plots <- quietly(PlotDimByFactors(obj, factors = c("cell_type", "batch"), sample_name = "small",
                                    reduction = "tsne", output_dir = out, width = 4, height = 4,
                                    dpi = 50, verbose = FALSE))

  expect_named(plots, c("cell_type", "batch"))
  expect_true(all(vapply(plots, inherits, logical(1), what = "ggplot")))
  expect_length(written_files(out, "^small_DimPlot_"), 2)
})

test_that("PlotDimByFactors validates its inputs", {
  obj <- small_seurat()
  out <- test_out_dir("dimplots_validate")
  expect_error(PlotDimByFactors(obj, factors = "cell_type", reduction = "nope", output_dir = out),
               "not found")
  expect_error(suppressWarnings(PlotDimByFactors(obj, factors = "not_a_column", reduction = "tsne",
                                                 output_dir = out, verbose = FALSE)),
               "None of the requested factors")
  expect_warning(quietly_warn <- PlotDimByFactors(obj, factors = c("cell_type", "not_a_column"),
                                                  reduction = "tsne", output_dir = out, dpi = 50,
                                                  verbose = FALSE),
                 "not found in metadata")
  expect_named(quietly_warn, "cell_type")
})

test_that("plot_combined pairs feature and density plots for genes and metadata", {
  skip_if_not_installed("Nebulosa")
  skip_if_not_installed("SingleCellExperiment")
  obj <- small_seurat()
  out <- test_out_dir("combined")

  res <- quietly(plot_combined(obj, cluster_col = "cell_type",
                               genes_list = list(B = c("MS4A1", "CD79A"), QC = "nCount_RNA",
                                                 Missing = "NOT_A_GENE"),
                               sample_name = "small", reduction = "tsne", output_dir = out,
                               width = 6, dpi = 30, verbose = FALSE))

  expect_named(res, c("combined_plots", "output_dir"))
  expect_named(res$combined_plots, c("B", "QC"))      # "Missing" is skipped
  expect_s3_class(res$combined_plots$B, "patchwork")
  expect_length(written_files(out, "CombinedPlot"), 2)
})

test_that("plot_violin_qc_single saves QC and gene violins", {
  obj <- small_seurat()
  out <- test_out_dir("violins")

  quietly(plot_violin_qc_single(obj, sample_name = "small",
                                qc_meta_features = c("nCount_RNA", "nFeature_RNA"),
                                gene_features = c("MS4A1", "LYZ", "NOT_A_GENE"),
                                res_col = "cell_type", output_dir = out, width = 6, dpi = 30))

  expect_length(written_files(out, "^small_ViolinQC_cell_type"), 1)
  expect_length(written_files(out, "^small_ViolinGene_cell_type"), 1)
  expect_warning(res <- plot_violin_qc_single(obj, res_col = "not_a_column", output_dir = out),
                 "not found")
  expect_null(res)
})

test_that("extract_binned_expression bins each gene into terciles", {
  obj <- small_seurat()

  binned <- extract_binned_expression(obj, gene_list = c("MS4A1", "LYZ", "GNLY"),
                                      group_col = "cell_type")

  expect_s3_class(binned, "data.frame")
  expect_named(binned, c("Cluster", "Gene", "AvgExpression", "PctExpress", "Expression_Level"))
  expect_equal(nrow(binned), 3 * 3)
  expect_equal(levels(binned$Expression_Level), c("Low", "Int", "High"))
})

test_that("plot_dot_dendro returns a patchwork and saves its companion files", {
  obj <- small_seurat()
  out <- test_out_dir("dot_dendro")

  p <- quietly(plot_dot_dendro(obj, meta_col = "cell_type", feature_df = small_feature_df(),
                               prefix = "Test", pct_threshold = 0, output_dir = out))

  expect_s3_class(p, "patchwork")
  expect_length(written_files(out, "^DotPlot_RowColDendro_cell_type"), 1)
  expect_length(written_files(out, "^OrderedGenes_cell_type"), 1)
  # an unknown prefix or column is handled gracefully
  expect_null(quietly(plot_dot_dendro(obj, "cell_type", small_feature_df(), prefix = "Nope",
                                      output_dir = out)))
  expect_error(plot_dot_dendro(obj, "not_a_column", small_feature_df(), output_dir = out),
               "does not exist")
})

test_that("plot_dot_dendro_split chunks long gene panels", {
  obj <- small_seurat()
  out <- test_out_dir("dot_split")

  chunks <- quietly(plot_dot_dendro_split(obj, meta_col = "cell_type", feature_df = small_feature_df(),
                                          prefix = "Test", pct_threshold = 0,
                                          max_genes_per_plot = 4, output_dir = out))

  expect_length(chunks, 3)          # 10 genes / 4 per chunk
  expect_true(all(vapply(chunks, inherits, logical(1), what = "patchwork")))
  expect_length(written_files(out, "_chunk"), 3)
  expect_length(written_files(out, "^DotPlot_GeneList_"), 1)
})

test_that("plot_dot_dendro_multi combines several metadata groupings", {
  obj <- small_seurat()
  out <- test_out_dir("dot_multi")

  df <- quietly(plot_dot_dendro_multi(obj, meta_cols = c("cell_type", "batch"),
                                      feature_df = small_feature_df(), prefix = "Test",
                                      pct_threshold = 0, output_dir = out,
                                      save_options = c("combined", "gene_list")))

  expect_s3_class(df, "data.frame")
  expect_true(all(c("C0_cell_type", "B1_batch") %in% as.character(df$cluster_var)))
  expect_length(written_files(out, "^DotPlot_RowColDendro_Test"), 1)
  expect_length(written_files(out, "^OrderedGenes_Test"), 1)
  expect_length(written_files(out, "Dendrogram"), 0)      # dendrograms not requested
})

test_that(".blend_feature_plot_v5 returns a four-panel patchwork", {
  obj <- small_seurat()
  p <- .blend_feature_plot_v5(obj, gene1 = "MS4A1", gene2 = "LYZ", reduction = "tsne")
  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 3)    # 3 stored patches + the last plot = 4 panels
})

test_that("plot_blend_nebulosa stacks one row per valid gene pair", {
  skip_if_not_installed("Nebulosa")
  obj <- small_seurat()
  out <- test_out_dir("blend")

  res <- quietly(plot_blend_nebulosa(obj, cluster_col = "cell_type",
                                     gene_pairs = list(c("MS4A1", "CD79A"), c("LYZ", "S100A9"),
                                                       c("LYZ", "NOT_A_GENE")),
                                     sample_name = "small", reduction = "tsne", width = 12,
                                     dpi = 30, output_dir = out, verbose = FALSE))

  expect_named(res, c("plot", "rows", "output_dir", "dimensions"))
  expect_length(res$rows, 2)
  expect_equal(res$dimensions$width, 12)
  expect_length(written_files(out, "BlendNebulosa"), 1)
  expect_error(plot_blend_nebulosa(obj, "cell_type", gene_pairs = list("MS4A1")), "length 2")
})

test_that("plot_gene_pair_correlations returns per-cluster correlations", {
  obj <- small_seurat()
  out <- test_out_dir("gene_pairs")

  res <- quietly(plot_gene_pair_correlations(obj, cluster_col = "cell_type",
                                             gene_pairs = list(c("MS4A1", "CD79A"), c("LYZ", "S100A9")),
                                             sample_name = "small", cor_method = "spearman",
                                             output_dir = out, dpi = 30, verbose = FALSE))

  expect_named(res, c("plot", "data", "output_dir"))
  expect_s3_class(res$plot, "ggplot")
  expect_true(all(c("Cluster", "Correlation", "N_Cells", "Pair") %in% names(res$data)))
  expect_setequal(unique(res$data$Pair), c("MS4A1 vs CD79A", "LYZ vs S100A9"))
  expect_true(all(abs(res$data$Correlation) <= 1, na.rm = TRUE))
  expect_error(plot_gene_pair_correlations(obj, "cell_type", list(c("A", "B", "C"))), "length 2")
  expect_error(plot_gene_pair_correlations(obj, "not_a_column", list(c("MS4A1", "LYZ"))), "not found")
})
