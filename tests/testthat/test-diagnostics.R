# PC diagnostics and technical-noise screening -----------------------------------

test_that("generate_dimheatmaps saves one file per PC window", {
  obj <- small_seurat()
  out <- test_out_dir("dimheatmaps")

  failed <- quietly(generate_dimheatmaps(obj, sample_name = "small", reduction = "pca",
                                         pc_windows = list(1:2, 3:4), output_dir = out,
                                         nfeatures = 5, cells = 30, width_in = 4,
                                         height_in = 4, res = 50))

  expect_length(failed, 0)
  expect_length(written_files(out, "^small_DimHeatmap_pca_PC"), 2)
})

test_that("generate_dimheatmaps reports a missing reduction", {
  obj <- small_seurat()
  failed <- quietly(generate_dimheatmaps(obj, "small", reduction = "nope",
                                         output_dir = test_out_dir("dimheatmaps_missing")))
  expect_equal(failed, "small")
})

test_that("plot_vizdimloadings saves loadings and returns no failures", {
  obj <- small_seurat()
  out <- test_out_dir("vizdim")

  failed <- quietly(plot_vizdimloadings(obj, sample_name = "small", reduction = "pca",
                                        pc_windows = list(1:4), output_dir = out,
                                        nfeatures = 5, width_in = 6, height_in = 4,
                                        res = 50, ncol = 2))

  expect_identical(failed, character(0))
  expect_length(written_files(out, "^small_VizDimLoadings_PC1-4"), 1)
  expect_equal(quietly(plot_vizdimloadings(obj, "small", reduction = "nope",
                                           output_dir = out)), "small")
})

test_that("plot_technical_contribution saves a plot and flags missing reductions", {
  obj <- small_seurat()
  out <- test_out_dir("tech_contrib")

  failed <- quietly(plot_technical_contribution(obj, sample_name = "small", reduction = "pca",
                                                output_dir = out, technical_keywords = c("^RP", "^MT-"),
                                                max_pcs = 5, n_top_genes = 10))

  expect_identical(failed, character(0))
  expect_length(written_files(out, "^tech_contrib_pca_posneg_split_small"), 1)
  expect_equal(quietly(plot_technical_contribution(obj, "small", reduction = "nope",
                                                   output_dir = out)), "small")
})

test_that("plot_stacked_technical_contribution saves the faceted plot", {
  obj <- small_seurat()
  out <- test_out_dir("tech_stacked")

  f <- quietly(plot_stacked_technical_contribution(obj, sample_name = "small", reduction = "pca",
                                                   technical_keywords = c("^RP", "^MT-"),
                                                   gene_depths = seq(0, 10, by = 5), max_pcs = 4,
                                                   output_dir = out, plot_width = 6,
                                                   plot_height = 4, dpi = 50))

  expect_true(file.exists(f))
  expect_match(basename(f), "^stacked_technical_contrib_posneg_small")
})

test_that("plot_pc_metadata_correlation saves a heatmap of PC-covariate correlations", {
  obj <- small_seurat()
  out <- test_out_dir("pc_cor")

  f <- quietly(plot_pc_metadata_correlation(obj, sample_name = "small",
                                            vars_to_test = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt"),
                                            reduction = "pca", n_pcs = 5, output_dir = out,
                                            plot_width_in = 5, plot_height_in = 3, res_dpi = 50))

  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)
  expect_error(
    quietly(plot_pc_metadata_correlation(obj, "small", vars_to_test = "nCount_RNA",
                                         reduction = "pca", output_dir = out)),
    "At least two"
  )
})
