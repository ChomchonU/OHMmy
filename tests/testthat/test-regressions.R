# Regression tests for bugs fixed in OHMmy 1.0.2

test_that("get_sample_names returns unnamed, cleaned sample IDs", {
  root <- file.path(tempdir(), "cr_names")
  dir.create(file.path(root, "S1_filtered_feature_bc_matrix"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "S2_filtered_feature_bc_matrix"), recursive = TRUE, showWarnings = FALSE)

  ids <- get_sample_names(root)
  expect_null(names(ids))
  expect_setequal(ids, c("S1", "S2"))

  unlink(root, recursive = TRUE)
})

test_that("score_consensus finds a GMM cutoff without mclust being attached", {
  skip_if_not_installed("mclust")
  if ("package:mclust" %in% search()) detach("package:mclust")

  set.seed(1)
  obj <- SeuratObject::pbmc_small
  n <- ncol(obj)
  cd8_group <- rep(c(TRUE, FALSE), length.out = n)
  obj$CD8_a <- ifelse(cd8_group, 2, -2) + rnorm(n, sd = 0.4)
  obj$CD8_b <- ifelse(cd8_group, 2, -2) + rnorm(n, sd = 0.4)
  obj$NK_Score <- ifelse(cd8_group, -2, 2) + rnorm(n, sd = 0.4)

  res <- score_consensus(obj, cd8_cols = c("CD8_a", "CD8_b"), method = "gmm")
  expect_false(any(res$cutoffs$fallback))
  expect_true(all(res$consensus[cd8_group] == "CD8"))
  expect_true(all(res$consensus[!cd8_group] == "NK"))
})

test_that("plot_pc_metadata_correlation writes a non-empty file while another device is open", {
  obj <- SeuratObject::pbmc_small
  out <- file.path(tempdir(), "pc_cor_test")

  grDevices::pdf(NULL)                         # simulate an open knitr / RStudio device
  on.exit(grDevices::dev.off(), add = TRUE)

  f <- suppressWarnings(plot_pc_metadata_correlation(
    obj, sample_name = "test", vars_to_test = c("nCount_RNA", "nFeature_RNA", "not_a_column"),
    reduction = "pca", n_pcs = 5, output_dir = out, res_dpi = 50
  ))
  expect_true(file.exists(f))
  expect_gt(file.size(f), 0)
  expect_warning(
    plot_pc_metadata_correlation(obj, "test", vars_to_test = c("nCount_RNA", "nFeature_RNA", "not_a_column"),
                                 reduction = "pca", n_pcs = 5, output_dir = out, res_dpi = 50),
    "not_a_column"
  )
})

test_that("plot_dot_dendro creates its output directory", {
  obj <- SeuratObject::pbmc_small
  out <- file.path(tempdir(), "dot_dendro_new_dir")
  unlink(out, recursive = TRUE)

  feature_df <- data.frame(
    group   = c("Test_A", "Test_A", "Test_B", "Test_B"),
    feature = c("MS4A1", "CD79A", "LYZ", "S100A9")
  )
  p <- suppressWarnings(plot_dot_dendro(obj, meta_col = "groups", feature_df = feature_df,
                                        prefix = "Test", pct_threshold = 0, output_dir = out))
  expect_s3_class(p, "patchwork")
  expect_true(dir.exists(out))
})

test_that("generate_and_save_heatmap accepts a plain matrix and creates out_dir", {
  set.seed(2)
  mat <- matrix(rnorm(20 * 6), nrow = 20,
                dimnames = list(paste0("G", 1:20), paste0("S", 1:6)))
  mat[1:5, 4:6] <- mat[1:5, 4:6] + 5
  res <- data.frame(log2FoldChange = c(rep(3, 5), rep(0, 15)),
                    padj = c(rep(1e-5, 5), rep(0.9, 15)),
                    row.names = rownames(mat))
  anno <- data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(mat))
  out <- file.path(tempdir(), "heatmap_new_dir")
  unlink(out, recursive = TRUE)

  f <- generate_and_save_heatmap(res, "B_vs_A", "B vs A", mat, anno, colnames(mat),
                                 n_padj = 5, n_lfc = 5, out_dir = out, ts = "test")
  expect_true(file.exists(f))
})
