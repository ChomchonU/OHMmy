# CD8 vs NK lineage resolution -------------------------------------------------------

# Cells alternate between a clear CD8 profile and a clear NK profile
consensus_fixture <- function() {
  obj <- small_seurat()
  set.seed(1)
  n <- ncol(obj)
  is_cd8 <- rep(c(TRUE, FALSE), length.out = n)
  obj$CD8_a    <- ifelse(is_cd8, 2, -2) + stats::rnorm(n, sd = 0.4)
  obj$CD8_b    <- ifelse(is_cd8, 2, -2) + stats::rnorm(n, sd = 0.4)
  obj$NK_Score <- ifelse(is_cd8, -2, 2) + stats::rnorm(n, sd = 0.4)
  list(obj = obj, is_cd8 = is_cd8)
}

test_that("score_consensus labels cells with a fixed cutoff and writes metadata", {
  fx <- consensus_fixture()

  res <- score_consensus(fx$obj, cd8_cols = c("CD8_a", "CD8_b"), method = "fixed")

  expect_named(res, c("obj", "per", "labels", "consensus", "agreement", "n_disagree", "cutoffs",
                      "bootstrap", "plots", "cells_used", "dropped", "params"))
  expect_equal(unname(res$consensus), ifelse(fx$is_cd8, "CD8", "NK"))
  expect_equal(res$agreement["CD8_a", "CD8_b"], 1)
  expect_equal(res$n_disagree, 0)
  expect_true(all(c("CD8vNK_consensus", "CD8vNK_nkscale", "CD8vNK_diff_CD8_a") %in%
                    colnames(res$obj@meta.data)))
  expect_null(res$bootstrap)
})

test_that("score_consensus finds data-driven cutoffs with trough and gmm", {
  fx <- consensus_fixture()

  trough <- score_consensus(fx$obj, cd8_cols = c("CD8_a", "CD8_b"), method = "trough")
  expect_false(any(trough$cutoffs$fallback))
  expect_equal(unname(trough$consensus), ifelse(fx$is_cd8, "CD8", "NK"))

  skip_if_not_installed("mclust")
  gmm <- score_consensus(fx$obj, cd8_cols = c("CD8_a", "CD8_b"), method = "gmm")
  expect_false(any(gmm$cutoffs$fallback))
  expect_equal(unname(gmm$consensus), ifelse(fx$is_cd8, "CD8", "NK"))
})

test_that("score_consensus flags disagreements between primary and confirm scores", {
  fx <- consensus_fixture()
  fx$obj$CD8_b[1:4] <- -10                     # confirm score disagrees on the first 4 cells

  res <- score_consensus(fx$obj, cd8_cols = c("CD8_a", "CD8_b"), method = "fixed")

  expect_true(all(res$consensus[c(1, 3)] == "Ambiguous_CD8_primary"))
  expect_gt(res$n_disagree, 0)
  expect_lt(res$agreement["CD8_a", "CD8_b"], 1)
})

test_that("score_consensus bootstraps, filters cells and draws plots", {
  fx <- consensus_fixture()
  out <- test_out_dir("consensus_plots")

  res <- quietly(score_consensus(fx$obj, cd8_cols = c("CD8_a", "CD8_b"), method = "trough",
                                 min_features = stats::median(fx$obj$nFeature_RNA),
                                 n_boot = 20, make_plots = TRUE, plot_n = 20, out_dir = out))

  expect_gt(length(res$dropped), 0)
  expect_equal(length(res$cells_used) + length(res$dropped), ncol(fx$obj))
  expect_true(all(is.na(res$obj$CD8vNK_consensus[res$dropped])))
  expect_s3_class(res$bootstrap$summary, "data.frame")
  expect_equal(nrow(res$bootstrap$summary), 2)
  expect_true(all(c("Violin_CD8_a", "Agreement_Matrix", "Primary_vs_Confirm",
                    "Consensus_Barplot", "Bootstrap_ForestPlot") %in% names(res$plots)))
  expect_gte(length(written_files(out, "\\.png$")), 5)
})

test_that("score_consensus validates its inputs", {
  fx <- consensus_fixture()
  expect_error(score_consensus(fx$obj, cd8_cols = "not_a_column"), "Missing meta.data")
  expect_error(score_consensus(fx$obj, cd8_cols = "CD8_a", primary = "CD8_b"), "'primary' not in")
  expect_warning(score_consensus(fx$obj, cd8_cols = "CD8_a", method = "fixed", n_boot = 5),
                 "Switching boot_target")
})
