# Donor-level statistics ------------------------------------------------------------

test_that("plot_cell_abundance tests two conditions with Mann-Whitney", {
  obj <- small_seurat()
  out <- test_out_dir("abundance_2")

  p <- quietly(plot_cell_abundance(obj, sample_col = "donor", condition_col = "condition",
                                   celltype_col = "cell_type", output_dir = out, base_size = 2,
                                   dpi = 30))

  expect_s3_class(p, "ggplot")
  # proportions are computed per donor and sum to 1
  props <- tapply(p$data$Proportion, p$data$Sample, sum)
  expect_equal(as.numeric(props), rep(1, 8))
  expect_length(written_files(out, "^cell_abundance_facet_by_cluster_condition_"), 1)
})

test_that("plot_cell_abundance handles three conditions and the unfaceted layout", {
  obj <- small_seurat()
  out <- test_out_dir("abundance_3")

  p_kw <- quietly(plot_cell_abundance(obj, "donor", "severity", "cell_type",
                                      global_test = "kruskal.test", output_dir = out,
                                      base_size = 2, dpi = 30))
  p_anova <- quietly(plot_cell_abundance(obj, "donor", "severity", "cell_type",
                                         global_test = "anova", facet_by_cluster = FALSE,
                                         output_dir = out, base_size = 2, dpi = 30))

  expect_s3_class(p_kw, "ggplot")
  expect_s3_class(p_anova, "ggplot")
  expect_setequal(unique(as.character(p_kw$data$Condition)), c("Healthy", "Mild", "Severe"))
  expect_length(written_files(out, "facet_by_celltype"), 1)
})

test_that("plot_metadata_stats routes continuous and categorical variables", {
  obj <- small_seurat()
  out <- test_out_dir("metadata_2")

  res <- quietly(plot_metadata_stats(obj, sample_col = "donor", condition_col = "condition",
                                     metadata_vars = c("age", "sex"), output_dir = out, dpi = 30))

  expect_named(res, c("plots", "stats"))
  expect_s3_class(res$plots$age, "ggplot")
  expect_s3_class(res$plots$sex, "ggplot")
  expect_s3_class(res$stats$age$Global, "data.frame")
  # deduplicated to one row per donor, not per cell
  expect_equal(res$stats$age$Global$n1 + res$stats$age$Global$n2, 8)
  expect_length(written_files(out, "^metadata_"), 2)
})

test_that("plot_metadata_stats supports 3+ groups, Fisher tests and skips bad columns", {
  obj <- small_seurat()

  res <- quietly(plot_metadata_stats(obj, sample_col = "donor", condition_col = "severity",
                                     metadata_vars = c("age", "sex", "not_a_column"),
                                     continuous_test_n3 = "anova", categorical_test = "fisher",
                                     output_dir = test_out_dir("metadata_3"), dpi = 30))

  expect_named(res$plots, c("age", "sex"))
  expect_equal(res$stats$sex$Global$Method, "fisher")
  expect_true(res$stats$sex$Global$p >= 0 && res$stats$sex$Global$p <= 1)
})

test_that("plot_cluster_distributions returns count and proportion plots", {
  obj <- small_seurat()
  out <- test_out_dir("distributions")

  res <- quietly(plot_cluster_distributions(obj, cluster_col = "cell_type", batch_col = "severity",
                                            output_dir = out, file_prefix = "small"))

  expect_named(res, c("proportions", "counts"))
  expect_s3_class(res$proportions, "patchwork")
  expect_length(written_files(out, "^cell_prop_plot_severity_small_"), 1)
  expect_length(written_files(out, "^cell_count_plot_severity_small_"), 1)
  expect_error(plot_cluster_distributions(obj, cluster_col = "nope", output_dir = out), "not found")
})
