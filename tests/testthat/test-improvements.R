# Robustness and returned plots ----------------------------------------------------

test_that("CleanSeuratReductions never overwrites an existing reduction", {
  obj <- SeuratObject::pbmc_small
  obj[["pca.log"]] <- obj[["pca"]]
  obj[["pcaLog"]] <- obj[["tsne"]]           # a different reduction already holds the clean name
  before <- SeuratObject::Embeddings(obj, "pcaLog")

  expect_warning(clean <- suppressMessages(CleanSeuratReductions(obj)), "already exists")
  expect_true(all(c("pca.log", "pcaLog") %in% SeuratObject::Reductions(clean)))
  expect_equal(SeuratObject::Embeddings(clean, "pcaLog"), before)
})

test_that("plot_violin_qc_single returns the QC and gene plots", {
  obj <- small_seurat()
  res <- quietly(plot_violin_qc_single(obj, "s", qc_meta_features = c("nCount_RNA", "nFeature_RNA"),
                                       gene_features = c("CD3E", "LYZ"), res_col = "cell_type",
                                       output_dir = test_out_dir("violin_return"), dpi = 30))
  expect_named(res, c("qc", "genes"))
  expect_s3_class(res$qc, "patchwork")
  expect_s3_class(res$genes, "patchwork")

  no_genes <- quietly(plot_violin_qc_single(obj, "s", qc_meta_features = "nCount_RNA", res_col = "cell_type",
                                            output_dir = test_out_dir("violin_return2"), dpi = 30))
  expect_null(no_genes$genes)
})

test_that("plot_dot_dendro_multi attaches the plots it draws", {
  obj <- small_seurat()
  res <- quietly(plot_dot_dendro_multi(obj, c("cell_type", "batch"), small_feature_df(), prefix = "Test",
                                       output_dir = test_out_dir("multi_return")))
  plots <- attr(res, "plots")
  expect_named(plots, c("combined", "dotplot", "row_dendrogram", "col_dendrogram"))
  expect_s3_class(plots$combined, "patchwork")
  expect_s3_class(plots$dotplot, "ggplot")

  only_list <- quietly(plot_dot_dendro_multi(obj, c("cell_type", "batch"), small_feature_df(), prefix = "Test",
                                             save_options = "gene_list",
                                             output_dir = test_out_dir("multi_return2")))
  expect_null(attr(only_list, "plots")$combined)
})

test_that("technical-noise functions attach their plots", {
  obj <- small_seurat()
  f <- quietly(plot_stacked_technical_contribution(obj, "s", reduction = "pca", gene_depths = c(0, 10),
                                                   max_pcs = 3, output_dir = test_out_dir("stack_attr"),
                                                   plot_width = 5, plot_height = 4, dpi = 30))
  expect_s3_class(attr(f, "plot"), "ggplot")

  g <- quietly(plot_pc_metadata_correlation(obj, "s", vars_to_test = c("nCount_RNA", "nFeature_RNA"),
                                            reduction = "pca", n_pcs = 4,
                                            output_dir = test_out_dir("pccor_attr"), res_dpi = 30))
  expect_s3_class(attr(g, "plot"), "pheatmap")
  expect_equal(dim(attr(g, "cor_matrix")), c(2, 4))
})

test_that("statistics plots compute their p-value labels", {
  obj <- small_seurat()
  for (facet in c(TRUE, FALSE)) {
    p <- quietly(plot_cell_abundance(obj, "donor", "severity", "cell_type", facet_by_cluster = facet,
                                     output_dir = test_out_dir("p_labels"), base_size = 2, dpi = 30))
    labels <- unlist(lapply(quietly(ggplot2::ggplot_build(p))$data, function(d) d$label))
    expect_true(any(grepl("^p [=<]", labels)))
  }
  res <- quietly(plot_metadata_stats(obj, "donor", "severity", "age",
                                     output_dir = test_out_dir("p_labels_meta"), dpi = 30))
  labels <- unlist(lapply(quietly(ggplot2::ggplot_build(res$plots$age))$data, function(d) d$label))
  expect_true(any(grepl("^Kruskal-Wallis, p = ", labels)))
})
