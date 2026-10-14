# Regression tests for the bugs fixed in the development version ---------------

# -- Ambient RNA (SoupX) ---------------------------------------------------------

test_that("create_final_seurat uses the counts corrected by estimate_contamination", {
  inp <- mock_contamination_inputs()
  soup <- lapply(inp$soup, SoupX::setContaminationFraction, contFrac = 0.05)
  # a recognisable "corrected" matrix: every count doubled
  corrected <- lapply(inp$toc, function(m) m * 2)

  default <- quietly(create_final_seurat(soup, inp$samples))
  given   <- quietly(create_final_seurat(soup, inp$samples, corrected_list = corrected))

  counts <- function(o) SeuratObject::LayerData(o, layer = "counts")
  expect_equal(sum(counts(given$S1)), sum(corrected$S1[rownames(given$S1), colnames(given$S1)]))
  expect_false(isTRUE(all.equal(sum(counts(default$S1)), sum(counts(given$S1)))))
})

test_that("run_soupx_post_clustering passes the method-specific corrected counts on", {
  seen <- NULL
  local_mocked_bindings(
    setClusters  = function(sc, clusters) sc,
    autoEstCont  = function(sc, ...) { sc$metaData$rho <- 0.05; sc },
    adjustCounts = function(sc, method = "subtraction", ...) {
      structure(sc$toc, method = method)
    },
    create_final_seurat = function(soup_list_all, sample_names, corrected_list = NULL) {
      seen <<- corrected_list
      stats::setNames(as.list(sample_names), sample_names)
    }
  )
  inp <- mock_contamination_inputs()
  prep <- list(sample_names = inp$samples, soup_list = inp$soup,
               seurat_for_clustering = inp$seur, counts_lists = list(filtered = inp$toc))

  quietly(run_soupx_post_clustering(prep, methods = "multinomial"))

  expect_named(seen, inp$samples)
  expect_equal(attr(seen$S1, "method"), "multinomial")
})

test_that("estimate_contamination keeps going when autoEstCont fails for a sample", {
  local_mocked_bindings(
    setClusters  = function(sc, clusters) sc,
    autoEstCont  = function(sc, ...) stop("No plausible marker genes found"),
    adjustCounts = function(sc, ...) sc$toc
  )
  inp <- mock_contamination_inputs()

  expect_warning(
    res <- quietly_messages(utils::capture.output(
      r <- estimate_contamination(inp$soup, inp$seur, inp$toc, inp$samples, mulFac = 2)
    )),
    "autoEstCont\\(\\) failed"
  )
  rho <- unname(sapply(r$soup_list_all, function(sc) unique(sc$metaData$rho)))
  expect_equal(rho, c(0.02, 0.02))
})

test_that("estimate_contamination applies a fraction above 0.5 with a warning and rejects > 1", {
  local_mocked_bindings(
    setClusters  = function(sc, clusters) sc,
    autoEstCont  = function(sc, ...) { sc$metaData$rho <- 0.05; sc },
    adjustCounts = function(sc, ...) sc$toc
  )
  inp <- mock_contamination_inputs()

  expect_warning(
    quietly_messages(utils::capture.output(
      r <- estimate_contamination(inp$soup, inp$seur, inp$toc, inp$samples, mulFac = 55)
    )),
    "above 0.5"
  )
  expect_equal(unname(sapply(r$soup_list_all, function(sc) unique(sc$metaData$rho))), c(0.6, 0.6))

  expect_error(
    quietly(estimate_contamination(inp$soup, inp$seur, inp$toc, inp$samples, mulFac = 100)),
    "exceeds 1"
  )
})

test_that("load_counts reads Cell Ranger v2 (genes.tsv, uncompressed) directories", {
  root <- test_out_dir("cellranger_v2")
  m <- synthetic_counts(n_cells = 20)
  for (kind in c("filtered", "raw")) {
    d <- file.path(root, kind, paste0("S1_", kind, "_feature_bc_matrix"))
    dir.create(d, recursive = TRUE)
    Matrix::writeMM(methods::as(m, "CsparseMatrix"), file.path(d, "matrix.mtx"))
    writeLines(paste(rownames(m), rownames(m), sep = "\t"), file.path(d, "genes.tsv"))
    writeLines(colnames(m), file.path(d, "barcodes.tsv"))
  }

  counts <- quietly(load_counts("S1", file.path(root, "filtered"), file.path(root, "raw")))

  expect_equal(dim(counts$filtered$S1), dim(m))
  expect_equal(rownames(counts$raw$S1), rownames(m))
  expect_error(OHMmy:::.find_10x_file(root, "nothing.tsv"), "nothing.tsv")
})

test_that("create_seurat_for_clustering works on samples too small for 50 PCs", {
  toc <- list(S1 = synthetic_counts(n_cells = 40, seed = 3))

  res <- quietly(create_seurat_for_clustering(toc, "S1"))

  expect_s4_class(res$S1, "Seurat")
  expect_lt(ncol(SeuratObject::Embeddings(res$S1, "pca")), 40)
  expect_true("umap" %in% SeuratObject::Reductions(res$S1))
})

# -- Clustering ----------------------------------------------------------------------

test_that("ClusterAndUMAP names mixed-precision resolutions as Seurat does", {
  res <- quietly(ClusterAndUMAP(SeuratObject::pbmc_small, dims = 1:3, reduction = "pca",
                                umap_name = "umap.test", cluster_resolutions = c(0.25, 0.5, 1),
                                final_resolution = 0.5, plot_dir = test_out_dir("clustree_mixed"),
                                dpi = 30))

  expect_true(all(c("RNA_snn_res.0.25", "RNA_snn_res.0.5", "RNA_snn_res.1") %in%
                    colnames(res$seurat@meta.data)))
  expect_false(any(c("RNA_snn_res.0.50", "RNA_snn_res.1.00") %in% colnames(res$seurat@meta.data)))
})

test_that("ClusterAndUMAP activates final_resolution when the UMAP already exists", {
  out <- test_out_dir("clustree_final")
  first <- quietly(ClusterAndUMAP(SeuratObject::pbmc_small, dims = 1:3, reduction = "pca",
                                  umap_name = "umap.test", cluster_resolutions = c(0.5, 1.5),
                                  final_resolution = 0.5, plot_dir = out, dpi = 30))
  # Re-run with another final resolution: the UMAP exists, the sweep ends at 1.5
  second <- quietly(ClusterAndUMAP(first$seurat, dims = 1:3, reduction = "pca",
                                   umap_name = "umap.test", cluster_resolutions = c(0.5, 1.5),
                                   final_resolution = 0.5, plot_dir = out, dpi = 30))

  md <- second$seurat@meta.data
  expect_equal(as.character(md$seurat_clusters), as.character(md$RNA_snn_res.0.5))
  expect_equal(as.character(SeuratObject::Idents(second$seurat)), as.character(md$RNA_snn_res.0.5))
})

test_that("ClusterAndUMAP passes algorithm to the final clustering", {
  algs <- c()
  local_mocked_bindings(
    FindClusters = function(object, resolution, graph.name, algorithm = 1, ...) {
      algs <<- c(algs, algorithm)
      object@meta.data[[paste0(graph.name, "_res.", resolution)]] <- factor(rep(0:1, length.out = ncol(object)))
      object$seurat_clusters <- object@meta.data[[paste0(graph.name, "_res.", resolution)]]
      object
    }
  )
  quietly(ClusterAndUMAP(SeuratObject::pbmc_small, dims = 1:3, reduction = "pca", algorithm = 4,
                         umap_name = "umap.test", cluster_resolutions = 0.5, final_resolution = 0.5,
                         plot_dir = test_out_dir("clustree_alg"), dpi = 30))
  expect_true(length(algs) >= 2)
  expect_true(all(algs == 4))
})

# -- Donor-level statistics ---------------------------------------------------------

test_that("plot_metadata_stats counts each sample once when its values vary", {
  obj <- small_seurat()
  obj$age[obj$donor == "D1"][1] <- 99   # one cell of D1 with a different age

  expect_warning(
    res <- quietly_messages(utils::capture.output(
      r <- plot_metadata_stats(obj, sample_col = "donor", condition_col = "condition",
                               metadata_vars = "age", output_dir = test_out_dir("metadata_dup"),
                               dpi = 30)
    )),
    "not constant within sample"
  )
  expect_equal(r$stats$age$Global$n1 + r$stats$age$Global$n2, 8)
})

test_that("plot_cell_abundance computes proportions within each sample x condition", {
  obj <- small_seurat()
  obj$condition[obj$donor == "D1"][1:3] <- "Case"   # D1 appears in both conditions

  expect_warning(
    quietly_messages(utils::capture.output(
      p <- plot_cell_abundance(obj, "donor", "condition", "cell_type",
                               output_dir = test_out_dir("abundance_multi"),
                               base_size = 2, dpi = 30)
    )),
    "more than one condition"
  )
  sums <- tapply(p$data$Proportion, paste(p$data$Sample, p$data$Condition), sum)
  expect_equal(as.numeric(sums), rep(1, length(sums)))
})

# -- General visualisation ----------------------------------------------------------

test_that("make_cluster_palette returns distinct colours and keeps small palettes unchanged", {
  expect_equal(anyDuplicated(make_cluster_palette(24, "kelly")), 0)
  expect_equal(anyDuplicated(make_cluster_palette(30, "polychrome")), 0)
  expect_false("#F2F3F4" %in% make_cluster_palette(24, "kelly"))
  expect_equal(make_cluster_palette(3, "kelly"), c("#F3C300", "#875692", "#F38400"))
  expect_equal(make_cluster_palette(3, "polychrome"), c("#5A5156", "#E4E1E3", "#F6222E"))
})

test_that("plot_gene_pair_correlations reports NA for clusters with too few cells", {
  obj <- small_seurat()
  expr <- SeuratObject::FetchData(obj, c("CD79A", "MS4A1"))
  both <- rownames(expr)[expr$CD79A > 0 & expr$MS4A1 > 0]
  obj$tiny <- ifelse(colnames(obj) %in% both[1:2], "Tiny", "Big")   # 2 co-expressing cells

  run <- function(...) quietly(plot_gene_pair_correlations(
    obj, "tiny", list(c("CD79A", "MS4A1")), quantile_thresh = 0,
    output_dir = test_out_dir("pair_min"), dpi = 30, add_timestamp = FALSE, verbose = FALSE, ...))
  res3 <- run()
  res1 <- run(min_cells = 1)

  tiny3 <- res3$data[res3$data$Cluster == "Tiny", ]
  tiny1 <- res1$data[res1$data$Cluster == "Tiny", ]
  expect_equal(tiny3$N_Cells, 2)
  expect_true(is.na(tiny3$Correlation))
  expect_equal(abs(tiny1$Correlation), 1)          # the spurious value min_cells removes
  big <- res3$data$Cluster == "Big"
  expect_equal(res3$data$Correlation[big], res1$data$Correlation[big])
})

test_that(".blend_feature_plot_v5 honours a numeric min_cutoff", {
  obj <- small_seurat()
  p_none <- OHMmy:::.blend_feature_plot_v5(obj, "CD79A", "MS4A1", "tsne", min_cutoff = NA)
  p_num  <- OHMmy:::.blend_feature_plot_v5(obj, "CD79A", "MS4A1", "tsne", min_cutoff = 5.5)   # about the median of expressing cells
  g_none <- p_none[[1]]$data$gene1
  g_num  <- p_num[[1]]$data$gene1
  expect_gt(sum(g_none > 0), sum(g_num > 0))
})

test_that(".blend_feature_plot_v5 works on split (unjoined) Seurat v5 layers", {
  obj <- small_seurat()
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)
  expect_gt(length(SeuratObject::Layers(obj[["RNA"]], search = "data")), 1)

  p <- OHMmy:::.blend_feature_plot_v5(obj, "CD79A", "MS4A1", "tsne")
  expect_s3_class(p, "patchwork")
  expect_equal(nrow(p[[1]]$data), ncol(obj))
})

# -- Markers -------------------------------------------------------------------------

test_that("FindTopMarkersAndHeatmap writes tables but no heatmap when no marker passes", {
  obj <- small_seurat()
  out <- test_out_dir("markers_none_pass")

  expect_warning(
    quietly_messages(utils::capture.output(
      res <- FindTopMarkersAndHeatmap(obj, sample_name = "S", marker_avg_log2FC_thresh = 1e6,
                                      output_dir_base = out, add_timestamp = FALSE, dpi = 30)
    )),
    "No markers passed the thresholds"
  )
  expect_null(res$heatmap)
  expect_equal(nrow(res$top_markers), 0)
  files <- written_files(out, recursive = TRUE)
  expect_true(any(grepl("_all_markers_", files)))
  expect_false(any(grepl("\\.jpg$", files)))
})

# -- Technical noise -----------------------------------------------------------------

test_that("technical contribution functions cap n_top_genes at the number of loaded genes", {
  obj <- SeuratObject::pbmc_small   # 20 genes in the PCA loadings
  out <- test_out_dir("tech_cap")

  expect_no_error(quietly(plot_technical_contribution(obj, "S", reduction = "pca", output_dir = out,
                                                      max_pcs = 3, n_top_genes = 500)))
  expect_no_error(quietly(plot_stacked_technical_contribution(obj, "S", reduction = "pca",
                                                              output_dir = out, max_pcs = 3,
                                                              gene_depths = c(20, 100), dpi = 30)))
})

test_that("plot_pc_metadata_correlation skips non-numeric metadata with a warning", {
  obj <- SeuratObject::pbmc_small
  expect_warning(
    f <- plot_pc_metadata_correlation(obj, "S", vars_to_test = c("nCount_RNA", "nFeature_RNA", "groups"),
                                      reduction = "pca", n_pcs = 5,
                                      output_dir = test_out_dir("pc_cor_nonnum"), res_dpi = 30),
    "non-numeric"
  )
  expect_true(file.exists(f))
})

# -- DESeq2 heatmap ------------------------------------------------------------------

test_that("generate_and_save_heatmap drops genes with constant expression", {
  skip_if_not_installed("pheatmap")
  set.seed(1)
  genes <- paste0("G", 1:6)
  res <- data.frame(padj = c(1e-5, 1e-4, 1e-3, 1e-2, 0.02, 0.03),
                    log2FoldChange = c(3, -2, 1.5, -1, 2, 0.5), row.names = genes)
  mat <- matrix(rnorm(36, 8), 6, dimnames = list(genes, paste0("s", 1:6)))
  mat["G1", ] <- 5   # constant gene
  anno <- data.frame(group = rep(c("A", "B"), each = 3), row.names = colnames(mat))
  out <- test_out_dir("heatmap_const")

  expect_message(
    f <- generate_and_save_heatmap(res, "AvsB", "A vs B", mat, anno, colnames(mat),
                                   n_padj = 6, n_lfc = 6, out_dir = out, ts = "t"),
    "constant expression"
  )
  expect_true(file.exists(f))
})

# -- Dot plot dendrograms ------------------------------------------------------------

test_that("plot_dot_dendro_split and plot_dot_dendro_multi skip a single gene instead of failing", {
  obj <- small_seurat()
  one_gene <- data.frame(group = "Test_B", feature = "MS4A1")

  expect_message(
    r1 <- plot_dot_dendro_split(obj, "cell_type", one_gene, pct_threshold = 0,
                                output_dir = test_out_dir("dd_split_one")),
    "Skipping plot"
  )
  expect_null(r1)
  expect_message(
    r2 <- plot_dot_dendro_multi(obj, c("cell_type", "batch"), one_gene, pct_threshold = 0,
                                prefix = "Test", output_dir = test_out_dir("dd_multi_one")),
    "Skipping plot"
  )
  expect_null(r2)
})
