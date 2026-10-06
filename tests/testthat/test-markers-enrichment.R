# Markers, enrichment and DESeq2 visualisation ---------------------------------------

test_that("FindTopMarkersAndHeatmap finds markers and writes heatmap + tables", {
  obj <- small_seurat()
  SeuratObject::Idents(obj) <- "cell_type"
  out <- test_out_dir("markers")

  res <- quietly(FindTopMarkersAndHeatmap(obj, sample_name = "small", top_n = 3,
                                          output_dir_base = out, width = 5, height = 5, dpi = 30,
                                          add_timestamp = FALSE))

  expect_named(res, c("top_markers", "heatmap", "markers", "output_dir"))
  expect_s3_class(res$heatmap, "ggplot")
  expect_true(all(c("gene", "cluster", "avg_log2FC", "p_val_adj") %in% names(res$markers)))
  expect_true(all(res$top_markers$avg_log2FC > 0))              # onlyPos = TRUE
  expect_lte(max(table(res$top_markers$cluster)), 3)              # top_n per cluster
  expect_true(file.exists(file.path(out, "small_Heatmap_topMarkers.jpg")))
  expect_length(written_files(out, "\\.csv$"), 3)
})

test_that("FindTopMarkersAndHeatmap supports pairwise comparisons in both directions", {
  obj <- small_seurat()
  SeuratObject::Idents(obj) <- "cell_type"

  res <- quietly(FindTopMarkersAndHeatmap(obj, sample_name = "pair", compare = c("C0", "C1"),
                                          onlyPos = FALSE, output_dir_base = test_out_dir("markers_pair"),
                                          width = 5, height = 5, dpi = 30))

  expect_true(all(c("gene", "cluster") %in% names(res$markers)))
  expect_setequal(unique(res$markers$cluster), c("C0", "C1"))
  expect_equal(res$markers$cluster == "C0", res$markers$avg_log2FC > 0)
})

# A marker table computed separately in two batches, with the `diff` column
multi_marker_table <- function() {
  obj <- small_seurat()
  SeuratObject::Idents(obj) <- "RNA_snn_res.1"
  do.call(rbind, lapply(c("B1", "B2"), function(b) {
    sub <- subset(obj, subset = batch == b)
    m <- suppressWarnings(Seurat::FindAllMarkers(sub, only.pos = TRUE, verbose = FALSE))
    m$diff <- m$pct.1 - m$pct.2
    m$batch <- b
    m
  }))
}

test_that("plot_gene_markers_with_dendro accepts gene vectors and regular expressions", {
  df <- multi_marker_table()
  out <- test_out_dir("marker_dendro")

  res <- quietly(plot_gene_markers_with_dendro(df, genes = c("MS4A1", "CD79A", "LYZ", "S100A9"),
                                               id_col = "batch", output_dir = out))
  expect_named(res, c("dotplot", "dendrogram", "combined"))
  expect_s3_class(res$combined, "patchwork")
  expect_length(written_files(out, "\\.png$"), 3)

  regex_res <- quietly(plot_gene_markers_with_dendro(df, genes = "^HLA-D", id_col = "batch",
                                                     output_dir = test_out_dir("marker_dendro_regex")))
  expect_true(all(grepl("^HLA-D", levels(regex_res$dotplot$data$gene))))
  expect_error(plot_gene_markers_with_dendro(df, genes = "LYZ", id_col = "nope"), "does not exist")
})

test_that("plot_split_dotplots_by_gene_cluster writes a dendrogram and split dot plots", {
  df <- multi_marker_table()
  out <- test_out_dir("split_dot")

  quietly(plot_split_dotplots_by_gene_cluster(df, gene_regex = "^HLA-D|^CD79", id_col = "batch",
                                              output_dir = out, n_splits = 2))

  expect_length(written_files(out, "^Dendrogram_ALL_"), 1)
  expect_length(written_files(out, "^DotPlot_Part"), 2)
  expect_null(quietly(plot_split_dotplots_by_gene_cluster(df, gene_list = "NOT_A_GENE",
                                                          id_col = "batch", output_dir = out)))
  expect_error(plot_split_dotplots_by_gene_cluster(df, gene_list = "LYZ", id_col = "nope"),
               "does not exist")
})

# Synthetic DE table with a planted signal: genes G1-G30 are strongly up in cluster A
synthetic_de <- function() {
  genes <- paste0("G", 1:200)
  set.seed(3)
  make <- function(cl, up) {
    lfc <- stats::rnorm(200, 0, 0.2)
    lfc[up] <- seq(3, 1, length.out = length(up))
    data.frame(gene = genes, cluster = cl, avg_log2FC = lfc,
               p_val_adj = ifelse(seq_along(genes) %in% up, 1e-10, 0.8))
  }
  rbind(make("A", 1:30), make("B", 101:130))
}

synthetic_t2g <- function() {
  data.frame(
    gs_name     = rep(c("SET_A", "SET_B", "SET_RANDOM"), each = 20),
    gene_symbol = c(paste0("G", 1:20), paste0("G", 101:120), paste0("G", seq(41, 98, by = 3)))
  )
}

test_that("run_global_ora recovers the planted pathway per cluster", {
  out <- test_out_dir("ora")

  res <- quietly(run_global_ora(synthetic_de(), synthetic_t2g(), output_dir = out,
                                title_prefix = "Toy", indiv_width = 5, indiv_height = 4,
                                global_width = 5, global_height = 4))

  expect_s3_class(res, "data.frame")
  top <- res[res$Direction == "Activated", c("cluster", "ID")]
  expect_true(any(top$cluster == "A" & top$ID == "SET_A"))
  expect_true(any(top$cluster == "B" & top$ID == "SET_B"))
  expect_length(written_files(out, "^Toy_Summary_AllClusters_.*csv$"), 1)
  expect_length(written_files(out, "^Toy_A_Dotplot_"), 1)
})

test_that("run_global_gsea recovers the planted pathway and accepts paths without a slash", {
  out <- test_out_dir("gsea")     # no trailing slash on purpose

  res <- quietly(run_global_gsea(synthetic_de(), synthetic_t2g(), output_dir = out,
                                 title_prefix = "Toy", dotplot_width = 5, dotplot_height = 4,
                                 gseaplot_width = 5, gseaplot_height = 4))

  expect_s3_class(res, "data.frame")
  expect_true(any(res$cluster == "A" & res$ID == "SET_A" & res$NES > 0))
  expect_true(any(res$cluster == "B" & res$ID == "SET_B" & res$NES > 0))
  expect_length(written_files(out, "^Toy_GSEA_Cluster_Summary_.*csv$"), 1)
  expect_length(written_files(out, "^Toy_A_GSEA_NES_Dotplot_"), 1)
})

# Synthetic DESeq2-style results: 5 genes up, 5 down, the rest unchanged
synthetic_deseq <- function() {
  set.seed(4)
  genes <- paste0("G", 1:100)
  data.frame(
    baseMean = stats::runif(100, 10, 1000),
    log2FoldChange = c(seq(4, 2, length.out = 5), seq(-4, -2, length.out = 5), stats::rnorm(90, 0, 0.2)),
    padj = c(10^-(10:6), 10^-(10:6), stats::runif(90, 0.2, 1)),
    row.names = genes
  )
}

test_that("get_top_mixed_genes unions the most significant and largest-effect genes", {
  res <- synthetic_deseq()

  top <- get_top_mixed_genes(res, n_padj = 3, n_lfc = 3)

  expect_type(top, "character")
  expect_true(all(top %in% paste0("G", 1:10)))          # only strictly significant genes
  expect_lte(length(top), 6)
  expect_length(get_top_mixed_genes(res, n_padj = 50, n_lfc = 50), 10)
})

test_that("generate_volcano_trio returns a three-panel patchwork", {
  skip_if_not_installed("EnhancedVolcano")

  p <- quietly(generate_volcano_trio(synthetic_deseq(), main_title = "Toy",
                                     target_genes = c("G1", "G6", "G50")))

  expect_s3_class(p, "patchwork")
  expect_length(p$patches$plots, 2)    # 2 stored patches + the last plot = 3 panels
})

test_that("generate_and_save_heatmap skips comparisons without enough genes", {
  res <- synthetic_deseq()
  res$padj <- 1
  mat <- matrix(stats::rnorm(100 * 4), 100, dimnames = list(rownames(res), paste0("S", 1:4)))
  anno <- data.frame(group = c("a", "a", "b", "b"), row.names = colnames(mat))

  expect_null(quietly(generate_and_save_heatmap(res, "none", "none", mat, anno, colnames(mat),
                                                n_padj = 5, n_lfc = 5,
                                                out_dir = test_out_dir("heatmap_skip"), ts = "x")))
})
