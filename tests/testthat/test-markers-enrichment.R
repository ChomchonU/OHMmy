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

test_that("FindTopMarkersAndHeatmap(use_sct = TRUE) works on a subset SCT object", {
  skip_if_not_installed("harmony")
  sct <- quietly(ProcessSeuratSCT(small_seurat(), batch_col = "batch", vars_to_regress = NULL,
                                  reduction_name = "pca.sct", integration_method = "HarmonyIntegration",
                                  integration_reduction = "harmony.sct", dims = 1:10, verbose = FALSE))
  sct <- quietly(Seurat::PrepSCTFindMarkers(sct))
  # Removing the deepest cells lowers every model's observed median UMI, the case in
  # which PrepSCTFindMarkers() skips re-correction and FindMarkers() would refuse to run
  sub <- subset(sct, cells = colnames(sct)[sct$nCount_RNA < stats::quantile(sct$nCount_RNA, 0.6)])
  SeuratObject::Idents(sub) <- "cell_type"
  expect_false(OHMmy:::.sct_models_consistent(quietly(Seurat::PrepSCTFindMarkers(sub))))

  expect_message(
    res <- suppressWarnings(FindTopMarkersAndHeatmap(sub, sample_name = "sct_subset", use_sct = TRUE,
                                                     output_dir_base = test_out_dir("markers_sct"),
                                                     width = 4, height = 4, dpi = 30)),
    "recorrect_umi = FALSE"
  )
  expect_gt(nrow(res$markers), 0)
  expect_s3_class(res$heatmap, "ggplot")
})

test_that("FindTopMarkersAndHeatmap warns instead of failing when no markers are found", {
  local_mocked_bindings(FindAllMarkers = function(...) data.frame())
  obj <- small_seurat()
  SeuratObject::Idents(obj) <- "cell_type"
  out <- test_out_dir("markers_empty")

  expect_warning(res <- quietly_messages(FindTopMarkersAndHeatmap(obj, sample_name = "empty",
                                                                  output_dir_base = out)),
                 "No markers were returned")
  expect_null(res$heatmap)
  expect_equal(nrow(res$markers), 0)
  expect_length(written_files(out, recursive = TRUE), 0)
})

# Edge cases of the two marker-table plots ---------------------------------------

toy_marker_table <- function() {
  set.seed(1)
  df <- expand.grid(gene = c("GZMB", "NKG7", "LYZ", "MS4A1"), cluster = c(0, 1, 2),
                    cond = c("A", "B", "C"), stringsAsFactors = FALSE)
  df$avg_log2FC <- stats::runif(nrow(df), 0.3, 3)
  df$diff <- stats::runif(nrow(df), 0, 0.5)
  df
}

capture_saved_plots <- function(expr) {
  plots <- list()
  testthat::with_mocked_bindings(
    force(expr),
    ggsave = function(filename, plot, ...) { plots[[basename(filename)]] <<- plot; invisible(filename) },
    .package = "OHMmy"
  )
  plots
}

test_that("split dot plots order genes and label clusters when cluster names are not numbers", {
  df <- toy_marker_table()
  df$cluster <- c("CD4 T", "NK", "B")[match(df$cluster, 0:2)]

  plots <- capture_saved_plots(expect_no_warning(suppressMessages(
    plot_split_dotplots_by_gene_cluster(df, id_col = "cond", n_splits = 1, output_dir = tempfile()))))
  labels <- levels(plots[[grep("^DotPlot_Part1", names(plots))]]$data$gene_cluster)

  expect_true(all(grepl("^(GZMB|NKG7|LYZ|MS4A1) \\((CD4 T|NK|B)\\)$", labels)))   # no "C" prefix
  gene_means <- sort(tapply(df$avg_log2FC, df$gene, mean))
  expect_equal(unique(sub(" \\(.*$", "", labels)), names(gene_means))                # ordered by mean log2FC
})

test_that("split dot plots never write more (or empty) plots than gene x cluster combinations", {
  df <- toy_marker_table()
  df <- df[df$gene == "GZMB" & df$cluster %in% 0:1, ]
  out <- test_out_dir("split_small")

  expect_message(plot_split_dotplots_by_gene_cluster(df, id_col = "cond", n_splits = 4, output_dir = out),
                 "instead of 4")
  expect_length(written_files(out, "^DotPlot_Part"), 2)
})

test_that("split dot plots handle one condition, duplicated rows and missing ids", {
  df <- toy_marker_table()
  out <- test_out_dir("split_edge")

  expect_message(plot_split_dotplots_by_gene_cluster(df[df$cond == "A", ], id_col = "cond", output_dir = out),
                 "Only one condition")
  expect_length(written_files(out, "^Dendrogram_ALL_"), 0)
  expect_gt(length(written_files(out, "^DotPlot_Part")), 0)

  dup <- rbind(df, transform(df, avg_log2FC = avg_log2FC + 1))
  plots <- capture_saved_plots(expect_message(
    plot_split_dotplots_by_gene_cluster(dup, id_col = "cond", n_splits = 1, output_dir = tempfile()),
    "Averaged"))
  part <- plots[[grep("^DotPlot_Part1", names(plots))]]$data
  expect_equal(nrow(part), nrow(df))                                       # one row per combination
  expect_equal(sort(part$avg_log2FC), sort(df$avg_log2FC + 0.5))           # duplicates averaged

  na_id <- df; na_id$cond[1:3] <- NA
  expect_message(plot_split_dotplots_by_gene_cluster(na_id, id_col = "cond", output_dir = tempfile()),
                 "Dropping 3 row")
  expect_error(plot_split_dotplots_by_gene_cluster(df[, names(df) != "diff"], id_col = "cond"),
               "missing: diff")
})

test_that("marker dendro plot handles no matches, one condition, missing ids and bad input", {
  df <- toy_marker_table()
  out <- test_out_dir("marker_edge")

  expect_message(res <- plot_gene_markers_with_dendro(df, genes = "NOT_A_GENE", id_col = "cond", output_dir = out),
                 "No genes matched")
  expect_null(res)

  expect_message(res <- plot_gene_markers_with_dendro(df[df$cond == "A", ], genes = c("GZMB", "LYZ"),
                                                      id_col = "cond", output_dir = out),
                 "Only one condition")
  expect_s3_class(res$dotplot, "ggplot")
  expect_null(res$dendrogram)
  expect_null(res$combined)

  na_id <- df; na_id$cond[1:3] <- NA
  expect_message(plot_gene_markers_with_dendro(na_id, genes = c("GZMB", "LYZ"), id_col = "cond", output_dir = out),
                 "Dropping 3 row")
  expect_error(plot_gene_markers_with_dendro(df[, names(df) != "diff"], genes = "LYZ", id_col = "cond"),
               "missing: diff")
  expect_error(plot_gene_markers_with_dendro(df, genes = "LYZ", id_col = "cond", dendro_side = "top"),
               "should be one of")
})

test_that("marker dendro plot treats a single non-gene string as a regular expression", {
  df <- toy_marker_table()
  out <- test_out_dir("marker_regex")
  genes_of <- function(genes) {
    levels(suppressMessages(plot_gene_markers_with_dendro(df, genes = genes, id_col = "cond",
                                                          output_dir = out))$dotplot$data$gene)
  }
  expect_setequal(genes_of("GZM|NKG"), c("GZMB", "NKG7"))       # regex with |
  expect_equal(genes_of("GZM"), "GZMB")                          # plain prefix
  expect_setequal(genes_of("[GL][ZY]"), c("GZMB", "LYZ"))        # character class
  expect_equal(genes_of("LYZ"), "LYZ")                           # exact gene name stays exact
})

test_that("dendrogram leaves line up with the dot-plot rows in the combined figure", {
  res <- suppressMessages(plot_gene_markers_with_dendro(toy_marker_table(), genes = c("GZMB", "LYZ", "NKG7"),
                                                        id_col = "cond", output_dir = tempfile()))
  dot_range <- ggplot2::ggplot_build(res$dotplot)$layout$panel_params[[1]]$y.range
  den_range <- ggplot2::ggplot_build(res$dendrogram)$layout$panel_params[[1]]$y.range
  expect_equal(den_range, dot_range)
})

# -- Analysis-quality options and returned plots ---------------------------------------

test_that("run_global_ora uses `universe` as the background and returns its plots", {
  de <- synthetic_de()
  # background genes that are annotated (enricher ignores unannotated universe genes)
  t2g <- rbind(synthetic_t2g(), data.frame(gs_name = "SET_BG", gene_symbol = paste0("BG", 1:200)))
  run <- function(...) quietly(run_global_ora(de, t2g, output_dir = test_out_dir("ora_universe"),
                                              title_prefix = "Toy", indiv_width = 5, indiv_height = 4,
                                              global_width = 5, global_height = 4, ...))
  default <- run()
  same    <- run(universe = unique(de$gene))
  wider   <- run(universe = c(unique(de$gene), paste0("BG", 1:200)))

  bg_size <- function(r) as.numeric(sub(".*/", "", r$BgRatio[r$cluster == "A" & r$ID == "SET_A"]))
  expect_equal(default$pvalue, same$pvalue)                 # NULL = genes of ORA_df
  expect_equal(bg_size(wider), bg_size(default) + 200)      # the universe is the background
  expect_s3_class(attr(default, "plots")$A_dotplot, "ggplot")
  expect_s3_class(attr(default, "plots")$global_activated, "ggplot")

  expect_warning(
    run_global_ora(de, t2g, output_dir = test_out_dir("ora_universe2"), universe = paste0("G", 1:150),
                   indiv_width = 5, indiv_height = 4, global_width = 5, global_height = 4) |>
      utils::capture.output() |> suppressMessages(),
    "not in `universe`"
  )
})

test_that("run_global_gsea can rank by a signed significance score", {
  de <- synthetic_de()
  de$p_val <- ifelse(de$p_val_adj < 0.05, 1e-12, 0.5) * stats::runif(nrow(de), 0.5, 1)
  run <- function(...) quietly(run_global_gsea(de, synthetic_t2g(), output_dir = test_out_dir("gsea_rank"),
                                               title_prefix = "Toy", dotplot_width = 5, dotplot_height = 4,
                                               gseaplot_width = 5, gseaplot_height = 4, ...))
  res <- run(rank_by = "signed_p")

  expect_true(any(res$cluster == "A" & res$ID == "SET_A" & res$NES > 0))
  expect_s3_class(attr(res, "plots")$A_nes_dotplot, "ggplot")
  expect_error(run_global_gsea(de[, c("gene", "cluster", "avg_log2FC")], synthetic_t2g(),
                               output_dir = test_out_dir("gsea_rank2"), rank_by = "signed_p"),
               "p_val")
  expect_error(run(rank_by = "nope"), "should be one of")
})

test_that("FindTopMarkersAndHeatmap scales top markers missing from scale.data", {
  obj <- small_seurat()
  SeuratObject::Idents(obj) <- "cell_type"
  obj <- suppressMessages(Seurat::ScaleData(obj, features = rownames(obj)[1:5], verbose = FALSE))

  expect_message(
    res <- quietly_output(FindTopMarkersAndHeatmap(obj, sample_name = "s", top_n = 3,
                                                   output_dir_base = test_out_dir("markers_scale"),
                                                   width = 5, height = 5, dpi = 30, add_timestamp = FALSE)),
    "missing from scale.data"
  )
  expect_true(all(unique(res$top_markers$gene) %in% as.character(res$heatmap$data$Feature)))
})

test_that("plot_split_dotplots_by_gene_cluster returns its plots", {
  df <- multi_marker_table()
  res <- quietly(plot_split_dotplots_by_gene_cluster(df, gene_regex = "^HLA-D|^CD79", id_col = "batch",
                                                     output_dir = test_out_dir("split_return"),
                                                     n_splits = 2))
  expect_named(res, c("dendrogram", "dotplots"))
  expect_s3_class(res$dendrogram, "ggplot")
  expect_length(res$dotplots, 2)
  expect_s3_class(res$dotplots$Part1, "ggplot")
})

test_that("generate_and_save_heatmap attaches the drawn heatmap", {
  res <- synthetic_deseq()
  set.seed(2)
  mat <- matrix(stats::rnorm(100 * 4), 100, dimnames = list(rownames(res), paste0("S", 1:4)))
  anno <- data.frame(group = c("a", "a", "b", "b"), row.names = colnames(mat))
  f <- quietly(generate_and_save_heatmap(res, "ab", "a vs b", mat, anno, colnames(mat), n_padj = 5,
                                         n_lfc = 5, out_dir = test_out_dir("heatmap_attr"), ts = "x"))
  expect_true(file.exists(f))
  expect_s3_class(attr(f, "plot"), "gtable")
})
