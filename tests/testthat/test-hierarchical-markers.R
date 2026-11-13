# find_hierarchical_markers(): resolution tree, sibling markers, patterns and keys ------

# A small object with a known resolution tree:
#   res 0.1: A, B
#   res 0.5: A1, A2, B            (B passes through unchanged)
#   res 1  : A1, A2, B1, B2       (B splits; 20 A2-like cells are labelled B2 -> B2 purity 0.75)
# Planted markers: LA (all of A), LA1, LA2, LB (all of B), LB1, LB2; the rest is noise.
hierarchy_seurat <- function(cells_per_group = 60, n_noise = 200) {
  set.seed(7)
  groups <- rep(c("A1", "A2", "B1", "B2"), each = cells_per_group)
  marker_rate <- function(on) ifelse(on, 8, 0.05)
  counts <- rbind(
    LA  = stats::rpois(length(groups), marker_rate(groups %in% c("A1", "A2"))),
    LA1 = stats::rpois(length(groups), marker_rate(groups == "A1")),
    LA2 = stats::rpois(length(groups), marker_rate(groups == "A2")),
    LB  = stats::rpois(length(groups), marker_rate(groups %in% c("B1", "B2"))),
    LB1 = stats::rpois(length(groups), marker_rate(groups == "B1")),
    LB2 = stats::rpois(length(groups), marker_rate(groups == "B2")),
    matrix(stats::rpois(n_noise * length(groups), 2), n_noise,
           dimnames = list(paste0("NOISE", seq_len(n_noise)), NULL))
  )
  colnames(counts) <- paste0("cell", seq_along(groups))
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(methods::as(counts, "dgCMatrix")))
  obj <- suppressMessages(Seurat::NormalizeData(obj, verbose = FALSE))

  obj$RNA_snn_res.0.1 <- ifelse(groups %in% c("A1", "A2"), "A", "B")
  obj$RNA_snn_res.0.5 <- ifelse(groups %in% c("B1", "B2"), "B", groups)
  res1 <- groups
  res1[which(groups == "A2")[seq_len(20)]] <- "B2"     # mixed cluster: cells from two parents
  obj$RNA_snn_res.1 <- res1
  obj$true_group <- groups
  obj
}

hier_obj <- hierarchy_seurat()
hier_res <- suppressMessages(find_hierarchical_markers(hier_obj, pipeline = "RNA_snn", verbose = FALSE))

test_that("the tree links each cluster to the parent holding most of its cells", {
  tree <- hier_res$tree
  parent_of <- stats::setNames(tree$parent, tree$node)
  expect_equal(unname(parent_of[c("res0.1_A", "res0.5_A1", "res0.5_B", "res1_B2", "res1_A2")]),
               c("root", "res0.1_A", "res0.1_B", "res0.5_B", "res0.5_A2"))
  expect_equal(tree$purity[tree$node == "res1_B2"], 60 / 80)
  expect_equal(tree$node[tree$low_purity], "res1_B2")
  # B passes through res 0.5 unchanged, so the splits are root, A and B(res 0.5)
  expect_equal(hier_res$splits$split, c("root", "res0.1_A", "res0.5_B"))
  expect_true(all(hier_res$splits$status == "ok"))
})

test_that("markers compare siblings within the parent only", {
  m <- hier_res$markers
  expect_equal(m$gene[m$split == "root" & m$child == "res0.1_A" & m$rank == 1], "LA")
  expect_equal(m$gene[m$split == "res0.1_A" & m$child == "res0.5_A1" & m$rank == 1], "LA1")
  # LA is shared by A1 and A2, so it does not separate them
  expect_false("LA" %in% m$gene[m$split == "res0.1_A"])
  # cells that moved from A2 to B2 are not used in the B split: B2 is all real B2 cells
  b2_cells <- colnames(hier_obj)[hier_obj$RNA_snn_res.0.5 == "B" & hier_obj$RNA_snn_res.1 == "B2"]
  expect_equal(hier_res$splits$n_cells[hier_res$splits$split == "res0.5_B"], 60 + length(b2_cells))
})

test_that("every reported marker passes the filters and is ranked by AUC", {
  m <- hier_res$markers
  expect_true(all(m$p_val_adj < 0.05 & m$avg_log2FC >= 0.5 & m$pct.1 >= 0.5 & m$pct_diff >= 0.3))
  for (grp in split(m, paste(m$split, m$child))) {
    expect_equal(grp$rank, seq_len(nrow(grp)))
    expect_false(is.unsorted(rev(grp$auc)))
  }
})

test_that(".auc equals the Mann-Whitney U divided by n1 * n2 (ties count half)", {
  expect_equal(OHMmy:::.auc(c(2, 3), c(1, 1)), 1)
  expect_equal(OHMmy:::.auc(1, 1), 0.5)
  set.seed(3)
  x1 <- c(0, 0, stats::rpois(20, 3))
  x2 <- c(0, 0, 0, stats::rpois(15, 2))
  u <- suppressWarnings(stats::wilcox.test(x1, x2, exact = FALSE)$statistic)
  expect_equal(OHMmy:::.auc(x1, x2), unname(u) / (length(x1) * length(x2)))
})

test_that("patterns: lineage, specific, unstable, and the trajectory columns", {
  m <- hier_res$markers
  row <- function(split, child, gene) m[m$split == split & m$child == child & m$gene == gene, ]
  expect_equal(row("root", "res0.1_A", "LA")$pattern, "lineage")
  expect_false(row("root", "res0.1_A", "LA")$also_marks_children)
  expect_equal(row("root", "res0.1_B", "LB")$pattern, "lineage")      # B splits two levels down
  expect_equal(row("res0.1_A", "res0.5_A1", "LA1")$pattern, "specific")
  expect_true(is.na(row("res0.1_A", "res0.5_A1", "LA1")$also_marks_children))
  expect_gte(row("res0.1_A", "res0.5_A1", "LA1")$auc_vs_all, 0.8)
  expect_equal(row("res0.5_B", "res1_B2", "LB2")$pattern, "unstable")  # purity 0.75 < 0.8
  expect_true(all(is.na(m$auc_vs_all[m$child %in% c("res0.1_A", "res0.1_B")])))
  # plain columns (no dim attributes), so the tables print and save cleanly
  expect_null(dim(m$pattern))
  expect_type(m$pattern, "character")
  expect_null(dim(hier_res$key$pattern))
})

test_that("the key picks a planted marker per split, with backup genes in the same columns", {
  key <- hier_res$key
  key_genes <- key$gene[key$role == "key"]
  expect_true(key_genes[key$split[key$role == "key"] == "root"] %in% c("LA", "LB"))
  expect_true(key_genes[key$split[key$role == "key"] == "res0.1_A"] %in% c("LA1", "LA2"))
  expect_true(any(key$role == "backup"))
  backup <- key[key$role == "backup", ]
  expect_false(anyNA(backup[, c("auc", "rank", "p_val_adj", "pct.1", "pattern", "agreement")]))
  expect_true(all(is.na(key$agreement[key$role == "key"])))
  expect_true(all(grepl(" >= ", key$rule)))
  expect_true(all(hier_res$splits$cv_accuracy > 0.9))
  expect_equal(sum(diag(hier_res$confusion$res0.1_A)), sum(hier_res$confusion$res0.1_A))
})

test_that("the high side of a rule is the child with high expression", {
  key <- hier_res$key
  a_split <- key[key$split == "res0.1_A", ]
  expect_equal(a_split$high_side[a_split$gene == "LA1"], "res0.5_A1")
  expect_equal(a_split$high_side[a_split$gene == "LA2"], "res0.5_A2")
  # LA1 is only expressed in A cells, so as a root backup gene its high side is A
  root_la1 <- key[key$split == "root" & key$gene == "LA1", ]
  if (nrow(root_la1)) expect_equal(root_la1$high_side, "res0.1_A")
})

test_that("resolution columns can be given directly and repeated partitions are skipped", {
  obj <- hier_obj
  obj$RNA_snn_res.0.3 <- obj$RNA_snn_res.0.1          # same clusters as res 0.1
  expect_message(cols <- OHMmy:::.resolve_res_cols(obj@meta.data, NULL, "RNA_snn"), "skipped")
  expect_equal(cols, c("RNA_snn_res.0.1", "RNA_snn_res.0.5", "RNA_snn_res.1"))
  expect_equal(OHMmy:::.resolve_res_cols(obj@meta.data, c("RNA_snn_res.1", "RNA_snn_res.0.1"), NULL),
               c("RNA_snn_res.0.1", "RNA_snn_res.1"))     # ordered by resolution
})

test_that("small children are skipped and inputs are validated", {
  res <- suppressMessages(find_hierarchical_markers(hier_obj, pipeline = "RNA_snn", min_cells = 70,
                                                    verbose = FALSE))
  expect_match(res$splits$status[res$splits$split == "res0.1_A"], "not tested")
  expect_match(res$splits$skipped_children[res$splits$split == "res0.1_A"], "res0.5_A1")

  expect_error(find_hierarchical_markers(hier_obj), "exactly one")
  expect_error(find_hierarchical_markers(hier_obj, res_cols = "x", pipeline = "RNA_snn"), "exactly one")
  expect_error(find_hierarchical_markers(hier_obj, res_cols = "not_a_column"), "not found")
  expect_error(find_hierarchical_markers(hier_obj, pipeline = "nope"), "No resolution columns")
  expect_error(find_hierarchical_markers(hier_obj, pipeline = "RNA_snn", assay = "ADT"), "not found")
})

test_that("the caller's random state is left unchanged", {
  set.seed(11)
  expected <- stats::runif(1)
  set.seed(11)
  suppressMessages(find_hierarchical_markers(hier_obj, pipeline = "RNA_snn", verbose = FALSE))
  expect_equal(stats::runif(1), expected)
})

test_that("the whole key walks every cell from the root to a final cluster", {
  full <- hier_res$confusion_full
  expect_equal(sum(full), ncol(hier_obj))
  expect_equal(rownames(full), c("res1_A1", "res1_A2", "res1_B1", "res1_B2"))
  # the 20 A2 cells clustered with B2 are sorted by their expression into A2
  expect_gte(full["res1_B2", "res1_A2"], 18)
  expect_equal(hier_res$key_accuracy, sum(diag(full)) / sum(full))
  expect_gt(hier_res$key_accuracy, 0.85)
})

test_that("per-split confusion uses every cell of the parent, with moved cells in their own row", {
  # res 0.5 where 20 B1 cells are clustered with A1: they leave parent B for a child of A
  obj <- hier_obj
  res05 <- obj$true_group
  res05[which(obj$true_group == "B1")[seq_len(20)]] <- "A1"
  obj$RNA_snn_res.0.5 <- res05
  res <- suppressMessages(find_hierarchical_markers(obj, res_cols = c("RNA_snn_res.0.1", "RNA_snn_res.0.5"),
                                                    verbose = FALSE))
  b_split <- res$confusion[["res0.1_B"]]
  expect_equal(rownames(b_split), c("res0.5_B1", "res0.5_B2", "moved elsewhere"))
  expect_equal(sum(b_split["moved elsewhere", ]), 20)
  expect_equal(sum(b_split), sum(obj$RNA_snn_res.0.1 == "B"))
  row_b <- res$splits[res$splits$split == "res0.1_B", ]
  expect_lte(row_b$accuracy_all_parent_cells, 100 / 120)        # moved cells always count as wrong
  # the confusion CSV rows add up to the cells of each parent
  long <- OHMmy:::.split_confusion_long(res$confusion)
  expect_equal(as.vector(tapply(long$n_cells, long$split, sum)[c("root", "res0.1_A", "res0.1_B")]),
               c(240, 120, 120))
})

test_that("saved image sizes grow with content and stay within the limit", {
  sizes <- OHMmy:::.image_sizes(hier_res)
  expect_true(all(unlist(sizes) >= 4))
  big <- hier_res
  extra <- data.frame(node = paste0("res1_x", 1:200), layer = 3L, resolution = 1, res_col = "RNA_snn_res.1",
                      cluster = paste0("x", 1:200), parent = "res0.5_A1", n_cells = 1L, purity = 1,
                      n_children = 0L, low_purity = FALSE)
  big$tree <- rbind(big$tree, extra)
  big_sizes <- OHMmy:::.image_sizes(big)
  expect_gt(big_sizes$dotplot[["height"]], sizes$dotplot[["height"]])
  expect_gt(big_sizes$tree[["width"]], sizes$tree[["width"]])
  expect_lte(big_sizes$tree[["width"]], 100)
})

test_that("plots are returned and tables/plots are saved", {
  expect_s3_class(hier_res$plots$tree, "ggplot")
  expect_s3_class(hier_res$plots$dotplot, "ggplot")
  expect_s3_class(hier_res$plots$confusion, "ggplot")
  expect_s3_class(hier_res$plots$confusion_full, "ggplot")
  out <- test_out_dir("hier_markers")
  suppressMessages(find_hierarchical_markers(hier_obj, pipeline = "RNA_snn", sample_name = "toy",
                                             output_dir = out, verbose = FALSE))
  expect_length(written_files(out, "^toy_HierMarkers_.*\\.csv$"), 6)
  expect_length(written_files(out, "^toy_HierMarkers_.*\\.jpg$"), 4)
  conf_csv <- utils::read.csv(written_files(out, "confusion_splits\\.csv$"))
  expect_named(conf_csv, c("split", "true_child", "predicted_child", "n_cells", "share_of_true"))
})
