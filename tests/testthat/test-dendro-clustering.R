# Cutting the dendrograms of plot_dot_dendro_multi() ---------------------------------

# A dot_df-like table built from 2-D points (2 "genes"): two tight groups close
# together (A, B) and one loose group far away (C). No single horizontal cut can
# separate A, B and C: C's internal spread is larger than the A-B distance.
variable_height_df <- function(seed = 1, c_center = 5, c_sd = 0.6) {
  set.seed(seed)
  ids <- c(paste0("A", 1:6), paste0("B", 1:6), paste0("C", 1:6))
  pts <- rbind(cbind(stats::rnorm(6, 0, 0.05), stats::rnorm(6, 0, 0.05)),
               cbind(stats::rnorm(6, 1, 0.05), stats::rnorm(6, 0, 0.05)),
               cbind(stats::rnorm(6, c_center, c_sd), stats::rnorm(6, 0, c_sd)))
  data.frame(cluster_var = rep(paste0(ids, "_grp"), 2), id = rep(ids, 2), var = "grp",
             features.plot = rep(c("G1", "G2"), each = 18), avg.exp.scaled = c(pts[, 1], pts[, 2]))
}

truth_of <- function(clusters) {
  rows <- clusters[clusters$type == "rows", ]
  stats::setNames(as.character(rows$cluster), substr(rows$label, 1, 1))
}

recovers_abc <- function(clusters) {
  lab <- truth_of(clusters)
  !anyNA(lab) && length(unique(lab)) == 3 &&
    all(tapply(lab, names(lab), function(x) length(unique(x))) == 1)
}

small_dot_df <- function() {
  quietly(plot_dot_dendro_multi(small_seurat(), meta_cols = c("cell_type", "batch", "donor"),
                                feature_df = small_feature_df(), prefix = "Test", pct_threshold = 0,
                                output_dir = test_out_dir("dot_multi_for_clusters"),
                                save_options = "gene_list"))
}

test_that("plot_dot_dendro_multi attaches the trees behind its dendrograms", {
  df <- small_dot_df()
  expect_true(all(c("dendro_matrix", "row_hclust", "col_hclust") %in% names(attributes(df))))
  expect_s3_class(attr(df, "row_hclust"), "hclust")
  expect_setequal(attr(df, "row_hclust")$labels, unique(as.character(df$cluster_var)))
})

test_that("trees rebuilt from the data frame match the attached ones exactly", {
  df <- small_dot_df()
  plain <- df
  attributes(plain)[c("dendro_matrix", "row_hclust", "col_hclust")] <- NULL

  attached <- OHMmy:::.dot_dendro_trees(df)
  rebuilt  <- OHMmy:::.dot_dendro_trees(plain)
  for (tp in c("rows", "genes")) {
    expect_identical(rebuilt$hclust[[tp]]$merge, attached$hclust[[tp]]$merge)
    expect_equal(rebuilt$hclust[[tp]]$height, attached$hclust[[tp]]$height)
    expect_identical(rebuilt$hclust[[tp]]$labels, attached$hclust[[tp]]$labels)
  }
})

test_that("the inconsistency method recovers clusters at different heights where a fixed cut cannot", {
  df <- variable_height_df()
  hc <- OHMmy:::.dot_dendro_trees(df)$hclust$rows
  truth <- substr(hc$labels, 1, 1)

  fixed_can <- vapply(2:8, function(k) {
    lab <- stats::cutree(hc, k)
    length(unique(lab)) == 3 && all(tapply(lab, truth, function(x) length(unique(x))) == 1)
  }, logical(1))
  expect_false(any(fixed_can))

  res <- cluster_dot_dendro(df)                                   # default: inconsistency
  expect_equal(res$parameters$method, "inconsistency")
  expect_true(recovers_abc(res$clusters))
  expect_equal(unname(res$k["rows"]), 3L)

  # the tight clusters are cut lower than the loose one
  h <- tapply(res$clusters$cut_height, truth_of(res$clusters), unique)
  lab <- truth_of(res$clusters)
  expect_lt(h[[lab[["A"]]]], h[[lab[["C"]]]])
  expect_lt(h[[lab[["B"]]]], h[[lab[["C"]]]])
})

test_that("cut heights lie between each cluster's top merge and the merge above it", {
  res <- cluster_dot_dendro(variable_height_df())
  hc <- res$hclust$rows
  cl <- res$clusters
  for (c_id in levels(cl$cluster)) {
    sub <- cl[cl$cluster %in% c_id, ]
    expect_true(all(sub$cut_height >= sub$cluster_height))
    expect_true(all(sub$cut_height <= max(hc$height)))
  }
})

test_that("auto settings pick the candidate with the best silhouette", {
  skip_if_not_installed("dynamicTreeCut")
  df <- variable_height_df(c_center = 3, c_sd = 0.3)

  dyn <- cluster_dot_dendro(df, method = "dynamic")
  sel <- dyn$selection$rows
  expect_named(sel, c("deepSplit", "n_clusters", "n_unassigned", "avg_silhouette", "selected"))
  expect_equal(sel$deepSplit[sel$selected], unname(dyn$parameters$setting$rows["deepSplit"]))
  expect_equal(sel$avg_silhouette[sel$selected], max(sel$avg_silhouette, na.rm = TRUE))
  expect_true(recovers_abc(dyn$clusters))

  inc <- cluster_dot_dendro(df, method = "inconsistency")
  sel_inc <- inc$selection$rows
  expect_equal(sel_inc$avg_silhouette[sel_inc$selected], max(sel_inc$avg_silhouette, na.rm = TRUE))
})

test_that("explicit deepSplit and cutoff values are used as given", {
  skip_if_not_installed("dynamicTreeCut")
  df <- variable_height_df()
  expect_null(cluster_dot_dendro(df, method = "dynamic", deepSplit = 1)$selection$rows)
  res <- cluster_dot_dendro(df, method = "inconsistency", cutoff = 100)   # nothing inconsistent
  expect_equal(unname(res$k["rows"]), 1L)
})

test_that("unclear leaves stay unassigned by default and are assigned with assign_all = TRUE", {
  df <- variable_height_df()
  # add one isolated leaf far from everything
  extra <- data.frame(cluster_var = "Z1_grp", id = "Z1", var = "grp", features.plot = c("G1", "G2"),
                      avg.exp.scaled = c(30, 30))
  df <- rbind(df, extra)

  res <- cluster_dot_dendro(df)
  z <- res$clusters[res$clusters$label == "Z1_grp", ]
  expect_true(is.na(z$cluster))
  expect_true(is.na(z$cut_height))

  res_all <- cluster_dot_dendro(df, assign_all = TRUE)
  expect_false(anyNA(res_all$clusters$cluster))
})

test_that("minClusterSize leaves small groups unassigned", {
  res <- cluster_dot_dendro(variable_height_df(), minClusterSize = 7)   # every true group has 6 leaves
  expect_true(all(is.na(res$clusters$cluster)) || unname(res$k["rows"]) < 3)
})

test_that("the fixed method cuts at k clusters or at height h", {
  df <- variable_height_df()
  by_k <- cluster_dot_dendro(df, method = "fixed", k = 2)
  expect_equal(unname(by_k$k["rows"]), 2L)
  expect_false(anyNA(by_k$clusters$cluster))

  hc <- by_k$hclust$rows
  by_h <- cluster_dot_dendro(df, method = "fixed", h = mean(range(hc$height)))
  expect_equal(unname(by_h$k["rows"]), length(unique(stats::cutree(hc, h = mean(range(hc$height))))))
})

test_that("genes, both trees, plots and saved files are returned", {
  df <- small_dot_df()
  out <- test_out_dir("dendro_clusters")

  res <- cluster_dot_dendro(df, what = "both", output_dir = out)

  expect_named(res, c("clusters", "k", "selection", "plot", "hclust", "parameters"))
  expect_setequal(unique(res$clusters$type), c("rows", "genes"))
  expect_named(res$k, c("rows", "genes"))
  expect_true(all(vapply(res$plot, inherits, logical(1), what = "ggplot")))
  expect_true(all(is.na(res$clusters$id[res$clusters$type == "genes"])))
  expect_false(anyNA(res$clusters$var[res$clusters$type == "rows"]))
  # every leaf appears once per tree, in dendrogram order
  rows <- res$clusters[res$clusters$type == "rows", ]
  expect_equal(rows$dendro_order, seq_len(nrow(rows)))
  expect_length(written_files(out, "\\.csv$"), 1)
  expect_length(written_files(out, "\\.jpg$"), 2)
})

test_that("row clusters are written back to the cells of the chosen metadata column", {
  obj <- small_seurat()
  df <- small_dot_df()

  res <- cluster_dot_dendro(df, method = "fixed", k = 3, seurat_obj = obj, meta_col = "cell_type",
                            new_col = "group_cluster")

  rows <- res$clusters[res$clusters$type == "rows" & res$clusters$var == "cell_type", ]
  expected <- as.character(rows$cluster)[match(obj$cell_type, rows$id)]
  expect_equal(as.character(res$seurat_obj$group_cluster), expected)
  expect_true(is.factor(res$seurat_obj$group_cluster))
})

test_that("cluster_dot_dendro validates its inputs", {
  df <- small_dot_df()
  obj <- small_seurat()
  expect_error(cluster_dot_dendro(data.frame(a = 1)), "plot_dot_dendro_multi")
  expect_error(cluster_dot_dendro(df, method = "fixed"), "exactly one of")
  expect_error(cluster_dot_dendro(df, method = "fixed", k = 2, h = 1), "exactly one of")
  expect_error(cluster_dot_dendro(df, deepSplit = 7), "deepSplit")
  expect_error(cluster_dot_dendro(df, seurat_obj = obj), "meta_col")
  expect_error(cluster_dot_dendro(df, what = "genes", seurat_obj = obj, meta_col = "cell_type"), "row clusters")
  expect_error(cluster_dot_dendro(df, seurat_obj = obj, meta_col = "not_used"), "must be one of")
})
