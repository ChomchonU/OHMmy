# track_cluster_markers(): gene paths, cluster profiles and suggested labels ------------
# Fixture: hierarchy_seurat() in helper-fixtures.R (A -> A1, A2; B -> B1, B2 with a
# mixed B2 of purity 0.75; planted markers LA, LA1, LA2, LB, LB1, LB2).

track_obj <- hierarchy_seurat()
track_hier <- suppressMessages(find_hierarchical_markers(track_obj, pipeline = "RNA_snn", verbose = FALSE))
track_res <- suppressMessages(track_cluster_markers(track_obj, hier = track_hier, verbose = FALSE))

test_that("gene paths: where each marker starts and where it stays on", {
  p <- track_res$gene_paths
  la <- p[p$gene == "LA", ]
  expect_equal(nrow(la), 1)
  expect_equal(la$start_node, "res0.1_A")
  expect_equal(la$start_resolution, 0.1)
  expect_equal(la$kind, "lineage")
  expect_equal(c(la$n_final, la$n_final_on), c(2, 2))     # on in both final A clusters
  expect_equal(la$stops_at, "")
  # LA1 is in half of all A cells, so it already marks A against B; it is marked again at
  # A1 and lost in A2
  la1 <- p[p$gene == "LA1", ]
  expect_equal(la1$start_node, "res0.1_A")
  expect_equal(la1$marked_again_below, "res0.5_A1")
  expect_match(la1$stops_at, "res0.5_A2")
  expect_equal(c(la1$n_final, la1$n_final_on), c(2, 1))
})

test_that("node profile: AUCs against all cells and against siblings", {
  prof <- track_res$node_profile
  la_a <- prof[prof$gene == "LA" & prof$node == "res0.1_A", ]
  expr <- SeuratObject::LayerData(track_obj, layer = "data")
  in_a <- track_obj$RNA_snn_res.0.1 == "A"
  expect_equal(la_a$auc_vs_all, OHMmy:::.auc(expr["LA", in_a], expr["LA", !in_a]))
  expect_true(la_a$is_marker)
  expect_true(is.na(prof$auc_vs_siblings[prof$node == "res0.5_B"][1]))   # passes through: no siblings
  expect_equal(prof$global, prof$auc_vs_all >= 0.7)
  expect_equal(nrow(prof), (nrow(track_hier$tree) - 1) * length(unique(prof$gene)))
})

test_that("labels: path of top markers, suggested or ambiguous with a reason", {
  a <- track_res$annotation
  expect_equal(a$path_label[a$node == "res1_A1"], "LA+ > LA1+")
  expect_equal(a$status[a$node == "res1_A1"], "suggested")
  expect_equal(a$path_label[a$node == "res1_B2"], "LB+ > LB2+")
  expect_equal(a$status[a$node == "res1_B2"], "ambiguous")
  expect_match(a$reason[a$node == "res1_B2"], "low purity")
  expect_match(a$lineage_markers[a$node == "res1_A1"], "LA")
  expect_match(a$own_markers[a$node == "res1_A1"], "LA1")
})

test_that("a label made of relative markers is ambiguous", {
  strict <- suppressMessages(track_cluster_markers(track_obj, hier = track_hier, min_global_auc = 1.01,   # nothing is global
                                                   verbose = FALSE))
  expect_true(all(strict$annotation$status == "ambiguous"))
  expect_match(strict$annotation$reason[strict$annotation$node == "res1_A1"], "only relative")
})

test_that("reusing a find_hierarchical_markers() result gives the same as running it", {
  standalone <- suppressMessages(track_cluster_markers(track_obj, pipeline = "RNA_snn", verbose = FALSE))
  expect_equal(standalone$annotation, track_res$annotation)
  expect_equal(standalone$gene_paths, track_res$gene_paths)
  expect_error(track_cluster_markers(track_obj[, 1:200], hier = track_hier, verbose = FALSE), "does not match")
  expect_warning(suppressMessages(track_cluster_markers(track_obj, hier = track_hier, pipeline = "RNA_snn",
                                                        verbose = FALSE)), "ignored")
})

test_that("plots are returned and tables/plots are saved", {
  expect_s3_class(track_res$plots$tree, "ggplot")
  expect_s3_class(track_res$plots$heatmap, "ggplot")
  expect_s3_class(track_res$plots$dotplot, "ggplot")
  out <- test_out_dir("track_markers")
  suppressMessages(track_cluster_markers(track_obj, hier = track_hier, sample_name = "toy",
                                         output_dir = out, verbose = FALSE))
  expect_length(written_files(out, "^toy_TrackMarkers_.*\\.csv$"), 3)
  expect_length(written_files(out, "^toy_TrackMarkers_.*\\.jpg$"), 3)
})
