test_that("CleanSeuratReductions fixes bad naming conventions", {
  library(Seurat)
  dummy <- pbmc_small

  # Intentionally break a reduction name
  dummy[["pca.bad_name"]] <- dummy[["pca"]]

  # Ensure it was broken
  expect_true("pca.bad_name" %in% Reductions(dummy))

  clean_dummy <- CleanSeuratReductions(dummy)

  # Ensure old name is gone and new standard name is present
  expect_false("pca.bad_name" %in% Reductions(clean_dummy))
  expect_true("pcaBadName" %in% Reductions(clean_dummy))
})

test_that("ClusterAndUMAP executes resolution sweeping", {
  library(Seurat)
  dummy <- pbmc_small

  temp_plot_dir <- file.path(tempdir(), "clustree_plots")

  suppressWarnings({
    res <- ClusterAndUMAP(
      seurat_obj = dummy,
      sample_name = "TestSample",
      dims = 1:3, # Tiny dims for tiny pbmc_small dataset
      reduction = "pca",
      umap_name = "umap.test",
      cluster_resolutions = c(0.2, 0.4),
      final_resolution = 0.4,
      plot_dir = temp_plot_dir
    )
  })

  # Check output structure
  expect_s4_class(res$seurat, "Seurat")
  expect_s3_class(res$clustree, "ggplot")

  # Check that UMAP was built
  expect_true("umap.test" %in% Reductions(res$seurat))

  # Check that metadata columns were created
  expect_true("RNA_snn_res.0.2" %in% colnames(res$seurat@meta.data))
  expect_true("RNA_snn_res.0.4" %in% colnames(res$seurat@meta.data))

  # Verify plot file exists
  expect_true(file.exists(res$plot_file))
})
