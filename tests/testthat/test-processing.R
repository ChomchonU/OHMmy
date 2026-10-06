# Normalisation, integration and clustering -------------------------------------

test_that("CleanSeuratReductions fixes bad naming conventions", {
  dummy <- SeuratObject::pbmc_small
  suppressWarnings(dummy[["pca.bad_name"]] <- dummy[["pca"]])
  expect_true("pca.bad_name" %in% SeuratObject::Reductions(dummy))

  clean_dummy <- suppressWarnings(quietly(CleanSeuratReductions(dummy)))

  expect_false("pca.bad_name" %in% SeuratObject::Reductions(clean_dummy))
  expect_true("pcaBadName" %in% SeuratObject::Reductions(clean_dummy))
  expect_equal(SeuratObject::Key(clean_dummy[["pcaBadName"]]), "pcaBadName_")
  # already well-formed names are untouched
  expect_true(all(c("pca", "tsne") %in% SeuratObject::Reductions(clean_dummy)))
})

test_that("ClusterAndUMAP sweeps resolutions, saves a clustree and builds a UMAP", {
  dummy <- SeuratObject::pbmc_small
  plot_dir <- test_out_dir("clustree")

  res <- quietly(ClusterAndUMAP(
    seurat_obj = dummy,
    sample_name = "TestSample",
    dims = 1:3,
    reduction = "pca",
    umap_name = "umap.test",
    cluster_resolutions = c(0.2, 0.4),
    final_resolution = 0.4,
    plot_dir = plot_dir,
    dpi = 50
  ))

  expect_named(res, c("seurat", "clustree", "plot_file"))
  expect_s4_class(res$seurat, "Seurat")
  expect_s3_class(res$clustree, "ggplot")
  expect_true("umap.test" %in% SeuratObject::Reductions(res$seurat))
  expect_true(all(c("RNA_snn_res.0.2", "RNA_snn_res.0.4") %in% colnames(res$seurat@meta.data)))
  expect_true(file.exists(res$plot_file))
})

test_that("ClusterAndUMAP reuses existing graphs and saves the object when asked", {
  dummy <- SeuratObject::pbmc_small
  out <- test_out_dir("clustree_reuse")
  rds <- file.path(out, "obj.rds")

  first <- quietly(ClusterAndUMAP(dummy, dims = 1:3, reduction = "pca", umap_name = "umap.test",
                                  cluster_resolutions = 0.5, final_resolution = 0.5,
                                  plot_dir = out, dpi = 50))
  expect_message(
    ClusterAndUMAP(first$seurat, dims = 1:3, reduction = "pca", umap_name = "umap.test",
                   cluster_resolutions = 0.5, final_resolution = 0.5, plot_dir = out,
                   dpi = 50, save_path = rds),
    "already exists"
  )
  expect_true(file.exists(rds))
})

test_that("ProcessSeuratLOG normalises, integrates with Harmony and re-joins layers", {
  skip_if_not_installed("harmony")
  obj <- small_seurat()
  elbow_dir <- test_out_dir("elbow_log")

  res <- quietly(ProcessSeuratLOG(
    obj,
    batch_col = "batch",
    vars_to_regress = "pct_counts_mt",
    reduction_name = "pca.log",
    integration_method = "HarmonyIntegration",
    integration_reduction = "harmony.log",
    dims = 1:10,
    elbow_plot_dir = elbow_dir,
    sample_name = "small",
    verbose = FALSE
  ))

  expect_s4_class(res, "Seurat")
  expect_true(all(c("pca.log", "harmony.log") %in% SeuratObject::Reductions(res)))
  expect_equal(SeuratObject::DefaultAssay(res), "RNA")
  expect_setequal(SeuratObject::Layers(res), c("counts", "data", "scale.data"))
  expect_length(written_files(elbow_dir, "^elbow_plot_small_"), 1)
})

test_that("ProcessSeuratLOG removes receptor genes and validates integration_method", {
  obj <- small_seurat()
  expect_error(
    quietly(ProcessSeuratLOG(obj, reduction_name = "pca.log", integration_method = "NotAMethod",
                             verbose = FALSE)),
    "Unknown integration method"
  )
  expect_error(
    quietly(ProcessSeuratLOG(obj, reduction_name = "pca.log", integration_method = Seurat::CCAIntegration,
                             verbose = FALSE)),
    "must be a single string"
  )
})

test_that("ProcessSeuratSCT normalises with SCTransform and integrates", {
  skip_if_not_installed("harmony")
  obj <- small_seurat()

  res <- quietly(ProcessSeuratSCT(
    obj,
    batch_col = "batch",
    vars_to_regress = NULL,
    reduction_name = "pca.sct",
    integration_method = "HarmonyIntegration",
    integration_reduction = "harmony.sct",
    dims = 1:10,
    verbose = FALSE
  ))

  expect_s4_class(res, "Seurat")
  expect_equal(SeuratObject::DefaultAssay(res), "SCT")
  expect_true(all(c("pca.sct", "harmony.sct") %in% SeuratObject::Reductions(res)))
})

# Internal helpers ---------------------------------------------------------------

test_that(".resolve_k_params shrinks k values to the smallest batch", {
  k <- quietly(OHMmy:::.resolve_k_params(NULL, NULL, NULL, NULL,
                                         interactive_mode = FALSE, min_batch_cells = 40))
  expect_equal(k, list(k.weight = 39, k.anchor = 5, k.filter = 39, k.score = 30))

  manual <- OHMmy:::.resolve_k_params(20, 3, 25, 10, interactive_mode = FALSE, min_batch_cells = 1000)
  expect_equal(manual, list(k.weight = 20, k.anchor = 3, k.filter = 25, k.score = 10))
})

test_that(".choose_dims caps dims at the number of available PCs", {
  expect_equal(OHMmy:::.choose_dims(1:10, interactive_mode = FALSE, max_pca = 30), 1:10)
  expect_warning(capped <- OHMmy:::.choose_dims(1:50, interactive_mode = FALSE, max_pca = 20),
                 "Adjusting down")
  expect_equal(capped, 1:20)
})

test_that(".remove_receptor_genes and .boost_scaled_genes work on variable features", {
  obj <- small_seurat()
  SeuratObject::VariableFeatures(obj) <- c("TRAV1-FAKE", "MS4A1", "CD79A")
  obj <- suppressWarnings(SeuratObject::SetAssayData(obj, layer = "scale.data",
          new.data = SeuratObject::GetAssayData(obj, layer = "scale.data")))

  cleaned <- quietly(OHMmy:::.remove_receptor_genes(obj, "^TR[ABDG]V|^IG[HKL]V"))
  expect_setequal(SeuratObject::VariableFeatures(cleaned), c("MS4A1", "CD79A"))

  scaled <- SeuratObject::GetAssayData(obj, layer = "scale.data")
  gene <- rownames(scaled)[1]
  boosted <- quietly(OHMmy:::.boost_scaled_genes(obj, "RNA", gene, 3))
  expect_equal(SeuratObject::GetAssayData(boosted, layer = "scale.data")[gene, ],
               scaled[gene, ] * 3)
  # multiplier of 1 leaves the object unchanged
  expect_identical(OHMmy:::.boost_scaled_genes(obj, "RNA", gene, 1), obj)
})

test_that(".suggest_n_pcs returns a value within the computed PCs", {
  obj <- SeuratObject::pbmc_small
  n <- OHMmy:::.suggest_n_pcs(obj, "pca", max_pca = 19)
  expect_true(n >= 1 && n <= 19)
})
