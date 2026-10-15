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

# Reusing existing normalisation / scaling (force_norm, force_scale) -------------------

test_that(".has_layer recognises joined and split layers and tells data from scale.data", {
  obj <- small_seurat()
  expect_true(OHMmy:::.has_layer(obj, "RNA", "data"))
  expect_true(OHMmy:::.has_layer(obj, "RNA", "scale.data"))
  expect_false(OHMmy:::.has_layer(obj, "SCT", "data"))           # assay does not exist

  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch)               # data -> data.B1, data.B2
  expect_true(OHMmy:::.has_layer(obj, "RNA", "data"))

  no_data <- SeuratObject::CreateSeuratObject(SeuratObject::LayerData(small_seurat(), layer = "counts"))
  expect_false(OHMmy:::.has_layer(no_data, "RNA", "data"))
  expect_false(OHMmy:::.has_layer(no_data, "RNA", "scale.data"))  # "data" must not match "scale.data"
})

# Run ProcessSeuratLOG/SCT with integration mocked out, capturing messages
run_without_integration <- function(fun, obj, ...) {
  msgs <- character(0)
  res <- withCallingHandlers(
    quietly_output(fun(obj, batch_col = "batch", reduction_name = "pca.x",
                       integration_reduction = "int.x", dims = 1:10, verbose = FALSE, ...)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  list(obj = res, messages = paste(msgs, collapse = "\n"))
}

test_that("ProcessSeuratLOG reuses or recomputes normalised data as requested", {
  local_mocked_bindings(.integrate_layers = function(seurat_obj, ...) seurat_obj)
  obj <- small_seurat()
  marked <- SeuratObject::LayerData(obj, layer = "data") * 2      # not what NormalizeData produces
  SeuratObject::LayerData(obj, layer = "data") <- marked

  reused <- run_without_integration(ProcessSeuratLOG, obj, force_norm = FALSE)
  expect_match(reused$messages, "Using the existing normalized data layer")
  expect_equal(as.matrix(SeuratObject::LayerData(reused$obj, layer = "data"))[rownames(obj), colnames(obj)],
               as.matrix(marked))

  forced <- run_without_integration(ProcessSeuratLOG, obj, force_norm = TRUE)
  expect_match(forced$messages, "Running LogNormalization")
  expect_false(isTRUE(all.equal(as.matrix(SeuratObject::LayerData(forced$obj, layer = "data"))[rownames(obj), colnames(obj)],
                                as.matrix(marked))))
})

test_that("ProcessSeuratLOG normalises anyway when force_norm = FALSE but no data layer exists", {
  local_mocked_bindings(.integrate_layers = function(seurat_obj, ...) seurat_obj)
  obj <- small_seurat()
  counts_only <- SeuratObject::CreateSeuratObject(SeuratObject::LayerData(obj, layer = "counts"))
  counts_only$batch <- obj$batch

  res <- run_without_integration(ProcessSeuratLOG, counts_only, force_norm = FALSE)
  expect_match(res$messages, "No existing data layer found")
  expect_match(res$messages, "Running LogNormalization")
  expect_true(OHMmy:::.has_layer(res$obj, "RNA", "data"))
})

test_that("ProcessSeuratLOG reuses scale.data only when it covers all variable features", {
  local_mocked_bindings(.integrate_layers = function(seurat_obj, ...) seurat_obj)
  obj <- small_seurat()     # scale.data holds only 20 genes; all 230 genes become variable

  partial <- run_without_integration(ProcessSeuratLOG, obj, force_scale = FALSE)
  expect_match(partial$messages, "does not contain all variable features")
  expect_match(partial$messages, "Scaling data and regressing variables")

  full <- quietly(Seurat::ScaleData(obj, features = rownames(obj), verbose = FALSE))
  reused <- run_without_integration(ProcessSeuratLOG, full, force_scale = FALSE)
  expect_match(reused$messages, "Using the existing scale.data layer")
  expect_false(grepl("Scaling data and regressing variables", reused$messages))

  forced <- run_without_integration(ProcessSeuratLOG, full, force_scale = TRUE)
  expect_match(forced$messages, "Scaling data and regressing variables")
})

test_that("ProcessSeuratSCT reuses an existing SCT assay when force_norm = FALSE", {
  local_mocked_bindings(.integrate_layers = function(seurat_obj, ...) seurat_obj)
  first <- run_without_integration(ProcessSeuratSCT, small_seurat(), vars_to_regress = NULL)
  expect_match(first$messages, "Running SCTransform")
  sct_obj <- first$obj

  # force_norm = FALSE: SCTransform must not run again
  local_mocked_bindings(SCTransform = function(...) stop("SCTransform should not be called"))
  recompute <- run_without_integration(ProcessSeuratSCT, sct_obj, force_norm = FALSE, force_scale = TRUE)
  expect_match(recompute$messages, "Using the existing SCT assay")
  expect_match(recompute$messages, "Recomputing SCT residuals")
  expect_equal(SeuratObject::DefaultAssay(recompute$obj), "SCT")

  # With one model per batch, scale.data lacks variable features that are missing from
  # one batch: force_scale = FALSE tries to add them and reports those it cannot add
  partial <- run_without_integration(ProcessSeuratSCT, sct_obj, force_norm = FALSE, force_scale = FALSE)
  expect_match(partial$messages, "Computing SCT residuals for [0-9]+ variable features")
  expect_match(partial$messages, "will not be used for PCA")

  # When every variable feature already has residuals, they are reused as they are
  covered <- sct_obj
  SeuratObject::DefaultAssay(covered) <- "SCT"
  SeuratObject::VariableFeatures(covered) <-
    rownames(SeuratObject::GetAssayData(covered, assay = "SCT", layer = "scale.data"))
  keep <- run_without_integration(ProcessSeuratSCT, covered, force_norm = FALSE, force_scale = FALSE)
  expect_match(keep$messages, "Using the existing SCT scale.data")
  expect_false(grepl("will not be used for PCA", keep$messages))
})

test_that("ProcessSeuratSCT runs SCTransform when force_norm = FALSE but no SCT assay exists", {
  local_mocked_bindings(.integrate_layers = function(seurat_obj, ...) seurat_obj)

  res <- run_without_integration(ProcessSeuratSCT, small_seurat(), vars_to_regress = NULL, force_norm = FALSE)

  expect_match(res$messages, "No existing SCT assay")
  expect_match(res$messages, "Running SCTransform")
  expect_true("SCT" %in% SeuratObject::Assays(res$obj))
})

# analyze_and_save_clusters() ------------------------------------------------------

# pbmc_small with two resolution sweeps in the OHMmy naming. Clusters are assigned
# in turn, so k clusters of 80 cells have a smallest size of floor(80 / k).
#   log.fastmnn: res 0.5 -> 4, 1 -> 6, 1.5 -> 8 clusters
#   sct.fastmnn: res 0.4 -> 4, 0.8 -> 6, 1.2 -> 6, 2 -> 10 clusters
# Shared cluster numbers: 4 and 6.
sweep_object <- function(drop_rna = TRUE) {
  obj <- SeuratObject::pbmc_small
  if (drop_rna) obj@meta.data <- obj@meta.data[, !grepl("_res\\.", colnames(obj@meta.data))]
  add <- function(col, k) obj@meta.data[[col]] <<- factor(rep(seq_len(k), length.out = ncol(obj)))
  add("log.fastmnn.log_snn_res.0.5", 4); add("log.fastmnn.log_snn_res.1", 6)
  add("log.fastmnn.log_snn_res.1.5", 8)
  add("sct.fastmnn.sct_snn_res.0.4", 4); add("sct.fastmnn.sct_snn_res.0.8", 6)
  add("sct.fastmnn.sct_snn_res.1.2", 6); add("sct.fastmnn.sct_snn_res.2", 10)
  obj
}

test_that("analyze_and_save_clusters compares methods at a shared number of clusters", {
  out <- test_out_dir("cluster_choice")
  res <- quietly(analyze_and_save_clusters(sweep_object(), output_dir = out, target_cluster_num = 6))

  expect_equal(res$Pipeline, c("log.fastmnn", "sct.fastmnn"))
  expect_equal(res$Resolution, c(1, 1.2))          # highest resolution giving 6 clusters
  expect_equal(res$Num_Clusters, c(6L, 6L))
  expect_equal(res$Min_Cluster_Size, c(13L, 13L))   # floor(80 / 6)
  expect_s3_class(attr(res, "plot"), "ggplot")
  expect_setequal(attr(res, "plot_data")$Num_Clusters, c(4, 6))
  expect_true(file.exists(attr(res, "file")))
  expect_length(written_files(out, "^min_cluster_size_plot_.*png$"), 1)
})

test_that("analyze_and_save_clusters defaults to the highest shared number when not interactive", {
  expect_message(
    res <- quietly_output(analyze_and_save_clusters(sweep_object(), output_dir = test_out_dir("cluster_choice2"))),
    "highest common cluster number, 6"
  )
  expect_equal(unique(res$Num_Clusters), 6L)
})

test_that("analyze_and_save_clusters names methods for any graph and can filter them", {
  obj <- sweep_object(drop_rna = FALSE)            # also has pbmc_small's RNA_snn_res.* columns
  plot_methods <- function(r) unique(attr(r, "plot_data")$Pipeline)

  only_rna <- quietly(analyze_and_save_clusters(obj, output_dir = test_out_dir("cluster_choice3"),
                                                pipelines = "RNA_snn"))
  expect_equal(plot_methods(only_rna), "RNA_snn")  # one method, not one per resolution

  two <- quietly(analyze_and_save_clusters(obj, output_dir = test_out_dir("cluster_choice4"),
                                           pipelines = c("log.fastmnn", "sct.fastmnn"), target_cluster_num = 4))
  expect_setequal(two$Pipeline, c("log.fastmnn", "sct.fastmnn"))
  expect_error(analyze_and_save_clusters(obj, pipelines = "nope"), "Unknown pipeline")
})

test_that("analyze_and_save_clusters ignores empty factor levels and non-resolution columns", {
  obj <- sweep_object()
  col <- "log.fastmnn.log_snn_res.1"
  obj@meta.data[[col]] <- factor(obj@meta.data[[col]], levels = 1:7)   # level 7 has no cells
  obj$log.fastmnn.log_snn_res.1_celltype <- "T"                       # not a resolution column

  res <- quietly(analyze_and_save_clusters(obj, output_dir = test_out_dir("cluster_choice5"),
                                           target_cluster_num = 6))
  expect_equal(res$Min_Cluster_Size, c(13L, 13L))   # not 0 from the empty level
  expect_false(any(grepl("celltype", attr(res, "plot_data")$Full_Column)))
})

test_that("analyze_and_save_clusters gives clear errors", {
  obj <- sweep_object()
  out <- test_out_dir("cluster_choice6")
  expect_error(quietly_output(analyze_and_save_clusters(obj, output_dir = out, target_cluster_num = 5)),
               "not among the common")
  expect_error(analyze_and_save_clusters(SeuratObject::pbmc_small@meta.data), "must be a Seurat object")

  none <- obj
  none@meta.data <- none@meta.data[, !grepl("_res\\.", colnames(none@meta.data))]
  expect_error(analyze_and_save_clusters(none, output_dir = out), "No cluster resolution columns")

  disjoint <- obj
  disjoint$sct.fastmnn.sct_snn_res.0.4 <- factor(rep(1:3, length.out = ncol(obj)))
  disjoint$sct.fastmnn.sct_snn_res.0.8 <- factor(rep(1:5, length.out = ncol(obj)))
  disjoint$sct.fastmnn.sct_snn_res.1.2 <- factor(rep(1:5, length.out = ncol(obj)))
  expect_error(quietly_output(analyze_and_save_clusters(disjoint, output_dir = out)),
               "no number of clusters that every pipeline produced")
})
