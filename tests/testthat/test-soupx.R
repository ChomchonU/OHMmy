# Ambient RNA decontamination (SoupX wrappers) ----------------------------------

test_that("get_sample_names parses Cell Ranger folder names", {
  fake_dir <- test_out_dir("cellranger_names")
  dir.create(file.path(fake_dir, "Patient1_filtered_feature_bc_matrix"))
  dir.create(file.path(fake_dir, "Patient2_filtered_feature_bc_matrix"))
  dir.create(file.path(fake_dir, "Control_feature_bc_matrix"))

  ids <- get_sample_names(fake_dir)

  expect_type(ids, "character")
  expect_null(names(ids))
  expect_setequal(ids, c("Patient1", "Patient2", "Control"))
})

test_that("load_counts reads filtered and raw matrices for every sample", {
  cr <- synthetic_cellranger(test_out_dir("cellranger_load"))
  samples <- get_sample_names(cr$filtered_dir)

  counts <- quietly(load_counts(samples, cr$filtered_dir, cr$raw_dir))

  expect_named(counts, c("filtered", "raw"))
  expect_named(counts$filtered, samples)
  expect_equal(unname(sapply(counts$filtered, ncol)), c(80, 80))
  expect_equal(unname(sapply(counts$raw, ncol)), c(380, 380))
  expect_true(all(colnames(counts$filtered$S1) %in% colnames(counts$raw$S1)))
})

test_that("create_soup_channels builds one SoupChannel per sample", {
  cr <- synthetic_cellranger(test_out_dir("cellranger_channels"))
  samples <- get_sample_names(cr$filtered_dir)
  counts <- quietly(load_counts(samples, cr$filtered_dir, cr$raw_dir))

  channels <- quietly(create_soup_channels(counts$raw, counts$filtered, samples))

  expect_named(channels, samples)
  expect_true(all(vapply(channels, inherits, logical(1), what = "SoupChannel")))
})

test_that("prepare_soupx_inputs bundles counts, channels and clustered objects", {
  cr <- synthetic_cellranger(test_out_dir("cellranger_prepare"))

  prep <- quietly(prepare_soupx_inputs(cr$filtered_dir, cr$raw_dir))

  expect_named(prep, c("sample_names", "counts_lists", "soup_list",
                       "seurat_for_clustering", "manual_contam"))
  expect_setequal(prep$sample_names, c("S1", "S2"))
  expect_s4_class(prep$seurat_for_clustering$S1, "Seurat")
  expect_true("seurat_clusters" %in% colnames(prep$seurat_for_clustering$S1@meta.data))
  expect_length(prep$manual_contam, 0)
})

# estimate_contamination() with SoupX's estimation mocked to a fixed rho = 0.05,
# so that only OHMmy's adjustment logic is tested.
mock_contamination_inputs <- function() {
  samples <- c("S1", "S2")
  toc <- lapply(1:2, function(i) synthetic_counts(seed = i))
  names(toc) <- samples
  tod <- lapply(toc, function(m) {
    set.seed(10)
    empty <- Matrix::Matrix(matrix(stats::rpois(nrow(m) * 300, 0.05), nrow(m), 300,
                                   dimnames = list(rownames(m), paste0("EMPTY", 1:300, "-1"))),
                            sparse = TRUE)
    cbind(m, empty)
  })
  soup <- lapply(samples, function(s) quietly(SoupX::SoupChannel(tod[[s]], toc[[s]])))
  names(soup) <- samples
  seur <- lapply(toc, function(m) {
    o <- suppressWarnings(SeuratObject::CreateSeuratObject(m))
    o$seurat_clusters <- factor(rep(0:1, length.out = ncol(o)))
    o
  })
  list(samples = samples, toc = toc, soup = soup, seur = seur)
}

test_that("estimate_contamination applies mulFac globally or selectively", {
  local_mocked_bindings(
    setClusters  = function(sc, clusters) sc,
    autoEstCont  = function(sc, ...) { sc$metaData$rho <- 0.05; sc },
    adjustCounts = function(sc, ...) sc$toc
  )
  inp <- mock_contamination_inputs()
  rho <- function(res) sapply(res$soup_list_all, function(sc) unique(sc$metaData$rho))

  auto <- quietly(estimate_contamination(inp$soup, inp$seur, inp$toc, inp$samples))
  global <- quietly(estimate_contamination(inp$soup, inp$seur, inp$toc, inp$samples, mulFac = 5))
  selective <- quietly(estimate_contamination(inp$soup, inp$seur, inp$toc, inp$samples,
                                              mulFac = 5, manual = TRUE, manual_contam = "S2"))

  expect_named(auto, c("soup_list_all", "corrected_list"))
  expect_equal(unname(rho(auto)), c(0.05, 0.05))
  expect_equal(unname(rho(global)), c(0.10, 0.10))
  expect_equal(unname(rho(selective)), c(0.05, 0.10))
})

test_that("run_soupx_post_clustering restricts multiFac to manual_contam samples", {
  local_mocked_bindings(
    setClusters  = function(sc, clusters) sc,
    autoEstCont  = function(sc, ...) { sc$metaData$rho <- 0.05; sc },
    adjustCounts = function(sc, ...) sc$toc,
    create_final_seurat = function(soup_list_all, sample_names) {
      stats::setNames(as.list(sample_names), sample_names)
    }
  )
  inp <- mock_contamination_inputs()
  prep <- list(sample_names = inp$samples, soup_list = inp$soup,
               seurat_for_clustering = inp$seur, counts_lists = list(filtered = inp$toc))
  rho <- function(res) unname(sapply(res$soup_objects, function(sc) unique(sc$metaData$rho)))

  all_samples <- quietly(run_soupx_post_clustering(prep, multiFac = 10))
  only_s1 <- quietly(run_soupx_post_clustering(prep, multiFac = 10, manual_contam = "S1"))

  expect_named(all_samples, c("final_seurat", "soup_objects", "seurat_clustering"))
  expect_equal(rho(all_samples), c(0.15, 0.15))
  expect_equal(rho(only_s1), c(0.15, 0.05))
})

test_that("process_soupx_samples chains preparation and correction", {
  calls <- list()
  local_mocked_bindings(
    prepare_soupx_inputs = function(filtered_dir, raw_dir, manual_contam = NULL) {
      calls$prepare <<- c(filtered_dir, raw_dir)
      list(tag = "prep")
    },
    run_soupx_post_clustering = function(prep_output, multiFac = 0, ...) {
      calls$run <<- list(prep_output, multiFac)
      list(final_seurat = "ok", soup_objects = NULL, seurat_clustering = NULL)
    }
  )

  res <- process_soupx_samples("filt", "raw", multiFac = 3)

  expect_equal(calls$prepare, c("filt", "raw"))
  expect_equal(calls$run, list(list(tag = "prep"), 3))
  expect_equal(res$final_seurat, "ok")
})

test_that("create_final_seurat builds Seurat objects with mitoPercent", {
  inp <- mock_contamination_inputs()
  soup <- lapply(inp$soup, SoupX::setContaminationFraction, contFrac = 0.05)

  res <- quietly(create_final_seurat(soup, inp$samples))

  expect_named(res, inp$samples)
  expect_s4_class(res$S1, "Seurat")
  expect_true("mitoPercent" %in% colnames(res$S1@meta.data))
})

test_that("addSoupXMetaToSeurat transfers corrected counts and QC metadata", {
  toc <- list(S1 = synthetic_counts(seed = 1), S2 = synthetic_counts(seed = 2))
  soupx_list <- lapply(names(toc), function(s) {
    m <- toc[[s]]
    m@x <- pmax(m@x - 1, 0)                      # pretend SoupX removed one UMI per entry
    o <- suppressWarnings(SeuratObject::CreateSeuratObject(m, project = s))
    o$mitoPercent <- Seurat::PercentageFeatureSet(o, pattern = "^MT-")
    o
  })
  names(soupx_list) <- names(toc)

  original_counts <- do.call(cbind, lapply(names(toc), function(s) {
    m <- toc[[s]]; colnames(m) <- paste0(s, "_", sub("-1$", "", colnames(m))); m
  }))
  original <- suppressWarnings(SeuratObject::CreateSeuratObject(original_counts))

  res <- quietly(addSoupXMetaToSeurat(original, soupx_list))

  expect_s4_class(res, "Seurat")
  expect_equal(ncol(res), ncol(original))
  expect_true(all(c("nCount_RNA", "nFeature_RNA", "pct_counts_mt") %in% colnames(res@meta.data)))
  expect_lt(sum(SeuratObject::LayerData(res, layer = "counts")), sum(original_counts))
})
