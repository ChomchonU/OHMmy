# =============================================================================
# Shared helpers for the OHMmy vignettes
# -----------------------------------------------------------------------------
# Every OHMmy article works on the same sample object: the public 10x Genomics
# "3k PBMCs from a healthy donor" dataset (pbmc3k), processed with OHMmy's own
# wrappers. To make the cohort-level functions (abundance testing, confounder
# checks, pseudo-bulk DE) demonstrable, the cells are distributed across a
# *simulated* cohort of 12 donors. The donor, severity, batch, age and sex
# columns are therefore synthetic and carry no biological meaning.
#
# The processed object is cached so it is only built once per machine. Set the
# environment variable OHMMY_VIGNETTE_CACHE to control where it is stored.
# =============================================================================

ohmmy_cache_dir <- function() {
  d <- Sys.getenv("OHMMY_VIGNETTE_CACHE", unset = tools::R_user_dir("OHMmy", "cache"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Locate (or download once) the pbmc3k filtered count matrix --------------------
pbmc3k_dir <- function() {
  candidates <- c(
    "filtered_gene_bc_matrices/hg19",
    file.path(ohmmy_cache_dir(), "filtered_gene_bc_matrices", "hg19")
  )
  hit <- candidates[file.exists(file.path(candidates, "matrix.mtx"))]
  if (length(hit) > 0) return(hit[1])

  tgz <- file.path(ohmmy_cache_dir(), "pbmc3k_filtered_gene_bc_matrices.tar.gz")
  utils::download.file(
    "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
    destfile = tgz, mode = "wb", quiet = TRUE
  )
  utils::untar(tgz, exdir = ohmmy_cache_dir())
  candidates[2]
}

pbmc3k_counts <- function() Seurat::Read10X(data.dir = pbmc3k_dir())

# Simulated cohort --------------------------------------------------------------
# 12 donors: 4 Healthy, 4 Mild, 4 Severe, split over two sequencing batches.
# Monocyte-like cells (high CD14/LYZ/S100A8/S100A9/FCGR3A/MS4A7) are assigned
# preferentially to Mild/Severe donors so that the abundance tests have a
# deliberate, known effect to recover.
mock_donor_table <- function() {
  data.frame(
    donor_id = sprintf("D%02d", 1:12),
    severity = factor(rep(c("Healthy", "Mild", "Severe"), each = 4),
                      levels = c("Healthy", "Mild", "Severe")),
    batch    = rep(c("Batch_1", "Batch_2"), times = 6),
    age      = c(34, 51, 28, 45, 39, 62, 47, 30, 55, 41, 66, 36),
    sex      = c("F", "M", "M", "F", "F", "M", "F", "M", "M", "F", "M", "F"),
    stringsAsFactors = FALSE
  )
}

add_mock_cohort <- function(seurat_obj, seed = 42) {
  set.seed(seed)
  donors <- mock_donor_table()

  counts   <- SeuratObject::LayerData(seurat_obj, assay = "RNA", layer = "counts")
  mono_g   <- intersect(c("CD14", "LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7"), rownames(counts))
  mono_scr <- log1p(Matrix::colSums(counts[mono_g, , drop = FALSE]) / Matrix::colSums(counts) * 1e4)
  is_mono  <- mono_scr > 4.5

  w_mono  <- c(Healthy = 1, Mild = 1.8, Severe = 3)[as.character(donors$severity)]
  assigned <- character(ncol(seurat_obj))
  assigned[is_mono]  <- sample(donors$donor_id, sum(is_mono),  replace = TRUE, prob = w_mono)
  assigned[!is_mono] <- sample(donors$donor_id, sum(!is_mono), replace = TRUE)

  md <- donors[match(assigned, donors$donor_id), ]
  rownames(md) <- colnames(seurat_obj)
  md$infection_status <- ifelse(md$severity == "Healthy", "Healthy", "Infected")
  SeuratObject::AddMetaData(seurat_obj, md)
}

# Marker-based annotation (canonical PBMC markers) ------------------------------
pbmc_marker_sets <- list(
  "CD4 T"       = c("IL7R", "CCR7", "LDHB", "CD3E", "MAL"),
  "CD8 T"       = c("CD8A", "CD8B", "CD3D", "GZMK"),
  "NK"          = c("GNLY", "NKG7", "GZMB", "PRF1", "FCGR3A", "KLRD1"),
  "B"           = c("MS4A1", "CD79A", "CD79B", "CD19"),
  "CD14 Mono"   = c("CD14", "LYZ", "S100A8", "S100A9"),
  "FCGR3A Mono" = c("FCGR3A", "MS4A7", "LST1", "IFITM2"),
  "DC"          = c("FCER1A", "CST3", "CLEC10A"),
  "Platelet"    = c("PPBP", "PF4", "GNG11")
)

annotate_pbmc <- function(seurat_obj, cluster_col) {
  sets <- lapply(pbmc_marker_sets, intersect, rownames(seurat_obj))
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = sets, name = "annot_", seed = 1)
  score_cols <- paste0("annot_", seq_along(sets))
  per_cluster <- aggregate(seurat_obj@meta.data[, score_cols],
                           by = list(cluster = seurat_obj@meta.data[[cluster_col]]), FUN = mean)
  winner <- names(sets)[apply(per_cluster[, score_cols], 1, which.max)]
  lookup <- setNames(winner, per_cluster$cluster)
  labels <- unname(lookup[as.character(seurat_obj@meta.data[[cluster_col]])])
  seurat_obj$cell_type <- droplevels(factor(labels, levels = names(pbmc_marker_sets)))
  seurat_obj@meta.data[score_cols] <- NULL
  seurat_obj
}

# Build the processed sample object with OHMmy -------------------------------------
# This mirrors, step by step, the code shown in the
# "Data Processing, Integration, and Clustering" vignette.
build_pbmc3k_ohmmy <- function(plot_dir = tempdir()) {
  pbmc <- Seurat::CreateSeuratObject(pbmc3k_counts(), project = "pbmc3k",
                                     min.cells = 3, min.features = 200)
  pbmc[["pct_counts_mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")
  pbmc[["percent.ribo"]]  <- Seurat::PercentageFeatureSet(pbmc, pattern = "^RP[SL]")
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & pct_counts_mt < 5)
  pbmc <- add_mock_cohort(pbmc)

  pbmc <- OHMmy::ProcessSeuratLOG(
    pbmc,
    batch_col             = "batch",
    vars_to_regress       = "pct_counts_mt",
    tcr_bcr_patterns      = "^TR[ABDG]V|^IG[HKL]V",
    reduction_name        = "pca.log",
    integration_method    = "HarmonyIntegration",
    integration_reduction = "harmony.log",
    dims                  = 1:30,
    verbose               = FALSE
  )

  pbmc <- Seurat::CellCycleScoring(pbmc,
                                   s.features   = Seurat::cc.genes.updated.2019$s.genes,
                                   g2m.features = Seurat::cc.genes.updated.2019$g2m.genes)

  clus <- OHMmy::ClusterAndUMAP(
    pbmc,
    sample_name         = "pbmc3k",
    reduction           = "harmony.log",
    dims                = 1:20,
    umap_name           = "umap.har",
    graph_name          = "har_snn",
    cluster_prefix      = "har_snn_res.",
    cluster_resolutions = seq(0.2, 1.2, by = 0.2),
    final_resolution    = 0.6,
    plot_dir            = plot_dir,
    dpi                 = 100
  )
  pbmc <- clus$seurat
  pbmc$seurat_clusters <- pbmc$har_snn_res.0.6
  annotate_pbmc(pbmc, cluster_col = "seurat_clusters")
}

load_pbmc3k_ohmmy <- function(rebuild = FALSE) {
  rds <- file.path(ohmmy_cache_dir(), sprintf("pbmc3k_ohmmy_v%s.rds", utils::packageVersion("OHMmy")))
  if (!rebuild && file.exists(rds)) return(readRDS(rds))
  # Building prints Seurat progress output; keep the articles clean
  utils::capture.output(obj <- suppressWarnings(suppressMessages(build_pbmc3k_ohmmy())))
  saveRDS(obj, rds)
  obj
}

# Simulated Cell Ranger output for the SoupX workflow ------------------------------
# pbmc3k only ships the *filtered* matrix, but SoupX also needs the *raw* matrix
# (all droplets) to learn the ambient-RNA profile. This helper splits pbmc3k into
# two pseudo-samples and, for each, writes:
#   <sample>_filtered_feature_bc_matrix/  the real cells
#   <sample>_raw_feature_bc_matrix/       the real cells + simulated empty droplets
# Empty droplets draw a small number of UMIs (median ~20) from the pooled
# expression profile of the sample, which is how ambient RNA behaves in practice.
write_mtx_dir <- function(mat, path, features) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  mtx <- file.path(path, "matrix.mtx")
  Matrix::writeMM(methods::as(mat, "CsparseMatrix"), mtx)
  R.utils::gzip(mtx, overwrite = TRUE)
  gz_lines <- function(x, f) { con <- gzfile(f, "w"); writeLines(x, con); close(con) }
  gz_lines(paste(features$id, features$symbol, "Gene Expression", sep = "\t"),
           file.path(path, "features.tsv.gz"))
  gz_lines(colnames(mat), file.path(path, "barcodes.tsv.gz"))
}

make_mock_cellranger <- function(root, n_empty = 6000, seed = 7) {
  set.seed(seed)
  filt_root <- file.path(root, "filtered")
  raw_root  <- file.path(root, "raw")

  genes <- utils::read.delim(file.path(pbmc3k_dir(), "genes.tsv"), header = FALSE,
                             col.names = c("id", "symbol"))
  counts <- Seurat::Read10X(pbmc3k_dir(), gene.column = 1)   # Ensembl IDs keep rows unique
  sample_of_cell <- sample(c("PBMC_A", "PBMC_B"), ncol(counts), replace = TRUE)

  for (s in c("PBMC_A", "PBMC_B")) {
    toc <- counts[, sample_of_cell == s]
    soup_prob <- Matrix::rowSums(toc) / sum(toc)
    lib <- pmax(1, round(stats::rlnorm(n_empty, meanlog = log(20), sdlog = 0.6)))
    empty <- Matrix::sparseMatrix(
      i    = sample.int(length(soup_prob), sum(lib), replace = TRUE, prob = soup_prob),
      j    = rep(seq_len(n_empty), lib),
      x    = 1,
      dims = c(length(soup_prob), n_empty)
    )
    empty_bc <- paste0(replicate(n_empty, paste(sample(c("A", "C", "G", "T"), 14, TRUE),
                                                collapse = "")), "-1")
    colnames(empty) <- make.unique(empty_bc)
    rownames(empty) <- rownames(toc)
    tod <- cbind(toc, empty)

    write_mtx_dir(toc, file.path(filt_root, paste0(s, "_filtered_feature_bc_matrix")), genes)
    write_mtx_dir(tod, file.path(raw_root,  paste0(s, "_raw_feature_bc_matrix")),      genes)
  }
  list(filtered_dir = filt_root, raw_dir = raw_root)
}

# Display helpers ---------------------------------------------------------------
# Many OHMmy functions write timestamped figures to disk. `show_saved()` finds the
# newest file matching `pattern` in `dir`, copies it next to the knitted figures
# (downsampling very large images when the 'jpeg'/'png' packages are installed)
# and embeds it.
.downsample_image <- function(src, dest, max_px = 2400) {
  ext <- tolower(tools::file_ext(src))
  pkg <- switch(ext, jpg = , jpeg = "jpeg", png = "png", NA_character_)
  if (is.na(pkg) || !requireNamespace(pkg, quietly = TRUE)) {
    return(file.copy(src, dest, overwrite = TRUE))  # optional: install jpeg/png to downsample
  }
  reader <- if (pkg == "jpeg") jpeg::readJPEG else png::readPNG
  img <- reader(src)
  f <- ceiling(max(dim(img)[1:2]) / max_px)
  if (f > 1) {
    h <- (dim(img)[1] %/% f) * f
    w <- (dim(img)[2] %/% f) * f
    img <- img[seq_len(h), seq_len(w), , drop = FALSE]
    small <- array(0, c(h / f, w / f, dim(img)[3]))
    for (i in seq_len(f)) for (j in seq_len(f)) {
      small <- small + img[seq(i, h, by = f), seq(j, w, by = f), , drop = FALSE]
    }
    img <- small / f^2
  }
  if (ext == "png") png::writePNG(img, dest) else jpeg::writeJPEG(img, dest, quality = 0.9)
}

saved_files <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  f[order(file.mtime(f), decreasing = TRUE)]
}

show_saved <- function(dir, pattern, n = 1, max_px = 2400) {
  files <- utils::head(saved_files(dir, pattern), n)
  if (length(files) == 0) stop("No saved figure matches '", pattern, "' in ", dir)
  fig_dir <- knitr::opts_chunk$get("fig.path")
  if (is.null(fig_dir) || !nzchar(fig_dir)) fig_dir <- "figure/"
  dir.create(dirname(paste0(fig_dir, "x")), recursive = TRUE, showWarnings = FALSE)
  dest <- paste0(fig_dir, "saved-", basename(files))
  for (i in seq_along(files)) .downsample_image(files[i], dest[i], max_px = max_px)
  knitr::include_graphics(dest)
}
