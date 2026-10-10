# Shared fixtures for the OHMmy test suite -------------------------------------
# Files named helper-*.R are sourced by testthat before the tests run.

# A small Seurat v5 object based on SeuratObject::pbmc_small (80 cells, 230 genes)
# with the metadata columns the OHMmy functions expect.
small_seurat <- function() {
  obj <- suppressWarnings(suppressMessages(SeuratObject::UpdateSeuratObject(SeuratObject::pbmc_small)))
  suppressWarnings(obj[["RNA"]] <- methods::as(obj[["RNA"]], "Assay5"))

  set.seed(42)
  n <- ncol(obj)
  obj$batch         <- rep(c("B1", "B2"), length.out = n)
  obj$pct_counts_mt <- stats::runif(n, 0, 5)
  obj$cell_type     <- paste0("C", obj$RNA_snn_res.1)
  obj$donor         <- rep(paste0("D", 1:8), length.out = n)
  obj$condition     <- ifelse(obj$donor %in% paste0("D", 1:4), "Ctrl", "Case")
  obj$severity      <- c(D1 = "Healthy", D2 = "Healthy", D3 = "Healthy",
                         D4 = "Mild", D5 = "Mild", D6 = "Mild",
                         D7 = "Severe", D8 = "Severe")[obj$donor] |> unname()
  obj$age           <- c(D1 = 30, D2 = 41, D3 = 52, D4 = 35, D5 = 47,
                         D6 = 58, D7 = 39, D8 = 61)[obj$donor] |> unname()
  obj$sex           <- c(D1 = "F", D2 = "M", D3 = "F", D4 = "M", D5 = "F",
                         D6 = "M", D7 = "F", D8 = "M")[obj$donor] |> unname()
  obj
}

# A fresh, empty output directory per test
test_out_dir <- function(name) {
  d <- file.path(tempdir(), "ohmmy-tests", name)
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

# Files written into a directory (non-recursive by default)
written_files <- function(dir, pattern = NULL, recursive = FALSE) {
  list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
}

# Run code while silencing console output, messages and progress bars
quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(suppressMessages(expr)))
  out
}

# Marker gene panels present in pbmc_small, in the prefix format used by the
# plot_dot_dendro* functions
small_feature_df <- function(prefix = "Test") {
  data.frame(
    group   = paste0(prefix, "_", c("B", "B", "B", "Mono", "Mono", "Mono", "T", "T", "NK", "NK")),
    feature = c("MS4A1", "CD79A", "CD79B", "LYZ", "S100A9", "CST3", "CD3E", "IL7R", "GNLY", "NKG7")
  )
}

# Write a 10x-style matrix directory (gzipped files with 3-column features.tsv)
write_10x_dir <- function(mat, path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  mtx <- file.path(path, "matrix.mtx")
  Matrix::writeMM(methods::as(mat, "CsparseMatrix"), mtx)
  gz_copy <- function(src, dest) {
    con <- gzfile(dest, "wb"); writeBin(readBin(src, "raw", file.size(src)), con); close(con)
  }
  gz_copy(mtx, paste0(mtx, ".gz")); unlink(mtx)
  con <- gzfile(file.path(path, "features.tsv.gz"), "w")
  writeLines(paste(rownames(mat), rownames(mat), "Gene Expression", sep = "\t"), con); close(con)
  con <- gzfile(file.path(path, "barcodes.tsv.gz"), "w")
  writeLines(colnames(mat), con); close(con)
}

# Synthetic dense count matrix: two cell groups, a few MT- genes, enough features
# per cell (> 200) to survive CreateSeuratObject(min.features = 200)
synthetic_counts <- function(n_cells = 80, n_genes = 300, seed = 1, prefix = "CELL") {
  set.seed(seed)
  genes <- c(paste0("MT-G", 1:5), paste0("GENE", seq_len(n_genes - 5)))
  group <- rep(1:2, length.out = n_cells)
  mu <- matrix(2, n_genes, n_cells)
  mu[6:40, group == 1] <- 8
  mu[41:75, group == 2] <- 8
  m <- matrix(stats::rpois(n_genes * n_cells, mu), n_genes, n_cells,
              dimnames = list(genes, paste0(prefix, seq_len(n_cells), "-1")))
  Matrix::Matrix(m, sparse = TRUE)
}

# Two-sample Cell Ranger-style layout (filtered + raw with empty droplets)
synthetic_cellranger <- function(root, n_empty = 300) {
  for (s in c("S1", "S2")) {
    toc <- synthetic_counts(seed = match(s, c("S1", "S2")))
    set.seed(10)
    empty <- Matrix::Matrix(matrix(stats::rpois(nrow(toc) * n_empty, 0.05), nrow(toc), n_empty,
                                   dimnames = list(rownames(toc), paste0("EMPTY", seq_len(n_empty), "-1"))),
                            sparse = TRUE)
    write_10x_dir(toc, file.path(root, "filtered", paste0(s, "_filtered_feature_bc_matrix")))
    write_10x_dir(cbind(toc, empty), file.path(root, "raw", paste0(s, "_raw_feature_bc_matrix")))
  }
  list(filtered_dir = file.path(root, "filtered"), raw_dir = file.path(root, "raw"))
}

# Silence messages only (keeps warnings visible to expect_warning)
quietly_messages <- function(expr) suppressMessages(expr)
