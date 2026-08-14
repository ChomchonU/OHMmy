#' Title
#'
#' @param filtered_dir
#'
#' @returns
#' @export
#'
#' @examples
get_sample_names <- function(filtered_dir) {
  dirs <- list.dirs(path = filtered_dir, recursive = FALSE, full.names = FALSE)
  sample_names <- sapply(dirs, function(x) gsub("(_filtered)?_feature_bc_matrix", "", x))
  return(sample_names)
}

#--------------------------------------------------

#' Title
#'
#' @param sample_names
#' @param filtered_dir
#' @param raw_dir
#'
#' @returns
#' @export
#'
#' @examples
load_counts <- function(sample_names, filtered_dir, raw_dir) {
  cts_filtered_list <- list()
  cts_raw_list <- list()

  cat("Loading counts:\n")
  pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    filtered_path <- file.path(filtered_dir, paste0(sample_name, "_filtered_feature_bc_matrix"))
    raw_path <- file.path(raw_dir, paste0(sample_name, "_raw_feature_bc_matrix"))

    cts_filtered <- ReadMtx(
      mtx = paste0(filtered_path, '/matrix.mtx.gz'),
      features = paste0(filtered_path, "/features.tsv.gz"),
      cells = paste0(filtered_path, "/barcodes.tsv.gz")
    )

    cts_raw <- ReadMtx(
      mtx = paste0(raw_path, '/matrix.mtx.gz'),
      features = paste0(raw_path, "/features.tsv.gz"),
      cells = paste0(raw_path, "/barcodes.tsv.gz")
    )

    cts_filtered_list[[sample_name]] <- cts_filtered
    cts_raw_list[[sample_name]] <- cts_raw

    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(list(filtered = cts_filtered_list, raw = cts_raw_list))
}

# ---------------------------------------------------------

#' Title
#'
#' @param cts_raw_list
#' @param cts_filtered_list
#' @param sample_names
#'
#' @returns
#' @export
#'
#' @examples
create_soup_channels <- function(cts_raw_list, cts_filtered_list, sample_names) {
  soup_list <- list()
  cat("Creating SoupChannel objects:\n")
  pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    sc <- SoupChannel(tod = cts_raw_list[[sample_name]], toc = cts_filtered_list[[sample_name]])
    soup_list[[sample_name]] <- sc
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(soup_list)
}


#---------------------------------------------------------

#' Title
#'
#' @param cts_filtered_list
#' @param sample_names
#'
#' @returns
#' @export
#'
#' @examples
create_seurat_for_clustering <- function(cts_filtered_list, sample_names) {
  seurat_list <- list()
  cat("Clustering for SoupX (Seurat pipeline):\n")
  pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    seurat_temp <- CreateSeuratObject(counts = cts_filtered_list[[sample_name]], project = sample_name)
    seurat_temp <- NormalizeData(seurat_temp)
    seurat_temp <- FindVariableFeatures(seurat_temp)
    seurat_temp <- ScaleData(seurat_temp)
    # seurat_temp <- SCTransform(seurat_temp)
    seurat_temp <- RunPCA(seurat_temp, verbose = FALSE)
    seurat_temp <- FindNeighbors(seurat_temp, dims = 1:50)
    seurat_temp <- FindClusters(seurat_temp, resolution = 1.8)
    seurat_temp <- RunUMAP(seurat_temp, dims = 1:50)

    seurat_list[[sample_name]] <- seurat_temp
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(seurat_list)
}

# -----------------------------------------------------

#' Title
#'
#' @param soup_list
#' @param seurat_list
#' @param cts_filtered_list
#' @param sample_names
#' @param mulFac
#' @param manual_contam
#' @param manual
#' @param methods
#' @param pCut
#'
#' @returns
#' @export
#'
#' @examples
estimate_contamination <- function(
    soup_list,
    seurat_list,
    cts_filtered_list,
    sample_names,
    mulFac = 0,
    manual_contam = list(),     # character vector of sample names to override
    manual = FALSE,             # if TRUE, apply mulFac to all
    methods = "subtraction",
    pCut = 0.01
) {
  soup_list_all <- list()
  corrected_list <- list()

  cat("Estimating contamination:\n")
  pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    sc <- soup_list[[sample_name]]
    seurat_temp <- seurat_list[[sample_name]]

    # Assign clusters
    clusters <- seurat_temp$seurat_clusters
    names(clusters) <- colnames(cts_filtered_list[[sample_name]])
    sc <- setClusters(sc, clusters)

    # Auto estimate contamination
    sc <- autoEstCont(sc)

    # Extract estimated rho
    estimated_rho <- unique(sc$metaData$rho)
    if (length(estimated_rho) != 1) {
      stop("Sample ", sample_name, ": Expected a single rho value but got: ", paste(estimated_rho, collapse = ", "))
    }
    estimated_rho <- as.numeric(estimated_rho)

    # Determine final contamination fraction
    if (!manual) {
      confrac <- estimated_rho + (mulFac / 100)
      message("Sample: ", sample_name, " | MANUAL mode active | Adjusted rho: ", round(confrac, 4))
    } else {
      if (sample_name %in% manual_contam) {
        confrac <- estimated_rho + (mulFac / 100)
        message("Sample: ", sample_name, " | SELECTIVE manual override | Adjusted rho: ", round(confrac, 4))
      } else {
        confrac <- estimated_rho
        message("Sample: ", sample_name, " | Auto-estimated rho: ", round(confrac, 4))
      }
    }

    # Apply contamination fraction
    sc <- setContaminationFraction(sc, confrac)
    soup_list_all[[sample_name]] <- sc

    # Announce method used
    msg <- paste0("Sample: ", sample_name, " | Adjusting counts using method = '", methods, "'")
    if (methods == "soupOnly") {
      msg <- paste0(msg, " (with pCut = ", pCut, ")")
    }
    message(msg)

    # Adjust counts
    corrected <- adjustCounts(sc, method = methods, pCut = pCut, roundToInt = TRUE)
    corrected_list[[sample_name]] <- corrected

    setTxtProgressBar(pb, i)
  }

  close(pb)

  return(list(
    soup_list_all = soup_list_all,
    corrected_list = corrected_list
  ))
}

# --------------------------------------------------------

#' Title
#'
#' @param soup_list_all
#' @param sample_names
#'
#' @returns
#' @export
#'
#' @examples
create_final_seurat <- function(soup_list_all, sample_names) {
  final_seurat_list <- list()
  cat("Creating final Seurat objects with corrected counts:\n")
  pb <- txtProgressBar(min = 0, max = length(sample_names), style = 3)

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    sc <- soup_list_all[[sample_name]]
    corrected_counts <- adjustCounts(sc, roundToInt = TRUE)

    seurat_ind <- CreateSeuratObject(counts = corrected_counts, project = sample_name, min.cells = 3, min.features = 200)
    seurat_ind$mitoPercent <- PercentageFeatureSet(seurat_ind, pattern = "^MT-")

    final_seurat_list[[sample_name]] <- seurat_ind
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(final_seurat_list)
}


# --------------------------------------------------------

#' Title
#'
#' @param filtered_dir
#' @param raw_dir
#' @param manual_contam
#'
#' @returns
#' @export
#'
#' @examples
prepare_soupx_inputs <- function(filtered_dir, raw_dir, manual_contam = NULL) {
  sample_names <- get_sample_names(filtered_dir)

  counts_lists <- load_counts(sample_names, filtered_dir, raw_dir)

  soup_list <- create_soup_channels(counts_lists$raw, counts_lists$filtered, sample_names)

  seurat_for_clustering <- create_seurat_for_clustering(counts_lists$filtered, sample_names)

  # If manual_contam is NULL, default to empty list
  if (is.null(manual_contam)) {
    manual_contam <- list()
  }

  return(list(
    sample_names = sample_names,
    counts_lists = counts_lists,
    soup_list = soup_list,
    seurat_for_clustering = seurat_for_clustering,
    manual_contam = manual_contam
  ))
}


# ----------------------------------------------------------

#' Title
#'
#' @param prep_output
#' @param multiFac
#' @param manual_contam
#' @param pCut
#' @param methods
#'
#' @returns
#' @export
#'
#' @examples
run_soupx_post_clustering <- function(
    prep_output,
    multiFac = 0,
    manual_contam = list(),
    pCut = 0.01,
    methods = "subtraction"       # optional: force ambient gene correction
) {
  sample_names <- prep_output$sample_names
  soup_list <- prep_output$soup_list
  seurat_for_clustering <- prep_output$seurat_for_clustering
  counts_filtered <- prep_output$counts_lists$filtered

  contamination_result <- estimate_contamination(
    soup_list = soup_list,
    seurat_list = seurat_for_clustering,
    cts_filtered_list = counts_filtered,
    sample_names = sample_names,
    mulFac = multiFac,
    manual_contam = manual_contam,
    methods = methods,
    pCut = pCut
  )

  final_seurat_list <- create_final_seurat(
    soup_list_all = contamination_result$soup_list_all,
    sample_names = sample_names
  )

  return(list(
    final_seurat = final_seurat_list,
    soup_objects = contamination_result$soup_list_all,
    seurat_clustering = seurat_for_clustering
  ))
}

# ---------------------------------------------------------------

#' Title
#'
#' @param filtered_dir
#' @param raw_dir
#' @param multiFac
#'
#' @returns
#' @export
#'
#' @examples
process_soupx_samples <- function(filtered_dir, raw_dir, multiFac=0) {
  sample_names <- get_sample_names(filtered_dir)

  counts_lists <- load_counts(sample_names, filtered_dir, raw_dir)

  soup_list <- create_soup_channels(counts_lists$raw, counts_lists$filtered, sample_names)

  seurat_for_clustering <- create_seurat_for_clustering(counts_lists$filtered, sample_names)

  soup_list_all <- estimate_contamination(soup_list, seurat_for_clustering, counts_lists$filtered, sample_names, multiFac)

  final_seurat_list <- create_final_seurat(soup_list_all, sample_names)

  return(list(final_seurat = final_seurat_list,
              soup_objects = soup_list_all,
              seurat_clustering = seurat_for_clustering))
}

# --------------------------------------------------------------

#' Title
#'
#' @param original_seurat
#' @param soupx_seurat_list
#' @param mito_colname
#' @param counts_assay
#' @param rename_barcodes
#' @param save_path
#'
#' @returns
#' @export
#'
#' @examples
addSoupXMetaToSeurat <- function(original_seurat,
                                 soupx_seurat_list,
                                 mito_colname = "mitoPercent",
                                 counts_assay = "RNA",
                                 rename_barcodes = TRUE,
                                 save_path = NULL) {
  library(Seurat)

  message("Merging SoupX Seurat objects...")
  merged_seurat <- merge(
    soupx_seurat_list[[1]],
    y = soupx_seurat_list[-1],
    add.cell.ids = names(soupx_seurat_list),
    project = "MergedSoupX"
  )

  if (rename_barcodes) {
    message("Renaming barcodes (removing '-1')...")
    old_barcodes <- colnames(merged_seurat)
    new_barcodes <- gsub("-1$", "", old_barcodes)
    colnames(merged_seurat) <- new_barcodes
  }

  # Standardize mito column name
  if (mito_colname %in% colnames(merged_seurat@meta.data)) {
    colnames(merged_seurat@meta.data)[colnames(merged_seurat@meta.data) == mito_colname] <- "pct_counts_mt"
  }

  message("Transferring metadata (nCount_RNA, nFeature_RNA, pct_counts_mt)...")
  required_cols <- c("nCount_RNA", "nFeature_RNA", "pct_counts_mt")
  missing_cols <- setdiff(required_cols, colnames(merged_seurat@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in SoupX metadata: ", paste(missing_cols, collapse = ", "))
  }

  # Match metadata by shared barcodes
  matched_cells <- intersect(colnames(original_seurat), rownames(merged_seurat@meta.data))
  if (length(matched_cells) == 0) {
    stop("No overlapping cell barcodes found between original and merged Seurat objects.")
  }

  matched_meta <- merged_seurat@meta.data[matched_cells, required_cols, drop = FALSE]
  updated_seurat <- AddMetaData(original_seurat, metadata = matched_meta)

  message("Joining RNA assay layers in merged object...")
  merged_seurat[[counts_assay]] <- JoinLayers(merged_seurat[[counts_assay]])

  message("Preparing raw count transfer...")
  cells_to_use <- colnames(original_seurat)
  genes_to_use <- rownames(original_seurat[[counts_assay]])

  missing_cells <- setdiff(cells_to_use, colnames(merged_seurat))
  missing_genes <- setdiff(genes_to_use, rownames(merged_seurat[[counts_assay]]))
  if (length(missing_cells) > 0) {
    message(" Warning: Some cells not found in merged object. Subsetting...")
  }

  if (length(missing_genes) > 0) {
    message(" Warning: Some genes not found in merged object. Subsetting...")
  }

  # Subset to only available ones (preserve order)
  valid_cells <- intersect(cells_to_use, colnames(merged_seurat))
  valid_genes <- intersect(genes_to_use, rownames(merged_seurat[[counts_assay]]))

  raw_counts <- GetAssayData(merged_seurat, assay = counts_assay, layer = "counts")[valid_genes, valid_cells, drop = FALSE]

  # Subset the original object to match valid dimensions (if necessary)
  updated_seurat <- subset(updated_seurat, cells = valid_cells, features = valid_genes)

  # Set counts
  updated_seurat[[counts_assay]] <- SetAssayData(
    updated_seurat[[counts_assay]],
    layer = "counts",
    new.data = raw_counts
  )

  if (!is.null(save_path)) {
    message("Saving updated Seurat object to: ", save_path)
    saveRDS(updated_seurat, file = save_path)
  }

  message("Done: Metadata and counts transferred successfully.")
  return(updated_seurat)
}
