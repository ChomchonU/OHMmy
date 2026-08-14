#' Title
#'
#' @param seurat_obj
#'
#' @returns
#' @export
#'
#' @examples
CleanSeuratReductions <- function(seurat_obj) {
  red_names <- Reductions(seurat_obj)

  for (old_name in red_names) {
    # 1. Generate the clean camelCase name
    words <- unlist(strsplit(old_name, "[._]"))
    if(length(words) > 1) {
      words[-1] <- paste0(toupper(substr(words[-1], 1, 1)), substring(words[-1], 2))
    }
    new_name <- paste(words, collapse = "")

    # Skip if already perfectly formatted
    if (old_name == new_name) next

    new_key <- paste0(new_name, "_")
    old_reduc <- seurat_obj[[old_name]]

    # 2. Extract the raw numerical matrices
    my_embeddings <- Embeddings(old_reduc)
    my_loadings <- Loadings(old_reduc)

    # 3. Rename the columns exactly how Seurat wants them
    colnames(my_embeddings) <- paste0(new_key, 1:ncol(my_embeddings))

    if (!is.null(my_loadings) && ncol(my_loadings) > 0) {
      colnames(my_loadings) <- paste0(new_key, 1:ncol(my_loadings))
    }

    # 4. Build a brand new, error-free DimReduc Object
    new_reduc <- CreateDimReducObject(
      embeddings = my_embeddings,
      loadings = my_loadings,
      assay = DefaultAssay(old_reduc),
      stdev = old_reduc@stdev,
      key = new_key,
      global = old_reduc@global
    )

    # 5. Insert the new one and delete the old broken one
    seurat_obj[[new_name]] <- new_reduc
    seurat_obj[[old_name]] <- NULL

    message(paste("Successfully rebuilt & renamed:", old_name, "->", new_name))
  }

  return(seurat_obj)
}

#------------------------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param dims
#' @param k.param
#' @param reduction
#' @param umap_name
#' @param cluster_prefix
#' @param graph_name
#' @param cluster_resolutions
#' @param final_resolution
#' @param save_path
#' @param force_neighbors
#' @param force_clustering
#' @param plot_dir
#' @param return.model
#' @param plot_format
#' @param width
#' @param height
#' @param dpi
#'
#' @returns
#' @export
#'
#' @examples
ClusterAndUMAP <- function(seurat_obj,
                           sample_name = "Sample",
                           dims = 1:50,
                           k.param = 20,
                           reduction = "integrated.har",
                           umap_name = "umap.har",
                           cluster_prefix = "RNA_snn_res.",
                           graph_name = "RNA_snn",
                           cluster_resolutions = seq(0.1, 2, by = 0.1),
                           final_resolution = 0.7,
                           save_path = NULL,
                           force_neighbors = FALSE,
                           force_clustering = FALSE,
                           plot_dir = "Plots_clustree",
                           return.model = TRUE,
                           plot_format = "jpg",
                           width = 10,
                           height = 15,
                           dpi = 300) {

  require(Seurat)
  require(clustree)
  require(ggplot2)
  require(lubridate)

  # Format and prepare metadata column names
  fac_names <- paste0(cluster_prefix, format(cluster_resolutions, nsmall = 1))
  fac_names <- sub("\\.0$", "", fac_names)
  existing_clusters <- fac_names %in% colnames(seurat_obj@meta.data)
  missing_clusters <- !existing_clusters

  # Step 1: Find Neighbors
  run_neighbors <- force_neighbors || !(graph_name %in% names(seurat_obj@graphs))
  if (run_neighbors) {
    message(" [", sample_name, "] Running FindNeighbors using graph: ", graph_name)
    seurat_obj <- FindNeighbors(
      seurat_obj,
      dims = dims,
      reduction = reduction,
      k.param = k.param,
      graph.name = graph_name
    )
  } else {
    message("[", sample_name, "] SNN graph '", graph_name, "' already exists.")
  }

  if (!(graph_name %in% names(seurat_obj@graphs))) {
    stop("SNN graph '", graph_name, "' not found after FindNeighbors.")
  }

  # Step 2: Find Clusters
  if (any(missing_clusters) || force_clustering) {
    message(" [", sample_name, "] Running FindClusters at resolutions: ", paste(cluster_resolutions, collapse = ", "))
    pb <- txtProgressBar(min = 0, max = length(cluster_resolutions), style = 3)

    for (i in seq_along(cluster_resolutions)) {
      res <- cluster_resolutions[i]
      meta_name <- paste0(cluster_prefix, format(res, nsmall = 1))
      meta_name <- sub("\\.0$", "", meta_name)

      if (force_clustering || !(meta_name %in% colnames(seurat_obj@meta.data))) {
        seurat_obj <- FindClusters(seurat_obj, resolution = res, graph.name = graph_name)
      }

      setTxtProgressBar(pb, i)
    }
    close(pb)
  } else {
    message(" [", sample_name, "] All cluster resolutions already present.")
  }

  # Step 3: Plot Clustree
  clustree_plot <- clustree(seurat_obj, prefix = cluster_prefix) +
    ggtitle(paste("Clustree:", sample_name))

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  plot_dir <- paste0(rtrim(plot_dir, "/"), "/")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_file <- paste0(plot_dir, sample_name, "_clustree_", timestamp, ".", plot_format)

  ggsave(plot_file, clustree_plot, width = width, height = height, dpi = dpi)
  message(" [", sample_name, "] Clustree plot saved to: ", plot_file)

  # Step 4: Final clustering & UMAP
  if (!umap_name %in% names(seurat_obj@reductions)) {
    message("🗺️ [", sample_name, "] Running UMAP + clustering at final resolution = ", final_resolution)
    seurat_obj <- seurat_obj %>%
      FindClusters(resolution = final_resolution, graph.name = graph_name) %>%
      RunUMAP(dims = dims, reduction = reduction, reduction.name = umap_name, return.model = return.model)
  } else {
    message("[", sample_name, "] UMAP '", umap_name, "' already exists.")
  }

  # Step 5: Save Seurat object
  if (!is.null(save_path)) {
    message(" [", sample_name, "] Saving Seurat object to: ", save_path)
    saveRDS(seurat_obj, file = save_path)
  }

  return(list(
    seurat = seurat_obj,
    clustree = clustree_plot,
    plot_file = plot_file
  ))
}

# -------------------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param batch_col
#' @param vars_to_regress
#' @param tcr_bcr_patterns
#' @param reduction_name
#' @param integration_method
#' @param integration_reduction
#' @param dims
#' @param interactive_mode
#' @param elbow_plot_dir
#' @param k.weight
#' @param k.anchor
#' @param k.filter
#' @param k.score
#' @param clustering_resolution
#' @param verbose
#' @param sample_name
#'
#' @returns
#' @export
#'
#' @examples
ProcessSeuratSCT <- function(
    seurat_obj,
    batch_col = "batch",
    vars_to_regress = "pct_counts_mt",
    tcr_bcr_patterns = "^TR[ABDG]|^IG[HKL]",
    reduction_name = "pca.SCT",
    integration_method = "HarmonyIntegration",    # Defaulted to string for if/else logic
    integration_reduction = "integrated.har.SCT",
    dims = 1:50,
    interactive_mode = FALSE,     # Covers PC and k-params
    elbow_plot_dir = NULL,
    k.weight = NULL,              # Default is NULL to trigger defaults or interactive prompts
    k.anchor = NULL,
    k.filter = NULL,
    k.score = NULL,
    clustering_resolution = 0.8,
    verbose = TRUE,
    sample_name = "seurat"
) {
  message("Splitting RNA layer by batch...")
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj[[batch_col]][, 1])

  message("Running SCTransform...")
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = vars_to_regress, verbose = verbose)

  message("Removing TCR/BCR genes from variable features...")
  tcr_bcr_genes <- grep(tcr_bcr_patterns, rownames(seurat_obj), value = TRUE)
  tcr_variable_genes <- intersect(tcr_bcr_genes, VariableFeatures(seurat_obj))
  VariableFeatures(seurat_obj) <- setdiff(VariableFeatures(seurat_obj), tcr_variable_genes)
  message("Removed ", length(tcr_variable_genes), " TCR & BCR genes.")

  # --- BATCH AWARENESS & PCA ---
  batch_counts <- table(seurat_obj[[batch_col]])
  min_batch_cells <- min(batch_counts)
  message("Smallest batch has ", min_batch_cells, " cells.")

  max_pca <- min(50, min_batch_cells - 1)

  message("Running PCA (calculating ", max_pca, " PCs)...")
  seurat_obj <- RunPCA(seurat_obj, reduction.name = reduction_name, npcs = max_pca, verbose = verbose)

  pca_stdev <- Seurat::Stdev(seurat_obj, reduction = reduction_name)
  prop_var <- (pca_stdev^2) / sum(pca_stdev^2)
  cumu_var <- cumsum(prop_var) * 100
  suggested_pcs <- which(cumu_var > 90 & (prop_var * 100) < 5)[1]

  if (is.na(suggested_pcs)) suggested_pcs <- max_pca
  suggested_pcs <- min(suggested_pcs, max_pca)

  # --- SAVE ELBOW PLOT ---
  if (!is.null(elbow_plot_dir)) {
    if (!dir.exists(elbow_plot_dir)) dir.create(elbow_plot_dir, recursive = TRUE)

    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    file_name <- paste0("elbow_plot_SCT_", sample_name, "_", timestamp, ".jpg")
    plot_path <- file.path(elbow_plot_dir, file_name)

    p <- ElbowPlot(seurat_obj, reduction = reduction_name, ndims = max_pca) +
      ggtitle(paste0(sample_name, " - SCT Elbow Plot (Min Batch Cells: ", min_batch_cells, ")"))

    ggsave(filename = plot_path, plot = p, width = 6, height = 4)
    message("Elbow plot saved to: ", plot_path)
  }

  # --- INTERACTIVE PC SELECTION ---
  message("--> Suggested number of PCs: ", suggested_pcs)

  if (interactive_mode && interactive()) {
    user_input <- readline(prompt = paste0("Enter the number of PCs to use (or press Enter to use default 'dims = 1:", max(dims), "'): "))
    if (user_input != "") {
      selected_pc <- as.integer(user_input)
      if (!is.na(selected_pc) && selected_pc > 0) {
        dims <- 1:selected_pc
        message("User override: Setting dims to 1:", selected_pc)
      } else {
        message("Invalid input. Proceeding with manually defined dims: 1:", max(dims))
      }
    } else {
      message("No input provided. Proceeding with manually defined dims: 1:", max(dims))
    }
  }

  if (max(dims) > max_pca) {
    warning("Requested dims (1:", max(dims), ") exceeds the maximum allowed by your smallest batch (", max_pca, "). Adjusting down to 1:", max_pca)
    dims <- 1:max_pca
  }

  # ====================================================================
  # --- INTERACTIVE & DYNAMIC K-PARAMETERS SELECTION ---
  # ====================================================================

  get_k_val <- function(manual_val, param_name, default_val) {
    if (!is.null(manual_val)) return(manual_val)

    if (interactive_mode && interactive()) {
      ans <- readline(prompt = paste0("Enter ", param_name, " (or press Enter for Seurat default ", default_val, "): "))
      if (ans != "") {
        parsed <- as.integer(ans)
        if (!is.na(parsed) && parsed > 0) return(parsed)
        message("Invalid input. Using default: ", default_val)
      }
    }
    return(default_val)
  }

  raw_k_weight <- get_k_val(k.weight, "k.weight", 100)
  raw_k_anchor <- get_k_val(k.anchor, "k.anchor", 5)
  raw_k_filter <- get_k_val(k.filter, "k.filter", 200)
  raw_k_score  <- get_k_val(k.score, "k.score", 30)

  safe_k_weight <- max(1, min(raw_k_weight, min_batch_cells - 1))
  safe_k_anchor <- max(1, min(raw_k_anchor, min_batch_cells - 1))
  safe_k_filter <- max(1, min(raw_k_filter, min_batch_cells - 1))
  safe_k_score  <- max(1, min(raw_k_score, min_batch_cells - 1))

  if (safe_k_weight < raw_k_weight) message(" k.weight dynamically reduced from ", raw_k_weight, " to ", safe_k_weight, " due to small batch size.")
  if (safe_k_anchor < raw_k_anchor) message(" k.anchor dynamically reduced from ", raw_k_anchor, " to ", safe_k_anchor, " due to small batch size.")
  if (safe_k_filter < raw_k_filter) message(" k.filter dynamically reduced from ", raw_k_filter, " to ", safe_k_filter, " due to small batch size.")
  if (safe_k_score < raw_k_score) message(" k.score dynamically reduced from ", raw_k_score, " to ", safe_k_score, " due to small batch size.")

  # ====================================================================

  message("Splitting SCT layer by batch (if not already split)...")
  # Ensure SCT layers are split appropriately based on batch
  if (length(Layers(seurat_obj, assay = "SCT")) == 1) {
    seurat_obj[["SCT"]] <- split(seurat_obj[["SCT"]], f = seurat_obj[[batch_col]][, 1])
  }

  message(" Running integration using ", integration_method,
          " → new reduction: ", integration_reduction)

  DefaultAssay(seurat_obj) <- "SCT"

  if (integration_method == "FastMNNIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = FastMNNIntegration,
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      batch = seurat_obj$batch,
      verbose = verbose
    )
  } else if (integration_method == "RPCAIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = RPCAIntegration,
      normalization.method = "SCT", # Essential for SCT data
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      k.weight = safe_k_weight,
      k.anchor = safe_k_anchor,
      k.filter = safe_k_filter,
      k.score = safe_k_score,
      dims = dims,
      verbose = verbose
    )
  } else if (integration_method == "HarmonyIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = HarmonyIntegration,
      normalization.method = "SCT", # Essential for SCT data
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      k.weight = safe_k_weight,
      verbose = verbose
    )
  } else if (integration_method == "CCAIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = CCAIntegration,
      normalization.method = "SCT", # Essential for SCT data
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      k.weight = safe_k_weight,
      k.anchor = safe_k_anchor,
      k.filter = safe_k_filter,
      k.score = safe_k_score,
      dims = dims,
      verbose = verbose
    )
  } else {
    stop("Unknown integration method. Please choose from: FastMNNIntegration, RPCAIntegration, HarmonyIntegration, or CCAIntegration.")
  }

  message("Joining layers...")
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])

  DefaultAssay(seurat_obj) <- "SCT"
  return(seurat_obj)
}

# -------------------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param batch_col
#' @param vars_to_regress
#' @param tcr_bcr_patterns
#' @param reduction_name
#' @param integration_method
#' @param integration_reduction
#' @param dims
#' @param interactive_mode
#' @param elbow_plot_dir
#' @param k.weight
#' @param k.anchor
#' @param k.filter
#' @param k.score
#' @param clustering_resolution
#' @param verbose
#' @param sample_name
#'
#' @returns
#' @export
#'
#' @examples
ProcessSeuratLOG <- function(
    seurat_obj,
    batch_col = "batch",
    vars_to_regress = NULL,
    tcr_bcr_patterns = "^TR[ABDG]|^IG[HKL]",
    reduction_name = "pca.SCT",
    integration_method = HarmonyIntegration,
    integration_reduction = "integrated.har.SCT",
    dims = 1:30,
    interactive_mode = FALSE,     # CHANGED: Renamed from interactive_pca to cover k-params as well
    elbow_plot_dir = NULL,
    k.weight = NULL,              # Default is NULL so it triggers either defaults or interactive prompts
    k.anchor = NULL,
    k.filter = NULL,
    k.score = NULL,
    clustering_resolution = 1,
    verbose = TRUE,
    sample_name = "seurat"
) {
  message("Splitting RNA layer by batch...")
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj[[batch_col]][, 1])

  message("Running LogNormalization...")
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", verbose = verbose)

  message("Finding variable features...")
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = verbose)

  message("Removing TCR/BCR genes from variable features...")
  tcr_bcr_genes <- grep(tcr_bcr_patterns, rownames(seurat_obj), value = TRUE)
  tcr_variable_genes <- intersect(tcr_bcr_genes, VariableFeatures(seurat_obj))
  VariableFeatures(seurat_obj) <- setdiff(VariableFeatures(seurat_obj), tcr_variable_genes)
  message("Removed ", length(tcr_variable_genes), " TCR & BCR genes.")

  message("Scaling data and regressing variables...")
  seurat_obj <- ScaleData(seurat_obj, vars.to.regress = vars_to_regress, features = VariableFeatures(seurat_obj), verbose = verbose)

  # --- BATCH AWARENESS & PCA ---
  batch_counts <- table(seurat_obj[[batch_col]])
  min_batch_cells <- min(batch_counts)
  message("Smallest batch has ", min_batch_cells, " cells.")

  max_pca <- min(50, min_batch_cells - 1)

  message("Running PCA (calculating ", max_pca, " PCs)...")
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = max_pca, reduction.name = reduction_name, verbose = verbose)

  pca_stdev <- Seurat::Stdev(seurat_obj, reduction = reduction_name)
  prop_var <- (pca_stdev^2) / sum(pca_stdev^2)
  cumu_var <- cumsum(prop_var) * 100
  suggested_pcs <- which(cumu_var > 90 & (prop_var * 100) < 5)[1]

  if (is.na(suggested_pcs)) suggested_pcs <- max_pca
  suggested_pcs <- min(suggested_pcs, max_pca)

  # --- SAVE ELBOW PLOT ---
  if (!is.null(elbow_plot_dir)) {
    if (!dir.exists(elbow_plot_dir)) dir.create(elbow_plot_dir, recursive = TRUE)

    # Generate a clean timestamp (Format: YYYYMMDD_HHMMSS)
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    # Construct the filename with sample name and timestamp
    # If sample_name is "Patient1", output is: elbow_plot_Patient1_20231027_143000.jpg
    file_name <- paste0("elbow_plot_", sample_name, "_", timestamp, ".jpg")
    plot_path <- file.path(elbow_plot_dir, file_name)

    # Create the plot (Added the sample name to the title as well)
    p <- ElbowPlot(seurat_obj, reduction = reduction_name, ndims = max_pca) +
      ggtitle(paste0(sample_name, " - Elbow Plot (Min Batch Cells: ", min_batch_cells, ")"))

    # Save the plot
    ggsave(filename = plot_path, plot = p, width = 6, height = 4)
    message("Elbow plot saved to: ", plot_path)
  }

  # --- INTERACTIVE PC SELECTION ---
  message("--> Suggested number of PCs: ", suggested_pcs)

  if (interactive_mode && interactive()) {
    user_input <- readline(prompt = paste0("Enter the number of PCs to use (or press Enter to use default 'dims = 1:", max(dims), "'): "))
    if (user_input != "") {
      selected_pc <- as.integer(user_input)
      if (!is.na(selected_pc) && selected_pc > 0) {
        dims <- 1:selected_pc
        message("User override: Setting dims to 1:", selected_pc)
      } else {
        message("Invalid input. Proceeding with manually defined dims: 1:", max(dims))
      }
    } else {
      message("No input provided. Proceeding with manually defined dims: 1:", max(dims))
    }
  }

  if (max(dims) > max_pca) {
    warning("Requested dims (1:", max(dims), ") exceeds the maximum allowed by your smallest batch (", max_pca, "). Adjusting down to 1:", max_pca)
    dims <- 1:max_pca
  }

  # ====================================================================
  # --- NEW: INTERACTIVE & DYNAMIC K-PARAMETERS SELECTION ---
  # ====================================================================

  # Helper function to get values (Manual -> Interactive -> Default)
  get_k_val <- function(manual_val, param_name, default_val) {
    if (!is.null(manual_val)) return(manual_val) # Manual input overrides everything

    if (interactive_mode && interactive()) {
      ans <- readline(prompt = paste0("Enter ", param_name, " (or press Enter for Seurat default ", default_val, "): "))
      if (ans != "") {
        parsed <- as.integer(ans)
        if (!is.na(parsed) && parsed > 0) return(parsed)
        message("Invalid input. Using default: ", default_val)
      }
    }
    return(default_val)
  }

  # 1. Gather the requested parameters
  raw_k_weight <- get_k_val(k.weight, "k.weight", 100)
  raw_k_anchor <- get_k_val(k.anchor, "k.anchor", 5)
  raw_k_filter <- get_k_val(k.filter, "k.filter", 200)
  raw_k_score  <- get_k_val(k.score, "k.score", 30)

  # 2. Scale them down safely if the batch size is too small
  safe_k_weight <- max(1, min(raw_k_weight, min_batch_cells - 1))
  safe_k_anchor <- max(1, min(raw_k_anchor, min_batch_cells - 1))
  safe_k_filter <- max(1, min(raw_k_filter, min_batch_cells - 1))
  safe_k_score  <- max(1, min(raw_k_score, min_batch_cells - 1))

  # 3. Warn the user if scaling occurred
  if (safe_k_weight < raw_k_weight) message(" k.weight dynamically reduced from ", raw_k_weight, " to ", safe_k_weight, " due to small batch size.")
  if (safe_k_anchor < raw_k_anchor) message(" k.anchor dynamically reduced from ", raw_k_anchor, " to ", safe_k_anchor, " due to small batch size.")
  if (safe_k_filter < raw_k_filter) message(" k.filter dynamically reduced from ", raw_k_filter, " to ", safe_k_filter, " due to small batch size.")
  if (safe_k_score < raw_k_score) message(" k.score dynamically reduced from ", raw_k_score, " to ", safe_k_score, " due to small batch size.")
  # ====================================================================

  message(" Running integration using ", integration_method,
          " → new reduction: ", integration_reduction)

  DefaultAssay(seurat_obj) <- "RNA"

  if (integration_method == "FastMNNIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = FastMNNIntegration,
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      batch = seurat_obj$batch,
      verbose = verbose
    )
  } else if (integration_method == "RPCAIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = RPCAIntegration,
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      k.weight = safe_k_weight, # Using safe variables
      k.anchor = safe_k_anchor,
      k.filter = safe_k_filter,
      k.score = safe_k_score,
      dims = dims,
      verbose = verbose
    )
  } else if (integration_method == "HarmonyIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = HarmonyIntegration,
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      k.weight = safe_k_weight, # Using safe variable
      verbose = verbose
    )
  } else if (integration_method == "CCAIntegration") {
    seurat_obj <- IntegrateLayers(
      object = seurat_obj,
      method = CCAIntegration,
      orig.reduction = reduction_name,
      new.reduction = integration_reduction,
      k.weight = safe_k_weight, # Using safe variables
      k.anchor = safe_k_anchor,
      k.filter = safe_k_filter,
      k.score = safe_k_score,
      dims = dims,
      verbose = verbose
    )
  } else {
    stop("Unknown integration method. Please choose from: FastMNN, RPCA, Harmony, or CCA.")
  }

  message("Joining layers...")
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])

  DefaultAssay(seurat_obj) <- "RNA"
  return(seurat_obj)
}
# -------------------------------------------------------------
