#' Title
#'
#' @param seurat_obj
#' @param cluster_col
#' @param genes_list
#' @param sample_name
#' @param reduction
#' @param viridis_palette
#' @param output_dir
#' @param save_format
#' @param width
#' @param dpi
#' @param add_timestamp
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
plot_combined <- function(seurat_obj,
                          cluster_col,
                          genes_list , # Changed from genes_list to reflect it takes metadata too
                          sample_name = "Sample",
                          reduction = "umap.cca",
                          viridis_palette = "viridis",
                          output_dir = "Plots_feat_dens",
                          save_format = "jpg",  # Options: "png", "pdf"
                          width = 20,
                          dpi = 300,
                          add_timestamp = TRUE,
                          verbose = TRUE) {
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  require(lubridate)

  chains <- names(genes_list )
  combined_plot_list <- list()

  Idents(seurat_obj) <- cluster_col

  # Create output folder
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), '%Y-%m-%d_%H-%M-%S') else NULL

  if (verbose) {
    message("Plotting feature & density for sample: ", sample_name)
    pb <- txtProgressBar(min = 0, max = length(chains), style = 3)
  }

  for (i in seq_along(chains)) {
    chain <- chains[i]
    features <- genes_list[[chain]]

    # MODIFIED: Check for features in both gene names AND metadata columns
    available_features <- c(rownames(seurat_obj), colnames(seurat_obj@meta.data))
    features_in_data <- intersect(features, available_features)

    if (length(features_in_data) == 0) {
      warning("No valid features (genes or metadata) found for '", chain, "' in Seurat object.")
      if (verbose) setTxtProgressBar(pb, i)
      next
    }

    # Dimensions for saving
    n_features <- length(features_in_data)

    # --- Build interleaved plots (feature | density) per feature
    interleaved_plots <- lapply(features_in_data, function(feat) {
      p_feat <- FeaturePlot(
        seurat_obj,
        reduction = reduction,
        features = feat,
        min.cutoff = "q10",
        label = TRUE,
        repel = TRUE,
        order = TRUE
      ) + ggtitle(paste(feat, "– Feature"))

      p_dens <- Plot_Density_Custom(
        seurat_obj,
        reduction = reduction,
        features = feat,
        viridis_palette = viridis_palette
      ) + ggtitle(paste(feat, "– Density"))

      list(p_feat, p_dens)
    })

    # Flatten: [feat1, dens1, feat2, dens2, feat3, dens3, ...]
    plot_list_flat <- unlist(interleaved_plots, recursive = FALSE)

    # 6 cols → 3 features per row; each feature takes 2 slots (feat + dens)
    ncol   <- 6
    nrow   <- ceiling(n_features / 3)
    height <- nrow * (width / ncol) * 2  # approximate square panels

    p_combined <- wrap_plots(plot_list_flat, ncol = ncol) +
      plot_annotation(
        title = paste(sample_name, chain, "Feature + Density")
      )

    # --- Save to file
    reduction_safe <- sanitize(reduction)
    file_combined <- file.path(output_dir, paste0(
      paste(
        sanitize(sample_name),
        "CombinedPlot",
        sanitize(chain),
        reduction_safe,
        if (!is.null(timestamp)) timestamp else NULL,
        sep = "_"
      ),
      ".", save_format
    ))

    ggsave(file_combined, p_combined, width = width * 2, height = height, dpi = dpi, limitsize = FALSE)

    # Store the combined plot
    combined_plot_list[[chain]] <- p_combined

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("\nFinished saving combined plots to: ", output_dir)
  }

  return(list(
    combined_plots = combined_plot_list,
    output_dir = output_dir
  ))
}

# Helper to clean filenames
sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

# ----------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param qc_meta_features
#' @param gene_features
#' @param res_col
#' @param assay
#' @param output_dir
#' @param save_format
#' @param width
#' @param dpi
#' @param add_timestamp
#'
#' @returns
#' @export
#'
#' @examples
plot_violin_qc_single <- function(seurat_obj,
                                  sample_name = "Sample",
                                  qc_meta_features = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt"),
                                  gene_features = NULL,
                                  res_col = "seurat_clusters",
                                  assay = "RNA",
                                  output_dir = "Plots_violin_qc",
                                  save_format = "jpg",
                                  width = 20,
                                  dpi = 300,
                                  add_timestamp = TRUE) {
  require(Seurat)
  require(ggplot2)
  require(patchwork)
  require(lubridate)

  # Helper to sanitize file components
  sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

  # Timestamp (optional)
  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), "%Y-%m-%d_%H-%M-%S") else NULL

  # Ensure output directory
  output_dir <- sub("/+$", "", output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (!(res_col %in% colnames(seurat_obj@meta.data))) {
    warning("Cluster resolution '", res_col, "' not found in sample ", sample_name)
    return(NULL)
  }

  Idents(seurat_obj) <- res_col
  message("Sample: ", sample_name, " | Cluster column: ", res_col)

  meta_cols <- colnames(seurat_obj@meta.data)
  message("Meta.data columns: ", paste(head(meta_cols, 10), collapse = ", "), " ...")

  # --- QC Feature Violin Plots ---
  qc_found <- qc_meta_features[qc_meta_features %in% meta_cols]
  if (length(qc_found) > 0) {
    message("QC features found: ", paste(qc_found, collapse = ", "))

    qc_plots <- lapply(qc_found, function(feature) {
      VlnPlot(seurat_obj, features = feature, pt.size = 0) +
        ggtitle(paste0("QC: ", feature, " | Clusters: ", res_col)) +
        theme(plot.title = element_text(size = 10))+
        geom_jitter(
          width = 0.2,       # Spread the dots out horizontally
          alpha = 0.25,      # 2. Fade the dots! (0 is invisible, 1 is solid)
          size = 0.5         # Make the dots smaller
        )
    })

    p_qc <- wrap_plots(qc_plots, ncol = 3)

    # File naming (structured like PlotDimByFactors)
    fname_parts <- c(sanitize(sample_name), "ViolinQC", sanitize(res_col))
    if (!is.null(timestamp)) fname_parts <- c(fname_parts, timestamp)
    file_name <- paste0(paste(fname_parts, collapse = "_"), ".", save_format)
    file_path <- file.path(output_dir, file_name)

    message("Saving QC plot to: ", file_path)

    n_qc <- length(qc_found)
    ncol <- 3
    nrow <- ceiling(n_qc / ncol)
    row_height <- 20 / 3
    height <- nrow * row_height

    tryCatch({
      ggsave(file_path, plot = p_qc, width = width, height = height, dpi = dpi, limitsize = FALSE)
      if (file.exists(file_path)) {
        message("QC plot saved: ", file_path)
      } else {
        warning("ggsave failed. Trying fallback PDF.")
        fallback_file <- sub(paste0("\\.", save_format, "$"), ".pdf", file_path)
        pdf(fallback_file, width = width, height = height)
        print(p_qc)
        dev.off()
        message("Fallback PDF saved: ", fallback_file)
      }
    }, error = function(e) {
      warning("Error saving QC plot: ", conditionMessage(e))
    })
  } else {
    warning("No QC metadata features found.")
  }

  # --- Gene Expression Violin Plots ---
  if (!is.null(gene_features) && length(gene_features) > 0) {
    DefaultAssay(seurat_obj) <- assay

    gene_plots <- list()
    for (gene in gene_features) {
      if (!(gene %in% rownames(seurat_obj))) {
        warning("Gene ", gene, " not found.")
        next
      }

      gene_plots[[gene]] <- VlnPlot(seurat_obj, features = gene, pt.size = 0) +
        ggtitle(paste0("Gene: ", gene, " | Clusters: ", res_col))+
        geom_jitter(
          width = 0.2,       # Spread the dots out horizontally
          alpha = 0.25,      # 2. Fade the dots! (0 is invisible, 1 is solid)
          size = 0.5         # Make the dots smaller
        )
    }

    if (length(gene_plots) > 0) {
      # Combine into 2-column layout
      p_gene <- wrap_plots(gene_plots, ncol = 3)

      fname_parts <- c(sanitize(sample_name), "ViolinGene", sanitize(res_col))
      if (!is.null(timestamp)) fname_parts <- c(fname_parts, timestamp)
      file_name <- paste0(paste(fname_parts, collapse = "_"), ".", save_format)
      file_path <- file.path(output_dir, file_name)

      message("Saving combined gene plots: ", file_path)

      n_gene <- length(gene_features)
      ncol <- 3
      nrow <- ceiling(n_gene / ncol)
      row_height <- 20 / 3
      height <- nrow * row_height

      tryCatch({
        ggsave(file_path, plot = p_gene,
               width = 20,  # fixed 2 cols (2 × 8)
               height = height,
               dpi = dpi, limitsize = FALSE)
        if (file.exists(file_path)) {
          message("Gene plots saved: ", file_path)
        } else {
          warning("ggsave failed. Trying fallback PDF.")
          fallback_file <- sub(paste0("\\.", save_format, "$"), ".pdf", file_path)
          pdf(fallback_file, width = 16, height = height * ceiling(length(gene_plots) / 2))
          print(p_gene)
          dev.off()
          message(" Fallback PDF saved: ", fallback_file)
        }
      }, error = function(e) {
        warning("Error saving gene plots: ", conditionMessage(e))
      })
    }
  } else {
    message(" No gene features provided.")
  }

  message("Done: ", sample_name)
}

# ------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param factors
#' @param sample_name
#' @param reduction
#' @param raster
#' @param label
#' @param repel
#' @param output_dir
#' @param plot_format
#' @param width
#' @param height
#' @param dpi
#' @param add_timestamp
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
PlotDimByFactors <- function(
    seurat_obj,
    factors,                        # <-- moved up: now 2nd arg
    sample_name = "Sample",
    reduction = "umap.har",
    raster = FALSE,
    label = TRUE,
    repel = TRUE,
    output_dir = "Plots_umap",
    plot_format = "jpg",            # "png" or "pdf"
    width = 10,
    height = 10,
    dpi = 300,
    add_timestamp = TRUE,           # new: turn timestamp on/off
    verbose = TRUE                  # new: message control
) {
  stopifnot(inherits(seurat_obj, "Seurat"))

  if (missing(factors) || length(factors) == 0) {
    stop("`factors` must be supplied (vector of metadata column names).")
  }

  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Reduction '", reduction, "' not found in seurat_obj@reductions. Available: ",
         paste(names(seurat_obj@reductions), collapse = ", "))
  }

  md_cols <- colnames(seurat_obj@meta.data)
  missing_fac <- setdiff(factors, md_cols)
  if (length(missing_fac) > 0) {
    warning("These factors not found in metadata and will be skipped: ",
            paste(missing_fac, collapse = ", "))
    factors <- intersect(factors, md_cols)
    if (length(factors) == 0) {
      stop("None of the requested factors are present in metadata.")
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), "%Y-%m-%d_%H-%M-%S") else NULL

  # filename-sanitizing helper
  sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

  plots <- list()
  n <- length(factors)

  if (verbose) {
    message("Generating ", n, " DimPlots for ", sample_name, " in ", output_dir)
    pb <- txtProgressBar(min = 0, max = n, style = 3)
  }

  for (i in seq_along(factors)) {
    fac <- factors[i]

    p <- DimPlot(
      seurat_obj,
      reduction = reduction,
      group.by = fac,
      raster = raster,
      label = label,
      repel = repel
    ) +
      ggtitle(paste0(sample_name, " – ", fac))

    plots[[fac]] <- p
    if (verbose) print(p)

    # Build filename
    fname_parts <- c(sanitize(sample_name), "DimPlot", sanitize(fac), sanitize(reduction))
    if (!is.null(timestamp)) fname_parts <- c(fname_parts, timestamp)
    file_name <- paste0(paste(fname_parts, collapse = "_"), ".", plot_format)
    file_path <- file.path(output_dir, file_name)

    # Save
    tryCatch(
      {
        ggsave(
          filename = file_path,
          plot = p,
          width = width,
          height = height,
          dpi = dpi,
          device = plot_format
        )
      },
      error = function(e) {
        warning("Failed to save plot for factor ", fac, ": ", e$message)
      }
    )

    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) {
    close(pb)
    message("UMAP plots for ", sample_name, " saved in: ", output_dir)
  }

  invisible(plots)
}

# -------------------------------------------

# ── Palette helper (add this alongside sanitize()) ────────────────────────────
#' Title
#'
#' @param n
#' @param palette
#'
#' @returns
#' @export
#'
#' @examples
make_cluster_palette <- function(n, palette = "auto") {

  # --- Curated 36-colour palette combining Kelly + Alphabet + extras ----------
  # Maximally distinct, works well on white and dark backgrounds
  kelly_cols <- c(
    "#F3C300","#875692","#F38400","#A1CAF1","#BE0032","#C2B280",
    "#848482","#008856","#E68FAC","#0067A5","#F99379","#604E97",
    "#F6A600","#B3446C","#DCD300","#882D17","#8DB600","#654522",
    "#E25822","#2B3D26","#F2F3F4","#222222","#F3C300","#875692"
  )

  discrete_palettes <- list(
    "kelly"      = kelly_cols,
    "alphabet"   = c("#AA0DFE","#3283FE","#85660D","#782AB6","#565656",
                     "#1C8356","#16FF32","#F7E1A0","#E2E2E2","#1CBE4F",
                     "#C4451C","#DEA0FD","#FE00FA","#325A9B","#FEAF16",
                     "#F8A19F","#90AD1C","#F6222E","#1CFFCE","#2ED9FF",
                     "#B10DA1","#C075A6","#FC1CBF","#B00068","#FBE426",
                     "#FA0087"),
    "polychrome" = c("#5A5156","#E4E1E3","#F6222E","#FE00FA","#16FF32",
                     "#3283FE","#FEAF16","#B00068","#1CFFCE","#90AD1C",
                     "#2ED9FF","#AA0DFE","#F8A19F","#325A9B","#C4451C",
                     "#1CBE4F","#DEA0FD","#F6222E","#1C8356","#85660D",
                     "#B10DA1","#FBE426","#FA0087","#C075A6","#782AB6",
                     "#565656","#FC1CBF","#3283FE","#FEAF16","#565656"),
    "auto"       = NULL  # triggers the hybrid approach below
  )

  if (palette == "auto" || is.null(discrete_palettes[[palette]])) {
    # Combine multiple RColorBrewer qualitative palettes then interpolate if needed
    base_cols <- unique(c(
      RColorBrewer::brewer.pal(8,  "Set2"),
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(8,  "Dark2"),
      RColorBrewer::brewer.pal(9,  "Paired")
    ))  # 37 unique colours before interpolation
    if (n > length(base_cols)) {
      base_cols <- colorRampPalette(base_cols)(n)
    }
    return(base_cols[seq_len(n)])
  }

  cols <- discrete_palettes[[palette]]
  if (n > length(cols)) {
    cols <- colorRampPalette(cols)(n)  # interpolate if still not enough
  }
  cols[seq_len(n)]
}

#----------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param cluster_col
#' @param gene_pairs
#' @param sample_name
#' @param cor_method
#' @param quantile_thresh
#' @param fill_palette
#' @param output_dir
#' @param save_format
#' @param width
#' @param height
#' @param dpi
#' @param add_timestamp
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
plot_gene_pair_correlations <- function(seurat_obj,
                                        cluster_col,
                                        gene_pairs,
                                        sample_name = "Sample",
                                        cor_method = "pearson",
                                        quantile_thresh = 0.01,
                                        fill_palette = "auto",   # ← was "Set2"
                                        output_dir = "Plots_gene_pair_cor",
                                        save_format = "jpg",
                                        width = 14,
                                        height = 8,
                                        dpi = 300,
                                        add_timestamp = TRUE,
                                        verbose = TRUE) {

  require(Seurat)
  require(ggplot2)
  require(dplyr)
  require(lubridate)

  # ── Validation ────────────────────────────────────────────────────────────────
  if (!is.list(gene_pairs) || !all(lengths(gene_pairs) == 2)) {
    stop("`gene_pairs` must be a list of character vectors, each of length 2.")
  }
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("`cluster_col` '", cluster_col, "' not found in Seurat object metadata.")
  }

  # ── Setup ─────────────────────────────────────────────────────────────────────
  Idents(seurat_obj) <- cluster_col
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), "%Y-%m-%d_%H-%M-%S") else NULL

  all_genes    <- unique(unlist(gene_pairs))
  genes_found  <- intersect(all_genes, rownames(seurat_obj))
  genes_missing <- setdiff(all_genes, genes_found)

  if (length(genes_missing) > 0) {
    warning("  The following genes were not found and will be skipped: ",
            paste(genes_missing, collapse = ", "))
  }

  # Drop pairs where either gene is missing
  valid_pairs <- Filter(function(p) all(p %in% genes_found), gene_pairs)

  if (length(valid_pairs) == 0) {
    stop("No valid gene pairs remaining after filtering for genes present in the object.")
  }

  if (verbose) {
    message(" Running gene-pair correlations for sample: ", sample_name)
    message("   Valid pairs  : ", length(valid_pairs), " / ", length(gene_pairs))
    message("   Cluster col  : ", cluster_col)
    message("   Cor method   : ", cor_method)
    message("   Quantile thr : ", quantile_thresh)
    pb <- txtProgressBar(min = 0, max = length(valid_pairs), style = 3)
  }

  # ── Fetch expression + cluster data once ──────────────────────────────────────
  fetch_genes <- unique(unlist(valid_pairs))
  exp_data    <- FetchData(seurat_obj, vars = c(fetch_genes, cluster_col))
  colnames(exp_data)[ncol(exp_data)] <- "Cluster"

  # ── Loop over pairs ───────────────────────────────────────────────────────────
  cor_results <- list()

  for (i in seq_along(valid_pairs)) {
    gene1     <- valid_pairs[[i]][1]
    gene2     <- valid_pairs[[i]][2]
    pair_name <- paste(gene1, "vs", gene2)

    thresh1 <- quantile(exp_data[[gene1]], quantile_thresh, na.rm = TRUE)
    thresh2 <- quantile(exp_data[[gene2]], quantile_thresh, na.rm = TRUE)

    cor_df <- exp_data %>%
      filter(.data[[gene1]] > thresh1 & .data[[gene2]] > thresh2) %>%
      group_by(Cluster) %>%
      summarize(
        Correlation = cor(.data[[gene1]], .data[[gene2]],
                          method = cor_method, use = "complete.obs"),
        N_Cells     = n(),
        Pair        = pair_name,
        .groups     = "drop"
      )

    cor_results[[i]] <- cor_df
    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  final_cor_df <- dplyr::bind_rows(cor_results)

  # ── Plot ───────────────────────────────────────────────────────────────────────
  n_clusters   <- length(unique(final_cor_df$Cluster))
  cluster_cols <- make_cluster_palette(n_clusters, palette = fill_palette)

  p <- ggplot(final_cor_df, aes(x = Cluster, y = Correlation, fill = Cluster)) +
    geom_bar(stat = "identity", color = "black", width = 0.6) +
    geom_text(
      aes(label = round(Correlation, 3),
          vjust = ifelse(Correlation >= 0, -0.5, 1.5)),
      size = 3.5
    ) +
    facet_wrap(~ Pair, scales = "free_y") +
    scale_fill_manual(values = cluster_cols) +
    theme_minimal() +
    labs(
      title    = paste(sample_name, "– Gene Pair Correlations by Cluster"),
      subtitle = paste0("Method: ", cor_method,
                        "  |  Quantile threshold: ", quantile_thresh),
      x = "Cluster",
      y = paste(tools::toTitleCase(cor_method), "Correlation Coefficient")
    ) +
    theme(
      legend.position  = "none",
      plot.title       = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle    = element_text(hjust = 0.5),
      strip.text       = element_text(size = 12, face = "bold"),
      axis.text        = element_text(size = 10),
      axis.title       = element_text(size = 12),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  # ── Save ───────────────────────────────────────────────────────────────────────
  file_out <- file.path(output_dir, paste0(
    paste(
      sanitize(sample_name),
      "GenePairCor",
      sanitize(cluster_col),
      if (!is.null(timestamp)) timestamp else NULL,
      sep = "_"
    ),
    ".", save_format
  ))

  ggsave(file_out, p, width = width, height = height, dpi = dpi)

  if (verbose) {
    message("Plot saved to: ", file_out)
  }

  return(list(          # ← must be AFTER ggsave and the verbose message
    plot       = p,
    data       = final_cor_df,
    output_dir = output_dir
  ))
}                       # ← closing brace of the function

#------------------------------------------

# ── Private helper: Seurat-v5-safe blend plot ──────────────────────────────────
# Returns a patchwork of 4 panels (gene1 | gene2 | blend | colour-key),
# matching the shape that FeaturePlot(blend=TRUE, combine=TRUE) used to give.
#' Title
#'
#' @param seurat_obj
#' @param gene1
#' @param gene2
#' @param reduction
#' @param cols
#' @param blend_threshold
#' @param min_cutoff
#' @param gamma
#' @param pt_size
#'
#' @returns
#' @export
#'
#' @examples
.blend_feature_plot_v5 <- function(seurat_obj, gene1, gene2, reduction,
                                   cols            = c("lightgrey", "#00ff00", "#ff0000"),
                                   blend_threshold = 0.1,
                                   min_cutoff      = "q10",
                                   gamma           = 1,
                                   pt_size         = 1) {

  # ── Coordinates ───────────────────────────────────────────────────────────────
  emb       <- Embeddings(seurat_obj, reduction = reduction)
  dim_names <- colnames(emb)[1:2]
  df        <- as.data.frame(emb[, 1:2])
  colnames(df) <- c("DIM1", "DIM2")

  # ── Expression (v5 uses layer=, v4 used slot=) ────────────────────────────────
  expr_mat <- tryCatch(
    GetAssayData(seurat_obj, layer = "data"),
    error = function(e) GetAssayData(seurat_obj, slot = "data")
  )

  e1 <- as.numeric(expr_mat[gene1, rownames(df)])
  e2 <- as.numeric(expr_mat[gene2, rownames(df)])

  # ── Apply quantile min-cutoff ─────────────────────────────────────────────────
  apply_cutoff <- function(x, co) {
    if (is.character(co) && grepl("^q", co)) {
      q      <- as.numeric(sub("q", "", co)) / 100
      thresh <- quantile(x, q, na.rm = TRUE)
      x[x < thresh] <- 0
    }
    x
  }
  e1 <- apply_cutoff(e1, min_cutoff)
  e2 <- apply_cutoff(e2, min_cutoff)

  # ── Scale each gene to [0, 1] ─────────────────────────────────────────────────
  scale01 <- function(x) { mx <- max(x, na.rm = TRUE); if (mx == 0) x else x / mx }
  s1 <- scale01(e1)
  s2 <- scale01(e2)

  # ── Remap: values below threshold → 0, threshold → 0, max → 1 ────────────────
  remap <- function(x) pmax(0, (x - blend_threshold) / (1 - blend_threshold))
  s1r <- remap(s1)
  s2r <- remap(s2)

  # ── Quantize remapped values to 10 steps ─────────────────────────────────────
  steps <- 10L
  snap  <- function(x) ceiling(x * steps) / steps
  s1q <- snap(s1r)
  s2q <- snap(s2r)

  # ── Colour setup ──────────────────────────────────────────────────────────────
  col_neither <- cols[1]
  col_g1      <- cols[2]
  col_g2      <- cols[3]
  rgb_neither <- col2rgb(col_neither) / 255
  rgb_g1      <- col2rgb(col_g1)      / 255
  rgb_g2      <- col2rgb(col_g2)      / 255

  # ── Blend function ────────────────────────────────────────────────────────────
  blend_one <- function(v1, v2) {
    if (v1 == 0 && v2 == 0) return(col_neither)

    mixed <- v1 * rgb_g1 + v2 * rgb_g2
    mx <- max(mixed)
    if (mx > 0) mixed <- mixed / mx

    brightness <- max(v1, v2) ^ gamma
    final <- rgb_neither * (1 - brightness) + mixed * brightness
    final <- pmax(0, pmin(1, final))
    rgb(final[1], final[2], final[3])
  }

  # ── Populate data frame ───────────────────────────────────────────────────────
  df$gene1     <- s1r
  df$gene2     <- s2r
  df$gene1q    <- s1q
  df$gene2q    <- s2q

  df$blend_col <- mapply(blend_one, s1q, s2q, SIMPLIFY = TRUE)
  df$expressed <- (s1q > 0) | (s2q > 0)
  df <- df[order(df$expressed, s1q + s2q), ]   # background cells plotted first

  # ── Shared theme ──────────────────────────────────────────────────────────────
  base_thm <- theme_classic(base_size = 10) + theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position   = "none",
    axis.title        = element_text(size = 8)
  )
  ax <- labs(x = dim_names[1], y = dim_names[2])

  # ── Panel 1: gene 1 ──────────────────────────────────────────────────────────
  p1 <- ggplot(df[order(df$gene1), ], aes(x = DIM1, y = DIM2, colour = gene1)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradient(low = col_neither, high = col_g1, limits = c(0, 1)) +
    labs(title = gene1) + ax + base_thm

  # ── Panel 2: gene 2 ──────────────────────────────────────────────────────────
  p2 <- ggplot(df[order(df$gene2), ], aes(x = DIM1, y = DIM2, colour = gene2)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradient(low = col_neither, high = col_g2, limits = c(0, 1)) +
    labs(title = gene2) + ax + base_thm

  # ── Panel 3: blend ────────────────────────────────────────────────────────────
  p_bl <- ggplot(df, aes(x = DIM1, y = DIM2)) +
    geom_point(colour = df$blend_col, size = pt_size, stroke = 0) +
    labs(title = paste(gene1, "+", gene2)) + ax + base_thm

  # ── Panel 4: Colour key (Matching full 0-1 scale with threshold mask) ─────────
  # Generate a high-res grid over the ORIGINAL 0-1 scale
  leg_vals <- seq(0, 1, length.out = 100)
  leg_df   <- expand.grid(g1 = leg_vals, g2 = leg_vals)

  # Push legend grid through the exact same remap & snap filters as the cells
  leg_df$g1_q <- snap(remap(leg_df$g1))
  leg_df$g2_q <- snap(remap(leg_df$g2))
  leg_df$col  <- mapply(blend_one, leg_df$g1_q, leg_df$g2_q, SIMPLIFY = TRUE)

  p_key <- ggplot(leg_df, aes(x = g1, y = g2)) +
    geom_tile(fill = leg_df$col) +
    # Add the crosshairs from your image
    geom_vline(xintercept = blend_threshold, linetype = "dashed", colour = "black", alpha = 0.4) +
    geom_hline(yintercept = blend_threshold, linetype = "dashed", colour = "black", alpha = 0.4) +
    scale_x_continuous(expand = c(0,0), breaks = c(0, blend_threshold, 1), labels = c("0", "Thr", "Max")) +
    scale_y_continuous(expand = c(0,0), breaks = c(0, blend_threshold, 1), labels = c("0", "Thr", "Max")) +
    labs(title = paste("Color Threshold:", blend_threshold), x = gene1, y = gene2) +
    theme_classic(base_size = 10) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.title   = element_text(size = 8),
      aspect.ratio = 1, # Keeps the legend perfectly square
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )

  wrap_plots(p1, p2, p_bl, p_key, nrow = 1)
}

#--------------------------------------------------------------------------------------

# ── Main function (only the FeaturePlot call is changed) ───────────────────────
#' Title
#'
#' @param seurat_obj
#' @param cluster_col
#' @param gene_pairs
#' @param sample_name
#' @param reduction
#' @param blend_cols
#' @param blend_threshold
#' @param min_cutoff
#' @param output_dir
#' @param save_format
#' @param width
#' @param dpi
#' @param add_timestamp
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
plot_blend_nebulosa <- function(seurat_obj,
                                cluster_col,
                                gene_pairs,
                                sample_name     = "Sample",
                                reduction       = "umap.cca",
                                blend_cols      = c("lightgrey", "#00ff00", "#ff0000"),
                                blend_threshold = 0.1,
                                min_cutoff      = "q10",
                                output_dir      = "Plots_blend_nebulosa",
                                save_format     = "jpg",
                                width           = 36,
                                dpi             = 300,
                                add_timestamp   = TRUE,
                                verbose         = TRUE) {

  require(Seurat)
  require(ggplot2)
  require(patchwork)
  require(Nebulosa)

  # ── Validation ────────────────────────────────────────────────────────────────
  if (!is.list(gene_pairs) || !all(lengths(gene_pairs) == 2))
    stop("`gene_pairs` must be a list of character vectors, each of length 2.")
  if (!cluster_col %in% colnames(seurat_obj@meta.data))
    stop("`cluster_col` '", cluster_col, "' not found in Seurat object metadata.")

  # ── Setup ─────────────────────────────────────────────────────────────────────
  Idents(seurat_obj) <- cluster_col
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), "%Y-%m-%d_%H-%M-%S") else NULL

  genes_missing <- setdiff(unique(unlist(gene_pairs)), rownames(seurat_obj))
  if (length(genes_missing) > 0)
    warning("  Genes not found and will be skipped: ",
            paste(genes_missing, collapse = ", "))

  valid_pairs <- Filter(function(p) all(p %in% rownames(seurat_obj)), gene_pairs)
  if (length(valid_pairs) == 0) stop("No valid gene pairs after filtering.")

  # ── Layout math ───────────────────────────────────────────────────────────────
  n_panels_per_row <- 6
  panel_width      <- width / n_panels_per_row
  panel_height     <- panel_width * 0.9
  total_height     <- panel_height * length(valid_pairs)

  if (verbose) {
    message(" Plotting all pairs in one figure for sample: ", sample_name)
    message("   Layout   : ", length(valid_pairs), " rows × ", n_panels_per_row, " panels")
    message("   Width    : ", width, "in  |  Total height: ", round(total_height, 2), "in")
    pb <- txtProgressBar(min = 0, max = length(valid_pairs), style = 3)
  }

  all_rows <- list()

  for (i in seq_along(valid_pairs)) {
    gene1     <- valid_pairs[[i]][1]
    gene2     <- valid_pairs[[i]][2]
    pair_name <- paste(gene1, "vs", gene2)

    # ── CHANGED: use v5-safe manual blend instead of FeaturePlot ─────────────
    p_blend <- .blend_feature_plot_v5(
      seurat_obj      = seurat_obj,
      gene1           = gene1,
      gene2           = gene2,
      reduction       = reduction,
      cols            = blend_cols,
      blend_threshold = blend_threshold,
      min_cutoff      = min_cutoff
    )

    # ── Nebulosa — gene 1 ────────────────────────────────────────────────────────
    p_neb1 <- plot_density(seurat_obj, features = gene1, reduction = reduction) +
      labs(title = paste(gene1, "density")) +
      theme(
        plot.title        = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position   = "right",
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.2, "cm"),
        legend.title      = element_text(size = 8),
        legend.text       = element_text(size = 7)
      )

    # ── Nebulosa — gene 2 ────────────────────────────────────────────────────────
    p_neb2 <- plot_density(seurat_obj, features = gene2, reduction = reduction) +
      labs(title = paste(gene2, "density")) +
      theme(
        plot.title        = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position   = "right",
        legend.key.height = unit(0.4, "cm"),
        legend.key.width  = unit(0.2, "cm"),
        legend.title      = element_text(size = 8),
        legend.text       = element_text(size = 7)
      )

    # ── Assemble row ──────────────────────────────────────────────────────────────
    p_row <- wrap_plots(
      p_blend, p_neb1, p_neb2,
      nrow   = 1,
      widths = c(4, 1, 1)
    ) + plot_annotation(
      title = pair_name,
      theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 12))
    )

    all_rows[[i]] <- p_row
    if (verbose) setTxtProgressBar(pb, i)
  }

  if (verbose) close(pb)

  # ── Stack all rows ─────────────────────────────────────────────────────────────
  p_final <- wrap_plots(all_rows, ncol = 1) +
    plot_annotation(
      title    = paste(sample_name, "— Blend + Density"),
      subtitle = paste0("Reduction: ", reduction,
                        "  |  Blend threshold: ", blend_threshold,
                        "  |  Min cutoff: ", min_cutoff),
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40")
      )
    )

  # ── Save ──────────────────────────────────────────────────────────────────────
  file_out <- file.path(output_dir, paste0(
    paste(
      sanitize(sample_name),
      "BlendNebulosa",
      sanitize(reduction),
      if (!is.null(timestamp)) timestamp else NULL,
      sep = "_"
    ),
    ".", save_format
  ))

  ggsave(file_out, p_final,
         width     = width,
         height    = total_height,
         dpi       = dpi,
         limitsize = FALSE)

  if (verbose) message("Saved to: ", file_out)

  return(list(
    plot       = p_final,
    rows       = all_rows,
    output_dir = output_dir,
    dimensions = list(width = width, height = total_height)
  ))
}

#---------------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param cluster_col
#' @param batch_col
#' @param output_dir
#' @param file_prefix
#'
#' @returns
#' @export
#'
#' @examples
plot_cluster_distributions <- function(
    seurat_obj,
    cluster_col = "label_T3_log_rpca",
    batch_col = "exist_GTS",
    output_dir = "C:/Users/ADMIN/Desktop/Dengue_summer/Output_R/Plots_count_dis/CD8_subset/rpca",
    file_prefix = "ind_label"
) {

  # Ensure columns exist in metadata
  meta_data <- seurat_obj@meta.data
  if (!(cluster_col %in% colnames(meta_data)) | !(batch_col %in% colnames(meta_data))) {
    stop("One or both specified columns not found in the Seurat object's metadata.")
  }

  # Step 1: Convert to data frame
  df_table <- as.data.frame(table(
    meta_data[[cluster_col]],
    meta_data[[batch_col]]
  ))
  # Standardize column names so ggplot can easily read them generically
  colnames(df_table) <- c("Cluster", "Batch", "Freq")

  # Step 2: Compute proportions
  df_table_prop <- df_table %>%
    group_by(Batch) %>%
    mutate(prop_batch = Freq / sum(Freq)) %>%
    ungroup() %>%
    group_by(Cluster) %>%
    mutate(prop_level = Freq / sum(Freq)) %>%
    ungroup()

  # Step 3: Create proportional plots
  p1_prop <- ggplot(df_table_prop, aes(x = Batch, y = prop_batch, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = scales::percent(prop_batch, accuracy = 1)),
              position = position_stack(vjust = 0.5), size = 3, color = "black") +
    theme_bw() +
    labs(
      title = "All Clusters (100%)",
      x = batch_col, y = "Percentage", fill = cluster_col
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2_prop <- ggplot(df_table_prop, aes(x = Cluster, y = prop_level, fill = Batch)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = scales::percent(prop_level, accuracy = 1)),
              position = position_stack(vjust = 0.5), size = 3, color = "black") +
    theme_bw() +
    labs(
      title = "Filtered Clusters (100%)",
      x = cluster_col, y = "Percentage", fill = batch_col
    ) +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_prop_combined <- p1_prop + p2_prop + plot_layout(ncol = 2)

  # Step 4: Create absolute count plots
  p1_count <- ggplot(df_table, aes(x = Batch, y = Freq, fill = Cluster)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 3, color = "black") +
    theme_bw() +
    labs(
      title = "Batch",
      x = batch_col, y = "Cell Count", fill = cluster_col
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2_count <- ggplot(df_table, aes(x = Cluster, y = Freq, fill = Batch)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 3, color = "black") +
    theme_bw() +
    labs(
      title = "Clusters",
      x = cluster_col, y = "Cell Count", fill = batch_col
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p_count_combined <- p1_count + p2_count + plot_layout(ncol = 2)

  # Step 5: Setup directory and timestamp
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Step 6: Save both combined plots
  filename_prop <- file.path(output_dir, paste0("cell_prop_plot_", batch_col, "_", file_prefix, "_", timestamp, ".jpg"))
  ggsave(filename_prop, plot = p_prop_combined, width = 20, height = 10, dpi = 300)
  message("Saved proportions to: ", filename_prop)

  filename_count <- file.path(output_dir, paste0("cell_count_plot_", batch_col, "_", file_prefix, "_", timestamp, ".jpg"))
  ggsave(filename_count, plot = p_count_combined, width = 20, height = 10, dpi = 300)
  message("Saved counts to: ", filename_count)

  # Return plots as a list just in case you want to view them in your R session
  return(invisible(list(proportions = p_prop_combined, counts = p_count_combined)))
}
