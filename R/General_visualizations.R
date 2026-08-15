#' Generate Paired Feature and Density Plots
#'
#' This function takes a Seurat object and a named list of features (which can
#' include both gene names and metadata columns). For each feature, it generates
#' a side-by-side composite of a standard Seurat `FeaturePlot` and a custom density
#' plot. The resulting interleaved plots are stitched together using `patchwork`
#' and automatically saved to disk.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param cluster_col Character. The name of the metadata column to set as the active identity (\code{Idents}).
#' @param genes_list A named list where each element is a character vector of features (e.g., genes or metadata columns). The names of the list (e.g., "T_cell_markers") are used for plot titles and file naming.
#' @param sample_name Character. The name of the biological sample, used for plot titles and file naming. Default is "Sample".
#' @param reduction Character. The dimensionality reduction to use for plotting (e.g., "umap", "umap.cca"). Default is "umap.cca".
#' @param viridis_palette Character. The viridis color palette to use for the density plots (e.g., "viridis", "magma", "plasma"). Default is "viridis".
#' @param output_dir Character. The directory path where the generated plots will be saved. Default is "Plots_feat_dens".
#' @param save_format Character. The file format for the saved plots. Options are "jpg", "png", or "pdf". Default is "jpg".
#' @param width Numeric. The base width multiplier for the saved plot dimensions. Default is 20.
#' @param dpi Numeric. The resolution of the saved plots. Default is 300.
#' @param add_timestamp Logical. Whether to append the current date and time to the saved filenames. Default is TRUE.
#' @param verbose Logical. Whether to display progress bars and console messages. Default is TRUE.
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{combined_plots}: A named list of the generated \code{patchwork} plot objects, where names correspond to the names of \code{genes_list}.
#'   \item \code{output_dir}: A character string of the path where the plots were saved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define a named list of features containing both genes and metadata
#' marker_list <- list(
#'   T_Cells = c("CD3E", "CD4", "CD8A", "percent.mt"),
#'   NK_Cells = c("NCAM1", "NKG7", "GNLY")
#' )
#'
#' # Run the plotting function
#' plot_results <- plot_combined(
#'   seurat_obj = my_seurat,
#'   cluster_col = "seurat_clusters",
#'   genes_list = marker_list,
#'   sample_name = "Donor_1",
#'   reduction = "umap",
#'   save_format = "png"
#' )
#'
#' # Access a specific plot directly in R without opening the saved file
#' plot_results$combined_plots$T_Cells
#' }
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
      ) + ggtitle(paste(feat, "- Feature"))

      p_dens <- Plot_Density_Custom(
        seurat_obj,
        reduction = reduction,
        features = feat,
        viridis_palette = viridis_palette
      ) + ggtitle(paste(feat, "- Density"))

      list(p_feat, p_dens)
    })

    # Flatten: [feat1, dens1, feat2, dens2, feat3, dens3, ...]
    plot_list_flat <- unlist(interleaved_plots, recursive = FALSE)

    # 6 cols -> 3 features per row; each feature takes 2 slots (feat + dens)
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

#' Generate Violin Plots for QC Metrics and Gene Expression
#'
#' This function generates highly customized violin plots for both quality control
#' (QC) metadata and specific gene expression markers across cell clusters.
#' It enhances standard Seurat `VlnPlot` outputs by adding faded, size-adjusted
#' jittered points to better visualize true single-cell distributions. Plots are
#' combined via `patchwork` and safely saved to disk with an automatic PDF fallback.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param sample_name Character. The name of the biological sample, used for plot titles and file naming. Default is "Sample".
#' @param qc_meta_features Character vector. Metadata columns to plot for quality control (e.g., "nCount_RNA", "nFeature_RNA"). Default is c("nCount_RNA", "nFeature_RNA", "pct_counts_mt").
#' @param gene_features Character vector. Specific genes to visualize across clusters. Default is NULL.
#' @param res_col Character. The metadata column containing cluster assignments or cell identities. Default is "seurat_clusters".
#' @param assay Character. The assay to pull gene expression data from. Default is "RNA".
#' @param output_dir Character. Directory path where the generated plots will be saved. Default is "Plots_violin_qc".
#' @param save_format Character. The file format for the saved plots ("jpg", "png", or "pdf"). Default is "jpg".
#' @param width Numeric. The base width of the saved plot. Default is 20.
#' @param dpi Numeric. The resolution of the saved plots. Default is 300.
#' @param add_timestamp Logical. Whether to append the current date and time to the saved filenames. Default is TRUE.
#'
#' @return Invisibly returns \code{NULL}. The primary purpose of this function is its side effect of saving combined plots to the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot basic QC metrics across clusters
#' plot_violin_qc_single(
#'   seurat_obj = my_seurat,
#'   sample_name = "Donor_1",
#'   res_col = "seurat_clusters"
#' )
#'
#' # Plot both QC metrics AND specific T/NK cell markers
#' plot_violin_qc_single(
#'   seurat_obj = my_seurat,
#'   sample_name = "Donor_1",
#'   qc_meta_features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
#'   gene_features = c("CD3E", "CD8A", "NKG7", "GNLY"),
#'   res_col = "CellType",
#'   save_format = "png"
#' )
#' }
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
               width = 20,  # fixed 2 cols (2  8)
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

#' Generate Dimensionality Reduction Plots by Metadata Factors
#'
#' This function iterates through a specified list of metadata columns (factors)
#' in a Seurat object and generates a dimensionality reduction plot (e.g., UMAP)
#' for each. It automatically handles filename sanitization, saves the plots to
#' a designated directory, and provides console progress updates.
#'
#' @param seurat_obj A Seurat object containing single-cell data and the specified reduction.
#' @param factors Character vector. The names of the metadata columns to group the cells by (e.g., c("seurat_clusters", "Phase")).
#' @param sample_name Character. The name of the biological sample, used for plot titles and file naming. Default is "Sample".
#' @param reduction Character. The dimensionality reduction to visualize (e.g., "umap", "umap.har" for Harmony). Default is "umap.har".
#' @param raster Logical. Whether to rasterize the points (useful for massive datasets to speed up plotting). Default is FALSE.
#' @param label Logical. Whether to label the clusters/groups directly on the plot. Default is TRUE.
#' @param repel Logical. Whether to repel the labels to prevent overlapping. Default is TRUE.
#' @param output_dir Character. Directory path where the generated plots will be saved. Default is "Plots_umap".
#' @param plot_format Character. The file format for the saved plots ("jpg", "png", or "pdf"). Default is "jpg".
#' @param width Numeric. The width of the saved plot. Default is 10.
#' @param height Numeric. The height of the saved plot. Default is 10.
#' @param dpi Numeric. The resolution of the saved plots. Default is 300.
#' @param add_timestamp Logical. Whether to append the current date and time to the saved filenames. Default is TRUE.
#' @param verbose Logical. Whether to print progress messages, display a progress bar, and print the plots to the console. Default is TRUE.
#'
#' @return Invisibly returns a named list of the generated \code{ggplot} objects, where names correspond to the provided \code{factors}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot UMAPs colored by multiple different metadata factors
#' my_plots <- PlotDimByFactors(
#'   seurat_obj = my_seurat,
#'   factors = c("seurat_clusters", "CellType", "Condition"),
#'   sample_name = "Donor_1",
#'   reduction = "umap",
#'   plot_format = "png"
#' )
#'
#' # The function saves the plots to disk, but you can also view one in R:
#' my_plots$CellType
#' }
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
      ggtitle(paste0(sample_name, " - ", fac))

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

# -- Palette helper (add this alongside sanitize()) ----------------------------
#' Generate a Discrete Color Palette for Clustering
#'
#' Creates a maximally distinct color palette tailored for high-dimensional single-cell
#' clustering visualizations. Users can choose from curated categorical palettes
#' (Kelly, Alphabet, Polychrome) or an "auto" mode that combines multiple \code{RColorBrewer}
#' sets. If the requested number of colors (\code{n}) exceeds the available base colors
#' in a given palette, the function automatically interpolates using \code{colorRampPalette}
#' to generate the required amount.
#'
#' @param n Integer. The number of distinct colors required (e.g., the number of Seurat clusters).
#' @param palette Character. The specific discrete palette to use. Options are \code{"auto"} (default), \code{"kelly"}, \code{"alphabet"}, or \code{"polychrome"}.
#'
#' @return A character vector of length \code{n} containing hexadecimal color codes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate a 12-color palette using the default hybrid RColorBrewer sets
#' auto_colors <- make_cluster_palette(n = 12, palette = "auto")
#'
#' # Generate a 20-color palette using the Kelly distinct color list
#' kelly_colors <- make_cluster_palette(n = 20, palette = "kelly")
#'
#' # Pass the colors directly into a Seurat plotting function
#' DimPlot(my_seurat, cols = make_cluster_palette(15, "alphabet"))
#'
#' # Requesting 50 colors will automatically trigger interpolation
#' massive_palette <- make_cluster_palette(n = 50, palette = "polychrome")
#' }
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

#' Calculate and Plot Gene-Pair Correlations Across Clusters
#'
#' Computes the expression correlation (e.g., Pearson or Spearman) between specified
#' pairs of genes within each cell cluster. To mitigate zero-inflation artifacts common
#' in scRNA-seq data, it filters cells based on a defined quantile expression threshold
#' before calculating the correlation. The results are visualized as faceted bar plots
#' and automatically saved to disk.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data.
#' @param cluster_col Character. The name of the metadata column defining cell clusters to group the correlations by.
#' @param gene_pairs A list of character vectors, where each vector contains exactly two gene names to correlate (e.g., \code{list(c("GeneA", "GeneB"))}).
#' @param sample_name Character. The name of the biological sample, used for plot titles and file naming. Default is "Sample".
#' @param cor_method Character. The correlation method to use, passed to \code{cor()}. Options include "pearson", "spearman", or "kendall". Default is "pearson".
#' @param quantile_thresh Numeric. The expression quantile threshold (0 to 1) used to filter out lowly expressing cells before correlation calculations. Default is 0.01.
#' @param fill_palette Character. The discrete color palette to pass to \code{make_cluster_palette()}. Default is "auto".
#' @param output_dir Character. Directory path where the generated plots will be saved. Default is "Plots_gene_pair_cor".
#' @param save_format Character. The file format for the saved plots ("jpg", "png", or "pdf"). Default is "jpg".
#' @param width Numeric. The width of the saved plot. Default is 14.
#' @param height Numeric. The height of the saved plot. Default is 8.
#' @param dpi Numeric. The resolution of the saved plots. Default is 300.
#' @param add_timestamp Logical. Whether to append the current date and time to the saved filenames. Default is TRUE.
#' @param verbose Logical. Whether to print progress messages and display a progress bar. Default is TRUE.
#'
#' @return A list containing three elements:
#' \itemize{
#'   \item \code{plot}: The generated \code{ggplot} object containing the faceted bar charts.
#'   \item \code{data}: A \code{tibble} (data frame) containing the calculated correlation values and cell counts per cluster for each valid gene pair.
#'   \item \code{output_dir}: A character string of the path where the plot was saved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define pairs of genes expected to co-express in specific subsets
#' t_nk_pairs <- list(
#'   c("CD8A", "GZMB"),   # Cytotoxic T cell markers
#'   c("NKG7", "GNLY"),   # NK cell effectors
#'   c("CD4", "IL7R")     # Naive/Memory CD4 T cell markers
#' )
#'
#' # Run the correlation analysis across predefined clusters
#' cor_results <- plot_gene_pair_correlations(
#'   seurat_obj = my_seurat,
#'   cluster_col = "CellType_Subset",
#'   gene_pairs = t_nk_pairs,
#'   sample_name = "Donor_1",
#'   cor_method = "spearman",
#'   quantile_thresh = 0.05,
#'   save_format = "png"
#' )
#'
#' # View the raw correlation dataframe directly
#' head(cor_results$data)
#' }
plot_gene_pair_correlations <- function(seurat_obj,
                                        cluster_col,
                                        gene_pairs,
                                        sample_name = "Sample",
                                        cor_method = "pearson",
                                        quantile_thresh = 0.01,
                                        fill_palette = "auto",   #  was "Set2"
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

  # -- Validation ----------------------------------------------------------------
  if (!is.list(gene_pairs) || !all(lengths(gene_pairs) == 2)) {
    stop("`gene_pairs` must be a list of character vectors, each of length 2.")
  }
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("`cluster_col` '", cluster_col, "' not found in Seurat object metadata.")
  }

  # -- Setup ---------------------------------------------------------------------
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

  # -- Fetch expression + cluster data once --------------------------------------
  fetch_genes <- unique(unlist(valid_pairs))
  exp_data    <- FetchData(seurat_obj, vars = c(fetch_genes, cluster_col))
  colnames(exp_data)[ncol(exp_data)] <- "Cluster"

  # -- Loop over pairs -----------------------------------------------------------
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

  # -- Plot -----------------------------------------------------------------------
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
      title    = paste(sample_name, "- Gene Pair Correlations by Cluster"),
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

  # -- Save -----------------------------------------------------------------------
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

  return(list(          #  must be AFTER ggsave and the verbose message
    plot       = p,
    data       = final_cor_df,
    output_dir = output_dir
  ))
}                       #  closing brace of the function

#------------------------------------------

# -- Private helper: Seurat-v5-safe blend plot ----------------------------------
# Returns a patchwork of 4 panels (gene1 | gene2 | blend | colour-key),
# matching the shape that FeaturePlot(blend=TRUE, combine=TRUE) used to give.
#' Generate a 4-Panel Blended Feature Plot (Seurat v5 Compatible)
#'
#' Replicates and enhances the behavior of Seurat's original \code{FeaturePlot(blend = TRUE, combine = TRUE)}
#' functionality, with built-in compatibility for Seurat v5's new \code{layer} architecture
#' (while maintaining a v4 \code{slot} fallback). It computes normalized, thresholded expression
#' and generates a 1x4 \code{patchwork} grid containing: individual expression of gene 1,
#' individual expression of gene 2, their blended co-expression on the embedding, and a customized
#' 2D color threshold key.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param gene1 Character. The name of the first feature/gene to plot (mapped to the x-axis of the color key).
#' @param gene2 Character. The name of the second feature/gene to plot (mapped to the y-axis of the color key).
#' @param reduction Character. The dimensionality reduction to visualize (e.g., "umap", "pca").
#' @param cols Character vector of length 3. Specifies the colors for the background (double-negative), gene 1, and gene 2. Default is c("lightgrey", "#00ff00", "#ff0000").
#' @param blend_threshold Numeric. A scaling threshold (0 to 1) below which normalized expression values are visually suppressed to 0. Default is 0.1.
#' @param min_cutoff Character or Numeric. The lowest expression cutoff (e.g., "q10" for the 10th quantile). Values below this are set to zero before scaling. Default is "q10".
#' @param gamma Numeric. The exponent used to adjust the brightness and saturation of the blended colors. Default is 1.
#' @param pt_size Numeric. The point size for cells in the scatter plots. Default is 1.
#'
#' @return A \code{patchwork} object containing four aligned \code{ggplot} panels.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Visualize co-expression of cytotoxic markers
#' blend_plot <- .blend_feature_plot_v5(
#'   seurat_obj = my_seurat,
#'   gene1 = "CD8A",
#'   gene2 = "GZMB",
#'   reduction = "umap",
#'   cols = c("lightgrey", "#ff0000", "#0000ff"), # Custom red and blue blend
#'   blend_threshold = 0.15,
#'   pt_size = 0.5
#' )
#'
#' # Display the 4-panel plot
#' print(blend_plot)
#' }
.blend_feature_plot_v5 <- function(seurat_obj, gene1, gene2, reduction,
                                   cols            = c("lightgrey", "#00ff00", "#ff0000"),
                                   blend_threshold = 0.1,
                                   min_cutoff      = "q10",
                                   gamma           = 1,
                                   pt_size         = 1) {

  # -- Coordinates ---------------------------------------------------------------
  emb       <- Embeddings(seurat_obj, reduction = reduction)
  dim_names <- colnames(emb)[1:2]
  df        <- as.data.frame(emb[, 1:2])
  colnames(df) <- c("DIM1", "DIM2")

  # -- Expression (v5 uses layer=, v4 used slot=) --------------------------------
  expr_mat <- tryCatch(
    GetAssayData(seurat_obj, layer = "data"),
    error = function(e) GetAssayData(seurat_obj, slot = "data")
  )

  e1 <- as.numeric(expr_mat[gene1, rownames(df)])
  e2 <- as.numeric(expr_mat[gene2, rownames(df)])

  # -- Apply quantile min-cutoff -------------------------------------------------
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

  # -- Scale each gene to [0, 1] -------------------------------------------------
  scale01 <- function(x) { mx <- max(x, na.rm = TRUE); if (mx == 0) x else x / mx }
  s1 <- scale01(e1)
  s2 <- scale01(e2)

  # -- Remap: values below threshold -> 0, threshold -> 0, max -> 1 ----------------
  remap <- function(x) pmax(0, (x - blend_threshold) / (1 - blend_threshold))
  s1r <- remap(s1)
  s2r <- remap(s2)

  # -- Quantize remapped values to 10 steps -------------------------------------
  steps <- 10L
  snap  <- function(x) ceiling(x * steps) / steps
  s1q <- snap(s1r)
  s2q <- snap(s2r)

  # -- Colour setup --------------------------------------------------------------
  col_neither <- cols[1]
  col_g1      <- cols[2]
  col_g2      <- cols[3]
  rgb_neither <- col2rgb(col_neither) / 255
  rgb_g1      <- col2rgb(col_g1)      / 255
  rgb_g2      <- col2rgb(col_g2)      / 255

  # -- Blend function ------------------------------------------------------------
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

  # -- Populate data frame -------------------------------------------------------
  df$gene1     <- s1r
  df$gene2     <- s2r
  df$gene1q    <- s1q
  df$gene2q    <- s2q

  df$blend_col <- mapply(blend_one, s1q, s2q, SIMPLIFY = TRUE)
  df$expressed <- (s1q > 0) | (s2q > 0)
  df <- df[order(df$expressed, s1q + s2q), ]   # background cells plotted first

  # -- Shared theme --------------------------------------------------------------
  base_thm <- theme_classic(base_size = 10) + theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position   = "none",
    axis.title        = element_text(size = 8)
  )
  ax <- labs(x = dim_names[1], y = dim_names[2])

  # -- Panel 1: gene 1 ----------------------------------------------------------
  p1 <- ggplot(df[order(df$gene1), ], aes(x = DIM1, y = DIM2, colour = gene1)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradient(low = col_neither, high = col_g1, limits = c(0, 1)) +
    labs(title = gene1) + ax + base_thm

  # -- Panel 2: gene 2 ----------------------------------------------------------
  p2 <- ggplot(df[order(df$gene2), ], aes(x = DIM1, y = DIM2, colour = gene2)) +
    geom_point(size = pt_size, stroke = 0) +
    scale_colour_gradient(low = col_neither, high = col_g2, limits = c(0, 1)) +
    labs(title = gene2) + ax + base_thm

  # -- Panel 3: blend ------------------------------------------------------------
  p_bl <- ggplot(df, aes(x = DIM1, y = DIM2)) +
    geom_point(colour = df$blend_col, size = pt_size, stroke = 0) +
    labs(title = paste(gene1, "+", gene2)) + ax + base_thm

  # -- Panel 4: Colour key (Matching full 0-1 scale with threshold mask) ---------
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

# -- Main function (only the FeaturePlot call is changed) -----------------------
#' Generate Combined Blend and Nebulosa Density Plots
#'
#' Iterates over a list of gene pairs to generate a comprehensive 6-panel visualization
#' per pair. For each pair, it creates a 4-panel Seurat v5-compatible blended feature plot
#' (showing individual and co-expression) alongside two \code{Nebulosa} density plots.
#' All gene pairs are stacked vertically into a single massive \code{patchwork} figure
#' and automatically saved to disk with dynamically calculated dimensions.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param cluster_col Character. The name of the metadata column to set as the active identity (\code{Idents}).
#' @param gene_pairs A list of character vectors, where each vector contains exactly two gene names to visualize (e.g., \code{list(c("GeneA", "GeneB"))}).
#' @param sample_name Character. The name of the biological sample, used for plot titles and file naming. Default is "Sample".
#' @param reduction Character. The dimensionality reduction to visualize (e.g., "umap.cca"). Default is "umap.cca".
#' @param blend_cols Character vector of length 3. Colors for the background, gene 1, and gene 2 in the blend plot. Default is \code{c("lightgrey", "#00ff00", "#ff0000")}.
#' @param blend_threshold Numeric. The scaling threshold (0 to 1) for the blended feature plot. Default is 0.1.
#' @param min_cutoff Character or Numeric. The lowest expression cutoff before scaling (e.g., "q10"). Default is "q10".
#' @param output_dir Character. Directory path where the generated plots will be saved. Default is "Plots_blend_nebulosa".
#' @param save_format Character. The file format for the saved plots ("jpg", "png", or "pdf"). Default is "jpg".
#' @param width Numeric. The base width of the saved plot. Note: Default is 36 inches to accommodate 6 horizontal panels.
#' @param dpi Numeric. The resolution of the saved plots. Default is 300.
#' @param add_timestamp Logical. Whether to append the current date and time to the saved filenames. Default is TRUE.
#' @param verbose Logical. Whether to print progress messages and display a progress bar. Default is TRUE.
#'
#' @return A list containing four elements:
#' \itemize{
#'   \item \code{plot}: The final combined \code{patchwork} object containing all stacked rows.
#'   \item \code{rows}: A list of the individual \code{patchwork} row plots for each gene pair.
#'   \item \code{output_dir}: A character string of the path where the plot was saved.
#'   \item \code{dimensions}: A list containing the final \code{width} and calculated \code{height} of the saved plot.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define pairs of genes to evaluate for co-expression
#' target_pairs <- list(
#'   c("CD8A", "GZMB"),
#'   c("NKG7", "GNLY")
#' )
#'
#' # Generate the master plot
#' master_blend <- plot_blend_nebulosa(
#'   seurat_obj = my_seurat,
#'   cluster_col = "seurat_clusters",
#'   gene_pairs = target_pairs,
#'   sample_name = "Donor_1",
#'   reduction = "umap",
#'   save_format = "png"
#' )
#'
#' # Extract just the first row (CD8A vs GZMB) if you want to view it directly in R
#' master_blend$rows[[1]]
#' }
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

  # -- Validation ----------------------------------------------------------------
  if (!is.list(gene_pairs) || !all(lengths(gene_pairs) == 2))
    stop("`gene_pairs` must be a list of character vectors, each of length 2.")
  if (!cluster_col %in% colnames(seurat_obj@meta.data))
    stop("`cluster_col` '", cluster_col, "' not found in Seurat object metadata.")

  # -- Setup ---------------------------------------------------------------------
  Idents(seurat_obj) <- cluster_col
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), "%Y-%m-%d_%H-%M-%S") else NULL

  genes_missing <- setdiff(unique(unlist(gene_pairs)), rownames(seurat_obj))
  if (length(genes_missing) > 0)
    warning("  Genes not found and will be skipped: ",
            paste(genes_missing, collapse = ", "))

  valid_pairs <- Filter(function(p) all(p %in% rownames(seurat_obj)), gene_pairs)
  if (length(valid_pairs) == 0) stop("No valid gene pairs after filtering.")

  # -- Layout math ---------------------------------------------------------------
  n_panels_per_row <- 6
  panel_width      <- width / n_panels_per_row
  panel_height     <- panel_width * 0.9
  total_height     <- panel_height * length(valid_pairs)

  if (verbose) {
    message(" Plotting all pairs in one figure for sample: ", sample_name)
    message("   Layout   : ", length(valid_pairs), " rows  ", n_panels_per_row, " panels")
    message("   Width    : ", width, "in  |  Total height: ", round(total_height, 2), "in")
    pb <- txtProgressBar(min = 0, max = length(valid_pairs), style = 3)
  }

  all_rows <- list()

  for (i in seq_along(valid_pairs)) {
    gene1     <- valid_pairs[[i]][1]
    gene2     <- valid_pairs[[i]][2]
    pair_name <- paste(gene1, "vs", gene2)

    # -- CHANGED: use v5-safe manual blend instead of FeaturePlot -------------
    p_blend <- .blend_feature_plot_v5(
      seurat_obj      = seurat_obj,
      gene1           = gene1,
      gene2           = gene2,
      reduction       = reduction,
      cols            = blend_cols,
      blend_threshold = blend_threshold,
      min_cutoff      = min_cutoff
    )

    # -- Nebulosa - gene 1 --------------------------------------------------------
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

    # -- Nebulosa - gene 2 --------------------------------------------------------
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

    # -- Assemble row --------------------------------------------------------------
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

  # -- Stack all rows -------------------------------------------------------------
  p_final <- wrap_plots(all_rows, ncol = 1) +
    plot_annotation(
      title    = paste(sample_name, "- Blend + Density"),
      subtitle = paste0("Reduction: ", reduction,
                        "  |  Blend threshold: ", blend_threshold,
                        "  |  Min cutoff: ", min_cutoff),
      theme = theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40")
      )
    )

  # -- Save ----------------------------------------------------------------------
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

#' Plot Cell Counts and Proportions Across Clusters and Batches
#'
#' Generates stacked bar charts to visualize the distribution of cells between two
#' metadata variables (typically cell clusters and experimental batches/conditions).
#' It calculates and plots both absolute cell counts and relative percentages (proportions),
#' assembling them into two separate dual-panel figures using \code{patchwork}.
#' The resulting plots are automatically saved to the specified output directory.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param cluster_col Character. The name of the metadata column representing cell clusters or identity classes. Default is "label_T3_log_rpca".
#' @param batch_col Character. The name of the metadata column representing batches, samples, or conditions (e.g., "exist_GTS"). Default is "exist_GTS".
#' @param output_dir Character. The directory path where the generated JPEGs will be saved. Default is a specific local directory path.
#' @param file_prefix Character. A prefix string to append to the saved filenames to help identify the analysis run. Default is "ind_label".
#'
#' @return Invisibly returns a list containing two \code{patchwork} plot objects:
#' \itemize{
#'   \item \code{proportions}: The combined plot showing relative percentages.
#'   \item \code{counts}: The combined plot showing absolute cell counts.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare CD8 T cell distributions across different disease states or batches
#' distribution_plots <- plot_cluster_distributions(
#'   seurat_obj = my_seurat,
#'   cluster_col = "CD8_Subsets",
#'   batch_col = "Disease_Status",
#'   output_dir = "Plots_Distributions",
#'   file_prefix = "CD8_analysis"
#' )
#'
#' # The function saves the JPEGs automatically, but you can also view them in R:
#' distribution_plots$proportions
#' distribution_plots$counts
#' }
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
