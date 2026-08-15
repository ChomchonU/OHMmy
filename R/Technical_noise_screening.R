#' Generate Dimensionality Reduction Heatmaps
#'
#' A robust wrapper around Seurat's \code{DimHeatmap} function designed to systematically
#' evaluate sources of heterogeneity across multiple dimensions. It chunks specified
#' principal components (or other dimensional reductions) into distinct windows, generates
#' a balanced heatmap for each window, and safely exports them as high-resolution JPEGs.
#' It features a built-in progress bar and robust error handling (\code{tryCatch}) to
#' prevent graphics device hangs during batch processing on remote compute nodes.
#'
#' @param seurat_obj A Seurat object containing single-cell data and the specified dimensional reduction.
#' @param sample_name Character. The name of the biological sample, used for console messages and prefixing the saved filenames.
#' @param reduction Character. The dimensional reduction to use (e.g., "pca", "ica"). Default is "pca".
#' @param pc_windows A list of numeric vectors. Each vector defines a chunk/window of dimensions to plot together in a single file. Default is \code{list(1:10, 11:20, 21:30)}.
#' @param output_dir Character. Directory path where the generated JPEGs will be saved. Default is "DimHeatmap_plots".
#' @param nfeatures Integer. The number of top features (genes) to display per dimension. Default is 60.
#' @param cells Integer. The number of cells to plot. If \code{cells} is a single number, it plots the most extreme cells from both ends of the spectrum. Default is 500.
#' @param width_in Numeric. The width of the saved JPEG in inches. Default is 15.
#' @param height_in Numeric. The height of the saved JPEG in inches. Default is 45.
#' @param res Numeric. The resolution (dpi) of the saved JPEG. Default is 300.
#'
#' @return Invisibly returns \code{NULL}. The primary purpose of this function is its side effect of saving JPEG plots to the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Evaluate the first 15 principal components in 3 chunks of 5
#' generate_dimheatmaps(
#'   seurat_obj = my_seurat,
#'   sample_name = "Donor_1",
#'   reduction = "pca",
#'   pc_windows = list(1:5, 6:10, 11:15),
#'   output_dir = "QC_Plots/PCA_Heatmaps",
#'   nfeatures = 30,
#'   cells = 500,
#'   width_in = 10,
#'   height_in = 20
#' )
#' }
generate_dimheatmaps <- function(seurat_obj,
                                 sample_name,
                                 reduction = "pca",
                                 pc_windows = list(1:10, 11:20, 21:30),
                                 output_dir = "DimHeatmap_plots",
                                 nfeatures = 60,
                                 cells = 500,
                                 width_in = 15,
                                 height_in = 45,
                                 res = 300) {
  # Load progress bar package
  if (!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  library(progress)

  # Create output directory if missing
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  failed_heatmap <- character(0)

  message("Generating DimHeatmap for: ", sample_name)

  # Check if reduction exists
  if (!(reduction %in% names(seurat_obj@reductions))) {
    message("  Reduction '", reduction, "' not found in ", sample_name)
    return(invisible(NULL))
  }

  # Initialize progress bar
  pb <- progress::progress_bar$new(
    total = length(pc_windows),
    format = "  [:bar] :percent - PC range :current/:total"
  )

  # Loop through PC windows
  for (dims_range in pc_windows) {
    pb$tick()
    range_label <- paste0(min(dims_range), "-", max(dims_range))
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    file_path <- file.path(
      output_dir,
      paste0(sample_name, "_DimHeatmap_", reduction, "_PC", range_label, "_", timestamp, ".jpg")
    )

    tryCatch({
      jpeg(file_path, width = width_in, height = height_in, units = "in", res = res)
      DimHeatmap(seurat_obj, dims = dims_range, nfeatures = nfeatures, cells = cells,
                 balanced = TRUE, reduction = reduction)
      dev.off()
      message("\nSaved heatmap: ", file_path)
    }, error = function(e) {
      dev.off()  # In case jpeg device is open
      message("\nFailed for ", sample_name, " (PC ", range_label, ") - ", e$message)
      failed_heatmap <<- c(failed_heatmap, paste0(sample_name, "_PC", range_label))
    })
  }

  # Summary
  if (length(failed_heatmap) > 0) {
    message("\nThe following samples or ranges failed:")
    message(paste(failed_heatmap, collapse = ", "))
  } else {
    message("\nAll DimHeatmaps generated successfully.")
  }
}

# ---------------------------------------------------

#' Generate and Save Dimensional Reduction Loading Plots
#'
#' A robust wrapper around Seurat's \code{VizDimLoadings} function designed for batch
#' processing in HPC environments. It chunks specified principal components (or other
#' dimensional reductions) into windows, extracts the top features (genes) driving each
#' dimension, and stitches the individual plots together using \code{patchwork}. The
#' combined layouts are safely exported as high-resolution JPEGs. It features a built-in
#' progress bar and strict error handling (\code{tryCatch}) to prevent graphics device
#' hangs during automated runs.
#'
#' @param seurat_obj A Seurat object containing single-cell data and the calculated dimensional reduction.
#' @param sample_name Character. The identifier for the biological sample, used for console messaging and prefixing the saved filenames.
#' @param reduction Character. The dimensional reduction to extract feature loadings from (e.g., "pca", "ica"). Default is "pca".
#' @param pc_windows A list of numeric vectors. Each vector defines a chunk/window of dimensions to plot together in a single stitched file. Default is \code{list(1:10, 11:20, 21:30)}.
#' @param output_dir Character. Directory path where the generated JPEGs will be saved. Default is "VizDimLoadings_plots".
#' @param nfeatures Integer. The number of top and bottom features (genes) with the highest absolute loadings to display per dimension. Default is 60.
#' @param width_in Numeric. The width of the saved JPEG in inches. Default is 20.
#' @param height_in Numeric. The height of the saved JPEG in inches. Default is 45.
#' @param res Numeric. The resolution (dpi) of the saved JPEG. Default is 300.
#' @param ncol Integer. The number of columns to use when wrapping the individual dimension plots together via \code{patchwork}. Default is 5.
#'
#' @return A character vector containing the names of any samples or PC ranges that failed to plot (e.g., due to missing reductions or graphics errors). If all ranges succeed, returns an empty character vector.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Evaluate the feature loadings for the first 12 PCs, grouped in chunks of 4
#' failed_plots <- plot_vizdimloadings(
#'   seurat_obj = my_seurat,
#'   sample_name = "Donor_1",
#'   reduction = "pca",
#'   pc_windows = list(1:4, 5:8, 9:12),
#'   output_dir = "QC_Plots/PCA_Loadings",
#'   nfeatures = 30,
#'   width_in = 16,
#'   height_in = 12,
#'   ncol = 2
#' )
#'
#' # Check if any plots failed during the run
#' if (length(failed_plots) > 0) {
#'   print("Some plots failed:")
#'   print(failed_plots)
#' }
#' }
plot_vizdimloadings <- function(seurat_obj,
                                sample_name,
                                reduction = "pca",
                                pc_windows = list(1:10, 11:20, 21:30),
                                output_dir = "VizDimLoadings_plots",
                                nfeatures = 60,
                                width_in = 20,
                                height_in = 45,
                                res = 300,
                                ncol = 5) {
  # Load progress bar
  if (!requireNamespace("progress", quietly = TRUE)) {
    install.packages("progress")
  }
  library(progress)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  failed_samples <- character(0)

  message("Plotting VizDimLoadings for: ", sample_name)

  if (!(reduction %in% names(seurat_obj@reductions))) {
    message(" Reduction '", reduction, "' not found in sample: ", sample_name)
    failed_samples <- c(failed_samples, sample_name)
    return(failed_samples)
  }

  # Initialize progress bar
  pb <- progress_bar$new(
    total = length(pc_windows),
    format = "  [:bar] :percent - Window :current/:total"
  )

  for (dims_range in pc_windows) {
    pb$tick()
    range_label <- paste0(min(dims_range), "-", max(dims_range))
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    file_path <- file.path(
      output_dir,
      paste0(sample_name, "_VizDimLoadings_PC", range_label, "_", timestamp, ".jpg")
    )

    tryCatch({
      plots <- lapply(dims_range, function(pc) {
        VizDimLoadings(seurat_obj, dims = pc, reduction = reduction, balanced = TRUE, nfeatures = nfeatures) +
          ggtitle(paste("PC", pc))
      })

      combined_plot <- patchwork::wrap_plots(plots, ncol = ncol)

      jpeg(file_path, width = width_in, height = height_in, units = "in", res = res)
      print(combined_plot)
      dev.off()

      message("\nSaved VizDimLoadings plot: ", file_path)

    }, error = function(e) {
      message("\nFailed plotting for ", sample_name, " PC range ", range_label, ": ", e$message)
      failed_samples <<- c(failed_samples, paste0(sample_name, "_PC", range_label))
      if (dev.cur() != 1) dev.off()
    })
  }

  return(failed_samples)
}

# -------------------------------------------------

#' Evaluate Technical Gene Contributions to Principal Components
#'
#' Extracts PCA loadings from a Seurat object and calculates the weighted percentage
#' of technical or nuisance genes (e.g., mitochondrial, ribosomal, immunoglobulins,
#' or specific lncRNAs) driving each principal component. It evaluates the top positive
#' and negative feature loadings independently. The output is a split bar chart that
#' visually highlights which PCs are overwhelmed by technical noise, helping determine
#' which components to exclude or which variables require regression (\code{vars.to.regress}).
#'
#' @note This function explicitly looks for a dimensionality reduction named \code{"pca.log"}.
#' Ensure your Seurat object's PCA was saved under this name, or modify the function
#' to accept a custom reduction name.
#'
#' @param seurat_obj A Seurat object containing single-cell data. Must have a reduction named "pca.log".
#' @param sample_name Character. The identifier for the biological sample, used for plot titles and file naming.
#' @param output_dir Character. Directory path where the generated JPEG will be saved. Default is "Output_R/Find_vartoregress".
#' @param technical_keywords Character vector. Regular expressions defining the "technical" genes to track. Default includes prefixes for mitochondrial ("^MT-"), ribosomal ("^RPL", "^RPS"), immunoglobulin (\code{"^IG[HKL]"}), and common lncRNAs ("MALAT1", "NEAT1", "XIST").
#' @param max_pcs Integer. The maximum number of principal components to evaluate. Default is 40.
#' @param n_top_genes Integer. The total number of top-loading genes to evaluate per PC. This is split evenly, meaning a value of 500 evaluates the top 250 positive and top 250 negative features. Default is 500.
#' @param cutoff Numeric. The threshold percentage used to draw a dashed horizontal warning line on the plot. Default is 15.
#'
#' @return A character vector containing the names of any samples that failed to process (e.g., due to missing reductions). Returns an empty character vector (\code{character(0)}) if successful.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Evaluate how much mitochondrial and ribosomal genes are driving your PCs
#' failed_qc <- plot_technical_contribution(
#'   seurat_obj = my_seurat,
#'   sample_name = "Viral_Infection_Cohort",
#'   technical_keywords = c("^MT-", "^RPL", "^RPS"),
#'   max_pcs = 30,
#'   n_top_genes = 400,
#'   cutoff = 10
#' )
#'
#' if (length(failed_qc) == 0) {
#'   print("Technical contribution plot generated successfully!")
#' }
#' }
plot_technical_contribution <- function(seurat_obj,
                                        sample_name,
                                        output_dir = "Output_R/Find_vartoregress",
                                        technical_keywords = c("^MT-", "^RPL", "^RPS", "^IG[HKL]", "MALAT1", "NEAT1", "XIST"),
                                        max_pcs = 40,
                                        n_top_genes = 500,
                                        cutoff = 15) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  failed_samples <- character(0)
  message("Processing sample: ", sample_name)

  # Check if PCA reduction exists
  if (!("pca.log" %in% names(seurat_obj@reductions))) {
    message(" PCA reduction 'pca.log' not found for ", sample_name)
    return(c(sample_name))
  }

  # Load PCA loadings
  pca_loadings <- tryCatch({
    Loadings(seurat_obj[["pca.log"]])
  }, error = function(e) {
    message("Failed to load PCA loadings for ", sample_name, ": ", e$message)
    return(NULL)
  })

  if (is.null(pca_loadings)) return(c(sample_name))

  results <- list()
  num_pcs <- min(max_pcs, ncol(pca_loadings))

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = num_pcs, style = 3)

  for (i in 1:num_pcs) {
    pc_name <- colnames(pca_loadings)[i]
    pc_load <- pca_loadings[, i]

    pos_genes <- names(sort(pc_load, decreasing = TRUE)[1:(n_top_genes / 2)])
    neg_genes <- names(sort(pc_load, decreasing = FALSE)[1:(n_top_genes / 2)])

    for (direction in c("Positive", "Negative")) {
      genes <- if (direction == "Positive") pos_genes else neg_genes
      vals <- abs(pc_load[genes])

      tech_genes <- genes[grepl(paste(technical_keywords, collapse = "|"), genes, ignore.case = TRUE)]
      tech_vals <- abs(pc_load[tech_genes])

      percent_tech <- if (sum(vals) == 0) 0 else 100 * sum(tech_vals) / sum(vals)

      results[[length(results) + 1]] <- data.frame(
        PC = pc_name,
        Direction = direction,
        PercentTechnical = percent_tech
      )
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)  # Close progress bar

  plot_df <- do.call(rbind, results)
  plot_df$PC <- factor(plot_df$PC, levels = paste0("PC_", 1:num_pcs))
  plot_df$Direction <- factor(plot_df$Direction, levels = c("Positive", "Negative"))

  p <- ggplot(plot_df, aes(x = PC, y = PercentTechnical, fill = Direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "orange")) +
    labs(
      title = paste0("Technical Gene Contribution (Top ", n_top_genes / 2, " Pos & Neg Loadings) - ", sample_name),
      y = "Weighted % Technical Contribution",
      x = "PC",
      fill = "Loading Direction"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- file.path(output_dir, paste0("tech_contrib_posneg_split_", sample_name, "_", timestamp, ".jpg"))

  tryCatch({
    ggsave(filename, plot = p, width = 15, height = 10, dpi = 300)
    message("Saved plot: ", filename)
    return(character(0))  # no failures
  }, error = function(e) {
    message("Failed to save plot for ", sample_name, ": ", e$message)
    return(c(sample_name))
  })
}

# --------------------------------------------

#' Evaluate Depth-Dependent Technical Contributions to Principal Components
#'
#' An advanced quality control tool that tracks how the influence of technical/nuisance
#' genes changes at varying depths of principal component loadings. Instead of checking
#' a single fixed number of top genes, it iterates across a sequence of depths (e.g.,
#' top 20, 40, ..., 500), calculating the weighted percentage of technical genes at each
#' step for both positive and negative directions. The output is a massive faceted
#' bar chart (one facet per PC) with color-coded flags indicating when a loading
#' direction exceeds the defined technical cutoff threshold.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param sample_name Character. The identifier for the biological sample, used for plot titles and file naming.
#' @param reduction Character. The dimensional reduction to extract feature loadings from. Default is "pca.log".
#' @param technical_keywords Character vector. Regular expressions defining the "technical" genes to track. Default includes mitochondrial ("^MT-"), ribosomal ("^RPL", "^RPS"), immunoglobulin (\code{"^IG[HKL]"}), and specific lncRNAs ("MALAT1", "NEAT1", "XIST").
#' @param gene_depths Numeric vector. A sequence defining the varying numbers of top/bottom genes to evaluate sequentially. Default is \code{seq(0, 500, by = 20)}.
#' @param max_pcs Integer. The maximum number of principal components to evaluate and facet. Default is 40.
#' @param cutoff Numeric. The threshold percentage used to draw a warning line and trigger the color-coded flag (Above/Below). Default is 15.
#' @param output_dir Character. Directory path where the generated JPEG will be saved. Default is "Output_R/Find_vartoregress".
#' @param plot_width Numeric. The width of the saved JPEG in inches. Because this plot contains many facets, the default is large (30).
#' @param plot_height Numeric. The height of the saved JPEG in inches. Default is 15.
#' @param dpi Numeric. The resolution (dpi) of the saved JPEG. Default is 300.
#'
#' @return Invisibly returns a character string containing the exact file path where the plot was saved.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Track technical noise saturation across the first 30 PCs at intervals of 25 genes
#' plot_stacked_technical_contribution(
#'   seurat_obj = my_seurat,
#'   sample_name = "T_Cell_Subset",
#'   reduction = "pca",
#'   gene_depths = seq(0, 400, by = 25),
#'   max_pcs = 30,
#'   cutoff = 10,
#'   plot_width = 24,
#'   plot_height = 12
#' )
#' }
plot_stacked_technical_contribution <- function(
    seurat_obj,
    sample_name,
    reduction = "pca.log",
    technical_keywords = c("^MT-", "^RPL", "^RPS", "^IG[HKL]", "MALAT1", "NEAT1", "XIST"),
    gene_depths = seq(0, 500, by = 20),
    max_pcs = 40,
    cutoff = 15,
    output_dir = "Output_R/Find_vartoregress",
    plot_width = 30,
    plot_height = 15,
    dpi = 300
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  pca_loadings <- tryCatch({
    Loadings(seurat_obj[[reduction]])
  }, error = function(e) {
    stop("Failed to load PCA loadings: ", e$message)
  })

  is_technical_gene <- function(genes) {
    grepl(paste(technical_keywords, collapse = "|"), genes, ignore.case = TRUE)
  }

  results_list <- list()
  num_pcs <- min(max_pcs, ncol(pca_loadings))

  total_steps <- length(gene_depths[gene_depths > 0]) * num_pcs  # count only non-zero gene_depths
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0

  for (n_top in gene_depths) {
    half_n <- floor(n_top / 2)
    if (half_n == 0) next

    for (i in 1:num_pcs) {
      pc_name <- colnames(pca_loadings)[i]
      pc_load <- pca_loadings[, i]

      pos_genes <- names(sort(pc_load, decreasing = TRUE)[1:half_n])
      pos_vals <- abs(pc_load[pos_genes])
      pos_tech_genes <- pos_genes[is_technical_gene(pos_genes)]
      pos_tech_vals <- abs(pc_load[pos_tech_genes])
      pos_pct <- if (sum(pos_vals) == 0) 0 else 100 * sum(pos_tech_vals) / sum(pos_vals)

      neg_genes <- names(sort(pc_load, decreasing = FALSE)[1:half_n])
      neg_vals <- abs(pc_load[neg_genes])
      neg_tech_genes <- neg_genes[is_technical_gene(neg_genes)]
      neg_tech_vals <- abs(pc_load[neg_tech_genes])
      neg_pct <- if (sum(neg_vals) == 0) 0 else 100 * sum(neg_tech_vals) / sum(neg_vals)

      results_list[[length(results_list) + 1]] <- data.frame(
        PC = pc_name,
        NTopGenes = n_top,
        Direction = "Positive",
        WeightedPercentTechnical = pos_pct
      )
      results_list[[length(results_list) + 1]] <- data.frame(
        PC = pc_name,
        NTopGenes = n_top,
        Direction = "Negative",
        WeightedPercentTechnical = neg_pct
      )

      step <- step + 1
      setTxtProgressBar(pb, step)
    }
  }
  close(pb)

  tech_pc_split_df <- do.call(rbind, results_list)
  tech_pc_split_df$color_flag <- ifelse(
    tech_pc_split_df$WeightedPercentTechnical > cutoff,
    "Above", "Below"
  )

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 to use this function.")
  }
  library(ggplot2)

  p <- ggplot(tech_pc_split_df, aes(x = factor(NTopGenes), y = WeightedPercentTechnical, fill = interaction(Direction, color_flag))) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ PC, scales = "fixed", ncol = 5) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "black", size = 0.5) +
    scale_fill_manual(
      values = c(
        "Positive.Above" = "red2",
        "Positive.Below" = "green2",
        "Negative.Above" = "pink2",
        "Negative.Below" = "blue2"
      ),
      labels = c(
        "Positive.Below" = "Pos <= cutoff",
        "Positive.Above" = "Pos > cutoff",
        "Negative.Below" = "Neg <= cutoff",
        "Negative.Above" = "Neg > cutoff"
      ),
      name = "Direction / Flag"
    ) +
    labs(
      title = paste0("Weighted % Technical Gene Contribution per PC (Positive & Negative) - ", sample_name),
      x = "Top + Bottom Genes per PC",
      y = "Weighted % Technical Contribution"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- file.path(output_dir, paste0("stacked_technical_contrib_posneg_", sample_name, "_", timestamp, ".jpg"))

  tryCatch({
    ggsave(filename, plot = p, width = plot_width, height = plot_height, dpi = dpi)
    message("Plot saved: ", filename)
  }, error = function(e) {
    message("Failed to save plot for ", sample_name, ": ", e$message)
  })

  invisible(filename)
}

# -----------------------------------------------

#' Calculate and Plot Correlation Between Principal Components and Metadata
#'
#' Extracts dimensional reduction embeddings (e.g., PCA cell scores) from a Seurat object
#' and computes the Spearman correlation against a specified list of continuous metadata
#' covariates (such as sequencing depth, mitochondrial percentage, or cell cycle scores).
#' It generates a hierarchically clustered heatmap of the correlation coefficients using
#' \code{pheatmap} and automatically saves it to disk. This is a critical QC step for
#' identifying whether specific principal components are capturing technical artifacts
#' rather than true biological variance.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param sample_name Character. The identifier for the biological sample, used for the plot title and file naming.
#' @param vars_to_test Character vector. The names of the metadata columns (or specific features in the active assay) to correlate against the PCs. Default includes common technical and cell-cycle covariates: \code{c("pct_counts_mt", "nCount_RNA", "percent.ribo", "percent.ig", "nuclear_rna", "S.Score", "G2M.Score", "MALAT1", "NEAT1", "XIST")}.
#' @param reduction Character. The dimensional reduction to extract cell embeddings from. Default is "pca.log".
#' @param n_pcs Integer. The number of principal components to evaluate, starting from PC 1. Default is 40.
#' @param output_dir Character. Directory path where the generated JPEG heatmap will be saved. Default is "Output_R/Find_vartoregress".
#' @param plot_width_in Numeric. The width of the saved JPEG in inches. Default is 15.
#' @param plot_height_in Numeric. The height of the saved JPEG in inches. Default is 7.
#' @param res_dpi Numeric. The resolution (dpi) of the saved JPEG. Default is 300.
#'
#' @return Invisibly returns a character string containing the exact file path where the correlation heatmap was saved.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming your Seurat object has cell cycle scores (S.Score, G2M.Score)
#' # and mitochondrial percentages (percent.mt) calculated in the metadata
#'
#' plot_pc_metadata_correlation(
#'   seurat_obj = my_seurat,
#'   sample_name = "PBMC_Condition_A",
#'   vars_to_test = c("nCount_RNA", "nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"),
#'   reduction = "pca",
#'   n_pcs = 30,
#'   plot_width_in = 12,
#'   plot_height_in = 6
#' )
#' }
plot_pc_metadata_correlation <- function(
    seurat_obj,
    sample_name,
    vars_to_test = c("pct_counts_mt", "nCount_RNA", "percent.ribo", "percent.ig", "nuclear_rna", "S.Score", "G2M.Score", "MALAT1", "NEAT1", "XIST"),
    reduction = "pca.log",
    n_pcs = 40,
    output_dir = "Output_R/Find_vartoregress",
    plot_width_in = 15,
    plot_height_in = 7,
    res_dpi = 300
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Extract embeddings for PCs
  pc_embeddings <- Embeddings(seurat_obj, reduction)[, 1:n_pcs]

  # Compute correlation matrix (Spearman)
  cor_matrix <- sapply(vars_to_test, function(varname) {
    apply(pc_embeddings, 2, function(pc) cor(pc, seurat_obj[[varname]][, 1], method = "spearman"))
  })

  cor_matrix <- t(cor_matrix)  # transpose: rows=PCs, cols=variables

  width_px <- plot_width_in * res_dpi
  height_px <- plot_height_in * res_dpi

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- file.path(output_dir, paste0("corr_pc_", sample_name, "_", timestamp, ".jpg"))

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' required but not installed.")
  }

  tryCatch({
    jpeg(filename, width = width_px, height = height_px, res = res_dpi)
    pheatmap::pheatmap(
      cor_matrix,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      display_numbers = TRUE,
      number_format = "%.2f",
      main = paste0("Spearman correlation between PCs and technical covariates - ", sample_name)
    )
    dev.off()
    message("Correlation heatmap saved to: ", filename)
  }, error = function(e) {
    message("Failed to save correlation heatmap for ", sample_name, ": ", e$message)
  })

  invisible(filename)
}
