#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param reduction
#' @param pc_windows
#' @param output_dir
#' @param nfeatures
#' @param cells
#' @param width_in
#' @param height_in
#' @param res
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param reduction
#' @param pc_windows
#' @param output_dir
#' @param nfeatures
#' @param width_in
#' @param height_in
#' @param res
#' @param ncol
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param output_dir
#' @param technical_keywords
#' @param max_pcs
#' @param n_top_genes
#' @param cutoff
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param reduction
#' @param technical_keywords
#' @param gene_depths
#' @param max_pcs
#' @param cutoff
#' @param output_dir
#' @param plot_width
#' @param plot_height
#' @param dpi
#'
#' @returns
#' @export
#'
#' @examples
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
        "Positive.Below" = "Pos ≤ cutoff",
        "Positive.Above" = "Pos > cutoff",
        "Negative.Below" = "Neg ≤ cutoff",
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

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param vars_to_test
#' @param reduction
#' @param n_pcs
#' @param output_dir
#' @param plot_width_in
#' @param plot_height_in
#' @param res_dpi
#'
#' @returns
#' @export
#'
#' @examples
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
