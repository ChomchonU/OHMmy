#' Generate a Bidirectionally Clustered DotPlot with Dendrograms
#'
#' Extracts expression data via Seurat's \code{DotPlot} function, filters out lowly
#' expressed genes based on a minimum percentage threshold, and performs two-way
#' hierarchical clustering on both the features (genes) and the identities (clusters).
#' It assembles a comprehensive \code{patchwork} layout featuring the central dot plot
#' flanked by column and row dendrograms. The function automatically saves the combined
#' plot, standalone dendrograms, and a text file of the final clustered gene order.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param meta_col Character. The metadata column used to group cells (e.g., "seurat_clusters" or "CellType").
#' @param feature_df A data frame mapping genes to categories. Must contain a \code{group} column (to filter by prefix) and a \code{feature} column (containing the gene names).
#' @param scale Logical. Whether to scale the average expression values (z-score) across groups before plotting and clustering. Default is TRUE.
#' @param pct_threshold Numeric. The minimum percentage of cells expressing the gene in at least one cluster required to retain the gene. Default is 15.
#' @param output_dir Character. Directory path where the generated JPEGs and text file will be saved. Default is "Output_R".
#' @param prefix Character. A string used to filter the \code{group} column in \code{feature_df}. If \code{NULL}, it automatically attempts to derive the prefix by removing "label_" from \code{meta_col}. Default is NULL.
#'
#' @return A \code{patchwork} object representing the final combined layout (dot plot + column dendrogram + row dendrogram).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a mock mapping dataframe for subset markers
#' marker_mapping <- data.frame(
#'   group = c("Effector_T", "Effector_NK", "Naive_T"),
#'   feature = c("GZMB", "NKG7", "CCR7")
#' )
#'
#' # Generate the bi-clustered dotplot
#' clustered_plot <- plot_dot_dendro(
#'   seurat_obj = my_seurat,
#'   meta_col = "Subset_Labels",
#'   feature_df = marker_mapping,
#'   scale = TRUE,
#'   pct_threshold = 10,
#'   output_dir = "Results/DotPlots",
#'   prefix = "Effector" # Will selectively plot and cluster the first two genes
#' )
#'
#' # Display the assembled plot in the R viewer
#' print(clustered_plot)
#' }
plot_dot_dendro <- function(
    seurat_obj,
    meta_col,
    feature_df,
    scale = TRUE,
    pct_threshold = 15,
    output_dir = "Output_R",
    prefix = NULL
) {
  # Ensure meta_col exists
  if (!meta_col %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column '", meta_col, "' does not exist in the Seurat object.")
  }

  if (is.null(prefix)) {
    prefix <- sub("label_", "", meta_col)
  }

  # Filter features by prefix
  genes_to_plot <- feature_df %>%
    filter(grepl(paste0("^", prefix, "_"), group)) %>%
    pull(feature) %>%
    unique()

  if (length(genes_to_plot) == 0) {
    message(" No genes found for prefix '", prefix, "'. Skipping plot.")
    return(NULL)
  }

  # --- Get DotPlot data ---

  if (scale) {
    dp_data <- DotPlot(seurat_obj, features = genes_to_plot, group.by = meta_col)$data
  } else {
    dp_data <- DotPlot(seurat_obj, features = genes_to_plot, group.by = meta_col, scale = F)$data
  }

  if (scale) {
    mat_avg <- dp_data %>%
      select(id, features.plot, avg.exp.scaled) %>%
      pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) %>%
      column_to_rownames("id") %>%
      as.matrix()
  } else {
    mat_avg <- dp_data %>%
      select(id, features.plot, avg.exp) %>%
      pivot_wider(names_from = features.plot, values_from = avg.exp) %>%
      column_to_rownames("id") %>%
      as.matrix()
  }

  mat_pct <- dp_data %>%
    select(id, features.plot, pct.exp) %>%
    pivot_wider(names_from = features.plot, values_from = pct.exp) %>%
    column_to_rownames("id") %>%
    as.matrix()

  # --- Prepare plot data ---
  plot_df <- expand.grid(
    cluster = rownames(mat_avg),
    feature = colnames(mat_avg)
  ) %>%
    mutate(
      avg_exp = as.vector(mat_avg),
      pct_exp = as.vector(mat_pct),
      cluster = factor(cluster, levels = rownames(mat_avg)),
      feature = factor(feature, levels = colnames(mat_avg))
    )

  # Keep genes expressed above threshold in at least one cluster
  filtered_plot_df <- plot_df %>%
    group_by(feature) %>%
    filter(any(pct_exp > pct_threshold)) %>%
    ungroup()

  if (nrow(filtered_plot_df) == 0) {
    message(" No features pass the pct_threshold filter. Skipping plot.")
    return(NULL)
  }

  # --- Build filtered matrices for clustering ---
  mat_avg_filt <- filtered_plot_df %>%
    select(cluster, feature, avg_exp) %>%
    pivot_wider(names_from = feature, values_from = avg_exp) %>%
    column_to_rownames("cluster") %>%
    as.matrix()

  mat_pct_filt <- filtered_plot_df %>%
    select(cluster, feature, pct_exp) %>%
    pivot_wider(names_from = feature, values_from = pct_exp) %>%
    column_to_rownames("cluster") %>%
    as.matrix()

  # --- Hierarchical clustering on filtered features ---

  if (nrow(mat_avg_filt) < 2 || ncol(mat_avg_filt) < 2) {
    message(" Not enough clusters or genes left after filtering. Skipping plot.")
    return(NULL)
  }

  row_hc <- hclust(dist(mat_avg_filt))
  col_hc <- hclust(dist(t(mat_avg_filt)))

  # Reorder matrices
  mat_avg_filt <- mat_avg_filt[row_hc$order, col_hc$order]
  mat_pct_filt <- mat_pct_filt[row_hc$order, col_hc$order]

  # Reorder filtered_plot_df
  filtered_plot_df <- filtered_plot_df %>%
    mutate(
      cluster = factor(cluster, levels = rownames(mat_avg_filt)),
      feature = factor(feature, levels = colnames(mat_avg_filt))
    )

  # --- Gene & cluster dendrogram ---
  col_dd <- ggdendro::dendro_data(as.dendrogram(col_hc))
  row_dd <- ggdendro::dendro_data(as.dendrogram(row_hc))

  # Dot plot
  if (0.3*length(unique(filtered_plot_df$feature)) > 75) {
    if (scale) {
      dot_plot <- ggplot(filtered_plot_df, aes(x = feature, y = cluster)) +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        scale_color_gradient2(high = "#fde725", mid = "#21918c", low = "#440154") +
        scale_radius(range = c(0.05, 3)) +
        theme_bw(base_size = 6) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank())
    } else {
      dot_plot <- ggplot(filtered_plot_df, aes(x = feature, y = cluster)) +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        scale_color_gradient(low = "gray90", high = "darkblue") +
        scale_radius(range = c(0.05, 3)) +
        theme_bw(base_size = 6) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank())
    }
  } else {
    if (scale) {
      dot_plot <- ggplot(filtered_plot_df, aes(x = feature, y = cluster)) +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        scale_color_gradient2(high = "#fde725", mid = "#21918c", low = "#440154") +
        scale_radius(range = c(0.25, 6)) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank())
    } else {
      dot_plot <- ggplot(filtered_plot_df, aes(x = feature, y = cluster)) +
        geom_point(aes(size = pct_exp, color = avg_exp)) +
        scale_color_gradient(low = "gray90", high = "darkblue") +
        scale_radius(range = c(0.25, 6)) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_blank())
    }
  }

  # Row dendrogram (rotated)
  row_dend_rot <- ggplot() +
    geom_segment(data = row_dd$segments, aes(x = y, y = x, xend = yend, yend = xend)) +
    geom_label(data = row_dd$labels, aes(x = y, y = x, label = label),
               hjust = -0.2, size = 3, color = "white", fill = "steelblue") +
    theme_void() +
    scale_y_continuous(expand = expansion(mult = c(0, 0),
                                          add = c(0.5, 0.5)))

  # Combine
  combined_plot <- dot_plot + row_dend_rot + patchwork::plot_layout(widths = c(8, 1))

  # --- Save outputs ---
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe_meta_col <- gsub("[^a-zA-Z0-9_]", "_", meta_col)
  safe_prefix <- gsub("[^a-zA-Z0-9_]", "_", paste0(prefix, "_pct", pct_threshold))

  # Main combined dot + dendrogram plot
  filename <- file.path(output_dir,
                        paste0("DotPlot_Dendro_ggplot_", safe_meta_col, "_", safe_prefix, "_", timestamp, ".jpg"))
  ggsave(filename, combined_plot,
         width = min(75, max(15, 0.3*length(unique(filtered_plot_df$feature)))),
         height = 0.75*nrow(mat_avg_filt),
         dpi = 300, limitsize = FALSE)

  # Save gene dendrogram (column dendrogram only)
  col_dend <- ggplot() +
    geom_segment(data = col_dd$segments,
                 aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = col_dd$labels,
              aes(x = x, y = y, label = label),
              angle = 90, hjust = 1, vjust = 0.5, size = 3) +
    theme_void() +
    scale_x_continuous(expand = expansion(mult = c(0, 0),
                                          add = c(0.5, 2.5))) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02),
                                          add = c(0.75, 0))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())

  col_filename <- file.path(output_dir,
                            paste0("GeneDendrogram_", safe_meta_col, "_", safe_prefix, "_", timestamp, ".jpg"))
  ggsave(col_filename, col_dend,
         width = min(75, max(15, 0.3*length(unique(filtered_plot_df$feature)))),
         height = 20,
         dpi = 300, limitsize = FALSE)

  # ADDED: Save row dendrogram (cluster dendrogram only)
  row_dend_filename <- file.path(output_dir,
                                 paste0("ClusterDendrogram_", safe_meta_col, "_", safe_prefix, "_", timestamp, ".jpg"))
  ggsave(row_dend_filename, row_dend_rot,
         width = 10,
         height = 0.75*nrow(mat_avg_filt),
         dpi = 300, limitsize = FALSE)

  # Save ordered gene list as text file
  ordered_genes <- colnames(mat_avg_filt)
  txt_filename <- file.path(output_dir,
                            paste0("OrderedGenes_", safe_meta_col, "_", safe_prefix, "_", timestamp, ".txt"))
  writeLines(ordered_genes, txt_filename)

  # Empty placeholder for alignment
  empty_plot <- ggplot() + theme_void()

  # Bottom: dot plot + row dendrogram
  dot_row <- dot_plot + row_dend_rot + patchwork::plot_layout(widths = c(8, 1))

  # Top: col dendrogram + empty plot to align with row dendrogram
  col_row <- col_dend + empty_plot + patchwork::plot_layout(widths = c(8, 1))

  # Final combined layout: stack top and bottom
  combined_plot <- col_row / dot_row + patchwork::plot_layout(heights = c(4, 10))

  filename <- file.path(output_dir,
                        paste0("DotPlot_RowColDendro_", safe_meta_col, "_", safe_prefix, "_", timestamp, ".jpg"))

  ggsave(filename, combined_plot,
         width = min(75, max(15, 0.3*length(unique(filtered_plot_df$feature)))),
         height = 0.75*nrow(mat_avg_filt) + 5,
         dpi = 300, limitsize = FALSE)

  message("Saved DotPlot with dendrogram: ", filename)
  message("Saved Gene dendrogram: ", col_filename)
  message("Saved Cluster dendrogram: ", row_dend_filename)
  message("Saved ordered genes: ", txt_filename)

  return(combined_plot)
}

# --------------------------------------

#' Generate Split Clustered DotPlots with Row Dendrograms
#'
#' An extension of the bidirectionally clustered DotPlot, specifically designed to handle
#' massive gene lists without overcrowding the x-axis. It performs global hierarchical
#' clustering on both genes and clusters to establish the true biological ordering.
#' If the number of filtered genes exceeds \code{max_genes_per_plot}, it chunks the x-axis
#' into multiple separate plots. The global cluster dendrogram (y-axis) is identically
#' attached to each chunk for easy cross-referencing.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param meta_col Character. The metadata column used to group cells (e.g., "seurat_clusters").
#' @param feature_df A data frame mapping genes to categories. Must contain a \code{group} column (to filter by prefix) and a \code{feature} column (containing the gene names).
#' @param scale Logical. Whether to scale the average expression values (z-score) across groups before plotting. Default is TRUE.
#' @param pct_threshold Numeric. The minimum percentage of cells expressing the gene in at least one cluster required to retain the gene. Default is 15.
#' @param output_dir Character. Directory path where the generated JPEGs and text file will be saved. Default is "Output_R".
#' @param prefix Character. A string used to filter the \code{group} column in \code{feature_df}. If \code{NULL}, it automatically derives the prefix by removing "label_" from \code{meta_col}. Default is NULL.
#' @param max_genes_per_plot Integer. The maximum number of genes to display per plot. If the filtered gene list exceeds this number, the output is split into multiple chunked plots. Default is NULL (no splitting).
#'
#' @return A list of \code{patchwork} objects, where each element represents a specific chunk of the x-axis (dot plot) combined with the global row (cluster) dendrogram.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming marker_mapping is a dataframe with 'group' and 'feature' columns
#'
#' # Generate chunked dotplots, forcing a maximum of 25 genes per plot
#' chunked_plots <- plot_dot_dendro_split(
#'   seurat_obj = my_seurat,
#'   meta_col = "T_Cell_States",
#'   feature_df = marker_mapping,
#'   scale = TRUE,
#'   pct_threshold = 10,
#'   output_dir = "Results/Split_DotPlots",
#'   prefix = "Exhaustion",
#'   max_genes_per_plot = 25
#' )
#'
#' # View the first chunk directly in R
#' chunked_plots[[1]]
#'
#' # View the second chunk
#' chunked_plots[[2]]
#' }
plot_dot_dendro_split <- function(
    seurat_obj,
    meta_col,
    feature_df,
    scale = TRUE,
    pct_threshold = 15,
    output_dir = "Output_R",
    prefix = NULL,
    max_genes_per_plot = NULL
) {
  # Ensure meta_col exists
  if (!meta_col %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column '", meta_col, "' does not exist in the Seurat object.")
  }

  if (is.null(prefix)) {
    prefix <- sub("label_", "", meta_col)
  }

  # Filter features by prefix
  genes_to_plot <- feature_df %>%
    filter(grepl(paste0("^", prefix, "_"), group)) %>%
    pull(feature) %>%
    unique()

  if (length(genes_to_plot) == 0) {
    message(" No genes found for prefix '", prefix, "'. Skipping plot.")
    return(NULL)
  }

  # --- Get DotPlot data ---
  if (scale) {
    dp_data <- DotPlot(seurat_obj, features = genes_to_plot, group.by = meta_col)$data
  } else {
    dp_data <- DotPlot(seurat_obj, features = genes_to_plot, group.by = meta_col, scale = F)$data
  }

  if (scale) {
    mat_avg <- dp_data %>%
      select(id, features.plot, avg.exp.scaled) %>%
      pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) %>%
      column_to_rownames("id") %>%
      as.matrix()
  } else {
    mat_avg <- dp_data %>%
      select(id, features.plot, avg.exp) %>%
      pivot_wider(names_from = features.plot, values_from = avg.exp) %>%
      column_to_rownames("id") %>%
      as.matrix()
  }

  mat_pct <- dp_data %>%
    select(id, features.plot, pct.exp) %>%
    pivot_wider(names_from = features.plot, values_from = pct.exp) %>%
    column_to_rownames("id") %>%
    as.matrix()

  # --- Prepare plot data ---
  plot_df <- expand.grid(
    cluster = rownames(mat_avg),
    feature = colnames(mat_avg)
  ) %>%
    mutate(
      avg_exp = as.vector(mat_avg),
      pct_exp = as.vector(mat_pct),
      cluster = factor(cluster, levels = rownames(mat_avg)),
      feature = factor(feature, levels = colnames(mat_avg))
    )

  # Keep genes expressed above threshold in at least one cluster
  filtered_plot_df <- plot_df %>%
    group_by(feature) %>%
    filter(any(pct_exp > pct_threshold)) %>%
    ungroup()

  if (nrow(filtered_plot_df) == 0) {
    message(" No features pass the pct_threshold filter. Skipping plot.")
    return(NULL)
  }

  # --- Build filtered matrices for clustering ---
  mat_avg_filt <- filtered_plot_df %>%
    select(cluster, feature, avg_exp) %>%
    pivot_wider(names_from = feature, values_from = avg_exp) %>%
    column_to_rownames("cluster") %>%
    as.matrix()

  mat_pct_filt <- filtered_plot_df %>%
    select(cluster, feature, pct_exp) %>%
    pivot_wider(names_from = feature, values_from = pct_exp) %>%
    column_to_rownames("cluster") %>%
    as.matrix()

  # --- Hierarchical clustering ---
  row_hc <- hclust(dist(mat_avg_filt))
  col_hc <- hclust(dist(t(mat_avg_filt)))

  # Reorder matrices
  mat_avg_filt <- mat_avg_filt[row_hc$order, col_hc$order]
  mat_pct_filt <- mat_pct_filt[row_hc$order, col_hc$order]

  # Reorder filtered_plot_df
  filtered_plot_df <- filtered_plot_df %>%
    mutate(
      cluster = factor(cluster, levels = rownames(mat_avg_filt)),
      feature = factor(feature, levels = colnames(mat_avg_filt))
    )

  # --- Dendrogram data ---
  col_dd <- ggdendro::dendro_data(as.dendrogram(col_hc))
  row_dd <- ggdendro::dendro_data(as.dendrogram(row_hc))

  # --- Build row dendrogram ONCE (reuse for all chunks) ---
  row_dend_rot <- ggplot() +
    geom_segment(data = row_dd$segments, aes(x = y, y = x, xend = yend, yend = xend)) +
    geom_label(data = row_dd$labels, aes(x = y, y = x, label = label),
               hjust = -0.2, size = 3, color = "white", fill = "steelblue") +
    theme_void() +
    scale_y_continuous(expand = expansion(mult = c(0, 0),
                                          add = c(0.5, 0.5)))

  # --- Split into chunks if too many genes ---
  all_features <- levels(filtered_plot_df$feature)
  if (!is.null(max_genes_per_plot) && length(all_features) > max_genes_per_plot) {
    feature_chunks <- split(all_features, ceiling(seq_along(all_features) / max_genes_per_plot))
  } else {
    feature_chunks <- list(all_features)
  }

  plots_out <- list()
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe_meta_col <- gsub("[^a-zA-Z0-9_]", "_", meta_col)
  safe_prefix <- gsub("[^a-zA-Z0-9_]", "_", paste0(prefix, "_pct", pct_threshold))

  for (i in seq_along(feature_chunks)) {
    chunk_features <- feature_chunks[[i]]
    chunk_df <- filtered_plot_df %>% filter(feature %in% chunk_features)

    # Dot plot
    dot_plot <- ggplot(chunk_df, aes(x = feature, y = cluster)) +
      geom_point(aes(size = pct_exp, color = avg_exp)) +
      {if (scale) scale_color_gradient2(high = "#fde725", mid = "#21918c", low = "#440154")
        else scale_color_gradient(low = "gray90", high = "darkblue")} +
      scale_radius(range = c(0.25, 6)) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank())

    # Combine with row dendrogram (reuse)
    dot_row <- dot_plot + row_dend_rot + patchwork::plot_layout(widths = c(8, 1))

    filename <- file.path(output_dir,
                          paste0("DotPlot_RowColDendro_", safe_meta_col, "_", safe_prefix,
                                 "_chunk", i, "_", timestamp, ".jpg"))
    ggsave(filename, dot_row,
           width = min(40, max(10, 0.3*length(chunk_features))),
           height = 0.75*nrow(mat_avg_filt),
           dpi = 300, limitsize = FALSE)

    message("Saved chunk ", i, " -> ", filename)
    plots_out[[i]] <- dot_row
  }

  # --- Save gene lists for each chunk ---
  gene_list_file <- file.path(output_dir,
                              paste0("DotPlot_GeneList_", safe_meta_col, "_", safe_prefix,
                                     "_", timestamp, ".txt"))

  sink(gene_list_file)
  for (i in seq_along(feature_chunks)) {
    cat(paste0("### Chunk ", i, " ###\n"))
    cat(paste(feature_chunks[[i]], collapse = "\n"))
    cat("\n\n") # blank line between chunks
  }
  sink()

  message(" Gene list saved -> ", gene_list_file)

  return(plots_out)
}

# -------------------------------------

#' Generate a Multi-Metadata Bidirectionally Clustered DotPlot
#'
#' An advanced extension of the clustered DotPlot function designed to handle and
#' concatenate multiple metadata groupings simultaneously. It iterates over a vector of
#' metadata columns, aggregates the expression data using Seurat's \code{DotPlot} engine,
#' and binds them into a single comprehensive matrix. It then performs 2D hierarchical
#' clustering across all combined groups and genes. The function features a modular
#' saving system, allowing users to independently export the combined plot, individual
#' dendrograms, or the clustered gene list.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param meta_cols Character vector. A list of metadata columns to evaluate and combine (e.g., \code{c("Condition_A", "Condition_B")} or \code{c("res.0.5", "res.1.0")}).
#' @param feature_df A data frame mapping genes to categories. Must contain a \code{group} column and a \code{feature} column.
#' @param scale Logical. Whether to use scaled average expression values (z-scores) for plotting and clustering. Default is TRUE.
#' @param pct_threshold Numeric. The minimum percentage of cells expressing the gene in at least one cluster/group combination required to retain the gene. Default is 15.
#' @param output_dir Character. Directory path where the generated outputs will be saved. Default is "Output_R".
#' @param prefix Character. A string used to filter the \code{group} column in \code{feature_df}. If \code{NULL}, it attempts to derive it from the metadata column names. Default is NULL.
#' @param base_size Numeric. The base font size for the \code{ggplot2} theme. Default is 12.
#' @param range Numeric vector of length 2. The minimum and maximum point sizes for the dots (\code{scale_radius}). Default is \code{c(0.25, 6)}.
#' @param save_options Character vector. Defines which specific outputs to generate and save. Options include "combined", "col_dend", "row_dend", and "gene_list". Default includes all four.
#'
#' @return Invisibly returns a \code{tibble} containing the filtered, long-format data used to generate the plots (useful for custom downstream plotting in \code{ggplot2}).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming marker_mapping contains your 'group' and 'feature' columns
#'
#' # Compare exhaustion markers across three different metadata cluster resolutions
#' plot_data <- plot_dot_dendro_multi(
#'   seurat_obj = my_seurat,
#'   meta_cols = c("SCT_snn_res.0.4", "SCT_snn_res.0.8", "SCT_snn_res.1.2"),
#'   feature_df = marker_mapping,
#'   pct_threshold = 10,
#'   prefix = "Exhaustion",
#'   output_dir = "Results/Multi_Resolution",
#'   save_options = c("combined", "gene_list") # Skip saving standalone dendrograms
#' )
#'
#' # The function invisibly returns the data, which you can inspect:
#' head(plot_data)
#' }
plot_dot_dendro_multi <- function(
    seurat_obj,
    meta_cols,
    feature_df,
    scale = TRUE,
    pct_threshold = 15,
    output_dir = "Output_R",
    prefix = NULL,
    base_size = 12,
    range = c(0.25, 6),
    save_options = c("combined", "col_dend", "row_dend", "gene_list")
) {
  # Load only necessary namespaces to avoid masking
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(ggdendro)
  require(patchwork)
  require(Seurat)

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # 1. Efficient Data Gathering
  dp_all <- lapply(meta_cols, function(mc) {
    if (!mc %in% colnames(seurat_obj@meta.data)) {
      message(" Metadata column '", mc, "' not found. Skipping.")
      return(NULL)
    }

    pfx <- if (!is.null(prefix)) prefix else sub("label_", "", mc)

    genes_to_plot <- feature_df %>%
      filter(grepl(paste0("^", pfx, "_"), group)) %>%
      pull(feature) %>% unique()

    if (length(genes_to_plot) == 0) return(NULL)

    # Extract data once
    DotPlot(seurat_obj, features = genes_to_plot, group.by = mc)$data %>%
      mutate(id = as.character(id), var = mc)
  }) %>% bind_rows()

  if (nrow(dp_all) == 0) return(NULL)

  # 2. Reshape and Filter (Combined Logic)
  # Instead of creating separate matrices first, we create one wide data frame
  # and filter genes based on the pct_threshold before clustering.

  plot_df_wide <- dp_all %>%
    tidyr::unite("cluster_var", id, var, sep = "_", remove = FALSE) %>%
    select(cluster_var, features.plot, avg.exp.scaled, pct.exp)

  # Find features to keep
  keep_features <- plot_df_wide %>%
    group_by(features.plot) %>%
    filter(any(pct.exp > pct_threshold)) %>%
    pull(features.plot) %>% unique()

  if (length(keep_features) == 0) {
    message(" No features pass pct_threshold filter."); return(NULL)
  }

  # 3. Perform Clustering Once
  mat_avg <- plot_df_wide %>%
    filter(features.plot %in% keep_features) %>%
    select(cluster_var, features.plot, avg.exp.scaled) %>%
    pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) %>%
    tibble::column_to_rownames("cluster_var")

  row_hc <- hclust(dist(mat_avg))
  col_hc <- hclust(dist(t(mat_avg)))

  # Prepare Final Long Data for ggplot
  filtered_plot_df <- dp_all %>%
    filter(features.plot %in% keep_features) %>%
    mutate(
      cluster_var = factor(paste(id, var, sep = "_"), levels = rownames(mat_avg)[row_hc$order]),
      feature = factor(features.plot, levels = colnames(mat_avg)[col_hc$order])
    )

  # 4. Reusable Plotting Components
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe_prefix <- gsub("[^a-zA-Z0-9_]", "_", paste0(prefix, "_pct", pct_threshold))

  # Pre-calculate Dendrogram data if needed by any option
  needs_dendro <- any(c("combined", "col_dend", "row_dend") %in% save_options)

  if (needs_dendro) {
    # Row Dendro (Rotated)
    row_dd <- dendro_data(as.dendrogram(row_hc))
    p_row_dend <- ggplot() +
      geom_segment(data = row_dd$segments, aes(x = y, y = x, xend = yend, yend = xend)) +
      geom_label(data = row_dd$labels, aes(x = y, y = x, label = label),
                 hjust = -0.2, size = 3, color = "white", fill = "steelblue") +
      scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
      theme_void()

    # Column Dendro
    col_dd <- dendro_data(as.dendrogram(col_hc))
    p_col_dend <- ggplot() +
      geom_segment(data = col_dd$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_text(data = col_dd$labels, aes(x = x, y = y, label = label),
                angle = 90, hjust = 1, vjust = 0.5, size = base_size/2) +
      scale_x_continuous(expand = expansion(add = c(0.5, 2.5))) +
      scale_y_continuous(expand = expansion(add = c(0.75, 0))) +
      theme_void()
  }

  # 5. Execute Save Operations
  if ("col_dend" %in% save_options) {
    ggsave(file.path(output_dir, paste0("GeneDendrogram_", safe_prefix, "_", timestamp, ".jpg")),
           p_col_dend, width = min(75, max(15, 0.3 * length(keep_features))), height = 20)
  }

  if ("row_dend" %in% save_options) {
    ggsave(file.path(output_dir, paste0("ClusterDendrogram_", safe_prefix, "_", timestamp, ".jpg")),
           p_row_dend, width = 10, height = max(10, 0.4 * nrow(mat_avg)))
  }

  if ("gene_list" %in% save_options) {
    writeLines(levels(filtered_plot_df$feature),
               file.path(output_dir, paste0("OrderedGenes_", safe_prefix, "_", timestamp, ".txt")))
  }

  if ("combined" %in% save_options) {
    dot_plot <- ggplot(filtered_plot_df, aes(x = feature, y = cluster_var)) +
      geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
      scale_color_gradient2(high = "#fde725", mid = "#21918c", low = "#440154") +
      scale_radius(range = range) +
      theme_bw(base_size = base_size) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())

    combined_plot <- (p_col_dend + plot_spacer() + plot_layout(widths = c(8, 1))) /
      (dot_plot + p_row_dend + plot_layout(widths = c(8, 1))) +
      plot_layout(heights = c(4, 10))

    ggsave(file.path(output_dir, paste0("DotPlot_RowColDendro_", safe_prefix, "_", timestamp, ".jpg")),
           combined_plot, width = min(75, max(15, 0.3 * length(keep_features))),
           height = 0.75 * nrow(mat_avg) + 5, limitsize = FALSE)
  }

  return(invisible(filtered_plot_df))
}
