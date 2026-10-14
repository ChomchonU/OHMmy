# Internal: validate a combined marker table and prepare the condition id column
#' @keywords internal
#' @noRd
.prepare_marker_table <- function(df, id_col, fun_name) {
  if (length(id_col) != 1 || !(id_col %in% colnames(df))) {
    stop("The specified id_col '", paste(id_col, collapse = ", "), "' does not exist in the data.", call. = FALSE)
  }
  required <- c("gene", "cluster", "avg_log2FC", "diff")
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop(fun_name, "() needs the columns ", paste(required, collapse = ", "), " and '", id_col,
         "'; missing: ", paste(missing_cols, collapse = ", "), ".", call. = FALSE)
  }
  df$gene <- as.character(df$gene)
  df$y_id <- as.character(df[[id_col]])
  n_na <- sum(is.na(df$y_id))
  if (n_na > 0) {
    message(" Dropping ", n_na, " row(s) with a missing '", id_col, "' value.")
    df <- df[!is.na(df$y_id), , drop = FALSE]
  }
  df
}

#' Generate Split DotPlots and Global Dendrogram for Marker Genes
#'
#' Filters a differential expression results data frame for specific genes (via list or regex),
#' performs hierarchical clustering on the specified identity column based on average log2FC,
#' and generates a global dendrogram. To prevent unreadable, overcrowded x-axes, the target
#' gene-cluster combinations are automatically split into multiple chunked dotplots.
#' All plots are automatically sized and saved directly to the specified output directory.
#'
#' @param df A data frame of differential expression results. Must contain the columns \code{gene}, \code{cluster}, \code{avg_log2FC}, \code{diff}, and the column specified in \code{id_col}. Rows with a missing \code{id_col} value are dropped, and duplicated gene x cluster x condition rows are averaged.
#' @param gene_list Character vector. A specific list of genes to filter and plot. Default is NULL.
#' @param gene_regex Character string. A regular expression to match target genes (e.g., \code{"^TR[AB][VD]"} for T cell receptor genes or "^KLR" for Killer-cell receptors). Default is NULL.
#' @param id_col Character string. The name of the column in \code{df} to map to the y-axis and use for hierarchical clustering. Default is "integration_method".
#' @param output_dir Character. Directory path where the generated PNGs will be saved. Default is "Markers_Plots/Split".
#' @param n_splits Integer. The number of parts to split the horizontal axis (gene-cluster pairs) into. If there are fewer gene-cluster pairs than \code{n_splits}, one plot per pair is made instead. Default is 4.
#' @param plot_title_prefix Character. The base string used for the title of each split plot. Default is "Gene DotPlot".
#' @param width_scale Numeric. A scaling multiplier to dynamically calculate plot width based on the number of x-axis items. Default is 0.3.
#' @param height_scale Numeric. A scaling multiplier to dynamically calculate plot height based on the number of y-axis items. Default is 0.5.
#' @param min_width Numeric. The absolute minimum width (in inches) for the saved plots. Default is 5.
#' @param min_height Numeric. The absolute minimum height (in inches) for the saved plots. Default is 5.
#' @param label_space Numeric. The spatial padding multiplier for the text labels on the global dendrogram. Default is 1.
#'
#' @return Invisibly returns \code{NULL}. Plots are directly written to disk. When no gene matches, a message is shown and nothing is written. With a single condition, the dot plots are saved but no dendrogram is drawn.
#'
#' @details Numeric cluster ids are labelled \code{"gene (C<cluster>)"} and ordered numerically; other cluster labels (e.g. cell type names) are shown as \code{"gene (<cluster>)"} and ordered alphabetically. Genes are ordered by their mean \code{avg_log2FC}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'marker_df' contains combined FindMarkers output across multiple
#' # integration methods or conditions.
#'
#' # Example 1: Plot specific T/NK cell markers, split into 2 plots
#' plot_split_dotplots_by_gene_cluster(
#'   df = marker_df,
#'   gene_list = c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "GZMB", "PRF1"),
#'   id_col = "Condition",
#'   output_dir = "Results/DotPlots",
#'   n_splits = 2,
#'   plot_title_prefix = "Cytotoxic Markers"
#' )
#'
#' # Example 2: Plot all genes starting with "KLR" using Regex
#' plot_split_dotplots_by_gene_cluster(
#'   df = marker_df,
#'   gene_regex = "^KLR",
#'   id_col = "integration_method",
#'   n_splits = 3
#' )
#' }
plot_split_dotplots_by_gene_cluster <- function(df,
                                                gene_list = NULL,
                                                gene_regex = NULL,
                                                id_col = "integration_method",
                                                output_dir = "Markers_Plots/Split",
                                                n_splits = 4,
                                                plot_title_prefix = "Gene DotPlot",
                                                width_scale = 0.3,
                                                height_scale = 0.5,
                                                min_width = 5,
                                                min_height = 5,
                                                label_space = 1) {
  df <- .prepare_marker_table(df, id_col, "plot_split_dotplots_by_gene_cluster")

  # 1. Filter genes by list or regex
  df_filtered <- df %>%
    filter(
      (if (!is.null(gene_list)) gene %in% gene_list else TRUE) &
        (if (!is.null(gene_regex)) str_detect(gene, stringr::regex(gene_regex, ignore_case = TRUE)) else TRUE)
    )

  if (nrow(df_filtered) == 0) {
    message(" No genes matched.")
    return(invisible(NULL))
  }

  # One value per gene x cluster x condition (duplicated rows are averaged)
  n_before <- nrow(df_filtered)
  df_filtered <- df_filtered %>%
    group_by(gene, cluster, y_id) %>%
    summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
              diff = mean(diff, na.rm = TRUE), .groups = "drop")
  if (nrow(df_filtered) < n_before) {
    message(" Averaged ", n_before - nrow(df_filtered),
            " duplicated gene x cluster x condition row(s).")
  }

  # 2. Order genes by mean log2FC, then clusters (numerically when they are numbers)
  gene_order_only <- df_filtered %>%
    group_by(gene) %>%
    summarize(avg = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
    arrange(avg) %>%
    pull(gene)

  cluster_chr <- as.character(df_filtered$cluster)
  numeric_clusters <- !anyNA(suppressWarnings(as.numeric(cluster_chr)))
  cluster_rank <- if (numeric_clusters) as.numeric(cluster_chr) else match(cluster_chr, sort(unique(cluster_chr)))

  df_filtered <- df_filtered %>%
    mutate(
      gene_cluster = if (numeric_clusters) paste0(gene, " (C", cluster_chr, ")") else paste0(gene, " (", cluster_chr, ")"),
      gene_clean   = factor(gene, levels = gene_order_only),
      cluster_rank = cluster_rank
    ) %>%
    arrange(gene_clean, cluster_rank)

  df_filtered$gene_cluster <- factor(df_filtered$gene_cluster, levels = unique(df_filtered$gene_cluster))

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 3. Global dendrogram (needs at least two conditions)
  conditions <- unique(df_filtered$y_id)
  if (length(conditions) < 2) {
    message(" Only one condition ('", conditions, "') present: no dendrogram is drawn.")
    y_order <- conditions
  } else {
    mat_all <- df_filtered %>%
      dplyr::select(y_id, gene_cluster, avg_log2FC) %>%
      pivot_wider(names_from = gene_cluster, values_from = avg_log2FC, values_fill = 0) %>%
      column_to_rownames("y_id") %>%
      as.matrix()

    dendro <- stats::hclust(dist(mat_all), method = "ward.D2")
    dendro_data <- ggdendro::dendro_data(as.dendrogram(dendro), type = "rectangle")
    y_order <- dendro$labels[dendro$order]

    # Save global dendrogram
    dendro_path <- file.path(output_dir, paste0("Dendrogram_ALL_", timestamp, ".png"))

    p_dendro <- ggplot() +
      geom_segment(
        data = ggdendro::segment(dendro_data),
        aes(x = y, xend = yend, y = x, yend = xend)
      ) +
      geom_text(
        data = ggdendro::label(dendro_data),
        aes(x = -label_space * 4, y = x, label = label),
        hjust = 0, size = 2.5
      ) +
      scale_x_continuous(expand = expansion(add = c(label_space, 0.5))) +
      scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
      theme_bw() +
      theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank()
      )

    ggsave(dendro_path, p_dendro, width = 10, height = max(5, height_scale * length(y_order)), dpi = 300)
    message("Global dendrogram saved: ", dendro_path)
  }

  # 4. Split gene-cluster labels into chunks (never more chunks than labels, no empty chunks)
  gene_cluster_levels <- levels(df_filtered$gene_cluster)
  n_chunks <- max(1, min(n_splits, length(gene_cluster_levels)))
  if (n_chunks < n_splits) {
    message(" Only ", length(gene_cluster_levels), " gene x cluster combination(s): using ",
            n_chunks, " plot(s) instead of ", n_splits, ".")
  }
  chunks <- if (n_chunks == 1) {
    list(Part1 = gene_cluster_levels)             # cut() cannot make a single interval
  } else {
    split(gene_cluster_levels,
          cut(seq_along(gene_cluster_levels), breaks = n_chunks, labels = paste0("Part", seq_len(n_chunks))))
  }
  chunks <- chunks[lengths(chunks) > 0]

  # 5. Loop and plot each chunk
  for (part_name in names(chunks)) {
    selected_genes <- chunks[[part_name]]

    df_part <- df_filtered %>%
      filter(as.character(gene_cluster) %in% selected_genes) %>%
      mutate(
        gene_cluster = factor(as.character(gene_cluster), levels = selected_genes),
        y_id = factor(y_id, levels = y_order)
      )

    n_x <- length(unique(df_part$gene_cluster))
    n_y <- length(unique(df_part$y_id))
    width <- max(width_scale * n_x, min_width)
    height <- max(height_scale * n_y, min_height)

    p <- ggplot(df_part, aes(
      x = gene_cluster,
      y = y_id,
      size = avg_log2FC,
      color = diff
    )) +
      geom_point(alpha = 0.9) +
      scale_size_continuous(name = "avg_log2FC") +
      scale_color_viridis_c(name = "diff", option = "C") +
      labs(
        title = paste(plot_title_prefix, "-", part_name),
        x = "Gene (Cluster)",
        y = id_col
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(face = "bold", size = 12)
      )

    # Save dotplot
    plot_path <- file.path(output_dir, paste0("DotPlot_", part_name, "_", timestamp, ".png"))
    ggsave(plot_path, p, width = width, height = height, dpi = 300, limitsize = FALSE)
    message("DotPlot saved: ", plot_path)
  }

  invisible(NULL)
}

# -------------------------------------------

#' Find and Visualize Top Cluster Markers
#'
#' Computes differential gene expression for a Seurat object-either across all clusters
#' (\code{FindAllMarkers}) or as a specific pairwise comparison. It calculates a custom
#' biological relevance score (combining log2 fold change, percentage difference, and
#' adjusted p-value) to strictly filter and rank the top \code{n} markers per cluster.
#' The function automatically generates a Seurat \code{DoHeatmap} and exports multiple
#' CSV tables containing both the filtered and unfiltered marker lists (ready for GSEA/ORA).
#'
#' @param seurat_obj A Seurat object containing single-cell data with active identities (\code{Idents}) set.
#' @param sample_name Character. The name of the biological sample, used for plot titles, filenames, and sub-directory routing. Default is "Sample".
#' @param use_sct Logical. If \code{TRUE}, sets the default assay to \code{"SCT"}, runs \code{\link[Seurat]{PrepSCTFindMarkers}} to ensure model comparability across samples/batches, and performs marker identification on the SCT assay. If the SCT models still have unequal library sizes afterwards (which happens when an SCT object has been subset), the already-corrected counts are used (\code{recorrect_umi = FALSE}) instead of failing. Defaults to \code{FALSE}.
#' @param marker_diff_thresh Numeric. The minimum required absolute difference in the percentage of expressing cells between groups (\code{abs(pct.1 - pct.2)}). Default is 0.1.
#' @param marker_pval_adj Numeric. The maximum adjusted p-value allowed for a gene to be considered a significant marker. Default is 0.05.
#' @param marker_avg_log2FC_thresh Numeric. The minimum absolute log2 fold change required. Default is 0.5.
#' @param top_n Integer. The number of top-scoring marker genes to select per cluster for plotting on the heatmap. Default is 20.
#' @param output_dir_base Character. The base directory path where outputs (plots and CSVs) will be saved. Default is "Plots_heatmap".
#' @param plot_format Character. The file format for the saved heatmap ("jpg", "png", or "pdf"). Default is "jpg".
#' @param width Numeric. The width of the saved heatmap. Default is 10.
#' @param height Numeric. The height of the saved heatmap. Default is 20.
#' @param dpi Numeric. The resolution of the saved heatmap. Default is 300.
#' @param compare Character vector of length 2. If provided, performs a pairwise \code{FindMarkers} test between these two specific identities (e.g., \code{c("Effector_T", "Naive_T")}). If \code{NULL}, runs \code{FindAllMarkers} across all clusters. Default is NULL.
#' @param add_timestamp Logical. Whether to append the current date and time to the saved filenames. Default is TRUE.
#' @param onlyPos Logical. If TRUE, only identifies positive markers (upregulated genes). Passed to the \code{only.pos} argument in Seurat. Default is TRUE.
#'
#' @return A list containing four elements (returned invisibly, with \code{heatmap} and \code{output_dir} set to \code{NULL} and a warning, when Seurat finds no markers at all):
#' \itemize{
#'   \item \code{top_markers}: A \code{tibble} of the highly filtered, top-scoring marker genes used for the heatmap.
#'   \item \code{heatmap}: The \code{ggplot} object of the generated \code{DoHeatmap}.
#'   \item \code{markers}: A data frame of the full, unfiltered differential expression results.
#'   \item \code{output_dir}: A character string of the specific path where the files were saved.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Find top 10 positive markers for ALL clusters in the object
#' all_markers_res <- FindTopMarkersAndHeatmap(
#'   seurat_obj = my_seurat,
#'   sample_name = "PBMC_Global",
#'   top_n = 10,
#'   onlyPos = TRUE
#' )
#'
#' # Example 2: Perform a specific pairwise comparison (e.g., CD8+ Effector vs Naive)
#' # Capturing both up and downregulated genes
#' pairwise_res <- FindTopMarkersAndHeatmap(
#'   seurat_obj = my_seurat,
#'   sample_name = "CD8_Subset",
#'   compare = c("CD8_Effector", "CD8_Naive"),
#'   onlyPos = FALSE,
#'   marker_avg_log2FC_thresh = 0.25,
#'   plot_format = "png"
#' )
#' }
FindTopMarkersAndHeatmap <- function(
    seurat_obj,
    sample_name = "Sample",
    use_sct = FALSE,
    marker_diff_thresh = 0.1,
    marker_pval_adj = 0.05,
    marker_avg_log2FC_thresh = 0.5,
    top_n = 20,
    output_dir_base = "Plots_heatmap",
    plot_format = "jpg",
    width = 10,
    height = 20,
    dpi = 300,
    compare = NULL,
    add_timestamp = TRUE,
    onlyPos = TRUE
) {
  # --- SCT preparation ---
  recorrect_umi <- TRUE
  if (isTRUE(use_sct)) {
    message("[", sample_name, "] Preparing SCT models for marker discovery...")
    DefaultAssay(seurat_obj) <- "SCT"
    seurat_obj <- PrepSCTFindMarkers(seurat_obj)
    active_assay <- "SCT"

    # PrepSCTFindMarkers() skips re-correction whenever every stored median UMI is
    # above the smallest observed one, which is typical after subsetting an SCT
    # object. FindMarkers() would then refuse every comparison, so fall back to the
    # already-corrected counts, as Seurat recommends for subset objects.
    if (!.sct_models_consistent(seurat_obj, "SCT")) {
      message("[", sample_name, "] SCT models still have unequal library sizes (common after ",
              "subsetting an SCT object). Using the existing corrected counts (recorrect_umi = FALSE).")
      recorrect_umi <- FALSE
    }
  } else {
    active_assay <- DefaultAssay(seurat_obj) # Will fallback to RNA or the current default
  }

  message("[", sample_name, "] Finding markers using ", active_assay, " assay...")

  # recorrect_umi is only understood by the SCT method of FindMarkers()
  if (is.null(compare)) {
    markers <- if (active_assay == "SCT") {
      FindAllMarkers(seurat_obj, assay = active_assay, only.pos = onlyPos, logfc.threshold = 0,
                     min.pct = 0, recorrect_umi = recorrect_umi)
    } else {
      FindAllMarkers(seurat_obj, assay = active_assay, only.pos = onlyPos, logfc.threshold = 0,
                     min.pct = 0)
    }
  } else {
    markers <- if (active_assay == "SCT") {
      FindMarkers(seurat_obj, assay = active_assay, ident.1 = compare[1], ident.2 = compare[2],
                  only.pos = onlyPos, logfc.threshold = 0, min.pct = 0, recorrect_umi = recorrect_umi)
    } else {
      FindMarkers(seurat_obj, assay = active_assay, ident.1 = compare[1], ident.2 = compare[2],
                  only.pos = onlyPos, logfc.threshold = 0, min.pct = 0)
    }

    # --- CRITICAL PATCH FOR PAIRWISE COMPARISON ---
    # FindMarkers lacks 'gene' and 'cluster' columns. We must create them.
    markers$gene <- rownames(markers)

    # Assign the cluster name based on the direction of the fold change
    if("avg_log2FC" %in% colnames(markers)) {
      markers$cluster <- ifelse(markers$avg_log2FC > 0, compare[1], compare[2])
    }
  }

  # Seurat returns an empty table (with warnings) when every comparison failed
  required_cols <- c("pct.1", "pct.2", "avg_log2FC", "p_val_adj")
  if (nrow(markers) == 0 || !all(required_cols %in% colnames(markers))) {
    warning("[", sample_name, "] No markers were returned, so no heatmap or tables were written. ",
            "See the warnings from Seurat above for the reason.", call. = FALSE)
    return(invisible(list(top_markers = markers, heatmap = NULL, markers = markers, output_dir = NULL)))
  }

  # --- PATCH FOR NEGATIVE FILTERING ---
  # Calculate absolute differences so genes upregulated in ident.2 (which have negative
  # log2FC and diff) are not accidentally filtered out by the positive thresholds.
  top_markers <- markers %>%
    mutate(
      abs_diff = abs(pct.1 - pct.2),
      abs_log2FC = abs(avg_log2FC),
      score = abs_log2FC * abs_diff / (p_val_adj + 1e-10)
    ) %>%
    filter(
      p_val_adj <= marker_pval_adj,
      abs_diff >= marker_diff_thresh,
      abs_log2FC >= marker_avg_log2FC_thresh
    ) %>%
    group_by(cluster) %>%
    slice_max(order_by = score, n = top_n)

  message(" [", sample_name, "] Plotting heatmap...")

  heatmap <- DoHeatmap(seurat_obj, features = top_markers$gene, assay = active_assay) +
    ggtitle(paste0("Top Markers - ", sample_name))

  # Prepare output path
  sample_path <- gsub("\\|", "/", sample_name)
  parent_path <- dirname(sample_path)
  output_dir <- file.path(output_dir_base, parent_path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  timestamp <- if (isTRUE(add_timestamp)) format(Sys.time(), "%Y-%m-%d_%H-%M-%S") else NULL

  # --- Generate standardized file names ---
  fname_parts <- c(
    sanitize(sample_name),
    "Heatmap",
    "topMarkers"
  )
  if (!is.null(timestamp)) fname_parts <- c(fname_parts, timestamp)
  file_base <- paste(fname_parts, collapse = "_")

  heatmap_file <- file.path(output_dir, paste0(file_base, ".", plot_format))
  markers_file <- file.path(output_dir, paste0(file_base, "_all_markers_",onlyPos,".csv"))
  markers_GSEA_file <- file.path(output_dir, paste0(file_base, "_for_GSEA&ORA_",onlyPos,".csv"))
  top_markers_file <- file.path(output_dir, paste0(file_base, "_top_markers_",onlyPos,".csv"))

  # --- Save outputs ---
  ggsave(heatmap_file, plot = heatmap, width = width, height = height, dpi = dpi)
  write.csv(markers, markers_file, row.names = FALSE)
  write.csv(markers, markers_GSEA_file, row.names = FALSE)
  write.csv(top_markers, top_markers_file, row.names = FALSE)

  message("[", sample_name, "] Heatmap saved to: ", heatmap_file)
  message("[", sample_name, "] Marker tables saved to: ", output_dir)

  return(list(
    top_markers = top_markers,
    heatmap = heatmap,
    markers = markers,
    output_dir = output_dir
  ))
}

# ------------------------------------------------

#' Generate a Gene Marker DotPlot with Aligned Dendrogram
#'
#' Filters a differential expression data frame for target genes (using either a specific
#' vector of genes or a regular expression), calculates average expression metrics, and
#' performs Ward.D2 hierarchical clustering on the specified condition/identity column.
#' It generates a customized dot plot (with clusters labeled directly on the points)
#' and aligns it side-by-side with the hierarchical dendrogram using \code{patchwork}.
#' The individual components and the combined layout are automatically saved to disk.
#'
#' @param df A data frame containing differential expression results. Must contain \code{gene}, \code{cluster}, \code{avg_log2FC}, \code{diff}, and the column specified in \code{id_col}. Rows with a missing \code{id_col} value are dropped.
#' @param genes Character vector or string. Either a vector of gene names, matched exactly (e.g., \code{c("CD4", "CD8A")}), or a single regular expression, matched case-insensitively (e.g., \code{"^TR[AB][VD]"}, \code{"GZM"} or \code{"HLA-D[PQR]"}). A single string is used as a regular expression unless it is itself a gene name in \code{df}.
#' @param id_col Character. The name of the column in \code{df} representing the experimental condition or identity to cluster on the y-axis. Default is "short_id".
#' @param output_dir Character. Directory path where the generated PNGs will be saved. Default is "Markers_Plots".
#' @param plot_title Character. The title displayed at the top of the dot plot. Default is "Marker Genes Plot".
#' @param dendro_side Character. Which side of the dot plot to attach the dendrogram: \code{"right"} (default) or \code{"left"}.
#' @param label_space Numeric. The horizontal padding multiplier for the text labels on the dendrogram to prevent truncation. Default is 1.
#' @param width_scale Numeric. A scaling multiplier to dynamically calculate plot width based on the number of genes. Default is 0.5.
#' @param height_scale Numeric. A scaling multiplier to dynamically calculate plot height based on the number of conditions (\code{id_col}). Default is 0.5.
#' @param min_width Numeric. The absolute minimum width (in inches) for the saved plots. Default is 5.
#' @param min_height Numeric. The absolute minimum height (in inches) for the saved plots. Default is 4.
#'
#' @return Invisibly returns a list containing three \code{patchwork}/\code{ggplot} objects (\code{NULL}, with a message, when no gene matches; with a single condition, only the dot plot is drawn and saved, and \code{dendrogram} and \code{combined} are \code{NULL}):
#' \itemize{
#'   \item \code{dotplot}: The gene expression dot plot.
#'   \item \code{dendrogram}: The hierarchical clustering dendrogram.
#'   \item \code{combined}: The side-by-side aligned layout of both plots.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'marker_results' is a dataframe combining FindMarkers output
#' # across multiple datasets or integration methods.
#'
#' # Example 1: Using a specific list of T cell exhaustion and activation markers
#' plot_gene_markers_with_dendro(
#'   df = marker_results,
#'   genes = c("PDCD1", "HAVCR2", "LAG3", "TOX", "IFNG", "GZMB"),
#'   id_col = "integration_method",
#'   plot_title = "Exhaustion Markers Across Integrations",
#'   dendro_side = "left"
#' )
#'
#' # Example 2: Using a Regex search to pull all Killer-cell immunoglobulin-like receptors
#' plot_gene_markers_with_dendro(
#'   df = marker_results,
#'   genes = "^KIR",
#'   id_col = "Condition",
#'   output_dir = "NK_Cell_Analysis"
#' )
#' }
plot_gene_markers_with_dendro <- function(df, genes,
                                          id_col = "short_id",
                                          output_dir = "Markers_Plots",
                                          plot_title = "Marker Genes Plot",
                                          dendro_side = "right",
                                          label_space = 1,
                                          width_scale = 0.5,
                                          height_scale = 0.5,
                                          min_width = 5,
                                          min_height = 4) {
  dendro_side <- match.arg(tolower(dendro_side), c("right", "left"))
  df <- .prepare_marker_table(df, id_col, "plot_gene_markers_with_dendro")

  # A single string that is not itself a gene name is used as a regular expression
  # (e.g. "^KIR", "GZM", "HLA-D[PQR]"); otherwise genes are matched exactly.
  is_regex <- is.character(genes) && length(genes) == 1 && !(genes %in% df$gene)

  # 1. Filter and prepare
  df_igh <- if (is_regex) {
    df %>% filter(grepl(genes, gene, ignore.case = TRUE))
  } else {
    df %>% filter(gene %in% genes)
  }
  df_igh <- df_igh %>% dplyr::select(gene, cluster, y_id, avg_log2FC, diff)

  if (nrow(df_igh) == 0) {
    message(" No genes matched.")
    return(invisible(NULL))
  }

  # 2. Aggregate per gene and y_id
  df_igh <- df_igh %>%
    group_by(gene, y_id) %>%
    summarise(
      avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
      diff = mean(diff, na.rm = TRUE),
      cluster = paste(sort(unique(cluster)), collapse = ","),
      .groups = "drop"
    )

  # 3. Hierarchical clustering on y_id (needs at least two conditions)
  conditions <- unique(df_igh$y_id)
  has_dendro <- length(conditions) >= 2
  if (has_dendro) {
    mat <- df_igh %>%
      dplyr::select(y_id, gene, avg_log2FC) %>%
      pivot_wider(names_from = gene, values_from = avg_log2FC, values_fill = 0) %>%
      column_to_rownames("y_id") %>%
      as.matrix()

    dendro <- stats::hclust(dist(mat), method = "ward.D2")
    ordered_y <- dendro$labels[dendro$order]
  } else {
    message(" Only one condition ('", conditions, "') present: no dendrogram is drawn.")
    ordered_y <- conditions
  }
  df_igh$y_id <- factor(df_igh$y_id, levels = ordered_y)

  # 4. Order genes by avg expression
  gene_order <- df_igh %>%
    group_by(gene) %>%
    summarise(avg = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
    arrange(avg) %>%
    pull(gene)
  df_igh$gene <- factor(df_igh$gene, levels = gene_order)

  # 5. Dynamic sizing
  n_x <- length(unique(df_igh$gene))
  n_y <- length(unique(df_igh$y_id))
  plot_width <- max(width_scale * n_x, min_width)
  plot_height <- max(height_scale * n_y, min_height)

  # 6. Dotplot
  p_dot <- ggplot(df_igh, aes(
    x = gene,
    y = y_id,
    size = avg_log2FC,
    color = diff
  )) +
    geom_point(alpha = 0.9) +
    geom_text(aes(label = cluster), color = "black", size = 1.5, vjust = -3, hjust = 0.5) +
    scale_size_continuous(name = "avg_log2FC") +
    scale_color_viridis_c(name = "diff", option = "C") +
    labs(
      title = plot_title,
      x = "Gene",
      y = paste("Condition (", id_col, ")", sep = "")
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  # 9. Save the dot plot
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  dot_path <- file.path(output_dir, paste0("DotPlot_", timestamp, ".png"))
  ggsave(dot_path, p_dot, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
  message("Dot plot saved: ", dot_path)

  if (!has_dendro) {
    return(invisible(list(dotplot = p_dot, dendrogram = NULL, combined = NULL)))
  }

  # 7. Dendrogram
  dendro_data <- ggdendro::dendro_data(as.dendrogram(dendro), type = "rectangle")

  p_dendro <- ggplot() +
    geom_segment(
      data = ggdendro::segment(dendro_data),
      aes(x = y, xend = yend, y = x, yend = xend)
    ) +
    geom_text(
      data = ggdendro::label(dendro_data),
      aes(x = -label_space * 4, y = x, label = label),
      hjust = 0,
      size = 2.5
    ) +
    scale_x_continuous(expand = expansion(add = c(label_space, 0.5))) +
    # same span as the discrete y axis of the dot plot, so leaves line up with its rows
    scale_y_continuous(limits = c(1, length(ordered_y)), expand = expansion(add = 0.6)) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  # 8. Combine
  if (dendro_side == "left") {
    combined_plot <- p_dendro + p_dot + patchwork::plot_layout(widths = c(1, 5))
  } else {
    combined_plot <- p_dot + p_dendro + patchwork::plot_layout(widths = c(5, 1))
  }

  # 10. Save the dendrogram and the combined figure
  dendro_path <- file.path(output_dir, paste0("Dendrogram_", timestamp, ".png"))
  combined_path <- file.path(output_dir, paste0("DotPlot_with_Dendrogram_", timestamp, ".png"))

  ggsave(dendro_path, p_dendro, width = 10, height = plot_height, dpi = 300)
  ggsave(combined_path, combined_plot, width = plot_width + 10, height = plot_height, dpi = 300, limitsize = FALSE)

  message("Dendrogram saved: ", dendro_path)
  message("Combined plot saved: ", combined_path)

  invisible(list(dotplot = p_dot, dendrogram = p_dendro, combined = combined_plot))
}

# ---------------------------------------------------------

#' Extract and Bin Average Gene Expression into Terciles
#'
#' Leverages the internal data extraction of Seurat's \code{DotPlot} function to calculate
#' the unscaled average expression and percent expressed for a specified list of genes
#' across cell groups. It then evaluates each gene independently, calculates its expression
#' terciles (33rd and 67th percentiles), and categorizes each cluster's expression into
#' "Low", "Int" (Intermediate), or "High" bins. This is highly useful for converting
#' continuous transcriptomic data into discrete categories for simplified metadata assignment
#' or categorical plotting.
#'
#' @param seurat_obj A Seurat object containing single-cell data.
#' @param gene_list Character vector. A list of specific gene names to extract and bin (e.g., \code{c("CD3E", "CD8A")}).
#' @param group_col Character. The metadata column name defining the cell groups or clusters to aggregate the expression by. Default is "seurat_clusters".
#'
#' @return A data frame (\code{tibble}) containing the following columns:
#' \itemize{
#'   \item \code{Cluster}: The identity class / group.
#'   \item \code{Gene}: The feature name.
#'   \item \code{AvgExpression}: The raw average expression within the group.
#'   \item \code{PctExpress}: The percentage of cells in the group expressing the gene.
#'   \item \code{Expression_Level}: An ordered factor (\code{"Low", "Int", "High"}) representing the binned category.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define a few key functional markers
#' target_genes <- c("CD8A", "GZMB", "PRF1", "IFNG")
#'
#' # Extract and bin the expression across fine-resolution clusters
#' binned_data <- extract_binned_expression(
#'   seurat_obj = my_seurat,
#'   gene_list = target_genes,
#'   group_col = "T_Cell_Subsets"
#' )
#'
#' # View the resulting data frame
#' head(binned_data)
#'
#' # Filter to find which clusters have "High" expression of GZMB
#' subset(binned_data, Gene == "GZMB" & Expression_Level == "High")
#' }
extract_binned_expression <- function(seurat_obj, gene_list, group_col = "seurat_clusters") {

  # 1. Use Seurat's internal DotPlot extraction to safely grab both Average Expression and Percent Expressed
  dp_info <- DotPlot(seurat_obj, features = gene_list, group.by = group_col, scale = F)

  # 2. Extract the matrix into a data frame format
  df <- dp_info$data %>%
    select(Cluster = id, Gene = features.plot, AvgExpression = avg.exp, PctExpress = pct.exp)

  # 3. Categorize into High, Int, and Low based on quantiles PER GENE
  df_categorized <- df %>%
    group_by(Gene) %>%
    mutate(
      lower_bound = quantile(AvgExpression, probs = 0.3333, na.rm = TRUE),
      upper_bound = quantile(AvgExpression, probs = 0.6667, na.rm = TRUE),

      Expression_Level = case_when(
        AvgExpression >= upper_bound ~ "High",
        AvgExpression >= lower_bound & AvgExpression < upper_bound ~ "Int",
        AvgExpression < lower_bound ~ "Low"
      ),
      # Convert to a factor so the legend displays in the correct order
      Expression_Level = factor(Expression_Level, levels = c("Low", "Int", "High"))
    ) %>%
    select(-lower_bound, -upper_bound) %>%
    ungroup()

  return(df_categorized)
}

#------------------------------------------------------------------

#' Generate a Trio of Volcano Plots for Differential Expression
#'
#' Takes a differential expression results data frame (e.g., from DESeq2 or Seurat)
#' and generates a three-panel volcano plot using \code{EnhancedVolcano} and \code{patchwork}.
#' It applies strict significance thresholds (\code{padj < 0.05} and \code{abs(log2FoldChange) >= 1.0})
#' to automatically extract and label:
#' 1) The top 10 significantly upregulated and downregulated genes.
#' 2) A user-provided list of target functional genes (filtered to only show significant ones).
#' 3) A combined view of both the top 20 genes and the significant target genes.
#'
#' @param res A data frame (or DESeq2 \code{DESeqResults} object) containing differential expression results. Must have gene symbols as row names and the columns \code{padj} (adjusted p-value) and \code{log2FoldChange}. Requires the \pkg{EnhancedVolcano} package.
#' @param main_title Character. The overarching title displayed at the very top of the combined plot.
#' @param target_genes Character vector. A specific list of genes of interest to highlight (e.g., functional pathway markers). Only genes in this list that meet the strict significance thresholds will be labeled.
#'
#' @return A \code{patchwork} object containing three horizontally aligned \code{ggplot}/volcano plot panels.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'de_results' is a data frame from DESeq2 or Seurat FindMarkers
#' # (Make sure your Seurat output column 'avg_log2FC' is renamed to 'log2FoldChange'
#' # and 'p_val_adj' to 'padj' before passing it to this function!)
#'
#' # Define a list of functional genes you care about for your specific subset
#' nk_functional_genes <- c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG", "KLRK1")
#'
#' # Generate the 3-panel plot
#' volcano_trio <- generate_volcano_trio(
#'   res = de_results,
#'   main_title = "Differential Expression: NK Cells (Infected vs Control)",
#'   target_genes = nk_functional_genes
#' )
#'
#' # Display the plot
#' print(volcano_trio)
#'
#' # Save the wide format plot
#' ggsave("Volcano_Trio.png", plot = volcano_trio, width = 24, height = 8, dpi = 300)
#' }
generate_volcano_trio <- function(res, main_title, target_genes) {
  .check_suggested("EnhancedVolcano", "generate_volcano_trio()")
  res <- as.data.frame(res)

  # 1. STRICT FILTER: Keep only genes significant in BOTH padj AND Fold Change
  sig_strict <- res[!is.na(res$padj) &
                      res$padj < 0.05 &
                      abs(res$log2FoldChange) >= 1.0, ]

  # Order by adjusted p-value
  sig_strict <- sig_strict[order(sig_strict$padj), ]

  # Extract Top 10 Up and Down from this strictly significant pool
  top_10_up <- head(rownames(sig_strict[sig_strict$log2FoldChange > 0, ]), 10)
  top_10_down <- head(rownames(sig_strict[sig_strict$log2FoldChange < 0, ]), 10)
  top_20_genes <- c(top_10_up, top_10_down)

  # 2. FILTER TARGET GENES: Only keep target genes that exist in the sig_strict list
  target_genes_sig <- intersect(target_genes, rownames(sig_strict))

  # 3. Combine the Top 20 with your filtered functional list for the 3rd plot
  all_genes_to_label <- unique(c(top_20_genes, target_genes_sig))

  # ---------------------------------------------------------
  # PLOT 1: Top 20 Genes Only
  # ---------------------------------------------------------
  p1 <- EnhancedVolcano::EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'log2FoldChange', y = 'padj',
                        title = 'Top 10 Up & Down',
                        subtitle = '',
                        pCutoff = 0.05, FCcutoff = 1.0,
                        selectLab = top_20_genes,
                        drawConnectors = TRUE, widthConnectors = 0.5, max.overlaps = Inf,
                        pointSize = 2.0, labSize = 3.5,
                        legendPosition = 'none') # Hide legend to save space

  # ---------------------------------------------------------
  # PLOT 2: Functional List Only (Filtered for Significance)
  # ---------------------------------------------------------
  p2 <- EnhancedVolcano::EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'log2FoldChange', y = 'padj',
                        title = 'Target Functional Genes',
                        subtitle = '',
                        pCutoff = 0.05, FCcutoff = 1.0,
                        selectLab = target_genes_sig,       # <-- Uses the newly filtered list
                        drawConnectors = TRUE, widthConnectors = 0.5, max.overlaps = Inf,
                        pointSize = 2.0, labSize = 3.5,
                        legendPosition = 'none')

  # ---------------------------------------------------------
  # PLOT 3: Combined (Top 20 + Filtered Functional)
  # ---------------------------------------------------------
  p3 <- EnhancedVolcano::EnhancedVolcano(res,
                        lab = rownames(res),
                        x = 'log2FoldChange', y = 'padj',
                        title = 'Combined Labels',
                        subtitle = '',
                        pCutoff = 0.05, FCcutoff = 1.0,
                        selectLab = all_genes_to_label,     # <-- Uses the combined filtered lists
                        drawConnectors = TRUE, widthConnectors = 0.5, max.overlaps = Inf,
                        pointSize = 2.0, labSize = 3.5,
                        legendPosition = 'right') # Keep legend on the far right

  # ---------------------------------------------------------
  # Stitch them together using patchwork
  # ---------------------------------------------------------
  combined_plot <- (p1 | p2 | p3) +
    plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
    )

  return(combined_plot)
}
