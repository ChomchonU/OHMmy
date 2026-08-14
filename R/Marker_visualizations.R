#' Title
#'
#' @param df
#' @param gene_list
#' @param gene_regex
#' @param id_col
#' @param output_dir
#' @param n_splits
#' @param plot_title_prefix
#' @param width_scale
#' @param height_scale
#' @param min_width
#' @param min_height
#' @param label_space
#'
#' @returns
#' @export
#'
#' @examples
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
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggdendro)
  library(patchwork)
  library(rlang)  # For tidy evaluation
  library(tibble)

  # Convert id_col to symbol for tidy evaluation
  id_sym <- sym(id_col)

  # 1. Filter genes by list or regex
  df_filtered <- df %>%
    mutate(gene = as.character(gene)) %>%
    filter(
      (if (!is.null(gene_list)) gene %in% gene_list else TRUE) &
        (if (!is.null(gene_regex)) str_detect(gene, regex(gene_regex, ignore_case = TRUE)) else TRUE)
    ) %>%
    dplyr::select(gene, cluster, {{id_sym}}, avg_log2FC, diff) %>%
    mutate(
      across(c(gene), as.character),
      y_id = as.character(!!id_sym),
      gene_cluster = paste0(gene, " (C", cluster, ")")
    ) %>%
    distinct(gene, cluster, gene_cluster, y_id, avg_log2FC, diff)

  if (nrow(df_filtered) == 0) {
    message(" No genes matched.")
    return(invisible(NULL))
  }

  # 2. Order genes
  gene_order_only <- df_filtered %>%
    group_by(gene) %>%
    summarize(avg = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
    arrange(avg) %>%
    pull(gene)

  df_filtered <- df_filtered %>%
    mutate(
      gene_cluster_label = gene_cluster,
      gene_clean = gsub(" \\(C\\d+\\)$", "", gene_cluster_label),
      cluster_num = as.integer(gsub("^.*\\(C(\\d+)\\)$", "\\1", gene_cluster_label)),
      gene_clean = factor(gene_clean, levels = gene_order_only)
    ) %>%
    arrange(gene_clean, cluster_num)

  df_filtered$gene_cluster <- factor(df_filtered$gene_cluster_label, levels = unique(df_filtered$gene_cluster_label))

  # 3. Global dendrogram
  mat_all <- df_filtered %>%
    dplyr::select(y_id, gene_cluster, avg_log2FC) %>%
    pivot_wider(names_from = gene_cluster, values_from = avg_log2FC, values_fill = 0) %>%
    column_to_rownames("y_id") %>%
    as.matrix()

  dendro <- stats::hclust(dist(mat_all), method = "ward.D2")
  dendro_data <- ggdendro::dendro_data(as.dendrogram(dendro), type = "rectangle")
  y_order <- dendro$labels[dendro$order]

  # Save global dendrogram
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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

  # 4. Split gene-cluster labels into chunks
  gene_cluster_levels <- levels(df_filtered$gene_cluster)
  chunks <- split(gene_cluster_levels, cut(seq_along(gene_cluster_levels), breaks = n_splits, labels = paste0("Part", 1:n_splits)))

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
}

# -------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param sample_name
#' @param marker_diff_thresh
#' @param marker_pval_adj
#' @param marker_avg_log2FC_thresh
#' @param top_n
#' @param output_dir_base
#' @param plot_format
#' @param width
#' @param height
#' @param dpi
#' @param compare
#' @param add_timestamp
#' @param onlyPos
#'
#' @returns
#' @export
#'
#' @examples
FindTopMarkersAndHeatmap <- function(
    seurat_obj,
    sample_name = "Sample",
    marker_diff_thresh = 0.1,
    marker_pval_adj = 0.05,
    marker_avg_log2FC_thresh = 0.5,
    top_n = 20,
    output_dir_base = "Plots_heatmap",
    plot_format = "jpg",  # "png" or "pdf"
    width = 10,
    height = 20,
    dpi = 300,
    compare = NULL,
    add_timestamp = TRUE,
    onlyPos = TRUE

) {
  require(Seurat)
  require(dplyr)
  require(ggplot2)
  require(lubridate)

  # Helper to sanitize filenames
  sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

  message("[", sample_name, "] Finding markers...")

  if(is.null(compare)) {
    markers <- FindAllMarkers(seurat_obj, only.pos = onlyPos, logfc.threshold = 0, min.pct = 0)
  } else {
    markers <- FindMarkers(seurat_obj, ident.1 = compare[1], ident.2 = compare[2], only.pos = onlyPos, logfc.threshold = 0, min.pct = 0)

    # --- CRITICAL PATCH FOR PAIRWISE COMPARISON ---
    # FindMarkers lacks 'gene' and 'cluster' columns. We must create them.
    markers$gene <- rownames(markers)

    # Assign the cluster name based on the direction of the fold change
    if("avg_log2FC" %in% colnames(markers)) {
      markers$cluster <- ifelse(markers$avg_log2FC > 0, compare[1], compare[2])
    }
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
  heatmap <- DoHeatmap(seurat_obj, features = top_markers$gene) +
    ggtitle(paste0("Top Markers – ", sample_name))

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

#' Title
#'
#' @param df
#' @param genes
#' @param id_col
#' @param output_dir
#' @param plot_title
#' @param dendro_side
#' @param label_space
#' @param width_scale
#' @param height_scale
#' @param min_width
#' @param min_height
#'
#' @returns
#' @export
#'
#' @examples
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
  # Check if id_col exists
  if (!(id_col %in% colnames(df))) {
    stop("The specified id_col '", id_col, "' does not exist in the data.")
  }

  # Determine if genes is a regex
  is_regex <- is.character(genes) && length(genes) == 1 && grepl("[\\^\\$|]", genes)

  # 1. Filter and prepare
  df_igh <- if (is_regex) {
    df %>%
      filter(grepl(genes, gene, ignore.case = TRUE)) %>%
      dplyr::select(gene, cluster, !!sym(id_col), avg_log2FC, diff)
  } else {
    df %>%
      filter(gene %in% genes) %>%
      dplyr::select(gene, cluster, !!sym(id_col), avg_log2FC, diff)
  }

  df_igh <- df_igh %>%
    mutate(across(c(gene, !!sym(id_col)), as.character)) %>%
    rename(y_id = !!sym(id_col))

  # 2. Aggregate per gene and y_id
  df_igh <- df_igh %>%
    group_by(gene, y_id) %>%
    summarise(
      avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
      diff = mean(diff, na.rm = TRUE),
      cluster = paste(sort(unique(cluster)), collapse = ","),
      .groups = "drop"
    )

  # 3. Hierarchical clustering on y_id
  mat <- df_igh %>%
    dplyr::select(y_id, gene, avg_log2FC) %>%
    pivot_wider(names_from = gene, values_from = avg_log2FC, values_fill = 0) %>%
    column_to_rownames("y_id") %>%
    as.matrix()

  dendro <- stats::hclust(dist(mat), method = "ward.D2")
  ordered_y <- dendro$labels[dendro$order]
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
    scale_y_continuous(expand = expansion(add = c(0.5, 0.5))) +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  # 8. Combine
  if (tolower(dendro_side) == "left") {
    combined_plot <- p_dendro + p_dot + patchwork::plot_layout(widths = c(1, 5))
  } else {
    combined_plot <- p_dot + p_dendro + patchwork::plot_layout(widths = c(5, 1))
  }

  # 9. Save
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  dot_path <- file.path(output_dir, paste0("DotPlot_", timestamp, ".png"))
  dendro_path <- file.path(output_dir, paste0("Dendrogram_", timestamp, ".png"))
  combined_path <- file.path(output_dir, paste0("DotPlot_with_Dendrogram_", timestamp, ".png"))

  ggsave(dot_path, p_dot, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
  ggsave(dendro_path, p_dendro, width = 10, height = plot_height, dpi = 300)
  ggsave(combined_path, combined_plot, width = plot_width + 10, height = plot_height, dpi = 300, limitsize = FALSE)

  message("Dot plot saved: ", dot_path)
  message("Dendrogram saved: ", dendro_path)
  message("Combined plot saved: ", combined_path)

  invisible(list(dotplot = p_dot, dendrogram = p_dendro, combined = combined_plot))
}

# ---------------------------------------------------------

# Function to extract, average, and bin gene expression
#' Title
#'
#' @param seurat_obj
#' @param gene_list
#' @param group_col
#'
#' @returns
#' @export
#'
#' @examples
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

#' Title
#'
#' @param res
#' @param main_title
#' @param target_genes
#'
#' @returns
#' @export
#'
#' @examples
generate_volcano_trio <- function(res, main_title, target_genes) {

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
  p1 <- EnhancedVolcano(res,
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
  p2 <- EnhancedVolcano(res,
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
  p3 <- EnhancedVolcano(res,
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
