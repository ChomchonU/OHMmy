#' Title
#'
#' @param ORA_df
#' @param m_t2g
#' @param output_dir
#' @param title_prefix
#' @param top_n_global
#' @param top_n_per_cluster
#' @param log2fc_cutoff
#' @param padj_cutoff
#' @param indiv_width
#' @param indiv_height
#' @param global_width
#' @param global_height
#' @param variable_per_clus
#'
#' @returns
#' @export
#'
#' @examples
run_global_ora <- function(ORA_df,
                           m_t2g,
                           output_dir,
                           title_prefix       = "Hallmark",
                           top_n_global       = 40,
                           top_n_per_cluster  = 10,
                           log2fc_cutoff      = 0.25,
                           padj_cutoff        = 0.05,
                           indiv_width        = 12,   #Width for individual  plots
                           indiv_height       = 8,    #Height for individual  plots
                           global_width       = 15,   #Width for the summary plots
                           global_height      = 15,   #Height for the summary plots
                           variable_per_clus  = FALSE) # <-- Added boolean parameter
{

  # =============================================================
  # STEP 1: Setup
  # =============================================================
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")  # Global timestamp applied to ALL outputs

  universe_genes <- unique(ORA_df$gene)
  all_clusters   <- unique(ORA_df$cluster)

  # Detect cell_type column for later per-cell-type grouping
  ct_col <- grep("cell_type|celltype", colnames(ORA_df), ignore.case = TRUE, value = TRUE)[1]
  has_celltype <- !is.na(ct_col)

  print(paste("=== Starting ORA for:", title_prefix, "==="))
  print(paste("Universe genes:", length(universe_genes)))
  print(paste("Total clusters:", length(all_clusters)))
  if (has_celltype) print(paste("Detected cell type column:", ct_col, "- Will generate per-cell-type plots!"))

  # =============================================================
  # STEP 2: Run ORA Loop (Generates Per-Cluster Plots & Summary)
  # =============================================================
  all_ora_results <- list()
  summary_list    <- list()

  for (current_cluster in all_clusters) {
    print(paste("--------------------------------------------------"))
    print(paste("Processing Cluster:", current_cluster))

    cluster_df <- ORA_df %>%
      filter(cluster == current_cluster) %>%
      filter(!is.na(avg_log2FC))

    print(paste("  Total genes in cluster:", nrow(cluster_df)))

    # ── Gene lists ──
    up_genes <- cluster_df %>%
      filter(avg_log2FC > log2fc_cutoff, p_val_adj < padj_cutoff) %>%
      pull(gene) %>% unique()

    dn_genes <- cluster_df %>%
      filter(avg_log2FC < -log2fc_cutoff, p_val_adj < padj_cutoff) %>%
      pull(gene) %>% unique()

    print(paste("  UP genes:", length(up_genes), "| DN genes:", length(dn_genes)))

    # ── ORA Upregulated ──
    ora_up <- if (length(up_genes) >= 5) {
      tryCatch(
        enricher(
          gene          = up_genes,
          TERM2GENE     = m_t2g,
          universe      = universe_genes,
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          minGSSize     = 10,
          maxGSSize     = 500,
          qvalueCutoff  = 0.25
        ),
        error = function(e) { NULL }
      )
    } else NULL

    # ── ORA Downregulated ──
    ora_dn <- if (length(dn_genes) >= 5) {
      tryCatch(
        enricher(
          gene          = dn_genes,
          TERM2GENE     = m_t2g,
          universe      = universe_genes,
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          minGSSize     = 10,
          maxGSSize     = 500,
          qvalueCutoff  = 0.25
        ),
        error = function(e) { NULL }
      )
    } else NULL

    n_up <- if (!is.null(ora_up)) nrow(as.data.frame(ora_up)) else 0
    n_dn <- if (!is.null(ora_dn)) nrow(as.data.frame(ora_dn)) else 0

    print(paste("  Enriched UP pathways:", n_up))
    print(paste("  Enriched DN pathways:", n_dn))

    # ── Collect results ──
    result_rows <- list()
    if (n_up > 0) result_rows[["up"]] <- as.data.frame(ora_up) %>% mutate(cluster = current_cluster, Direction = "Activated")
    if (n_dn > 0) result_rows[["dn"]] <- as.data.frame(ora_dn) %>% mutate(cluster = current_cluster, Direction = "Suppressed")

    # ── G. Store Summary Data ──
    summary_list[[current_cluster]] <- data.frame(
      cluster       = current_cluster,
      status        = ifelse(n_up == 0 && n_dn == 0, "no_pathways", "success"),
      n_up_pathways = n_up,
      n_dn_pathways = n_dn,
      top_up        = ifelse(n_up > 0, result_rows[["up"]]$Description[1], "none"),
      top_dn        = ifelse(n_dn > 0, result_rows[["dn"]]$Description[1], "none")
    )

    if (n_up == 0 && n_dn == 0) {
      print(paste("  [SKIP] No enriched pathways to plot for:", current_cluster))
      next
    }

    # ── Combine & Calculate Stats for the Cluster ──
    cluster_combined <- bind_rows(result_rows) %>%
      mutate(
        GeneRatio = sapply(GeneRatio, function(x) {
          v <- as.numeric(strsplit(x, "/")[[1]])
          v[1] / v[2]
        }),
        neg_log10_padj = -log10(p.adjust),
        Count          = as.numeric(Count)
      )

    # Map cell type into the results if available
    if (has_celltype) {
      c_type_val <- ORA_df[[ct_col]][ORA_df$cluster == current_cluster][1]
      cluster_combined[[ct_col]] <- c_type_val
    }

    all_ora_results[[current_cluster]] <- cluster_combined

    # =============================================================
    # PER-CLUSTER PLOTTING (Top N per direction)
    # =============================================================
    plot_df_cluster <- cluster_combined %>%
      group_by(Direction) %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n_per_cluster) %>%
      ungroup() %>%
      arrange(desc(Direction), desc(neg_log10_padj)) %>%
      mutate(Description = ifelse(
        duplicated(Description) | duplicated(Description, fromLast = TRUE),
        paste0(Description, " (", Direction, ")"),
        Description
      )) %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))

    # ── Per-Cluster Dotplot ──
    tryCatch({
      p_dot <- ggplot(plot_df_cluster, aes(x = GeneRatio, y = Description, color = neg_log10_padj, size = Count)) +
        geom_point() +
        facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
        scale_color_viridis_c(option = "plasma", name = "-log10(p.adjust)", direction = 1) +
        scale_size_continuous(range = c(3, 10), name = "Count") +
        labs(
          title    = paste(title_prefix, "ORA Dotplot - Cluster:", current_cluster),
          subtitle = paste("Activated:", n_up, "| Suppressed:", n_dn),
          x        = "GeneRatio", y = NULL
        ) +
        theme_bw() +
        theme(
          axis.text.y      = element_text(size = 8,  face = "bold"),
          axis.text.x      = element_text(size = 9),
          axis.title.x     = element_text(size = 10, face = "bold"),
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
          plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
          strip.text       = element_text(face = "bold", size = 11),
          strip.background = element_rect(fill = "grey90", color = "grey70"),
          panel.grid.major = element_line(color = "grey92"),
          panel.grid.minor = element_blank(),
          plot.margin      = ggplot2::margin(t = 15, r = 20, b = 10, l = 10)
        )
      ggsave(paste0(output_dir, title_prefix, "_", current_cluster, "_Dotplot_", ts, ".jpg"),
             plot = p_dot, width = indiv_width, height = indiv_height, dpi = 300, limitsize = FALSE)
    }, error = function(e) { print(paste("  [ERROR] Cluster Dotplot failed:", e$message)) })

    # ── Per-Cluster Barplot ──
    tryCatch({
      p_bar <- ggplot(plot_df_cluster, aes(x = neg_log10_padj, y = Description, fill = Direction)) +
        geom_bar(stat = "identity", width = 0.7) +
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey30", linewidth = 0.7) +
        facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
        scale_fill_manual(values = c("Activated" = "firebrick", "Suppressed" = "steelblue"), guide = "none") +
        labs(
          title    = paste(title_prefix, "ORA Barplot - Cluster:", current_cluster),
          subtitle = paste("Activated:", n_up, "| Suppressed:", n_dn, "| Dashed = p.adjust 0.05"),
          x        = "-log10(p.adjust)", y = NULL
        ) +
        theme_bw() +
        theme(
          axis.text.y      = element_text(size = 8,  face = "bold"),
          axis.text.x      = element_text(size = 9),
          axis.title.x     = element_text(size = 10, face = "bold"),
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
          plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
          strip.text       = element_text(face = "bold", size = 11),
          strip.background = element_rect(fill = "grey90", color = "grey70"),
          panel.grid.major = element_line(color = "grey92"),
          panel.grid.minor = element_blank(),
          plot.margin      = ggplot2::margin(t = 15, r = 20, b = 10, l = 10)
        )
      ggsave(paste0(output_dir, title_prefix, "_", current_cluster, "_Barplot_", ts, ".jpg"),
             plot = p_bar, width = indiv_width, height = indiv_height, dpi = 300, limitsize = FALSE)
    }, error = function(e) { print(paste("  [ERROR] Cluster Barplot failed:", e$message)) })

  }

  # =============================================================
  # STEP 3: Save Summary Table
  # =============================================================
  summary_df <- bind_rows(summary_list)
  write.csv(summary_df, file = paste0(output_dir, title_prefix, "_Summary_AllClusters_", ts, ".csv"), row.names = FALSE)
  print(paste("Saved Summary table to:", output_dir))

  # Check if any results were found at all for global analysis
  if (length(all_ora_results) == 0) {
    print(paste("No significant pathways found overall for", title_prefix))
    return(NULL)
  }

  # =============================================================
  # STEP 4: Combine All Results for Global Plots
  # =============================================================
  combined_df <- bind_rows(all_ora_results)
  print(paste("Total global pathway-cluster combinations:", nrow(combined_df)))

  # =============================================================
  # STEP 5: Select Global Pathways to Display
  # =============================================================
  df_up <- combined_df %>% filter(Direction == "Activated")
  df_dn <- combined_df %>% filter(Direction == "Suppressed")

  get_top_pathways <- function(df, n_top) {
    if (nrow(df) == 0) return(character(0))
    df %>%
      group_by(Description) %>%
      summarise(
        n_clusters    = n_distinct(cluster),
        mean_neg_logp = mean(neg_log10_padj),
        .groups       = "drop"
      ) %>%
      arrange(desc(n_clusters), desc(mean_neg_logp)) %>%
      head(n_top) %>%
      pull(Description)
  }

  top_pathways_up <- get_top_pathways(df_up, top_n_global)
  top_pathways_dn <- get_top_pathways(df_dn, top_n_global)

  plot_df_up <- df_up %>% filter(Description %in% top_pathways_up)
  plot_df_dn <- df_dn %>% filter(Description %in% top_pathways_dn)

  cluster_order <- unique(ORA_df$cluster) %>%
    .[order(as.numeric(str_extract(., "\\d+")))]

  order_pathways <- function(df) {
    if (nrow(df) == 0) return(character(0))
    df %>%
      group_by(Description) %>%
      summarise(mean_neg_logp = mean(neg_log10_padj), .groups = "drop") %>%
      arrange(mean_neg_logp) %>%
      pull(Description) %>%
      unique()
  }

  pathway_order_up <- order_pathways(plot_df_up)
  pathway_order_dn <- order_pathways(plot_df_dn)

  plot_df_up <- plot_df_up %>%
    mutate(cluster     = factor(cluster,     levels = cluster_order),
           Description = factor(Description, levels = pathway_order_up))

  plot_df_dn <- plot_df_dn %>%
    mutate(cluster     = factor(cluster,     levels = cluster_order),
           Description = factor(Description, levels = pathway_order_dn))

  # =============================================================
  # STEP 6: Plot & Save Global Views
  # =============================================================
  global_theme <- theme_bw() + theme(
    axis.text.x      = element_text(size = 9, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.title.x     = element_text(size = 10, face = "bold"),
    axis.text.y      = element_text(size = 8, face = "bold"),
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "grey40"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.margin      = ggplot2::margin(t = 15, r = 10, b = 10, l = 10)
  )

  # ── Save Global Activated Plot ──
  if (nrow(plot_df_up) > 0) {
    p_up <- ggplot(plot_df_up, aes(x = cluster, y = Description, color = neg_log10_padj, size = GeneRatio)) +
      geom_point() +
      scale_color_viridis_c(option = "plasma", name = "-log10(p.adjust)", direction = 1,
                            limits = c(0, max(plot_df_up$neg_log10_padj, na.rm = TRUE))) +
      scale_size_continuous(range = c(1.5, 8), name = "GeneRatio") +
      labs(title    = paste(title_prefix, "Global ORA - Activated Pathways"),
           subtitle = paste("Top", length(top_pathways_up), "pathways | Color = -log10(p.adjust) | Size = GeneRatio"),
           x = "Cluster", y = NULL) +
      global_theme

    path_up <- paste0(output_dir, title_prefix, "_Global_Activated_Dotplot_", ts, ".jpg")
    ggsave(filename = path_up, plot = p_up, width = global_width, height = global_height, dpi = 300, limitsize = FALSE)
  }

  # ── Save Global Suppressed Plot ──
  if (nrow(plot_df_dn) > 0) {
    p_dn <- ggplot(plot_df_dn, aes(x = cluster, y = Description, color = neg_log10_padj, size = GeneRatio)) +
      geom_point() +
      scale_color_viridis_c(option = "plasma", name = "-log10(p.adjust)", direction = 1,
                            limits = c(0, max(plot_df_dn$neg_log10_padj, na.rm = TRUE))) +
      scale_size_continuous(range = c(1.5, 8), name = "GeneRatio") +
      labs(title    = paste(title_prefix, "Global ORA - Suppressed Pathways"),
           subtitle = paste("Top", length(top_pathways_dn), "pathways | Color = -log10(p.adjust) | Size = GeneRatio"),
           x = "Cluster", y = NULL) +
      global_theme

    path_dn <- paste0(output_dir, title_prefix, "_Global_Suppressed_Dotplot_", ts, ".jpg")
    ggsave(filename = path_dn, plot = p_dn, width = global_width, height = global_height, dpi = 300, limitsize = FALSE)
  }

  # =============================================================
  # STEP 7: Clustered Global Plots + Dendrogram
  # =============================================================
  print("--- Generating Clustered Global Plots with Dendrograms ---")
  library(patchwork)

  # Helper function to generate clustered plots
  generate_clustered_plot <- function(plot_df, title_label, filename_suffix) {
    if(nrow(plot_df) == 0 || length(unique(plot_df$cluster)) < 2 || length(unique(plot_df$Description)) < 2) {
      print(paste("  [SKIP] Not enough clusters or pathways to build dendrogram for", title_label))
      return(NULL)
    }

    tryCatch({
      # 1. Pivot to matrix (fill missing with 0 significance)
      mat_wide <- plot_df %>%
        dplyr::select(cluster, Description, neg_log10_padj) %>%
        tidyr::pivot_wider(names_from = cluster, values_from = neg_log10_padj, values_fill = 0)

      mat <- as.matrix(mat_wide[, -1])
      rownames(mat) <- mat_wide$Description

      # 2. Perform hierarchical clustering
      hc_cols   <- hclust(dist(t(mat)))
      hc_rows   <- hclust(dist(mat))

      col_order <- hc_cols$labels[hc_cols$order]
      row_order <- hc_rows$labels[hc_rows$order]

      # 3. Apply clustering orders to dataframe
      clustered_df <- plot_df %>%
        mutate(cluster     = factor(cluster, levels = col_order),
               Description = factor(Description, levels = row_order))

      # 4. Main Plot
      p_main <- ggplot(clustered_df, aes(x = cluster, y = Description)) +
        geom_point(aes(size = GeneRatio, fill = neg_log10_padj), shape = 21, color = "black", stroke = 0.5) +
        scale_fill_viridis_c(option = "plasma", name = "-log10(p.adj)", direction = 1) +
        scale_size_continuous(range = c(2, 8), name = "GeneRatio") +
        theme_bw() +
        labs(x = "Cluster", y = NULL) +
        theme(
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
          axis.text.y      = element_text(size = 9, face = "bold"),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          legend.position  = "right",
          plot.margin      = ggplot2::margin(t = 0, r = 10, b = 10, l = 10) # Connect to dendrogram
        )

      # 5. Dendrogram (X-axis clusters)
      dendro_data <- ggdendro::dendro_data(hc_cols, type = "rectangle")
      p_dendro <- ggplot(ggdendro::segment(dendro_data)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        scale_x_continuous(limits = c(0.5, length(col_order) + 0.5), expand = c(0, 0)) +
        theme_void() +
        theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 10))

      # 6. Combine
      p_combined <- p_dendro / p_main +
        plot_layout(heights = c(1, 10)) +
        plot_annotation(
          title = paste(title_prefix, "- Clustered ORA:", title_label),
          theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
        )

      save_path <- paste0(output_dir, title_prefix, filename_suffix, "_", ts, ".jpg")
      # Added +1 to height to give the dendrogram space without squashing the y-axis text
      ggsave(filename = save_path, plot = p_combined, width = global_width, height = global_height + 1, dpi = 300, limitsize = FALSE)
      print(paste("  [OK] Clustered dotplot saved for", title_label))

    }, error = function(e) { print(paste("  [ERROR] Clustered dotplot failed for", title_label, ":", e$message)) })
  }

  # Run function for both up and down pathways
  generate_clustered_plot(plot_df_up, "Activated Pathways", "_Global_Activated_Clustered_Dendro")
  generate_clustered_plot(plot_df_dn, "Suppressed Pathways", "_Global_Suppressed_Clustered_Dendro")


  # =============================================================
  # STEP 8: Signed Significance Dotplot (GLOBAL)
  # =============================================================
  print("--- Generating Signed Significance Dotplot (Global) ---")
  plot_df_signed <- bind_rows(plot_df_up, plot_df_dn)

  tryCatch({
    if (nrow(plot_df_signed) > 0) {

      # Since ORA has no NES, we force a "sign" based on whether the gene set was run on Up or Down genes
      plot_df_signed <- plot_df_signed %>%
        mutate(
          Signed_log10_padj = ifelse(Direction == "Activated", neg_log10_padj, -neg_log10_padj)
        )

      # Order pathways by their average signed significance to cascade top to bottom
      pathway_order_signed <- plot_df_signed %>%
        group_by(Description) %>%
        summarise(mean_sig = mean(Signed_log10_padj, na.rm = TRUE), .groups = "drop") %>%
        arrange(mean_sig) %>%
        pull(Description)

      plot_df_signed <- plot_df_signed %>%
        mutate(Description = factor(Description, levels = pathway_order_signed))

      p_signed <- ggplot(plot_df_signed, aes(x = Signed_log10_padj, y = Description)) +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
        geom_vline(xintercept = -1.301, color = "gray50", linewidth = 0.6) + # p = 0.05
        geom_vline(xintercept = 1.301, color = "gray50", linewidth = 0.6) +  # p = 0.05

        geom_point(aes(fill = cluster, size = GeneRatio),
                   shape = 21, color = "black", stroke = 0.8, alpha = 0.6) +

        labs(
          title = paste(title_prefix, "- Signed Significance (All Clusters)"),
          x = expression("Signed " * -log[10] * " adjusted P value"),
          y = NULL,
          size = "GeneRatio",
          fill = "Cluster"
        ) +
        theme_bw() +
        theme(
          axis.text.y        = element_text(size = 10, face = "bold"),
          axis.text.x        = element_text(size = 10),
          axis.title.x       = element_text(size = 11, face = "bold"),
          plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "bottom",
          legend.box         = "horizontal"
        ) +
        guides(fill = guide_legend(override.aes = list(size = 4)))

      signed_plot_path <- paste0(output_dir, title_prefix, "_ORA_Signed_Significance_Dotplot_", ts, ".jpg")
      dyn_height <- max(8, length(unique(plot_df_signed$Description)) * 0.25 + 3)
      ggsave(filename = signed_plot_path, plot = p_signed, width = global_width, height = dyn_height, dpi = 300, limitsize = FALSE)

      if (file.exists(signed_plot_path)) print(paste("  [OK] Signed Significance dotplot saved"))
    } else {
      print("  [SKIP] Not enough data for Signed Significance plot.")
    }
  }, error = function(e) { print(paste("  [ERROR] Signed Significance dotplot failed:", e$message)) })


  # =============================================================
  # STEP 8b: Signed Significance Dotplot (PER CELL TYPE)
  # =============================================================
  if (has_celltype && variable_per_clus && nrow(plot_df_signed) > 0) {
    print("--- Generating Signed Significance Dotplot Per Cell Type ---")

    unique_ctypes <- unique(plot_df_signed[[ct_col]])
    unique_ctypes <- unique_ctypes[!is.na(unique_ctypes)]

    for (c_type in unique_ctypes) {
      tryCatch({
        # Subset data for just this cell type
        df_ct <- plot_df_signed[plot_df_signed[[ct_col]] == c_type, ]

        if (nrow(df_ct) == 0) next

        # Re-order pathways just for this cell type subset
        ct_pathway_order <- df_ct %>%
          group_by(Description) %>%
          summarise(mean_sig = mean(Signed_log10_padj, na.rm = TRUE), .groups = "drop") %>%
          arrange(mean_sig) %>%
          pull(Description)

        df_ct <- df_ct %>% mutate(Description = factor(Description, levels = ct_pathway_order))

        p_ct <- ggplot(df_ct, aes(x = Signed_log10_padj, y = Description)) +
          geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
          geom_vline(xintercept = -1.301, color = "gray50", linewidth = 0.6) +
          geom_vline(xintercept = 1.301, color = "gray50", linewidth = 0.6) +
          geom_point(aes(fill = cluster, size = GeneRatio),
                     shape = 21, color = "black", stroke = 0.8, alpha = 0.6) +
          labs(
            title = paste(title_prefix, "- Signed Significance ORA (", c_type, ")"),
            x = expression("Signed " * -log[10] * " adjusted P value"),
            y = NULL,
            size = "GeneRatio",
            fill = "Cluster"
          ) +
          theme_bw() +
          theme(
            axis.text.y        = element_text(size = 10, face = "bold"),
            axis.text.x        = element_text(size = 10),
            axis.title.x       = element_text(size = 11, face = "bold"),
            plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
            panel.grid.major.y = element_blank(),
            panel.grid.minor   = element_blank(),
            legend.position    = "bottom",
            legend.box         = "horizontal"
          ) +
          guides(fill = guide_legend(override.aes = list(size = 4)))

        # Scale height dynamically for this specific plot
        ct_dyn_height <- max(5, length(unique(df_ct$Description)) * 0.25 + 3)
        safe_c_type <- gsub("[^A-Za-z0-9_]", "_", c_type)

        ct_plot_path <- paste0(output_dir, title_prefix, "_ORA_Signed_Sig_CT_", safe_c_type, "_", ts, ".jpg")
        ggsave(filename = ct_plot_path, plot = p_ct, width = 12, height = ct_dyn_height, dpi = 300, limitsize = FALSE)

        if (file.exists(ct_plot_path)) print(paste("  [OK] Per-cell-type Signed Significance dotplot saved for", c_type))

      }, error = function(e) { print(paste("  [ERROR] Per-cell-type Signed Significance plot failed for", c_type, ":", e$message)) })
    }
  }


  # ── Save Combined Data CSV ──
  write.csv(combined_df,
            file      = paste0(output_dir, title_prefix, "_Global_Combined_Data_", ts, ".csv"),
            row.names = FALSE)

  print("==================================================")
  print(paste("Finished!", title_prefix, "results saved to:", output_dir))
  return(combined_df)
}

# -------------------------------------------------------------

#' Title
#'
#' @param GSEA_df
#' @param m_t2g
#' @param output_dir
#' @param title_prefix
#' @param top_n_per_direction
#' @param padj_cutoff
#' @param dotplot_width
#' @param dotplot_height
#' @param gseaplot_width
#' @param gseaplot_height
#' @param top_n_overall
#' @param variable_per_clus
#'
#' @returns
#' @export
#'
#' @examples
run_global_gsea <- function(GSEA_df,
                            m_t2g,
                            output_dir,
                            title_prefix        = "Hallmark",
                            top_n_per_direction = 10,
                            padj_cutoff         = 0.05,
                            dotplot_width       = 12,
                            dotplot_height      = 8,
                            gseaplot_width      = 8,
                            gseaplot_height     = 6,
                            top_n_overall       = 5,
                            variable_per_clus   = FALSE) {

  # ── 1. Fix Windows parallel connection error ───────────────────
  BiocParallel::register(BiocParallel::SerialParam())

  # ── 2. Setup ───────────────────────────────────────────────────
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S") # Global timestamp for this run

  all_clusters <- unique(GSEA_df$cluster)

  # Detect cell_type column for later per-cell-type grouping
  ct_col <- grep("cell_type|celltype", colnames(GSEA_df), ignore.case = TRUE, value = TRUE)[1]
  has_celltype <- !is.na(ct_col)

  print(paste("=== Starting GSEA for:", title_prefix, "==="))
  print(paste("Output folder:", output_dir))
  print(paste("Total clusters to process:", length(all_clusters)))
  if (has_celltype) print(paste("Detected cell type column:", ct_col, "- Will generate per-cell-type plots!"))

  # =============================================================
  # PRE-PASS: Run ALL GSEA first to find global NES range
  # This ensures fixed x-axis is identical across every cluster plot
  # =============================================================
  print("--- PRE-PASS: Calculating global NES range ---")

  all_gsea_results <- list()   # store raw GSEA objects
  summary_list     <- list()   # store summary data

  for (current_cluster in all_clusters) {

    cluster_df <- GSEA_df %>%
      filter(cluster == current_cluster) %>%
      filter(!is.na(avg_log2FC)) %>%
      arrange(desc(abs(avg_log2FC))) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(desc(avg_log2FC))

    ranked_genes <- sort(setNames(cluster_df$avg_log2FC, cluster_df$gene),
                         decreasing = TRUE)

    gsea_res <- tryCatch({
      GSEA(
        geneList      = ranked_genes,
        TERM2GENE     = m_t2g,
        pvalueCutoff  = padj_cutoff,
        pAdjustMethod = "BH",
        minGSSize     = 10,
        maxGSSize     = 500,
        eps           = 1e-10,
        verbose       = FALSE
      )
    }, error = function(e) {
      print(paste("  [ERROR]", current_cluster, ":", e$message))
      return(NULL)
    })

    all_gsea_results[[current_cluster]] <- gsea_res

    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      print(paste("  [OK]", current_cluster,
                  "| pathways:", nrow(as.data.frame(gsea_res)),
                  "| NES range:", round(min(gsea_res@result$NES), 2),
                  "to", round(max(gsea_res@result$NES), 2)))
    } else {
      print(paste("  [SKIP]", current_cluster, "- no significant pathways"))
    }
  }

  # ── Calculate global NES min/max across ALL clusters ──────────
  all_nes_values <- unlist(lapply(all_gsea_results, function(res) {
    if (!is.null(res) && nrow(as.data.frame(res)) > 0) res@result$NES
  }))

  if (length(all_nes_values) == 0) {
    print(paste("No significant GSEA pathways found overall for", title_prefix))
    return(NULL)
  }

  # Add 10% padding on each side for breathing room
  nes_padding    <- 0.1 * diff(range(all_nes_values))
  global_nes_min <- floor((min(all_nes_values) - nes_padding) * 10) / 10
  global_nes_max <- ceiling((max(all_nes_values) + nes_padding) * 10) / 10

  print(paste("--- Global NES range: [", global_nes_min, ",", global_nes_max, "] ---"))

  # =============================================================
  # MAIN LOOP: Plot with fixed x-axis
  # =============================================================
  combined_df_list <- list()

  for (current_cluster in all_clusters) {

    print(paste("--------------------------------------------------"))
    print(paste("Plotting Cluster:", current_cluster))

    gsea_results <- all_gsea_results[[current_cluster]]

    # ── Skip if no results ──────────────────────────────────────
    if (is.null(gsea_results) || nrow(as.data.frame(gsea_results)) == 0) {
      print(paste("  [SKIP] No results for:", current_cluster))
      summary_list[[current_cluster]] <- data.frame(cluster = current_cluster, status = "no_pathways")
      next
    }

    # ── Calculate plotting columns ──────────────────────────────
    gsea_results@result <- gsea_results@result %>%
      mutate(
        Count          = str_count(core_enrichment, "/") + 1,
        GeneRatio      = Count / setSize,
        neg_log10_padj = -log10(p.adjust),
        Direction      = ifelse(NES > 0, "Activated", "Suppressed")
      )

    n_results    <- nrow(gsea_results@result)
    n_activated  <- sum(gsea_results@result$Direction == "Activated")
    n_suppressed <- sum(gsea_results@result$Direction == "Suppressed")

    print(paste("  Significant pathways:", n_results))
    print(paste("  Activated:", n_activated, "| Suppressed:", n_suppressed))

    # Store for combined CSV and add cell_type if it exists
    res_df <- as.data.frame(gsea_results) %>% mutate(cluster = current_cluster)
    if (has_celltype) {
      c_type_val <- GSEA_df[[ct_col]][GSEA_df$cluster == current_cluster][1]
      res_df[[ct_col]] <- c_type_val
    }
    combined_df_list[[current_cluster]] <- res_df

    summary_list[[current_cluster]] <- data.frame(
      cluster      = current_cluster,
      status       = "success",
      total_path   = n_results,
      activated    = n_activated,
      suppressed   = n_suppressed,
      top_pathway  = res_df$ID[1]
    )

    # ── Select Top N per direction ──────────────────────────────
    plot_df <- bind_rows(
      res_df %>% filter(NES > 0) %>% arrange(p.adjust) %>% head(top_n_per_direction),
      res_df %>% filter(NES < 0) %>% arrange(p.adjust) %>% head(top_n_per_direction)
    )

    if (nrow(plot_df) == 0) plot_df <- res_df %>% arrange(p.adjust) %>% head(top_n_per_direction * 2)

    # ── Order pathways by NES (negative at bottom, positive at top)
    plot_df <- plot_df %>%
      arrange(NES) %>%
      mutate(Description = factor(Description, levels = unique(Description)))

    # ── Build Single Combined Dotplot ───────────────────────────
    tryCatch({

      p_dot <- ggplot(plot_df,
                      aes(x     = NES,
                          y     = Description,
                          color = neg_log10_padj,
                          size  = GeneRatio)) +
        geom_point() +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +

        annotate("rect", xmin = global_nes_min, xmax = 0, ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.04) +
        annotate("rect", xmin = 0, xmax = global_nes_max, ymin = -Inf, ymax = Inf, fill = "firebrick", alpha = 0.04) +

        annotate("text", x = global_nes_min + 0.05, y = Inf, label = "← Suppressed", hjust = 0, vjust = 1.5, color = "steelblue", fontface = "bold", size = 3.5) +
        annotate("text", x = global_nes_max - 0.05, y = Inf, label = "Activated →", hjust = 1, vjust = 1.5, color = "firebrick", fontface = "bold", size = 3.5) +

        scale_x_continuous(
          limits = c(global_nes_min, global_nes_max),
          breaks = seq(round(global_nes_min, 0), round(global_nes_max, 0), by = 0.5)
        ) +
        scale_color_viridis_c(option = "plasma", name = "-log10(p.adjust)", direction = 1) +
        scale_size_continuous(range = c(3, 10), name = "GeneRatio") +
        labs(
          title = paste(title_prefix, "GSEA - Cluster:", current_cluster),
          x     = "NES (Normalized Enrichment Score)",
          y     = NULL
        ) +
        theme_bw() +
        theme(
          axis.text.y      = element_text(size = 8,  face = "bold"),
          axis.text.x      = element_text(size = 9),
          axis.title.x     = element_text(size = 10, face = "bold"),
          plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
          panel.grid.major = element_line(color = "grey92"),
          panel.grid.minor = element_blank(),
          legend.position  = "right",
          legend.title     = element_text(size = 9, face = "bold"),
          legend.text      = element_text(size = 8),
          plot.margin      = ggplot2::margin(t = 15, r = 20, b = 10, l = 10)
        )

      dotplot_path <- paste0(output_dir, title_prefix, "_", current_cluster, "_GSEA_NES_Dotplot_", ts, ".jpg")
      ggsave(filename = dotplot_path, plot = p_dot, width = dotplot_width, height = dotplot_height, dpi = 300 , limitsize = FALSE)

    }, error = function(e) { print(paste("  [ERROR] Dotplot failed:", e$message)) })

    # ── Save Top Pathway GSEA Enrichment Plot ───────────────────
    tryCatch({

      top_pathway <- gsea_results$ID[1]

      p_classic <- gseaplot2(
        gsea_results,
        geneSetID    = top_pathway,
        title        = paste(current_cluster, "-", top_pathway),
        pvalue_table = TRUE
      )

      safe_pathway_name <- gsub("[^A-Za-z0-9_]", "_", top_pathway)
      gseaplot_path <- paste0(output_dir, title_prefix, "_", current_cluster, "_", safe_pathway_name, "_Gseaplot_", ts, ".jpg")
      ggsave(filename = gseaplot_path, plot = p_classic, width = gseaplot_width, height = gseaplot_height, dpi = 300, limitsize = FALSE)

    }, error = function(e) { print(paste("  [ERROR] GSEAplot failed:", e$message)) })
  }

  # =============================================================
  # STEP 4: Save Summary Dataframes
  # =============================================================
  summary_df <- bind_rows(summary_list)
  write.csv(summary_df, file = paste0(output_dir, title_prefix, "_GSEA_Cluster_Summary_", ts, ".csv"), row.names = FALSE)

  combined_df <- bind_rows(combined_df_list)
  if (nrow(combined_df) > 0) {
    write.csv(combined_df, file = paste0(output_dir, title_prefix, "_GSEA_Global_Combined_Data_", ts, ".csv"), row.names = FALSE)
  }

  # =============================================================
  # DATA PREP FOR OVERALL PLOTS (Steps 5, 6, 7 & 7b)
  # =============================================================
  if (nrow(combined_df) > 0) {

    # Extract top N activated and suppressed pathways per cluster
    top_pathways_overall <- combined_df %>%
      group_by(cluster, Direction) %>%
      slice_min(order_by = p.adjust, n = top_n_overall, with_ties = FALSE) %>%
      pull(Description) %>%
      unique()

    # Filter main dataframe to only these top pathways
    overall_plot_df <- combined_df %>%
      filter(Description %in% top_pathways_overall)

    # Dynamically calculate plot height and width
    n_pathways <- length(unique(overall_plot_df$Description))
    n_clusters <- length(unique(overall_plot_df$cluster))
    dyn_width  <- max(10, (n_clusters * 0.8) + 6)
    dyn_height <- max(8, (n_pathways * 0.25) + 3)

    # =============================================================
    # STEP 5: Alphabetical Overall Dotplot
    # =============================================================
    print("--- Generating Overall Combined Plot (Alphabetical) ---")
    tryCatch({
      plot_df_alpha <- overall_plot_df %>%
        mutate(Description = factor(Description, levels = rev(sort(unique(Description)))))

      p_overall <- ggplot(plot_df_alpha, aes(x = cluster, y = Description)) +
        geom_point(aes(size = neg_log10_padj, fill = NES), shape = 21, color = "black", stroke = 0.5) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
        scale_size_continuous(range = c(2, 8), name = "-log10(p.adj)") +
        theme_bw() +
        labs(title = paste(title_prefix, "- Top GSEA Pathways Across All Clusters"), x = "Cluster", y = NULL) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
          axis.text.y = element_text(size = 9, face = "bold"),
          plot.title  = element_text(hjust = 0.5, face = "bold", size = 14)
        )

      overall_plot_path <- paste0(output_dir, title_prefix, "_GSEA_Overall_Dotplot_", ts, ".jpg")
      ggsave(filename = overall_plot_path, plot = p_overall, width = dyn_width, height = dyn_height, dpi = 300, limitsize = FALSE)
    }, error = function(e) { print(paste("  [ERROR] Overall dotplot failed:", e$message)) })

    # =============================================================
    # STEP 6: Hierarchically Clustered Plot + Dendrogram
    # =============================================================
    print("--- Generating Clustered Overall Plot with Dendrogram ---")
    tryCatch({
      nes_wide <- overall_plot_df %>%
        dplyr::select(cluster, Description, NES) %>%
        tidyr::pivot_wider(names_from = cluster, values_from = NES, values_fill = 0)

      nes_mat <- as.matrix(nes_wide[, -1])
      rownames(nes_mat) <- nes_wide$Description

      hc_cols   <- hclust(dist(t(nes_mat)))
      col_order <- hc_cols$labels[hc_cols$order]

      hc_rows   <- hclust(dist(nes_mat))
      row_order <- hc_rows$labels[hc_rows$order]

      clustered_df <- overall_plot_df %>%
        mutate(cluster     = factor(cluster, levels = col_order),
               Description = factor(Description, levels = row_order))

      p_main_clustered <- ggplot(clustered_df, aes(x = cluster, y = Description)) +
        geom_point(aes(size = neg_log10_padj, fill = NES), shape = 21, color = "black", stroke = 0.5) +
        scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0, name = "NES") +
        scale_size_continuous(range = c(2, 8), name = "-log10(p.adj)") +
        theme_bw() +
        labs(x = "Cluster", y = NULL) +
        theme(
          axis.text.x      = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
          axis.text.y      = element_text(size = 9, face = "bold"),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          legend.position  = "right",
          plot.margin      = ggplot2::margin(t = 0, r = 10, b = 10, l = 10)
        )

      dendro_data <- ggdendro::dendro_data(hc_cols, type = "rectangle")
      p_dendro <- ggplot(ggdendro::segment(dendro_data)) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        scale_x_continuous(limits = c(0.5, length(col_order) + 0.5), expand = c(0, 0)) +
        theme_void() +
        theme(plot.margin = ggplot2::margin(t = 10, r = 10, b = 0, l = 10))

      library(patchwork)
      p_combined <- p_dendro / p_main_clustered +
        plot_layout(heights = c(1, 10)) +
        plot_annotation(
          title = paste(title_prefix, "- Hierarchically Clustered GSEA Pathways"),
          theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
        )

      clustered_plot_path <- paste0(output_dir, title_prefix, "_GSEA_Clustered_Dendro_Dotplot_", ts, ".jpg")
      ggsave(filename = clustered_plot_path, plot = p_combined, width = dyn_width, height = dyn_height + 1, dpi = 300, limitsize = FALSE)

    }, error = function(e) { print(paste("  [ERROR] Clustered dotplot failed:", e$message)) })

    # =============================================================
    # STEP 7: Signed Significance Dotplot (GLOBAL Reference Match)
    # =============================================================
    print("--- Generating Signed Significance Dotplot (Global) ---")
    tryCatch({

      plot_df_signed <- overall_plot_df %>%
        mutate(
          Signed_log10_padj = sign(NES) * neg_log10_padj,
          Abs_NES = abs(NES)
        )

      # Order pathways by their average signed significance to cascade top to bottom
      pathway_order_signed <- plot_df_signed %>%
        group_by(Description) %>%
        summarise(mean_sig = mean(Signed_log10_padj, na.rm = TRUE), .groups = "drop") %>%
        arrange(mean_sig) %>%
        pull(Description)

      plot_df_signed <- plot_df_signed %>%
        mutate(Description = factor(Description, levels = pathway_order_signed))

      p_signed <- ggplot(plot_df_signed, aes(x = Signed_log10_padj, y = Description)) +
        geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
        geom_vline(xintercept = -1.301, color = "gray50", linewidth = 0.6) +
        geom_vline(xintercept = 1.301, color = "gray50", linewidth = 0.6) +
        geom_point(aes(fill = cluster, size = Abs_NES),
                   shape = 21, color = "black", stroke = 0.8, alpha = 0.5) +
        labs(
          title = paste(title_prefix, "- Signed Significance (All Clusters)"),
          x = expression("Signed " * -log[10] * " adjusted P value"),
          y = NULL,
          size = "Abs. (NES)",
          fill = "Cluster"
        ) +
        theme_bw() +
        theme(
          axis.text.y        = element_text(size = 10, face = "bold"),
          axis.text.x        = element_text(size = 10),
          axis.title.x       = element_text(size = 11, face = "bold"),
          plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "bottom",
          legend.box         = "horizontal"
        ) +
        guides(fill = guide_legend(override.aes = list(size = 4)))

      signed_plot_path <- paste0(output_dir, title_prefix, "_GSEA_Signed_Significance_Dotplot_", ts, ".jpg")
      ggsave(filename = signed_plot_path, plot = p_signed, width = dyn_width, height = dyn_height + 0.5, dpi = 300, limitsize = FALSE)

      if (file.exists(signed_plot_path)) print(paste("  [OK] Signed Significance dotplot saved"))

    }, error = function(e) { print(paste("  [ERROR] Signed Significance dotplot failed:", e$message)) })

    # =============================================================
    # STEP 7b: Signed Significance Dotplot (PER CELL TYPE)
    # =============================================================
    if (has_celltype) {
      print("--- Generating Signed Significance Dotplot Per Cell Type ---")

      # Extract unique cell types
      unique_ctypes <- unique(plot_df_signed[[ct_col]])
      unique_ctypes <- unique_ctypes[!is.na(unique_ctypes)]

      for (c_type in unique_ctypes) {
        tryCatch({
          # Subset data for just this cell type
          df_ct <- plot_df_signed[plot_df_signed[[ct_col]] == c_type, ]

          if (nrow(df_ct) == 0) next

          # Re-order pathways just for this cell type subset to cascade properly
          ct_pathway_order <- df_ct %>%
            group_by(Description) %>%
            summarise(mean_sig = mean(Signed_log10_padj, na.rm = TRUE), .groups = "drop") %>%
            arrange(mean_sig) %>%
            pull(Description)

          df_ct <- df_ct %>% mutate(Description = factor(Description, levels = ct_pathway_order))

          p_ct <- ggplot(df_ct, aes(x = Signed_log10_padj, y = Description)) +
            geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
            geom_vline(xintercept = -1.301, color = "gray50", linewidth = 0.6) +
            geom_vline(xintercept = 1.301, color = "gray50", linewidth = 0.6) +
            geom_point(aes(fill = cluster, size = Abs_NES),
                       shape = 21, color = "black", stroke = 0.8, alpha = 0.6) +
            labs(
              title = paste(title_prefix, "- Signed Significance (", c_type, ")"),
              x = expression("Signed " * -log[10] * " adjusted P value"),
              y = NULL,
              size = "Abs. (NES)",
              fill = "Cluster"
            ) +
            theme_bw() +
            theme(
              axis.text.y        = element_text(size = 10, face = "bold"),
              axis.text.x        = element_text(size = 10),
              axis.title.x       = element_text(size = 11, face = "bold"),
              plot.title         = element_text(hjust = 0.5, face = "bold", size = 14),
              panel.grid.major.y = element_blank(),
              panel.grid.minor   = element_blank(),
              legend.position    = "bottom",
              legend.box         = "horizontal"
            ) +
            guides(fill = guide_legend(override.aes = list(size = 4)))

          # Scale height dynamically for this specific plot
          ct_dyn_height <- max(5, length(unique(df_ct$Description)) * 0.25 + 3)
          safe_c_type <- gsub("[^A-Za-z0-9_]", "_", c_type)

          ct_plot_path <- paste0(output_dir, title_prefix, "_GSEA_Signed_Sig_CT_", safe_c_type, "_", ts, ".jpg")
          ggsave(filename = ct_plot_path, plot = p_ct, width = 12, height = ct_dyn_height, dpi = 300, limitsize = FALSE)

          if (file.exists(ct_plot_path)) print(paste("  [OK] Per-cell-type Signed Significance dotplot saved for", c_type))

        }, error = function(e) { print(paste("  [ERROR] Per-cell-type Signed Significance plot failed for", c_type, ":", e$message)) })
      }
    }
  }

  print("==================================================")
  print(paste("Finished!", title_prefix, "GSEA completed. Output saved to:", output_dir))

  return(combined_df)
}
# -------------------------------------------------------------

# Filter and extract top genes
#' Title
#'
#' @param res_obj
#' @param n_padj
#' @param n_lfc
#'
#' @returns
#' @export
#'
#' @examples
get_top_mixed_genes <- function(res_obj, n_padj = 40, n_lfc = 40) {
  sig <- res_obj[!is.na(res_obj$padj) &
                   res_obj$padj < 0.05 &
                   !is.na(res_obj$log2FoldChange) &
                   abs(res_obj$log2FoldChange) >= 1, ]

  sig_by_padj <- sig[order(sig$padj), ]
  top_padj_genes <- rownames(head(sig_by_padj, n_padj))

  sig_by_lfc <- sig[order(abs(sig$log2FoldChange), decreasing = TRUE), ]
  top_lfc_genes <- rownames(head(sig_by_lfc, n_lfc))

  return(unique(c(top_padj_genes, top_lfc_genes)))
}

# --------------------------------------------------------------

# Generate and save the side-by-side heatmap
#' Title
#'
#' @param res_obj
#' @param comp_name
#' @param comp_title
#' @param vsd_data
#' @param anno_col
#' @param ordered_samps
#' @param n_padj
#' @param n_lfc
#' @param out_dir
#' @param ts
#' @param clus
#'
#' @returns
#' @export
#'
#' @examples
generate_and_save_heatmap <- function(res_obj, comp_name, comp_title, vsd_data, anno_col,
                                      ordered_samps, n_padj, n_lfc, out_dir, ts, clus = "all") {

  # 1. Extract genes specifically for this comparison
  target_genes <- get_top_mixed_genes(res_obj, n_padj, n_lfc)

  if(length(target_genes) < 2) {
    message("   Not enough significant genes for ", comp_name, ". Skipping.")
    return(NULL)
  }

  # 2. Extract and order the normalized data
  heatmap_data <- assay(vsd_data)[target_genes, ordered_samps, drop = FALSE]

  # 3. Setup Labels & Titles
  custom_row_labels <- rownames(heatmap_data)
  title_A <- paste0(comp_title, "\nTop Padj = ", n_padj, " genes | Top LFC = ", n_lfc, " genes")

  # 4. Generate Plots
  p1 <- pheatmap(heatmap_data,
                 cluster_rows = TRUE,
                 cluster_cols = FALSE,
                 show_rownames = TRUE,
                 labels_row = custom_row_labels,
                 show_colnames = TRUE,
                 annotation_col = anno_col,
                 scale = "row",
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                 main = title_A,
                 silent = TRUE)

  p2 <- pheatmap(heatmap_data,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_rownames = TRUE,
                 labels_row = custom_row_labels,
                 show_colnames = TRUE,
                 annotation_col = anno_col,
                 scale = "row",
                 color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                 main = "Clustered by Gene Expression",
                 silent = TRUE)

  # 5. Stitch and Save
  combined_heatmaps <- arrangeGrob(p1$gtable, p2$gtable, ncol = 2)
  heatmap_filename <- file.path(out_dir, paste0("Heatmap_", comp_name,"_", clus, "_topP_", n_padj, "_topLFC_", n_lfc, "_", ts, ".png"))

  # DYNAMIC HEIGHT: Roughly 0.15 inches per gene, min 10, max 100
  # This prevents 2000 genes from turning into a black smudge
  plot_height <- max(10, min(100, length(target_genes) * 0.15))

  ggsave(filename = heatmap_filename,
         plot = combined_heatmaps,
         width = 16,
         height = plot_height,
         dpi = 300,
         limitsize = FALSE) # Required if height exceeds standard ggplot limits

  message("  Saved Heatmap: ", heatmap_filename)
}
