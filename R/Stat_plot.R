#' Title
#'
#' @param seurat_obj
#' @param sample_col
#' @param condition_col
#' @param celltype_col
#' @param global_test
#' @param strict_posthoc
#' @param pairwise_test
#' @param p_adjust
#' @param facet_by_cluster
#' @param facet_scales
#' @param pairwise_label
#' @param y_expand
#' @param output_dir
#' @param base_size
#' @param dpi
#'
#' @returns
#' @export
#'
#' @examples
plot_cell_abundance <- function(seurat_obj, sample_col, condition_col, celltype_col,
                                global_test = "kruskal.test",
                                strict_posthoc = TRUE, # NEW: Only calculate post-hoc if global is sig (p < 0.05)
                                pairwise_test = "mann_whitney", # Only used if n_conditions == 2
                                p_adjust = "BH",
                                facet_by_cluster = TRUE,
                                facet_scales = "free_y",
                                pairwise_label = "p.signif", # Changed default to 'p.signif' for cleaner bracket labels
                                y_expand = 0.2, # Slightly increased to make room for brackets
                                output_dir = ".",
                                base_size = 3.5,
                                dpi = 300) {

  # =======================================================================
  # AVAILABLE STATISTICAL PARAMETERS:
  #
  # global_test:   "kruskal.test" -> Uses Dunn's Test for post-hoc
  #                "anova"        -> Uses Tukey's HSD for post-hoc
  #
  # strict_posthoc: TRUE limits pairwise tests to celltypes where global test p < 0.05.
  #
  # pairwise_label:"p.adj"        (Numeric adjusted p-value)
  #                "p.format"     (Numeric raw p-value)
  #                "p.signif"     (Significance stars: ns, *, **, ***, ****)
  # =======================================================================

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # Extract and standardize metadata
  md <- seurat_obj@meta.data %>%
    select(
      Sample    = all_of(sample_col),
      Condition = all_of(condition_col),
      CellType  = all_of(celltype_col)
    ) %>%
    tidyr::drop_na()

  # Calculate proportions
  abundance_df <- md %>%
    group_by(Sample, Condition, CellType) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    tidyr::complete(tidyr::nesting(Sample, Condition), CellType, fill = list(Count = 0)) %>%
    group_by(Sample) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()

  # ---------------------------------------------------------------------------
  # Condition-aware statistics setup using rstatix
  # ---------------------------------------------------------------------------
  conditions       <- unique(as.character(abundance_df$Condition))
  celltypes        <- unique(as.character(abundance_df$CellType))
  n_conditions     <- length(conditions)

  run_global_test  <- n_conditions > 2
  effective_adjust <- if (n_conditions > 2) p_adjust else "none"

  message(sprintf(
    "Detected %d condition(s): %s\n  Global test: %s | Strict Post-hoc: %s",
    n_conditions, paste(conditions, collapse = ", "),
    if (run_global_test) global_test else "skipped",
    if (run_global_test) strict_posthoc else "N/A"
  ))

  # Compute statistics
  if (run_global_test) {
    if (global_test == "kruskal.test") {
      # Kruskal-Wallis -> Dunn's test
      global_res <- abundance_df %>% group_by(CellType) %>% kruskal_test(Proportion ~ Condition)
      posthoc_res <- abundance_df %>% group_by(CellType) %>% dunn_test(Proportion ~ Condition, p.adjust.method = effective_adjust)
    } else {
      # ANOVA -> Tukey's HSD (Tukey inherently adjusts for multiple comparisons)
      global_res <- abundance_df %>% group_by(CellType) %>% anova_test(Proportion ~ Condition)
      posthoc_res <- abundance_df %>% group_by(CellType) %>% tukey_hsd(Proportion ~ Condition)
    }

    # Filter post-hocs if strict_posthoc is TRUE
    if (strict_posthoc) {
      sig_cells <- global_res %>% filter(p < 0.05) %>% pull(CellType)
      posthoc_res <- posthoc_res %>% filter(CellType %in% sig_cells)
    }
  } else {
    # 2 Conditions -> Wilcoxon or T-test
    if (pairwise_test == "mann_whitney") {
      # Independent non-parametric (Mann-Whitney U)
      posthoc_res <- abundance_df %>%
        group_by(CellType) %>%
        wilcox_test(Proportion ~ Condition, paired = FALSE)

    } else if (pairwise_test == "wilcoxon_paired") {
      # Paired non-parametric (Wilcoxon Signed-Rank)
      posthoc_res <- abundance_df %>%
        group_by(CellType) %>%
        wilcox_test(Proportion ~ Condition, paired = TRUE)

    } else if (pairwise_test == "t_test") {
      # Independent parametric (Student's t-test)
      posthoc_res <- abundance_df %>%
        group_by(CellType) %>%
        t_test(Proportion ~ Condition, paired = FALSE)
    }
  }

  # Calculate proper Y-positions for the post-hoc brackets
  if (nrow(posthoc_res) > 0) {

    # --- ADDED FIX: Force rstatix to generate significance stars ---
    if ("p.adj" %in% colnames(posthoc_res) && !("p.adj.signif" %in% colnames(posthoc_res))) {
      posthoc_res <- posthoc_res %>% add_significance(p.col = "p.adj", output.col = "p.adj.signif")
    }
    if ("p" %in% colnames(posthoc_res) && !("p.signif" %in% colnames(posthoc_res))) {
      posthoc_res <- posthoc_res %>% add_significance(p.col = "p", output.col = "p.signif")
    }
    # ---------------------------------------------------------------

    # --- ADDED FIX: Format numeric p-values to exactly 3 decimal places ---
    format_p <- function(pval) {
      ifelse(pval < 0.001, "< 0.001", sprintf("%.3f", pval))
    }

    if ("p" %in% colnames(posthoc_res)) {
      posthoc_res$p <- format_p(as.numeric(posthoc_res$p))
    }
    if ("p.adj" %in% colnames(posthoc_res)) {
      posthoc_res$p.adj <- format_p(as.numeric(posthoc_res$p.adj))
    }
    # ----------------------------------------------------------------------

    if (facet_by_cluster) {
      # 1. Let rstatix calculate the X-axis bounds
      posthoc_res <- posthoc_res %>% add_xy_position(x = "Condition", data = abundance_df, formula = Proportion ~ Condition)

      # 2. OVERRIDE the Y-axis position dynamically per facet
      facet_maxes <- abundance_df %>%
        group_by(CellType) %>%
        summarise(facet_max = max(Proportion, na.rm = TRUE), .groups = "drop")

      posthoc_res <- posthoc_res %>%
        left_join(facet_maxes, by = "CellType") %>%
        group_by(CellType) %>%
        mutate(
          # Start 3% above the highest data point, and step up 8% for every extra bracket
          y.position = facet_max + (facet_max * 0.1) + (facet_max * 0.05 * (row_number() - 1))
        ) %>%
        ungroup()

    } else {
      posthoc_res <- posthoc_res %>% add_xy_position(x = "CellType", dodge = 0.8, data = abundance_df, formula = Proportion ~ Condition, step.increase = 0.05)
    }
  }

  # Safely map user's requested label to rstatix column names
  if (pairwise_label == "p.adj") {
    p_label <- if ("p.adj" %in% colnames(posthoc_res)) "{p.adj}" else "{p}"
  } else if (pairwise_label == "p.format") {
    p_label <- "{p}"
  } else if (pairwise_label == "p.signif") {
    p_label <- if ("p.adj.signif" %in% colnames(posthoc_res)) "{p.adj.signif}" else "{p.signif}"
  } else {
    p_label <- pairwise_label
  }

  # ---------------------------------------------------------------------------
  # Dynamic plot dimensions
  # ---------------------------------------------------------------------------
  if (facet_by_cluster) {
    n_panels    <- length(celltypes)
    n_cols      <- ceiling(sqrt(n_panels))
    n_rows      <- ceiling(n_panels / n_cols)
    plot_width  <- n_cols * base_size + 1
    plot_height <- n_rows * base_size + 1
  } else {
    n_x         <- length(celltypes)
    plot_width  <- max(n_x * base_size + 2, 6)
    plot_height <- base_size * 2 + 2
  }

  # ---------------------------------------------------------------------------
  # Build shared stat layers
  # ---------------------------------------------------------------------------
  stat_layers <- list()

  if (run_global_test) {
    # Add the global p-value text to the top
    if (facet_by_cluster) {
      stat_layers[["global"]] <- stat_compare_means(
        method   = global_test,
        aes(label = ifelse(..p.. < 0.001, "p < 0.001", sprintf("p = %.3f", ..p..))), # Removed duplicate label = "p.format"
        vjust    = -2.5,
        hide.ns  = FALSE
      )
    } else {
      stat_layers[["global"]] <- stat_compare_means(
        aes(group = Condition, label = ifelse(..p.. < 0.001, "p < 0.001", sprintf("p = %.3f", ..p..))), # Removed duplicate label = "p.format"
        method   = global_test,
        vjust    = -2.5,
        hide.ns  = FALSE
      )
    }
  }

  # Add the manual post-hoc brackets if any exist
  if (nrow(posthoc_res) > 0) {
    stat_layers[["pairwise"]] <- stat_pvalue_manual(
      posthoc_res,
      label = p_label,
      tip.length = 0.02,
      vjust    = -0.25,
      hide.ns = FALSE # Set to TRUE if you only want to see significant brackets
    )
  }

  # ---------------------------------------------------------------------------
  # Generate the plot
  # ---------------------------------------------------------------------------
  if (facet_by_cluster) {
    p <- ggplot(abundance_df, aes(x = Condition, y = Proportion)) +
      geom_boxplot(aes(fill = Condition), outlier.shape = NA, alpha = 0.4) +
      geom_jitter(aes(color = Condition),
                  position = position_jitter(width = 0.2),
                  size = 1.5, alpha = 0.8) +
      facet_wrap(~ CellType, scales = facet_scales) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            legend.position = "none") +
      labs(title = paste("Cell Type Abundance Across Conditions\n",
                         "(Sample:", sample_col,
                         "| Condition:", condition_col,
                         "| Cell Type:", celltype_col, ")"),
           y = "Fraction of Total Cells",
           x = "") +
      scale_y_continuous(expand = expansion(mult = c(0.05, y_expand))) +
      stat_layers

  } else {
    p <- ggplot(abundance_df, aes(x = CellType, y = Proportion)) +
      geom_boxplot(aes(fill = Condition), outlier.shape = NA, alpha = 0.4, position = position_dodge(0.8)) +
      geom_jitter(aes(color = Condition),
                  position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.8) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),
            legend.position = "right") +
      labs(title = paste("Cell Type Abundance Across Conditions\n",
                         "(Sample:", sample_col,
                         "| Condition:", condition_col,
                         "| Cell Type:", celltype_col, ")"),
           y = "Fraction of Total Cells",
           x = "") +
      scale_y_continuous(expand = expansion(mult = c(0.05, y_expand))) +
      stat_layers
  }

  # ---------------------------------------------------------------------------
  # Save plot
  # ---------------------------------------------------------------------------
  date_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  layout_tag <- if (facet_by_cluster) "facet_by_cluster" else "facet_by_celltype"
  file_name  <- sprintf("cell_abundance_%s_%s_%s.jpg", layout_tag, condition_col, date_stamp)
  file_path  <- file.path(output_dir, file_name)

  ggsave(
    filename = file_path,
    plot     = p,
    width    = plot_width,
    height   = plot_height,
    units    = "in",
    dpi      = dpi,
    device   = "jpeg"
  )

  message(sprintf("Plot saved: %s  [%.1f x %.1f in @ %d dpi]",
                  file_path, plot_width, plot_height, dpi))

  return(p)
}

#--------------------------------------------------

#' Title
#'
#' @param seurat_obj
#' @param sample_col
#' @param condition_col
#' @param metadata_vars
#' @param continuous_test_n2
#' @param continuous_test_n3
#' @param categorical_test
#' @param strict_posthoc
#' @param p_adjust
#' @param add_facet
#' @param output_dir
#' @param plot_width
#' @param plot_height
#' @param dpi
#'
#' @returns
#' @export
#'
#' @examples
plot_metadata_stats <- function(seurat_obj,
                                sample_col = "Sample",
                                condition_col = "Severity",
                                metadata_vars = c("Age", "Gender"),
                                continuous_test_n2 = "mann_whitney",
                                continuous_test_n3 = "kruskal.test",
                                categorical_test = "chisq",
                                strict_posthoc = TRUE, # Only run pairwise if global test p < 0.05
                                p_adjust = "BH",       # Multiple testing correction method
                                add_facet = NULL,
                                output_dir = "metadata_plots",
                                plot_width = 6,
                                plot_height = 5,
                                dpi = 300) {

  # =======================================================================
  # AVAILABLE STATISTICAL PARAMETERS:
  #
  # continuous_test_n2 (For 2 conditions):
  #   - "mann_whitney" -> Non-parametric (Wilcoxon rank-sum)
  #   - "t_test"       -> Parametric (Student's t-test)
  #
  # continuous_test_n3 (For 3+ conditions):
  #   - "kruskal.test" -> Non-parametric global test
  #   - "anova"        -> Parametric global test
  #
  # categorical_test (For character/factor variables):
  #   - "fisher"       -> Fisher's Exact Test (Best for small/imbalanced n)
  #   - "chisq"        -> Pearson's Chi-Square Test (Needs expected counts >= 5)
  # =======================================================================

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  # 1. Extract and deduplicate metadata
  cols_to_keep <- c(sample_col, condition_col, metadata_vars)
  if (!is.null(add_facet)) cols_to_keep <- c(cols_to_keep, add_facet)

  patient_md <- seurat_obj@meta.data %>%
    select(all_of(cols_to_keep)) %>%
    distinct()

  plot_list <- list()
  stat_results <- list()

  # 2. Loop through variables
  for (var in metadata_vars) {

    if (!var %in% colnames(patient_md) || all(is.na(patient_md[[var]]))) {
      message("Skipping ", var, ": Missing or all NAs")
      next
    }

    # Remove NAs for accurate group counting
    clean_md <- patient_md %>% filter(!is.na(.data[[var]]) & !is.na(.data[[condition_col]]))
    n_groups <- length(unique(clean_md[[condition_col]]))

    if (n_groups < 2) {
      message("Skipping ", var, ": Less than 2 conditions found")
      next
    }

    is_continuous <- is.numeric(clean_md[[var]])

    # Grouping logic for facet-aware stats
    if (!is.null(add_facet)) {
      stat_md <- clean_md %>% group_by(.data[[add_facet]])
    } else {
      stat_md <- clean_md
    }

    # =========================================================================
    # CONTINUOUS VARIABLES
    # =========================================================================
    if (is_continuous) {
      test_formula <- as.formula(paste(var, "~", condition_col))

      # Determine global and pairwise methods
      if (n_groups == 2) {
        if (continuous_test_n2 == "mann_whitney") {
          global_method <- "wilcox.test"
          posthoc_res <- stat_md %>% wilcox_test(test_formula) %>% add_significance()
        } else {
          global_method <- "t.test"
          posthoc_res <- stat_md %>% t_test(test_formula) %>% add_significance()
        }
        global_res <- posthoc_res # For 2 groups, global is the pairwise

      } else {
        if (continuous_test_n3 == "kruskal.test") {
          global_method <- "kruskal.test"
          global_res <- stat_md %>% kruskal_test(test_formula)
          posthoc_res <- stat_md %>% dunn_test(test_formula, p.adjust.method = p_adjust)
        } else {
          global_method <- "anova"
          global_res <- stat_md %>% anova_test(test_formula)
          posthoc_res <- stat_md %>% tukey_hsd(test_formula)
        }
      }

      # Filter post-hocs if strict
      if (n_groups > 2 && strict_posthoc) {
        if (!is.null(add_facet)) {
          sig_facets <- global_res %>% filter(p < 0.05) %>% pull(add_facet)
          posthoc_res <- posthoc_res %>% filter(.data[[add_facet]] %in% sig_facets)
        } else {
          if (global_res$p >= 0.05) posthoc_res <- posthoc_res[0, ] # Empty if not sig
        }
      }

      # Ensure rstatix generates significance stars for adjusted P
      if (nrow(posthoc_res) > 0) {
        if ("p.adj" %in% colnames(posthoc_res) && !("p.adj.signif" %in% colnames(posthoc_res))) {
          posthoc_res <- posthoc_res %>% add_significance(p.col = "p.adj", output.col = "p.adj.signif")
        }

        # Calculate Y-positions for brackets automatically
        if (!is.null(add_facet)) {
          posthoc_res <- posthoc_res %>% add_xy_position(x = condition_col, dodge = 0.8, data = clean_md, formula = test_formula, step.increase = 0.1)
        } else {
          posthoc_res <- posthoc_res %>% add_xy_position(x = condition_col, data = clean_md, formula = test_formula, step.increase = 0.1)
        }
      }

      stat_results[[var]] <- list(Global = global_res, PostHoc = posthoc_res)

      # Create continuous plot
      p <- ggboxplot(
        clean_md,
        x = condition_col,
        y = var,
        fill = condition_col,
        palette = "igv",
        add = "jitter",
        title = paste("Patient", var, "by", condition_col),
        ylab = var,
        xlab = condition_col,
        facet.by = add_facet
      ) +
        theme_bw() +
        theme(legend.position = "none") +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) # Extra headroom for brackets

      # --- FIX: Force the correct test name on the plot ---
      display_method <- if (global_method == "wilcox.test") {
        "Mann-Whitney U"
      } else if (global_method == "kruskal.test") {
        "Kruskal-Wallis"
      } else if (global_method == "t.test") {
        "Student's t-test"
      } else {
        "ANOVA"
      }

      # Notice the !! before display_method
      p <- p + stat_compare_means(
        method = global_method,
        aes(label = paste0(!!display_method, ", p = ", ..p.format..)),
        vjust = -2
      )
      # ----------------------------------------------------

      # Add brackets if any post-hocs exist
      if (nrow(posthoc_res) > 0) {
        p_label <- if ("p.adj.signif" %in% colnames(posthoc_res)) "{p.adj.signif}" else "{p.signif}"
        p <- p + stat_pvalue_manual(posthoc_res, label = p_label, tip.length = 0.02, hide.ns = TRUE)
      }

      # =========================================================================
      # CATEGORICAL VARIABLES
      # =========================================================================
    } else {

      run_cat_test <- function(tab, test_type) {
        if (test_type == "fisher") {
          tryCatch(fisher.test(tab), error = function(e) fisher.test(tab, simulate.p.value = TRUE))
        } else {
          suppressWarnings(chisq.test(tab))
        }
      }

      # If faceted, we must loop through facets to compute categorical stats
      facet_levels <- if (!is.null(add_facet)) unique(clean_md[[add_facet]]) else "All"

      cat_global_list <- list()
      cat_posthoc_list <- list()

      for (fl in facet_levels) {
        sub_md <- if (fl == "All") clean_md else clean_md %>% filter(.data[[add_facet]] == fl)
        cont_tab <- table(sub_md[[condition_col]], sub_md[[var]])

        # Global Test
        g_test <- run_cat_test(cont_tab, categorical_test)
        g_df <- data.frame(Facet = fl, Method = categorical_test, p = g_test$p.value)
        cat_global_list[[fl]] <- g_df

        # Post-hoc Test (if > 2 groups)
        if (n_groups > 2) {
          if (!strict_posthoc || g_test$p.value < 0.05) {
            pairs <- combn(rownames(cont_tab), 2, simplify = FALSE)
            pw_res <- lapply(pairs, function(pair) {
              pair_tab <- cont_tab[pair, ]
              pw_test <- run_cat_test(pair_tab, categorical_test)
              data.frame(Facet = fl, Group1 = pair[1], Group2 = pair[2], p = pw_test$p.value)
            })
            pw_df <- bind_rows(pw_res)
            pw_df$p.adj <- p.adjust(pw_df$p, method = p_adjust)
            cat_posthoc_list[[fl]] <- pw_df
          }
        }
      }

      global_res <- bind_rows(cat_global_list)
      posthoc_res <- if (length(cat_posthoc_list) > 0) bind_rows(cat_posthoc_list) else data.frame()

      stat_results[[var]] <- list(Global = global_res, PostHoc = posthoc_res)

      # Determine subtitle p-value
      if (length(facet_levels) == 1) {
        p_val <- global_res$p[1]
        p_label <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
        sub_title <- paste(ifelse(categorical_test=="fisher", "Fisher's Exact", "Chi-Square"), "Test (Global):", p_label)
      } else {
        sub_title <- "See stats output for per-facet p-values"
      }

      # Create categorical plot
      p <- ggplot(clean_md, aes(x = .data[[condition_col]], fill = .data[[var]])) +
        geom_bar(position = "fill", color = "black", alpha = 0.8) +
        theme_bw() +
        labs(
          title = paste("Patient", var, "by", condition_col),
          x = condition_col,
          y = "Proportion",
          fill = var,
          subtitle = sub_title
        ) +
        scale_y_continuous(labels = percent_format()) +
        scale_fill_npg()

      if (!is.null(add_facet)) {
        p <- p + facet_wrap(as.formula(paste("~", add_facet)))
      }
    }

    plot_list[[var]] <- p

    # =========================================================================
    # SAVE PLOTS
    # =========================================================================
    date_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    layout_tag <- if (!is.null(add_facet)) paste0("faceted_by_", add_facet) else "unfaceted"

    file_name  <- sprintf("metadata_%s_by_%s_%s_%s.jpg", var, condition_col, layout_tag, date_stamp)
    file_path  <- file.path(output_dir, file_name)

    ggsave(
      filename = file_path,
      plot     = p,
      width    = plot_width,
      height   = plot_height,
      units    = "in",
      dpi      = dpi,
      device   = "jpeg"
    )

    message(sprintf("Plot saved: %s", file_path))
  }

  return(list(plots = plot_list, stats = stat_results))
}
