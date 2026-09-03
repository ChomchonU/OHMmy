# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.find_trough <- function(x, min_n = 50, adjust = 1, n_grid = 512) {
  x <- x[is.finite(x)]
  if (length(x) < min_n) return(list(cut = NA_real_, fallback = TRUE))
  dens <- try(stats::density(x, adjust = adjust, n = n_grid), silent = TRUE)
  if (inherits(dens, "try-error")) return(list(cut = NA_real_, fallback = TRUE))
  xs <- dens$x
  ys <- dens$y
  pk <- which(diff(sign(diff(ys))) == -2L) + 1L
  if (length(pk) < 2L) return(list(cut = NA_real_, fallback = TRUE))
  top2 <- sort(pk[order(ys[pk], decreasing = TRUE)][1:2])
  seg  <- seq.int(top2[1], top2[2])
  list(cut = xs[seg][which.min(ys[seg])], fallback = FALSE)
}

#' @keywords internal
#' @noRd
.find_gmm <- function(x, min_n = 50, n_grid = 1001) {
  fail <- list(cut = NA_real_, fallback = TRUE)
  if (!requireNamespace("mclust", quietly = TRUE)) return(fail)
  x <- x[is.finite(x)]
  if (length(x) < min_n) return(fail)
  fit <- try(mclust::Mclust(x, G = 2, verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) return(fail)
  if (is.null(fit$G) || fit$G != 2L)              return(fail)
  mu <- fit$parameters$mean
  if (length(mu) != 2L || !all(is.finite(mu)) || diff(range(mu)) == 0) return(fail)
  grid <- seq(min(mu), max(mu), length.out = n_grid)
  z    <- try(stats::predict(fit, newdata = grid)$z, silent = TRUE)
  if (inherits(z, "try-error") || is.null(z)) return(fail)
  p   <- z[, which.max(mu)]
  idx <- which(diff(sign(p - 0.5)) != 0)[1]
  if (is.na(idx)) return(fail)
  list(cut = grid[idx], fallback = FALSE)
}

#' @keywords internal
#' @noRd
.zs <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(x - mean(x, na.rm = TRUE))
  (x - mean(x, na.rm = TRUE)) / s
}

# ---------------------------------------------------------------------------
# Main Function
# ---------------------------------------------------------------------------

#' Calculate Consensus Labels for CD8 vs NK Cells
#'
#' Evaluates competing CD8 scoring variants against a baseline NK score to assign discrete
#' cell labels (CD8 vs NK). It establishes a cutoff threshold (using fixed, density trough,
#' or Gaussian Mixture Model methods), calculates consensus between a primary and confirmatory
#' variant, optionally bootstraps for stability, and generates diagnostic plots.
#'
#' @param obj A single-cell object (e.g., Seurat) that contains a `@meta.data` slot.
#' @param cd8_cols Character vector; names of the columns in `obj@meta.data` containing the CD8 scores to test.
#' @param nk_col Character string; name of the column containing the baseline NK score. Default: `"NK_Score"`.
#' @param cells_use Character vector; specific cell barcodes to include. Defaults to all cells in `obj`.
#' @param method Character; the method to determine the CD8 vs NK decision boundary. One of `"fixed"`, `"trough"` (density valley), or `"gmm"` (Gaussian Mixture Model).
#' @param fixed_cut Numeric; the threshold to use if `method = "fixed"`. Default: `0`.
#' @param zscore Logical; whether to standardize (Z-score) the CD8 and NK columns independently before subtracting them. Default: `TRUE`.
#' @param primary Character string; the CD8 variant to act as the primary label for the consensus calculation. Defaults to `cd8_cols[1]`.
#' @param confirm Character string; the CD8 variant to act as the confirming label for the consensus calculation. Defaults to `cd8_cols[2]`.
#' @param min_features Numeric; optional threshold to filter out low-quality cells prior to scoring.
#' @param feature_col Character string; the metadata column used for `min_features` filtering. Default: `"nFeature_RNA"`.
#' @param n_boot Numeric; number of bootstrap iterations to run to test cutoff/proportion stability. Set to `0` to skip.
#' @param boot_target Character; what metric to track during bootstrapping. Either `"cutoff"` (the threshold value) or `"proportion"` (fraction of cells labeled CD8).
#' @param add_meta Logical; whether to append the computed differences, scaled scores, and labels back into `obj@meta.data`.
#' @param prefix Character string; the prefix applied to the new columns added to `meta.data`. Default: `"CD8vNK"`.
#' @param seed Numeric; the random seed used for bootstrapping and plot downsampling. Default: `42`.
#' @param make_plots Logical; whether to generate diagnostic plots (Violin distributions, Agreement Heatmaps, etc.).
#' @param plot_group_by Character; how to group/color cells in the violin plots. Either `"consensus"` or by individual `"variant_label"`.
#' @param plot_n Numeric; maximum number of cells per class to sample for violin plots to ensure fast rendering. Default: `1000`.
#' @param out_dir Character string; optional directory path to save the generated plots (e.g., `"./plots/"`). If `NULL`, plots are only returned in the output list.
#'
#' @returns A list containing:
#' \itemize{
#'   \item \code{obj}: The updated single-cell object with new metadata columns (if `add_meta = TRUE`).
#'   \item \code{per}: A list of data frames with detailed, per-cell scoring and labeling for each variant.
#'   \item \code{labels}: A matrix of the assigned labels (CD8/NK) for all cells across all variants.
#'   \item \code{consensus}: A named character vector of the final consensus labels per cell.
#'   \item \code{agreement}: A pairwise matrix showing the proportion of label agreement between variants.
#'   \item \code{n_disagree}: Integer; the number of cells that had conflicting labels among any variants.
#'   \item \code{cutoffs}: A data frame summarizing the calculated thresholds and statistics for each variant.
#'   \item \code{bootstrap}: A list containing bootstrap draws and summary statistics (if `n_boot > 0`).
#'   \item \code{plots}: A list of `ggplot2` objects generated during the run (if `make_plots = TRUE`).
#'   \item \code{cells_used}: Character vector of valid cell IDs retained after filtering.
#'   \item \code{dropped}: Character vector of cell IDs dropped during filtering.
#'   \item \code{params}: A list recording the input parameters used for this run.
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with a fixed cutoff of 0 (on Z-scored data)
#' res <- score_consensus(
#'   obj = seurat_obj,
#'   cd8_cols = c("CD8_Score_TCR", "CD8_Score_GEX"),
#'   nk_col = "NK_Score",
#'   method = "fixed",
#'   fixed_cut = 0,
#'   make_plots = TRUE,
#'   out_dir = "./scoring_plots"
#' )
#'
#' # Extract the updated object for downstream Seurat workflows
#' seurat_obj <- res$obj
#' }
score_consensus <- function(obj,
                            cd8_cols,
                            nk_col        = "NK_Score",
                            cells_use     = colnames(obj),
                            method        = c("fixed", "trough", "gmm"),
                            fixed_cut     = 0,
                            zscore        = TRUE,
                            primary       = NULL,
                            confirm       = NULL,
                            min_features  = NULL,
                            feature_col   = "nFeature_RNA",
                            n_boot        = 0,
                            boot_target   = c("cutoff", "proportion"),
                            add_meta      = TRUE,
                            prefix        = "CD8vNK",
                            seed          = 42,
                            # Plotting arguments
                            make_plots    = FALSE,
                            plot_group_by = c("consensus", "variant_label"),
                            plot_n        = 1000,
                            out_dir       = NULL) {

  ## ---- 1. Validate -------------------------------------------------------
  method        <- match.arg(method)
  boot_target   <- match.arg(boot_target)
  plot_group_by <- match.arg(plot_group_by)

  md <- obj@meta.data

  stopifnot(is.character(cd8_cols), length(cd8_cols) >= 1L,
            is.character(nk_col),   length(nk_col)   == 1L)

  miss <- setdiff(c(cd8_cols, nk_col), colnames(md))
  if (length(miss))
    stop("Missing meta.data column(s): ", paste(miss, collapse = ", "), call. = FALSE)

  if (is.null(primary)) primary <- cd8_cols[1]
  if (is.null(confirm))
    confirm <- if (length(cd8_cols) >= 2L) cd8_cols[2] else primary

  if (!primary %in% cd8_cols) stop("'primary' not in cd8_cols.", call. = FALSE)
  if (!confirm %in% cd8_cols) stop("'confirm' not in cd8_cols.", call. = FALSE)

  if (n_boot > 0 && method == "fixed" && boot_target == "cutoff") {
    warning("method = 'fixed': cutoff is not data-estimated, so its bootstrap SD is 0. Switching boot_target to 'proportion'.", call. = FALSE)
    boot_target <- "proportion"
  }

  ## ---- 2. Select cells (filter BEFORE standardising) ---------------------
  cells   <- intersect(cells_use, rownames(md))
  dropped <- character(0)

  if (!is.null(min_features)) {
    if (!feature_col %in% colnames(md)) stop("feature_col '", feature_col, "' not found.", call. = FALSE)
    val     <- md[cells, feature_col]
    keep    <- cells[!is.na(val) & val >= min_features]
    dropped <- setdiff(cells, keep)
    cells   <- keep
  }

  if (length(cells) < 3L) stop("Fewer than 3 cells remain after filtering.", call. = FALSE)

  ## ---- 3. Cutoff dispatcher ----------------------------------------------
  .get_cut <- function(d) {
    r <- switch(method,
                fixed  = list(cut = fixed_cut, fallback = FALSE),
                trough = .find_trough(d),
                gmm    = .find_gmm(d))
    if (!is.finite(r$cut)) r <- list(cut = fixed_cut, fallback = TRUE)
    r
  }

  ## ---- 4. Per-variant scoring --------------------------------------------
  nk_raw <- as.numeric(md[cells, nk_col])
  nk_use <- if (zscore) .zs(nk_raw) else nk_raw

  per <- lapply(cd8_cols, function(cc) {
    cd8_raw <- as.numeric(md[cells, cc])
    cd8_use <- if (zscore) .zs(cd8_raw) else cd8_raw

    d  <- cd8_use - nk_use
    ct <- .get_cut(d)

    data.frame(cell      = cells,
               variant   = cc,
               cd8       = cd8_raw,
               nk        = nk_raw,
               cd8_scale = cd8_use,
               nk_scale  = nk_use,
               diff      = d,
               cutoff    = ct$cut,
               fallback  = ct$fallback,
               label     = ifelse(d > ct$cut, "CD8", "NK"),
               stringsAsFactors = FALSE)
  })
  names(per) <- cd8_cols

  ## ---- 5. Label matrix, consensus, agreement -----------------------------
  labels <- vapply(per, `[[`, character(length(cells)), "label")
  labels <- matrix(labels, nrow = length(cells), dimnames = list(cells, cd8_cols))

  lp <- labels[, primary]
  lc <- labels[, confirm]

  consensus <- ifelse(
    is.na(lp) | is.na(lc),            NA_character_,
    ifelse(lp == lc,                  lp,
           ifelse(lp == "CD8",               "Ambiguous_CD8_primary", "Ambiguous_CD8_confirm"))
  )
  names(consensus) <- cells

  disagree_dir <- factor(
    ifelse(is.na(lp) | is.na(lc), NA_character_,
           ifelse(lp == lc, "agree", ifelse(lp == "CD8", "primary_only", "confirm_only"))),
    levels = c("agree", "primary_only", "confirm_only")
  )
  names(disagree_dir) <- cells

  primary_vs_confirm <- table(primary = lp, confirm = lc, useNA = "ifany")

  k <- length(cd8_cols)
  agreement <- matrix(NA_real_, k, k, dimnames = list(cd8_cols, cd8_cols))
  for (i in seq_len(k)) {
    agreement[i, i] <- 1
    for (j in seq_len(i - 1L)) {
      ok <- !is.na(labels[, i]) & !is.na(labels[, j])
      a  <- if (any(ok)) mean(labels[ok, i] == labels[ok, j]) else NA_real_
      agreement[i, j] <- agreement[j, i] <- a
    }
  }

  n_disagree <- sum(apply(labels, 1L, function(r) {
    r <- r[!is.na(r)]
    length(r) > 1L && length(unique(r)) > 1L
  }))

  ## ---- 6. Cutoff summary --------------------------------------------------
  cut_tab <- do.call(rbind, lapply(cd8_cols, function(cc) {
    p <- per[[cc]]
    data.frame(variant   = cc,
               cutoff    = p$cutoff[1],
               fallback  = p$fallback[1],
               med_diff  = stats::median(p$diff, na.rm = TRUE),
               q05_diff  = unname(stats::quantile(p$diff, 0.05, na.rm = TRUE)),
               q95_diff  = unname(stats::quantile(p$diff, 0.95, na.rm = TRUE)),
               frac_cd8  = mean(p$label == "CD8"),
               stringsAsFactors = FALSE)
  }))
  rownames(cut_tab) <- NULL

  ## ---- 7. Bootstrap -------------------------------------------------------
  bootstrap <- NULL
  if (n_boot > 0) {
    set.seed(seed)
    boot_list <- lapply(cd8_cols, function(cc) {
      d       <- per[[cc]]$diff
      d       <- d[is.finite(d)]
      cut_obs <- cut_tab$cutoff[cut_tab$variant == cc]

      draws_list <- lapply(seq_len(n_boot), function(i) {
        db <- sample(d, length(d), replace = TRUE)
        if (method == "fixed") {
          list(cut = cut_obs, fb = FALSE, prop = mean(db > cut_obs))
        } else {
          res <- .get_cut(db)
          list(cut = res$cut, fb = res$fallback, prop = mean(db > res$cut))
        }
      })

      if (boot_target == "cutoff") {
        draws <- vapply(draws_list, `[[`, numeric(1), "cut")
        obs_val <- cut_obs
      } else {
        draws <- vapply(draws_list, `[[`, numeric(1), "prop")
        obs_val <- mean(d > cut_obs)
      }

      fb_rate <- if (method != "fixed") mean(vapply(draws_list, `[[`, logical(1), "fb")) else NA_real_

      list(draws = draws,
           summary = data.frame(variant = cc, target = boot_target, observed = obs_val,
                                mean = mean(draws), sd = stats::sd(draws),
                                ci_lo = unname(stats::quantile(draws, 0.025, na.rm = TRUE)),
                                ci_hi = unname(stats::quantile(draws, 0.975, na.rm = TRUE)),
                                fb_rate = fb_rate, stringsAsFactors = FALSE))
    })
    names(boot_list) <- cd8_cols
    bootstrap <- list(target = boot_target, draws = lapply(boot_list, `[[`, "draws"),
                      summary = do.call(rbind, lapply(boot_list, `[[`, "summary")))
    rownames(bootstrap$summary) <- NULL
  }

  ## ---- 8. Write back to meta.data ----------------------------------------
  if (add_meta) {
    for (cc in cd8_cols) {
      obj@meta.data[[paste0(prefix, "_diff_", cc)]]      <- NA_real_
      obj@meta.data[[paste0(prefix, "_label_", cc)]]     <- NA_character_
      obj@meta.data[[paste0(prefix, "_cd8scale_", cc)]]  <- NA_real_

      obj@meta.data[cells, paste0(prefix, "_diff_", cc)]     <- per[[cc]]$diff
      obj@meta.data[cells, paste0(prefix, "_label_", cc)]    <- per[[cc]]$label
      obj@meta.data[cells, paste0(prefix, "_cd8scale_", cc)] <- per[[cc]]$cd8_scale
    }
    nkcol <- paste0(prefix, "_nkscale")
    obj@meta.data[[nkcol]] <- NA_real_
    obj@meta.data[cells, nkcol] <- nk_use

    ccol <- paste0(prefix, "_consensus")
    obj@meta.data[[ccol]] <- NA_character_
    obj@meta.data[cells, ccol] <- unname(consensus)
  }

  ## ---- 9. Plotting (Violins + Matrices/Tables) ---------------------------
  out_plots <- list()
  if (make_plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("tidyr", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
      warning("Packages 'ggplot2', 'dplyr', and 'tidyr' are required for plots. Skipping plotting.", call. = FALSE)
    } else {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      if (!is.null(out_dir) && !dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

      # --- A. Variant Violin Plots ---
      for (cc in cd8_cols) {
        plot_data <- data.frame(
          cell_id          = cells,
          CD8_Score_Scaled = per[[cc]]$cd8_scale,
          NK_Score_Scaled  = nk_use,
          sc_Class         = if (plot_group_by == "consensus") unname(consensus) else per[[cc]]$label,
          stringsAsFactors = FALSE
        )

        plot_data <- plot_data[!is.na(plot_data$sc_Class), ]

        set.seed(seed)
        plot_data_sampled <- plot_data |>
          dplyr::group_by(sc_Class) |>
          dplyr::slice_sample(n = plot_n, replace = FALSE) |>
          dplyr::ungroup()

        plot_data_long <- plot_data_sampled |>
          tidyr::pivot_longer(cols = c("CD8_Score_Scaled", "NK_Score_Scaled"),
                              names_to = "Score_Type", values_to = "Score_Value")

        p_violin <- ggplot2::ggplot(plot_data_long, ggplot2::aes(x = Score_Type, y = Score_Value, color = sc_Class)) +
          ggplot2::geom_violin(ggplot2::aes(group = Score_Type, fill = sc_Class), alpha = 0.2, color = "gray30", linewidth = 0.5) +
          ggplot2::geom_line(ggplot2::aes(group = cell_id), alpha = 0.3, linewidth = 0.25) +
          ggplot2::theme_bw() +
          ggplot2::labs(
            title = paste("CD8 vs NK Scaled Scores:", cc),
            subtitle = sprintf("Assigned by %s | up to %d cells/class | Method: %s", plot_group_by, plot_n, toupper(method)),
            x = "Module", y = "Scaled UCell Score", color = "Assigned Class", fill = "Assigned Class"
          ) +
          ggplot2::facet_wrap(~ sc_Class) +
          ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold"), panel.grid.minor = ggplot2::element_blank())

        out_plots[[paste0("Violin_", cc)]] <- p_violin

        if (!is.null(out_dir)) {
          file_name <- paste0("CD8vNK_Violin_", method, "_", cc, "_", timestamp, ".png")
          ggplot2::ggsave(filename = file.path(out_dir, file_name), plot = p_violin, width = 12, height = 6)
        }
      }

      # --- B. Agreement Matrix Heatmap ---
      ag_df <- as.data.frame(as.table(agreement))
      colnames(ag_df) <- c("Variant_1", "Variant_2", "Agreement")
      p_ag <- ggplot2::ggplot(ag_df, ggplot2::aes(x = Variant_1, y = Variant_2, fill = Agreement)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = round(Agreement, 3)), color = "black") +
        ggplot2::scale_fill_gradient(low = "white", high = "dodgerblue", na.value = "grey90") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Pairwise Agreement Matrix | Method:", toupper(method)), x = "", y = "") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      out_plots[["Agreement_Matrix"]] <- p_ag
      if (!is.null(out_dir)) ggplot2::ggsave(filename = file.path(out_dir, paste0("CD8vNK_Agreement_", method, "_", timestamp, ".png")), plot = p_ag, width = 7, height = 6)

      # --- C. Primary vs Confirm Heatmap ---
      pvc_df <- as.data.frame(primary_vs_confirm)
      # ensure names match if missing
      if(ncol(pvc_df) == 3) colnames(pvc_df) <- c("Primary", "Confirm", "Freq")
      p_pvc <- ggplot2::ggplot(pvc_df, ggplot2::aes(x = Confirm, y = Primary, fill = Freq)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = Freq), color = "black", fontface = "bold") +
        ggplot2::scale_fill_gradient(low = "white", high = "firebrick") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = paste("Primary vs Confirm | Method:", toupper(method)),
                      subtitle = paste("Primary:", primary, "| Confirm:", confirm),
                      x = "Confirm Label", y = "Primary Label", fill = "Cell Count")

      out_plots[["Primary_vs_Confirm"]] <- p_pvc
      if (!is.null(out_dir)) ggplot2::ggsave(filename = file.path(out_dir, paste0("CD8vNK_PrimaryVsConfirm_", method, "_", timestamp, ".png")), plot = p_pvc, width = 6, height = 5)

      # --- D. Consensus Distribution Barplot ---
      cons_df <- as.data.frame(table(consensus, useNA = "always"))
      cons_df$consensus <- as.character(cons_df$consensus)
      cons_df$consensus[is.na(cons_df$consensus)] <- "Unassigned/NA"

      p_cons <- ggplot2::ggplot(cons_df, ggplot2::aes(x = consensus, y = Freq, fill = consensus)) +
        ggplot2::geom_col(color = "black", alpha = 0.8) +
        ggplot2::geom_text(ggplot2::aes(label = Freq), vjust = -0.5, fontface = "bold") +
        ggplot2::theme_bw() +
        ggplot2::labs(title = paste("Consensus Label Distribution | Method:", toupper(method)), x = "Consensus Label", y = "Total Cells") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")

      out_plots[["Consensus_Barplot"]] <- p_cons
      if (!is.null(out_dir)) ggplot2::ggsave(filename = file.path(out_dir, paste0("CD8vNK_Consensus_", method, "_", timestamp, ".png")), plot = p_cons, width = 7, height = 6)

      # --- E. Bootstrap Summary Forest Plot (if applicable) ---
      if (n_boot > 0 && !is.null(bootstrap)) {
        bs_df <- bootstrap$summary
        p_boot <- ggplot2::ggplot(bs_df, ggplot2::aes(x = variant, y = mean, ymin = ci_lo, ymax = ci_hi, color = variant)) +
          ggplot2::geom_pointrange(size = 1) +
          ggplot2::geom_point(ggplot2::aes(y = observed), color = "black", shape = 4, size = 3, stroke = 1.5) +
          ggplot2::theme_bw() +
          ggplot2::labs(title = paste("Bootstrap Summary (", boot_target, ") | Method:", toupper(method)),
                        subtitle = "Crosses (x) represent the original observed value on the full dataset.",
                        y = boot_target, x = "Variant") +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")

        out_plots[["Bootstrap_ForestPlot"]] <- p_boot
        if (!is.null(out_dir)) ggplot2::ggsave(filename = file.path(out_dir, paste0("CD8vNK_Bootstrap_", method, "_", timestamp, ".png")), plot = p_boot, width = 8, height = 5)
      }
    }
  }

  ## ---- 10. Return ---------------------------------------------------------
  list(obj        = obj,
       per        = per,
       labels     = labels,
       consensus  = consensus,
       agreement  = agreement,
       n_disagree = n_disagree,
       cutoffs    = cut_tab,
       bootstrap  = bootstrap,
       plots      = out_plots,
       cells_used = cells,
       dropped    = dropped,
       params     = list(cd8_cols     = cd8_cols,
                         nk_col       = nk_col,
                         method       = method,
                         fixed_cut    = fixed_cut,
                         zscore       = zscore,
                         primary      = primary,
                         confirm      = confirm,
                         min_features = min_features,
                         feature_col  = feature_col,
                         n_boot       = n_boot,
                         boot_target  = boot_target,
                         prefix       = prefix,
                         seed         = seed))
}
