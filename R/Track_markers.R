# Track markers down the resolution tree and suggest a label for every cluster.
#
# track_cluster_markers() takes the tree and sibling markers of
# find_hierarchical_markers() (or runs it), follows each marker from the cluster where
# it first appears down its branch, profiles the tracked genes at every cluster, and
# builds a suggested path label per cluster from the top marker of each split above it.

#' Track Markers Down the Resolution Tree to Annotate Clusters
#'
#' Follows the markers found by \code{\link{find_hierarchical_markers}()} from low to
#' high resolution, to help name clusters. For every gene it shows where it first
#' marks a cluster, how far down that branch it stays detected and where it is lost.
#' For every cluster it suggests a label made of the top marker at each split on the
#' way from the root, e.g. \code{"CST3+ > FCGR3A+"}, which you translate into a cell
#' type. Labels are suggestions only: you make the final call.
#'
#' @details
#' \strong{Tracked genes.} The top \code{n_track} markers (by AUC against the
#' siblings) of every child at every split.
#'
#' \strong{Profile.} For each tracked gene at every cluster of every resolution: mean
#' expression, share of cells detected (\code{pct}) and outside the cluster
#' (\code{pct_rest}), log2 fold change and AUC against all other cells, and AUC against
#' its siblings (the other children of its parent; all their cells). \code{global} is
#' \code{auc_vs_all >= min_global_auc}: the gene is high in the cluster compared with the
#' whole object, not only compared with its siblings ("relative").
#'
#' \strong{Gene paths.} A gene's path starts at a cluster where it is a marker and no
#' cluster above it on the branch is marked by it. Below the start, a cluster is "on"
#' when the gene is detected in at least \code{min_pct} of its cells (the detection
#' filter of \code{find_hierarchical_markers()}); \code{stops_at} lists the first
#' clusters where it is no longer on.
#'
#' \strong{Labels.} For each split on the way from the root, the label gene is the
#' highest-ranked global marker among the top \code{n_track} of that child (if none is
#' global, the top marker, flagged relative). A label is \code{"ambiguous"} when, at any
#' split on its path, there is no marker (or the child was too small to test), the label gene's AUC against the siblings is
#' below \code{min_label_auc}, the label gene is only relative, it is lower in the
#' cluster than in the rest of the object (AUC against all < 0.5), or a cluster on the
#' path has low purity. \code{reason} says which.
#'
#' Statistics use the same cells and clusters, so they describe the data rather than
#' prove the labels; check the labels against known biology.
#'
#' @param seurat_obj A \code{Seurat} object (the one used for \code{hier}).
#' @param hier Optional result of \code{\link{find_hierarchical_markers}()} on
#'   \code{seurat_obj}. If \code{NULL} (default), it is run with \code{res_cols} or
#'   \code{pipeline} and \code{...}.
#' @param res_cols,pipeline Resolution columns or integration method, as in
#'   \code{find_hierarchical_markers()}. Used only when \code{hier} is \code{NULL}.
#' @param n_track Integer. Top markers per child and split to track. Default 10.
#' @param min_label_auc Numeric. Minimum AUC against the siblings for a label gene.
#'   Default 0.8.
#' @param min_global_auc Numeric. Minimum AUC against all other cells for a gene to be
#'   \code{global}. Default 0.7.
#' @param n_own Integer. Own markers listed per cluster. Default 5.
#' @param sample_name Character. Used in plot titles and file names. Default "Sample".
#' @param output_dir Character or \code{NULL}. If given, the tables (CSV) and plots
#'   (JPEG, sized from their content) are saved there. Default \code{NULL}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#' @param ... Passed to \code{find_hierarchical_markers()} when \code{hier} is
#'   \code{NULL} (e.g. \code{assay}, \code{min_cells}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{annotation}: one row per cluster: \code{node}, \code{layer},
#'     \code{resolution}, \code{cluster}, \code{n_cells}, \code{is_final},
#'     \code{path_label}, \code{status} (\code{"suggested"} / \code{"ambiguous"}),
#'     \code{reason}, \code{confidence} (lowest AUC vs siblings of the label genes on
#'     the path), \code{min_auc_vs_all} (lowest AUC vs all of the label genes),
#'     \code{lineage_markers} (lineage markers inherited from the splits above),
#'     \code{own_markers} (top global markers at the cluster's own split),
#'     \code{relative_markers} (top markers there that are only relative),
#'     \code{lowest_purity} (on the path).
#'   \item \code{gene_paths}: one row per gene and start: \code{gene}, \code{start_node},
#'     \code{start_resolution}, \code{start_split}, \code{kind} (the marker pattern at
#'     the start), \code{auc_at_start}, \code{global_at_start}, \code{marked_again_below},
#'     \code{n_final}, \code{n_final_on}, \code{final_on}, \code{stops_at}.
#'   \item \code{node_profile}: one row per cluster and tracked gene: \code{node},
#'     \code{layer}, \code{resolution}, \code{gene}, \code{mean}, \code{pct},
#'     \code{pct_rest}, \code{log2fc_vs_all}, \code{auc_vs_all}, \code{global},
#'     \code{auc_vs_siblings}, \code{is_marker}.
#'   \item \code{plots}: \code{tree} (tree with the label gene of every split),
#'     \code{heatmap} (AUC against all cells of the tracked genes at every cluster,
#'     markers dotted) and \code{dotplot} (label genes across the final clusters in tree
#'     order).
#'   \item \code{hier}: the \code{find_hierarchical_markers()} result used.
#'   \item \code{parameters}.
#' }
#'
#' @seealso \code{\link{find_hierarchical_markers}()}
#' @export
#'
#' @examples
#' \dontrun{
#' hier <- find_hierarchical_markers(my_seurat, pipeline = "log.rpca")
#' ann <- track_cluster_markers(my_seurat, hier = hier, output_dir = "Track_markers")
#' ann$annotation[, c("node", "path_label", "status", "reason")]
#' subset(ann$gene_paths, gene == "GZMK")
#' ann$plots$heatmap
#' }
track_cluster_markers <- function(seurat_obj,
                                  hier = NULL,
                                  res_cols = NULL,
                                  pipeline = NULL,
                                  n_track = 10,
                                  min_label_auc = 0.8,
                                  min_global_auc = 0.7,
                                  n_own = 5,
                                  sample_name = "Sample",
                                  output_dir = NULL,
                                  verbose = TRUE,
                                  ...) {
  if (!inherits(seurat_obj, "Seurat")) stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  say <- function(...) if (verbose) message(...)

  if (is.null(hier)) {
    hier <- find_hierarchical_markers(seurat_obj, res_cols = res_cols, pipeline = pipeline,
                                      sample_name = sample_name, verbose = verbose,
                                      min_global_auc = min_global_auc, ...)
  } else if (!is.null(res_cols) || !is.null(pipeline)) {
    warning("`hier` is given, so `res_cols` / `pipeline` are ignored.", call. = FALSE)
  }
  tree <- hier$tree
  markers <- hier$markers
  assay <- hier$parameters$assay
  min_pct <- hier$parameters$min_pct

  # Cluster membership of every cell at every resolution, rebuilt from the metadata
  md <- seurat_obj@meta.data[, hier$res_cols, drop = FALSE]
  complete <- stats::complete.cases(md)
  seurat_obj <- seurat_obj[, complete]
  built <- .build_resolution_tree(md[complete, , drop = FALSE], hier$res_cols, hier$parameters$min_purity)
  if (!identical(built$tree$node, tree$node) || !identical(built$tree$n_cells, tree$n_cells)) {
    stop("`hier` does not match `seurat_obj`: run find_hierarchical_markers() on this object.", call. = FALSE)
  }
  seurat_obj <- .join_data_layers(seurat_obj, assay)
  expr <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = "data")

  tracked <- unique(markers$gene[markers$rank <= n_track])
  say("Tracking ", length(tracked), " genes over ", nrow(tree) - 1, " clusters.")

  node_profile <- .profile_nodes(tree, built$cell_nodes, expr, tracked, markers, min_global_auc)
  gene_paths <- .gene_paths(tree, markers[markers$gene %in% tracked, , drop = FALSE], node_profile, min_pct)
  annotation <- .annotate_nodes(tree, markers, hier$splits, node_profile, n_track, n_own, min_label_auc,
                                min_global_auc, hier$parameters$min_purity)
  n_ambiguous <- sum(annotation$status == "ambiguous")
  say("Suggested labels: ", sum(annotation$status == "suggested"), " clusters; ambiguous: ", n_ambiguous, ".")

  plots <- list(
    tree    = .plot_annotation_tree(tree, annotation, markers, n_track, sample_name),
    heatmap = .plot_tracking_heatmap(tree, node_profile, gene_paths, sample_name),
    dotplot = .plot_label_dotplot(seurat_obj, tree, annotation, hier$res_cols, assay, sample_name)
  )

  out <- list(annotation = annotation, gene_paths = gene_paths, node_profile = node_profile,
              plots = plots, hier = hier,
              parameters = list(n_track = n_track, min_label_auc = min_label_auc,
                                min_global_auc = min_global_auc, n_own = n_own, min_pct = min_pct))
  if (!is.null(output_dir)) .save_tracking(out, output_dir, sample_name, say)
  out
}

# -------------------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------------------

#' Statistics of every tracked gene at every cluster (all resolutions)
#'
#' Each gene is ranked once over all cells (AUC vs all) and once within every family
#' of siblings (AUC vs siblings); AUCs then come from rank sums.
#' @keywords internal
#' @noRd
.profile_nodes <- function(tree, cell_nodes, expr, genes, markers, min_global_auc) {
  nodes <- tree[tree$node != "root", , drop = FALSE]
  positions <- lapply(seq_len(nrow(nodes)), function(i) which(cell_nodes[, nodes$layer[i] + 1] == nodes$node[i]))
  names(positions) <- nodes$node
  families <- split(nodes$node, nodes$parent)
  families <- families[lengths(families) >= 2]          # a single child has no siblings
  n_all <- ncol(expr)
  marker_pairs <- paste(markers$child, markers$gene)

  rows <- lapply(genes, function(gene) {
    x <- as.numeric(expr[gene, ])
    ranks_all <- rank(x)
    linear_all <- expm1(x)
    total_linear <- sum(linear_all)
    total_detected <- sum(x > 0)

    stats <- t(vapply(positions, function(pos) {
      n_in <- length(pos)
      n_out <- n_all - n_in
      mean_in_linear <- mean(linear_all[pos])
      mean_out_linear <- if (n_out > 0) (total_linear - sum(linear_all[pos])) / n_out else NA_real_
      detected_in <- sum(x[pos] > 0)
      c(mean = mean(x[pos]), pct = detected_in / n_in,
        pct_rest = if (n_out > 0) (total_detected - detected_in) / n_out else NA_real_,
        log2fc_vs_all = log2(mean_in_linear + 1) - log2(mean_out_linear + 1),   # Seurat-style
        auc_vs_all = .auc_from_ranks(ranks_all, pos))
    }, numeric(5)))

    auc_siblings <- stats::setNames(rep(NA_real_, nrow(nodes)), nodes$node)
    for (family in families) {
      family_pos <- unlist(positions[family], use.names = FALSE)
      family_ranks <- rank(x[family_pos])
      offset <- 0
      for (child in family) {
        n_child <- length(positions[[child]])
        in_child <- offset + seq_len(n_child)
        n_other <- length(family_pos) - n_child
        auc_siblings[[child]] <- if (n_other > 0) {
          (sum(family_ranks[in_child]) - n_child * (n_child + 1) / 2) / (n_child * n_other)
        } else NA_real_
        offset <- offset + n_child
      }
    }

    data.frame(node = nodes$node, layer = nodes$layer, resolution = nodes$resolution, gene = gene,
               stats, global = stats[, "auc_vs_all"] >= min_global_auc,
               auc_vs_siblings = unname(auc_siblings[nodes$node]),
               is_marker = paste(nodes$node, gene) %in% marker_pairs,
               stringsAsFactors = FALSE, row.names = NULL)
  })
  if (length(rows) == 0) {
    return(data.frame(node = character(), layer = integer(), resolution = numeric(), gene = character(),
                      mean = numeric(), pct = numeric(), pct_rest = numeric(), log2fc_vs_all = numeric(),
                      auc_vs_all = numeric(), global = logical(), auc_vs_siblings = numeric(),
                      is_marker = logical(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' All nodes below a node (not including it)
#' @keywords internal
#' @noRd
.descendants <- function(node, children_of) {
  kids <- children_of[[node]]
  if (length(kids) == 0) return(character())
  c(kids, unlist(lapply(kids, .descendants, children_of = children_of), use.names = FALSE))
}

#' Ancestors of a node from its parent up to (not including) root
#' @keywords internal
#' @noRd
.ancestors <- function(node, parent_of) {
  out <- character()
  current <- parent_of[[node]]
  while (!is.na(current) && current != "root") {
    out <- c(out, current)
    current <- parent_of[[current]]
  }
  out
}

#' Where each tracked gene starts to mark a branch, and where it is lost below
#' @keywords internal
#' @noRd
.gene_paths <- function(tree, markers, node_profile, min_pct) {
  empty <- data.frame(gene = character(), start_node = character(), start_resolution = numeric(),
                      start_split = character(), kind = character(), auc_at_start = numeric(),
                      global_at_start = logical(), marked_again_below = character(), n_final = integer(),
                      n_final_on = integer(), final_on = character(), stops_at = character(),
                      stringsAsFactors = FALSE)
  if (nrow(markers) == 0) return(empty)
  children_of <- split(tree$node, tree$parent)
  parent_of <- stats::setNames(tree$parent, tree$node)
  is_final <- stats::setNames(tree$n_children == 0, tree$node)
  profile_key <- paste(node_profile$node, node_profile$gene)

  rows <- lapply(split(markers, markers$gene), function(m) {
    gene <- m$gene[1]
    pct_of <- stats::setNames(node_profile$pct[node_profile$gene == gene], node_profile$node[node_profile$gene == gene])
    marked <- unique(m$child)
    starts <- marked[vapply(marked, function(n) !any(.ancestors(n, parent_of) %in% marked), logical(1))]

    do.call(rbind, lapply(starts, function(start) {
      row <- m[m$child == start, , drop = FALSE][1, ]
      below <- .descendants(start, children_of)
      branch <- c(start, below)
      finals <- branch[is_final[branch]]
      on <- pct_of[branch] >= min_pct
      first_off <- below[!on[below] & on[parent_of[below]]]
      data.frame(
        gene = gene, start_node = start, start_resolution = tree$resolution[tree$node == start],
        start_split = row$split, kind = row$pattern, auc_at_start = row$auc,
        global_at_start = node_profile$global[match(paste(start, gene), profile_key)],
        marked_again_below = paste(intersect(below, marked), collapse = ", "),
        n_final = length(finals), n_final_on = sum(on[finals]),
        final_on = paste(finals[on[finals]], collapse = ", "),
        stops_at = paste(first_off, collapse = ", "),
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  out <- out[order(match(out$start_node, tree$node), out$auc_at_start * -1), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' The splits on the way from the root to a node: (split, child) pairs, top first
#'
#' The child of each split is the start of the unchanged chain the node belongs to.
#' @keywords internal
#' @noRd
.path_steps <- function(node, tops, parent_of) {
  steps <- list()
  current <- node
  repeat {
    top <- tops[[current]]
    if (top == "root") break
    split_node <- parent_of[[top]]
    steps <- c(list(c(split = split_node, child = top)), steps)
    current <- split_node
  }
  steps
}

#' The label gene of one split for one child: the best global marker among the top
#' n_track, otherwise the top marker (relative)
#' @keywords internal
#' @noRd
.label_gene <- function(step, markers, auc_vs_all_of, n_track, min_global_auc) {
  m <- markers[markers$split == step[["split"]] & markers$child == step[["child"]] & markers$rank <= n_track, ,
               drop = FALSE]
  m <- m[order(m$rank), , drop = FALSE]
  if (nrow(m) == 0) {
    return(list(gene = NA_character_, auc = NA_real_, auc_vs_all = NA_real_, global = NA,
                markers = cbind(m, global_here = logical())))
  }
  vs_all <- auc_vs_all_of[paste(step[["child"]], m$gene)]
  global <- !is.na(vs_all) & vs_all >= min_global_auc
  pick <- if (any(global)) which(global)[1] else 1
  list(gene = m$gene[pick], auc = m$auc[pick], auc_vs_all = unname(vs_all[pick]), global = global[pick],
       markers = cbind(m, global_here = global))
}

#' Suggested path label, evidence and status for every cluster
#' @keywords internal
#' @noRd
.annotate_nodes <- function(tree, markers, splits, node_profile, n_track, n_own, min_label_auc,
                            min_global_auc, min_purity) {
  tops <- .chain_tops(tree)
  parent_of <- stats::setNames(tree$parent, tree$node)
  purity_of <- stats::setNames(tree$purity, tree$node)
  auc_vs_all_of <- stats::setNames(node_profile$auc_vs_all, paste(node_profile$node, node_profile$gene))
  # Children left out of a split for having too few cells (find_hierarchical_markers' min_cells)
  skipped_pairs <- unlist(lapply(seq_len(nrow(splits)), function(i) {
    kids <- strsplit(splits$skipped_children[i], ", ", fixed = TRUE)[[1]]
    if (length(kids) == 0) character() else paste(splits$split[i], kids)
  }))

  rows <- lapply(tree$node[tree$node != "root"], function(node) {
    steps <- .path_steps(node, tops, parent_of)
    picks <- lapply(steps, .label_gene, markers = markers, auc_vs_all_of = auc_vs_all_of, n_track = n_track,
                    min_global_auc = min_global_auc)
    genes <- vapply(picks, `[[`, character(1), "gene")
    aucs <- vapply(picks, `[[`, numeric(1), "auc")
    vs_all <- vapply(picks, `[[`, numeric(1), "auc_vs_all")
    global <- vapply(picks, function(p) isTRUE(p$global), logical(1))
    step_splits <- vapply(steps, `[[`, character(1), "split")
    too_small <- vapply(steps, function(st) paste(st[["split"]], st[["child"]]) %in% skipped_pairs, logical(1))
    no_markers <- is.na(genes) & !too_small

    path_nodes <- c(node, .ancestors(node, parent_of))
    lowest_purity <- min(purity_of[path_nodes])
    reasons <- c(
      if (any(too_small)) paste0("too few cells to test at ", paste(step_splits[too_small], collapse = ", ")),
      if (any(no_markers)) paste0("no markers at ", paste(step_splits[no_markers], collapse = ", ")),
      if (any(!is.na(aucs) & aucs < min_label_auc))
        paste0("weak marker at ", paste(step_splits[!is.na(aucs) & aucs < min_label_auc], collapse = ", ")),
      if (any(!is.na(genes) & !global))
        paste0("only relative markers at ", paste(step_splits[!is.na(genes) & !global], collapse = ", ")),
      if (any(!is.na(vs_all) & vs_all < 0.5))
        paste0("lower than the rest: ", paste(genes[!is.na(vs_all) & vs_all < 0.5], collapse = ", ")),
      if (lowest_purity < min_purity)
        paste0("low purity: ", paste(path_nodes[purity_of[path_nodes] < min_purity], collapse = ", ")),
      if (length(steps) == 0) "no split above this cluster"
    )

    own <- if (length(picks) > 0) picks[[length(picks)]]$markers else NULL
    lineage <- unlist(lapply(picks[-length(picks)], function(p) {
      if (is.null(p$markers) || nrow(p$markers) == 0) return(character())
      utils::head(p$markers$gene[p$markers$pattern == "lineage"], 2)
    }), use.names = FALSE)

    data.frame(
      node = node, layer = tree$layer[tree$node == node], resolution = tree$resolution[tree$node == node],
      cluster = tree$cluster[tree$node == node], n_cells = tree$n_cells[tree$node == node],
      is_final = tree$n_children[tree$node == node] == 0,
      path_label = paste(ifelse(is.na(genes), "?", paste0(genes, "+")), collapse = " > "),
      status = if (length(reasons) == 0) "suggested" else "ambiguous",
      reason = paste(reasons, collapse = "; "),
      confidence = if (all(is.na(aucs))) NA_real_ else min(aucs, na.rm = TRUE),
      min_auc_vs_all = if (all(is.na(vs_all))) NA_real_ else min(vs_all, na.rm = TRUE),
      lineage_markers = paste(unique(lineage), collapse = ", "),
      own_markers = if (is.null(own)) "" else paste(utils::head(own$gene[own$global_here], n_own), collapse = ", "),
      relative_markers = if (is.null(own)) "" else paste(utils::head(own$gene[!own$global_here], n_own), collapse = ", "),
      lowest_purity = lowest_purity,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Resolution tree with the label gene of every split (red: relative or weak)
#' @keywords internal
#' @noRd
.plot_annotation_tree <- function(tree, annotation, markers, n_track, sample_name) {
  x <- .tree_x_positions(tree)
  nodes <- data.frame(node = tree$node, x = x[tree$node], y = -tree$layer, n_cells = tree$n_cells,
                      label = ifelse(tree$node == "root", "all cells", tree$cluster), stringsAsFactors = FALSE)
  child_rows <- tree$node != "root"
  edges <- data.frame(x = x[tree$parent[child_rows]], xend = x[tree$node[child_rows]],
                      y = -tree$layer[child_rows] + 1, yend = -tree$layer[child_rows])

  # The label gene of each chain start is the last part of its path label
  starts <- annotation[annotation$node %in% unique(.chain_tops(tree)), , drop = FALSE]
  gene_labels <- data.frame(
    x = x[starts$node], y = -starts$layer,
    label = sub(".* > ", "", starts$path_label),
    status = starts$status, stringsAsFactors = FALSE
  )
  layer_labels <- unique(data.frame(y = -tree$layer[child_rows],
                                    label = ifelse(is.na(tree$resolution[child_rows]), tree$res_col[child_rows],
                                                   paste("res", tree$resolution[child_rows]))))

  ggplot() +
    geom_segment(data = edges, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y), colour = "grey60") +
    geom_segment(data = edges, aes(x = .data$xend, xend = .data$xend, y = .data$y, yend = .data$yend),
                 colour = "grey60") +
    geom_point(data = nodes, aes(x = .data$x, y = .data$y, size = .data$n_cells), colour = "grey20") +
    geom_text(data = nodes, aes(x = .data$x, y = .data$y, label = .data$label), vjust = -1.5, size = 3) +
    geom_text(data = gene_labels, aes(x = .data$x, y = .data$y, label = .data$label, colour = .data$status),
              vjust = 1.8, size = 3, fontface = "italic") +
    scale_colour_manual(values = c(suggested = "steelblue4", ambiguous = "red3"), name = "Label") +
    guides(colour = guide_legend(override.aes = list(label = "gene+", size = 3.5))) +
    scale_size_area(max_size = 6, name = "Cells") +
    scale_y_continuous(breaks = c(0, layer_labels$y), labels = c("root", layer_labels$label),
                       expand = expansion(mult = 0.15)) +
    labs(title = paste(sample_name, "- label gene of every split"), x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank())
}

#' Heatmap of tracked genes x clusters: colour = AUC against all other cells, dot = marker
#'
#' Genes are ordered by where their path starts; clusters by resolution, then tree order.
#' @keywords internal
#' @noRd
.plot_tracking_heatmap <- function(tree, node_profile, gene_paths, sample_name) {
  if (nrow(node_profile) == 0) return(NULL)
  x <- .tree_x_positions(tree)
  node_order <- tree$node[tree$node != "root"]
  node_order <- node_order[order(tree$layer[match(node_order, tree$node)], x[node_order])]
  gene_order <- unique(c(gene_paths$gene, node_profile$gene))

  df <- node_profile
  df$node <- factor(df$node, levels = node_order)
  df$gene <- factor(df$gene, levels = rev(gene_order))
  df$resolution_label <- factor(paste("res", df$resolution),
                                levels = unique(paste("res", tree$resolution[match(node_order, tree$node)])))
  df$cluster <- tree$cluster[match(as.character(df$node), tree$node)]

  ggplot(df, aes(x = .data$node, y = .data$gene)) +
    geom_tile(aes(fill = .data$auc_vs_all)) +
    geom_point(data = df[df$is_marker, , drop = FALSE], size = 0.8) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "red3", midpoint = 0.5, limits = c(0, 1),
                         name = "AUC vs\nall cells") +
    scale_x_discrete(labels = function(n) tree$cluster[match(n, tree$node)]) +
    facet_grid(. ~ resolution_label, scales = "free_x", space = "free_x") +
    labs(title = paste(sample_name, "- tracked markers across resolutions"),
         subtitle = "Dot: marker of that cluster against its siblings", x = "Cluster", y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 7, face = "italic"))
}

#' Dot plot of all label genes across the final clusters, in tree order
#' @keywords internal
#' @noRd
.plot_label_dotplot <- function(seurat_obj, tree, annotation, res_cols, assay, sample_name) {
  genes <- unique(unlist(strsplit(gsub("\\+", "", annotation$path_label[annotation$is_final]), " > ")))
  genes <- setdiff(genes, "?")
  if (length(genes) == 0) return(NULL)
  top_col <- res_cols[length(res_cols)]
  x <- .tree_x_positions(tree)
  top_layer <- max(tree$layer)
  top_nodes <- tree$node[tree$layer == top_layer]
  cluster_order <- tree$cluster[match(top_nodes[order(x[top_nodes])], tree$node)]
  seurat_obj$.tracked_cluster <- factor(as.character(seurat_obj@meta.data[[top_col]]), levels = cluster_order)
  suppressWarnings(suppressMessages(
    Seurat::DotPlot(seurat_obj, features = genes, group.by = ".tracked_cluster", assay = assay)
  )) +
    Seurat::RotatedAxis() +
    labs(title = paste(sample_name, "- label genes"), subtitle = paste(top_col, "(tree order)"), y = NULL)
}

#' Save the tables (CSV) and plots (JPEG, sized from their content)
#' @keywords internal
#' @noRd
.save_tracking <- function(out, output_dir, sample_name, say) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(output_dir, paste0(sanitize(sample_name), "_TrackMarkers_",
                                       format(Sys.time(), "%Y%m%d_%H%M%S")))
  for (table_name in c("annotation", "gene_paths", "node_profile")) {
    utils::write.csv(out[[table_name]], paste0(stem, "_", table_name, ".csv"), row.names = FALSE)
  }
  sizes <- .tracking_image_sizes(out)
  for (plot_name in names(sizes)) {
    if (is.null(out$plots[[plot_name]])) next
    ggsave(paste0(stem, "_", plot_name, ".jpg"), out$plots[[plot_name]],
           width = sizes[[plot_name]][["width"]], height = sizes[[plot_name]][["height"]],
           dpi = 300, limitsize = FALSE, bg = "white")
  }
  say("Saved tracked markers to: ", output_dir)
}

#' Width and height (inches) of each saved tracking image
#' @keywords internal
#' @noRd
.tracking_image_sizes <- function(out) {
  tree <- out$hier$tree
  n_layers <- max(tree$layer)
  n_final <- sum(tree$n_children == 0)
  n_nodes <- nrow(tree) - 1
  n_genes <- length(unique(out$node_profile$gene))
  n_label_genes <- length(setdiff(unique(unlist(strsplit(
    gsub("\\+", "", out$annotation$path_label[out$annotation$is_final]), " > "))), "?"))
  n_top <- sum(tree$layer == n_layers)
  list(
    # gene labels sit side by side, so each final cluster gets more room than in the key tree
    tree    = .bounded_size(1.1 * n_final, 2 + HIER_IMAGE$per_layer * n_layers, 8, 4),
    heatmap = .bounded_size(3 + 0.25 * n_nodes, 2 + 0.15 * n_genes, 6, 4),
    dotplot = .bounded_size(3 + HIER_IMAGE$per_gene * n_label_genes, 2 + HIER_IMAGE$per_dot_row * n_top, 6, 4)
  )
}
