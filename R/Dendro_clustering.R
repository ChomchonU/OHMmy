#' Cut the Dendrograms of \code{plot_dot_dendro_multi()} into Clusters
#'
#' Rebuilds the exact hierarchical clustering behind the dendrograms drawn by
#' \code{\link{plot_dot_dendro_multi}()} and cuts it into clusters of groups (rows),
#' genes (columns), or both. Because natural clusters often sit at different heights
#' of a dendrogram, the default methods cut each branch at its own height instead of
#' using a single horizontal cut, so both the number of clusters and the cut heights
#' are derived from the tree:
#' \itemize{
#'   \item \code{"inconsistency"} (default): a branch is kept as one cluster while none of
#'     its merges is unusually high compared with the merges directly below it
#'     (inconsistency coefficient, depth 2). Each branch is therefore cut relative to its
#'     own spread, which suits clusters at clearly different heights.
#'   \item \code{"dynamic"}: the hybrid dynamic tree cut of \pkg{dynamicTreeCut}
#'     (Langfelder et al., 2008), which uses both the tree shape and the distance matrix.
#'     It works best when clusters have broadly similar heights.
#'   \item \code{"fixed"}: an ordinary horizontal cut at \code{k} clusters or height \code{h}.
#' }
#'
#' @details
#' The trees are taken from the \code{"row_hclust"} and \code{"col_hclust"} attributes that
#' \code{plot_dot_dendro_multi()} attaches to its output. For older results without these
#' attributes, the same trees are rebuilt from the \code{avg.exp.scaled} values
#' (Euclidean distance, complete linkage).
#'
#' With \code{deepSplit = "auto"} (dynamic method) or \code{cutoff = "auto"}
#' (inconsistency method), every candidate setting is tried and the one with the highest
#' average silhouette width is chosen. Unassigned leaves count as a silhouette of 0, so
#' leaving many leaves out is not rewarded. All candidates are reported in
#' \code{selection}.
#'
#' Note that \code{DotPlot()} z-scores expression separately within each metadata column,
#' so when several columns are combined, groups are compared by expression pattern rather
#' than on one absolute scale. This is the same basis as the plotted dendrograms.
#'
#' @param dot_df The data frame returned by \code{\link{plot_dot_dendro_multi}()}.
#' @param what Character. Cut the tree of the groups (\code{"rows"}, default), of the genes (\code{"genes"}), or both (\code{"both"}).
#' @param method Character. The cutting method: \code{"inconsistency"} (default), \code{"dynamic"} or \code{"fixed"}. See Description.
#' @param deepSplit Integer from 0 to 4, or \code{"auto"} (default). Dynamic method only: how readily branches are split into sub-clusters (higher gives more, smaller clusters). \code{"auto"} picks the value with the best average silhouette width.
#' @param cutoff Numeric or \code{"auto"} (default). Inconsistency method only: branches whose merges all have an inconsistency coefficient at or below \code{cutoff} form one cluster. \code{"auto"} picks the value with the best average silhouette width.
#' @param minClusterSize Integer. The smallest number of leaves that can form a cluster with the dynamic and inconsistency methods; smaller groups are left unassigned. Default is 2. Not used by the fixed method.
#' @param assign_all Logical. If \code{FALSE} (default), leaves that do not clearly belong to any cluster stay unassigned (\code{NA}). If \code{TRUE}, each of them is assigned to the cluster with the smallest average distance.
#' @param k,h Numeric. Fixed method only: the number of clusters (\code{k}) or the height at which to cut (\code{h}). Give exactly one.
#' @param seurat_obj Optional Seurat object. If given (together with \code{meta_col}), the row clusters are written back to the cells as a new metadata column.
#' @param meta_col Character. The metadata column (one of the \code{meta_cols} used in \code{plot_dot_dendro_multi()}) whose groups are mapped to cells. Each cell receives the cluster of its group's row; cells of unassigned groups get \code{NA}.
#' @param new_col Character. The name of the new metadata column. Default is \code{"dendro_cluster"}.
#' @param output_dir Character. Optional directory in which the cluster table (CSV) and the coloured dendrogram(s) (JPEG) are saved. Default is \code{NULL} (nothing is saved).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{clusters}: a data frame with one row per leaf: \code{type} (\code{"rows"} or \code{"genes"}), \code{label}, \code{id} and \code{var} (rows only), \code{cluster} (factor, \code{NA} = unassigned), \code{dendro_order} (position in the dendrogram), \code{cluster_height} (height of the merge joining the cluster's leaves) and \code{cut_height} (the height at which the cluster's branch was cut: midway between \code{cluster_height} and the merge joining it to the rest of the tree).
#'   \item \code{k}: a named integer vector with the number of clusters found per tree.
#'   \item \code{selection}: a named list with, per tree, the candidate settings and their scores (for \code{"auto"} settings), or \code{NULL}.
#'   \item \code{plot}: a named list of \code{ggplot} dendrograms coloured by cluster, with a dashed mark at each cluster's cut height; unassigned leaves are grey.
#'   \item \code{hclust}: a named list of the \code{hclust} trees that were cut.
#'   \item \code{parameters}: the method and the settings that were used.
#'   \item \code{seurat_obj}: the updated Seurat object (only when \code{seurat_obj} was given).
#' }
#'
#' @references Langfelder P, Zhang B, Horvath S (2008). Defining clusters from a hierarchical cluster tree: the Dynamic Tree Cut package for R. \emph{Bioinformatics} 24(5):719-720.
#'
#' @seealso \code{\link{plot_dot_dendro_multi}()}
#' @export
#'
#' @examples
#' \dontrun{
#' dot_df <- plot_dot_dendro_multi(seurat_obj, meta_cols = c("cell_type", "seurat_clusters"),
#'                                 feature_df = marker_df, prefix = "PBMC")
#'
#' # Cut the group dendrogram at variable heights
#' res <- cluster_dot_dendro(dot_df)
#' res$clusters
#' res$plot$rows
#'
#' # Genes and groups with the dynamic tree cut, and write the group clusters to the cells
#' res <- cluster_dot_dendro(dot_df, what = "both", method = "dynamic",
#'                           seurat_obj = seurat_obj, meta_col = "seurat_clusters")
#' table(res$seurat_obj$dendro_cluster)
#' }
cluster_dot_dendro <- function(dot_df,
                               what = c("rows", "genes", "both"),
                               method = c("inconsistency", "dynamic", "fixed"),
                               deepSplit = "auto",
                               cutoff = "auto",
                               minClusterSize = 2,
                               assign_all = FALSE,
                               k = NULL,
                               h = NULL,
                               seurat_obj = NULL,
                               meta_col = NULL,
                               new_col = "dendro_cluster",
                               output_dir = NULL) {
  what   <- match.arg(what)
  method <- match.arg(method)
  types  <- switch(what, rows = "rows", genes = "genes", both = c("rows", "genes"))

  # --- Validate arguments ------------------------------------------------------------
  if (method == "dynamic") .check_suggested("dynamicTreeCut", "cluster_dot_dendro(method = \"dynamic\")")
  if (method == "fixed" && (is.null(k) == is.null(h))) {
    stop("method = \"fixed\" needs exactly one of `k` or `h`.", call. = FALSE)
  }
  if (!identical(deepSplit, "auto") && !(is.numeric(deepSplit) && length(deepSplit) == 1 && deepSplit %in% 0:4)) {
    stop("`deepSplit` must be \"auto\" or a single integer from 0 to 4.", call. = FALSE)
  }
  if (!identical(cutoff, "auto") && !(is.numeric(cutoff) && length(cutoff) == 1)) {
    stop("`cutoff` must be \"auto\" or a single number.", call. = FALSE)
  }
  if (!is.null(seurat_obj)) {
    if (!"rows" %in% types) stop("Writing clusters to cells needs the row clusters: use what = \"rows\" or \"both\".", call. = FALSE)
    if (is.null(meta_col)) stop("`meta_col` is required when `seurat_obj` is given.", call. = FALSE)
  }

  trees <- .dot_dendro_trees(dot_df)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # --- Cut each requested tree ----------------------------------------------------------
  results <- lapply(types, function(tp) {
    hc <- trees$hclust[[tp]]
    dm <- as.matrix(trees$dist[[tp]])[hc$labels, hc$labels]
    if (length(hc$labels) < 3) stop("The ", tp, " dendrogram has fewer than 3 leaves; nothing to cluster.", call. = FALSE)

    cut <- switch(method,
                  dynamic       = .cut_dynamic(hc, dm, deepSplit, minClusterSize),
                  inconsistency = .cut_inconsistency(hc, dm, cutoff, minClusterSize),
                  fixed         = list(labels = if (!is.null(k)) stats::cutree(hc, k = k) else stats::cutree(hc, h = h),
                                       selection = NULL, setting = if (!is.null(k)) c(k = k) else c(h = h)))

    labels <- cut$labels
    if (assign_all) labels <- .assign_leftovers(labels, dm)
    labels <- .relabel_in_dendro_order(labels, hc)

    heights <- .cluster_heights(hc, labels)
    df <- data.frame(
      type           = tp,
      label          = hc$labels,
      cluster        = factor(ifelse(is.na(labels), NA, paste0("C", labels)),
                              levels = paste0("C", sort(unique(stats::na.omit(labels))))),
      dendro_order   = match(seq_along(hc$labels), hc$order),
      cluster_height = heights$cluster_height[labels],
      cut_height     = heights$cut_height[labels],
      stringsAsFactors = FALSE
    )
    df <- df[order(df$dendro_order), ]

    list(clusters = df, k = nlevels(df$cluster), selection = cut$selection, setting = cut$setting,
         plot = .plot_dendro_clusters(hc, labels, heights, title = paste0(
           "Dendrogram clusters (", tp, ") - ", method, ": ", nlevels(df$cluster), " clusters")))
  })
  names(results) <- types

  clusters <- do.call(rbind, lapply(results, `[[`, "clusters"))
  rownames(clusters) <- NULL

  # Add the group id and metadata column for rows
  key <- unique(data.frame(label = as.character(dot_df$cluster_var), id = as.character(dot_df$id),
                           var = as.character(dot_df$var), stringsAsFactors = FALSE))
  clusters$id  <- key$id[match(ifelse(clusters$type == "rows", clusters$label, NA), key$label)]
  clusters$var <- key$var[match(ifelse(clusters$type == "rows", clusters$label, NA), key$label)]
  clusters <- clusters[, c("type", "label", "id", "var", "cluster", "dendro_order", "cluster_height", "cut_height")]

  out <- list(
    clusters   = clusters,
    k          = vapply(results, `[[`, integer(1), "k"),
    selection  = lapply(results, `[[`, "selection"),
    plot       = lapply(results, `[[`, "plot"),
    hclust     = trees$hclust[types],
    parameters = list(method = method, what = what, minClusterSize = minClusterSize,
                      assign_all = assign_all, setting = lapply(results, `[[`, "setting"))
  )

  # --- Write the row clusters back to the cells -------------------------------------
  if (!is.null(seurat_obj)) {
    out$seurat_obj <- .write_dendro_clusters(seurat_obj, clusters, meta_col, new_col)
  }

  # --- Save ---------------------------------------------------------------------------
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(clusters, file.path(output_dir, paste0("DendroClusters_", what, "_", method, "_", ts, ".csv")),
                     row.names = FALSE)
    for (tp in types) {
      n_leaves <- length(trees$hclust[[tp]]$labels)
      ggsave(file.path(output_dir, paste0("DendroClusters_", tp, "_", method, "_", ts, ".jpg")),
             out$plot[[tp]], width = min(40, max(8, 0.25 * n_leaves)), height = 7, dpi = 300, limitsize = FALSE)
    }
    message("Saved dendrogram clusters to: ", output_dir)
  }

  out
}

# -------------------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------------------

#' Matrix, trees and distances behind plot_dot_dendro_multi()
#' @keywords internal
#' @noRd
.dot_dendro_trees <- function(dot_df) {
  needed <- c("cluster_var", "features.plot", "avg.exp.scaled", "id", "var")
  missing <- setdiff(needed, colnames(dot_df))
  if (!is.data.frame(dot_df) || length(missing) > 0) {
    stop("`dot_df` must be the data frame returned by plot_dot_dendro_multi() ",
         "(missing columns: ", paste(missing, collapse = ", "), ").", call. = FALSE)
  }

  mat <- attr(dot_df, "dendro_matrix")
  row_hc <- attr(dot_df, "row_hclust")
  col_hc <- attr(dot_df, "col_hclust")

  if (is.null(mat) || is.null(row_hc) || is.null(col_hc)) {
    # Rebuild exactly as plot_dot_dendro_multi() does
    long <- data.frame(cluster_var = as.character(dot_df$cluster_var),
                       feature = as.character(dot_df$features.plot),
                       value = dot_df$avg.exp.scaled, stringsAsFactors = FALSE)
    wide <- tidyr::pivot_wider(long, names_from = "feature", values_from = "value")
    mat <- as.matrix(tibble::column_to_rownames(as.data.frame(wide), "cluster_var"))
    if (anyNA(mat)) {
      stop("The scaled expression matrix contains missing values; the dendrograms cannot be rebuilt.", call. = FALSE)
    }
    row_hc <- stats::hclust(stats::dist(mat))
    col_hc <- stats::hclust(stats::dist(t(mat)))
  }
  mat <- as.matrix(mat)

  list(
    matrix = mat,
    hclust = list(rows = row_hc, genes = col_hc),
    dist   = list(rows = stats::dist(mat), genes = stats::dist(t(mat)))
  )
}

#' Average silhouette width; unassigned leaves (NA) count as 0
#' @keywords internal
#' @noRd
.avg_silhouette <- function(labels, dm) {
  ok <- !is.na(labels)
  if (length(unique(labels[ok])) < 2) return(NA_real_)
  lab <- labels[ok]
  d <- dm[ok, ok, drop = FALSE]
  idx <- seq_along(lab)
  s <- vapply(idx, function(i) {
    own <- lab == lab[i]
    if (sum(own) == 1) return(0)
    a <- mean(d[i, own & idx != i])
    b <- min(vapply(setdiff(unique(lab), lab[i]), function(cl) mean(d[i, lab == cl]), numeric(1)))
    if (max(a, b) == 0) 0 else (b - a) / max(a, b)
  }, numeric(1))
  sum(s) / length(labels)
}

#' Choose the candidate with the best average silhouette (ties: first candidate)
#' @keywords internal
#' @noRd
.pick_best <- function(selection, fallback) {
  if (all(is.na(selection$avg_silhouette))) return(fallback)
  which.max(selection$avg_silhouette)
}

#' Hybrid dynamic tree cut (dynamicTreeCut)
#' @keywords internal
#' @noRd
.cut_dynamic <- function(hc, dm, deepSplit, minClusterSize) {
  run <- function(ds) {
    lab <- tryCatch(
      dynamicTreeCut::cutreeDynamic(dendro = hc, distM = dm, method = "hybrid", deepSplit = ds,
                                    minClusterSize = minClusterSize, pamStage = FALSE, verbose = 0),
      error = function(e) rep(0, length(hc$labels))
    )
    lab <- as.integer(lab)
    lab[lab == 0] <- NA
    lab
  }

  if (!identical(deepSplit, "auto")) {
    return(list(labels = run(deepSplit), selection = NULL, setting = c(deepSplit = deepSplit)))
  }

  candidates <- 0:4
  runs <- lapply(candidates, run)
  selection <- data.frame(
    deepSplit      = candidates,
    n_clusters     = vapply(runs, function(l) length(unique(stats::na.omit(l))), integer(1)),
    n_unassigned   = vapply(runs, function(l) sum(is.na(l)), integer(1)),
    avg_silhouette = vapply(runs, .avg_silhouette, numeric(1), dm = dm)
  )
  best <- .pick_best(selection, fallback = 3)                     # deepSplit = 2 (package default)
  selection$selected <- seq_len(nrow(selection)) == best
  list(labels = runs[[best]], selection = selection, setting = c(deepSplit = candidates[best]))
}

#' Inconsistency coefficient of every merge (depth 2: the merge and the merges below it)
#' @keywords internal
#' @noRd
.inconsistency <- function(hc) {
  m <- hc$merge
  ht <- hc$height
  vapply(seq_len(nrow(m)), function(i) {
    below <- m[i, m[i, ] > 0]
    hs <- c(ht[i], ht[below])
    if (length(hs) < 2 || stats::sd(hs) == 0) 0 else (ht[i] - mean(hs)) / stats::sd(hs)
  }, numeric(1))
}

#' Members (leaf indices) and parent of every merge node
#' @keywords internal
#' @noRd
.tree_nodes <- function(hc) {
  m <- hc$merge
  n <- nrow(m)
  members <- vector("list", n)
  parent <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    kids <- lapply(m[i, ], function(x) if (x < 0) -x else members[[x]])
    members[[i]] <- unlist(kids)
    for (x in m[i, ]) if (x > 0) parent[x] <- i
  }
  list(members = members, parent = parent)
}

#' Cut where a subtree contains no inconsistent merge
#' @keywords internal
#' @noRd
.cut_inconsistency <- function(hc, dm, cutoff, minClusterSize) {
  coef <- .inconsistency(hc)
  m <- hc$merge
  n_leaves <- length(hc$labels)

  # largest coefficient anywhere in each subtree
  sub_max <- numeric(nrow(m))
  for (i in seq_len(nrow(m))) {
    kids <- m[i, m[i, ] > 0]
    sub_max[i] <- max(c(coef[i], sub_max[kids]))
  }

  run <- function(cut_value) {
    labels <- integer(n_leaves)
    next_id <- 0L
    assign_node <- function(node) {
      if (node < 0) {                                  # a leaf on its own
        next_id <<- next_id + 1L
        labels[-node] <<- next_id
      } else if (sub_max[node] <= cut_value) {         # consistent subtree -> one cluster
        next_id <<- next_id + 1L
        labels[.leaves_of(m, node)] <<- next_id
      } else {
        assign_node(m[node, 1])
        assign_node(m[node, 2])
      }
    }
    assign_node(nrow(m))
    sizes <- table(labels)
    labels[labels %in% as.integer(names(sizes)[sizes < minClusterSize])] <- NA
    labels
  }

  if (!identical(cutoff, "auto")) {
    return(list(labels = run(cutoff), selection = NULL, setting = c(cutoff = cutoff)))
  }

  candidates <- sort(unique(round(c(coef, max(coef) + 1), 10)))
  runs <- lapply(candidates, run)
  selection <- data.frame(
    cutoff         = candidates,
    n_clusters     = vapply(runs, function(l) length(unique(stats::na.omit(l))), integer(1)),
    n_unassigned   = vapply(runs, function(l) sum(is.na(l)), integer(1)),
    avg_silhouette = vapply(runs, .avg_silhouette, numeric(1), dm = dm)
  )
  best <- .pick_best(selection, fallback = length(candidates))
  selection$selected <- seq_len(nrow(selection)) == best
  list(labels = runs[[best]], selection = selection, setting = c(cutoff = candidates[best]))
}

#' Leaf indices below a merge node
#' @keywords internal
#' @noRd
.leaves_of <- function(m, node) {
  if (node < 0) return(-node)
  c(.leaves_of(m, m[node, 1]), .leaves_of(m, m[node, 2]))
}

#' Give unassigned leaves to the cluster with the smallest average distance
#' @keywords internal
#' @noRd
.assign_leftovers <- function(labels, dm) {
  todo <- which(is.na(labels))
  clusters <- unique(stats::na.omit(labels))
  if (length(todo) == 0 || length(clusters) == 0) return(labels)
  for (i in todo) {
    avg <- vapply(clusters, function(cl) mean(dm[i, which(labels == cl)]), numeric(1))
    labels[i] <- clusters[which.min(avg)]
  }
  labels
}

#' Renumber clusters 1..k from left to right along the dendrogram
#' @keywords internal
#' @noRd
.relabel_in_dendro_order <- function(labels, hc) {
  ordered <- labels[hc$order]
  first_seen <- unique(stats::na.omit(ordered))
  out <- match(labels, first_seen)
  out[is.na(labels)] <- NA
  as.integer(out)
}

#' Height of the merge joining each cluster, and the height at which it was cut
#' @keywords internal
#' @noRd
.cluster_heights <- function(hc, labels) {
  nodes <- .tree_nodes(hc)
  ids <- sort(unique(stats::na.omit(labels)))
  cluster_height <- cut_height <- stats::setNames(numeric(length(ids)), ids)

  for (cl in ids) {
    leaves <- which(labels == cl)
    # smallest merge containing every leaf of the cluster
    node <- which(vapply(nodes$members, function(mm) all(leaves %in% mm), logical(1)))[1]
    if (length(leaves) == 1) {
      top <- 0
      parent_node <- which(vapply(nodes$members, function(mm) leaves %in% mm, logical(1)))[1]
    } else {
      top <- hc$height[node]
      parent_node <- nodes$parent[node]
    }
    parent_height <- if (is.na(parent_node)) top else hc$height[parent_node]
    cluster_height[as.character(cl)] <- top
    cut_height[as.character(cl)] <- (top + parent_height) / 2
  }
  list(cluster_height = unname(cluster_height), cut_height = unname(cut_height))
}

#' Dendrogram with branches coloured by cluster and a mark at each cut height
#' @keywords internal
#' @noRd
.plot_dendro_clusters <- function(hc, labels, heights, title) {
  m <- hc$merge
  n_merge <- nrow(m)
  nodes <- .tree_nodes(hc)
  leaf_x <- match(seq_along(hc$labels), hc$order)

  node_x <- numeric(n_merge)
  coord <- function(x) if (x < 0) c(leaf_x[-x], 0) else c(node_x[x], hc$height[x])
  group_of <- function(leaves) {
    cl <- unique(labels[leaves])
    if (length(cl) == 1 && !is.na(cl)) paste0("C", cl) else if (all(is.na(cl))) "Unassigned" else "Between clusters"
  }

  segs <- list()
  for (i in seq_len(n_merge)) {
    a <- coord(m[i, 1]); b <- coord(m[i, 2])
    node_x[i] <- (a[1] + b[1]) / 2
    y <- hc$height[i]
    kids_group <- lapply(m[i, ], function(x) group_of(if (x < 0) -x else nodes$members[[x]]))
    segs[[length(segs) + 1]] <- data.frame(x = a[1], xend = a[1], y = a[2], yend = y, group = kids_group[[1]])
    segs[[length(segs) + 1]] <- data.frame(x = b[1], xend = b[1], y = b[2], yend = y, group = kids_group[[2]])
    segs[[length(segs) + 1]] <- data.frame(x = a[1], xend = b[1], y = y, yend = y,
                                           group = group_of(nodes$members[[i]]))
  }
  segs <- do.call(rbind, segs)

  ids <- sort(unique(stats::na.omit(labels)))
  cluster_levels <- paste0("C", ids)
  colours <- c(stats::setNames(make_cluster_palette(max(1, length(ids)), "kelly")[seq_along(ids)], cluster_levels),
               "Unassigned" = "grey70", "Between clusters" = "grey20")
  segs$group <- factor(segs$group, levels = names(colours))

  leaves <- data.frame(x = leaf_x, label = hc$labels,
                       group = factor(ifelse(is.na(labels), "Unassigned", paste0("C", labels)), levels = names(colours)))

  cuts <- do.call(rbind, lapply(ids, function(cl) {
    xs <- leaf_x[labels == cl & !is.na(labels)]
    data.frame(x = min(xs) - 0.45, xend = max(xs) + 0.45, y = heights$cut_height[cl], group = paste0("C", cl))
  }))

  ymax <- max(hc$height)
  p <- ggplot() +
    geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend, colour = group), linewidth = 0.6) +
    geom_text(data = leaves, aes(x = x, y = -0.02 * ymax, label = label, colour = group),
              angle = 90, hjust = 1, size = 3, show.legend = FALSE) +
    scale_colour_manual(values = colours, drop = TRUE, name = "Cluster") +
    guides(colour = guide_legend(override.aes = list(linewidth = 2))) +
    # leave room below the tree for the (rotated) leaf labels
    scale_y_continuous(expand = expansion(mult = c(max(0.3, 0.03 * max(nchar(hc$labels))), 0.05)),
                       breaks = function(lim) { b <- pretty(c(0, lim[2])); b[b >= 0] }) +
    labs(title = title, x = NULL, y = "Height") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())

  if (!is.null(cuts) && nrow(cuts) > 0) {
    p <- p + geom_segment(data = cuts, aes(x = x, xend = xend, y = y, yend = y),
                          linetype = "dashed", colour = "red3", linewidth = 0.5)
  }
  p
}

#' Map row clusters of one metadata column back to the cells
#' @keywords internal
#' @noRd
.write_dendro_clusters <- function(seurat_obj, clusters, meta_col, new_col) {
  rows <- clusters[clusters$type == "rows", ]
  row_levels <- levels(droplevels(rows$cluster))
  if (!meta_col %in% rows$var) {
    stop("`meta_col` must be one of the metadata columns used in plot_dot_dendro_multi(): ",
         paste(unique(rows$var), collapse = ", "), ".", call. = FALSE)
  }
  if (!meta_col %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column '", meta_col, "' is not in `seurat_obj`.", call. = FALSE)
  }
  rows <- rows[rows$var == meta_col, ]
  lookup <- stats::setNames(as.character(rows$cluster), rows$id)
  cell_clusters <- unname(lookup[as.character(seurat_obj@meta.data[[meta_col]])])
  seurat_obj[[new_col]] <- factor(cell_clusters, levels = row_levels)
  seurat_obj
}
