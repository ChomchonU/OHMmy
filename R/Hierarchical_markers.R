# Markers along a resolution tree, and a dichotomous key per split.
#
# find_hierarchical_markers() links clusters across increasing resolutions into a
# tree (each cluster -> the lower-resolution cluster holding most of its cells),
# finds markers between sibling clusters at every split, ranks them by AUC and
# builds a short decision-tree key (rpart) per split.

#' Find Markers Along a Resolution Tree and Build a Dichotomous Key
#'
#' Tracks clusters from low to high resolution and finds, at every point where a
#' cluster splits, the genes that distinguish its sub-clusters. Low resolutions give
#' broad markers of major populations; higher resolutions give finer markers of
#' sub-populations. For every split, a short decision tree (\pkg{rpart}) picks the
#' genes and thresholds that sort cells into the sub-clusters: a dichotomous key.
#'
#' @details
#' \strong{1. Tree.} Resolutions are ordered from low to high. Each cluster is linked
#' to the cluster of the previous resolution that holds most of its cells (its
#' parent); that share is its \code{purity}. Clusterings at different resolutions are
#' not always nested, so a low purity means the cluster takes cells from several
#' parents. The root holds all cells, so its children are the clusters of the lowest
#' resolution. A cluster with a single child (unchanged at the next resolution) is
#' not a split.
#'
#' \strong{2. Markers per split.} At each split, every child is compared with its
#' siblings only, using the parent's cells (\code{Seurat::FindMarkers()}, Wilcoxon,
#' genes higher in the child). Cells of the parent that moved to a child of another
#' parent are left out. Children with fewer than \code{min_cells} cells are skipped.
#' Genes are kept when \code{p_val_adj < padj_max}, \code{avg_log2FC >= min_log2fc},
#' \code{pct.1 >= min_pct} and \code{pct.1 - pct.2 >= min_pct_diff}, then ranked by
#' AUC (the probability that a random child cell expresses the gene more highly than
#' a random sibling cell; 0.5 = no separation, 1 = perfect), with the detection
#' difference as tie-break. Clusters were defined on the same data, so p-values are
#' optimistic: use them as a filter, not as evidence on their own.
#'
#' \strong{3. Patterns.} Each marker row gets \code{also_marks_parent} (the gene also
#' marks the parent against the parent's siblings), \code{also_marks_children} (it
#' also separates the child's own sub-clusters; \code{NA} for a child that does not
#' split further) and a \code{pattern}:
#' \itemize{
#'   \item \code{"unstable"}: the child (or a cluster it passes through unchanged) has a
#'     purity below \code{min_purity}, so the split is not clean across resolutions;
#'   \item \code{"specific"}: the child does not split further and the gene also
#'     separates it from all other cells (\code{auc_vs_all >= specific_auc});
#'   \item \code{"lineage"}: the child splits further, but the gene does not separate
#'     its sub-clusters (shared by all of them);
#'   \item \code{"split"}: any other marker of this split.
#' }
#'
#' \strong{4. Key.} For each split, \code{rpart} is given the top
#' \code{max_candidate_genes} ranked markers of every child, grown to
#' \code{max_depth} questions and pruned by cross-validation (1-SE rule: the smallest
#' tree whose error is within one standard error of the best). The data therefore
#' decide how many genes each split needs. For every question, the chosen gene
#' (\code{role = "key"}) and up to \code{n_backup} surrogate genes that sort the cells
#' almost the same way (\code{role = "backup"}) are reported with all ranking values.
#'
#' Expression is read from the \code{data} layer of \code{assay}; split Seurat v5
#' layers are joined on a local copy.
#'
#' @param seurat_obj A \code{Seurat} object with normalised data and cluster columns
#'   at several resolutions.
#' @param res_cols Character vector of metadata columns holding the clusters, e.g.
#'   \code{c("RNA_snn_res.0.2", "RNA_snn_res.0.5", "RNA_snn_res.1")}. Ordered by the
#'   resolution in their name when it can be read, otherwise used as given. Give this
#'   or \code{pipeline}.
#' @param pipeline Character. Instead of \code{res_cols}: use every
#'   \code{<graph>_res.<resolution>} column of this method, named as in
#'   \code{\link{analyze_and_save_clusters}()} (e.g. \code{"log.fastmnn"} or
#'   \code{"RNA_snn"}).
#' @param assay Character. Assay used for expression. Default is \code{"RNA"}.
#' @param padj_max,min_log2fc,min_pct,min_pct_diff Numeric marker filters: adjusted
#'   p-value below \code{padj_max} (default 0.05), \code{avg_log2FC} at least
#'   \code{min_log2fc} (default 0.5), detected in at least \code{min_pct} of the child's
#'   cells (default 0.5), and detection difference \code{pct.1 - pct.2} at least
#'   \code{min_pct_diff} (default 0.3).
#' @param min_purity Numeric. Children below this purity are flagged and their markers
#'   labelled \code{"unstable"} (they are still tested). Default is 0.8.
#' @param min_cells Integer. Children with fewer cells in a split are skipped. Default 20.
#' @param specific_auc Numeric. AUC against all other cells needed for
#'   \code{"specific"}. Default is 0.8.
#' @param max_candidate_genes Integer. Top-ranked markers per child offered to
#'   \code{rpart}. Default is 100.
#' @param max_depth Integer or \code{NULL}. Maximum number of questions on any path of a
#'   split's key before pruning. \code{NULL} (default) uses the number of children of
#'   the split (at most 30).
#' @param n_backup Integer. Backup (surrogate) genes reported per question. Default 10.
#' @param cv_folds Integer. Cross-validation folds for pruning. Default is 10.
#' @param seed Integer. Seed for the cross-validation folds (the global random state
#'   is restored afterwards). Default is 42.
#' @param sample_name Character. Used in plot titles and file names. Default "Sample".
#' @param output_dir Character or \code{NULL}. If given, the tables (CSV) and plots
#'   (JPEG) are saved there. Default is \code{NULL} (nothing saved).
#' @param verbose Logical. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{tree}: one row per cluster (node): \code{node}, \code{layer},
#'     \code{resolution}, \code{res_col}, \code{cluster}, \code{parent},
#'     \code{n_cells}, \code{purity}, \code{n_children}, \code{low_purity}.
#'   \item \code{splits}: one row per split: \code{split}, \code{layer},
#'     \code{children}, \code{skipped_children} (too few cells), \code{n_cells},
#'     \code{n_markers}, \code{n_key_genes}, \code{cv_accuracy},
#'     \code{train_accuracy}, \code{status}.
#'   \item \code{markers}: every gene passing the filters, per split and child:
#'     \code{split}, \code{child}, \code{gene}, \code{rank}, \code{auc},
#'     \code{auc_vs_all}, \code{p_val}, \code{p_val_adj}, \code{avg_log2FC},
#'     \code{pct.1}, \code{pct.2}, \code{pct_diff}, \code{mean_child},
#'     \code{mean_siblings}, \code{also_marks_parent}, \code{also_marks_children},
#'     \code{pattern}.
#'   \item \code{key}: the key and backup genes of every split: \code{split},
#'     \code{rpart_node}, \code{role}, \code{gene}, \code{threshold},
#'     \code{rule} (e.g. \code{"GZMK >= 1.2 -> res0.5_3"}), \code{high_side} (the child
#'     predicted for cells at or above the threshold), \code{agreement} (backup genes:
#'     rpart's adjusted agreement with the key gene), \code{importance}, then
#'     \code{marker_of} and the ranking columns of \code{markers} for that gene.
#'   \item \code{confusion}: per split, a table of true vs key-predicted child.
#'   \item \code{plots}: \code{tree} (the resolution tree with key genes at each split)
#'     and \code{dotplot} (key genes across the clusters of the highest resolution).
#'   \item \code{res_cols} and \code{parameters}.
#' }
#'
#' @seealso \code{\link{analyze_and_save_clusters}()} to choose resolutions,
#'   \code{\link{plot_dot_dendro_multi}()} to plot the key genes.
#' @export
#'
#' @examples
#' \dontrun{
#' res <- find_hierarchical_markers(
#'   seurat_obj = my_seurat,
#'   res_cols   = c("RNA_snn_res.0.1", "RNA_snn_res.0.3", "RNA_snn_res.0.8"),
#'   output_dir = "Hierarchical_markers"
#' )
#' res$key                                  # key + backup genes with their ranking values
#' subset(res$markers, pattern == "lineage")
#' res$plots$tree
#'
#' # All resolutions of one integration method
#' res <- find_hierarchical_markers(my_seurat, pipeline = "log.fastmnn")
#' }
find_hierarchical_markers <- function(seurat_obj,
                                      res_cols = NULL,
                                      pipeline = NULL,
                                      assay = "RNA",
                                      padj_max = 0.05,
                                      min_log2fc = 0.5,
                                      min_pct = 0.5,
                                      min_pct_diff = 0.3,
                                      min_purity = 0.8,
                                      min_cells = 20,
                                      specific_auc = 0.8,
                                      max_candidate_genes = 100,
                                      max_depth = NULL,
                                      n_backup = 10,
                                      cv_folds = 10,
                                      seed = 42,
                                      sample_name = "Sample",
                                      output_dir = NULL,
                                      verbose = TRUE) {
  if (!inherits(seurat_obj, "Seurat")) stop("`seurat_obj` must be a Seurat object.", call. = FALSE)
  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop("Assay '", assay, "' not found in `seurat_obj`.", call. = FALSE)
  }
  res_cols <- .resolve_res_cols(seurat_obj@meta.data, res_cols, pipeline)
  say <- function(...) if (verbose) message(...)

  # Cells without a cluster at any resolution cannot be placed in the tree
  md <- seurat_obj@meta.data[, res_cols, drop = FALSE]
  complete <- stats::complete.cases(md)
  if (!all(complete)) {
    warning(sum(!complete), " cells with a missing cluster label are left out.", call. = FALSE)
    seurat_obj <- seurat_obj[, complete]
    md <- md[complete, , drop = FALSE]
  }
  seurat_obj <- .join_data_layers(seurat_obj, assay)
  expr <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = "data")

  # 1. Tree across resolutions
  built <- .build_resolution_tree(md, res_cols, min_purity)
  tree <- built$tree
  cell_nodes <- built$cell_nodes
  split_nodes <- tree$node[tree$n_children >= 2]
  say("Resolution tree: ", nrow(tree) - 1, " clusters over ", length(res_cols), " resolutions, ",
      length(split_nodes), " splits.")

  # 2. Markers per split (child vs siblings within the parent)
  filters <- list(padj_max = padj_max, min_log2fc = min_log2fc, min_pct = min_pct,
                  min_pct_diff = min_pct_diff)
  per_split <- lapply(split_nodes, function(split_node) {
    say("  Split ", split_node, " ...")
    .markers_for_split(seurat_obj, expr, tree, cell_nodes, split_node, assay, filters, min_cells)
  })
  names(per_split) <- split_nodes
  markers <- do.call(rbind, lapply(per_split, `[[`, "markers"))
  if (is.null(markers)) markers <- .empty_marker_table()
  rownames(markers) <- NULL

  # 3. Patterns along the tree
  markers <- .add_auc_vs_all(markers, per_split, expr, tree)
  markers <- .label_marker_patterns(markers, tree, min_purity, specific_auc)

  # 4. Dichotomous key per split
  keys <- lapply(split_nodes, function(split_node) {
    .key_for_split(per_split[[split_node]], markers, expr, split_node, max_candidate_genes,
                   max_depth, n_backup, cv_folds, seed)
  })
  names(keys) <- split_nodes
  key <- do.call(rbind, lapply(keys, `[[`, "key"))
  if (is.null(key)) key <- .empty_key_table(markers)
  rownames(key) <- NULL

  splits <- .summarise_splits(split_nodes, tree, per_split, keys, markers)

  # 5. Plots
  plots <- list(
    tree    = .plot_resolution_tree(tree, key, sample_name),
    dotplot = .plot_key_dotplot(seurat_obj, key, res_cols[length(res_cols)], assay, sample_name)
  )

  out <- list(
    tree = tree, splits = splits, markers = markers, key = key,
    confusion = lapply(keys, `[[`, "confusion"), plots = plots, res_cols = res_cols,
    parameters = list(assay = assay, padj_max = padj_max, min_log2fc = min_log2fc, min_pct = min_pct,
                      min_pct_diff = min_pct_diff, min_purity = min_purity, min_cells = min_cells,
                      specific_auc = specific_auc, max_candidate_genes = max_candidate_genes,
                      max_depth = max_depth, n_backup = n_backup, cv_folds = cv_folds, seed = seed)
  )

  if (!is.null(output_dir)) .save_hierarchical_markers(out, output_dir, sample_name, say)
  out
}

# -------------------------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------------------------

#' Resolution columns from `res_cols` or `pipeline`, ordered from low to high
#' @keywords internal
#' @noRd
.resolve_res_cols <- function(md, res_cols, pipeline) {
  if (is.null(res_cols) == is.null(pipeline)) {
    stop("Give exactly one of `res_cols` or `pipeline`.", call. = FALSE)
  }
  if (!is.null(pipeline)) {
    all_res <- grep("_res\\.[0-9.]+$", colnames(md), value = TRUE)
    res_cols <- all_res[.pipeline_from_res_col(all_res) == pipeline]
    if (length(res_cols) == 0) {
      stop("No resolution columns found for pipeline '", pipeline, "'. Available: ",
           paste(unique(.pipeline_from_res_col(all_res)), collapse = ", "), call. = FALSE)
    }
  }
  missing <- setdiff(res_cols, colnames(md))
  if (length(missing) > 0) {
    stop("Column(s) not found in the metadata: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  resolutions <- .resolution_of(res_cols)
  if (!anyNA(resolutions)) res_cols <- res_cols[order(resolutions)]
  .drop_repeated_partitions(md, res_cols)
}

#' Resolution value in a `<graph>_res.<resolution>` column name (NA if absent)
#' @keywords internal
#' @noRd
.resolution_of <- function(res_cols) {
  value <- rep(NA_real_, length(res_cols))
  named <- grepl("_res\\.[0-9.]+$", res_cols)
  value[named] <- suppressWarnings(as.numeric(sub(".*_res\\.", "", res_cols[named])))
  value
}

#' Drop a resolution whose clusters are identical to the previous one (no new split)
#' @keywords internal
#' @noRd
.drop_repeated_partitions <- function(md, res_cols) {
  keep <- rep(TRUE, length(res_cols))
  for (i in seq_along(res_cols)[-1]) {
    previous <- as.character(md[[res_cols[max(which(keep[seq_len(i - 1)]))]]])
    current <- as.character(md[[res_cols[i]]])
    same <- length(unique(paste(previous, current))) == length(unique(previous)) &&
      length(unique(previous)) == length(unique(current))
    if (same) keep[i] <- FALSE
  }
  if (any(!keep)) {
    message("Same clusters as the previous resolution, skipped: ", paste(res_cols[!keep], collapse = ", "))
  }
  res_cols[keep]
}

#' Join split Seurat v5 data layers on this local copy
#' @keywords internal
#' @noRd
.join_data_layers <- function(seurat_obj, assay) {
  if (inherits(seurat_obj[[assay]], "Assay5") &&
      length(SeuratObject::Layers(seurat_obj[[assay]], search = "data")) > 1) {
    seurat_obj[[assay]] <- SeuratObject::JoinLayers(seurat_obj[[assay]])
  }
  seurat_obj
}

#' Link clusters across resolutions: parent = the previous-resolution cluster holding
#' most of the cluster's cells
#'
#' Returns the node table and a cells x (1 + resolutions) matrix of node ids
#' (column 1 is the root).
#' @keywords internal
#' @noRd
.build_resolution_tree <- function(md, res_cols, min_purity) {
  resolutions <- .resolution_of(res_cols)
  layer_prefix <- ifelse(is.na(resolutions), paste0("L", seq_along(res_cols)), paste0("res", resolutions))

  cell_nodes <- matrix("root", nrow = nrow(md), ncol = length(res_cols) + 1,
                       dimnames = list(rownames(md), c("root", res_cols)))
  rows <- list(data.frame(node = "root", layer = 0L, resolution = NA_real_, res_col = NA_character_,
                          cluster = NA_character_, parent = NA_character_, n_cells = nrow(md),
                          purity = 1, stringsAsFactors = FALSE))

  for (layer in seq_along(res_cols)) {
    clusters <- as.character(md[[res_cols[layer]]])
    nodes <- paste0(layer_prefix[layer], "_", clusters)
    cell_nodes[, layer + 1] <- nodes

    overlap <- table(nodes, cell_nodes[, layer])        # rows: clusters, columns: candidate parents
    node_ids <- rownames(overlap)
    rows[[layer + 1]] <- data.frame(
      node       = node_ids,
      layer      = layer,
      resolution = resolutions[layer],
      res_col    = res_cols[layer],
      cluster    = sub(paste0("^", layer_prefix[layer], "_"), "", node_ids),
      parent     = colnames(overlap)[apply(overlap, 1, which.max)],
      n_cells    = as.integer(rowSums(overlap)),
      purity     = apply(overlap, 1, max) / rowSums(overlap),
      stringsAsFactors = FALSE
    )
    # Numeric cluster names in numeric order (2 before 10)
    layer_rows <- rows[[layer + 1]]
    rows[[layer + 1]] <- layer_rows[order(suppressWarnings(as.numeric(layer_rows$cluster)), layer_rows$cluster), ]
  }

  tree <- do.call(rbind, rows)
  rownames(tree) <- NULL
  tree$n_children <- vapply(tree$node, function(n) sum(tree$parent == n, na.rm = TRUE), integer(1),
                            USE.NAMES = FALSE)
  tree$low_purity <- tree$purity < min_purity
  list(tree = tree, cell_nodes = cell_nodes)
}

#' The node at which each node's unchanged chain starts (a child of a split, or root)
#'
#' A cluster that is its parent's only child is the same population seen at a higher
#' resolution, so it belongs to its parent's chain.
#' @keywords internal
#' @noRd
.chain_tops <- function(tree) {
  top <- stats::setNames(tree$node, tree$node)
  n_children <- stats::setNames(tree$n_children, tree$node)
  for (i in order(tree$layer)) {
    parent <- tree$parent[i]
    if (!is.na(parent) && n_children[[parent]] == 1) top[[tree$node[i]]] <- top[[parent]]
  }
  top
}

#' Children of one split: their cells, sizes and the child-vs-sibling markers
#' @keywords internal
#' @noRd
.markers_for_split <- function(seurat_obj, expr, tree, cell_nodes, split_node, assay, filters, min_cells) {
  split_layer <- tree$layer[tree$node == split_node]
  children <- tree$node[!is.na(tree$parent) & tree$parent == split_node]

  # Parent's cells that stayed in one of its own children
  in_split <- cell_nodes[, split_layer + 1] == split_node & cell_nodes[, split_layer + 2] %in% children
  child_of_cell <- cell_nodes[in_split, split_layer + 2]
  sizes <- table(factor(child_of_cell, levels = children))
  tested <- names(sizes)[sizes >= min_cells]
  skipped <- setdiff(children, tested)

  result <- list(children = children, tested = tested, skipped = skipped,
                 cells = split(names(child_of_cell), child_of_cell)[tested], markers = NULL)
  if (length(tested) < 2) return(result)

  result$markers <- do.call(rbind, lapply(tested, function(child) {
    .child_vs_siblings(seurat_obj, expr, result$cells, child, split_node, assay, filters)
  }))
  result
}

#' Markers higher in one child than in its siblings, filtered and ranked by AUC
#' @keywords internal
#' @noRd
.child_vs_siblings <- function(seurat_obj, expr, cells_by_child, child, split_node, assay, filters) {
  child_cells <- cells_by_child[[child]]
  sibling_cells <- unlist(cells_by_child[names(cells_by_child) != child], use.names = FALSE)

  de <- tryCatch(
    suppressWarnings(suppressMessages(Seurat::FindMarkers(
      seurat_obj, ident.1 = child_cells, ident.2 = sibling_cells, assay = assay,
      test.use = "wilcox", only.pos = TRUE, logfc.threshold = filters$min_log2fc,
      min.pct = filters$min_pct, min.diff.pct = filters$min_pct_diff, verbose = FALSE
    ))),
    error = function(e) NULL
  )
  if (is.null(de) || nrow(de) == 0) return(NULL)

  de$gene <- rownames(de)
  de$pct_diff <- de$pct.1 - de$pct.2
  de <- de[de$p_val_adj < filters$padj_max & de$avg_log2FC >= filters$min_log2fc &
             de$pct.1 >= filters$min_pct & de$pct_diff >= filters$min_pct_diff, , drop = FALSE]
  if (nrow(de) == 0) return(NULL)

  values <- as.matrix(expr[de$gene, c(child_cells, sibling_cells), drop = FALSE])
  in_child <- seq_along(child_cells)
  de$auc <- apply(values, 1, function(x) .auc(x[in_child], x[-in_child]))
  de$mean_child <- rowMeans(values[, in_child, drop = FALSE])
  de$mean_siblings <- rowMeans(values[, -in_child, drop = FALSE])

  de <- de[order(-de$auc, -de$pct_diff), , drop = FALSE]
  data.frame(split = split_node, child = child, gene = de$gene, rank = seq_len(nrow(de)),
             auc = de$auc, auc_vs_all = NA_real_, p_val = de$p_val, p_val_adj = de$p_val_adj,
             avg_log2FC = de$avg_log2FC, pct.1 = de$pct.1, pct.2 = de$pct.2, pct_diff = de$pct_diff,
             mean_child = de$mean_child, mean_siblings = de$mean_siblings,
             stringsAsFactors = FALSE)
}

#' AUC: probability that a random cell of group 1 has a higher value than a random
#' cell of group 2 (ties count half); the Mann-Whitney U divided by n1 * n2
#' @keywords internal
#' @noRd
.auc <- function(x1, x2) {
  n1 <- length(x1)
  n2 <- length(x2)
  ranks <- rank(c(x1, x2))
  (sum(ranks[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}

#' Empty marker table with the final columns (no split had markers)
#' @keywords internal
#' @noRd
.empty_marker_table <- function() {
  data.frame(split = character(), child = character(), gene = character(), rank = integer(),
             auc = numeric(), auc_vs_all = numeric(), p_val = numeric(), p_val_adj = numeric(),
             avg_log2FC = numeric(), pct.1 = numeric(), pct.2 = numeric(), pct_diff = numeric(),
             mean_child = numeric(), mean_siblings = numeric(), stringsAsFactors = FALSE)
}

#' AUC of each marker of a child that does not split further, against all other cells
#' @keywords internal
#' @noRd
.add_auc_vs_all <- function(markers, per_split, expr, tree) {
  if (nrow(markers) == 0) return(markers)
  tops <- .chain_tops(tree)
  has_own_split <- unique(tops[tree$node[tree$n_children >= 2]])
  leaf_rows <- which(!markers$child %in% has_own_split)

  for (child in unique(markers$child[leaf_rows])) {
    rows <- leaf_rows[markers$child[leaf_rows] == child]
    child_cells <- per_split[[markers$split[rows[1]]]]$cells[[child]]
    in_child <- colnames(expr) %in% child_cells
    values <- as.matrix(expr[markers$gene[rows], , drop = FALSE])
    markers$auc_vs_all[rows] <- apply(values, 1, function(x) .auc(x[in_child], x[!in_child]))
  }
  markers
}

#' Trajectory columns and pattern label of every marker row
#' @keywords internal
#' @noRd
.label_marker_patterns <- function(markers, tree, min_purity, specific_auc) {
  if (nrow(markers) == 0) {
    markers$also_marks_parent <- logical()
    markers$also_marks_children <- logical()
    markers$pattern <- character()
    return(markers)
  }
  tops <- .chain_tops(tree)
  split_nodes <- tree$node[tree$n_children >= 2]
  own_split <- stats::setNames(split_nodes, tops[split_nodes])     # chain top -> its split
  # Lowest purity along each chain (a plain named vector: a tapply() array would leak
  # its dim attribute through ifelse() into the pattern column)
  chain_purity <- vapply(split(tree$purity, tops[tree$node]), min, numeric(1))

  marker_key <- paste(markers$child, markers$gene)
  # The parent's chain starts at a child of the split above; root has no parent marker
  parent_top <- tops[markers$split]
  markers$also_marks_parent <- parent_top != "root" & paste(parent_top, markers$gene) %in% marker_key

  child_split <- unname(own_split[markers$child])
  split_genes <- paste(markers$split, markers$gene)
  markers$also_marks_children <- ifelse(is.na(child_split), NA,
                                        paste(child_split, markers$gene) %in% split_genes)

  unstable <- unname(chain_purity[markers$child] < min_purity)
  specific <- is.na(child_split) & !is.na(markers$auc_vs_all) & markers$auc_vs_all >= specific_auc
  lineage <- !is.na(child_split) & !markers$also_marks_children
  markers$pattern <- ifelse(unstable, "unstable",
                     ifelse(specific, "specific",
                     ifelse(lineage, "lineage", "split")))
  markers
}

#' rpart key for one split: grow, prune by the 1-SE rule, extract key + backup genes
#' @keywords internal
#' @noRd
.key_for_split <- function(split_info, markers, expr, split_node, max_candidate_genes, max_depth,
                           n_backup, cv_folds, seed) {
  empty <- list(key = NULL, confusion = NULL, cv_accuracy = NA_real_, train_accuracy = NA_real_)
  split_markers <- markers[markers$split == split_node, , drop = FALSE]
  if (length(split_info$tested) < 2 || nrow(split_markers) == 0) return(empty)

  candidates <- unique(split_markers$gene[split_markers$rank <= max_candidate_genes])
  cells <- unlist(split_info$cells, use.names = FALSE)
  child <- factor(rep(names(split_info$cells), lengths(split_info$cells)), levels = split_info$tested)

  # rpart needs syntactic column names (e.g. "HLA-DRA" is not); map them back afterwards
  values <- as.data.frame(t(as.matrix(expr[candidates, cells, drop = FALSE])))
  safe_names <- make.names(candidates, unique = TRUE)
  gene_of <- stats::setNames(candidates, safe_names)
  colnames(values) <- safe_names
  values$.child <- child

  depth <- min(30L, if (is.null(max_depth)) length(split_info$tested) else as.integer(max_depth))
  fit <- .with_seed(seed, rpart::rpart(
    .child ~ ., data = values, method = "class",
    control = rpart::rpart.control(maxdepth = depth, cp = 0, xval = cv_folds,
                                   maxcompete = 0, maxsurrogate = n_backup)
  ))
  pruned <- .prune_one_se(fit)

  predicted <- stats::predict(pruned, type = "class")
  confusion <- table(true = child, predicted = predicted)
  key <- .key_rows(pruned, gene_of, levels(child), split_node, split_markers)
  list(key = key, confusion = confusion, cv_accuracy = attr(pruned, "cv_accuracy"),
       train_accuracy = mean(as.character(predicted) == as.character(child)))
}

#' Run an expression with a fixed seed, then restore the caller's random state
#' @keywords internal
#' @noRd
.with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

#' Prune to the smallest tree whose cross-validated error is within 1 SE of the best
#'
#' Attaches the cross-validated accuracy of the chosen tree as attribute "cv_accuracy".
#' @keywords internal
#' @noRd
.prune_one_se <- function(fit) {
  cp_table <- fit$cptable
  root_error <- fit$frame$dev[1] / fit$frame$n[1]       # share of cells misclassified with no question
  if (!"xerror" %in% colnames(cp_table)) {
    attr(fit, "cv_accuracy") <- NA_real_
    return(fit)
  }
  best <- which.min(cp_table[, "xerror"])
  limit <- cp_table[best, "xerror"] + cp_table[best, "xstd"]
  chosen <- which(cp_table[, "xerror"] <= limit)[1]
  pruned <- rpart::prune(fit, cp = cp_table[chosen, "CP"])
  attr(pruned, "cv_accuracy") <- 1 - cp_table[chosen, "xerror"] * root_error
  pruned
}

#' Key (primary) and backup (surrogate) genes of a pruned rpart tree, with the
#' ranking values of each gene
#'
#' rpart stores, for every question (internal node, in frame order), one primary
#' split followed by its surrogates. A continuous split with ncat = -1 sends
#' values below the threshold left, ncat = +1 sends them right; the child predicted
#' on the high side is the majority class of that branch.
#' @keywords internal
#' @noRd
.key_rows <- function(fit, gene_of, classes, split_node, split_markers) {
  frame <- fit$frame
  questions <- which(frame$var != "<leaf>")
  if (length(questions) == 0) return(NULL)
  node_ids <- as.integer(rownames(frame))
  importance <- fit$variable.importance

  rows <- list()
  split_row <- 1
  for (q in questions) {
    n_rows <- 1 + frame$ncompete[q] + frame$nsurrogate[q]
    block <- fit$splits[split_row:(split_row + n_rows - 1), , drop = FALSE]
    split_row <- split_row + n_rows
    if (rownames(block)[1] != as.character(frame$var[q])) {
      stop("Internal error: rpart split table does not match its tree.", call. = FALSE)
    }
    keep <- c(1, if (frame$nsurrogate[q] > 0) seq(frame$ncompete[q] + 2, n_rows))

    for (i in keep) {
      safe_name <- rownames(block)[i]
      high_node <- if (block[i, "ncat"] < 0) 2 * node_ids[q] + 1 else 2 * node_ids[q]
      high_side <- classes[frame$yval[match(high_node, node_ids)]]
      gene <- gene_of[[safe_name]]
      rows[[length(rows) + 1]] <- data.frame(
        split = split_node, rpart_node = node_ids[q], role = if (i == 1) "key" else "backup",
        gene = gene, threshold = signif(block[i, "index"], 4),
        rule = paste0(gene, " >= ", signif(block[i, "index"], 3), " -> ", high_side),
        high_side = high_side, agreement = if (i == 1) NA_real_ else block[i, "adj"],
        importance = if (safe_name %in% names(importance)) unname(importance[[safe_name]]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  key <- do.call(rbind, rows)
  cbind(key, .ranking_values(key, split_markers))
}

#' Ranking values of each key/backup gene: from its marker row for the high-side child,
#' or (if it is not a marker of that child) its best-AUC marker row in the split
#' @keywords internal
#' @noRd
.ranking_values <- function(key, split_markers) {
  value_cols <- c("rank", "auc", "auc_vs_all", "p_val", "p_val_adj", "avg_log2FC", "pct.1", "pct.2",
                  "pct_diff", "mean_child", "mean_siblings", "also_marks_parent",
                  "also_marks_children", "pattern")
  picked <- lapply(seq_len(nrow(key)), function(i) {
    own <- split_markers[split_markers$gene == key$gene[i], , drop = FALSE]
    match_row <- own[own$child == key$high_side[i], , drop = FALSE]
    if (nrow(match_row) == 0) match_row <- own[which.max(own$auc), , drop = FALSE]
    data.frame(marker_of = match_row$child, match_row[, value_cols], stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, picked)
  rownames(out) <- NULL
  out
}

#' Empty key table with the final columns
#' @keywords internal
#' @noRd
.empty_key_table <- function(markers) {
  cbind(data.frame(split = character(), rpart_node = integer(), role = character(), gene = character(),
                   threshold = numeric(), rule = character(), high_side = character(),
                   agreement = numeric(), importance = numeric(), marker_of = character(),
                   stringsAsFactors = FALSE),
        markers[0, c("rank", "auc", "auc_vs_all", "p_val", "p_val_adj", "avg_log2FC", "pct.1", "pct.2",
                     "pct_diff", "mean_child", "mean_siblings", "also_marks_parent",
                     "also_marks_children", "pattern")])
}

#' One summary row per split
#' @keywords internal
#' @noRd
.summarise_splits <- function(split_nodes, tree, per_split, keys, markers) {
  rows <- lapply(split_nodes, function(s) {
    info <- per_split[[s]]
    k <- keys[[s]]
    n_key <- if (is.null(k$key)) 0L else length(unique(k$key$gene[k$key$role == "key"]))
    status <- if (length(info$tested) < 2) "not tested (fewer than 2 children with enough cells)"
              else if (sum(markers$split == s) == 0) "no markers passed the filters"
              else if (n_key == 0) "no key (no question improved on chance)"
              else "ok"
    data.frame(split = s, layer = tree$layer[tree$node == s],
               children = paste(info$children, collapse = ", "),
               skipped_children = paste(info$skipped, collapse = ", "),
               n_cells = sum(lengths(info$cells)), n_markers = sum(markers$split == s),
               n_key_genes = n_key, cv_accuracy = k$cv_accuracy, train_accuracy = k$train_accuracy,
               status = status, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    out <- data.frame(split = character(), layer = integer(), children = character(),
                      skipped_children = character(), n_cells = integer(), n_markers = integer(),
                      n_key_genes = integer(), cv_accuracy = numeric(), train_accuracy = numeric(),
                      status = character(), stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  out
}

#' Resolution tree drawn top-down, with each split's key genes under it
#' @keywords internal
#' @noRd
.plot_resolution_tree <- function(tree, key, sample_name) {
  children_of <- split(tree$node, tree$parent)
  x <- stats::setNames(numeric(nrow(tree)), tree$node)
  next_x <- 0
  place <- function(node) {                     # leaves left to right, parents centred above
    kids <- children_of[[node]]
    if (length(kids) == 0) {
      next_x <<- next_x + 1
      x[[node]] <<- next_x
    } else {
      for (kid in kids) place(kid)
      x[[node]] <<- mean(x[kids])
    }
  }
  place("root")

  nodes <- data.frame(node = tree$node, x = x[tree$node], y = -tree$layer,
                      n_cells = tree$n_cells, label = ifelse(tree$node == "root", "all cells", tree$cluster),
                      purity = ifelse(tree$low_purity, "low purity", "ok"), stringsAsFactors = FALSE)
  child_rows <- tree$node != "root"
  edges <- data.frame(x = x[tree$parent[child_rows]], xend = x[tree$node[child_rows]],
                      y = -tree$layer[child_rows] + 1, yend = -tree$layer[child_rows])

  key_genes <- if (nrow(key) > 0) tapply(key$gene[key$role == "key"], key$split[key$role == "key"],
                                         function(g) paste(unique(g), collapse = ", ")) else character()
  key_labels <- data.frame(x = x[names(key_genes)], y = -tree$layer[match(names(key_genes), tree$node)],
                           label = unname(key_genes), stringsAsFactors = FALSE)

  layer_labels <- unique(data.frame(y = -tree$layer[child_rows],
                                    label = ifelse(is.na(tree$resolution[child_rows]),
                                                   tree$res_col[child_rows],
                                                   paste("res", tree$resolution[child_rows]))))

  p <- ggplot() +
    geom_segment(data = edges, aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
                 colour = "grey60") +
    geom_segment(data = edges, aes(x = .data$xend, xend = .data$xend, y = .data$y, yend = .data$yend),
                 colour = "grey60") +
    geom_point(data = nodes, aes(x = .data$x, y = .data$y, size = .data$n_cells, colour = .data$purity)) +
    geom_text(data = nodes, aes(x = .data$x, y = .data$y, label = .data$label), vjust = -1.5, size = 3) +
    scale_colour_manual(values = c("ok" = "grey20", "low purity" = "red3"), name = "Purity") +
    scale_size_area(max_size = 6, name = "Cells") +
    scale_y_continuous(breaks = c(0, layer_labels$y), labels = c("root", layer_labels$label),
                       expand = expansion(mult = 0.15)) +
    labs(title = paste(sample_name, "- markers along the resolution tree"),
         subtitle = "Blue: key genes of each split", x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid = element_blank())

  if (nrow(key_labels) > 0) {
    p <- p + geom_text(data = key_labels, aes(x = .data$x, y = .data$y, label = .data$label),
                       vjust = 1.8, size = 3, colour = "steelblue4", fontface = "italic")
  }
  p
}

#' Dot plot of all key genes across the clusters of the highest resolution
#' @keywords internal
#' @noRd
.plot_key_dotplot <- function(seurat_obj, key, group_col, assay, sample_name) {
  genes <- unique(key$gene[key$role == "key"])
  if (length(genes) == 0) return(NULL)
  suppressWarnings(suppressMessages(
    Seurat::DotPlot(seurat_obj, features = genes, group.by = group_col, assay = assay)
  )) +
    Seurat::RotatedAxis() +
    labs(title = paste(sample_name, "- key genes"), subtitle = group_col)
}

#' Save the tables (CSV) and plots (JPEG)
#' @keywords internal
#' @noRd
.save_hierarchical_markers <- function(out, output_dir, sample_name, say) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  stem <- file.path(output_dir, paste0(sanitize(sample_name), "_HierMarkers_",
                                       format(Sys.time(), "%Y%m%d_%H%M%S")))
  for (table_name in c("tree", "splits", "markers", "key")) {
    utils::write.csv(out[[table_name]], paste0(stem, "_", table_name, ".csv"), row.names = FALSE)
  }
  n_leaves <- sum(out$tree$n_children == 0)
  ggsave(paste0(stem, "_tree.jpg"), out$plots$tree, width = min(40, max(8, 0.6 * n_leaves)),
         height = 2 + 1.5 * max(out$tree$layer), dpi = 300, limitsize = FALSE, bg = "white")
  if (!is.null(out$plots$dotplot)) {
    n_genes <- length(unique(out$key$gene[out$key$role == "key"]))
    ggsave(paste0(stem, "_dotplot.jpg"), out$plots$dotplot, width = min(40, max(6, 0.35 * n_genes + 3)),
           height = 6, dpi = 300, limitsize = FALSE, bg = "white")
  }
  say("Saved hierarchical markers to: ", output_dir)
}
