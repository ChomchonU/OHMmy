# Track Markers Down the Resolution Tree to Annotate Clusters

Follows the markers found by
[`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md)
from low to high resolution, to help name clusters. For every gene it
shows where it first marks a cluster, how far down that branch it stays
detected and where it is lost. For every cluster it suggests a label
made of the top marker at each split on the way from the root, e.g.
`"CST3+ > FCGR3A+"`, which you translate into a cell type. Labels are
suggestions only: you make the final call.

## Usage

``` r
track_cluster_markers(
  seurat_obj,
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
  ...
)
```

## Arguments

- seurat_obj:

  A `Seurat` object (the one used for `hier`).

- hier:

  Optional result of
  [`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md)
  on `seurat_obj`. If `NULL` (default), it is run with `res_cols` or
  `pipeline` and `...`.

- res_cols, pipeline:

  Resolution columns or integration method, as in
  [`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md).
  Used only when `hier` is `NULL`.

- n_track:

  Integer. Top markers per child and split to track. Default 10.

- min_label_auc:

  Numeric. Minimum AUC against the siblings for a label gene. Default
  0.8.

- min_global_auc:

  Numeric. Minimum AUC against all other cells for a gene to be
  `global`. Default 0.7.

- n_own:

  Integer. Own markers listed per cluster. Default 5.

- sample_name:

  Character. Used in plot titles and file names. Default "Sample".

- output_dir:

  Character or `NULL`. If given, the tables (CSV) and plots (JPEG, sized
  from their content) are saved there. Default `NULL`.

- verbose:

  Logical. Print progress messages. Default `TRUE`.

- ...:

  Passed to
  [`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md)
  when `hier` is `NULL` (e.g. `assay`, `min_cells`).

## Value

A list with:

- `annotation`: one row per cluster: `node`, `layer`, `resolution`,
  `cluster`, `n_cells`, `is_final`, `path_label`, `status`
  (`"suggested"` / `"ambiguous"`), `reason`, `confidence` (lowest AUC vs
  siblings of the label genes on the path), `min_auc_vs_all` (lowest AUC
  vs all of the label genes), `lineage_markers` (lineage markers
  inherited from the splits above), `own_markers` (top global markers at
  the cluster's own split), `relative_markers` (top markers there that
  are only relative), `lowest_purity` (on the path).

- `gene_paths`: one row per gene and start: `gene`, `start_node`,
  `start_resolution`, `start_split`, `kind` (the marker pattern at the
  start), `auc_at_start`, `global_at_start`, `marked_again_below`,
  `n_final`, `n_final_on`, `final_on`, `stops_at`.

- `node_profile`: one row per cluster and tracked gene: `node`, `layer`,
  `resolution`, `gene`, `mean`, `pct`, `pct_rest`, `log2fc_vs_all`,
  `auc_vs_all`, `global`, `auc_vs_siblings`, `is_marker`.

- `plots`: `tree` (tree with the label gene of every split), `heatmap`
  (AUC against all cells of the tracked genes at every cluster, markers
  dotted) and `dotplot` (label genes across the final clusters in tree
  order).

- `hier`: the
  [`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md)
  result used.

- `parameters`.

## Details

**Tracked genes.** The top `n_track` markers (by AUC against the
siblings) of every child at every split.

**Profile.** For each tracked gene at every cluster of every resolution:
mean expression, share of cells detected (`pct`) and outside the cluster
(`pct_rest`), log2 fold change and AUC against all other cells, and AUC
against its siblings (the other children of its parent; all their
cells). `global` is `auc_vs_all >= min_global_auc`: the gene is high in
the cluster compared with the whole object, not only compared with its
siblings ("relative").

**Gene paths.** A gene's path starts at a cluster where it is a marker
and no cluster above it on the branch is marked by it. Below the start,
a cluster is "on" when the gene is detected in at least `min_pct` of its
cells (the detection filter of
[`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md));
`stops_at` lists the first clusters where it is no longer on.

**Labels.** For each split on the way from the root, the label gene is
the highest-ranked global marker among the top `n_track` of that child
(if none is global, the top marker, flagged relative). A label is
`"ambiguous"` when, at any split on its path, there is no marker (or the
child was too small to test), the label gene's AUC against the siblings
is below `min_label_auc`, the label gene is only relative, it is lower
in the cluster than in the rest of the object (AUC against all \< 0.5),
or a cluster on the path has low purity. `reason` says which.

Statistics use the same cells and clusters, so they describe the data
rather than prove the labels; check the labels against known biology.

## See also

[`find_hierarchical_markers()`](https://chomchonu.github.io/OHMmy/reference/find_hierarchical_markers.md)

## Examples

``` r
if (FALSE) { # \dontrun{
hier <- find_hierarchical_markers(my_seurat, pipeline = "log.rpca")
ann <- track_cluster_markers(my_seurat, hier = hier, output_dir = "Track_markers")
ann$annotation[, c("node", "path_label", "status", "reason")]
subset(ann$gene_paths, gene == "GZMK")
ann$plots$heatmap
} # }
```
