# Cut the Dendrograms of `plot_dot_dendro_multi()` into Clusters

Rebuilds the exact hierarchical clustering behind the dendrograms drawn
by
[`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
and cuts it into clusters of groups (rows), genes (columns), or both.
Because natural clusters often sit at different heights of a dendrogram,
the default methods cut each branch at its own height instead of using a
single horizontal cut, so both the number of clusters and the cut
heights are derived from the tree:

- `"inconsistency"` (default): a branch is kept as one cluster while
  none of its merges is unusually high compared with the merges directly
  below it (inconsistency coefficient, depth 2). Each branch is
  therefore cut relative to its own spread, which suits clusters at
  clearly different heights.

- `"dynamic"`: the hybrid dynamic tree cut of dynamicTreeCut (Langfelder
  et al., 2008), which uses both the tree shape and the distance matrix.
  It works best when clusters have broadly similar heights.

- `"fixed"`: an ordinary horizontal cut at `k` clusters or height `h`.

## Usage

``` r
cluster_dot_dendro(
  dot_df,
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
  output_dir = NULL
)
```

## Arguments

- dot_df:

  The data frame returned by
  [`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md).

- what:

  Character. Cut the tree of the groups (`"rows"`, default), of the
  genes (`"genes"`), or both (`"both"`).

- method:

  Character. The cutting method: `"inconsistency"` (default),
  `"dynamic"` or `"fixed"`. See Description.

- deepSplit:

  Integer from 0 to 4, or `"auto"` (default). Dynamic method only: how
  readily branches are split into sub-clusters (higher gives more,
  smaller clusters). `"auto"` picks the value with the best average
  silhouette width.

- cutoff:

  Numeric or `"auto"` (default). Inconsistency method only: branches
  whose merges all have an inconsistency coefficient at or below
  `cutoff` form one cluster. `"auto"` picks the value with the best
  average silhouette width.

- minClusterSize:

  Integer. The smallest number of leaves that can form a cluster with
  the dynamic and inconsistency methods; smaller groups are left
  unassigned. Default is 2. Not used by the fixed method.

- assign_all:

  Logical. If `FALSE` (default), leaves that do not clearly belong to
  any cluster stay unassigned (`NA`). If `TRUE`, each of them is
  assigned to the cluster with the smallest average distance.

- k, h:

  Numeric. Fixed method only: the number of clusters (`k`) or the height
  at which to cut (`h`). Give exactly one.

- seurat_obj:

  Optional Seurat object. If given (together with `meta_col`), the row
  clusters are written back to the cells as a new metadata column.

- meta_col:

  Character. The metadata column (one of the `meta_cols` used in
  [`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md))
  whose groups are mapped to cells. Each cell receives the cluster of
  its group's row; cells of unassigned groups get `NA`.

- new_col:

  Character. The name of the new metadata column. Default is
  `"dendro_cluster"`.

- output_dir:

  Character. Optional directory in which the cluster table (CSV) and the
  coloured dendrogram(s) (JPEG) are saved. Default is `NULL` (nothing is
  saved).

## Value

A list with:

- `clusters`: a data frame with one row per leaf: `type` (`"rows"` or
  `"genes"`), `label`, `id` and `var` (rows only), `cluster` (factor,
  `NA` = unassigned), `dendro_order` (position in the dendrogram),
  `cluster_height` (height of the merge joining the cluster's leaves)
  and `cut_height` (the height at which the cluster's branch was cut:
  midway between `cluster_height` and the merge joining it to the rest
  of the tree).

- `k`: a named integer vector with the number of clusters found per
  tree.

- `selection`: a named list with, per tree, the candidate settings and
  their scores (for `"auto"` settings), or `NULL`.

- `plot`: a named list of `ggplot` dendrograms coloured by cluster, with
  a dashed mark at each cluster's cut height; unassigned leaves are
  grey.

- `hclust`: a named list of the `hclust` trees that were cut.

- `parameters`: the method and the settings that were used.

- `seurat_obj`: the updated Seurat object (only when `seurat_obj` was
  given).

## Details

The trees are taken from the `"row_hclust"` and `"col_hclust"`
attributes that
[`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
attaches to its output. For older results without these attributes, the
same trees are rebuilt from the `avg.exp.scaled` values (Euclidean
distance, complete linkage).

With `deepSplit = "auto"` (dynamic method) or `cutoff = "auto"`
(inconsistency method), every candidate setting is tried and the one
with the highest average silhouette width is chosen. Unassigned leaves
count as a silhouette of 0, so leaving many leaves out is not rewarded.
All candidates are reported in `selection`.

Note that
[`DotPlot()`](https://satijalab.org/seurat/reference/DotPlot.html)
z-scores expression separately within each metadata column, so when
several columns are combined, groups are compared by expression pattern
rather than on one absolute scale. This is the same basis as the plotted
dendrograms.

## References

Langfelder P, Zhang B, Horvath S (2008). Defining clusters from a
hierarchical cluster tree: the Dynamic Tree Cut package for R.
*Bioinformatics* 24(5):719-720.

## See also

[`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dot_df <- plot_dot_dendro_multi(seurat_obj, meta_cols = c("cell_type", "seurat_clusters"),
                                feature_df = marker_df, prefix = "PBMC")

# Cut the group dendrogram at variable heights
res <- cluster_dot_dendro(dot_df)
res$clusters
res$plot$rows

# Genes and groups with the dynamic tree cut, and write the group clusters to the cells
res <- cluster_dot_dendro(dot_df, what = "both", method = "dynamic",
                          seurat_obj = seurat_obj, meta_col = "seurat_clusters")
table(res$seurat_obj$dendro_cluster)
} # }
```
