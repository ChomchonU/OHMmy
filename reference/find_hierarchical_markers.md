# Find Markers Along a Resolution Tree and Build a Dichotomous Key

Tracks clusters from low to high resolution and finds, at every point
where a cluster splits, the genes that distinguish its sub-clusters. Low
resolutions give broad markers of major populations; higher resolutions
give finer markers of sub-populations. For every split, a short decision
tree (rpart) picks the genes and thresholds that sort cells into the
sub-clusters: a dichotomous key.

## Usage

``` r
find_hierarchical_markers(
  seurat_obj,
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
  max_candidate_genes = 50,
  max_depth = NULL,
  n_backup = 3,
  cv_folds = 10,
  seed = 42,
  sample_name = "Sample",
  output_dir = NULL,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A `Seurat` object with normalised data and cluster columns at several
  resolutions.

- res_cols:

  Character vector of metadata columns holding the clusters, e.g.
  `c("RNA_snn_res.0.2", "RNA_snn_res.0.5", "RNA_snn_res.1")`. Ordered by
  the resolution in their name when it can be read, otherwise used as
  given. Give this or `pipeline`.

- pipeline:

  Character. Instead of `res_cols`: use every `<graph>_res.<resolution>`
  column of this method, named as in
  [`analyze_and_save_clusters()`](https://chomchonu.github.io/OHMmy/reference/analyze_and_save_clusters.md)
  (e.g. `"log.fastmnn"` or `"RNA_snn"`).

- assay:

  Character. Assay used for expression. Default is `"RNA"`.

- padj_max, min_log2fc, min_pct, min_pct_diff:

  Numeric marker filters: adjusted p-value below `padj_max` (default
  0.05), `avg_log2FC` at least `min_log2fc` (default 0.5), detected in
  at least `min_pct` of the child's cells (default 0.5), and detection
  difference `pct.1 - pct.2` at least `min_pct_diff` (default 0.3).

- min_purity:

  Numeric. Children below this purity are flagged and their markers
  labelled `"unstable"` (they are still tested). Default is 0.8.

- min_cells:

  Integer. Children with fewer cells in a split are skipped. Default 20.

- specific_auc:

  Numeric. AUC against all other cells needed for `"specific"`. Default
  is 0.8.

- max_candidate_genes:

  Integer. Top-ranked markers per child offered to `rpart`. Default is
  50.

- max_depth:

  Integer or `NULL`. Maximum number of questions on any path of a
  split's key before pruning. `NULL` (default) uses the number of
  children of the split (at most 30).

- n_backup:

  Integer. Backup (surrogate) genes reported per question. Default 3.

- cv_folds:

  Integer. Cross-validation folds for pruning. Default is 10.

- seed:

  Integer. Seed for the cross-validation folds (the global random state
  is restored afterwards). Default is 42.

- sample_name:

  Character. Used in plot titles and file names. Default "Sample".

- output_dir:

  Character or `NULL`. If given, the tables (CSV) and plots (JPEG) are
  saved there. Default is `NULL` (nothing saved).

- verbose:

  Logical. Print progress messages. Default is `TRUE`.

## Value

A list with:

- `tree`: one row per cluster (node): `node`, `layer`, `resolution`,
  `res_col`, `cluster`, `parent`, `n_cells`, `purity`, `n_children`,
  `low_purity`.

- `splits`: one row per split: `split`, `layer`, `children`,
  `skipped_children` (too few cells), `n_cells`, `n_markers`,
  `n_key_genes`, `cv_accuracy`, `train_accuracy`, `status`.

- `markers`: every gene passing the filters, per split and child:
  `split`, `child`, `gene`, `rank`, `auc`, `auc_vs_all`, `p_val`,
  `p_val_adj`, `avg_log2FC`, `pct.1`, `pct.2`, `pct_diff`, `mean_child`,
  `mean_siblings`, `also_marks_parent`, `also_marks_children`,
  `pattern`.

- `key`: the key and backup genes of every split: `split`, `rpart_node`,
  `role`, `gene`, `threshold`, `rule` (e.g.
  `"GZMK >= 1.2 -> res0.5_3"`), `high_side` (the child predicted for
  cells at or above the threshold), `agreement` (backup genes: rpart's
  adjusted agreement with the key gene), `importance`, then `marker_of`
  and the ranking columns of `markers` for that gene.

- `confusion`: per split, a table of true vs key-predicted child.

- `plots`: `tree` (the resolution tree with key genes at each split) and
  `dotplot` (key genes across the clusters of the highest resolution).

- `res_cols` and `parameters`.

## Details

**1. Tree.** Resolutions are ordered from low to high. Each cluster is
linked to the cluster of the previous resolution that holds most of its
cells (its parent); that share is its `purity`. Clusterings at different
resolutions are not always nested, so a low purity means the cluster
takes cells from several parents. The root holds all cells, so its
children are the clusters of the lowest resolution. A cluster with a
single child (unchanged at the next resolution) is not a split.

**2. Markers per split.** At each split, every child is compared with
its siblings only, using the parent's cells
([`Seurat::FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html),
Wilcoxon, genes higher in the child). Cells of the parent that moved to
a child of another parent are left out. Children with fewer than
`min_cells` cells are skipped. Genes are kept when
`p_val_adj < padj_max`, `avg_log2FC >= min_log2fc`, `pct.1 >= min_pct`
and `pct.1 - pct.2 >= min_pct_diff`, then ranked by AUC (the probability
that a random child cell expresses the gene more highly than a random
sibling cell; 0.5 = no separation, 1 = perfect), with the detection
difference as tie-break. Clusters were defined on the same data, so
p-values are optimistic: use them as a filter, not as evidence on their
own.

**3. Patterns.** Each marker row gets `also_marks_parent` (the gene also
marks the parent against the parent's siblings), `also_marks_children`
(it also separates the child's own sub-clusters; `NA` for a child that
does not split further) and a `pattern`:

- `"unstable"`: the child (or a cluster it passes through unchanged) has
  a purity below `min_purity`, so the split is not clean across
  resolutions;

- `"specific"`: the child does not split further and the gene also
  separates it from all other cells (`auc_vs_all >= specific_auc`);

- `"lineage"`: the child splits further, but the gene does not separate
  its sub-clusters (shared by all of them);

- `"split"`: any other marker of this split.

**4. Key.** For each split, `rpart` is given the top
`max_candidate_genes` ranked markers of every child, grown to
`max_depth` questions and pruned by cross-validation (1-SE rule: the
smallest tree whose error is within one standard error of the best). The
data therefore decide how many genes each split needs. For every
question, the chosen gene (`role = "key"`) and up to `n_backup`
surrogate genes that sort the cells almost the same way
(`role = "backup"`) are reported with all ranking values.

Expression is read from the `data` layer of `assay`; split Seurat v5
layers are joined on a local copy.

## See also

[`analyze_and_save_clusters()`](https://chomchonu.github.io/OHMmy/reference/analyze_and_save_clusters.md)
to choose resolutions,
[`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
to plot the key genes.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- find_hierarchical_markers(
  seurat_obj = my_seurat,
  res_cols   = c("RNA_snn_res.0.1", "RNA_snn_res.0.3", "RNA_snn_res.0.8"),
  output_dir = "Hierarchical_markers"
)
res$key                                  # key + backup genes with their ranking values
subset(res$markers, pattern == "lineage")
res$plots$tree

# All resolutions of one integration method
res <- find_hierarchical_markers(my_seurat, pipeline = "log.fastmnn")
} # }
```
