# Compare Clustering Methods at a Shared Number of Clusters

Finds the numbers of clusters that *every* clustering method (pipeline)
produced somewhere in its resolution sweep, so methods can be compared
at the same number of clusters rather than at the same resolution value.
For each method and each shared cluster number it takes the highest
resolution giving that number, plots the size of the smallest cluster
against the number of clusters, saves the plot, and returns the
resolutions for one target number.

## Usage

``` r
analyze_and_save_clusters(
  seurat_obj,
  output_dir = ".",
  plot_prefix = "min_cluster_size_plot",
  target_cluster_num = NULL,
  pipelines = NULL,
  min_cells_line = 30
)
```

## Arguments

- seurat_obj:

  A `Seurat` object whose metadata holds the cluster columns of one or
  more resolution sweeps.

- output_dir:

  Character. Directory where the plot is saved (created if missing).
  Default is `"."`.

- plot_prefix:

  Character. File-name prefix of the saved PNG; a timestamp is appended.
  Default is `"min_cluster_size_plot"`.

- target_cluster_num:

  Numeric or `NULL`. The number of clusters to report. Must be one of
  the shared cluster numbers. If `NULL` (default), you are asked in an
  interactive session; otherwise (or on invalid input) the highest
  shared number is used.

- pipelines:

  Character vector or `NULL`. Only compare these methods (names as
  described above). `NULL` (default) uses every method found.

- min_cells_line:

  Numeric. Height of the dashed red reference line on the plot, a
  warning level for clusters that are too small to interpret. Default is
  30.

## Value

Invisibly, a tibble with one row per method for the target number of
clusters: `Pipeline`, `Resolution` (the highest resolution giving that
number), `Num_Clusters` and `Min_Cluster_Size`. It is also printed.
Attributes: `"plot"` (the `ggplot`), `"plot_data"` (the values for every
shared cluster number) and `"file"` (path of the saved PNG). The plot is
also printed in interactive sessions.

## Details

Resolution columns are recognised by Seurat's naming,
`<graph>_res.<resolution>` (e.g. from
[`ClusterAndUMAP()`](https://chomchonu.github.io/OHMmy/reference/ClusterAndUMAP.md)
or
[`Seurat::FindClusters()`](https://satijalab.org/seurat/reference/FindClusters.html)).
The method name is the part before `.log_snn`/`.sct_snn` for the OHMmy
naming (`"log.fastmnn.log_snn_res.1"` -\> `"log.fastmnn"`,
`"sct.fastmnn.sct_snn_res.1"` -\> `"sct.fastmnn"`), and the graph name
otherwise (`"RNA_snn_res.0.5"` -\> `"RNA_snn"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Two resolution sweeps, e.g. LogNormalize and SCTransform with FastMNN
res <- analyze_and_save_clusters(
  seurat_obj         = my_seurat,     # has log.fastmnn.log_snn_res.* and sct.fastmnn.sct_snn_res.*
  output_dir         = "Plots_cluster_choice",
  target_cluster_num = 12
)
res                       # resolution of each method that gives 12 clusters
attr(res, "plot")         # smallest cluster size vs number of clusters

# Compare only some of the methods
analyze_and_save_clusters(my_seurat, pipelines = c("log.fastmnn", "log.har"))
} # }
```
