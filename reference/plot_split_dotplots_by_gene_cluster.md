# Generate Split DotPlots and Global Dendrogram for Marker Genes

Filters a differential expression results data frame for specific genes
(via list or regex), performs hierarchical clustering on the specified
identity column based on average log2FC, and generates a global
dendrogram. To prevent unreadable, overcrowded x-axes, the target
gene-cluster combinations are automatically split into multiple chunked
dotplots. All plots are automatically sized and saved directly to the
specified output directory.

## Usage

``` r
plot_split_dotplots_by_gene_cluster(
  df,
  gene_list = NULL,
  gene_regex = NULL,
  id_col = "integration_method",
  output_dir = "Markers_Plots/Split",
  n_splits = 4,
  plot_title_prefix = "Gene DotPlot",
  width_scale = 0.3,
  height_scale = 0.5,
  min_width = 5,
  min_height = 5,
  label_space = 1
)
```

## Arguments

- df:

  A data frame of differential expression results. Must contain the
  columns `gene`, `cluster`, `avg_log2FC`, `diff`, and the column
  specified in `id_col`.

- gene_list:

  Character vector. A specific list of genes to filter and plot. Default
  is NULL.

- gene_regex:

  Character string. A regular expression to match target genes (e.g.,
  `"^TR[AB][VD]"` for T cell receptor genes or "^KLR" for Killer-cell
  receptors). Default is NULL.

- id_col:

  Character string. The name of the column in `df` to map to the y-axis
  and use for hierarchical clustering. Default is "integration_method".

- output_dir:

  Character. Directory path where the generated PNGs will be saved.
  Default is "Markers_Plots/Split".

- n_splits:

  Integer. The number of parts to split the horizontal axis
  (gene-cluster pairs) into. Default is 4.

- plot_title_prefix:

  Character. The base string used for the title of each split plot.
  Default is "Gene DotPlot".

- width_scale:

  Numeric. A scaling multiplier to dynamically calculate plot width
  based on the number of x-axis items. Default is 0.3.

- height_scale:

  Numeric. A scaling multiplier to dynamically calculate plot height
  based on the number of y-axis items. Default is 0.5.

- min_width:

  Numeric. The absolute minimum width (in inches) for the saved plots.
  Default is 5.

- min_height:

  Numeric. The absolute minimum height (in inches) for the saved plots.
  Default is 5.

- label_space:

  Numeric. The spatial padding multiplier for the text labels on the
  global dendrogram. Default is 1.

## Value

Invisibly returns `NULL`. Plots are directly written to disk.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'marker_df' contains combined FindMarkers output across multiple
# integration methods or conditions.

# Example 1: Plot specific T/NK cell markers, split into 2 plots
plot_split_dotplots_by_gene_cluster(
  df = marker_df,
  gene_list = c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "GZMB", "PRF1"),
  id_col = "Condition",
  output_dir = "Results/DotPlots",
  n_splits = 2,
  plot_title_prefix = "Cytotoxic Markers"
)

# Example 2: Plot all genes starting with "KLR" using Regex
plot_split_dotplots_by_gene_cluster(
  df = marker_df,
  gene_regex = "^KLR",
  id_col = "integration_method",
  n_splits = 3
)
} # }
```
