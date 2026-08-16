# Generate Split Clustered DotPlots with Row Dendrograms

An extension of the bidirectionally clustered DotPlot, specifically
designed to handle massive gene lists without overcrowding the x-axis.
It performs global hierarchical clustering on both genes and clusters to
establish the true biological ordering. If the number of filtered genes
exceeds `max_genes_per_plot`, it chunks the x-axis into multiple
separate plots. The global cluster dendrogram (y-axis) is identically
attached to each chunk for easy cross-referencing.

## Usage

``` r
plot_dot_dendro_split(
  seurat_obj,
  meta_col,
  feature_df,
  scale = TRUE,
  pct_threshold = 15,
  output_dir = "Output_R",
  prefix = NULL,
  max_genes_per_plot = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- meta_col:

  Character. The metadata column used to group cells (e.g.,
  "seurat_clusters").

- feature_df:

  A data frame mapping genes to categories. Must contain a `group`
  column (to filter by prefix) and a `feature` column (containing the
  gene names).

- scale:

  Logical. Whether to scale the average expression values (z-score)
  across groups before plotting. Default is TRUE.

- pct_threshold:

  Numeric. The minimum percentage of cells expressing the gene in at
  least one cluster required to retain the gene. Default is 15.

- output_dir:

  Character. Directory path where the generated JPEGs and text file will
  be saved. Default is "Output_R".

- prefix:

  Character. A string used to filter the `group` column in `feature_df`.
  If `NULL`, it automatically derives the prefix by removing "label\_"
  from `meta_col`. Default is NULL.

- max_genes_per_plot:

  Integer. The maximum number of genes to display per plot. If the
  filtered gene list exceeds this number, the output is split into
  multiple chunked plots. Default is NULL (no splitting).

## Value

A list of `patchwork` objects, where each element represents a specific
chunk of the x-axis (dot plot) combined with the global row (cluster)
dendrogram.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming marker_mapping is a dataframe with 'group' and 'feature' columns

# Generate chunked dotplots, forcing a maximum of 25 genes per plot
chunked_plots <- plot_dot_dendro_split(
  seurat_obj = my_seurat,
  meta_col = "T_Cell_States",
  feature_df = marker_mapping,
  scale = TRUE,
  pct_threshold = 10,
  output_dir = "Results/Split_DotPlots",
  prefix = "Exhaustion",
  max_genes_per_plot = 25
)

# View the first chunk directly in R
chunked_plots[[1]]

# View the second chunk
chunked_plots[[2]]
} # }
```
