# Generate a Bidirectionally Clustered DotPlot with Dendrograms

Extracts expression data via Seurat's `DotPlot` function, filters out
lowly expressed genes based on a minimum percentage threshold, and
performs two-way hierarchical clustering on both the features (genes)
and the identities (clusters). It assembles a comprehensive `patchwork`
layout featuring the central dot plot flanked by column and row
dendrograms. The function automatically saves the combined plot,
standalone dendrograms, and a text file of the final clustered gene
order.

## Usage

``` r
plot_dot_dendro(
  seurat_obj,
  meta_col,
  feature_df,
  scale = TRUE,
  pct_threshold = 15,
  output_dir = "Output_R",
  prefix = NULL
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- meta_col:

  Character. The metadata column used to group cells (e.g.,
  "seurat_clusters" or "CellType").

- feature_df:

  A data frame mapping genes to categories. Must contain a `group`
  column (to filter by prefix) and a `feature` column (containing the
  gene names).

- scale:

  Logical. Whether to scale the average expression values (z-score)
  across groups before plotting and clustering. Default is TRUE.

- pct_threshold:

  Numeric. The minimum percentage of cells expressing the gene in at
  least one cluster required to retain the gene. Default is 15.

- output_dir:

  Character. Directory path where the generated JPEGs and text file will
  be saved. Default is "Output_R".

- prefix:

  Character. A string used to filter the `group` column in `feature_df`.
  If `NULL`, it automatically attempts to derive the prefix by removing
  "label\_" from `meta_col`. Default is NULL.

## Value

A `patchwork` object representing the final combined layout (dot plot +
column dendrogram + row dendrogram).

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a mock mapping dataframe for subset markers
marker_mapping <- data.frame(
  group = c("Effector_T", "Effector_NK", "Naive_T"),
  feature = c("GZMB", "NKG7", "CCR7")
)

# Generate the bi-clustered dotplot
clustered_plot <- plot_dot_dendro(
  seurat_obj = my_seurat,
  meta_col = "Subset_Labels",
  feature_df = marker_mapping,
  scale = TRUE,
  pct_threshold = 10,
  output_dir = "Results/DotPlots",
  prefix = "Effector" # Will selectively plot and cluster the first two genes
)

# Display the assembled plot in the R viewer
print(clustered_plot)
} # }
```
