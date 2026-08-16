# Generate a Multi-Metadata Bidirectionally Clustered DotPlot

An advanced extension of the clustered DotPlot function designed to
handle and concatenate multiple metadata groupings simultaneously. It
iterates over a vector of metadata columns, aggregates the expression
data using Seurat's `DotPlot` engine, and binds them into a single
comprehensive matrix. It then performs 2D hierarchical clustering across
all combined groups and genes. The function features a modular saving
system, allowing users to independently export the combined plot,
individual dendrograms, or the clustered gene list.

## Usage

``` r
plot_dot_dendro_multi(
  seurat_obj,
  meta_cols,
  feature_df,
  scale = TRUE,
  pct_threshold = 15,
  output_dir = "Output_R",
  prefix = NULL,
  base_size = 12,
  range = c(0.25, 6),
  save_options = c("combined", "col_dend", "row_dend", "gene_list")
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- meta_cols:

  Character vector. A list of metadata columns to evaluate and combine
  (e.g., `c("Condition_A", "Condition_B")` or
  `c("res.0.5", "res.1.0")`).

- feature_df:

  A data frame mapping genes to categories. Must contain a `group`
  column and a `feature` column.

- scale:

  Logical. Whether to use scaled average expression values (z-scores)
  for plotting and clustering. Default is TRUE.

- pct_threshold:

  Numeric. The minimum percentage of cells expressing the gene in at
  least one cluster/group combination required to retain the gene.
  Default is 15.

- output_dir:

  Character. Directory path where the generated outputs will be saved.
  Default is "Output_R".

- prefix:

  Character. A string used to filter the `group` column in `feature_df`.
  If `NULL`, it attempts to derive it from the metadata column names.
  Default is NULL.

- base_size:

  Numeric. The base font size for the `ggplot2` theme. Default is 12.

- range:

  Numeric vector of length 2. The minimum and maximum point sizes for
  the dots (`scale_radius`). Default is `c(0.25, 6)`.

- save_options:

  Character vector. Defines which specific outputs to generate and save.
  Options include "combined", "col_dend", "row_dend", and "gene_list".
  Default includes all four.

## Value

Invisibly returns a `tibble` containing the filtered, long-format data
used to generate the plots (useful for custom downstream plotting in
`ggplot2`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming marker_mapping contains your 'group' and 'feature' columns

# Compare exhaustion markers across three different metadata cluster resolutions
plot_data <- plot_dot_dendro_multi(
  seurat_obj = my_seurat,
  meta_cols = c("SCT_snn_res.0.4", "SCT_snn_res.0.8", "SCT_snn_res.1.2"),
  feature_df = marker_mapping,
  pct_threshold = 10,
  prefix = "Exhaustion",
  output_dir = "Results/Multi_Resolution",
  save_options = c("combined", "gene_list") # Skip saving standalone dendrograms
)

# The function invisibly returns the data, which you can inspect:
head(plot_data)
} # }
```
