# Generate a Gene Marker DotPlot with Aligned Dendrogram

Filters a differential expression data frame for target genes (using
either a specific vector of genes or a regular expression), calculates
average expression metrics, and performs Ward.D2 hierarchical clustering
on the specified condition/identity column. It generates a customized
dot plot (with clusters labeled directly on the points) and aligns it
side-by-side with the hierarchical dendrogram using `patchwork`. The
individual components and the combined layout are automatically saved to
disk.

## Usage

``` r
plot_gene_markers_with_dendro(
  df,
  genes,
  id_col = "short_id",
  output_dir = "Markers_Plots",
  plot_title = "Marker Genes Plot",
  dendro_side = "right",
  label_space = 1,
  width_scale = 0.5,
  height_scale = 0.5,
  min_width = 5,
  min_height = 4
)
```

## Arguments

- df:

  A data frame containing differential expression results. Must contain
  `gene`, `cluster`, `avg_log2FC`, `diff`, and the column specified in
  `id_col`.

- genes:

  Character vector or string. Can be a specific vector of gene names
  (e.g., `c("CD4", "CD8A")`) or a single regular expression string
  (e.g., `"^TR[AB][VD]"`).

- id_col:

  Character. The name of the column in `df` representing the
  experimental condition or identity to cluster on the y-axis. Default
  is "short_id".

- output_dir:

  Character. Directory path where the generated PNGs will be saved.
  Default is "Markers_Plots".

- plot_title:

  Character. The title displayed at the top of the dot plot. Default is
  "Marker Genes Plot".

- dendro_side:

  Character. Which side of the dot plot to attach the dendrogram.
  Options are "left" or "right". Default is "right".

- label_space:

  Numeric. The horizontal padding multiplier for the text labels on the
  dendrogram to prevent truncation. Default is 1.

- width_scale:

  Numeric. A scaling multiplier to dynamically calculate plot width
  based on the number of genes. Default is 0.5.

- height_scale:

  Numeric. A scaling multiplier to dynamically calculate plot height
  based on the number of conditions (`id_col`). Default is 0.5.

- min_width:

  Numeric. The absolute minimum width (in inches) for the saved plots.
  Default is 5.

- min_height:

  Numeric. The absolute minimum height (in inches) for the saved plots.
  Default is 4.

## Value

Invisibly returns a list containing three `patchwork`/`ggplot` objects:

- `dotplot`: The gene expression dot plot.

- `dendrogram`: The hierarchical clustering dendrogram.

- `combined`: The side-by-side aligned layout of both plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'marker_results' is a dataframe combining FindMarkers output
# across multiple datasets or integration methods.

# Example 1: Using a specific list of T cell exhaustion and activation markers
plot_gene_markers_with_dendro(
  df = marker_results,
  genes = c("PDCD1", "HAVCR2", "LAG3", "TOX", "IFNG", "GZMB"),
  id_col = "integration_method",
  plot_title = "Exhaustion Markers Across Integrations",
  dendro_side = "left"
)

# Example 2: Using a Regex search to pull all Killer-cell immunoglobulin-like receptors
plot_gene_markers_with_dendro(
  df = marker_results,
  genes = "^KIR",
  id_col = "Condition",
  output_dir = "NK_Cell_Analysis"
)
} # }
```
