# Plot Cell Counts and Proportions Across Clusters and Batches

Generates stacked bar charts to visualize the distribution of cells
between two metadata variables (typically cell clusters and experimental
batches/conditions). It calculates and plots both absolute cell counts
and relative percentages (proportions), assembling them into two
separate dual-panel figures using `patchwork`. The resulting plots are
automatically saved to the specified output directory.

## Usage

``` r
plot_cluster_distributions(
  seurat_obj,
  cluster_col = "label_T3_log_rpca",
  batch_col = "exist_GTS",
  output_dir =
    "C:/Users/ADMIN/Desktop/Dengue_summer/Output_R/Plots_count_dis/CD8_subset/rpca",
  file_prefix = "ind_label"
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- cluster_col:

  Character. The name of the metadata column representing cell clusters
  or identity classes. Default is "label_T3_log_rpca".

- batch_col:

  Character. The name of the metadata column representing batches,
  samples, or conditions (e.g., "exist_GTS"). Default is "exist_GTS".

- output_dir:

  Character. The directory path where the generated JPEGs will be saved.
  Default is a specific local directory path.

- file_prefix:

  Character. A prefix string to append to the saved filenames to help
  identify the analysis run. Default is "ind_label".

## Value

Invisibly returns a list containing two `patchwork` plot objects:

- `proportions`: The combined plot showing relative percentages.

- `counts`: The combined plot showing absolute cell counts.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare CD8 T cell distributions across different disease states or batches
distribution_plots <- plot_cluster_distributions(
  seurat_obj = my_seurat,
  cluster_col = "CD8_Subsets",
  batch_col = "Disease_Status",
  output_dir = "Plots_Distributions",
  file_prefix = "CD8_analysis"
)

# The function saves the JPEGs automatically, but you can also view them in R:
distribution_plots$proportions
distribution_plots$counts
} # }
```
