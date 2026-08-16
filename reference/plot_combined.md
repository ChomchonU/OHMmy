# Generate Paired Feature and Density Plots

This function takes a Seurat object and a named list of features (which
can include both gene names and metadata columns). For each feature, it
generates a side-by-side composite of a standard Seurat `FeaturePlot`
and a custom density plot. The resulting interleaved plots are stitched
together using `patchwork` and automatically saved to disk.

## Usage

``` r
plot_combined(
  seurat_obj,
  cluster_col,
  genes_list,
  sample_name = "Sample",
  reduction = "umap.cca",
  viridis_palette = "viridis",
  output_dir = "Plots_feat_dens",
  save_format = "jpg",
  width = 20,
  dpi = 300,
  add_timestamp = TRUE,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- cluster_col:

  Character. The name of the metadata column to set as the active
  identity (`Idents`).

- genes_list:

  A named list where each element is a character vector of features
  (e.g., genes or metadata columns). The names of the list (e.g.,
  "T_cell_markers") are used for plot titles and file naming.

- sample_name:

  Character. The name of the biological sample, used for plot titles and
  file naming. Default is "Sample".

- reduction:

  Character. The dimensionality reduction to use for plotting (e.g.,
  "umap", "umap.cca"). Default is "umap.cca".

- viridis_palette:

  Character. The viridis color palette to use for the density plots
  (e.g., "viridis", "magma", "plasma"). Default is "viridis".

- output_dir:

  Character. The directory path where the generated plots will be saved.
  Default is "Plots_feat_dens".

- save_format:

  Character. The file format for the saved plots. Options are "jpg",
  "png", or "pdf". Default is "jpg".

- width:

  Numeric. The base width multiplier for the saved plot dimensions.
  Default is 20.

- dpi:

  Numeric. The resolution of the saved plots. Default is 300.

- add_timestamp:

  Logical. Whether to append the current date and time to the saved
  filenames. Default is TRUE.

- verbose:

  Logical. Whether to display progress bars and console messages.
  Default is TRUE.

## Value

A list containing two elements:

- `combined_plots`: A named list of the generated `patchwork` plot
  objects, where names correspond to the names of `genes_list`.

- `output_dir`: A character string of the path where the plots were
  saved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define a named list of features containing both genes and metadata
marker_list <- list(
  T_Cells = c("CD3E", "CD4", "CD8A", "percent.mt"),
  NK_Cells = c("NCAM1", "NKG7", "GNLY")
)

# Run the plotting function
plot_results <- plot_combined(
  seurat_obj = my_seurat,
  cluster_col = "seurat_clusters",
  genes_list = marker_list,
  sample_name = "Donor_1",
  reduction = "umap",
  save_format = "png"
)

# Access a specific plot directly in R without opening the saved file
plot_results$combined_plots$T_Cells
} # }
```
