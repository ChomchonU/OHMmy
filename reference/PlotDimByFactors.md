# Generate Dimensionality Reduction Plots by Metadata Factors

This function iterates through a specified list of metadata columns
(factors) in a Seurat object and generates a dimensionality reduction
plot (e.g., UMAP) for each. It automatically handles filename
sanitization, saves the plots to a designated directory, and provides
console progress updates.

## Usage

``` r
PlotDimByFactors(
  seurat_obj,
  factors,
  sample_name = "Sample",
  reduction = "umap.har",
  raster = FALSE,
  label = TRUE,
  repel = TRUE,
  output_dir = "Plots_umap",
  plot_format = "jpg",
  width = 10,
  height = 10,
  dpi = 300,
  add_timestamp = TRUE,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data and the specified
  reduction.

- factors:

  Character vector. The names of the metadata columns to group the cells
  by (e.g., c("seurat_clusters", "Phase")).

- sample_name:

  Character. The name of the biological sample, used for plot titles and
  file naming. Default is "Sample".

- reduction:

  Character. The dimensionality reduction to visualize (e.g., "umap",
  "umap.har" for Harmony). Default is "umap.har".

- raster:

  Logical. Whether to rasterize the points (useful for massive datasets
  to speed up plotting). Default is FALSE.

- label:

  Logical. Whether to label the clusters/groups directly on the plot.
  Default is TRUE.

- repel:

  Logical. Whether to repel the labels to prevent overlapping. Default
  is TRUE.

- output_dir:

  Character. Directory path where the generated plots will be saved.
  Default is "Plots_umap".

- plot_format:

  Character. The file format for the saved plots ("jpg", "png", or
  "pdf"). Default is "jpg".

- width:

  Numeric. The width of the saved plot. Default is 10.

- height:

  Numeric. The height of the saved plot. Default is 10.

- dpi:

  Numeric. The resolution of the saved plots. Default is 300.

- add_timestamp:

  Logical. Whether to append the current date and time to the saved
  filenames. Default is TRUE.

- verbose:

  Logical. Whether to print progress messages, display a progress bar,
  and print the plots to the console. Default is TRUE.

## Value

Invisibly returns a named list of the generated `ggplot` objects, where
names correspond to the provided `factors`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot UMAPs colored by multiple different metadata factors
my_plots <- PlotDimByFactors(
  seurat_obj = my_seurat,
  factors = c("seurat_clusters", "CellType", "Condition"),
  sample_name = "Donor_1",
  reduction = "umap",
  plot_format = "png"
)

# The function saves the plots to disk, but you can also view one in R:
my_plots$CellType
} # }
```
