# Generate Dimensionality Reduction Heatmaps

A robust wrapper around Seurat's `DimHeatmap` function designed to
systematically evaluate sources of heterogeneity across multiple
dimensions. It chunks specified principal components (or other
dimensional reductions) into distinct windows, generates a balanced
heatmap for each window, and safely exports them as high-resolution
JPEGs. It features a built-in progress bar and robust error handling
(`tryCatch`) to prevent graphics device hangs during batch processing on
remote compute nodes.

## Usage

``` r
generate_dimheatmaps(
  seurat_obj,
  sample_name,
  reduction = "pca",
  pc_windows = list(1:10, 11:20, 21:30),
  output_dir = "DimHeatmap_plots",
  nfeatures = 60,
  cells = 500,
  width_in = 15,
  height_in = 45,
  res = 300
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data and the specified
  dimensional reduction.

- sample_name:

  Character. The name of the biological sample, used for console
  messages and prefixing the saved filenames.

- reduction:

  Character. The dimensional reduction to use (e.g., "pca", "ica").
  Default is "pca".

- pc_windows:

  A list of numeric vectors. Each vector defines a chunk/window of
  dimensions to plot together in a single file. Default is
  `list(1:10, 11:20, 21:30)`.

- output_dir:

  Character. Directory path where the generated JPEGs will be saved.
  Default is "DimHeatmap_plots".

- nfeatures:

  Integer. The number of top features (genes) to display per dimension.
  Default is 60.

- cells:

  Integer. The number of cells to plot. If `cells` is a single number,
  it plots the most extreme cells from both ends of the spectrum.
  Default is 500.

- width_in:

  Numeric. The width of the saved JPEG in inches. Default is 15.

- height_in:

  Numeric. The height of the saved JPEG in inches. Default is 45.

- res:

  Numeric. The resolution (dpi) of the saved JPEG. Default is 300.

## Value

Invisibly returns `NULL`. The primary purpose of this function is its
side effect of saving JPEG plots to the specified output directory.

## Examples

``` r
if (FALSE) { # \dontrun{
# Evaluate the first 15 principal components in 3 chunks of 5
generate_dimheatmaps(
  seurat_obj = my_seurat,
  sample_name = "Donor_1",
  reduction = "pca",
  pc_windows = list(1:5, 6:10, 11:15),
  output_dir = "QC_Plots/PCA_Heatmaps",
  nfeatures = 30,
  cells = 500,
  width_in = 10,
  height_in = 20
)
} # }
```
