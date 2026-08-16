# Generate a 4-Panel Blended Feature Plot (Seurat v5 Compatible)

Replicates and enhances the behavior of Seurat's original
`FeaturePlot(blend = TRUE, combine = TRUE)` functionality, with built-in
compatibility for Seurat v5's new `layer` architecture (while
maintaining a v4 `slot` fallback). It computes normalized, thresholded
expression and generates a 1x4 `patchwork` grid containing: individual
expression of gene 1, individual expression of gene 2, their blended
co-expression on the embedding, and a customized 2D color threshold key.

## Usage

``` r
.blend_feature_plot_v5(
  seurat_obj,
  gene1,
  gene2,
  reduction,
  cols = c("lightgrey", "#00ff00", "#ff0000"),
  blend_threshold = 0.1,
  min_cutoff = "q10",
  gamma = 1,
  pt_size = 1
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- gene1:

  Character. The name of the first feature/gene to plot (mapped to the
  x-axis of the color key).

- gene2:

  Character. The name of the second feature/gene to plot (mapped to the
  y-axis of the color key).

- reduction:

  Character. The dimensionality reduction to visualize (e.g., "umap",
  "pca").

- cols:

  Character vector of length 3. Specifies the colors for the background
  (double-negative), gene 1, and gene 2. Default is c("lightgrey",
  "#00ff00", "#ff0000").

- blend_threshold:

  Numeric. A scaling threshold (0 to 1) below which normalized
  expression values are visually suppressed to 0. Default is 0.1.

- min_cutoff:

  Character or Numeric. The lowest expression cutoff (e.g., "q10" for
  the 10th quantile). Values below this are set to zero before scaling.
  Default is "q10".

- gamma:

  Numeric. The exponent used to adjust the brightness and saturation of
  the blended colors. Default is 1.

- pt_size:

  Numeric. The point size for cells in the scatter plots. Default is 1.

## Value

A `patchwork` object containing four aligned `ggplot` panels.

## Examples

``` r
if (FALSE) { # \dontrun{
# Visualize co-expression of cytotoxic markers
blend_plot <- .blend_feature_plot_v5(
  seurat_obj = my_seurat,
  gene1 = "CD8A",
  gene2 = "GZMB",
  reduction = "umap",
  cols = c("lightgrey", "#ff0000", "#0000ff"), # Custom red and blue blend
  blend_threshold = 0.15,
  pt_size = 0.5
)

# Display the 4-panel plot
print(blend_plot)
} # }
```
