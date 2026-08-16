# Generate a Discrete Color Palette for Clustering

Creates a maximally distinct color palette tailored for high-dimensional
single-cell clustering visualizations. Users can choose from curated
categorical palettes (Kelly, Alphabet, Polychrome) or an "auto" mode
that combines multiple `RColorBrewer` sets. If the requested number of
colors (`n`) exceeds the available base colors in a given palette, the
function automatically interpolates using `colorRampPalette` to generate
the required amount.

## Usage

``` r
make_cluster_palette(n, palette = "auto")
```

## Arguments

- n:

  Integer. The number of distinct colors required (e.g., the number of
  Seurat clusters).

- palette:

  Character. The specific discrete palette to use. Options are `"auto"`
  (default), `"kelly"`, `"alphabet"`, or `"polychrome"`.

## Value

A character vector of length `n` containing hexadecimal color codes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate a 12-color palette using the default hybrid RColorBrewer sets
auto_colors <- make_cluster_palette(n = 12, palette = "auto")

# Generate a 20-color palette using the Kelly distinct color list
kelly_colors <- make_cluster_palette(n = 20, palette = "kelly")

# Pass the colors directly into a Seurat plotting function
DimPlot(my_seurat, cols = make_cluster_palette(15, "alphabet"))

# Requesting 50 colors will automatically trigger interpolation
massive_palette <- make_cluster_palette(n = 50, palette = "polychrome")
} # }
```
