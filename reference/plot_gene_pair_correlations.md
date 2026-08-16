# Calculate and Plot Gene-Pair Correlations Across Clusters

Computes the expression correlation (e.g., Pearson or Spearman) between
specified pairs of genes within each cell cluster. To mitigate
zero-inflation artifacts common in scRNA-seq data, it filters cells
based on a defined quantile expression threshold before calculating the
correlation. The results are visualized as faceted bar plots and
automatically saved to disk.

## Usage

``` r
plot_gene_pair_correlations(
  seurat_obj,
  cluster_col,
  gene_pairs,
  sample_name = "Sample",
  cor_method = "pearson",
  quantile_thresh = 0.01,
  fill_palette = "auto",
  output_dir = "Plots_gene_pair_cor",
  save_format = "jpg",
  width = 14,
  height = 8,
  dpi = 300,
  add_timestamp = TRUE,
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell expression data.

- cluster_col:

  Character. The name of the metadata column defining cell clusters to
  group the correlations by.

- gene_pairs:

  A list of character vectors, where each vector contains exactly two
  gene names to correlate (e.g., `list(c("GeneA", "GeneB"))`).

- sample_name:

  Character. The name of the biological sample, used for plot titles and
  file naming. Default is "Sample".

- cor_method:

  Character. The correlation method to use, passed to
  [`cor()`](https://rdrr.io/r/stats/cor.html). Options include
  "pearson", "spearman", or "kendall". Default is "pearson".

- quantile_thresh:

  Numeric. The expression quantile threshold (0 to 1) used to filter out
  lowly expressing cells before correlation calculations. Default is
  0.01.

- fill_palette:

  Character. The discrete color palette to pass to
  [`make_cluster_palette()`](https://chomchonu.github.io/OHMmy/reference/make_cluster_palette.md).
  Default is "auto".

- output_dir:

  Character. Directory path where the generated plots will be saved.
  Default is "Plots_gene_pair_cor".

- save_format:

  Character. The file format for the saved plots ("jpg", "png", or
  "pdf"). Default is "jpg".

- width:

  Numeric. The width of the saved plot. Default is 14.

- height:

  Numeric. The height of the saved plot. Default is 8.

- dpi:

  Numeric. The resolution of the saved plots. Default is 300.

- add_timestamp:

  Logical. Whether to append the current date and time to the saved
  filenames. Default is TRUE.

- verbose:

  Logical. Whether to print progress messages and display a progress
  bar. Default is TRUE.

## Value

A list containing three elements:

- `plot`: The generated `ggplot` object containing the faceted bar
  charts.

- `data`: A `tibble` (data frame) containing the calculated correlation
  values and cell counts per cluster for each valid gene pair.

- `output_dir`: A character string of the path where the plot was saved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Define pairs of genes expected to co-express in specific subsets
t_nk_pairs <- list(
  c("CD8A", "GZMB"),   # Cytotoxic T cell markers
  c("NKG7", "GNLY"),   # NK cell effectors
  c("CD4", "IL7R")     # Naive/Memory CD4 T cell markers
)

# Run the correlation analysis across predefined clusters
cor_results <- plot_gene_pair_correlations(
  seurat_obj = my_seurat,
  cluster_col = "CellType_Subset",
  gene_pairs = t_nk_pairs,
  sample_name = "Donor_1",
  cor_method = "spearman",
  quantile_thresh = 0.05,
  save_format = "png"
)

# View the raw correlation dataframe directly
head(cor_results$data)
} # }
```
