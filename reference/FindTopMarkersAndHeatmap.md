# Find and Visualize Top Cluster Markers

Computes differential gene expression for a Seurat object-either across
all clusters (`FindAllMarkers`) or as a specific pairwise comparison. It
calculates a custom biological relevance score (combining log2 fold
change, percentage difference, and adjusted p-value) to strictly filter
and rank the top `n` markers per cluster. The function automatically
generates a Seurat `DoHeatmap` and exports multiple CSV tables
containing both the filtered and unfiltered marker lists (ready for
GSEA/ORA).

## Usage

``` r
FindTopMarkersAndHeatmap(
  seurat_obj,
  sample_name = "Sample",
  use_sct = FALSE,
  marker_diff_thresh = 0.1,
  marker_pval_adj = 0.05,
  marker_avg_log2FC_thresh = 0.5,
  top_n = 20,
  output_dir_base = "Plots_heatmap",
  plot_format = "jpg",
  width = 10,
  height = 20,
  dpi = 300,
  compare = NULL,
  add_timestamp = TRUE,
  onlyPos = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data with active identities
  (`Idents`) set.

- sample_name:

  Character. The name of the biological sample, used for plot titles,
  filenames, and sub-directory routing. Default is "Sample".

- use_sct:

  Logical. If `TRUE`, sets the default assay to `"SCT"`, runs
  [`PrepSCTFindMarkers`](https://satijalab.org/seurat/reference/PrepSCTFindMarkers.html)
  to ensure model comparability across samples/batches, and performs
  marker identification on the SCT assay. Defaults to `FALSE`.

- marker_diff_thresh:

  Numeric. The minimum required absolute difference in the percentage of
  expressing cells between groups (`abs(pct.1 - pct.2)`). Default is
  0.1.

- marker_pval_adj:

  Numeric. The maximum adjusted p-value allowed for a gene to be
  considered a significant marker. Default is 0.05.

- marker_avg_log2FC_thresh:

  Numeric. The minimum absolute log2 fold change required. Default is
  0.5.

- top_n:

  Integer. The number of top-scoring marker genes to select per cluster
  for plotting on the heatmap. Default is 20.

- output_dir_base:

  Character. The base directory path where outputs (plots and CSVs) will
  be saved. Default is "Plots_heatmap".

- plot_format:

  Character. The file format for the saved heatmap ("jpg", "png", or
  "pdf"). Default is "jpg".

- width:

  Numeric. The width of the saved heatmap. Default is 10.

- height:

  Numeric. The height of the saved heatmap. Default is 20.

- dpi:

  Numeric. The resolution of the saved heatmap. Default is 300.

- compare:

  Character vector of length 2. If provided, performs a pairwise
  `FindMarkers` test between these two specific identities (e.g.,
  `c("Effector_T", "Naive_T")`). If `NULL`, runs `FindAllMarkers` across
  all clusters. Default is NULL.

- add_timestamp:

  Logical. Whether to append the current date and time to the saved
  filenames. Default is TRUE.

- onlyPos:

  Logical. If TRUE, only identifies positive markers (upregulated
  genes). Passed to the `only.pos` argument in Seurat. Default is TRUE.

## Value

A list containing four elements:

- `top_markers`: A `tibble` of the highly filtered, top-scoring marker
  genes used for the heatmap.

- `heatmap`: The `ggplot` object of the generated `DoHeatmap`.

- `markers`: A data frame of the full, unfiltered differential
  expression results.

- `output_dir`: A character string of the specific path where the files
  were saved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Find top 10 positive markers for ALL clusters in the object
all_markers_res <- FindTopMarkersAndHeatmap(
  seurat_obj = my_seurat,
  sample_name = "PBMC_Global",
  top_n = 10,
  onlyPos = TRUE
)

# Example 2: Perform a specific pairwise comparison (e.g., CD8+ Effector vs Naive)
# Capturing both up and downregulated genes
pairwise_res <- FindTopMarkersAndHeatmap(
  seurat_obj = my_seurat,
  sample_name = "CD8_Subset",
  compare = c("CD8_Effector", "CD8_Naive"),
  onlyPos = FALSE,
  marker_avg_log2FC_thresh = 0.25,
  save_format = "png"
)
} # }
```
