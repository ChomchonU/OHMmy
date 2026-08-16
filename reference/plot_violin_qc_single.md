# Generate Violin Plots for QC Metrics and Gene Expression

This function generates highly customized violin plots for both quality
control (QC) metadata and specific gene expression markers across cell
clusters. It enhances standard Seurat `VlnPlot` outputs by adding faded,
size-adjusted jittered points to better visualize true single-cell
distributions. Plots are combined via `patchwork` and safely saved to
disk with an automatic PDF fallback.

## Usage

``` r
plot_violin_qc_single(
  seurat_obj,
  sample_name = "Sample",
  qc_meta_features = c("nCount_RNA", "nFeature_RNA", "pct_counts_mt"),
  gene_features = NULL,
  res_col = "seurat_clusters",
  assay = "RNA",
  output_dir = "Plots_violin_qc",
  save_format = "jpg",
  width = 20,
  dpi = 300,
  add_timestamp = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- sample_name:

  Character. The name of the biological sample, used for plot titles and
  file naming. Default is "Sample".

- qc_meta_features:

  Character vector. Metadata columns to plot for quality control (e.g.,
  "nCount_RNA", "nFeature_RNA"). Default is c("nCount_RNA",
  "nFeature_RNA", "pct_counts_mt").

- gene_features:

  Character vector. Specific genes to visualize across clusters. Default
  is NULL.

- res_col:

  Character. The metadata column containing cluster assignments or cell
  identities. Default is "seurat_clusters".

- assay:

  Character. The assay to pull gene expression data from. Default is
  "RNA".

- output_dir:

  Character. Directory path where the generated plots will be saved.
  Default is "Plots_violin_qc".

- save_format:

  Character. The file format for the saved plots ("jpg", "png", or
  "pdf"). Default is "jpg".

- width:

  Numeric. The base width of the saved plot. Default is 20.

- dpi:

  Numeric. The resolution of the saved plots. Default is 300.

- add_timestamp:

  Logical. Whether to append the current date and time to the saved
  filenames. Default is TRUE.

## Value

Invisibly returns `NULL`. The primary purpose of this function is its
side effect of saving combined plots to the specified output directory.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot basic QC metrics across clusters
plot_violin_qc_single(
  seurat_obj = my_seurat,
  sample_name = "Donor_1",
  res_col = "seurat_clusters"
)

# Plot both QC metrics AND specific T/NK cell markers
plot_violin_qc_single(
  seurat_obj = my_seurat,
  sample_name = "Donor_1",
  qc_meta_features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  gene_features = c("CD3E", "CD8A", "NKG7", "GNLY"),
  res_col = "CellType",
  save_format = "png"
)
} # }
```
