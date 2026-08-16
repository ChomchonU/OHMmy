# Automate Multi-Resolution Clustering, Clustree Visualization, and UMAP

A comprehensive wrapper for Seurat's graph-based clustering and UMAP
workflow. It intelligently checks for existing Nearest Neighbor graphs
and cluster resolutions to avoid redundant computations unless
explicitly forced. The function sweeps through a provided sequence of
clustering resolutions, generates and saves a `clustree` plot to help
visualize cluster stability, and computes a final UMAP embedding at a
specified resolution.

## Usage

``` r
ClusterAndUMAP(
  seurat_obj,
  sample_name = "Sample",
  dims = 1:50,
  k.param = 20,
  reduction = "integrated.har",
  umap_name = "umap.har",
  cluster_prefix = "RNA_snn_res.",
  graph_name = "RNA_snn",
  cluster_resolutions = seq(0.1, 2, by = 0.1),
  final_resolution = 0.7,
  save_path = NULL,
  force_neighbors = FALSE,
  force_clustering = FALSE,
  plot_dir = "Plots_clustree",
  return.model = TRUE,
  plot_format = "jpg",
  width = 10,
  height = 15,
  dpi = 300
)
```

## Arguments

- seurat_obj:

  A Seurat object containing the dimensional reduction specified in
  `reduction`.

- sample_name:

  Character. The identifier for the sample or dataset, used for console
  logging, plot titles, and file naming. Default is "Sample".

- dims:

  Numeric vector. The dimensions of the reduction to use as input for
  constructing the neighbor graph and UMAP (e.g., `1:50`). Default is
  `1:50`.

- k.param:

  Integer. The number of nearest neighbors to compute during
  `FindNeighbors`. Default is 20.

- reduction:

  Character. The name of the dimensional reduction to use (e.g., "pca",
  "integrated.har", "harmony"). Default is "integrated.har".

- umap_name:

  Character. The name to assign to the generated UMAP reduction. Default
  is "umap.har".

- cluster_prefix:

  Character. The prefix to use for naming the cluster metadata columns.
  Default is "RNA_snn_res.".

- graph_name:

  Character. The name to assign to the generated Shared Nearest Neighbor
  (SNN) graph. Default is "RNA_snn".

- cluster_resolutions:

  Numeric vector. A sequence of resolutions to sweep through for
  clustering. Default is `seq(0.1, 2, by = 0.1)`.

- final_resolution:

  Numeric. The specific resolution to use for the final clustering step
  immediately prior to calculating the UMAP. Default is 0.7.

- save_path:

  Character. An optional file path (.rds) to save the updated Seurat
  object. Default is `NULL` (does not save).

- force_neighbors:

  Logical. If `TRUE`, forces recalculation of the neighbor graph even if
  `graph_name` already exists. Default is `FALSE`.

- force_clustering:

  Logical. If `TRUE`, forces recalculation of clusters even if the
  resolution columns already exist in the metadata. Default is `FALSE`.

- plot_dir:

  Character. Directory path where the generated `clustree` plot will be
  saved. Default is "Plots_clustree".

- return.model:

  Logical. If `TRUE`, retains the UMAP model in the Seurat object
  (useful for projecting new data later). Default is `TRUE`.

- plot_format:

  Character. The file format for the saved `clustree` plot (e.g., "jpg",
  "png", "pdf"). Default is "jpg".

- width:

  Numeric. The width of the saved plot in inches. Default is 10.

- height:

  Numeric. The height of the saved plot in inches. Default is 15.

- dpi:

  Numeric. The resolution (dpi) of the saved plot. Default is 300.

## Value

A named list containing three elements:

- `seurat`: The updated `Seurat` object containing the new graph,
  clusters, and UMAP reduction.

- `clustree`: The generated `ggplot` object containing the clustering
  tree.

- `plot_file`: A character string of the exact file path where the plot
  was saved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'my_seurat' has already undergone PCA or Harmony integration

clustering_results <- ClusterAndUMAP(
  seurat_obj = my_seurat,
  sample_name = "PBMC_Donor_1",
  dims = 1:30,
  reduction = "pca",
  umap_name = "umap",
  cluster_resolutions = seq(0.2, 1.2, by = 0.2),
  final_resolution = 0.6,
  plot_dir = "Results/Clustree"
)

# Extract the updated Seurat object for downstream analysis
updated_seurat <- clustering_results$seurat

# View the clustree plot in R
print(clustering_results$clustree)
} # }
```
