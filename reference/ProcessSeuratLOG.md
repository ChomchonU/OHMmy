# Execute Comprehensive LogNormalization and Integration Pipeline (Seurat v5)

A unified wrapper for preprocessing and integrating multiple batches or
samples using Seurat v5's `IntegrateLayers` framework with standard
LogNormalization. The pipeline automatically splits the RNA assay by
batch, normalizes, and removes TCR/BCR genes from the variable feature
list to prevent clustering driven by clonotype. It then scales the data,
computes PCA, and executes the chosen integration algorithm (Harmony,
RPCA, CCA, or FastMNN). Crucially, it dynamically scales integration `k`
parameters (e.g., `k.weight`) downwards to accommodate the smallest
batch size, preventing common integration failures.

## Usage

``` r
ProcessSeuratLOG(
  seurat_obj,
  batch_col = "batch",
  vars_to_regress = NULL,
  tcr_bcr_patterns = "^TR[ABDG]|^IG[HKL]",
  reduction_name = "pca.SCT",
  integration_method = "HarmonyIntegration",
  integration_reduction = "integrated.har.SCT",
  dims = 1:30,
  interactive_mode = FALSE,
  elbow_plot_dir = NULL,
  k.weight = NULL,
  k.anchor = NULL,
  k.filter = NULL,
  k.score = NULL,
  clustering_resolution = 1,
  verbose = TRUE,
  sample_name = "seurat",
  boost_genes = c("CD3D", "CD3E", "CD3G", "TYROBP", "FCGR3A", "NCAM1"),
  boost_multiplier = 1
)
```

## Arguments

- seurat_obj:

  A Seurat object containing raw count data in the "RNA" assay.

- batch_col:

  Character. The name of the metadata column defining the biological
  batches or samples to split and integrate across. Default is "batch".

- vars_to_regress:

  Character vector. Variables to regress out during `ScaleData` (e.g.,
  cell cycle scores or mitochondrial percentage). Default is `NULL`.

- tcr_bcr_patterns:

  Character. A regular expression matching TCR and BCR gene segments to
  exclude from the variable features list. Default is
  `"^TR[ABDG]|^IG[HKL]"`. Note that this default also matches
  non-receptor genes such as *TRAF1* or *IGHMBP2*;
  `"^TR[ABDG]V|^IG[HKL]V"` restricts the filter to variable segments.

- reduction_name:

  Character. The name to assign to the pre-integration PCA reduction.
  Default is "pca.SCT" (kept for backward compatibility; a name such as
  "pca.log" is clearer for this workflow).

- integration_method:

  Character. The integration algorithm to use in `IntegrateLayers`.
  Options: "HarmonyIntegration", "RPCAIntegration", "CCAIntegration", or
  "FastMNNIntegration" (the latter requires the SeuratWrappers package
  to be attached). Default is "HarmonyIntegration".

- integration_reduction:

  Character. The name to assign to the final integrated dimensional
  reduction. Default is "integrated.har.SCT" (kept for backward
  compatibility; a name such as "integrated.har.log" is clearer for this
  workflow).

- dims:

  Numeric vector. The dimensions (PCs) to use for the integration step
  (RPCA/CCA). Default is `1:30`.

- interactive_mode:

  Logical. If `TRUE` and running in an interactive session, pauses to
  prompt the user for the optimal number of PCs and k-parameters after
  computing the initial PCA. Default is `FALSE`.

- elbow_plot_dir:

  Character. An optional directory path to save a JPG of the PCA elbow
  plot. Default is `NULL` (does not save).

- k.weight:

  Integer. The number of neighbors to consider when weighting anchors.
  If `NULL`, defaults to 100 or the size of the smallest batch minus 1.

- k.anchor:

  Integer. The number of neighbors to use for picking anchors
  (RPCA/CCA). If `NULL`, defaults to 5.

- k.filter:

  Integer. The number of neighbors to use for filtering anchors
  (RPCA/CCA). If `NULL`, defaults to 200.

- k.score:

  Integer. The number of neighbors to use for scoring anchors
  (RPCA/CCA). If `NULL`, defaults to 30.

- clustering_resolution:

  Numeric. Retained for pipeline compatibility; clustering itself is
  performed downstream (e.g., by
  [`ClusterAndUMAP()`](https://chomchonu.github.io/OHMmy/reference/ClusterAndUMAP.md)).
  Default is 1.

- verbose:

  Logical. If `TRUE`, outputs progress messages and Seurat logs to the
  console. Default is `TRUE`.

- sample_name:

  Character. A prefix used for saving the elbow plot file. Default is
  "seurat".

- boost_genes:

  A character vector of gene names whose variance should be artificially
  increased prior to PCA. This forces the dimensionality reduction to
  prioritize these specific lineage markers, which is highly useful for
  cleanly separating biologically distinct but transcriptomically
  similar populations (e.g., NK cells vs. CD8+ T cells). Only used when
  `boost_multiplier > 1`; set to `NULL` to disable. Default is
  `c("CD3D", "CD3E", "CD3G", "TYROBP", "FCGR3A", "NCAM1")`.

- boost_multiplier:

  A numeric value indicating the weight factor applied to the
  `boost_genes`. The scaled expression data for these genes will be
  multiplied by this number. Default is `1` (no boosting).

## Value

An integrated `Seurat` object with the `DefaultAssay` set to "RNA". The
split RNA layers are automatically re-joined at the end of the pipeline.

## See also

[`ProcessSeuratSCT()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratSCT.md)
for the SCTransform equivalent and
[`ClusterAndUMAP()`](https://chomchonu.github.io/OHMmy/reference/ClusterAndUMAP.md)
for the next step.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard Harmony integration across batches using LogNormalization
integrated_seurat <- ProcessSeuratLOG(
  seurat_obj = raw_seurat,
  batch_col = "Batch_ID",
  vars_to_regress = "pct_counts_mt",
  reduction_name = "pca",
  integration_method = "HarmonyIntegration",
  integration_reduction = "integrated.har",
  dims = 1:30,
  interactive_mode = TRUE,
  elbow_plot_dir = "QC_Plots/Integrations/"
)

# The returned object is ready for downstream UMAP and Clustering:
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30, reduction = "integrated.har")
} # }
```
