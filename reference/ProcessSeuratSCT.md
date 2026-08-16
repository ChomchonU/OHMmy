# Execute Comprehensive SCTransform and Integration Pipeline (Seurat v5)

A unified wrapper for preprocessing and integrating multiple batches or
samples using Seurat v5's `IntegrateLayers` framework. The pipeline
automatically splits the RNA assay by batch, normalizes via
`SCTransform`, and surgically removes TCR/BCR genes from the variable
feature list to prevent clustering based on clonotype. It then computes
PCA and executes the chosen integration algorithm (Harmony, RPCA, CCA,
or FastMNN). Crucially, it dynamically scales integration `k` parameters
(e.g., `k.weight`) downwards to accommodate the smallest batch size,
preventing common integration failures. It also features an interactive
mode to manually select optimal PC dimensions via an elbow plot.

## Usage

``` r
ProcessSeuratSCT(
  seurat_obj,
  batch_col = "batch",
  vars_to_regress = "pct_counts_mt",
  tcr_bcr_patterns = "^TR[ABDG]|^IG[HKL]",
  reduction_name = "pca.SCT",
  integration_method = "HarmonyIntegration",
  integration_reduction = "integrated.har.SCT",
  dims = 1:50,
  interactive_mode = FALSE,
  elbow_plot_dir = NULL,
  k.weight = NULL,
  k.anchor = NULL,
  k.filter = NULL,
  k.score = NULL,
  clustering_resolution = 0.8,
  verbose = TRUE,
  sample_name = "seurat"
)
```

## Arguments

- seurat_obj:

  A Seurat object containing raw count data in the "RNA" assay.

- batch_col:

  Character. The name of the metadata column defining the biological
  batches or samples to split and integrate across. Default is "batch".

- vars_to_regress:

  Character vector. Variables to regress out during `SCTransform` (e.g.,
  cell cycle scores or mitochondrial percentage). Default is
  "pct_counts_mt".

- tcr_bcr_patterns:

  Character. A regular expression matching TCR and BCR gene segments
  (e.g., TRAV, TRBV, IGHV) to exclude them from the variable features
  list. Default is `"^TR[ABDG]|^IG[HKL]"`.

- reduction_name:

  Character. The name to assign to the pre-integration PCA reduction.
  Default is "pca.SCT".

- integration_method:

  Character. The integration algorithm to use in `IntegrateLayers`.
  Options: "HarmonyIntegration", "RPCAIntegration", "CCAIntegration", or
  "FastMNNIntegration". Default is "HarmonyIntegration".

- integration_reduction:

  Character. The name to assign to the final integrated dimensional
  reduction. Default is "integrated.har.SCT".

- dims:

  Numeric vector. The dimensions (PCs) to use for the integration step.
  Default is `1:50`.

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

  Numeric. Included for pipeline compatibility; sets the target
  resolution (though actual clustering is typically handled downstream).
  Default is 0.8.

- verbose:

  Logical. If `TRUE`, outputs progress messages and Seurat logs to the
  console. Default is `TRUE`.

- sample_name:

  Character. A prefix used for saving the elbow plot file. Default is
  "seurat".

## Value

An integrated `Seurat` object with the `DefaultAssay` set to "SCT". The
split RNA and SCT layers are automatically re-joined at the end of the
pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard Harmony integration across patients
integrated_seurat <- ProcessSeuratSCT(
  seurat_obj = raw_seurat,
  batch_col = "Patient_ID",
  vars_to_regress = c("pct_counts_mt", "S.Score", "G2M.Score"),
  integration_method = "HarmonyIntegration",
  dims = 1:30,
  elbow_plot_dir = "QC_Plots/Integrations/"
)

# The returned object is ready for downstream UMAP and Clustering:
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:30, reduction = "integrated.har.SCT")
} # }
```
