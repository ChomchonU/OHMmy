# Calculate Consensus Labels for CD8 vs NK Cells

Evaluates competing CD8 scoring variants against a baseline NK score to
assign discrete cell labels (CD8 vs NK). It establishes a cutoff
threshold (using fixed, density trough, or Gaussian Mixture Model
methods), calculates consensus between a primary and confirmatory
variant, optionally bootstraps for stability, and generates diagnostic
plots.

## Usage

``` r
score_consensus(
  obj,
  cd8_cols,
  nk_col = "NK_Score",
  cells_use = colnames(obj),
  method = c("fixed", "trough", "gmm"),
  fixed_cut = 0,
  zscore = TRUE,
  primary = NULL,
  confirm = NULL,
  min_features = NULL,
  feature_col = "nFeature_RNA",
  n_boot = 0,
  boot_target = c("cutoff", "proportion"),
  add_meta = TRUE,
  prefix = "CD8vNK",
  seed = 42,
  make_plots = FALSE,
  plot_group_by = c("consensus", "variant_label"),
  plot_n = 1000,
  out_dir = NULL
)
```

## Arguments

- obj:

  A single-cell object (e.g., Seurat) that contains a `@meta.data` slot.

- cd8_cols:

  Character vector; names of the columns in `obj@meta.data` containing
  the CD8 scores to test.

- nk_col:

  Character string; name of the column containing the baseline NK score.
  Default: `"NK_Score"`.

- cells_use:

  Character vector; specific cell barcodes to include. Defaults to all
  cells in `obj`.

- method:

  Character; the method to determine the CD8 vs NK decision boundary.
  One of `"fixed"`, `"trough"` (density valley), or `"gmm"` (Gaussian
  Mixture Model).

- fixed_cut:

  Numeric; the threshold to use if `method = "fixed"`. Default: `0`.

- zscore:

  Logical; whether to standardize (Z-score) the CD8 and NK columns
  independently before subtracting them. Default: `TRUE`.

- primary:

  Character string; the CD8 variant to act as the primary label for the
  consensus calculation. Defaults to `cd8_cols[1]`.

- confirm:

  Character string; the CD8 variant to act as the confirming label for
  the consensus calculation. Defaults to `cd8_cols[2]`.

- min_features:

  Numeric; optional threshold to filter out low-quality cells prior to
  scoring.

- feature_col:

  Character string; the metadata column used for `min_features`
  filtering. Default: `"nFeature_RNA"`.

- n_boot:

  Numeric; number of bootstrap iterations to run to test
  cutoff/proportion stability. Set to `0` to skip.

- boot_target:

  Character; what metric to track during bootstrapping. Either
  `"cutoff"` (the threshold value) or `"proportion"` (fraction of cells
  labeled CD8).

- add_meta:

  Logical; whether to append the computed differences, scaled scores,
  and labels back into `obj@meta.data`.

- prefix:

  Character string; the prefix applied to the new columns added to
  `meta.data`. Default: `"CD8vNK"`.

- seed:

  Numeric; the random seed used for bootstrapping and plot downsampling.
  Default: `42`.

- make_plots:

  Logical; whether to generate diagnostic plots (Violin distributions,
  Agreement Heatmaps, etc.).

- plot_group_by:

  Character; how to group/color cells in the violin plots. Either
  `"consensus"` or by individual `"variant_label"`.

- plot_n:

  Numeric; maximum number of cells per class to sample for violin plots
  to ensure fast rendering. Default: `1000`.

- out_dir:

  Character string; optional directory path to save the generated plots
  (e.g., `"./plots/"`). If `NULL`, plots are only returned in the output
  list.

## Value

A list containing:

- `obj`: The updated single-cell object with new metadata columns (if
  `add_meta = TRUE`).

- `per`: A list of data frames with detailed, per-cell scoring and
  labeling for each variant.

- `labels`: A matrix of the assigned labels (CD8/NK) for all cells
  across all variants.

- `consensus`: A named character vector of the final consensus labels
  per cell.

- `agreement`: A pairwise matrix showing the proportion of label
  agreement between variants.

- `n_disagree`: Integer; the number of cells that had conflicting labels
  among any variants.

- `cutoffs`: A data frame summarizing the calculated thresholds and
  statistics for each variant.

- `bootstrap`: A list containing bootstrap draws and summary statistics
  (if `n_boot > 0`).

- `plots`: A list of `ggplot2` objects generated during the run (if
  `make_plots = TRUE`).

- `cells_used`: Character vector of valid cell IDs retained after
  filtering.

- `dropped`: Character vector of cell IDs dropped during filtering.

- `params`: A list recording the input parameters used for this run.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with a fixed cutoff of 0 (on Z-scored data)
res <- score_consensus(
  obj = seurat_obj,
  cd8_cols = c("CD8_Score_TCR", "CD8_Score_GEX"),
  nk_col = "NK_Score",
  method = "fixed",
  fixed_cut = 0,
  make_plots = TRUE,
  out_dir = "./scoring_plots"
)

# Extract the updated object for downstream Seurat workflows
seurat_obj <- res$obj
} # }
```
