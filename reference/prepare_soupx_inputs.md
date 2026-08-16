# Run End-to-End Preparation Pipeline for SoupX Decontamination

Acts as a master wrapper to execute the entire data loading and
preparation phase for ambient RNA removal. It automatically extracts
sample names from the directory structure, loads both raw and filtered
Cell Ranger matrices, initializes the `SoupChannel` objects, and runs a
rapid Seurat clustering pipeline. The output is a comprehensive bundled
list containing all the necessary objects to feed directly into the
final contamination estimation and correction steps.

## Usage

``` r
prepare_soupx_inputs(filtered_dir, raw_dir, manual_contam = NULL)
```

## Arguments

- filtered_dir:

  Character. The file path to the parent directory containing the Cell
  Ranger `"*_filtered_feature_bc_matrix"` folders.

- raw_dir:

  Character. The file path to the parent directory containing the Cell
  Ranger `"*_raw_feature_bc_matrix"` folders.

- manual_contam:

  Character vector. An optional list of specific sample names that will
  later require a manual contamination override (e.g., adding an
  artificial penalty to the estimated rho). Default is `NULL` (which
  evaluates to an empty list).

## Value

A structured, named list containing five top-level elements:

- `sample_names`: A character vector of the extracted sample
  identifiers.

- `counts_lists`: A list containing both `raw` and `filtered` sparse
  count matrices.

- `soup_list`: A list of initialized `SoupChannel` objects.

- `seurat_for_clustering`: A list of processed `Seurat` objects with
  high-resolution clusters stored in their active identities.

- `manual_contam`: The forwarded list of samples designated for manual
  overrides.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run the entire preparation pipeline in one line
soupx_prep_data <- prepare_soupx_inputs(
  filtered_dir = "data/cellranger/filtered",
  raw_dir = "data/cellranger/raw",
  manual_contam = c("Patient_3") # Flag this patient for manual adjustment later
)

# Feed the bundled output directly into your estimation function
clean_data <- estimate_contamination(
  soup_list = soupx_prep_data$soup_list,
  seurat_list = soupx_prep_data$seurat_for_clustering,
  cts_filtered_list = soupx_prep_data$counts_lists$filtered,
  sample_names = soupx_prep_data$sample_names,
  manual_contam = soupx_prep_data$manual_contam,
  mulFac = 5,
  manual = TRUE
)
} # }
```
