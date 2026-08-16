# Update Seurat Object with SoupX-Corrected Counts and Metadata

Merges a list of SoupX-processed Seurat objects and seamlessly maps
their background-corrected count matrices and updated QC metrics back
onto a pre-existing Seurat object. It handles the tedious formatting
steps, such as barcode stripping, standardizing mitochondrial column
names, and intersecting cells/genes to prevent dimension mismatches.
Finally, it replaces the raw counts layer in the original object with
the decontaminated data.

## Usage

``` r
addSoupXMetaToSeurat(
  original_seurat,
  soupx_seurat_list,
  mito_colname = "mitoPercent",
  counts_assay = "RNA",
  rename_barcodes = TRUE,
  save_path = NULL
)
```

## Arguments

- original_seurat:

  A pre-existing `Seurat` object that you wish to update with the
  cleaned counts.

- soupx_seurat_list:

  A named list of `Seurat` objects containing SoupX-corrected counts,
  typically the output of
  [`create_final_seurat()`](https://chomchonu.github.io/OHMmy/reference/create_final_seurat.md).

- mito_colname:

  Character. The name of the mitochondrial percentage column in the
  SoupX objects. It will be automatically standardized to
  `"pct_counts_mt"` during the transfer. Default is `"mitoPercent"`.

- counts_assay:

  Character. The name of the assay in the original object where the
  counts layer should be updated. Default is `"RNA"`.

- rename_barcodes:

  Logical. If `TRUE`, strips the common `"-1"` suffix from the barcodes
  in the merged SoupX object to ensure they match standard Seurat
  formatting. Default is `TRUE`.

- save_path:

  Character. An optional file path (.rds) to automatically save the
  updated Seurat object to disk. Default is `NULL` (does not save).

## Value

An updated `Seurat` object containing the decontaminated count matrix
and matching QC metadata (`nCount_RNA`, `nFeature_RNA`,
`pct_counts_mt`). The object is strictly subsetted to contain only the
overlapping valid cells and genes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'my_seurat' is an existing object you want to clean,
# and 'pipeline_out' is the list returned by process_soupx_samples()

cleaned_seurat <- addSoupXMetaToSeurat(
  original_seurat = my_seurat,
  soupx_seurat_list = pipeline_out$final_seurat,
  mito_colname = "mitoPercent",
  counts_assay = "RNA",
  rename_barcodes = TRUE,
  save_path = "Results/Seurat_SoupX_Cleaned.rds"
)

# The 'cleaned_seurat' object now contains the background-corrected counts
# in its RNA assay, ready for standard downstream normalization!
} # }
```
