# Initialize Final Seurat Objects from Decontaminated Counts

Takes a list of fully processed `SoupChannel` objects (where the
contamination fraction has already been estimated and set), extracts the
background-corrected count matrices, and initializes fresh `Seurat`
objects. During initialization, it applies standard baseline filtering
(`min.cells = 3`, `min.features = 200`) and automatically calculates the
mitochondrial read percentage (`mitoPercent`) for immediate use in
downstream quality control.

## Usage

``` r
create_final_seurat(soup_list_all, sample_names)
```

## Arguments

- soup_list_all:

  A named list of `SoupChannel` objects that have their contamination
  fraction (`rho`) set. Typically the `soup_list_all` element returned
  by
  [`estimate_contamination()`](https://chomchonu.github.io/OHMmy/reference/estimate_contamination.md).

- sample_names:

  Character vector. A list of sample identifiers used to match the list
  elements and assign the `project` name to each Seurat object.

## Value

A named list of initialized `Seurat` objects containing the
decontaminated count matrices and a newly appended `mitoPercent`
metadata column.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming clean_data is the list returned by estimate_contamination()

final_seurat_list <- create_final_seurat(
  soup_list_all = clean_data$soup_list_all,
  sample_names = my_samples
)

# The objects are now ready for standard Seurat QC filtering!
patient1_seurat <- final_seurat_list[["Patient_1"]]
head(patient1_seurat@meta.data)
} # }
```
