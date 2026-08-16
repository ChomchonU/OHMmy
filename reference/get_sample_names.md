# Extract Sample Names from Cell Ranger Output Directories

Scans a specified parent directory for 10x Genomics Cell Ranger output
folders and extracts clean sample names. It automatically formats the
names by stripping standard Cell Ranger suffixes (such as
`"_filtered_feature_bc_matrix"` or `"_feature_bc_matrix"`) from the
sub-directory folder names.

## Usage

``` r
get_sample_names(filtered_dir)
```

## Arguments

- filtered_dir:

  Character. The file path to the parent directory containing the
  individual sample output folders.

## Value

A character vector of cleaned sample names, derived from the
sub-directory names.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have a directory "data/cellranger_outputs" containing folders like:
# "Patient1_filtered_feature_bc_matrix" and "Patient2_filtered_feature_bc_matrix"

my_samples <- get_sample_names(filtered_dir = "data/cellranger_outputs")

print(my_samples)
# [1] "Patient1" "Patient2"
} # }
```
