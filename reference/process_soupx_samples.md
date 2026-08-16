# Execute the Complete One-Step SoupX Decontamination Pipeline

A unified master wrapper that executes the entire ambient RNA removal
workflow in a single command. It automatically extracts sample names
from the directory structure, loads the raw and filtered Cell Ranger
matrices, initializes `SoupChannel` objects, runs a quick clustering
pipeline to define groups, estimates and removes ambient RNA
contamination, and immediately wraps the cleaned count matrices into a
list of fresh `Seurat` objects ready for downstream quality control.

## Usage

``` r
process_soupx_samples(filtered_dir, raw_dir, multiFac = 0)
```

## Arguments

- filtered_dir:

  Character. The file path to the parent directory containing the Cell
  Ranger `"*_filtered_feature_bc_matrix"` folders.

- raw_dir:

  Character. The file path to the parent directory containing the Cell
  Ranger `"*_raw_feature_bc_matrix"` folders.

- multiFac:

  Numeric. A percentage to artificially add to the auto-estimated
  contamination fraction (`rho`) across all samples (e.g., passing 5
  adds 0.05 to the estimated rho). Default is 0.

## Value

A structured, named list containing three top-level elements:

- `final_seurat`: A list of the final, initialized `Seurat` objects
  containing the background-corrected counts and baseline metadata.

- `soup_objects`: The output list from the contamination estimation
  step, containing the updated `SoupChannel` objects and corrected
  matrices.

- `seurat_clustering`: A list of the preliminary `Seurat` objects used
  to generate the high-resolution clusters for ambient RNA estimation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run the entire start-to-finish pipeline in a single line
final_pipeline_output <- process_soupx_samples(
  filtered_dir = "data/cellranger/filtered",
  raw_dir = "data/cellranger/raw",
  multiFac = 2 # Globally add a 2% contamination penalty to all samples
)

# Extract your final cleaned Seurat objects for downstream integration
my_clean_seurats <- final_pipeline_output$final_seurat
} # }
```
