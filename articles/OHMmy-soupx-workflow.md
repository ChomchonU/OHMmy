# Ambient RNA Decontamination with SoupX

## Introduction

In droplet-based single-cell RNA sequencing, cell lysis before
partitioning can lead to a “soup” of ambient mRNA. This cell-free RNA
gets incorporated into the droplets and sequenced alongside the actual
cellular RNA. This is especially problematic in immune profiling (e.g.,
T cell and NK cell subsets), where highly expressed transcripts (like
viral genes or hemoglobin) can bleed into unrelated cell clusters.

The `OHMmy` package provides a streamlined, iterative workflow built on
top of `SoupX` to estimate and systematically remove this contamination.

## Part 1: The Standard Workflow

If you want to run a standard decontamination using a default
contamination multiplier (`multiFac`), you can use our one-step wrapper.

``` r


library(Seurat)
library(OHMmy)

# 1. Process SoupX across all samples in your directories
soup_results <- process_soupx_samples(
  filtered_dir = "data/filtered/", 
  raw_dir = "data/raw/", 
  multiFac = 1
)

# 2. Load your raw (pre-SoupX) Seurat object
seurat_raw <- readRDS("data/seurat_raw.rds")

# 3. Transfer the cleaned counts and metadata back to your Seurat object
seurat_clean <- addSoupXMetaToSeurat(
  original_seurat = seurat_raw,
  soupx_seurat_list = soup_results$final_seurat,
  save_path = "data/seurat_soupX_cleaned.rds"
)

# Inspect the new metadata
head(seurat_clean@meta.data)
```

## Part 2: Advanced Iterative Workflow (Testing Parameters)

Often, a default contamination estimate is not aggressive enough for
highly contaminated datasets. `OHMmy` allows you to efficiently loop
through different `multiFac` (multiplier) values to find the optimal
decontamination threshold.

To save computing time, we first prepare the inputs *once*, and then
iterate the SoupX algorithm over those pre-computed inputs.

``` r


# 1. Prepare inputs once (loads matrices and clusters)
prep_data <- prepare_soupx_inputs(
  filtered_dir = "data/filtered/", 
  raw_dir = "data/raw/"
)

# 2. Define the multiplier values you want to test
multiFac_values <- c(5, 10, 15)

# 3. Load the raw Seurat object
seurat_raw <- readRDS("data/seurat_raw.rds")

# 4. Loop through parameters and save a Seurat object for each
for (mf in multiFac_values) {
  message("\nProcessing multiFac = ", mf)
  
  # Run SoupX pipeline with current multiplier
  soup_iter <- run_soupx_post_clustering(
    prep = prep_data, 
    multiFac = mf, 
    methods = "subtraction"
  )
  
  # Define output filename dynamically
  output_name <- paste0("data/seurat_soupX_mf", mf, ".rds")
  
  # Update and save the Seurat object
  seurat_updated <- addSoupXMetaToSeurat(
    original_seurat = seurat_raw,
    soupx_seurat_list = soup_iter$final_seurat,
    save_path = output_name
  )
  
  message("Successfully saved: ", output_name)
}
```

## Part 3: Handling Manual Sample Overrides

In large cohorts, specific samples (e.g., highly necrotic tissue samples
or specific donor IDs) may require forced manual contamination
thresholds because the automated algorithm fails to accurately estimate
their background soup.

You can pass a list of these problematic sample IDs directly to the
pipeline.

``` r


# Identify problematic samples requiring manual intervention
problem_samples <- c("Donor_075", "Donor_089", "Donor_160")

# Run the pipeline applying manual thresholds to these specific samples
soup_manual <- run_soupx_post_clustering(
  prep = prep_data, 
  multiFac = 10, 
  manual_contam = problem_samples
)

# Save the final manually-adjusted object
seurat_final <- addSoupXMetaToSeurat(
  original_seurat = seurat_raw,
  soupx_seurat_list = soup_manual$final_seurat,
  save_path = "data/seurat_soupX_manual_adjusted.rds"
)
```

## Part 4: Evaluating the Results

Once your loops have finished running, you can easily load the different
output files to compare the metadata and ensure the contamination
metrics were applied correctly across your iterations.

``` r


# Load a few of the generated objects
seurat_mf5  <- readRDS("data/seurat_soupX_mf5.rds")
seurat_mf10 <- readRDS("data/seurat_soupX_mf10.rds")

# Compare SoupX metadata columns (e.g., nUMI_soupX, contamination_fraction)
head(seurat_mf5@meta.data)
head(seurat_mf10@meta.data)
```
