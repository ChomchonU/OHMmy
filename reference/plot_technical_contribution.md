# Evaluate Technical Gene Contributions to Principal Components

Extracts PCA loadings from a Seurat object and calculates the weighted
percentage of technical or nuisance genes (e.g., mitochondrial,
ribosomal, immunoglobulins, or specific lncRNAs) driving each principal
component. It evaluates the top positive and negative feature loadings
independently. The output is a split bar chart that visually highlights
which PCs are overwhelmed by technical noise, helping determine which
components to exclude or which variables require regression
(`vars.to.regress`).

## Usage

``` r
plot_technical_contribution(
  seurat_obj,
  sample_name,
  output_dir = "Output_R/Find_vartoregress",
  technical_keywords = c("^MT-", "^RPL", "^RPS", "^IG[HKL]", "MALAT1", "NEAT1", "XIST"),
  max_pcs = 40,
  n_top_genes = 500,
  cutoff = 15
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data. Must have a reduction
  named "pca.log".

- sample_name:

  Character. The identifier for the biological sample, used for plot
  titles and file naming.

- output_dir:

  Character. Directory path where the generated JPEG will be saved.
  Default is "Output_R/Find_vartoregress".

- technical_keywords:

  Character vector. Regular expressions defining the "technical" genes
  to track. Default includes prefixes for mitochondrial ("^MT-"),
  ribosomal ("^RPL", "^RPS"), immunoglobulin (`"^IG[HKL]"`), and common
  lncRNAs ("MALAT1", "NEAT1", "XIST").

- max_pcs:

  Integer. The maximum number of principal components to evaluate.
  Default is 40.

- n_top_genes:

  Integer. The total number of top-loading genes to evaluate per PC.
  This is split evenly, meaning a value of 500 evaluates the top 250
  positive and top 250 negative features. Default is 500.

- cutoff:

  Numeric. The threshold percentage used to draw a dashed horizontal
  warning line on the plot. Default is 15.

## Value

A character vector containing the names of any samples that failed to
process (e.g., due to missing reductions). Returns an empty character
vector (`character(0)`) if successful.

## Note

This function explicitly looks for a dimensionality reduction named
`"pca.log"`. Ensure your Seurat object's PCA was saved under this name,
or modify the function to accept a custom reduction name.

## Examples

``` r
if (FALSE) { # \dontrun{
# Evaluate how much mitochondrial and ribosomal genes are driving your PCs
failed_qc <- plot_technical_contribution(
  seurat_obj = my_seurat,
  sample_name = "Viral_Infection_Cohort",
  technical_keywords = c("^MT-", "^RPL", "^RPS"),
  max_pcs = 30,
  n_top_genes = 400,
  cutoff = 10
)

if (length(failed_qc) == 0) {
  print("Technical contribution plot generated successfully!")
}
} # }
```
