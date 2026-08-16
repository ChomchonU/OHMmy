# Evaluate Depth-Dependent Technical Contributions to Principal Components

An advanced quality control tool that tracks how the influence of
technical/nuisance genes changes at varying depths of principal
component loadings. Instead of checking a single fixed number of top
genes, it iterates across a sequence of depths (e.g., top 20, 40, ...,
500), calculating the weighted percentage of technical genes at each
step for both positive and negative directions. The output is a massive
faceted bar chart (one facet per PC) with color-coded flags indicating
when a loading direction exceeds the defined technical cutoff threshold.

## Usage

``` r
plot_stacked_technical_contribution(
  seurat_obj,
  sample_name,
  reduction = "pca.log",
  technical_keywords = c("^MT-", "^RPL", "^RPS", "^IG[HKL]", "MALAT1", "NEAT1", "XIST"),
  gene_depths = seq(0, 500, by = 20),
  max_pcs = 40,
  cutoff = 15,
  output_dir = "Output_R/Find_vartoregress",
  plot_width = 30,
  plot_height = 15,
  dpi = 300
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data.

- sample_name:

  Character. The identifier for the biological sample, used for plot
  titles and file naming.

- reduction:

  Character. The dimensional reduction to extract feature loadings from.
  Default is "pca.log".

- technical_keywords:

  Character vector. Regular expressions defining the "technical" genes
  to track. Default includes mitochondrial ("^MT-"), ribosomal ("^RPL",
  "^RPS"), immunoglobulin (`"^IG[HKL]"`), and specific lncRNAs
  ("MALAT1", "NEAT1", "XIST").

- gene_depths:

  Numeric vector. A sequence defining the varying numbers of top/bottom
  genes to evaluate sequentially. Default is `seq(0, 500, by = 20)`.

- max_pcs:

  Integer. The maximum number of principal components to evaluate and
  facet. Default is 40.

- cutoff:

  Numeric. The threshold percentage used to draw a warning line and
  trigger the color-coded flag (Above/Below). Default is 15.

- output_dir:

  Character. Directory path where the generated JPEG will be saved.
  Default is "Output_R/Find_vartoregress".

- plot_width:

  Numeric. The width of the saved JPEG in inches. Because this plot
  contains many facets, the default is large (30).

- plot_height:

  Numeric. The height of the saved JPEG in inches. Default is 15.

- dpi:

  Numeric. The resolution (dpi) of the saved JPEG. Default is 300.

## Value

Invisibly returns a character string containing the exact file path
where the plot was saved.

## Examples

``` r
if (FALSE) { # \dontrun{
# Track technical noise saturation across the first 30 PCs at intervals of 25 genes
plot_stacked_technical_contribution(
  seurat_obj = my_seurat,
  sample_name = "T_Cell_Subset",
  reduction = "pca",
  gene_depths = seq(0, 400, by = 25),
  max_pcs = 30,
  cutoff = 10,
  plot_width = 24,
  plot_height = 12
)
} # }
```
