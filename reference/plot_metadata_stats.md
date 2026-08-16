# Evaluate and Visualize Sample-Level Metadata Statistics

Automatically computes statistical differences and generates plots for
sample-level metadata covariates (e.g., Age, Gender, Clinical Scores)
across experimental conditions. Crucially, it deduplicates the Seurat
metadata down to the unique sample level prior to testing, preventing
false discoveries caused by single-cell pseudoreplication. The function
automatically routes continuous variables to boxplots (using parametric
or non-parametric tests like ANOVA/Kruskal-Wallis) and categorical
variables to stacked proportional bar charts (using Fisher's Exact or
Chi-Square tests).

## Usage

``` r
plot_metadata_stats(
  seurat_obj,
  sample_col = "Sample",
  condition_col = "Severity",
  metadata_vars = c("Age", "Gender"),
  continuous_test_n2 = "mann_whitney",
  continuous_test_n3 = "kruskal.test",
  categorical_test = "chisq",
  strict_posthoc = TRUE,
  p_adjust = "BH",
  add_facet = NULL,
  output_dir = "metadata_plots",
  plot_width = 6,
  plot_height = 5,
  dpi = 300
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data with populated metadata.

- sample_col:

  Character. The metadata column identifying unique biological samples.
  Used to deduplicate the data so that N equals the number of
  patients/samples, not the number of cells. Default is "Sample".

- condition_col:

  Character. The metadata column defining the experimental groups to
  compare (e.g., "Severity", "Treatment"). Default is "Severity".

- metadata_vars:

  Character vector. The specific metadata columns to evaluate (e.g.,
  `c("Age", "Gender")`). The function will automatically detect if each
  is continuous or categorical.

- continuous_test_n2:

  Character. The statistical test for continuous variables when there
  are exactly 2 conditions. Options: "mann_whitney" or "t_test". Default
  is "mann_whitney".

- continuous_test_n3:

  Character. The global statistical test for continuous variables when
  there are \>2 conditions. Options: "kruskal.test" or "anova". Default
  is "kruskal.test".

- categorical_test:

  Character. The statistical test for categorical variables. Options:
  "fisher" (recommended for imbalanced/small N) or "chisq". Default is
  "chisq".

- strict_posthoc:

  Logical. If TRUE, pairwise post-hoc tests (and plot brackets) are only
  executed if the global test (e.g., ANOVA/Kruskal) is significant
  (`p < 0.05`). Default is TRUE.

- p_adjust:

  Character. The multiple testing correction method for post-hoc
  pairwise comparisons (e.g., "BH", "bonferroni"). Default is "BH".

- add_facet:

  Character. An optional metadata column to facet the plots and stratify
  the statistical tests by (e.g., faceting by "Tissue"). Default is
  NULL.

- output_dir:

  Character. Directory path where the generated JPEGs will be saved.
  Default is "metadata_plots".

- plot_width:

  Numeric. The width of the saved JPEGs in inches. Default is 6.

- plot_height:

  Numeric. The height of the saved JPEGs in inches. Default is 5.

- dpi:

  Numeric. The resolution of the saved JPEGs. Default is 300.

## Value

A list containing two elements:

- `plots`: A named list of the generated `ggplot` objects.

- `stats`: A named list of statistical results, where each element
  contains data frames for both the `Global` test and the `PostHoc`
  pairwise tests.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming your metadata has sample-level metrics like Age (numeric)
# and Sex (categorical/character)

metadata_results <- plot_metadata_stats(
  seurat_obj = my_seurat,
  sample_col = "PatientID",
  condition_col = "Disease_Status",
  metadata_vars = c("Age", "Sex", "BMI", "Smoking_Status"),
  continuous_test_n3 = "anova",  # Use parametric ANOVA for continuous
  categorical_test = "fisher",   # Use Fisher's exact for categorical
  strict_posthoc = TRUE,
  output_dir = "Results/Demographics"
)

# Extract and view the exact statistical p-values for Age
print(metadata_results$stats[["Age"]]$Global)
print(metadata_results$stats[["Age"]]$PostHoc)

# View the generated plot for Sex within R
print(metadata_results$plots[["Sex"]])
} # }
```
