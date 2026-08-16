# Calculate, Test, and Visualize Cell Type Abundances

Computes the fractional proportion of each cell type per sample from a
Seurat object's metadata. It automatically detects the number of
experimental conditions and applies the appropriate statistical
framework using `rstatix` (e.g., Kruskal-Wallis + Dunn's test for \>2
conditions, or Mann-Whitney for 2 conditions). It generates highly
customized, publication-ready boxplots with jittered points and
automatically positions significance brackets (`ggpubr`) dynamically per
facet to avoid overlapping.

## Usage

``` r
plot_cell_abundance(
  seurat_obj,
  sample_col,
  condition_col,
  celltype_col,
  global_test = "kruskal.test",
  strict_posthoc = TRUE,
  pairwise_test = "mann_whitney",
  p_adjust = "BH",
  facet_by_cluster = TRUE,
  facet_scales = "free_y",
  pairwise_label = "p.signif",
  y_expand = 0.2,
  output_dir = ".",
  base_size = 3.5,
  dpi = 300
)
```

## Arguments

- seurat_obj:

  A Seurat object containing single-cell data with populated metadata.

- sample_col:

  Character. The metadata column containing unique sample IDs (used to
  calculate independent proportions).

- condition_col:

  Character. The metadata column containing the experimental conditions
  or groups to compare.

- celltype_col:

  Character. The metadata column containing cell type or cluster labels.

- global_test:

  Character. The overarching statistical test to use if there are \>2
  conditions. Options: "kruskal.test" (uses Dunn's for post-hoc) or
  "anova" (uses Tukey's HSD). Default is "kruskal.test".

- strict_posthoc:

  Logical. If TRUE, post-hoc pairwise tests and brackets are ONLY
  generated for cell types that pass the global test significance
  threshold (p \< 0.05). Default is TRUE.

- pairwise_test:

  Character. The statistical test to use if there are exactly 2
  conditions. Options: "mann_whitney", "wilcoxon_paired", or "t_test".
  Default is "mann_whitney".

- p_adjust:

  Character. The multiple testing correction method to pass to `rstatix`
  (e.g., "BH", "bonferroni"). Default is "BH".

- facet_by_cluster:

  Logical. If TRUE, creates a faceted plot where each panel is a cell
  type. If FALSE, plots all cell types grouped on the x-axis. Default is
  TRUE.

- facet_scales:

  Character. The `scales` argument passed to `facet_wrap`. Default is
  "free_y".

- pairwise_label:

  Character. What to display on the significance brackets. Options:
  "p.adj" (numeric), "p.format" (raw p-value), or "p.signif" (stars).
  Default is "p.signif".

- y_expand:

  Numeric. A multiplier used to expand the top of the Y-axis to ensure
  significance brackets are not cut off. Default is 0.2.

- output_dir:

  Character. Directory path where the generated JPEG will be saved.
  Default is the current working directory (".").

- base_size:

  Numeric. A baseline metric used to dynamically calculate the width and
  height of the saved plot based on the number of panels. Default is
  3.5.

- dpi:

  Numeric. The resolution of the saved JPEG. Default is 300.

## Value

Returns the generated `ggplot` object. The function also saves the plot
to disk as a JPEG as a side effect.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming your Seurat object 'my_seurat' has metadata columns:
# "Patient_ID", "Disease_State" (e.g., Healthy, Mild, Severe), and "Cell_Subset"

abundance_plot <- plot_cell_abundance(
  seurat_obj = my_seurat,
  sample_col = "Patient_ID",
  condition_col = "Disease_State",
  celltype_col = "Cell_Subset",
  global_test = "kruskal.test",
  strict_posthoc = TRUE,
  p_adjust = "BH",
  pairwise_label = "p.signif",
  output_dir = "Results/Abundance"
)

# Display the plot in the R viewer
print(abundance_plot)
} # }
```
