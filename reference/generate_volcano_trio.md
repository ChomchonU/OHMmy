# Generate a Trio of Volcano Plots for Differential Expression

Takes a differential expression results data frame (e.g., from DESeq2 or
Seurat) and generates a three-panel volcano plot using `EnhancedVolcano`
and `patchwork`. It applies strict significance thresholds
(`padj < 0.05` and `abs(log2FoldChange) >= 1.0`) to automatically
extract and label:

1.  The top 10 significantly upregulated and downregulated genes.

2.  A user-provided list of target functional genes (filtered to only
    show significant ones).

3.  A combined view of both the top 20 genes and the significant target
    genes.

## Usage

``` r
generate_volcano_trio(res, main_title, target_genes)
```

## Arguments

- res:

  A data frame containing differential expression results. Must contain
  row names as gene symbols, and the columns `padj` (adjusted p-value)
  and `log2FoldChange`.

- main_title:

  Character. The overarching title displayed at the very top of the
  combined plot.

- target_genes:

  Character vector. A specific list of genes of interest to highlight
  (e.g., functional pathway markers). Only genes in this list that meet
  the strict significance thresholds will be labeled.

## Value

A `patchwork` object containing three horizontally aligned
`ggplot`/volcano plot panels.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'de_results' is a data frame from DESeq2 or Seurat FindMarkers
# (Make sure your Seurat output column 'avg_log2FC' is renamed to 'log2FoldChange'
# and 'p_val_adj' to 'padj' before passing it to this function!)

# Define a list of functional genes you care about for your specific subset
nk_functional_genes <- c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG", "KLRK1")

# Generate the 3-panel plot
volcano_trio <- generate_volcano_trio(
  res = de_results,
  main_title = "Differential Expression: NK Cells (Infected vs Control)",
  target_genes = nk_functional_genes
)

# Display the plot
print(volcano_trio)

# Save the wide format plot
ggsave("Volcano_Trio.png", plot = volcano_trio, width = 24, height = 8, dpi = 300)
} # }
```
