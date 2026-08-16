# Extract Top Differentially Expressed Genes by Mixed Criteria

Filters a differential expression results data frame for strictly
significant genes (`padj < 0.05` and `abs(log2FoldChange) >= 1`). It
then sorts the remaining genes using two different metrics: statistical
confidence (lowest adjusted p-value) and effect size (highest absolute
log2 fold change). Finally, it returns the unique union of the top genes
from both sorting methods.

## Usage

``` r
get_top_mixed_genes(res_obj, n_padj = 40, n_lfc = 40)
```

## Arguments

- res_obj:

  A data frame containing differential expression results (e.g., from
  DESeq2). Must contain `padj` and `log2FoldChange` columns, with gene
  names set as the row names.

- n_padj:

  Integer. The number of top genes to extract based on statistical
  significance (lowest `padj`). Default is 40.

- n_lfc:

  Integer. The number of top genes to extract based on absolute effect
  size (highest `abs(log2FoldChange)`). Default is 40.

## Value

A character vector containing the unique combined list of top gene
names. Note that the total length of this vector may be less than
`n_padj + n_lfc` due to overlapping genes between the two lists.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'deseq_res' is a DESeqResults object or dataframe

# Extract a robust list of up to 80 top markers (40 by p-value, 40 by fold change)
top_mixed_markers <- get_top_mixed_genes(
  res_obj = deseq_res,
  n_padj = 40,
  n_lfc = 40
)

# Use this curated list for a targeted heatmap
pheatmap(assay(vsd)[top_mixed_markers, ])
} # }
```
