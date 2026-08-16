# Generate and Save Dual-Layout Differential Expression Heatmaps

Automatically extracts the most biologically and statistically relevant
genes for a given comparison using
[`get_top_mixed_genes()`](https://chomchonu.github.io/OHMmy/reference/get_top_mixed_genes.md),
pulls their normalized expression data, and generates two side-by-side
`pheatmap` visualizations. The left heatmap forces a specific,
user-defined sample order (ideal for gradients or chronological data),
while the right heatmap allows free hierarchical clustering of samples
to reveal natural groupings. The function dynamically scales the height
of the output PNG based on the number of genes to ensure row labels
remain readable.

## Usage

``` r
generate_and_save_heatmap(
  res_obj,
  comp_name,
  comp_title,
  vsd_data,
  anno_col,
  ordered_samps,
  n_padj,
  n_lfc,
  out_dir,
  ts,
  clus = "all"
)
```

## Arguments

- res_obj:

  A data frame containing differential expression results (must contain
  `padj` and `log2FoldChange` columns).

- comp_name:

  Character. A filesystem-safe string representing the comparison, used
  for the output filename (e.g., "Infected_vs_Mock").

- comp_title:

  Character. A human-readable title displayed at the top of the left
  heatmap.

- vsd_data:

  A `SummarizedExperiment` object or matrix containing normalized
  expression data (e.g., the output of DESeq2's `vst()` or `rlog()`).

- anno_col:

  A data frame containing sample metadata for the heatmap annotations.
  Row names must match the column names of `vsd_data`.

- ordered_samps:

  Character vector. The exact order of sample IDs (column names) to be
  plotted in the unclustered (left) heatmap.

- n_padj:

  Integer. The number of top significant genes to extract based on
  lowest adjusted p-value.

- n_lfc:

  Integer. The number of top significant genes to extract based on
  highest absolute log2 fold change.

- out_dir:

  Character. Directory path where the generated PNG will be saved.

- ts:

  Character. A timestamp string appended to the filename for version
  control.

- clus:

  Character. An optional identifier (e.g., "CD8_T_Cells") used in the
  filename if looping across multiple subsets or clusters. Default is
  "all".

## Value

Invisibly returns `NULL`. The function is called for its side effect of
saving the combined plot to disk.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'deseq_res' is your results dataframe and 'vsd' is a DESeqTransform object

# Define annotation colors and metadata
sample_anno <- data.frame(Condition = colData(vsd)$Condition)
rownames(sample_anno) <- colnames(vsd)

# Specify the exact order you want samples to appear in the fixed heatmap
ordered_samples <- c("Mock_1", "Mock_2", "Infected_1", "Infected_2")

generate_and_save_heatmap(
  res_obj = deseq_res,
  comp_name = "Infection_Effect",
  comp_title = "Infected vs Mock (Global)",
  vsd_data = vsd,
  anno_col = sample_anno,
  ordered_samps = ordered_samples,
  n_padj = 40,
  n_lfc = 40,
  out_dir = "Results/Heatmaps",
  ts = format(Sys.time(), "%Y%m%d_%H%M%S")
)
} # }
```
