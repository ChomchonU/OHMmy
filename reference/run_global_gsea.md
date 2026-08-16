# Run Comprehensive Global Gene Set Enrichment Analysis (GSEA)

An end-to-end pipeline for performing Gene Set Enrichment Analysis
(GSEA) on differential expression results across multiple clusters or
cell types. Unlike ORA, which relies on strict significance thresholds,
this function ranks all available genes by their average log2 fold
change and evaluates the entire spectrum using
[`clusterProfiler::GSEA`](https://rdrr.io/pkg/clusterProfiler/man/GSEA.html).
It performs a "pre-pass" calculation to establish a global Normalized
Enrichment Score (NES) range, ensuring all subsequent per-cluster
dotplots share a unified, directly comparable x-axis. It automatically
generates standard NES dotplots, classic GSEA enrichment ridge plots,
hierarchically clustered global dotplots, and signed significance plots.

## Usage

``` r
run_global_gsea(
  GSEA_df,
  m_t2g,
  output_dir,
  title_prefix = "Hallmark",
  top_n_per_direction = 10,
  padj_cutoff = 0.05,
  dotplot_width = 12,
  dotplot_height = 8,
  gseaplot_width = 8,
  gseaplot_height = 6,
  top_n_overall = 5,
  variable_per_clus = FALSE
)
```

## Arguments

- GSEA_df:

  A data frame containing differential expression results. Must contain
  columns `gene`, `cluster`, and `avg_log2FC`. Optionally, a column
  containing "cell_type" or "celltype" can be included for nested
  plotting.

- m_t2g:

  A two-column data frame mapping pathways/terms to genes (TERM2GENE
  format), typically sourced from MSigDB.

- output_dir:

  Character. Directory path where all plots and CSV summaries will be
  saved.

- title_prefix:

  Character. A prefix used for plot titles and file naming to identify
  the pathway database (e.g., "Hallmark", "KEGG"). Default is
  "Hallmark".

- top_n_per_direction:

  Integer. The maximum number of top activated and top suppressed
  pathways to display on individual cluster dotplots. Default is 10.

- padj_cutoff:

  Numeric. The adjusted p-value cutoff for statistical significance in
  the GSEA test. Default is 0.05.

- dotplot_width:

  Numeric. The width (in inches) of the per-cluster NES dotplots.
  Default is 12.

- dotplot_height:

  Numeric. The height (in inches) of the per-cluster NES dotplots.
  Default is 8.

- gseaplot_width:

  Numeric. The width (in inches) of the classic GSEA enrichment plot for
  the top pathway. Default is 8.

- gseaplot_height:

  Numeric. The height (in inches) of the classic GSEA enrichment plot
  for the top pathway. Default is 6.

- top_n_overall:

  Integer. The number of top pathways to extract from *each* cluster to
  build the combined global summary plots. Default is 5.

- variable_per_clus:

  Logical. A toggle to dictate specific nested behavior (retained for
  pipeline compatibility with the ORA function framework). Default is
  FALSE.

## Value

Invisibly returns a `tibble` (`combined_df`) containing the concatenated
GSEA results across all evaluated clusters. Outputs multiple JPEG plots
and CSV tables as side effects to the specified `output_dir`.

## Note

This function automatically registers a serial `BiocParallel` parameter
at the start of the run to prevent common parallel backend connection
errors on Windows environments.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fetch MSigDB Hallmark gene sets
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
m_t2g <- m_df[, c("gs_name", "gene_symbol")]

# Assume 'find_markers_output' is your DE results containing 'gene', 'cluster',
# and 'avg_log2FC'. (Do not pre-filter this dataframe for p-value thresholds,
# as GSEA requires the full ranked list of background genes!)

global_gsea_results <- run_global_gsea(
  GSEA_df = find_markers_output,
  m_t2g = m_t2g,
  output_dir = "Results/Pathways_GSEA/",
  title_prefix = "Hallmark",
  top_n_per_direction = 10,
  top_n_overall = 5
)
} # }
```
