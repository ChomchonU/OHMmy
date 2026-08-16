# Run Comprehensive Global Over-Representation Analysis (ORA)

This function acts as an end-to-end pipeline for performing
Over-Representation Analysis (ORA) on differential expression results
across multiple clusters. It automatically splits genes into upregulated
(Activated) and downregulated (Suppressed) sets based on user-defined
log2FC and p-value cutoffs. Using
[`clusterProfiler::enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html),
it tests against a provided term-to-gene database. The function
automatically generates and saves an extensive suite of visualizations:
per-cluster dotplots and barplots, global summary dotplots,
hierarchically clustered pathway dendrograms, and a unique "signed
significance" dotplot. If a cell type column is detected, it can also
split the signed significance plots by cell type.

## Usage

``` r
run_global_ora(
  ORA_df,
  m_t2g,
  output_dir,
  title_prefix = "Hallmark",
  top_n_global = 40,
  top_n_per_cluster = 10,
  log2fc_cutoff = 0.25,
  padj_cutoff = 0.05,
  indiv_width = 12,
  indiv_height = 8,
  global_width = 15,
  global_height = 15,
  variable_per_clus = FALSE
)
```

## Arguments

- ORA_df:

  A data frame containing differential expression results. Must contain
  columns `gene`, `cluster`, `avg_log2FC`, and `p_val_adj`. Optionally,
  a column containing "cell_type" or "celltype" in its name can be
  included for nested plotting.

- m_t2g:

  A two-column data frame mapping pathways/terms to genes (TERM2GENE
  format), typically sourced from MSigDB via the `msigdbr` package.

- output_dir:

  Character. Directory path where all plots and CSV summaries will be
  saved.

- title_prefix:

  Character. A prefix used for plot titles and file naming to identify
  the pathway database (e.g., "Hallmark", "KEGG"). Default is
  "Hallmark".

- top_n_global:

  Integer. The maximum number of top pathways to display on the global
  summary plots. Default is 40.

- top_n_per_cluster:

  Integer. The maximum number of top pathways to display per direction
  on the individual cluster plots. Default is 10.

- log2fc_cutoff:

  Numeric. The minimum absolute log2 fold change required to include a
  gene in the ORA test. Default is 0.25.

- padj_cutoff:

  Numeric. The maximum adjusted p-value required to include a gene in
  the ORA test. Default is 0.05.

- indiv_width:

  Numeric. The width (in inches) of the per-cluster plots. Default is
  12.

- indiv_height:

  Numeric. The height (in inches) of the per-cluster plots. Default is
  8.

- global_width:

  Numeric. The width (in inches) of the global summary plots. Default is
  15.

- global_height:

  Numeric. The height (in inches) of the global summary plots. Default
  is 15.

- variable_per_clus:

  Logical. If `TRUE` and a cell type column is detected in `ORA_df`, it
  generates separate signed significance dotplots for each cell type.
  Default is FALSE.

## Value

Invisibly returns a `tibble` (`combined_df`) containing the concatenated
enrichment results across all clusters and directions. Outputs multiple
JPEG plots and CSV tables as side effects to the specified `output_dir`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fetch MSigDB Hallmark gene sets
library(msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
m_t2g <- m_df[, c("gs_name", "gene_symbol")]

# Assume 'find_markers_output' is a dataframe resulting from Seurat's FindAllMarkers
# Ensure it has 'gene', 'cluster', 'avg_log2FC', and 'p_val_adj' columns.

global_ora_results <- run_global_ora(
  ORA_df = find_markers_output,
  m_t2g = m_t2g,
  output_dir = "Results/Pathways_ORA/",
  title_prefix = "Hallmark",
  top_n_global = 30,
  top_n_per_cluster = 15,
  log2fc_cutoff = 0.5,
  variable_per_clus = TRUE
)
} # }
```
