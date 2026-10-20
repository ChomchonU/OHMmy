# Changelog

## OHMmy (development version)

### New features

- New
  [`analyze_and_save_clusters()`](https://chomchonu.github.io/OHMmy/reference/analyze_and_save_clusters.md)
  compares clustering methods at the same number of clusters instead of
  the same resolution value. It finds the cluster numbers that every
  method produced in its resolution sweep, plots the size of the
  smallest cluster against the number of clusters (saved as PNG), and
  returns, for a target number, the highest resolution of each method
  that gives it. Compared with the first draft of this function:
  - methods are named for any graph (`RNA_snn_res.0.5` -\> `RNA_snn`);
    before, only the `.log_snn` / `.sct_snn` naming worked and default
    Seurat names always failed with “no single number of clusters”;
  - empty factor levels (e.g. after
    [`subset()`](https://rdrr.io/r/base/subset.html)) no longer give a
    smallest cluster of 0, and columns such as `x_res.1_celltype` are no
    longer taken as resolutions;
  - new optional `pipelines` (compare only some methods) and
    `min_cells_line` (the reference line, default 30);
  - the plot is printed only in interactive sessions (no stray
    `Rplots.pdf`) and is returned as `attr(res, "plot")`, with the
    plotted values in `attr(res, "plot_data")` and the file path in
    `attr(res, "file")`;
  - plain-ASCII messages, documentation and tests. For the OHMmy naming
    (`log.fastmnn.log_snn_res.1`, `sct.fastmnn.sct_snn_res.1`), the
    returned table is unchanged.
- New
  [`cluster_dot_dendro()`](https://chomchonu.github.io/OHMmy/reference/cluster_dot_dendro.md)
  cuts the group and/or gene dendrograms of
  [`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
  into clusters. The default `"inconsistency"` method and the
  `"dynamic"` method (dynamicTreeCut) cut each branch at its own height,
  so clusters at different heights are found without choosing a single
  `k` or `h`; `"fixed"` performs an ordinary cut. Settings can be chosen
  automatically by silhouette width, unclear leaves stay unassigned (or
  are assigned with `assign_all = TRUE`), and group clusters can be
  written back to the cells of a Seurat object.
- [`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
  attaches the matrix and `hclust` trees behind its dendrograms to its
  returned data frame as attributes.
- [`ProcessSeuratLOG()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratLOG.md)
  and
  [`ProcessSeuratSCT()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratSCT.md)
  gain `force_norm` and `force_scale` (both `TRUE` by default, which
  keeps the previous behaviour). With `force_norm = FALSE`, an existing
  normalised `data` layer (LOG) or SCT assay (SCT) is reused instead of
  being recomputed; with `force_scale = FALSE`, an existing `scale.data`
  layer is reused when it covers the variable features. Missing layers
  are always computed.

### Bug fixes

- [`plot_split_dotplots_by_gene_cluster()`](https://chomchonu.github.io/OHMmy/reference/plot_split_dotplots_by_gene_cluster.md):
  - non-numeric cluster labels (e.g. cell type names) no longer lose the
    gene ordering or produce coercion warnings; they are labelled
    `gene (<cluster>)`;
  - `n_splits = 1` no longer fails, and no empty plots are written when
    there are fewer gene-cluster combinations than `n_splits`;
  - a single condition no longer fails (dot plots are saved without a
    dendrogram);
  - duplicated gene x cluster x condition rows are averaged instead of
    failing;
  - rows with a missing `id_col` value are dropped instead of failing.
- [`plot_gene_markers_with_dendro()`](https://chomchonu.github.io/OHMmy/reference/plot_gene_markers_with_dendro.md):
  - returns `NULL` with a message when no gene matches, instead of
    failing;
  - with a single condition, saves and returns the dot plot without a
    dendrogram;
  - a single string that is not a gene name is now always used as a
    regular expression (previously only strings containing `^`, `$` or
    `|` were, so e.g. `"GZM"` or `"HLA-D[PQR]"` matched nothing);
  - the dendrogram leaves now line up exactly with the dot-plot rows in
    the combined figure;
  - an invalid `dendro_side` gives an error instead of silently using
    `"right"`.
- Both functions give a clear error naming any missing required column.
  Output for inputs that already worked is unchanged.
- `FindTopMarkersAndHeatmap(use_sct = TRUE)` now works on SCT objects
  that were subset after
  [`SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html)/[`PrepSCTFindMarkers()`](https://satijalab.org/seurat/reference/PrepSCTFindMarkers.html).
  In that situation
  [`PrepSCTFindMarkers()`](https://satijalab.org/seurat/reference/PrepSCTFindMarkers.html)
  skips re-correction and every
  [`FindMarkers()`](https://satijalab.org/seurat/reference/FindMarkers.html)
  test failed with “Object contains multiple models with unequal library
  sizes”, followed by “object ‘pct.1’ not found”. The function now
  detects this and uses the existing corrected counts
  (`recorrect_umi = FALSE`), as Seurat recommends for subsets.
- [`FindTopMarkersAndHeatmap()`](https://chomchonu.github.io/OHMmy/reference/FindTopMarkersAndHeatmap.md)
  warns and returns without writing files when Seurat finds no markers
  at all, instead of failing with “object ‘pct.1’ not found”. When
  markers are found but none passes the thresholds, it now warns,
  returns `heatmap = NULL` and still writes the marker tables
  (previously it failed with “‘names’ attribute \[3\] must be the same
  length as the vector \[0\]”).

#### Ambient RNA (SoupX)

- [`run_soupx_post_clustering()`](https://chomchonu.github.io/OHMmy/reference/run_soupx_post_clustering.md)
  now builds the final Seurat objects from the counts corrected with the
  requested `methods` and `pCut`. Before, these were computed and then
  discarded, and the final objects always used SoupX’s default method.
  [`create_final_seurat()`](https://chomchonu.github.io/OHMmy/reference/create_final_seurat.md)
  gains an optional `corrected_list` argument for this (default `NULL`,
  which keeps the previous behaviour).
- [`estimate_contamination()`](https://chomchonu.github.io/OHMmy/reference/estimate_contamination.md):
  if
  [`SoupX::autoEstCont()`](https://rdrr.io/pkg/SoupX/man/autoEstCont.html)
  fails for one sample, a warning is issued and that sample continues
  with an estimated `rho` of 0, instead of stopping the whole run. A
  final contamination fraction above 0.5 (from a large `mulFac`) is
  applied with a warning instead of being refused by SoupX; above 1 is a
  clear error.
- [`load_counts()`](https://chomchonu.github.io/OHMmy/reference/load_counts.md)
  also reads Cell Ranger v2 directories (`genes.tsv`) and uncompressed
  files.
- [`create_seurat_for_clustering()`](https://chomchonu.github.io/OHMmy/reference/create_seurat_for_clustering.md)
  (and so
  [`prepare_soupx_inputs()`](https://chomchonu.github.io/OHMmy/reference/prepare_soupx_inputs.md))
  no longer fails on samples with 50 or fewer cells or variable
  features: the number of PCs and UMAP neighbours is capped to what the
  sample supports.

#### Clustering

- [`ClusterAndUMAP()`](https://chomchonu.github.io/OHMmy/reference/ClusterAndUMAP.md):
  - `algorithm` is now also used for the final clustering (it was only
    used for the resolution sweep);
  - when the UMAP already exists, `final_resolution` is now made the
    active clustering (`seurat_clusters` and `Idents`); before,
    whichever resolution the sweep ran last stayed active;
  - resolutions with different numbers of decimals
    (e.g. `c(0.25, 0.5, 1)`) are now matched to the column names Seurat
    creates (`RNA_snn_res.0.5`, not `RNA_snn_res.0.50`), so they are no
    longer re-computed on every run.

#### Statistics

- [`plot_metadata_stats()`](https://chomchonu.github.io/OHMmy/reference/plot_metadata_stats.md)
  counts every sample once. When a metadata value varies within a
  sample, the sample was counted once per distinct value; it now uses
  the first value and warns.
- [`plot_cell_abundance()`](https://chomchonu.github.io/OHMmy/reference/plot_cell_abundance.md):
  when a sample appears under more than one condition, its proportions
  are computed within each sample x condition (so they sum to 1 in
  each), with a warning. Before, they were split across conditions.

#### Visualisation and diagnostics

- [`make_cluster_palette()`](https://chomchonu.github.io/OHMmy/reference/make_cluster_palette.md):
  the `"kelly"` palette no longer repeats colours or includes a
  near-white colour above 20 clusters, and `"polychrome"` no longer
  repeats colours above 17 clusters. Colours for smaller numbers of
  clusters are unchanged.
- [`plot_gene_pair_correlations()`](https://chomchonu.github.io/OHMmy/reference/plot_gene_pair_correlations.md)
  gains `min_cells` (default 3). Clusters with fewer cells get an `NA`
  correlation; before, a cluster with two cells always showed a
  correlation of +/-1.
- [`plot_blend_nebulosa()`](https://chomchonu.github.io/OHMmy/reference/plot_blend_nebulosa.md)
  works on Seurat v5 objects with split layers (it failed with
  “GetAssayData doesn’t work for multiple layers”), and a numeric
  `min_cutoff` is now applied as an expression value (it was silently
  ignored).
- [`plot_blend_nebulosa()`](https://chomchonu.github.io/OHMmy/reference/plot_blend_nebulosa.md)
  no longer fails with “The `slot` argument of
  [`FetchData()`](https://satijalab.github.io/seurat-object/reference/FetchData.html)
  was deprecated in SeuratObject 5.0.0 and is now defunct” when an older
  Nebulosa (e.g. 1.16) is used with a newer SeuratObject (e.g. 5.4). Its
  density panels now use the same version-independent route as
  [`plot_combined()`](https://chomchonu.github.io/OHMmy/reference/plot_combined.md)
  (expression fetched with `layer = "data"`, plotted through Nebulosa’s
  SingleCellExperiment method), so it now also needs the suggested
  package SingleCellExperiment. Plots are unchanged.
- [`plot_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_technical_contribution.md)
  and
  [`plot_stacked_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_stacked_technical_contribution.md)
  no longer fail when `n_top_genes` / `gene_depths` ask for more genes
  than the PCA loadings contain.
- [`plot_pc_metadata_correlation()`](https://chomchonu.github.io/OHMmy/reference/plot_pc_metadata_correlation.md)
  skips non-numeric metadata columns with a warning (it failed with “‘y’
  must be numeric”).
- [`generate_and_save_heatmap()`](https://chomchonu.github.io/OHMmy/reference/generate_and_save_heatmap.md)
  drops genes with constant expression across the samples (with a
  message) instead of failing during row scaling.
- [`plot_dot_dendro_split()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_split.md)
  and
  [`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
  return `NULL` with a message when fewer than two genes or groups
  remain, instead of failing in
  [`hclust()`](https://rdrr.io/r/stats/hclust.html).
- For inputs that already worked, output is unchanged.

### Improvements

#### Analysis quality

- [`run_global_gsea()`](https://chomchonu.github.io/OHMmy/reference/run_global_gsea.md)
  gains `rank_by`. The default `"avg_log2FC"` ranks genes by fold change
  as before; `"signed_p"` ranks by `sign(avg_log2FC) * -log10(p)` (using
  `p_val`, or `p_val_adj` if absent), which combines direction and
  significance and gives far fewer tied ranks.
- [`run_global_ora()`](https://chomchonu.github.io/OHMmy/reference/run_global_ora.md)
  gains `universe`, the background gene set for the test. The default
  `NULL` keeps using the genes in `ORA_df`, which is only a valid
  background when the table contains every tested gene, not just the
  significant ones.
- [`FindTopMarkersAndHeatmap()`](https://chomchonu.github.io/OHMmy/reference/FindTopMarkersAndHeatmap.md)
  now shows every top marker in the heatmap. Markers missing from
  `scale.data` (non-variable genes) were silently dropped by
  [`DoHeatmap()`](https://satijalab.org/seurat/reference/DoHeatmap.html);
  they are now scaled for the heatmap only (Pearson residuals for SCT),
  with a message. Genes that were already shown keep the same values.

#### Robustness

- [`CleanSeuratReductions()`](https://chomchonu.github.io/OHMmy/reference/CleanSeuratReductions.md)
  no longer overwrites a reduction that already has the clean name
  (e.g. renaming `pca.log` when `pcaLog` exists). It keeps both and
  warns.

#### Plots are returned consistently

Every plotting function now gives you the plots it saves, without
changing what it returned before:

- [`plot_split_dotplots_by_gene_cluster()`](https://chomchonu.github.io/OHMmy/reference/plot_split_dotplots_by_gene_cluster.md)
  returns `list(dendrogram, dotplots)`, and
  [`plot_violin_qc_single()`](https://chomchonu.github.io/OHMmy/reference/plot_violin_qc_single.md)
  returns `list(qc, genes)`. Both returned `NULL` before.
- Functions that return data or file paths attach their plots as an
  attribute:
  - `attr(x, "plots")` on the results of
    [`run_global_ora()`](https://chomchonu.github.io/OHMmy/reference/run_global_ora.md)
    and
    [`run_global_gsea()`](https://chomchonu.github.io/OHMmy/reference/run_global_gsea.md)
    (a named list of every saved plot) and of
    [`plot_dot_dendro_multi()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_multi.md)
    (`combined`, `dotplot`, `row_dendrogram`, `col_dendrogram`);
  - `attr(x, "plot")` on the results of
    [`plot_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_technical_contribution.md),
    [`plot_stacked_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_stacked_technical_contribution.md),
    [`plot_pc_metadata_correlation()`](https://chomchonu.github.io/OHMmy/reference/plot_pc_metadata_correlation.md)
    (a `pheatmap` object, plus the correlation matrix as
    `attr(x, "cor_matrix")`) and
    [`generate_and_save_heatmap()`](https://chomchonu.github.io/OHMmy/reference/generate_and_save_heatmap.md)
    (a `gtable`).
- [`plot_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_technical_contribution.md)
  now returns its value invisibly. Because the value now carries an
  attribute, compare it with `length(x) == 0` rather than
  `identical(x, character(0))`.

#### Maintenance

- [`plot_metadata_stats()`](https://chomchonu.github.io/OHMmy/reference/plot_metadata_stats.md)
  uses
  [`after_stat()`](https://ggplot2.tidyverse.org/reference/aes_eval.html)
  instead of the deprecated `..p.format..`; the labels are unchanged.
  [`plot_cell_abundance()`](https://chomchonu.github.io/OHMmy/reference/plot_cell_abundance.md)
  keeps `..p..`, which ggpubr (\>= 1.0) converts to
  [`after_stat()`](https://ggplot2.tidyverse.org/reference/aes_eval.html)
  itself (an
  [`after_stat()`](https://ggplot2.tidyverse.org/reference/aes_eval.html)
  version of that label is mangled by ggpubr and fails).
- [`.blend_feature_plot_v5()`](https://chomchonu.github.io/OHMmy/reference/dot-blend_feature_plot_v5.md)
  no longer falls back to Seurat’s deprecated `slot` argument. OHMmy now
  requires SeuratObject (\>= 5.0.0) and ggplot2 (\>= 3.3.0).

## OHMmy 1.0.2

### Bug fixes

- [`process_soupx_samples()`](https://chomchonu.github.io/OHMmy/reference/process_soupx_samples.md)
  no longer fails with “sc must be an object of type SoupChannel”; it
  now chains
  [`prepare_soupx_inputs()`](https://chomchonu.github.io/OHMmy/reference/prepare_soupx_inputs.md)
  and
  [`run_soupx_post_clustering()`](https://chomchonu.github.io/OHMmy/reference/run_soupx_post_clustering.md)
  and returns the same structure.
- [`run_soupx_post_clustering()`](https://chomchonu.github.io/OHMmy/reference/run_soupx_post_clustering.md)
  now honours `manual_contam`: when sample names are given, `multiFac`
  is added only to those samples (previously it was applied to all
  samples).
  [`estimate_contamination()`](https://chomchonu.github.io/OHMmy/reference/estimate_contamination.md)
  log messages now describe the mode correctly.
- [`plot_pc_metadata_correlation()`](https://chomchonu.github.io/OHMmy/reference/plot_pc_metadata_correlation.md)
  no longer writes an empty JPEG when another graphics device is open
  (knitr, RStudio plot pane). Metadata columns that do not exist are
  skipped with a warning.
- `score_consensus(method = "gmm")` now works without attaching mclust;
  before, the fit failed silently and fell back to `fixed_cut`.
- [`plot_dot_dendro()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro.md)
  and
  [`plot_dot_dendro_split()`](https://chomchonu.github.io/OHMmy/reference/plot_dot_dendro_split.md)
  create `output_dir` if it does not exist.
- [`ProcessSeuratLOG()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratLOG.md)’s
  default `integration_method` is now the string `"HarmonyIntegration"`
  (the previous default was a function and errored). FastMNN integration
  now uses `batch_col` instead of a hard-coded `batch` column.
- [`generate_volcano_trio()`](https://chomchonu.github.io/OHMmy/reference/generate_volcano_trio.md),
  [`generate_and_save_heatmap()`](https://chomchonu.github.io/OHMmy/reference/generate_and_save_heatmap.md),
  [`run_global_gsea()`](https://chomchonu.github.io/OHMmy/reference/run_global_gsea.md),
  [`plot_blend_nebulosa()`](https://chomchonu.github.io/OHMmy/reference/plot_blend_nebulosa.md)
  and
  [`plot_combined()`](https://chomchonu.github.io/OHMmy/reference/plot_combined.md)
  no longer depend on packages being attached by the user. Missing
  optional packages give an informative error.
- [`run_global_gsea()`](https://chomchonu.github.io/OHMmy/reference/run_global_gsea.md)
  and
  [`run_global_ora()`](https://chomchonu.github.io/OHMmy/reference/run_global_ora.md)
  accept `output_dir` with or without a trailing slash.
- `plot_cell_abundance(global_test = "anova")` and
  [`plot_metadata_stats()`](https://chomchonu.github.io/OHMmy/reference/plot_metadata_stats.md)
  with `continuous_test_n3 = "anova"` no longer fail when filtering the
  ANOVA results.
- [`plot_metadata_stats()`](https://chomchonu.github.io/OHMmy/reference/plot_metadata_stats.md)
  skips `metadata_vars` that are not in the metadata (previously it
  errored before reaching its skip check).
- [`ProcessSeuratLOG()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratLOG.md)
  and
  [`ProcessSeuratSCT()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratSCT.md)
  check `integration_method` before any computation and give a clear
  error for unsupported values.

### Improvements

- [`plot_cluster_distributions()`](https://chomchonu.github.io/OHMmy/reference/plot_cluster_distributions.md)
  defaults no longer point to a personal path:
  `cluster_col = "seurat_clusters"`, `batch_col = "orig.ident"`,
  `output_dir = "Plots_cluster_distributions"`.
- [`generate_and_save_heatmap()`](https://chomchonu.github.io/OHMmy/reference/generate_and_save_heatmap.md)
  also accepts a plain expression matrix and returns the saved file path
  invisibly.
  [`generate_dimheatmaps()`](https://chomchonu.github.io/OHMmy/reference/generate_dimheatmaps.md)
  returns the failed PC windows invisibly.
- [`plot_stacked_technical_contribution()`](https://chomchonu.github.io/OHMmy/reference/plot_stacked_technical_contribution.md)
  orders facets numerically (PC_1, PC_2, …).
- [`plot_violin_qc_single()`](https://chomchonu.github.io/OHMmy/reference/plot_violin_qc_single.md)
  uses `width` for both saved figures.
- [`ProcessSeuratLOG()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratLOG.md)
  and
  [`ProcessSeuratSCT()`](https://chomchonu.github.io/OHMmy/reference/ProcessSeuratSCT.md)
  share internal helpers; functions no longer call
  [`library()`](https://rdrr.io/r/base/library.html),
  [`require()`](https://rdrr.io/r/base/library.html) or
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  internally.
- Documentation corrected and expanded throughout; all vignettes
  rewritten to run end to end on the pbmc3k dataset.
