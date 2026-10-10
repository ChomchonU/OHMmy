# Changelog

## OHMmy (development version)

### Bug fixes

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
  at all, instead of failing with “object ‘pct.1’ not found”.

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
