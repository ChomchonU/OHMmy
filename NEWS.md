# OHMmy (development version)

## Bug fixes

* `FindTopMarkersAndHeatmap(use_sct = TRUE)` now works on SCT objects that were
  subset after `SCTransform()`/`PrepSCTFindMarkers()`. In that situation
  `PrepSCTFindMarkers()` skips re-correction and every `FindMarkers()` test failed
  with "Object contains multiple models with unequal library sizes", followed by
  "object 'pct.1' not found". The function now detects this and uses the existing
  corrected counts (`recorrect_umi = FALSE`), as Seurat recommends for subsets.
* `FindTopMarkersAndHeatmap()` warns and returns without writing files when Seurat
  finds no markers at all, instead of failing with "object 'pct.1' not found".

# OHMmy 1.0.2

## Bug fixes

* `process_soupx_samples()` no longer fails with "sc must be an object of type
  SoupChannel"; it now chains `prepare_soupx_inputs()` and
  `run_soupx_post_clustering()` and returns the same structure.
* `run_soupx_post_clustering()` now honours `manual_contam`: when sample names are
  given, `multiFac` is added only to those samples (previously it was applied to
  all samples). `estimate_contamination()` log messages now describe the mode
  correctly.
* `plot_pc_metadata_correlation()` no longer writes an empty JPEG when another
  graphics device is open (knitr, RStudio plot pane). Metadata columns that do
  not exist are skipped with a warning.
* `score_consensus(method = "gmm")` now works without attaching mclust; before,
  the fit failed silently and fell back to `fixed_cut`.
* `plot_dot_dendro()` and `plot_dot_dendro_split()` create `output_dir` if it
  does not exist.
* `ProcessSeuratLOG()`'s default `integration_method` is now the string
  `"HarmonyIntegration"` (the previous default was a function and errored).
  FastMNN integration now uses `batch_col` instead of a hard-coded `batch` column.
* `generate_volcano_trio()`, `generate_and_save_heatmap()`, `run_global_gsea()`,
  `plot_blend_nebulosa()` and `plot_combined()` no longer depend on packages being
  attached by the user. Missing optional packages give an informative error.
* `run_global_gsea()` and `run_global_ora()` accept `output_dir` with or without
  a trailing slash.
* `plot_cell_abundance(global_test = "anova")` and `plot_metadata_stats()` with
  `continuous_test_n3 = "anova"` no longer fail when filtering the ANOVA results.
* `plot_metadata_stats()` skips `metadata_vars` that are not in the metadata
  (previously it errored before reaching its skip check).
* `ProcessSeuratLOG()` and `ProcessSeuratSCT()` check `integration_method` before
  any computation and give a clear error for unsupported values.

## Improvements

* `plot_cluster_distributions()` defaults no longer point to a personal path:
  `cluster_col = "seurat_clusters"`, `batch_col = "orig.ident"`,
  `output_dir = "Plots_cluster_distributions"`.
* `generate_and_save_heatmap()` also accepts a plain expression matrix and returns
  the saved file path invisibly. `generate_dimheatmaps()` returns the failed PC
  windows invisibly.
* `plot_stacked_technical_contribution()` orders facets numerically (PC_1, PC_2, ...).
* `plot_violin_qc_single()` uses `width` for both saved figures.
* `ProcessSeuratLOG()` and `ProcessSeuratSCT()` share internal helpers; functions no
  longer call `library()`, `require()` or `install.packages()` internally.
* Documentation corrected and expanded throughout; all vignettes rewritten to run
  end to end on the pbmc3k dataset.
