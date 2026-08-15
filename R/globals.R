#' @import ggplot2
#' @import dplyr
#' @import Seurat
NULL

utils::globalVariables(c(
  "pct.1", "pct.2", "avg_log2FC", "abs_log2FC", "abs_diff", "p_val_adj",
  "cluster", "score", "FastMNNIntegration", "RPCAIntegration",
  "CCAIntegration", "HarmonyIntegration", "id", "features.plot", "avg.exp",
  "pct.exp", "Gene", "AvgExpression", "Expression_Level", "lower_bound",
  "upper_bound", "Sample", "Condition", "CellType", "Count", "Proportion",
  "facet_max", "..p..", "Batch", "Freq", "prop_batch", "prop_level",
  "group", "feature", "avg.exp.scaled", "pct_exp", "avg_exp", "y", "x",
  "yend", "xend", "label", "var", "cluster_var", "gene", "y_id", "avg",
  ".data", "Correlation", "..p.format..", "gene_cluster",
  "gene_cluster_label", "gene_clean", "cluster_num", "NTopGenes",
  "WeightedPercentTechnical", "Direction", "color_flag", "PC",
  "PercentTechnical", "progress_bar", "core_enrichment", "setSize",
  "NES", "Description", "neg_log10_padj", "GeneRatio", "Signed_log10_padj",
  "mean_sig", "Abs_NES", ".", "n_clusters", "mean_neg_logp", "DIM1", "DIM2",
  "g1", "g2", "Cluster"
))
